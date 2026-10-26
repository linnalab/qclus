"""
fraction_unspliced_from_bam with a GTF file, which counts molecules with velocyto's rules.

Each rule has a barcode of its own in one small BAM file. test_counts_match_velocyto compares the counts on a
larger random BAM file (synthetic_bam.py) with those of velocyto 0.17.17's own code.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import pytest

import qclus.utils as utils
from qclus import velocyto_rules
from synthetic_bam import make_dataset

DATA = Path(__file__).parent / "data"

# Gene models (1-based, inclusive). GENE1's two transcripts differ in their second exon, so that 3201-3400 is
# exonic in T2 and intronic in T1. GENE3 lies inside GENE1's first intron, on the same strand. velocyto never
# searches the last feature of a chromosome strand, so the LAST genes keep that quirk away from the others.
GTF = [
    ("chr1", "+", "GENE1", "T1", [(1001, 1200), (3001, 3200), (5001, 5200)]),
    ("chr1", "+", "GENE1", "T2", [(1001, 1200), (3001, 3400), (5001, 5200)]),
    ("chr1", "-", "GENE2", "T3", [(6001, 6200), (8001, 8300)]),
    ("chr1", "+", "GENE3", "T4", [(2001, 2100), (2501, 2600)]),
    ("chr1", "+", "LAST1", "T5", [(15001, 15100)]),
    ("chr1", "-", "LAST2", "T6", [(15001, 15100)]),
    ("chrM", "+", "MT-GENE", "T7", [(1, 1000)]),
    ("chrM", "+", "MT-LAST", "T8", [(2001, 2100)]),
]
# Repeats: one on each strand of GENE1's first intron, at 1501-1700 (+) and 2701-2900 (-). The others are there
# because of how velocyto joins repeats: the first repeat of a chromosome strand is joined to the last interval of
# the one before, the last repeat moves on to the next, the last interval is never searched, and the repeats of
# the last chromosome strand are never stored.
MASK = [
    ("chr1", "+", 1, 10), ("chr1", "+", 1501, 1700), ("chr1", "+", 18001, 18100), ("chr1", "+", 19001, 19100),
    ("chr1", "-", 1, 10), ("chr1", "-", 2701, 2900), ("chr1", "-", 18001, 18100), ("chr1", "-", 19001, 19100),
    ("chrX", "+", 1, 100),
]


def barcode(i):
    return "".join("ACGT"[(i >> (2 * k)) & 3] for k in range(16))


# One barcode per rule: name -> list of reads (start 1-based, CIGAR, reverse, tags other than CB, contig)
RULES = {
    "one_molecule_three_reads": [(1021, "90M", False, {}), (1051, "90M", False, {}), (1081, "90M", False, {})],
    "two_molecules": [(1021, "90M", False, {"UB": "UMI2"}), (1051, "90M", False, {})],
    "exonic_and_intronic_reads": [(1021, "90M", False, {}), (1301, "90M", False, {})],
    "intronic": [(1801, "90M", False, {})],
    "across_exon_end": [(1151, "90M", False, {})],
    "ambiguous": [(3251, "90M", False, {})],
    "spliced_read": [(1151, "50M1800N40M", False, {})],
    "multimapper": [(1021, "90M", False, {"NH": 2})],
    "antisense": [(1021, "90M", True, {})],
    "without_umi": [(1021, "90M", False, {"UB": None})],
    "in_repeat_on_its_strand": [(1551, "90M", False, {})],
    "in_repeat_on_other_strand": [(2751, "90M", False, {})],
    "two_genes_one_strand": [(2001, "90M", False, {})],
    "mitochondrial": [(101, "90M", False, {"contig": "chrM"})],
    # spliced 1, unspliced 1, ambiguous 2: the fraction is 1 / 4
    "all_classes": [
        (1021, "90M", False, {"UB": "UMI1"}),
        (1801, "90M", False, {"UB": "UMI2"}),
        (3251, "90M", False, {"UB": "UMI3"}),
        (3301, "90M", False, {"UB": "UMI4"}),
    ],
}
BARCODES = {name: barcode(i + 1) for i, name in enumerate(RULES)}


def write_case(directory, genes, repeats, reads, contigs, barcodes):
    """
    Write a BAM file with its gene models, repeats and barcodes.

    genes: (contig, strand, gene, transcript, exons); repeats: (contig, strand, start, end); reads: (contig, start
    (1-based), CIGAR, reverse, tags); contigs: (name, length).
    """
    directory = Path(directory)
    lines = []
    for contig, strand, gene, transcript, exons in genes:
        for number, (start, end) in enumerate(exons if strand == "+" else exons[::-1], start=1):
            tags = f'gene_id "{gene}"; transcript_id "{transcript}"; exon_number "{number}";'
            lines.append(f"{contig}\ttest\texon\t{start}\t{end}\t.\t{strand}\t.\t{tags}")
    (directory / "genes.gtf").write_text("\n".join(lines) + "\n")
    (directory / "rmsk.gtf").write_text("".join(
        f'{contig}\trmsk\texon\t{start}\t{end}\t.\t{strand}\t.\tgene_id "R{i}"; transcript_id "R{i}";\n'
        for i, (contig, strand, start, end) in enumerate(repeats)
    ))

    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": length} for name, length in contigs],
    })
    order = {name: i for i, (name, _) in enumerate(contigs)}
    path = directory / "possorted_genome_bam.bam"
    with pysam.AlignmentFile(str(path), "wb", header=header) as handle:
        for i, (contig, start, cigar, reverse, tags) in enumerate(sorted(reads, key=lambda r: (order[r[0]], r[1]))):
            read = pysam.AlignedSegment(header)
            read.query_name = f"read{i}"
            read.reference_id = header.get_tid(contig)
            read.reference_start = start - 1
            read.cigarstring = cigar
            length = read.query_length or read.infer_query_length()
            read.query_sequence = "A" * length
            read.query_qualities = pysam.qualitystring_to_array("I" * length)
            read.flag = 16 if reverse else 0
            read.mapping_quality = 255
            for tag, value in tags.items():
                if value is not None:
                    read.set_tag(tag, value)
            handle.write(read)
    pysam.index(str(path))
    (directory / "barcodes.tsv").write_text("".join(f"{b}-1\n" for b in barcodes))
    return {name: str(directory / file) for name, file in (
        ("bam", "possorted_genome_bam.bam"), ("barcodes", "barcodes.tsv"), ("gtf", "genes.gtf"), ("mask", "rmsk.gtf"))}


@pytest.fixture(scope="module")
def rules_files(tmp_path_factory):
    reads = []
    for name, rule_reads in RULES.items():
        for start, cigar, reverse, extra in rule_reads:
            tags = {"CB": BARCODES[name] + "-1", "UB": "UMI1", "NH": 1}
            tags.update({key: value for key, value in extra.items() if key != "contig"})
            reads.append((extra.get("contig", "chr1"), start, cigar, reverse, tags))
    contigs = [("chr1", 20_000), ("chrM", 16_569), ("chrX", 1_000)]
    return write_case(tmp_path_factory.mktemp("rules"), GTF, MASK, reads, contigs, BARCODES.values())


@pytest.fixture(scope="module")
def counts(rules_files):
    return utils.fraction_unspliced_from_bam(
        bam_path=rules_files["bam"], barcodes_path=rules_files["barcodes"], gtf_path=rules_files["gtf"],
        repeat_mask_path=rules_files["mask"], cores=1,
    )


def molecules(counts, rule):
    """(spliced, unspliced, ambiguous) of a rule's barcode, or None if the barcode has no read velocyto reads."""
    if BARCODES[rule] not in counts.index:
        return None
    return tuple(int(x) for x in counts.loc[BARCODES[rule], list(velocyto_rules.LAYERS)])


def test_columns(counts):
    assert list(counts.columns) == ["spliced", "unspliced", "ambiguous", "fraction_unspliced"]


def test_reads_of_one_molecule_count_once(counts):
    assert molecules(counts, "one_molecule_three_reads") == (1, 0, 0)
    assert molecules(counts, "two_molecules") == (2, 0, 0)


def test_molecule_with_exonic_and_intronic_reads_is_unspliced(counts):
    assert molecules(counts, "exonic_and_intronic_reads") == (0, 1, 0)
    assert molecules(counts, "intronic") == (0, 1, 0)


def test_read_across_an_exon_end_is_unspliced(counts):
    assert molecules(counts, "across_exon_end") == (0, 1, 0)


def test_exonic_in_one_transcript_and_intronic_in_another_is_ambiguous(counts):
    assert molecules(counts, "ambiguous") == (0, 0, 1)


def test_spliced_read(counts):
    assert molecules(counts, "spliced_read") == (1, 0, 0)


def test_multimapper_is_not_read(counts):
    assert molecules(counts, "multimapper") is None


def test_read_without_umi_is_not_read(counts):
    assert molecules(counts, "without_umi") is None


def test_antisense_read_is_not_counted(counts):
    assert molecules(counts, "antisense") == (0, 0, 0)


def test_repeat_masks_reads_on_its_own_strand_only(counts):
    assert molecules(counts, "in_repeat_on_its_strand") == (0, 0, 0)
    assert molecules(counts, "in_repeat_on_other_strand") == (0, 1, 0)


def test_molecule_compatible_with_two_genes_is_not_counted(counts):
    assert molecules(counts, "two_genes_one_strand") == (0, 0, 0)


def test_mitochondrial_reads_are_not_counted(counts):
    # velocyto names chrM 'MT' for reads and 'M' for genes, so the two never meet
    assert molecules(counts, "mitochondrial") == (0, 0, 0)
    assert counts.loc[BARCODES["mitochondrial"], "fraction_unspliced"] == 0


def test_fraction_has_ambiguous_molecules_in_the_denominator(counts):
    assert molecules(counts, "all_classes") == (1, 1, 2)
    assert counts.loc[BARCODES["all_classes"], "fraction_unspliced"] == pytest.approx(0.25)


def test_without_the_mask_reads_in_repeats_count(rules_files):
    result = utils.fraction_unspliced_from_bam(
        bam_path=rules_files["bam"], barcodes_path=rules_files["barcodes"], gtf_path=rules_files["gtf"], cores=1)
    assert molecules(result, "in_repeat_on_its_strand") == (0, 1, 0)


@pytest.fixture(scope="module")
def synthetic(tmp_path_factory):
    return make_dataset(tmp_path_factory.mktemp("synthetic"))


@pytest.mark.parametrize("masked", [True, False], ids=["mask", "no mask"])
def test_counts_match_velocyto(synthetic, masked):
    """velocyto 0.17.17's own counting code (run10x's steps, with its default logic) gave these counts."""
    reference = pd.read_csv(DATA / f"velocyto_synthetic_{'mask' if masked else 'nomask'}.csv", index_col=0)
    result = utils.fraction_unspliced_from_bam(
        bam_path=synthetic["bam"], barcodes_path=synthetic["barcodes"], gtf_path=synthetic["gtf"],
        repeat_mask_path=synthetic["mask"] if masked else None, tiles=20, cores=1,
    )
    pd.testing.assert_frame_equal(
        result[list(velocyto_rules.LAYERS)].sort_index(), reference.sort_index(), check_names=False, check_dtype=False)


def test_tiles_and_cores_do_not_change_the_counts(synthetic):
    arguments = dict(bam_path=synthetic["bam"], barcodes_path=synthetic["barcodes"], gtf_path=synthetic["gtf"],
                     repeat_mask_path=synthetic["mask"])
    one = utils.fraction_unspliced_from_bam(tiles=1, cores=1, **arguments)
    many = utils.fraction_unspliced_from_bam(tiles=200, cores=2, **arguments)
    pd.testing.assert_frame_equal(one.sort_index(), many.sort_index())


def test_user_regions_count_overlapping_reads(rules_files):
    result = utils.fraction_unspliced_from_bam(
        bam_path=rules_files["bam"], barcodes_path=rules_files["barcodes"], gtf_path=rules_files["gtf"],
        regions=[("chr1", 1100, 1300)], cores=1,
    )
    # The read from 1021 overlaps the region although it starts before it; the one from 1801 does not overlap
    assert molecules(result, "one_molecule_three_reads") == (1, 0, 0)
    assert molecules(result, "intronic") is None


def test_tile_boundaries_lie_between_genes(synthetic):
    models = velocyto_rules.read_gene_models(synthetic["gtf"])
    with pysam.AlignmentFile(synthetic["bam"]) as handle:
        regions = velocyto_rules.tiles_between_genes(handle.references, handle.lengths, models["blocks"], 500)
    assert len(regions) > 20
    half_gap = velocyto_rules.MIN_TILE_GAP // 2
    for contig, start, _ in regions:
        if start == 0:
            continue  # the start of a contig, not a cut
        starts, ends = models["blocks"].get(velocyto_rules.bam_chromosome(contig), ([], []))
        assert not any(a - half_gap < start < b + half_gap for a, b in zip(starts, ends))


def test_parse_cigar_follows_velocyto():
    # A soft clip moves the position as if the bases were aligned
    assert velocyto_rules.parse_cigar([(4, 5), (0, 85)], 1001) == ([(1006, 1090)], False)
    # Deletions and insertions of up to 3 bases join the segments around them
    assert velocyto_rules.parse_cigar([(0, 40), (2, 3), (0, 50)], 1001) == ([(1001, 1093)], False)
    assert velocyto_rules.parse_cigar([(0, 40), (1, 2), (0, 48)], 1001) == ([(1001, 1088)], False)
    assert velocyto_rules.parse_cigar([(0, 40), (2, 4), (0, 50)], 1001) == ([(1001, 1040), (1045, 1094)], False)
    assert velocyto_rules.parse_cigar([(0, 40), (3, 100), (0, 50)], 1001) == ([(1001, 1040), (1141, 1190)], True)


def test_chromosome_names_follow_velocyto():
    assert [velocyto_rules.bam_chromosome(name) for name in ("chr1", "chrM", "chr1_KI270706v1_random", "KI270728.1")] \
        == ["1", "MT", "KI270706v1", "KI270728.1"]
    assert [velocyto_rules.gtf_chromosome(name) for name in ("chr1", "chrM", "MT")] == ["1", "M", "MT"]


def test_nested_repeat_shortens_the_masked_interval(tmp_path):
    path = tmp_path / "rmsk.gtf"
    path.write_text("".join(f"chr1\tr\texon\t{start}\t{end}\t.\t+\t.\tx\n" for start, end in
                            [(1001, 5000), (2001, 2300), (2303, 2400), (8001, 8100), (9001, 9100)])
                    + "chr2\tr\texon\t1\t10\t.\t+\t.\tx\n")
    masks = velocyto_rules.read_repeat_mask(str(path))
    starts, ends = masks["1+"]
    # (1001, 5000) ends at 2400 once the nested repeats are joined to it; the last interval is never searched
    assert list(zip(starts.tolist(), ends.tolist())) == [(1001, 2400)]
    assert "2+" not in masks


def test_validated_intron_changes_only_the_spliced_ambiguous_case():
    # Model 1 crosses from an exon into intron 7, model 2 is exonic only
    record = ((1, 0, True, frozenset({(7, True)})), (2, 0, True, frozenset()))
    assert velocyto_rules.classify(record, validated={7}) == velocyto_rules.SPLICED
    assert velocyto_rules.classify(record, validated=set()) == velocyto_rules.AMBIGUOUS
    # A model of another gene: not counted
    assert velocyto_rules.classify(record + ((3, 1, True, frozenset()),), validated=set()) is None


def test_errors(rules_files, tmp_path):
    with pytest.raises(ValueError, match="needs gtf_path"):
        utils.fraction_unspliced_from_bam(bam_path=rules_files["bam"], barcodes_path=rules_files["barcodes"],
                                          repeat_mask_path=rules_files["mask"])
    with pytest.raises(FileNotFoundError, match="GTF"):
        utils.fraction_unspliced_from_bam(bam_path=rules_files["bam"], barcodes_path=rules_files["barcodes"],
                                          gtf_path=str(tmp_path / "missing.gtf"))
    other_names = tmp_path / "genes.gtf"
    other_names.write_text(
        Path(rules_files["gtf"]).read_text().replace("chr1\t", "NC_000001\t").replace("chrM\t", "NC_012920\t"))
    with pytest.raises(ValueError, match="could be assigned to a gene"):
        utils.fraction_unspliced_from_bam(bam_path=rules_files["bam"], barcodes_path=rules_files["barcodes"],
                                          gtf_path=str(other_names), cores=1)
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("TTTTTTTTTTTTTTTT-1\n")
    with pytest.raises(ValueError, match="matched a barcode"):
        utils.fraction_unspliced_from_bam(bam_path=rules_files["bam"], barcodes_path=str(barcodes),
                                          gtf_path=rules_files["gtf"], cores=1)


def test_command_line(rules_files, tmp_path, capsys):
    from qclus.cli import main

    output = tmp_path / "fraction_unspliced.csv"
    arguments = ["splicing-from-bam", "--bam", rules_files["bam"], "--barcodes", rules_files["barcodes"],
                 "--gtf", rules_files["gtf"], "--repeat-mask", rules_files["mask"], "-o", str(output), "--cores", "1"]
    assert main(arguments) == 0
    written = pd.read_csv(output, index_col=0)
    assert written.loc[BARCODES["all_classes"], "fraction_unspliced"] == pytest.approx(0.25)
    assert list(written.columns) == ["spliced", "unspliced", "ambiguous", "fraction_unspliced"]
    assert capsys.readouterr().out.strip() == f"barcodes\t{len(written)}"

    with pytest.raises(SystemExit) as exit_info:
        main(["splicing-from-bam", "--bam", rules_files["bam"], "--barcodes", rules_files["barcodes"],
              "--repeat-mask", rules_files["mask"], "-o", str(tmp_path / "other.csv")])
    assert exit_info.value.code == 2
    assert "--repeat-mask needs --gtf" in capsys.readouterr().err

    # The table is read by `qclus run` like any other
    fraction = utils.prepare_fraction_unspliced(written)
    assert np.isclose(fraction[BARCODES["all_classes"]], 0.25)


def case_counts(paths, **arguments):
    result = utils.fraction_unspliced_from_bam(
        bam_path=paths["bam"], barcodes_path=paths["barcodes"], gtf_path=paths["gtf"], cores=1, **arguments)
    return tuple(int(x) for x in result.loc[barcode(1), list(velocyto_rules.LAYERS)])


def test_molecule_spliced_across_a_tile_boundary_counts_once(tmp_path):
    # GENE-B starts at 5001; with two tiles the boundary falls at 4501, in the gap before it. The first read starts
    # in the gap and is spliced into GENE-B's first exon; the second, of the same molecule, crosses that exon's end.
    genes = [
        ("chr1", "+", "GENE-A", "TA", [(1001, 3000)]),
        ("chr1", "+", "GENE-B", "TB", [(5001, 5100), (5501, 6000)]),
        ("chr1", "+", "LAST", "TL", [(9001, 9100)]),
    ]
    tags = {"CB": barcode(1) + "-1", "UB": "UMI1", "NH": 1}
    reads = [("chr1", 4001, "50M950N40M", False, tags), ("chr1", 5051, "90M", False, tags)]
    paths = write_case(tmp_path, genes, [], reads, [("chr1", 10_000)], [barcode(1)])
    with pysam.AlignmentFile(paths["bam"]) as handle:
        models = velocyto_rules.read_gene_models(paths["gtf"])
        regions = velocyto_rules.tiles_between_genes(handle.references, handle.lengths, models["blocks"], 2)
    assert regions == [("chr1", 0, 4501), ("chr1", 4501, 10_000)]
    # velocyto: one molecule, unspliced
    assert case_counts(paths, tiles=1) == case_counts(paths, tiles=2) == (0, 1, 0)


def test_transcripts_that_start_together_are_ordered_by_their_gtf_lines(tmp_path):
    # T-A and T-B share their first and last exons; velocyto never searches the last feature of a chromosome
    # strand, which is the last exon of the transcript that sorts last by its GTF line, T-B. The molecule's first
    # read is exonic in T-A and intronic in T-B, its second lies in the last exon: velocyto counts it spliced.
    genes = [
        ("chr1", "+", "GENE", "T-B", [(1001, 1100), (3001, 3100), (9001, 9100)]),
        ("chr1", "+", "GENE", "T-A", [(1001, 1100), (2001, 2100), (9001, 9100)]),
    ]
    tags = {"CB": barcode(1) + "-1", "UB": "UMI1", "NH": 1}
    reads = [("chr1", 2011, "90M", False, tags), ("chr1", 9011, "90M", False, tags)]
    paths = write_case(tmp_path, genes, [], reads, [("chr1", 10_000)], [barcode(1)])
    assert case_counts(paths) == (1, 0, 0)


def test_interval_inverted_by_a_change_of_chromosome_still_touches_reads(tmp_path):
    # The minus strand's last repeat (2000, 2100) moves on to the plus strand, where it is joined to the first
    # repeat and becomes (2000, 200). A read inside (1900, 2200) that crosses position 2000 then touches that
    # interval, and velocyto does not mask it.
    repeats = [
        ("chr1", "-", 2000, 2100),
        ("chr1", "+", 150, 200), ("chr1", "+", 300, 400), ("chr1", "+", 1900, 2200), ("chr1", "+", 3000, 3100),
        ("chr1", "+", 5000, 5100),
        ("chr2", "+", 1, 10),
    ]
    path = tmp_path / "rmsk.gtf"
    path.write_text("".join(f"{c}\tr\texon\t{s}\t{e}\t.\t{strand}\t.\tx\n" for c, strand, s, e in repeats))
    masks = velocyto_rules.read_repeat_mask(str(path))
    starts, ends = masks["1+"]
    assert list(zip(starts.tolist(), ends.tolist())) == [(300, 400), (1900, 2200), (2000, 200)]
    strand = velocyto_rules._Strand(None, masks["1+"], 0, {})
    assert not strand.masked(1980, [(1980, 2069)])
    assert strand.masked(2051, [(2051, 2140)])
