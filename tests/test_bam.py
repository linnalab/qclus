import pandas as pd
import pysam
import pytest

import qclus.utils as utils

CONTIG_LENGTH = 1000
INSIDE, SPANNING, INTERGENIC, UNLISTED, MIXED = (
    "AAACCCAAGAAACACT-1",
    "AAACCCAAGAAACCAT-1",
    "AAACCCAAGAAACCCA-1",
    "AAACCCAAGAAACCCG-1",
    "AAACCCAAGAAACCTG-1",
)
# (barcode, start, region type). Reads are 50 bases long.
READS = [
    (SPANNING, 10, "E"),
    (INTERGENIC, 20, "I"),
    (UNLISTED, 30, "N"),
    (SPANNING, 90, "N"),  # covers 90 to 140, across the boundary between the first two of ten tiles
    (INTERGENIC, 250, "I"),
    (SPANNING, 300, "N"),
    (MIXED, 500, "E"),
    (MIXED, 520, "E"),
    (MIXED, 540, "N"),
    (INSIDE, 700, "N"),
]
LISTED = [INSIDE, SPANNING, INTERGENIC, MIXED]


def name(barcode):
    return barcode.split("-")[0]


@pytest.fixture(scope="module")
def bam(tmp_path_factory):
    directory = tmp_path_factory.mktemp("bam")
    path = directory / "possorted_genome_bam.bam"
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": CONTIG_LENGTH}]})
    with pysam.AlignmentFile(path, "wb", header=header) as handle:
        for i, (barcode, start, region) in enumerate(READS):
            read = pysam.AlignedSegment(header)
            read.query_name = f"read{i}"
            read.query_sequence = "A" * 50
            read.query_qualities = pysam.qualitystring_to_array("I" * 50)
            read.flag = 0
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = 60
            read.cigartuples = [(0, 50)]
            read.set_tag("CB", barcode)
            read.set_tag("RE", region)
            handle.write(read)
    pysam.index(str(path))

    barcodes = directory / "barcodes.tsv"
    barcodes.write_text("\n".join(LISTED) + "\n")
    return str(path), str(barcodes)


def test_fractions_from_generated_tiles(bam):
    bam_path, barcodes_path = bam
    result = utils.fraction_unspliced_from_bam(bam_path=bam_path, barcodes_path=barcodes_path, tiles=10, cores=1)

    expected = pd.Series(
        {
            name(INSIDE): 1.0,
            # One exonic and two intronic reads; the read across the tile boundary is counted once
            name(SPANNING): 2 / 3,
            # Only intergenic reads
            name(INTERGENIC): 0.0,
            name(MIXED): 1 / 3,
        },
        name="fraction_unspliced",
    )
    pd.testing.assert_series_equal(result["fraction_unspliced"].sort_index(), expected.sort_index())
    assert list(result.columns) == ["fraction_unspliced"]


def test_reads_overlapping_a_supplied_region_are_counted(bam):
    bam_path, barcodes_path = bam
    result = utils.fraction_unspliced_from_bam(
        bam_path=bam_path, barcodes_path=barcodes_path, regions=[("chr1", 100, 200)], cores=1
    )
    # The read starting at 90 overlaps the region although it starts before it
    assert result["fraction_unspliced"].to_dict() == {name(SPANNING): 1.0}


def test_index_under_another_name_reaches_the_workers(bam, tmp_path):
    bam_path, barcodes_path = bam
    copy = tmp_path / "copy.bam"
    copy.write_bytes(open(bam_path, "rb").read())
    index = tmp_path / "index_elsewhere.bai"
    pysam.index(str(copy), str(index))

    result = utils.fraction_unspliced_from_bam(
        bam_path=str(copy), bam_index_path=str(index), barcodes_path=barcodes_path, tiles=10, cores=2
    )
    assert result.loc[name(SPANNING), "fraction_unspliced"] == pytest.approx(2 / 3)


def test_no_matching_barcodes(bam, tmp_path):
    bam_path, _ = bam
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("\n".join(barcode.split("-")[0] for barcode in LISTED) + "\n")
    with pytest.raises(ValueError, match="carry the same suffix"):
        utils.fraction_unspliced_from_bam(bam_path=bam_path, barcodes_path=str(barcodes), tiles=10, cores=1)


def test_more_tiles_than_reference_bases(bam):
    bam_path, barcodes_path = bam
    result = utils.fraction_unspliced_from_bam(
        bam_path=bam_path, barcodes_path=barcodes_path, tiles=CONTIG_LENGTH * 5, cores=1
    )
    assert result.loc[name(SPANNING), "fraction_unspliced"] == pytest.approx(2 / 3)


def test_invalid_arguments(bam, tmp_path):
    bam_path, barcodes_path = bam
    with pytest.raises(ValueError, match="tiles must be at least 1"):
        utils.fraction_unspliced_from_bam(bam_path=bam_path, barcodes_path=barcodes_path, tiles=0)
    with pytest.raises(ValueError, match="Please provide"):
        utils.fraction_unspliced_from_bam(bam_path=bam_path)
    with pytest.raises(FileNotFoundError):
        utils.fraction_unspliced_from_bam(bam_path=str(tmp_path / "missing.bam"), barcodes_path=barcodes_path)
    with pytest.raises(FileNotFoundError):
        utils.fraction_unspliced_from_bam(
            bam_path=bam_path, bam_index_path=str(tmp_path / "missing.bai"), barcodes_path=barcodes_path
        )


def test_command_line(bam, tmp_path, capsys):
    from qclus.cli import main

    bam_path, barcodes_path = bam
    output = tmp_path / "fraction_unspliced.csv"
    arguments = ["splicing-from-bam", "--bam", bam_path, "--barcodes", barcodes_path, "-o", str(output), "--tiles", "10", "--cores", "1"]

    assert main(arguments) == 0
    written = pd.read_csv(output, index_col=0)
    assert written.loc[name(SPANNING), "fraction_unspliced"] == pytest.approx(2 / 3)
    assert capsys.readouterr().out.strip() == f"barcodes\t{len(LISTED)}"

    with pytest.raises(SystemExit) as exit_info:
        main(arguments)
    assert exit_info.value.code == 2


@pytest.mark.parametrize("suffix", [".bam.bai", ".bai", ".bam.csi", ".csi"])
def test_command_line_protects_the_index_it_finds(bam, tmp_path, capsys, suffix):
    from qclus.cli import main

    bam_path, barcodes_path = bam
    copy = tmp_path / "sample.bam"
    copy.write_bytes(open(bam_path, "rb").read())
    index = tmp_path / ("sample" + suffix)
    if suffix.endswith(".csi"):
        pysam.index("-c", str(copy), str(index))
    else:
        pysam.index(str(copy), str(index))
    before = index.read_bytes()

    # Each of these names is one that pysam finds by itself, so the reader does use the file
    with pysam.AlignmentFile(str(copy), "rb") as handle:
        assert handle.has_index()

    with pytest.raises(SystemExit) as exit_info:
        main(["splicing-from-bam", "--bam", str(copy), "--barcodes", barcodes_path, "-o", str(index), "--overwrite", "--cores", "1"])

    assert exit_info.value.code == 2
    assert "inputs are never overwritten" in capsys.readouterr().err
    assert index.read_bytes() == before
