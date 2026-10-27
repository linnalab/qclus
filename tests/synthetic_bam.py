"""
A synthetic Cell Ranger BAM file, with its gene models and repeats, for the tests of qclus.velocyto_rules.

The reads are random but cover each rule of velocyto: molecules of several reads, exonic, intronic, ambiguous and
spliced reads, reads across exon-intron junctions, antisense reads, multimappers, reads without a UMI or a cell
barcode, reads in repeats, soft clips and short indels, genes that overlap on one strand, introns longer than
1 Mb, adjacent exons, and genes on chrM and on a contig without the chr prefix.

tests/data/velocyto_synthetic_*.csv hold velocyto 0.17.17's counts for these files (see make_dataset).
"""
from pathlib import Path

import numpy as np
import pysam

SEED = 7
READ_LENGTH = 90
CONTIGS = [("chr1", 2_600_000), ("chr2", 200_000), ("chrM", 16_569), ("KI270728.1", 50_000)]

# gene: (contig, strand, {transcript: [(start, end), ...]}); 1-based, inclusive, in order of position
GENES = {
    "GA": ("chr1", "+", {
        "GA-1": [(1001, 1200), (3001, 3200), (5001, 5600)],
        "GA-2": [(1001, 1200), (3001, 3400), (5001, 5600)],  # a longer second exon: 3201-3400 is ambiguous
        "GA-3": [(901, 1300), (3001, 3200)],
    }),
    "GB": ("chr1", "-", {  # opposite strand of GA's second intron and last exon
        "GB-1": [(4001, 4200), (6001, 6300)],
        "GB-2": [(4101, 4200), (6001, 6500)],
    }),
    "GC": ("chr1", "+", {"GC-1": [(2001, 2100), (2501, 2600)]}),  # inside GA's first intron, same strand
    "GF": ("chr1", "+", {"GF-1": [(20001, 20100), (20101, 20200), (20301, 20400)]}),  # adjacent exons
    "GD": ("chr1", "+", {  # an intron longer than 1 Mb
        "GD-1": [(10001, 10200), (1_100_001, 1_100_300), (1_200_001, 1_200_200)],
        "GD-2": [(1_100_001, 1_100_300), (1_150_001, 1_150_100)],
    }),
    "GE": ("chr1", "-", {"GE-1": [(1_300_001, 1_300_200), (1_301_001, 1_301_100), (2_400_001, 2_400_200)]}),
    "GH": ("chr2", "-", {
        "GH-1": [(1001, 1500), (2001, 2300), (3001, 3500)],
        "GH-2": [(1001, 1500), (3001, 3500)],
    }),
    "GI": ("chr2", "+", {"GI-1": [(50001, 51000)]}),
    "GM": ("chrM", "+", {"GM-1": [(1, 1000)]}),
    "GK": ("KI270728.1", "+", {"GK-1": [(1001, 2000), (3001, 3500)]}),
}

# (contig, strand, start, end)
REPEATS = [
    ("chr1", "+", 2201, 2400),
    ("chr1", "+", 3601, 3800),
    ("chr1", "+", 3803, 3900),  # within 5 bases of the previous one: joined
    ("chr1", "+", 4001, 4900),
    ("chr1", "+", 4101, 4300),  # nested: velocyto's joined interval then ends at 4300
    ("chr1", "+", 4950, 4990),
    ("chr1", "-", 4401, 4600),
    ("chr1", "-", 6101, 6200),
    ("chr1", "+", 1_100_101, 1_100_200),
    ("chr1", "+", 1_500_001, 1_500_100),
    ("chr2", "-", 1601, 1900),
    ("chr2", "-", 3101, 3300),
    ("chr2", "+", 50101, 50300),
    ("chrM", "+", 101, 200),
]

N_CELLS = 40
N_FOREIGN = 4  # barcodes with reads that are not in the barcode list


def _barcode(rng):
    return "".join(rng.choice(list("ACGT"), 16))


def _umi(rng):
    return "".join(rng.choice(list("ACGT"), 12))


def write_gtf(path):
    lines = []
    for gene, (contig, strand, transcripts) in GENES.items():
        for transcript, exons in transcripts.items():
            numbers = range(1, len(exons) + 1) if strand == "+" else range(len(exons), 0, -1)
            for (start, end), number in zip(exons, numbers):
                tags = f'gene_id "{gene}"; transcript_id "{transcript}"; exon_number "{number}"; gene_name "{gene}";'
                lines.append(f"{contig}\tsynthetic\texon\t{start}\t{end}\t.\t{strand}\t.\t{tags}")
    # Shuffled, so that the readers have to sort
    order = np.random.default_rng(SEED).permutation(len(lines))
    Path(path).write_text("\n".join(lines[i] for i in order) + "\n")


def write_mask(path):
    lines = [
        f"{contig}\thg38_rmsk\texon\t{start}\t{end}\t.\t{strand}\t.\t"
        f'gene_id "R{i}"; transcript_id "R{i}_{i}"; family_id "Alu"; class_id "SINE";'
        for i, (contig, strand, start, end) in enumerate(REPEATS)
    ]
    Path(path).write_text("\n".join(lines) + "\n")


def _place(rng, intervals, length=READ_LENGTH):
    """A 1-based start such that a read of `length` lies inside one of the intervals, if one is long enough."""
    usable = [(start, end) for start, end in intervals if end - start + 1 >= length]
    start, end = usable[rng.integers(len(usable))]
    return int(rng.integers(start, end - length + 2))


def _cigar(rng):
    """A CIGAR of READ_LENGTH query bases: mostly a plain match, sometimes with clips or indels."""
    kind = rng.random()
    if kind < 0.65:
        return [(0, READ_LENGTH)]
    if kind < 0.75:
        clip = int(rng.integers(2, 12))
        return [(4, clip), (0, READ_LENGTH - clip)]
    if kind < 0.82:
        clip = int(rng.integers(2, 12))
        return [(0, READ_LENGTH - clip), (4, clip)]
    if kind < 0.88:
        size = int(rng.integers(1, 4))
        return [(0, 40), (2, size), (0, READ_LENGTH - 40)]
    if kind < 0.92:
        size = int(rng.integers(1, 4))
        return [(0, 40), (1, size), (0, READ_LENGTH - 40 - size)]
    if kind < 0.96:
        size = int(rng.integers(4, 12))
        return [(0, 45), (2, size), (0, READ_LENGTH - 45)]
    gap = int(rng.integers(20, 400))  # a splice that matches no annotated junction
    return [(0, 50), (3, gap), (0, READ_LENGTH - 50)]


def _molecule_reads(rng, gene):
    """Alignments (1-based start, CIGAR, reverse) of one molecule of a gene."""
    contig, strand, transcripts = GENES[gene]
    names = list(transcripts)
    exons = transcripts[names[rng.integers(len(names))]]
    introns = [(a[1] + 1, b[0] - 1) for a, b in zip(exons, exons[1:]) if b[0] - a[1] - 1 >= READ_LENGTH]
    reverse = strand == "-"
    kind = rng.random()
    reads = []
    for _ in range(int(rng.integers(1, 5))):
        part = rng.random()
        if kind < 0.45 or not introns:
            choice = "exon" if part < 0.85 or not introns else "intron"
        elif kind < 0.75:
            choice = "intron" if part < 0.8 else "exon"
        else:
            choice = "junction" if part < 0.5 else "splice"
        if choice == "exon":
            reads.append((_place(rng, exons), _cigar(rng), reverse))
        elif choice == "intron":
            reads.append((_place(rng, introns), _cigar(rng), reverse))
        elif choice == "junction":
            # An unspliced read across the end of an exon into the next intron, or into an exon from an intron
            exon = exons[rng.integers(len(exons))]
            boundary = exon[1] if rng.random() < 0.5 else exon[0]
            start = boundary - int(rng.integers(10, READ_LENGTH - 10))
            reads.append((max(1, start), [(0, READ_LENGTH)], reverse))
        else:
            # Exon pairs separated by an intron; adjacent exons cannot be spliced
            pairs = [i for i in range(len(exons) - 1) if exons[i + 1][0] - exons[i][1] > 20]
            i = pairs[int(rng.integers(len(pairs)))]
            left = int(rng.integers(15, READ_LENGTH - 15))
            jitter = int(rng.integers(-8, 9)) if rng.random() < 0.3 else 0
            gap = exons[i + 1][0] - exons[i][1] - 1 + jitter
            start = exons[i][1] - left + 1
            reads.append((start, [(0, left), (3, gap), (0, READ_LENGTH - left)], reverse))
        if rng.random() < 0.15:
            reverse_read = reads[-1]
            reads[-1] = (reverse_read[0], reverse_read[1], not reverse_read[2])  # antisense
    return contig, reads


def make_reads(rng):
    """All alignments: (contig, 1-based start, CIGAR, reverse, tags, flag)."""
    barcodes = [_barcode(rng) for _ in range(N_CELLS + N_FOREIGN)]
    genes = [gene for gene in GENES if gene != "GM"]
    weights = np.array([6 if gene in ("GA", "GB", "GH") else 2 for gene in genes], dtype=float)
    weights /= weights.sum()
    intergenic_end = min(900_000, dict(CONTIGS)["chr1"] - 1_000)
    alignments = []
    for cell, barcode in enumerate(barcodes):
        if cell == 1:
            continue  # a listed barcode without reads
        n_molecules = 3 if cell == 2 else int(rng.integers(30, 90))
        for _ in range(n_molecules):
            umi = _umi(rng)
            if cell == 2:
                gene = "GM"  # a barcode with mitochondrial reads only
            elif rng.random() < 0.06:
                gene = "GM"
            else:
                gene = genes[rng.choice(len(genes), p=weights)]
            if rng.random() < 0.05:
                # Intergenic
                contig, reads = "chr1", [(int(rng.integers(30_000, intergenic_end)), _cigar(rng), bool(rng.random() < 0.5))]
            else:
                contig, reads = _molecule_reads(rng, gene)
            for start, cigar, reverse in reads:
                copies = 1 + int(rng.random() < 0.3)  # PCR duplicates
                for _ in range(copies):
                    tags = {"CB": f"{barcode}-1", "UB": umi, "NH": 1}
                    roll = rng.random()
                    if roll < 0.03:
                        del tags["UB"]
                    elif roll < 0.05:
                        del tags["CB"]
                        tags["CR"] = barcode
                    elif roll < 0.10:
                        tags["NH"] = 2
                    alignments.append((contig, start, cigar, reverse, tags, 0))
                    if tags["NH"] == 2:
                        # The other alignment of the multimapper
                        other = int(rng.integers(30_000, intergenic_end))
                        alignments.append(("chr1", other, [(0, READ_LENGTH)], reverse, dict(tags), 256))
        # A UMI used by two genes in one cell, of the same region and of two chromosomes: velocyto drops the molecule
        for shared_cell, pair in ((3, ("GA", "GB")), (4, ("GA", "GH"))):
            if cell == shared_cell and all(gene in GENES for gene in pair):
                umi = _umi(rng)
                for gene in pair:
                    contig, reads = _molecule_reads(rng, gene)
                    for start, cigar, reverse in reads:
                        alignments.append((contig, start, cigar, reverse, {"CB": f"{barcode}-1", "UB": umi, "NH": 1}, 0))
    return barcodes, alignments


def write_bam(path, alignments):
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": length} for name, length in CONTIGS],
    })
    contig_ids = {name: i for i, (name, _) in enumerate(CONTIGS)}
    rows = sorted(
        (contig_ids[contig], start, i, cigar, reverse, tags, flag)
        for i, (contig, start, cigar, reverse, tags, flag) in enumerate(alignments)
    )
    with pysam.AlignmentFile(str(path), "wb", header=header) as handle:
        for contig_id, start, i, cigar, reverse, tags, flag in rows:
            read = pysam.AlignedSegment(header)
            read.query_name = f"read{i}"
            length = sum(n for op, n in cigar if op in (0, 1, 4))
            read.query_sequence = "A" * length
            read.query_qualities = pysam.qualitystring_to_array("I" * length)
            read.flag = flag | (16 if reverse else 0)
            read.reference_id = contig_id
            read.reference_start = start - 1
            read.mapping_quality = 255 if tags["NH"] == 1 else 3
            read.cigartuples = cigar
            for tag, value in tags.items():
                read.set_tag(tag, value)
            handle.write(read)
    pysam.index(str(path))


def make_dataset(directory, seed=SEED):
    """Write the BAM file, barcodes, gene models and repeats; return their paths."""
    directory = Path(directory)
    rng = np.random.default_rng(seed)
    barcodes, alignments = make_reads(rng)
    paths = {
        "bam": directory / "possorted_genome_bam.bam",
        "barcodes": directory / "barcodes.tsv",
        "gtf": directory / "genes.gtf",
        "mask": directory / "rmsk.gtf",
    }
    write_bam(paths["bam"], alignments)
    paths["barcodes"].write_text("".join(f"{barcode}-1\n" for barcode in barcodes[:N_CELLS]))
    write_gtf(paths["gtf"])
    write_mask(paths["mask"])
    return {name: str(path) for name, path in paths.items()}
