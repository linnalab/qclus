"""
Count spliced, unspliced and ambiguous molecules per barcode from a BAM file with the rules of velocyto.

velocyto 0.17.17, run as `velocyto run10x` with its default logic (Permissive10X), is QClus's reference for the
fraction of unspliced reads. This module counts the same molecules from the position-sorted BAM file of Cell
Ranger, without the two steps that make velocyto slow and large: it does not sort the BAM file by barcode, and it
keeps three totals per barcode instead of a gene by barcode matrix per layer.

The rules, with the velocyto file each comes from:

- Reads (counter.py, iter_alignments): only unique alignments (NH = 1) that carry a cell barcode and a UMI, of the
  given barcodes. Barcodes are compared without their suffix ('-1'). Reads that span more than 3 Mb are dropped.
- Gene models (counter.py, read_transcriptmodels; transcript_model.py): the exons of each transcript of the GTF
  file, with an intron between consecutive exons. Transcript models are cut at introns longer than 1 Mb.
- Matching (indexes.py, find_overlapping_ivls): each read is matched to the exons and introns it overlaps, on its
  own strand only, by more than 5 bases. Of the transcript models it overlaps, a read keeps those with the most
  matching segments, and a spliced read keeps only those whose exon or intron ends lie within 6 bases of the
  read's segment ends.
- Molecules (molitem.py): the reads of one barcode and UMI form a molecule, which keeps the transcript models that
  all its reads are compatible with.
- Repeats (counter.py, read_repeats and _count_cell_batch_stranded): a read that lies entirely inside a masked
  repeat on its own strand is left out. Repeats closer than 5 bases are joined first.
- Introns are validated (counter.py, mark_up_introns) when an unspliced read crosses one of their ends into the
  neighbouring exon.
- Classes (logic.py, Permissive10X): a molecule compatible with exactly one gene is spliced if all its
  transcript models are exonic only, ambiguous if it is exonic in some models and intronic in others, and
  unspliced otherwise. Molecules compatible with no gene or with several are not counted.

Some details of velocyto's code change the counts, and are reproduced so that the result matches its loom files:

- velocyto names chromosomes differently in the BAM file and in the GTF file: 'chrM' becomes 'MT' in one and 'M' in
  the other, so reads on the mitochondrial genome are never counted when the reference names it 'chrM', as Cell
  Ranger's references do.
- A soft clip moves the position of the following segments, as if the clipped bases were aligned.
- When repeats are joined, the end of the joined interval becomes the end of the last repeat added, even if that
  repeat ends earlier, so a repeat nested inside a longer one shortens the masked interval. The first repeat of
  each chromosome strand is joined to the last interval of the one before it and lost, the last repeat of each
  chromosome strand is never stored, and neither are the repeats of the last chromosome strand of the file.
- The last feature (in order of position) of each chromosome strand, and the last masked interval, are never
  searched.

The genome is processed in tiles whose boundaries lie between genes, so that all reads of a molecule, and all reads
that can validate an intron, reach the same worker. What remains different from velocyto:

- velocyto keeps the reads of a barcode and UMI together across the whole genome, so a UMI that a barcode uses in two
  genes of different tiles makes one molecule, which it drops; here it makes two. With 12-base UMIs this is rare.
- A read without an NH tag is taken as unique; velocyto stops with an error.
- If any exon lacks an exon_number, all exons are numbered by position; velocyto looks at the first 500 lines only.
"""
import gzip
import logging
import re
from bisect import bisect_left, bisect_right
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from contextlib import ExitStack
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from tqdm import tqdm

logger = logging.getLogger("qclus")

# Constants of velocyto 0.17.17 (constants.py and counter.py)
MIN_FLANK = 5
PATCH_INDELS = 3
SPLIC_INACUR = 6
LONGEST_INTRON_ALLOWED = 1_000_000
LONGEST_READ_SPAN = 3_000_000
MASK_TOLERANCE = 5
MATCH_INSIDE = 1

# Tile boundaries are placed in gaps of at least this size between genes
MIN_TILE_GAP = 1_000
# Active features are pruned when a read is this far past the last pruning
PRUNE_STEP = 1_000

LAYERS = ("spliced", "unspliced", "ambiguous")
SPLICED, UNSPLICED, AMBIGUOUS = range(3)
# The mapping record of a read that matches no transcript model
EMPTY = frozenset()

_TRANSCRIPT_ID = re.compile('transcript_id "([^"]+)"')
_GENE_ID = re.compile('gene_id "([^"]+)"')
_EXON_NUMBER = re.compile(r'exon_number "*?([\w]+)')


def bam_chromosome(name: str) -> str:
    """Chromosome name of a read, as velocyto writes it (counter.py, iter_alignments)."""
    if name.startswith("chr"):
        if "_" in name:
            return name.split("_")[1]
        name = name[3:]
        return "MT" if name == "M" else name
    return name


def gtf_chromosome(name: str) -> str:
    """Chromosome name of a gene model, as velocyto writes it (counter.py, read_transcriptmodels)."""
    return name[3:] if "chr" in name[:4] else name


def mask_chromosome(name: str) -> str:
    """Chromosome name of a masked repeat, as velocyto writes it (counter.py, read_repeats)."""
    return name[3:] if name[:3].lower() == "chr" else name


def parse_cigar(cigartuples: Sequence[Tuple[int, int]], pos: int) -> Tuple[List[Tuple[int, int]], bool]:
    """
    Split an alignment into the reference segments velocyto uses (counter.py, parse_cigar_tuple).

    Parameters:
        cigartuples (Sequence[Tuple[int, int]]): CIGAR operations and lengths, as pysam gives them.
        pos (int): 1-based position of the alignment.

    Returns:
        Tuple[List[Tuple[int, int]], bool]: 1-based inclusive segments, and whether the read is spliced (has an N).
    """
    segments = []
    hole_to_remove = set()
    ref_skip = False
    p = pos
    for i, (operation, length) in enumerate(cigartuples):
        if operation == 0:  # M
            segments.append((p, p + length - 1))
            p += length
        elif operation == 3:  # N, a splice
            ref_skip = True
            p += length
        elif operation == 2:  # D: a short deletion between two matches joins them
            if length <= PATCH_INDELS:
                try:
                    if cigartuples[i + 1][0] == 0 and cigartuples[i - 1][0] == 0:
                        hole_to_remove.add(len(segments) - 1)
                except IndexError:
                    pass
            p += length
        elif operation == 4:  # S: velocyto moves the position over soft clipped bases too
            p += length
        elif operation == 1:  # I: a short insertion between two matches joins them
            if length <= PATCH_INDELS:
                try:
                    if cigartuples[i + 1][0] == 0 and cigartuples[i - 1][0] == 0:
                        hole_to_remove.add(len(segments) - 1)
                except IndexError:
                    pass
    for a, b in enumerate(sorted(hole_to_remove)):
        segments[b - a] = (segments.pop(b - a)[0], segments[b - a][1])
    return segments, ref_skip


def _open_text(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _chop_long_intron(features: List[tuple], plus: bool) -> List[tuple]:
    """
    Cut a transcript model at an intron longer than 1 Mb (transcript_model.py, chop_if_long_intron).

    The part upstream of the last such intron (plus strand) or downstream of the first (minus strand) is removed, so
    that a very long intron does not hide the genes inside it, and exons and introns are numbered again.
    """
    long_introns = [f for f in features if f[2] and f[1] - f[0] + 1 > LONGEST_INTRON_ALLOWED]
    if not long_introns:
        return features
    if plus:
        cut = long_introns[-1]
        kept = [f for f in features if (f[0], f[1]) > (cut[0], cut[1])]
    else:
        cut = long_introns[0]
        kept = [f for f in features if (f[0], f[1]) < (cut[0], cut[1])][::-1]
    numbered, exon_count, intron_count = [], 1, 1
    for start, end, is_intron, _ in kept:
        if is_intron:
            numbered.append((start, end, True, intron_count))
            intron_count += 1
        else:
            numbered.append((start, end, False, exon_count))
            exon_count += 1
    return numbered if plus else numbered[::-1]


def read_gene_models(gtf_path: str) -> dict:
    """
    Read the transcript models of a GTF file as velocyto does (counter.py, read_transcriptmodels).

    Parameters:
        gtf_path (str): GTF file, optionally gzipped, such as genes/genes.gtf of the Cell Ranger reference.

    Returns:
        dict: 'features' maps each chromosome strand (velocyto's name plus '+' or '-') to arrays of its exons and
            introns, ordered by position: 'start', 'end' (1-based, inclusive), 'intron', 'model', 'gene', 'down'
            (start of the next exon of an intron, or -1), 'up' (end of the previous exon, or -1) and 'reach'
            (the larger of start and end). 'blocks' maps each chromosome to the sorted (start, end) arrays of the
            regions covered by gene models on either strand, joined when closer than MIN_TILE_GAP.

    Raises:
        ValueError: If the file holds no exon with a transcript_id and a gene_id.
    """
    exons = []
    gene_index: Dict[str, int] = {}
    with _open_text(gtf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            # Exons without a strand cannot match a read, which is always on one
            if len(fields) < 9 or fields[2] != "exon" or fields[6] not in ("+", "-"):
                continue
            tags = fields[8]
            transcript, gene = _TRANSCRIPT_ID.search(tags), _GENE_ID.search(tags)
            if transcript is None or gene is None:
                continue
            number = _EXON_NUMBER.search(tags)
            exons.append((
                fields[0], fields[6], int(fields[3]), int(fields[4]), transcript.group(1),
                gene_index.setdefault(gene.group(1), len(gene_index)),
                int(number.group(1)) if number is not None else None,
            ))
    if not exons:
        raise ValueError(f"'{gtf_path}' holds no exon with a transcript_id and a gene_id.")

    # velocyto sorts the lines as `sort -k1,1 -k7,7 -k4,4n` does, minus strand first
    exons.sort(key=lambda exon: (exon[0], exon[1] == "+", exon[2], exon[3]))
    if any(exon[6] is None for exon in exons):
        # Without exon numbers velocyto numbers the exons of each transcript by position (peek_and_correct)
        position: Dict[tuple, int] = {}
        counts: Dict[tuple, int] = {}
        for exon in exons:
            counts[(exon[0], exon[1], exon[4])] = counts.get((exon[0], exon[1], exon[4]), 0) + 1
        numbered = []
        for exon in exons:
            key = (exon[0], exon[1], exon[4])
            position[key] = position.get(key, 0) + 1
            number = position[key] if exon[1] == "+" else counts[key] - position[key] + 1
            numbered.append(exon[:6] + (number,))
        exons = numbered

    # Transcript models in order of their first exon; each feature is (start, end, is_intron, number)
    models: Dict[Tuple[str, str], Tuple[int, List[tuple]]] = {}
    for chrom, strand, start, end, transcript, gene, number in exons:
        key = (gtf_chromosome(chrom) + strand, transcript)
        model = models.get(key)
        if model is None:
            model = models[key] = (gene, [])
        features = model[1]
        if features:
            previous = features[-1]
            intron_number = previous[3] if strand == "+" else previous[3] - 1
            features.append((previous[1] + 1, start - 1, True, intron_number))
        features.append((start, end, False, number))
    del exons

    by_chromstrand: Dict[str, list] = {}
    for (chromstrand, _), (gene, features) in models.items():
        features = _chop_long_intron(features, chromstrand.endswith("+"))
        by_chromstrand.setdefault(chromstrand, []).append((gene, features))
    del models

    result: Dict[str, dict] = {}
    extents: Dict[str, list] = {}
    model_id = 0
    for chromstrand, chromstrand_models in by_chromstrand.items():
        plus = chromstrand.endswith("+")
        # velocyto orders the models by their first position and then all features by (start, end); sorts are stable
        chromstrand_models.sort(key=lambda model: model[1][0][0])
        rows = []
        for gene, features in chromstrand_models:
            extents.setdefault(chromstrand[:-1], []).append((features[0][0], features[-1][1]))
            n = len(features)
            for start, end, is_intron, number in features:
                down = up = -1
                if is_intron:
                    # The neighbouring exons, found by velocyto from the intron's number (feature.py)
                    down_index = number * 2 if plus else n - 2 * number + 1
                    up_index = number * 2 - 2 if plus else n - 2 * number - 1
                    if -n <= down_index < n:
                        down = features[down_index][0]
                    if -n <= up_index < n:
                        up = features[up_index][1]
                rows.append((start, end, is_intron, model_id, gene, down, up))
            model_id += 1
        rows.sort(key=lambda row: (row[0], row[1]))
        # velocyto's search loops stop before the last feature (indexes.py, `i < self.maxiidx`)
        rows = rows[:-1]
        if not rows:
            continue
        columns = list(zip(*rows))
        arrays = {
            "start": np.array(columns[0], dtype=np.int64),
            "end": np.array(columns[1], dtype=np.int64),
            "intron": np.array(columns[2], dtype=bool),
            "model": np.array(columns[3], dtype=np.int64),
            "gene": np.array(columns[4], dtype=np.int64),
            "down": np.array(columns[5], dtype=np.int64),
            "up": np.array(columns[6], dtype=np.int64),
        }
        arrays["reach"] = np.maximum(arrays["start"], arrays["end"])
        result[chromstrand] = arrays

    blocks = {}
    for chromosome, spans in extents.items():
        spans.sort()
        merged = [list(spans[0])]
        for start, end in spans[1:]:
            if start <= merged[-1][1] + MIN_TILE_GAP:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        merged = np.array(merged, dtype=np.int64)
        blocks[chromosome] = (merged[:, 0], merged[:, 1])
    return {"features": result, "blocks": blocks}


def read_repeat_mask(mask_path: str) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """
    Read the repeats to mask from a GTF file and join them as velocyto does (counter.py, read_repeats).

    Parameters:
        mask_path (str): GTF file of repeats, optionally gzipped, such as RepeatMasker's track from UCSC.

    Returns:
        Dict[str, Tuple[np.ndarray, np.ndarray]]: Start and end (1-based, inclusive) of the masked intervals of each
            chromosome strand, ordered by position.

    Raises:
        ValueError: If the file holds no repeat.
    """
    chroms, strands, starts, ends, kinds = [], [], [], [], []
    chrom_codes: Dict[str, int] = {}
    kind_codes: Dict[str, int] = {}
    with _open_text(mask_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 9:
                continue
            chroms.append(chrom_codes.setdefault(fields[0], len(chrom_codes)))
            strands.append(fields[6])
            starts.append(int(fields[3]))
            ends.append(int(fields[4]))
            kinds.append(kind_codes.setdefault(fields[1] + "\t" + fields[2], len(kind_codes)))
    if not starts:
        raise ValueError(f"'{mask_path}' holds no repeat.")

    # velocyto sorts by chromosome, strand (minus first), start and then the whole line, in which the source,
    # the feature type and the end come next
    def ranks(codes: Dict[str, int]) -> np.ndarray:
        rank = np.empty(len(codes), dtype=np.int64)
        for position, name in enumerate(sorted(codes)):
            rank[codes[name]] = position
        return rank

    chroms = ranks(chrom_codes)[np.array(chroms, dtype=np.int64)]
    kinds = ranks(kind_codes)[np.array(kinds, dtype=np.int64)]
    starts = np.array(starts, dtype=np.int64)
    ends = np.array(ends, dtype=np.int64)
    plus = np.array([strand == "+" for strand in strands])
    order = np.lexsort((np.char.mod("%d", ends).astype("S"), kinds, starts, plus, chroms))

    names = [mask_chromosome(name) for name in sorted(chrom_codes)]
    strands = np.array(strands, dtype=object)[order].tolist()
    chroms, starts, ends = chroms[order].tolist(), starts[order].tolist(), ends[order].tolist()

    # The loop of velocyto, kept as it is: on a change of chromosome or strand, the interval being joined moves on
    # to the next list, and the list of the last chromosome and strand is never stored
    masks: Dict[str, list] = {}
    intervals: list = []
    current = names[chroms[0]] + strands[0]
    current_start, current_end = starts[0], ends[0]
    for chrom, strand, start, end in zip(chroms[1:], strands[1:], starts[1:], ends[1:]):
        chromstrand = names[chrom] + strand
        if chromstrand != current:
            masks[current] = intervals
            intervals = []
            current = chromstrand
        if start > current_end + MASK_TOLERANCE:
            intervals.append((current_start, current_end))
            current_start, current_end = start, end
        else:
            current_end = end

    result = {}
    for chromstrand, intervals in masks.items():
        # velocyto never searches the last interval; one whose start lies after its end can contain no read
        intervals = sorted(intervals)[:-1]
        intervals = [interval for interval in intervals if interval[0] <= interval[1]]
        if intervals:
            array = np.array(intervals, dtype=np.int64)
            result[chromstrand] = (array[:, 0].copy(), array[:, 1].copy())
    return result


def tiles_between_genes(
    references: Sequence[str],
    lengths: Sequence[int],
    blocks: Dict[str, Tuple[np.ndarray, np.ndarray]],
    tiles: int,
) -> List[Tuple[str, int, int]]:
    """
    Split the genome into about `tiles` regions, with boundaries in gaps between genes.

    Parameters:
        references (Sequence[str]): Reference sequence names of the BAM header.
        lengths (Sequence[int]): Their lengths.
        blocks (Dict[str, Tuple[np.ndarray, np.ndarray]]): Regions covered by genes, from read_gene_models.
        tiles (int): Approximate number of regions.

    Returns:
        List[Tuple[str, int, int]]: (contig, start, end) in 0-based, half-open coordinates.
    """
    tile_size = max(1, sum(lengths) // tiles)
    half_gap = MIN_TILE_GAP // 2
    regions = []
    for contig, length in zip(references, lengths):
        starts, ends = blocks.get(bam_chromosome(contig), (np.zeros(0, np.int64), np.zeros(0, np.int64)))
        bounds = [0]
        for target in range(tile_size, length, tile_size):
            i = bisect_right(starts, target) - 1
            if i >= 0 and target <= ends[i]:
                # Inside a gene: cut in the middle of the gap that follows
                if i + 1 >= len(starts):
                    continue
                cut = int(ends[i] + starts[i + 1]) // 2
            else:
                low = int(ends[i]) + half_gap if i >= 0 else 0
                high = int(starts[i + 1]) - half_gap if i + 1 < len(starts) else length
                cut = min(max(target, low), high)
            if bounds[-1] < cut < length:
                bounds.append(cut)
        bounds.append(length)
        regions.extend((contig, start, end) for start, end in zip(bounds, bounds[1:]) if start < end)
    return regions


def tile_inputs(
    models: dict,
    masks: Optional[Dict[str, Tuple[np.ndarray, np.ndarray]]],
    contig: str,
    start: int,
    end: int,
    count_by_start: bool,
) -> Tuple[dict, dict]:
    """
    The features and masked intervals that reads of one region can reach.

    Returns:
        Tuple[dict, dict]: For '+' and '-', the feature arrays and the (start, end) arrays of masked intervals.
    """
    upper = end + LONGEST_READ_SPAN
    lower = start if count_by_start else start - LONGEST_READ_SPAN
    chromosome = bam_chromosome(contig)
    features, intervals = {}, {}
    for strand in "+-":
        arrays = models["features"].get(chromosome + strand)
        if arrays is not None:
            keep = (arrays["start"] < upper) & (arrays["reach"] > lower)
            if keep.any():
                features[strand] = {name: values[keep] for name, values in arrays.items()}
        if masks is not None and chromosome + strand in masks:
            mask_starts, mask_ends = masks[chromosome + strand]
            keep = (mask_starts < upper) & (mask_ends > lower)
            if keep.any():
                intervals[strand] = (mask_starts[keep], mask_ends[keep])
    return features, intervals


class _Strand:
    """The exons, introns and masked intervals of one strand of a tile, searched by reads in order of position."""

    def __init__(self, features: Optional[dict], masks: Optional[tuple], offset: int, intern: dict):
        empty = features is None
        self.start = [] if empty else features["start"].tolist()
        self.end = [] if empty else features["end"].tolist()
        self.intron = [] if empty else features["intron"].tolist()
        self.model = [] if empty else features["model"].tolist()
        self.gene = [] if empty else features["gene"].tolist()
        self.down = [] if empty else features["down"].tolist()
        self.up = [] if empty else features["up"].tolist()
        self.reach = [] if empty else features["reach"].tolist()
        self.offset = offset
        self.next = 0
        self.active: List[int] = []
        self.pruned_at = 0
        # Mapping records by the matches that make them, and one copy of each record shared by both strands
        self.records: Dict[tuple, frozenset] = {}
        self.intern = intern

        self.mask_start = [] if masks is None else masks[0].tolist()
        self.mask_end = [] if masks is None else masks[1].tolist()
        # Largest start or end of the intervals so far, to stop searching backwards
        self.mask_reach = np.maximum.accumulate(np.maximum(masks[0], masks[1])).tolist() if masks is not None else []

    def candidates(self, pos: int, last: int) -> List[int]:
        """Features that a read from `pos` to `last` may overlap; also features that no later read can reach."""
        start, n = self.start, len(self.start)
        while self.next < n and start[self.next] < last:
            self.active.append(self.next)
            self.next += 1
        if pos - self.pruned_at > PRUNE_STEP:
            limit, reach = pos + MIN_FLANK, self.reach
            self.active = [f for f in self.active if reach[f] > limit]
            self.pruned_at = pos
        return self.active

    def masked(self, segments: List[Tuple[int, int]]) -> bool:
        """Whether every segment lies inside a masked interval without touching another (indexes.py, has_ivls_enclosing)."""
        mask_start, mask_end, mask_reach = self.mask_start, self.mask_end, self.mask_reach
        if not mask_start:
            return False
        for s, e in segments:
            found = 0
            # Intervals that start before the segment ends; those reaching no further than `stop` cannot matter
            k = bisect_left(mask_start, e) - 1
            stop = min(e - MIN_FLANK, s + MIN_FLANK + 1)
            while k >= 0 and mask_reach[k] >= stop:
                ms, me = mask_start[k], mask_end[k]
                if s + MIN_FLANK >= ms and e - MIN_FLANK <= me and e - s > MIN_FLANK:
                    found |= 1
                if s + MIN_FLANK < ms < e - MIN_FLANK:
                    found |= 2
                if s + MIN_FLANK < me < e - MIN_FLANK:
                    found |= 4
                k -= 1
            if found != MATCH_INSIDE:
                return False
        return True

    def match(self, pos: int, segments: List[Tuple[int, int]], spliced: bool, validated: set) -> frozenset:
        """
        Match a read to the transcript models of this strand.

        Introns whose end an unspliced read crosses into the neighbouring exon are added to `validated`
        (indexes.py, mark_overlapping_ivls). The read's mapping record (indexes.py, find_overlapping_ivls) is
        returned summarised per transcript model, as a frozenset of (model, gene, has an exon, frozenset of
        (intron, crosses the intron's end)).
        """
        start, end, intron, down, up = self.start, self.end, self.intron, self.down, self.up
        offset = self.offset
        # Each match of a segment and a feature as 4 * feature + 2 * crosses an intron end + agrees with a splice
        matches = []
        for f in self.candidates(pos, segments[-1][1]):
            fs, fe = start[f], end[f]
            for s, e in segments:
                if fs >= e:
                    continue
                inner_s, inner_e = s + MIN_FLANK, e - MIN_FLANK
                crosses = False
                if intron[f]:
                    crosses = (inner_s < fe < inner_e and inner_s < down[f] < inner_e) or (
                        inner_s < fs < inner_e and inner_s < up[f] < inner_e)
                    if crosses and not spliced:
                        validated.add(f + offset)
                if e - s > MIN_FLANK and inner_e > fs and inner_s < fe:
                    agrees = not spliced or abs(fs - s) <= SPLIC_INACUR or abs(fe - e) <= SPLIC_INACUR
                    matches.append(4 * f + 2 * crosses + agrees)
        if not matches:
            return EMPTY
        # Reads in the same exons and introns have the same record
        signature = tuple(matches)
        record = self.records.get(signature)
        if record is None:
            record = self.records[signature] = self._record(signature)
        return record

    def _record(self, signature: tuple) -> frozenset:
        model, gene, intron, offset = self.model, self.gene, self.intron, self.offset
        entries: Dict[int, list] = {}
        for code in signature:
            f = code >> 2
            entry = entries.get(model[f])
            if entry is None:
                # matches, gene, has an exon, introns, ends agree with a splice
                entry = entries[model[f]] = [0, gene[f], False, set(), True]
            entry[0] += 1
            if intron[f]:
                entry[3].add((f + offset, bool(code & 2)))
            else:
                entry[2] = True
            if not code & 1:
                entry[4] = False
        # Only the models with the most matching segments, and for a spliced read only those whose exon or intron
        # ends agree with the splice, are kept
        most = max(entry[0] for entry in entries.values())
        record = frozenset(
            (tm, entry[1], entry[2], frozenset(entry[3]))
            for tm, entry in entries.items() if entry[0] == most and entry[4]
        )
        return self.intern.setdefault(record, record)



def _merge(first: frozenset, second: frozenset) -> frozenset:
    """The transcript models two mapping records share, with their matches combined (molitem.py)."""
    other = {entry[0]: entry for entry in second}
    return frozenset(
        (tm, gene, has_exon or other[tm][2], introns | other[tm][3])
        for tm, gene, has_exon, introns in first if tm in other
    )


def classify(record: frozenset, validated: set) -> Optional[int]:
    """
    The class of a molecule by velocyto's Permissive10X logic (logic.py).

    Parameters:
        record (frozenset): The molecule's mapping record, as returned by _Strand.match and combined by _merge.
        validated (set): Introns crossed into their neighbouring exon by an unspliced read.

    Returns:
        Optional[int]: SPLICED, UNSPLICED or AMBIGUOUS, or None if the molecule is not counted.
    """
    if not record:
        return None
    gene = next(iter(record))[1]
    if any(entry[1] != gene for entry in record):
        return None

    onlyintron = onlyintron_valid = onlyexon = mixed = valid_mixed = invalid_mixed = False
    only_span = True
    for _, _, has_exons, introns in record:
        has_introns = bool(introns)
        has_valid = any(f in validated for f, _ in introns)
        has_span = any(crosses and f in validated for f, crosses in introns)
        if has_valid and not has_exons:
            onlyintron_valid = True
        if has_introns and not has_exons:
            onlyintron = True
        if has_exons and not has_introns:
            onlyexon = True
        if has_exons and has_introns and not has_valid and not has_span:
            invalid_mixed = mixed = True
        if has_exons and has_introns and has_valid and not has_span:
            valid_mixed = mixed = True
        if not has_span:
            only_span = False

    # The branches of Permissive10X.count, in its order
    if onlyexon and not onlyintron and not mixed:
        return SPLICED
    if only_span:
        return UNSPLICED
    if onlyintron_valid and not mixed and not onlyexon:
        return UNSPLICED
    if onlyintron and not onlyintron_valid and not mixed and not onlyexon:
        return UNSPLICED
    if invalid_mixed and not valid_mixed and not onlyintron and not onlyexon and not only_span:
        return UNSPLICED
    if valid_mixed and not onlyintron and not onlyexon and not only_span:
        return UNSPLICED
    if onlyintron and onlyexon and not mixed:
        return AMBIGUOUS
    if onlyintron and not onlyexon and mixed:
        return UNSPLICED
    if not onlyintron and onlyexon and mixed:
        return AMBIGUOUS
    if onlyintron and onlyexon and mixed:
        return AMBIGUOUS
    return None


def count_tile(
    interval: tuple,
    bam_path: str,
    barcodes: Dict[str, int],
    features: dict,
    masks: dict,
    CB_tag: str = "CB",
    UB_tag: str = "UB",
    bam_index_path: Optional[str] = None,
    count_by_start: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Count the spliced, unspliced and ambiguous molecules of each barcode in one genomic region.

    Parameters:
        interval (tuple): (contig, start, end), 0-based and half-open.
        bam_path (str): Position-sorted BAM file.
        barcodes (Dict[str, int]): Barcodes without suffix, mapped to their row in the result.
        features (dict): Feature arrays per strand, from tile_inputs.
        masks (dict): Masked intervals per strand, from tile_inputs.
        CB_tag (str, optional): Tag of the cell barcode.
        UB_tag (str, optional): Tag of the UMI.
        bam_index_path (str, optional): BAM index. Found automatically if omitted.
        count_by_start (bool, optional): If True, use only reads that start inside the region.

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray]: Rows of the barcodes with a read that velocyto reads; rows of
            the barcodes with counted molecules; and their counts (one column per class of LAYERS).
    """
    import pysam

    contig, start, end = interval
    plus_features = features.get("+")
    offset = 0 if plus_features is None else len(plus_features["start"])
    records: Dict[frozenset, frozenset] = {}
    strands = (
        _Strand(plus_features, masks.get("+"), 0, records),
        _Strand(features.get("-"), masks.get("-"), offset, records),
    )
    has_features = bool(features)

    seen = set()
    molecules: Dict[tuple, frozenset] = {}
    validated: set = set()
    merged_records: Dict[tuple, frozenset] = {}
    # Rows of the cell barcode tags seen so far, None for barcodes not in the list
    tag_rows: Dict[str, Optional[int]] = {}
    # Reads that start at the same position and align the same way get the same result
    cache: Dict[tuple, Optional[frozenset]] = {}
    cache_position = -1

    with pysam.AlignmentFile(bam_path, "rb", index_filename=bam_index_path) as bam_file:
        for read in bam_file.fetch(contig, start, end):
            reference_start = read.reference_start
            if count_by_start and reference_start < start:
                continue
            if read.is_unmapped:
                continue
            try:
                if read.get_tag("NH") != 1:
                    continue
            except KeyError:
                pass  # velocyto needs NH; a read without it is taken as unique
            try:
                barcode = read.get_tag(CB_tag)
                umi = read.get_tag(UB_tag)
            except KeyError:
                continue
            row = tag_rows.get(barcode, -1)
            if row == -1:
                row = tag_rows[barcode] = barcodes.get(barcode.split("-")[0])
            if row is None:
                continue

            if reference_start != cache_position:
                cache.clear()
                cache_position = reference_start
            is_reverse = read.is_reverse
            key = (is_reverse, read.cigarstring)
            if key in cache:
                record = cache[key]
            else:
                segments, spliced = parse_cigar(read.cigartuples, reference_start + 1)
                if not segments or segments[-1][1] - segments[0][0] + 1 > LONGEST_READ_SPAN:
                    record = None
                elif not has_features:
                    record = EMPTY
                else:
                    strand = strands[is_reverse]
                    record = strand.match(reference_start + 1, segments, spliced, validated)
                    if record and strand.masked(segments):
                        record = EMPTY
                cache[key] = record
            if record is None:
                continue
            seen.add(row)
            if not record:
                continue

            molecule = (row, umi)
            current = molecules.get(molecule)
            if current is None:
                molecules[molecule] = record
            elif current is not record:
                pair = (id(current), id(record))
                combined = merged_records.get(pair)
                if combined is None:
                    combined = _merge(current, record)
                    combined = merged_records[pair] = records.setdefault(combined, combined)
                molecules[molecule] = combined

    counts: Dict[int, List[int]] = {}
    classes: Dict[int, Optional[int]] = {}
    for (row, _), record in molecules.items():
        key = id(record)
        if key not in classes:
            classes[key] = classify(record, validated)
        layer = classes[key]
        if layer is not None:
            counts.setdefault(row, [0, 0, 0])[layer] += 1

    rows = np.fromiter(counts.keys(), dtype=np.int64, count=len(counts))
    values = np.array(list(counts.values()), dtype=np.int64).reshape(-1, 3)
    return np.fromiter(seen, dtype=np.int64, count=len(seen)), rows, values


def count_molecules(
    bam_path: str,
    barcodes: Sequence[str],
    gtf_path: str,
    repeat_mask_path: Optional[str] = None,
    regions: Optional[List[tuple]] = None,
    tiles: int = 100,
    cores: int = 1,
    CB_tag: str = "CB",
    UB_tag: str = "UB",
    bam_index_path: Optional[str] = None,
) -> pd.DataFrame:
    """
    Count the spliced, unspliced and ambiguous molecules of each barcode with velocyto's rules.

    By default the genome is split into tiles with boundaries between genes, and every read is used once, in the
    tile that contains its start. If regions are given, every read overlapping a region is used for that region,
    so a molecule overlapping two regions is counted in both.

    Parameters:
        bam_path (str): Position-sorted BAM file with cell barcode and UMI tags.
        barcodes (Sequence[str]): Barcodes to count, with or without a suffix such as '-1'.
        gtf_path (str): GTF file of the reference the reads were aligned to.
        repeat_mask_path (str, optional): GTF file of repeats to mask.
        regions (List[tuple], optional): Genomic intervals (contig, start, end) to process.
        tiles (int, optional): Number of genomic regions to process in parallel.
        cores (int, optional): Number of worker processes.
        CB_tag (str, optional): Tag of the cell barcode.
        UB_tag (str, optional): Tag of the UMI.
        bam_index_path (str, optional): BAM index. Found automatically if omitted.

    Returns:
        pd.DataFrame: Molecules per class (columns of LAYERS), indexed by barcode without suffix, for the barcodes
            with at least one read that velocyto reads (unique, with both tags), in the order given.
    """
    import pysam

    names = list(dict.fromkeys(str(barcode).split("-")[0] for barcode in barcodes))
    rows = {name: row for row, name in enumerate(names)}

    with ExitStack() as stack:
        executor = None
        if cores > 1:
            executor = stack.enter_context(ProcessPoolExecutor(max_workers=cores))
            # Where workers are forked, they start now, before the gene models and repeats fill this process's memory
            executor.submit(int).result()

        logger.info("Reading the gene models of %s.", gtf_path)
        models = read_gene_models(gtf_path)
        masks = None
        if repeat_mask_path is not None:
            logger.info("Reading the repeats of %s.", repeat_mask_path)
            masks = read_repeat_mask(repeat_mask_path)

        count_by_start = regions is None
        if regions is None:
            with pysam.AlignmentFile(bam_path, "rb", index_filename=bam_index_path) as bam_file:
                references, lengths = bam_file.references, bam_file.lengths
            if sum(lengths) <= 0:
                raise ValueError(f"The header of '{bam_path}' lists no reference sequences.")
            regions = tiles_between_genes(references, lengths, models["blocks"], tiles)

        def jobs():
            for region in regions:
                features, intervals = tile_inputs(models, masks, *region, count_by_start)
                yield (region, bam_path, rows, features, intervals, CB_tag, UB_tag, bam_index_path, count_by_start)

        seen = np.zeros(len(names), dtype=bool)
        counts = np.zeros((len(names), len(LAYERS)), dtype=np.int64)

        def add(result):
            seen_rows, count_rows, values = result
            seen[seen_rows] = True
            counts[count_rows] += values  # a barcode appears once per tile

        progress = tqdm(total=len(regions), desc="Processing BAM tiles", unit="tile")
        if executor is None:
            for job in jobs():
                add(count_tile(*job))
                progress.update()
        else:
            # Tiles are prepared as workers become free, so that their inputs are not all held at once
            pending = set()
            for job in jobs():
                pending.add(executor.submit(count_tile, *job))
                if len(pending) >= 2 * cores:
                    done, pending = wait(pending, return_when=FIRST_COMPLETED)
                    for future in done:
                        add(future.result())
                        progress.update()
            for future in wait(pending).done:
                add(future.result())
                progress.update()
        progress.close()

        return pd.DataFrame(counts[seen], index=pd.Index(np.array(names, dtype=object)[seen]), columns=list(LAYERS))
