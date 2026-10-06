"""Command line interface: `qclus run`, `qclus splicing-from-bam` and `qclus splicing-from-loom`."""
import argparse
import inspect
import json
import logging
import os
import sys

import pandas as pd

from qclus._version import __version__
from qclus.qclus import TISSUE_PRESETS, quickstart_qclus, run_qclus
from qclus.utils import fraction_unspliced_from_bam, fraction_unspliced_from_loom

# Order in which the number of barcodes per label is printed
LABELS = ["passed", "initial filter", "clustering filter", "outlier filter", "scrublet filter"]


def default_text(name: str) -> str:
    """Describe the default of a run_qclus argument, including what each --tissue sets."""
    def show(value):
        return " ".join(map(str, value)) if isinstance(value, (list, tuple)) else str(value)

    default = inspect.signature(run_qclus).parameters[name].default
    if any(name in settings for settings in TISSUE_PRESETS.values()):
        by_tissue = "; ".join(f"{tissue}: {show(settings.get(name, default))}" for tissue, settings in TISSUE_PRESETS.items())
        return f"Default set by --tissue ({by_tissue})."
    return f"Default: {show(default)}."


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="qclus",
        description="QClus: droplet filtering for single-nucleus RNA-seq data.",
    )
    parser.add_argument("--version", action="version", version=f"qclus {__version__}")
    commands = parser.add_subparsers(dest="command", required=True, metavar="COMMAND")

    # Options shared by every command
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--overwrite", action="store_true", help="Replace output files that already exist.")
    common.add_argument("--debug", action="store_true", help="Show the full traceback when an error occurs.")

    # ------------------------------------------------------------------ qclus run
    run = commands.add_parser(
        "run",
        parents=[common],
        help="Run QClus on a count matrix.",
        description="Run QClus on one sample and write the annotated result.",
    )
    run.set_defaults(handler=run_command)

    files = run.add_argument_group("input and output")
    files.add_argument("--counts", required=True, metavar="FILE", help="Count matrix, .h5 (10x Genomics) or .h5ad.")
    files.add_argument(
        "--fraction-unspliced", required=True, metavar="FILE",
        help="Table (.csv or .tsv, optionally gzipped) with barcodes in the first column and the fraction of "
             "unspliced reads in a column named fraction_unspliced, or in its only other column.",
    )
    files.add_argument(
        "-o", "--output", metavar="FILE",
        help="Annotated .h5ad file. It holds the barcodes that have splicing information.",
    )
    files.add_argument(
        "--obs-csv", metavar="FILE",
        help="Per-barcode table (.csv or .tsv) with the QClus label and all metrics. At least one of --output and "
             "--obs-csv is required.",
    )
    files.add_argument(
        "--passed-only", action="store_true",
        help="Write only barcodes labelled 'passed' to --output. The table of --obs-csv always holds every barcode.",
    )

    run.add_argument(
        "--tissue", choices=sorted(TISSUE_PRESETS), default="heart",
        help="Workflow preset. 'heart' uses the clustering settings of the publication; 'other' leaves out the "
             "cell type specific metrics. Default: heart.",
    )

    # One option per argument of run_qclus. None has a default here: what is not given comes from --tissue,
    # and then from run_qclus itself. tests/test_cli.py checks that this list matches the signature.
    options = run.add_argument_group("pipeline options")
    unset = argparse.SUPPRESS
    options.add_argument(
        "--nucl-gene-set", dest="nucl_gene_set", default=unset, metavar="FILE",
        help="Text file with one nuclear gene per line. Default: the built-in set of 30 genes.",
    )
    options.add_argument(
        "--celltype-gene-sets", dest="celltype_gene_set_dict", default=unset, metavar="FILE",
        help="JSON file mapping cell type names to lists of marker genes. Default: the built-in cardiac sets.",
    )
    options.add_argument(
        "--minimum-genes", type=int, default=unset, metavar="N",
        help=f"Fewest genes detected for a barcode to pass the initial filter. {default_text('minimum_genes')}",
    )
    options.add_argument(
        "--maximum-genes", type=int, default=unset, metavar="N",
        help=f"Most genes detected for a barcode to pass the initial filter. {default_text('maximum_genes')}",
    )
    options.add_argument(
        "--max-mito-perc", type=float, default=unset, metavar="PERCENT",
        help=f"Highest mitochondrial percentage to pass the initial filter. {default_text('max_mito_perc')}",
    )
    options.add_argument(
        "--clustering-features", nargs="+", default=unset, metavar="FEATURE",
        help=f"Features used for clustering; fraction_unspliced must be one. {default_text('clustering_features')}",
    )
    options.add_argument(
        "--clustering-k", type=int, default=unset, metavar="K",
        help=f"Number of k-means clusters. {default_text('clustering_k')}",
    )
    options.add_argument(
        "--clusters-to-select", nargs="+", default=unset, metavar="CLUSTER",
        help="Clusters that pass the clustering filter. Clusters are numbered from 0 by decreasing mean "
             f"fraction_unspliced. {default_text('clusters_to_select')}",
    )
    options.add_argument(
        "--scrublet-filter", action=argparse.BooleanOptionalAction, default=unset,
        help=f"Filter doublets with Scrublet. {default_text('scrublet_filter')}",
    )
    options.add_argument(
        "--scrublet-expected-rate", type=float, default=unset, metavar="RATE",
        help=f"Expected doublet rate. {default_text('scrublet_expected_rate')}",
    )
    options.add_argument(
        "--scrublet-minimum-counts", type=int, default=unset, metavar="N",
        help=f"Scrublet gene filter: minimum counts. {default_text('scrublet_minimum_counts')}",
    )
    options.add_argument(
        "--scrublet-minimum-cells", type=int, default=unset, metavar="N",
        help=f"Scrublet gene filter: minimum cells. {default_text('scrublet_minimum_cells')}",
    )
    options.add_argument(
        "--scrublet-minimum-gene-variability-pctl", type=float, default=unset, metavar="PERCENTILE",
        help=f"Scrublet gene filter: variability percentile. {default_text('scrublet_minimum_gene_variability_pctl')}",
    )
    options.add_argument(
        "--scrublet-n-pcs", type=int, default=unset, metavar="N",
        help=f"Principal components used by Scrublet. {default_text('scrublet_n_pcs')}",
    )
    options.add_argument(
        "--scrublet-thresh", type=float, default=unset, metavar="SCORE",
        help=f"Doublet score at or above which a barcode is removed. {default_text('scrublet_thresh')}",
    )
    options.add_argument(
        "--scrublet-approx-neighbors", action=argparse.BooleanOptionalAction, default=unset,
        help="Use approximate nearest neighbours (annoy) in Scrublet, as QClus 0.2.0 and earlier did. "
             f"{default_text('scrublet_approx_neighbors')}",
    )
    options.add_argument(
        "--outlier-filter", action=argparse.BooleanOptionalAction, default=unset,
        help=f"Filter outliers within the selected clusters. {default_text('outlier_filter')}",
    )
    options.add_argument(
        "--outlier-unspliced-diff", type=float, default=unset, metavar="DIFF",
        help="Subtracted from the 25th percentile of fraction_unspliced in the reference cluster to get the "
             f"lower threshold. {default_text('outlier_unspliced_diff')}",
    )
    options.add_argument(
        "--outlier-mito-diff", type=float, default=unset, metavar="DIFF",
        help="Added to the 75th percentile of the mitochondrial percentage in the reference cluster to get the "
             f"upper threshold. {default_text('outlier_mito_diff')}",
    )
    options.add_argument(
        "--kmeans-n-init", type=int, default=unset, metavar="N",
        help=f"Number of k-means restarts. {default_text('kmeans_n_init')}",
    )
    options.add_argument(
        "--compute-embedding", action=argparse.BooleanOptionalAction, default=unset,
        help="Compute the UMAP of the clustering features, which is used only for plotting. "
             f"{default_text('compute_embedding')}",
    )

    # ------------------------------------------------------------------ qclus splicing-from-bam
    bam_defaults = inspect.signature(fraction_unspliced_from_bam).parameters
    bam = commands.add_parser(
        "splicing-from-bam",
        parents=[common],
        help="Calculate fraction_unspliced from a 10X BAM file.",
        description="Calculate the fraction of unspliced reads per barcode from the region tags of a 10X BAM file.",
    )
    bam.set_defaults(handler=bam_command)
    bam.add_argument("--bam", required=True, metavar="FILE", help="Position-sorted BAM file from Cell Ranger.")
    bam.add_argument("--barcodes", required=True, metavar="FILE", help="Barcodes to count reads for, one per line (barcodes.tsv.gz).")
    bam.add_argument("-o", "--output", required=True, metavar="FILE", help="Table to write (.csv or .tsv), as read by 'qclus run'.")
    bam.add_argument("--bam-index", metavar="FILE", help="BAM index. Default: found next to the BAM file.")
    bam.add_argument("--tiles", type=int, default=bam_defaults["tiles"].default, metavar="N",
                     help="Number of genomic tiles processed in parallel. Default: %(default)s.")
    bam.add_argument("--cores", type=int, metavar="N", help="CPU cores to use. Default: all but one of those available.")
    bam.add_argument("--cb-tag", default=bam_defaults["CB_tag"].default, metavar="TAG", help="Cell barcode tag. Default: %(default)s.")
    bam.add_argument("--re-tag", default=bam_defaults["RE_tag"].default, metavar="TAG", help="Region type tag. Default: %(default)s.")
    bam.add_argument("--exon-tag", default=bam_defaults["EXON_tag"].default, metavar="VALUE",
                     help="Value of the region type tag for exonic reads. Default: %(default)s.")
    bam.add_argument("--intron-tag", default=bam_defaults["INTRON_tag"].default, metavar="VALUE",
                     help="Value of the region type tag for intronic reads. Default: %(default)s.")

    # ------------------------------------------------------------------ qclus splicing-from-loom
    loom = commands.add_parser(
        "splicing-from-loom",
        parents=[common],
        help="Calculate fraction_unspliced from a Velocyto loom file.",
        description="Calculate the fraction of unspliced reads per barcode from the layers of a Velocyto loom file.",
    )
    loom.set_defaults(handler=loom_command)
    loom.add_argument("--loom", required=True, metavar="FILE", help="Loom file written by Velocyto.")
    loom.add_argument("-o", "--output", required=True, metavar="FILE", help="Table to write (.csv or .tsv), as read by 'qclus run'.")

    return parser


def pipeline_options(args: argparse.Namespace) -> dict:
    """The run_qclus arguments that were given on the command line."""
    names = [name for name in inspect.signature(run_qclus).parameters if name not in ("counts_path", "fraction_unspliced")]
    return {name: getattr(args, name) for name in names if hasattr(args, name)}


def same_file(first: str, second: str) -> bool:
    """Whether two paths lead to the same file, also through symbolic or hard links."""
    if os.path.exists(first) and os.path.exists(second):
        return os.path.samefile(first, second)
    return os.path.realpath(first) == os.path.realpath(second)


def check_paths(parser: argparse.ArgumentParser, inputs: dict, outputs: dict, overwrite: bool) -> None:
    """Stop with a usage error, before anything is computed, if a file cannot be read or must not be written."""
    for option, path in inputs.items():
        if not os.path.isfile(path):
            parser.error(f"{option}: '{path}' does not exist")

    checked = {}
    for option, path in outputs.items():
        if path is None:
            continue
        for input_option, input_path in inputs.items():
            if same_file(path, input_path):
                parser.error(f"{option}: '{path}' is the file given as {input_option}; inputs are never overwritten")
        for other_option, other_path in checked.items():
            if same_file(path, other_path):
                parser.error(f"{option} and {other_option} point to the same file '{path}'")
        checked[option] = path

        target = os.path.realpath(path)
        if os.path.isdir(target):
            parser.error(f"{option}: '{path}' is a directory")
        if not os.path.isdir(os.path.dirname(target)):
            parser.error(f"{option}: the directory '{os.path.dirname(target)}' does not exist")
        if os.path.exists(target) and not overwrite:
            parser.error(f"{option}: '{path}' already exists; use --overwrite to replace it")


def separator(path: str) -> str:
    return "\t" if ".tsv" in os.path.basename(path).lower() else ","


def read_gene_sets(parser: argparse.ArgumentParser, path: str, writes_h5ad: bool) -> dict:
    with open(path) as handle:
        try:
            gene_sets = json.load(handle)
        except json.JSONDecodeError as error:
            parser.error(f"--celltype-gene-sets: '{path}' is not valid JSON ({error})")
    valid = (
        isinstance(gene_sets, dict)
        and len(gene_sets) > 0
        and all(isinstance(genes, list) and all(isinstance(gene, str) for gene in genes) for genes in gene_sets.values())
    )
    if not valid:
        parser.error(f"--celltype-gene-sets: '{path}' must hold a JSON object that maps names to lists of gene names")
    with_slash = [name for name in gene_sets if "/" in name]
    if writes_h5ad and with_slash:
        parser.error(
            f"--celltype-gene-sets: the names {with_slash} contain '/', which cannot be stored in an .h5ad file; "
            "rename them or write only --obs-csv"
        )
    return gene_sets


def run_command(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    if args.output is None and args.obs_csv is None:
        parser.error("at least one of --output and --obs-csv is required")
    if args.passed_only and args.output is None:
        parser.error("--passed-only applies to --output, which was not given")

    options = pipeline_options(args)
    inputs = {"--counts": args.counts, "--fraction-unspliced": args.fraction_unspliced}
    if "nucl_gene_set" in options:
        inputs["--nucl-gene-set"] = options["nucl_gene_set"]
    if "celltype_gene_set_dict" in options:
        inputs["--celltype-gene-sets"] = options["celltype_gene_set_dict"]
    check_paths(parser, inputs, {"--output": args.output, "--obs-csv": args.obs_csv}, args.overwrite)

    if "nucl_gene_set" in options:
        with open(options["nucl_gene_set"]) as handle:
            options["nucl_gene_set"] = [line.strip() for line in handle if line.strip()]
    if "celltype_gene_set_dict" in options:
        options["celltype_gene_set_dict"] = read_gene_sets(parser, options["celltype_gene_set_dict"], args.output is not None)

    fraction_unspliced = pd.read_csv(args.fraction_unspliced, index_col=0, sep=separator(args.fraction_unspliced))
    adata = quickstart_qclus(args.counts, fraction_unspliced, tissue=args.tissue, **options)

    if args.obs_csv is not None:
        adata.obs.to_csv(args.obs_csv, sep=separator(args.obs_csv), index_label="barcode")
    if args.output is not None:
        if args.passed_only:
            adata = adata[adata.obs["qclus"] == "passed"].copy()
            # This array has one row per barcode that reached clustering; obsm["QClus_umap"] stays aligned
            adata.uns.pop("QClus_umap", None)
        adata.write_h5ad(args.output)

    # Counts of the whole run, also when only the passed barcodes were written
    counts = adata.uns["qclus"]["n_barcodes"]
    print(f"barcodes\t{sum(counts.values())}")
    for label in LABELS:
        if label in counts:
            print(f"{label}\t{counts[label]}")
    return 0


def bam_command(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    inputs = {"--bam": args.bam, "--barcodes": args.barcodes}
    if args.bam_index is not None:
        inputs["--bam-index"] = args.bam_index
    else:
        # The index that pysam finds next to the BAM file is an input too, and is protected like one
        stem = os.path.splitext(args.bam)[0]
        for candidate in (args.bam + ".bai", args.bam + ".csi", stem + ".bai", stem + ".csi"):
            if os.path.isfile(candidate):
                inputs[f"the index of --bam ({os.path.basename(candidate)})"] = candidate
    check_paths(parser, inputs, {"--output": args.output}, args.overwrite)
    if args.tiles < 1:
        parser.error("--tiles must be at least 1")

    result = fraction_unspliced_from_bam(
        bam_path=args.bam,
        bam_index_path=args.bam_index,
        barcodes_path=args.barcodes,
        tiles=args.tiles,
        cores=args.cores,
        CB_tag=args.cb_tag,
        RE_tag=args.re_tag,
        EXON_tag=args.exon_tag,
        INTRON_tag=args.intron_tag,
    )
    result.to_csv(args.output, sep=separator(args.output))
    print(f"barcodes\t{len(result)}")
    return 0


def loom_command(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    check_paths(parser, {"--loom": args.loom}, {"--output": args.output}, args.overwrite)
    result = fraction_unspliced_from_loom(args.loom)
    result.to_csv(args.output, sep=separator(args.output))
    print(f"barcodes\t{len(result)}")
    return 0


def main(argv=None) -> int:
    """
    Run the command line interface.

    Returns the exit status: 0 on success and 1 on an error while running. Mistakes in the arguments,
    including files that are missing or would be overwritten, exit with status 2 before anything is computed.
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    # Progress messages of the library go to stderr; results go to stdout
    logger = logging.getLogger("qclus")
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("%(message)s"))
        logger.addHandler(handler)
    logger.setLevel(logging.INFO)

    try:
        return args.handler(args, parser)
    except Exception as error:
        if args.debug:
            raise
        print(f"qclus: error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
