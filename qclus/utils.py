from sklearn.cluster import KMeans
from sklearn.preprocessing import MinMaxScaler
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import os
import warnings
import scanpy as sc
from typing import Iterable, List, Optional, Sequence, Tuple, Union
import numpy as np
import pandas as pd
from anndata import AnnData

# loompy, pysam, scrublet and umap are imported inside the functions that use them:
# together they make up most of the import time, and most calls need only some of them.


def check_unique_barcodes(barcodes: pd.Index, source: str) -> None:
    """
    Raise an error if cell barcodes are not unique after truncation.

    Parameters:
        barcodes (pd.Index): Truncated cell barcodes.
        source (str): Where the barcodes come from, used in the error message.

    Raises:
        ValueError: If any barcode occurs more than once.
    """
    barcodes = pd.Index(barcodes)
    duplicated = barcodes[barcodes.duplicated()].unique()
    if len(duplicated) > 0:
        raise ValueError(
            f"{len(duplicated)} barcodes in {source} are not unique after truncation to 16 characters "
            f"(for example '{duplicated[0]}'). QClus runs one sample at a time."
        )


def prepare_fraction_unspliced(fraction_unspliced: Union[pd.Series, pd.DataFrame]) -> pd.Series:
    """
    Validate the fraction of unspliced reads and index it by truncated cell barcodes.

    Parameters:
        fraction_unspliced (pd.Series or pd.DataFrame): Fraction of unspliced reads per cell. A DataFrame
            must have a 'fraction_unspliced' column or exactly one column.

    Returns:
        pd.Series: Fraction of unspliced reads, indexed by cell barcodes truncated to 16 characters.

    Raises:
        TypeError: If the input is neither a Series nor a DataFrame.
        ValueError: If the column cannot be identified, barcodes are not unique after truncation,
            or values are missing, not numeric, or outside 0 to 1.
    """
    if isinstance(fraction_unspliced, pd.DataFrame):
        if "fraction_unspliced" in fraction_unspliced.columns:
            series = fraction_unspliced["fraction_unspliced"]
        elif fraction_unspliced.shape[1] == 1:
            series = fraction_unspliced.iloc[:, 0]
        else:
            raise ValueError(
                f"fraction_unspliced has {fraction_unspliced.shape[1]} columns and none is named "
                "'fraction_unspliced'."
            )
    elif isinstance(fraction_unspliced, pd.Series):
        series = fraction_unspliced
    else:
        raise TypeError(
            f"fraction_unspliced must be a pandas Series or DataFrame, got {type(fraction_unspliced)}."
        )

    series = series.rename("fraction_unspliced")
    series.index = pd.Index(create_new_index(series.index))
    check_unique_barcodes(series.index, "fraction_unspliced")

    if not pd.api.types.is_numeric_dtype(series):
        raise ValueError(f"fraction_unspliced must be numeric, got dtype {series.dtype}.")
    if series.isna().any():
        raise ValueError(f"fraction_unspliced has {int(series.isna().sum())} missing values.")
    if ((series < 0) | (series > 1)).any():
        raise ValueError(
            "fraction_unspliced must lie between 0 and 1, but its values range from "
            f"{series.min()} to {series.max()}."
        )
    return series


def add_fraction_unspliced(
        adata: AnnData,
        fraction_unspliced: pd.Series
) -> Optional[AnnData]:
    """
    Filters an AnnData object to keep only cells that have corresponding splicing information available.
    Adds the 'fraction_unspliced' annotation to the filtered AnnData object.

    Parameters:
        adata (AnnData): The input AnnData object.
        fraction_unspliced (pd.Series): Series containing the fraction of unspliced reads per cell.

    Returns:
        AnnData: Filtered AnnData object with the 'fraction_unspliced' annotation.

    Raises:
        ValueError: If no common barcodes are found between counts data and fraction_unspliced.
        Optional[AnnData]: None if no matching cells remain after filtering.
    """
    # Find common cell barcodes between the AnnData object and the `fraction_unspliced` series
    common_barcodes = adata.obs.index.intersection(fraction_unspliced.index)

    if len(common_barcodes) == 0:
        raise ValueError("No common barcodes found between counts data and fraction_unspliced.")

    if len(common_barcodes) < len(adata.obs.index):
        warnings.warn(
            f"Removing {len(adata.obs.index) - len(common_barcodes)} barcodes without splicing information."
        )

    # Filter adata
    adata = adata[common_barcodes].copy()

    # Add fraction_unspliced as an observation annotation
    adata.obs["fraction_unspliced"] = fraction_unspliced.loc[common_barcodes]

    return adata


def read_count_file(file_path: str) -> AnnData:
    """
    Load a counts file as an AnnData object. Supports .h5 and .h5ad file formats.

    Parameters:
        file_path (str): Path to the counts file.

    Returns:
        AnnData: Loaded AnnData object.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file format is not supported.
        IOError: If there is an error reading the file.
    """
    file_path = os.fspath(file_path)
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The counts file '{file_path}' does not exist.")

    if file_path.endswith('.h5'):
        reader = sc.read_10x_h5
    elif file_path.endswith('.h5ad'):
        reader = sc.read_h5ad
    else:
        raise ValueError(
            f"Unsupported file format for '{file_path}'. Only .h5 and .h5ad are supported."
        )

    try:
        adata = reader(file_path)
    except Exception as e:
        raise IOError(f"Failed to read counts file at '{file_path}': {e}") from e

    # Ensure unique variable names
    adata.var_names_make_unique()
    return adata


def warn_if_not_counts(adata: AnnData) -> None:
    """
    Warn if adata.X does not look like raw counts.

    Parameters:
        adata (AnnData): AnnData object whose X is checked.
    """
    # A sparse matrix keeps its non-zero values in .data
    values = adata.X.data if hasattr(adata.X, "nnz") else np.asarray(adata.X)
    sample = np.ravel(values)[:100_000]
    if sample.size > 0 and not np.all(np.mod(sample, 1) == 0):
        warnings.warn(
            "The count matrix holds values that are not whole numbers. QClus expects raw counts; "
            "normalized or transformed data gives misleading QC metrics and doublet scores."
        )



def get_qc_metrics(
    adata: sc.AnnData,
    gene_set: Sequence[str],
    key_name: str,
    normlog: bool = False,
    scale: bool = False,
) -> None:
    """
    Calculate QC metrics for a given gene set and add them to the AnnData object.

    Parameters:
        adata (sc.AnnData): AnnData object containing single-cell data.
        gene_set (Sequence[str]): Genes to calculate QC metrics for.
        key_name (str): Key name under which to store the metrics.
        normlog (bool, optional): Whether to normalize and log-transform the data.
        scale (bool, optional): Whether to scale the data.

    Raises:
        TypeError: If gene_set is a single string.
        ValueError: If none of the genes are in adata.var_names.
    """
    # A single string would be read as a sequence of one-letter gene names
    if isinstance(gene_set, str):
        raise TypeError("gene_set must be a sequence of gene names, not a single string.")
    gene_set = list(gene_set)

    # Check if genes are in adata.var_names
    missing_genes = [gene for gene in gene_set if gene not in adata.var_names]
    if len(missing_genes) == len(gene_set):
        raise ValueError(
            f"None of the {len(gene_set)} genes in gene set '{key_name}' are in adata.var_names. "
            "The gene sets and adata.var_names must use the same gene names; "
            "the built-in gene sets are human gene symbols."
        )
    if missing_genes:
        warnings.warn(
            f"The following genes of gene set '{key_name}' are not in adata.var_names and will be ignored: "
            f"{missing_genes}"
        )

    # Create a boolean mask for the genes in the gene_set
    adata.var[key_name] = adata.var_names.isin(gene_set)

    # Create a layer to store the modified data
    adata.layers[key_name] = adata.X.copy()

    # Use the new layer for normalization and scaling
    if normlog:
        sc.pp.normalize_total(adata, target_sum=1e4, layer=key_name)
        sc.pp.log1p(adata, layer=key_name)

    if scale:
        sc.pp.scale(adata, max_value=10, layer=key_name)

    # Calculate QC metrics on the specified layer
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=[key_name],
        percent_top=None,
        log1p=False,
        inplace=True,
        layer=key_name,
    )
    sc.tl.score_genes(adata, gene_list=gene_set, score_name=f"score_{key_name}")

    # Remove the layer to save memory
    del adata.layers[key_name]


def add_qclus_embedding(
    adata: sc.AnnData,
    features: List[str],
    random_state: int = 1,
    n_components: int = 2,
) -> np.ndarray:
    """
    Compute UMAP embedding using specified features and add it to the AnnData object.

    Parameters:
        adata (sc.AnnData): AnnData object containing single-cell data.
        features (List[str]): List of features to use for embedding.
        random_state (int, optional): Random state for reproducibility.
        n_components (int, optional): Number of UMAP components.

    Returns:
        np.ndarray: UMAP embedding coordinates.
    """
    # Check if features are available
    missing_features = [feat for feat in features if feat not in adata.obs.columns]
    if missing_features:
        raise ValueError(f"The following features are missing in adata.obs: {missing_features}")

    # Extract the features
    X_full = adata.obs[features]

    # Scale the features
    scaler = MinMaxScaler()
    X_scaled = scaler.fit_transform(X_full)

    # Compute the UMAP embedding
    from umap import UMAP

    umap_embedding = UMAP(random_state=random_state, n_components=n_components).fit_transform(X_scaled)

    return umap_embedding


def create_new_index(index: pd.Index) -> List[str]:
    """
    Truncate cell barcodes to the first 16 characters.

    Parameters:
        index (pd.Index): Index of cell barcodes.

    Returns:
        List[str]: List of truncated cell barcodes.
    """
    return [str(x)[:16] for x in index]


def fraction_unspliced_from_loom(loompy_path: str, batch_size: int = 512) -> pd.DataFrame:
    """
    Calculate the fraction of unspliced reads per cell from a Loom file.

    Parameters:
        loompy_path (str): Path to the Loom file.
        batch_size (int, optional): Number of cells read into memory at a time.

    Returns:
        pd.DataFrame: DataFrame containing the fraction of unspliced reads per cell.

    Raises:
        FileNotFoundError: If the Loom file does not exist.
        ValueError: If the cell IDs are not in Velocyto's 'sample:barcode' form.
    """
    import loompy

    if not os.path.exists(loompy_path):
        raise FileNotFoundError(f"The Loom file '{loompy_path}' does not exist.")

    # The file is only read, so it is opened read-only (loompy's default mode needs write access)
    with loompy.connect(loompy_path, mode='r') as loompy_con:
        cell_ids = [str(x) for x in loompy_con.ca['CellID']]
        malformed = [x for x in cell_ids if ':' not in x]
        if malformed:
            raise ValueError(
                f"{len(malformed)} cell IDs in '{loompy_path}' are not in Velocyto's 'sample:barcode' form "
                f"(for example '{malformed[0]}')."
            )
        barcodes = [x.split(':')[1][:16] for x in cell_ids]

        # Sum each layer over batches of cells instead of loading it whole
        n_cells = loompy_con.shape[1]
        counts = {}
        for layer in ('spliced', 'unspliced', 'ambiguous'):
            counts[layer] = np.concatenate([
                loompy_con.layers[layer][:, start:start + batch_size].sum(axis=0)
                for start in range(0, n_cells, batch_size)
            ])

        total_counts = counts['spliced'] + counts['unspliced'] + counts['ambiguous']
        with np.errstate(invalid='ignore'):
            fraction_unspliced = counts['unspliced'] / total_counts
        fraction_unspliced = np.nan_to_num(fraction_unspliced)  # Replace NaN with zero

        return pd.DataFrame({'fraction_unspliced': fraction_unspliced}, index=barcodes)


def annoy_works() -> bool:
    """
    Check that the annoy library returns sensible neighbours in this environment.

    Scrublet's approximate neighbour search relies on annoy. Some builds of annoy return wrong
    neighbours, in which case every cell gets the same doublet score.

    Returns:
        bool: True if random vectors are their own nearest neighbour, as they must be.
    """
    try:
        from annoy import AnnoyIndex
    except ImportError:
        return False

    vectors = np.random.default_rng(0).normal(size=(200, 30)).astype("float32")
    index = AnnoyIndex(vectors.shape[1], "euclidean")
    for i, vector in enumerate(vectors):
        index.add_item(i, vector.tolist())
    index.build(10)
    hits = sum(index.get_nns_by_item(i, 1)[0] == i for i in range(len(vectors)))
    return hits >= 0.9 * len(vectors)


def calculate_scrublet(
    adata: sc.AnnData,
    expected_rate: float = 0.06,
    minimum_counts: int = 2,
    minimum_cells: int = 3,
    minimum_gene_variability_pctl: float = 85.0,
    n_pcs: int = 30,
    thresh: float = 0.1,
    approx_neighbors: bool = False,
) -> np.ndarray:
    """
    Calculate doublet scores using Scrublet.

    Parameters:
        adata (sc.AnnData): AnnData object containing single-cell data.
        expected_rate (float, optional): Expected doublet rate.
        minimum_counts (int, optional): Minimum counts per cell.
        minimum_cells (int, optional): Minimum cells per gene.
        minimum_gene_variability_pctl (float, optional): Minimum gene variability percentile.
        n_pcs (int, optional): Number of principal components.
        thresh (float, optional): Not used here; the threshold is applied to the returned scores
            by run_qclus. Kept so that existing calls keep working.
        approx_neighbors (bool, optional): Whether to use approximate nearest neighbours (annoy)
            instead of exact ones.

    Returns:
        np.ndarray: Array of doublet scores.

    Raises:
        RuntimeError: If approximate neighbours are requested but annoy is broken in this
            environment, or if the scores are not finite or are the same for every cell.
        ValueError: If Scrublet itself fails, for example because there are too few cells or genes.
    """
    if approx_neighbors and not annoy_works():
        raise RuntimeError(
            "The annoy library returns wrong neighbours in this environment, so Scrublet's approximate "
            "neighbour search cannot be used. Use exact neighbours (scrublet_approx_neighbors=False), "
            "or install python-annoy from conda-forge."
        )

    import scrublet as scr

    scrub = scr.Scrublet(adata.X, expected_doublet_rate=expected_rate)
    try:
        scrub.scrub_doublets(
            min_counts=minimum_counts,
            min_cells=minimum_cells,
            min_gene_variability_pctl=minimum_gene_variability_pctl,
            n_prin_comps=n_pcs,
            use_approx_neighbors=approx_neighbors,
            verbose=False,
        )
    except Exception as e:
        raise ValueError(
            f"Scrublet failed on {adata.n_obs} barcodes with scrublet_n_pcs={n_pcs}: {e}. "
            "Scrublet needs more barcodes, and more genes passing its own gene filter, than principal "
            "components. Lower scrublet_n_pcs, or skip the doublet filter with scrublet_filter=False."
        ) from e

    doublet_scores = np.asarray(scrub.doublet_scores_obs_, dtype=float)
    if not np.all(np.isfinite(doublet_scores)):
        raise RuntimeError(
            f"Scrublet returned {int((~np.isfinite(doublet_scores)).sum())} doublet scores that are not "
            "finite numbers. Skip the doublet filter with scrublet_filter=False to continue without it."
        )
    if doublet_scores.size > 0 and np.all(doublet_scores == doublet_scores[0]):
        raise RuntimeError(
            f"Scrublet gave every barcode the same doublet score ({doublet_scores[0]:.3g}), so the doublet "
            "filter cannot separate anything on this input. Skip it with scrublet_filter=False to continue."
        )
    return doublet_scores


def do_kmeans(X_full: pd.DataFrame, k: int, n_init: int = 1) -> List[str]:
    """
    Perform k-means clustering and sort clusters by decreasing mean fraction_unspliced.

    Parameters:
        X_full (pd.DataFrame): DataFrame containing features for clustering.
        k (int): Number of clusters.
        n_init (int, optional): Number of k-means restarts; the best one is kept.

    Returns:
        List[str]: List of cluster labels as strings.
    """
    if not isinstance(X_full, pd.DataFrame):
        raise TypeError("X_full must be a pandas DataFrame.")

    # Scale the features
    scaler = MinMaxScaler()
    X_scaled = scaler.fit_transform(X_full)

    # Perform k-means clustering
    kmeans = KMeans(n_clusters=k, random_state=0, n_init=n_init)
    labels = pd.Series(kmeans.fit_predict(X_scaled).astype(str), index=X_full.index)

    # Sort clusters by decreasing mean fraction_unspliced
    clusters = []
    for cluster_label in map(str, range(k)):
        cluster_df = X_full[labels == cluster_label]
        clusters.append((cluster_label, cluster_df['fraction_unspliced'].mean()))

    # A cluster that k-means left empty has no mean; it is ordered last
    sorted_clusters = sorted(clusters, key=lambda x: -np.inf if np.isnan(x[1]) else x[1], reverse=True)
    cluster_order = {cluster_label: str(idx) for idx, (cluster_label, _) in enumerate(sorted_clusters)}

    # Reassign cluster labels based on sorted order
    return labels.map(cluster_order).tolist()


def outlier_thresholds(
    df: pd.DataFrame,
    unspliced_diff: float,
    mito_diff: float,
) -> Tuple[float, float]:
    """
    Compute the outlier thresholds for fraction_unspliced and pct_counts_MT.

    The reference cluster is determined as the cluster with the highest mean
    fraction_unspliced. Outlier thresholds are then computed from this cluster.

    Parameters:
        df (pd.DataFrame): DataFrame containing 'fraction_unspliced', 'pct_counts_MT', and 'kmeans'.
        unspliced_diff (float): Value to subtract from the 25th percentile of fraction_unspliced.
        mito_diff (float): Value to add to the 75th percentile of pct_counts_MT.

    Returns:
        Tuple[float, float]: Lower threshold for fraction_unspliced and upper threshold for pct_counts_MT.
    """
    # Find the cluster with the highest mean fraction_unspliced
    cluster_means = df.groupby('kmeans')['fraction_unspliced'].mean()
    max_unspliced_cluster = cluster_means.idxmax()

    # Use the selected cluster as the reference to calculate thresholds
    ref_cluster = df[df['kmeans'] == max_unspliced_cluster]
    unspliced_threshold = ref_cluster['fraction_unspliced'].quantile(0.25) - unspliced_diff
    mito_threshold = ref_cluster['pct_counts_MT'].quantile(0.75) + mito_diff
    return unspliced_threshold, mito_threshold


def annotate_outliers(
    df: pd.DataFrame,
    unspliced_diff: float,
    mito_diff: float,
) -> pd.Series:
    """
    Annotate outlier cells based on fraction_unspliced and pct_counts_MT.

    The reference cluster is determined as the cluster with the highest mean
    fraction_unspliced. Outlier thresholds are then computed from this cluster.

    Parameters:
        df (pd.DataFrame): DataFrame containing 'fraction_unspliced', 'pct_counts_MT', and 'kmeans'.
        unspliced_diff (float): Value to subtract from the 25th percentile of fraction_unspliced.
        mito_diff (float): Value to add to the 75th percentile of pct_counts_MT.

    Returns:
        pd.Series: Boolean Series indicating outlier cells (True means outlier).
    """
    unspliced_threshold, mito_threshold = outlier_thresholds(df, unspliced_diff, mito_diff)

    # Annotate cells as outliers:
    # A cell is considered "good" (not an outlier) if its fraction_unspliced exceeds the threshold
    # and its pct_counts_MT is below the threshold. We invert this logic for the outlier flag.
    is_outlier = ~(
        (df['fraction_unspliced'] > unspliced_threshold) &
        (df['pct_counts_MT'] < mito_threshold)
    )
    return is_outlier


def available_cpus() -> int:
    """
    Number of CPU cores this process may use.

    On a compute cluster this is the number of cores the job was given, not the size of the node.
    """
    if hasattr(os, "sched_getaffinity"):
        return len(os.sched_getaffinity(0))
    return os.cpu_count() or 1


def parse_bam_tags(
    interval: tuple,
    bam_path: str,
    barcodes: Iterable[str],
    CB_tag: str,
    RE_tag: str,
    EXON_tag: str,
    INTRON_tag: str,
    bam_index_path: Optional[str] = None,
    count_by_start: bool = False,
) -> Optional[pd.DataFrame]:
    """
    Parse BAM file tags for a given genomic interval.

    Parameters:
        interval (tuple): Tuple of (contig, start, end).
        bam_path (str): Path to the BAM file.
        barcodes (Iterable[str]): Cell barcodes to count reads for.
        CB_tag (str): Tag for cell barcode in BAM file.
        RE_tag (str): Tag for region type in BAM file.
        EXON_tag (str): Tag indicating exon region.
        INTRON_tag (str): Tag indicating intron region.
        bam_index_path (str, optional): Path to the BAM index file. Found automatically if omitted.
        count_by_start (bool, optional): If True, count only reads that start inside the interval,
            so that adjacent intervals never count the same read twice. If False, count every read
            that overlaps the interval.

    Returns:
        Optional[pd.DataFrame]: DataFrame with counts of exon and intron reads per barcode.
    """
    import pysam

    if not isinstance(barcodes, (set, frozenset)):
        barcodes = set(barcodes)
    contig, start, end = interval

    counts = Counter()
    with pysam.AlignmentFile(bam_path, "rb", index_filename=bam_index_path) as bam_file:
        for read in bam_file.fetch(contig, start, end):
            if count_by_start and read.reference_start < start:
                continue
            if not read.has_tag(CB_tag) or not read.has_tag(RE_tag):
                continue
            cb_tag = read.get_tag(CB_tag)
            if cb_tag in barcodes:
                counts[(cb_tag, read.get_tag(RE_tag))] += 1

    if not counts:
        return None

    # Rows are barcodes, columns are region types
    count_data = (
        pd.Series(counts)
        .unstack(fill_value=0)
        .rename_axis(index='CB', columns='RE')
        .reindex(columns=[EXON_tag, INTRON_tag], fill_value=0)
    )
    return count_data


def fraction_unspliced_from_bam(
    bam_path: Optional[str] = None,
    bam_index_path: Optional[str] = None,
    barcodes_path: Optional[str] = None,
    regions: Optional[List[tuple]] = None,
    tiles: int = 100,
    cores: Optional[int] = None,
    CB_tag: str = "CB",
    RE_tag: str = "RE",
    EXON_tag: str = "E",
    INTRON_tag: str = "N",
) -> Optional[pd.DataFrame]:
    """
    Calculate the fraction of unspliced reads per cell from a BAM file.

    By default the genome is split into tiles that do not overlap, and every read is counted once,
    in the tile that contains its start. If regions are given, every read overlapping a region is
    counted for that region, so a read overlapping two regions is counted in both.

    Parameters:
        bam_path (str, optional): Path to the BAM file.
        bam_index_path (str, optional): Path to the BAM index file. Found automatically if omitted.
        barcodes_path (str, optional): Path to cell barcodes.
        regions (List[tuple], optional): List of genomic intervals to process.
        tiles (int, optional): Number of genomic regions to process in parallel.
        cores (int, optional): Number of CPU cores to use.
        CB_tag (str, optional): Tag for cell barcode in BAM file.
        RE_tag (str, optional): Tag for region type in BAM file.
        EXON_tag (str, optional): Tag indicating exon region.
        INTRON_tag (str, optional): Tag indicating intron region.

    Returns:
        Optional[pd.DataFrame]: DataFrame containing the fraction of unspliced reads per cell.

    Raises:
        ValueError: If a required path is missing, tiles is less than 1, the BAM header lists no
            reference sequences, or no read matches the barcodes and tags.
        FileNotFoundError: If the BAM file or the given index file does not exist.
    """
    import pysam

    if bam_path is None or barcodes_path is None:
        raise ValueError("Please provide bam_path and barcodes_path.")

    if not os.path.exists(bam_path):
        raise FileNotFoundError(f"The BAM file '{bam_path}' does not exist.")

    if bam_index_path is not None and not os.path.exists(bam_index_path):
        raise FileNotFoundError(f"The BAM index file '{bam_index_path}' does not exist.")

    if tiles < 1:
        raise ValueError(f"tiles must be at least 1, got {tiles}.")

    barcode_df = pd.read_csv(barcodes_path,
                           header=None,
                           sep="\t")
    barcodes = set(barcode_df[0])

    if cores is None:
        cores = max(1, available_cpus() - 1)

    # Reads are assigned to generated tiles by their start; user-supplied regions count overlaps
    count_by_start = regions is None
    if regions is None:
        with pysam.AlignmentFile(bam_path, "rb", index_filename=bam_index_path) as bam_file:
            references, lengths = bam_file.references, bam_file.lengths

        # Split the genome into regions
        total_length = sum(lengths)
        if total_length <= 0:
            raise ValueError(f"The header of '{bam_path}' lists no reference sequences.")
        tile_size = max(1, total_length // tiles)
        regions = []
        for contig, length in zip(references, lengths):
            for start in range(0, length, tile_size):
                end = min(start + tile_size, length)
                regions.append((contig, start, end))

    worker_args = (barcodes, CB_tag, RE_tag, EXON_tag, INTRON_tag, bam_index_path, count_by_start)
    results = []
    if cores == 1:
        for region in tqdm(regions, desc="Processing BAM tiles", unit="tile"):
            result = parse_bam_tags(region, bam_path, *worker_args)
            if result is not None:
                results.append(result)
    else:
        with ProcessPoolExecutor(max_workers=cores) as executor:
            futures = {
                executor.submit(parse_bam_tags, region, bam_path, *worker_args): region
                for region in regions
            }

            for future in tqdm(as_completed(futures), total=len(futures), desc="Processing BAM tiles", unit="tile"):
                result = future.result()
                if result is not None:
                    results.append(result)

    if not results:
        raise ValueError(
            f"No read in '{bam_path}' matched a barcode from '{barcodes_path}'. Check that the barcodes "
            f"carry the same suffix as the {CB_tag} tag in the BAM file (for example '-1'), and that the "
            f"tag names are right (CB_tag='{CB_tag}', RE_tag='{RE_tag}')."
        )

    final_df = pd.concat(results).groupby(level=0).sum()
    total_counts = final_df[INTRON_tag] + final_df[EXON_tag]
    # A barcode with neither exonic nor intronic reads gets 0
    final_df['fraction_unspliced'] = (final_df[INTRON_tag] / total_counts).fillna(0)

    final_df.index.name = None
    final_df.columns.name = None
    final_df.index = [index.split('-')[0] for index in final_df.index]

    return final_df[['fraction_unspliced']]
