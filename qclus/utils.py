import loompy
import pandas as pd
import scrublet as scr
from sklearn.cluster import KMeans
from sklearn.preprocessing import MinMaxScaler
import pysam
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import os
import scanpy as sc
from typing import List, Optional
import numpy as np
from umap import UMAP


def get_qc_metrics(
    adata: sc.AnnData,
    gene_set: List[str],
    key_name: str,
    normlog: bool = False,
    scale: bool = False,
) -> None:
    """
    Calculate QC metrics for a given gene set and add them to the AnnData object.

    Parameters:
        adata (sc.AnnData): AnnData object containing single-cell data.
        gene_set (List[str]): List of genes to calculate QC metrics for.
        key_name (str): Key name under which to store the metrics.
        normlog (bool, optional): Whether to normalize and log-transform the data.
        scale (bool, optional): Whether to scale the data.
    """
    # Check if gene_set is a list
    if not isinstance(gene_set, list):
        raise TypeError(f"gene_set must be a list, got {type(gene_set)}.")

    # Check if genes are in adata.var_names
    missing_genes = [gene for gene in gene_set if gene not in adata.var_names]
    if missing_genes:
        print(f"Warning: The following genes are not in adata.var_names and will be ignored: {missing_genes}")

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


def fraction_unspliced_from_loom(loompy_path: str) -> pd.DataFrame:
    """
    Calculate the fraction of unspliced reads per cell from a Loom file.

    Parameters:
        loompy_path (str): Path to the Loom file.

    Returns:
        pd.DataFrame: DataFrame containing the fraction of unspliced reads per cell.
    """
    if not os.path.exists(loompy_path):
        raise FileNotFoundError(f"The Loom file '{loompy_path}' does not exist.")

    with loompy.connect(loompy_path) as loompy_con:
        barcodes = [x.split(':')[1][:16] for x in loompy_con.ca['CellID']]
        spliced_counts = loompy_con.layers['spliced'][:, :].sum(axis=0)
        unspliced_counts = loompy_con.layers['unspliced'][:, :].sum(axis=0)
        ambiguous_counts = loompy_con.layers['ambiguous'][:, :].sum(axis=0)

        total_counts = spliced_counts + unspliced_counts + ambiguous_counts
        fraction_unspliced = unspliced_counts / total_counts
        fraction_unspliced = np.nan_to_num(fraction_unspliced)  # Replace NaN with zero

        return pd.DataFrame({'fraction_unspliced': fraction_unspliced}, index=barcodes)


def calculate_scrublet(
    adata: sc.AnnData,
    expected_rate: float = 0.06,
    minimum_counts: int = 2,
    minimum_cells: int = 3,
    minimum_gene_variability_pctl: float = 85.0,
    n_pcs: int = 30,
    thresh: float = 0.1,
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
        thresh (float, optional): Threshold for doublet calling.

    Returns:
        np.ndarray: Array of doublet scores.
    """
    # if not isinstance(adata.X, (np.ndarray, np.matrix)):
    #     raise TypeError("adata.X must be a numpy array or matrix.")

    scrub = scr.Scrublet(adata.X, expected_doublet_rate=expected_rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        min_counts=minimum_counts,
        min_cells=minimum_cells,
        min_gene_variability_pctl=minimum_gene_variability_pctl,
        n_prin_comps=n_pcs,
        verbose=False,
    )
    scrub.call_doublets(threshold=thresh, verbose=False)
    return scrub.doublet_scores_obs_


def do_kmeans(X_full: pd.DataFrame, k: int) -> List[str]:
    """
    Perform k-means clustering and sort clusters by decreasing mean fraction_unspliced.

    Parameters:
        X_full (pd.DataFrame): DataFrame containing features for clustering.
        k (int): Number of clusters.

    Returns:
        List[str]: List of cluster labels as strings.
    """
    if not isinstance(X_full, pd.DataFrame):
        raise TypeError("X_full must be a pandas DataFrame.")

    # Scale the features
    scaler = MinMaxScaler()
    X_scaled = scaler.fit_transform(X_full)

    # Perform k-means clustering
    kmeans = KMeans(n_clusters=k, random_state=0)
    labels = kmeans.fit_predict(X_scaled)
    X_full['kmeans'] = labels.astype(str)

    # Sort clusters by decreasing mean fraction_unspliced
    clusters = []
    for cluster_label in map(str, range(k)):
        cluster_df = X_full[X_full['kmeans'] == cluster_label]
        clusters.append((cluster_label, cluster_df['fraction_unspliced'].mean()))

    sorted_clusters = sorted(clusters, key=lambda x: x[1], reverse=True)
    cluster_order = {cluster_label: str(idx) for idx, (cluster_label, _) in enumerate(sorted_clusters)}

    # Reassign cluster labels based on sorted order
    X_full['kmeans'] = X_full['kmeans'].map(cluster_order)
    return X_full['kmeans'].tolist()


def annotate_outliers(
    df: pd.DataFrame,
    unspliced_diff: float,
    mito_diff: float,
) -> pd.Series:
    """
    Annotate outlier cells based on fraction_unspliced and pct_counts_MT.

    Parameters:
        df (pd.DataFrame): DataFrame containing 'fraction_unspliced', 'pct_counts_MT', 'kmeans', 'pct_counts_nonCM'.
        unspliced_diff (float): Threshold for fraction_unspliced difference.
        mito_diff (float): Threshold for pct_counts_MT difference.

    Returns:
        pd.Series: Boolean Series indicating outlier cells.
    """
    # Find cluster with highest mean pct_counts_nonCM
    cluster_means = df.groupby('kmeans')['pct_counts_nonCM'].mean()
    max_noncm_cluster = cluster_means.idxmax()

    # Reference quantiles from the cluster with highest pct_counts_nonCM
    ref_cluster = df[df['kmeans'] == max_noncm_cluster]
    unspliced_threshold = ref_cluster['fraction_unspliced'].quantile(0.25) - unspliced_diff
    mito_threshold = ref_cluster['pct_counts_MT'].quantile(0.75) + mito_diff

    # Annotate outliers
    is_outlier = ~(
        (df['fraction_unspliced'] > unspliced_threshold) &
        (df['pct_counts_MT'] < mito_threshold)
    )
    return is_outlier


def parse_bam_tags(
    interval: tuple,
    bam_path: str,
    barcodes: List[str],
    CB_tag: str,
    RE_tag: str,
    EXON_tag: str,
    INTRON_tag: str,
) -> Optional[pd.DataFrame]:
    """
    Parse BAM file tags for a given genomic interval.

    Parameters:
        interval (tuple): Tuple of (contig, start, end).
        bam_path (str): Path to the BAM file.
        barcodes (List[str]): List of cell barcodes.
        CB_tag (str): Tag for cell barcode in BAM file.
        RE_tag (str): Tag for region type in BAM file.
        EXON_tag (str): Tag indicating exon region.
        INTRON_tag (str): Tag indicating intron region.

    Returns:
        Optional[pd.DataFrame]: DataFrame with counts of exon and intron reads per barcode.
    """
    bam_file = pysam.AlignmentFile(bam_path, "rb")
    tags = {'CB': [], 'RE': []}

    for read in bam_file.fetch(interval[0], interval[1], interval[2]):
        if not read.has_tag(CB_tag) or not read.has_tag(RE_tag):
            continue
        cb_tag = read.get_tag(CB_tag)
        re_tag = read.get_tag(RE_tag)
        if cb_tag in barcodes:
            tags['CB'].append(cb_tag)
            tags['RE'].append(re_tag)

    if not tags['CB']:
        return None

    tags_df = pd.DataFrame(tags)
    tags_df = tags_df.dropna()

    # Count occurrences
    count_data = (
        tags_df.groupby(['CB', 'RE'])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=[EXON_tag, INTRON_tag], fill_value=0)
    )
    return count_data


def fraction_unspliced_from_bam(
    bam_path: Optional[str] = None,
    bam_index: Optional[str] = None,
    barcodes: Optional[List[str]] = None,
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

    Parameters:
        bam_path (str, optional): Path to the BAM file.
        bam_index (str, optional): Path to the BAM index file.
        barcodes (List[str], optional): List of cell barcodes.
        regions (List[tuple], optional): List of genomic intervals to process.
        tiles (int, optional): Number of genomic regions to process in parallel.
        cores (int, optional): Number of CPU cores to use.
        CB_tag (str, optional): Tag for cell barcode in BAM file.
        RE_tag (str, optional): Tag for region type in BAM file.
        EXON_tag (str, optional): Tag indicating exon region.
        INTRON_tag (str, optional): Tag indicating intron region.

    Returns:
        Optional[pd.DataFrame]: DataFrame containing the fraction of unspliced reads per cell.
    """
    if bam_path is None or bam_index is None or barcodes is None:
        raise ValueError("Please provide bam_path, bam_index, and barcodes.")

    if not os.path.exists(bam_path):
        raise FileNotFoundError(f"The BAM file '{bam_path}' does not exist.")

    if not os.path.exists(bam_index):
        raise FileNotFoundError(f"The BAM index file '{bam_index}' does not exist.")

    # Remove any suffixes from barcodes
    barcodes = [barcode.split("-")[0] for barcode in barcodes]

    bam_file = pysam.AlignmentFile(bam_path, "rb", index_filename=bam_index)
    if cores is None:
        cores = max(1, os.cpu_count() - 1)

    if regions is None:
        # Split the genome into regions
        total_length = sum(bam_file.lengths)
        tile_size = total_length // tiles
        regions = []
        for contig, length in zip(bam_file.references, bam_file.lengths):
            for start in range(0, length, tile_size):
                end = min(start + tile_size, length)
                regions.append((contig, start, end))

    results = []
    with ProcessPoolExecutor(max_workers=cores) as executor:
        futures = {
            executor.submit(
                parse_bam_tags,
                region,
                bam_path,
                barcodes,
                CB_tag,
                RE_tag,
                EXON_tag,
                INTRON_tag,
            ): region
            for region in regions
        }

        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing BAM tiles", unit="tile"):
            result = future.result()
            if result is not None:
                results.append(result)

    if not results:
        return None

    final_df = pd.concat(results).groupby(level=0).sum()
    total_counts = final_df[INTRON_tag] + final_df[EXON_tag]
    final_df['fraction_unspliced'] = final_df[INTRON_tag] / total_counts
    final_df['fraction_unspliced'].fillna(0, inplace=True)

    return final_df[['fraction_unspliced']]