from qclus.utils import *
from qclus.gene_lists import *
import scanpy as sc
from typing import List, Dict

def run_qclus(
    counts_path: str,
    fraction_unspliced: pd.Series,
    nucl_gene_set: List[str] = nucl_30,
    celltype_gene_set_dict: Dict[str, List[str]] = celltype_gene_set_dict,
    minimum_genes: int = 500,
    maximum_genes: int = 6000,
    max_mito_perc: float = 40.0,
    clustering_features: List[str] = [
        'pct_counts_nonCM',
        'pct_counts_nuclear',
        'pct_counts_MT',
        'pct_counts_CM_cyto',
        'pct_counts_CM_nucl',
        'fraction_unspliced',
    ],
    clustering_k: int = 4,
    clusters_to_select: List[str] = ["0", "1", "2"],
    scrublet_filter: bool = True,
    scrublet_expected_rate: float = 0.06,
    scrublet_minimum_counts: int = 2,
    scrublet_minimum_cells: int = 3,
    scrublet_minimum_gene_variability_pctl: float = 85.0,
    scrublet_n_pcs: int = 30,
    scrublet_thresh: float = 0.1,
    outlier_filter: bool = True,
    outlier_unspliced_diff: float = 0.1,
    outlier_mito_diff: float = 5.0,
) -> sc.AnnData:
    """
    Run the QClus pipeline on single-cell RNA sequencing data.

    Parameters:
        counts_path (str): Path to the 10x Genomics counts .h5 file.
        fraction_unspliced (pd.Series): Series containing the fraction of unspliced reads per cell.
        nucl_gene_set (List[str], optional): List of nuclear genes for QC metrics.
        celltype_gene_set_dict (Dict[str, List[str]], optional): Dictionary of cell type-specific gene sets.
        minimum_genes (int, optional): Minimum number of genes expressed to pass initial filter.
        maximum_genes (int, optional): Maximum number of genes expressed to pass initial filter.
        max_mito_perc (float, optional): Maximum mitochondrial gene percentage to pass initial filter.
        clustering_features (List[str], optional): List of features used for clustering.
        clustering_k (int, optional): Number of clusters to use in k-means clustering.
        clusters_to_select (List[str], optional): List of cluster labels to select (pass filtering).
        scrublet_filter (bool, optional): Whether to perform doublet filtering using Scrublet.
        scrublet_expected_rate (float, optional): Expected doublet rate for Scrublet.
        scrublet_minimum_counts (int, optional): Minimum counts per cell for Scrublet.
        scrublet_minimum_cells (int, optional): Minimum cells per gene for Scrublet.
        scrublet_minimum_gene_variability_pctl (float, optional): Minimum gene variability percentile for Scrublet.
        scrublet_n_pcs (int, optional): Number of principal components for Scrublet.
        scrublet_thresh (float, optional): Threshold for Scrublet doublet calling.
        outlier_filter (bool, optional): Whether to perform outlier filtering.
        outlier_unspliced_diff (float, optional): Unspliced fraction difference threshold for outlier filtering.
        outlier_mito_diff (float, optional): Mitochondrial percentage difference threshold for outlier filtering.

    Returns:
        sc.AnnData: AnnData object containing the raw data with QClus annotations.
    """
    # Initialize AnnData object
    adata = read_count_file(counts_path)

    adata.obs.index = create_new_index(adata.obs.index)
    adata_raw = adata.copy()

    # Filter adata and adata_raw using new utility function
    adata = add_fraction_unspliced(adata, fraction_unspliced)
    adata_raw = add_fraction_unspliced(adata_raw, fraction_unspliced)

    # Calculate QC metrics
    get_qc_metrics(adata, nucl_gene_set, 'nuclear', normlog=True)

    get_qc_metrics(adata, MT_genes, 'MT')

    # Add CM-specific annotations if included in clustering features
    if 'pct_counts_CM_cyto' in clustering_features and 'pct_counts_CM_nucl' in clustering_features:
        for entry in CM_gene_set_dict:
            get_qc_metrics(adata, CM_gene_set_dict[entry], entry)

    if 'pct_counts_nonCM' in clustering_features:
        # Add cell type-specific annotations from given gene sets
        for entry in celltype_gene_set_dict:
            get_qc_metrics(adata, celltype_gene_set_dict[entry], entry)

        # Create nonCM annotations
        ct_spec_columns = ['pct_counts_' + ct for ct in celltype_gene_set_dict.keys()]

        missing_columns = [col for col in ct_spec_columns if col not in adata.obs.columns]
        if missing_columns:
            raise ValueError(f"The following required columns are missing in adata.obs: {missing_columns}")

        adata.obs["pct_counts_nonCM"] = adata.obs[ct_spec_columns].max(axis=1)

    # Synchronize observations in raw data
    adata_raw.obs = adata.obs.copy()

    # Initial filter based on gene counts and mitochondrial percentage
    adata.obs["initial_filter"] = (
        (adata.obs['n_genes_by_counts'] < minimum_genes) |
        (adata.obs['n_genes_by_counts'] > maximum_genes) |
        (adata.obs['pct_counts_MT'] > max_mito_perc)
    )
    initial_filter_list = adata.obs.index[adata.obs.initial_filter].tolist()
    adata = adata[~adata.obs.initial_filter]

    # Calculate Scrublet scores
    if scrublet_filter:
        adata.obs['score_scrublet'] = calculate_scrublet(
            adata,
            expected_rate=scrublet_expected_rate,
            minimum_counts=scrublet_minimum_counts,
            minimum_cells=scrublet_minimum_cells,
            minimum_gene_variability_pctl=scrublet_minimum_gene_variability_pctl,
            n_pcs=scrublet_n_pcs,
            thresh=scrublet_thresh,
        )

    # Normalize and logarithmize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Check if clustering features are available
    missing_features = [feat for feat in clustering_features if feat not in adata.obs.columns]
    if missing_features:
        raise ValueError(f"The following clustering features are missing in adata.obs: {missing_features}")

    # Add QClus embedding
    cluster_embedding = add_qclus_embedding(adata, clustering_features, random_state=1, n_components=2)
    adata_raw.uns["QClus_umap"] = cluster_embedding  # Store the embedding

    # Perform unsupervised clustering
    adata.obs["kmeans"] = do_kmeans(adata.obs.loc[:, clustering_features], k=clustering_k)

    # Add cluster results to adata_raw
    adata_raw.obs["kmeans"] = "initial filter"
    adata_raw.obs.loc[adata.obs.index, "kmeans"] = adata.obs.kmeans

    # Clustering filter
    adata.obs["clustering_filter"] = ~adata.obs.kmeans.isin(clusters_to_select)
    clustering_filter_list = adata.obs.index[adata.obs.clustering_filter].tolist()
    adata = adata[~adata.obs.clustering_filter]

    # Outlier filter
    outlier_filter_list = []
    if outlier_filter:
        adata.obs["outlier_filter"] = annotate_outliers(
            adata.obs[["fraction_unspliced", "pct_counts_MT", "kmeans"]],
            unspliced_diff=outlier_unspliced_diff,
            mito_diff=outlier_mito_diff,
        )
        outlier_filter_list = adata.obs.index[adata.obs.outlier_filter].tolist()
        adata = adata[~adata.obs.outlier_filter]

    # Scrublet filter
    scrublet_filter_list = []
    if scrublet_filter:
        adata.obs["scrublet_filter"] = adata.obs.score_scrublet >= scrublet_thresh
        scrublet_filter_list = adata.obs.index[adata.obs.scrublet_filter].tolist()
        adata = adata[~adata.obs.scrublet_filter]

    # Annotate raw counts with the results of QClus
    adata_raw.obs["qclus"] = "passed"
    adata_raw.obs.loc[initial_filter_list, "qclus"] = "initial filter"
    adata_raw.obs.loc[clustering_filter_list, "qclus"] = "clustering filter"
    if outlier_filter_list:
        adata_raw.obs.loc[outlier_filter_list, "qclus"] = "outlier filter"
    if scrublet_filter_list:
        adata_raw.obs.loc[scrublet_filter_list, "qclus"] = "scrublet filter"

    return adata_raw
