from qclus.utils import *
from qclus.gene_lists import *
import scanpy as sc
import warnings
from typing import List, Dict, Union

def run_qclus(
    counts_path: str,
    fraction_unspliced: Union[pd.Series, pd.DataFrame],
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
    kmeans_n_init: int = 1,
    scrublet_approx_neighbors: bool = False,
) -> sc.AnnData:
    """
    Run the QClus pipeline on single-cell RNA sequencing data.

    Parameters:
        counts_path (str): Path to the 10x Genomics counts .h5 file.
        fraction_unspliced (pd.Series or pd.DataFrame): Fraction of unspliced reads per cell, indexed by
            cell barcode. A DataFrame must have a 'fraction_unspliced' column or exactly one column.
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
        kmeans_n_init (int, optional): Number of k-means restarts; the best one is kept. The default of 1
            is what scikit-learn 1.4 and later run when it is not set, as in QClus 0.2.0 and earlier.
        scrublet_approx_neighbors (bool, optional): Whether Scrublet uses approximate nearest neighbours
            (annoy) instead of exact ones. QClus 0.2.0 and earlier used approximate neighbours, which
            give wrong results where annoy is broken.

    Returns:
        sc.AnnData: AnnData object containing the raw data with QClus annotations. It holds the barcodes
            that have splicing information, under their 16-character names; obs['original_barcode']
            keeps the names from the counts file.
    """
    # Check the settings before any data is read
    if 'fraction_unspliced' not in clustering_features:
        raise ValueError(
            "'fraction_unspliced' must be one of the clustering_features, because the clusters are "
            "ordered by it."
        )

    clusters_to_select = [str(cluster) for cluster in clusters_to_select]
    valid_clusters = [str(i) for i in range(clustering_k)]
    if not clusters_to_select or not set(clusters_to_select) <= set(valid_clusters):
        raise ValueError(
            f"clusters_to_select must be a non-empty subset of {valid_clusters} for "
            f"clustering_k={clustering_k}, got {clusters_to_select}."
        )

    fraction_unspliced = prepare_fraction_unspliced(fraction_unspliced)

    # Initialize AnnData object

    adata = read_count_file(counts_path)

    adata.obs["original_barcode"] = adata.obs.index.astype(str)
    adata.obs.index = create_new_index(adata.obs.index)
    check_unique_barcodes(adata.obs.index, "the counts file")
    adata_raw = adata.copy()

    # Filter adata and adata_raw using new utility function
    adata = add_fraction_unspliced(adata, fraction_unspliced)
    adata_raw = add_fraction_unspliced(adata_raw, fraction_unspliced)

    # Calculate QC metrics
    get_qc_metrics(adata, nucl_gene_set, 'nuclear', normlog=True)

    get_qc_metrics(adata, MT_genes, 'MT')

    # Add CM-specific annotations if included in clustering features
    for entry in CM_gene_set_dict:
        if 'pct_counts_' + entry in clustering_features:
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

    # Check that enough barcodes are left for the steps that follow
    minimum_barcodes = [(clustering_k, f"k-means with clustering_k={clustering_k}"), (4, "the UMAP embedding")]
    if scrublet_filter:
        minimum_barcodes.append((scrublet_n_pcs + 1, f"the doublet filter with scrublet_n_pcs={scrublet_n_pcs}"))
    for minimum, step in minimum_barcodes:
        if adata.n_obs < minimum:
            raise ValueError(
                f"Only {adata.n_obs} of {adata_raw.n_obs} barcodes remain after the initial filter, "
                f"but {step} needs at least {minimum}."
            )

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
            approx_neighbors=scrublet_approx_neighbors,
        )
        adata_raw.obs["score_scrublet"] = np.nan
        adata_raw.obs.loc[adata.obs.index, "score_scrublet"] = adata.obs.score_scrublet

    # Normalize and logarithmize
    # sc.pp.normalize_total(adata, target_sum=1e4)
    # sc.pp.log1p(adata)

    # Check if clustering features are available
    missing_features = [feat for feat in clustering_features if feat not in adata.obs.columns]
    if missing_features:
        raise ValueError(f"The following clustering features are missing in adata.obs: {missing_features}")

    # Add QClus embedding
    cluster_embedding = add_qclus_embedding(adata, clustering_features, random_state=1, n_components=2)
    adata_raw.uns["QClus_umap"] = cluster_embedding  # Store the embedding

    # Perform unsupervised clustering
    adata.obs["kmeans"] = do_kmeans(adata.obs.loc[:, clustering_features], k=clustering_k, n_init=kmeans_n_init)

    # Add cluster results to adata_raw
    adata_raw.obs["kmeans"] = "initial filter"
    adata_raw.obs.loc[adata.obs.index, "kmeans"] = adata.obs.kmeans

    # Clustering filter
    adata.obs["clustering_filter"] = ~adata.obs.kmeans.isin(clusters_to_select)
    clustering_filter_list = adata.obs.index[adata.obs.clustering_filter].tolist()
    adata = adata[~adata.obs.clustering_filter]

    # K-means can fill fewer than clustering_k clusters when feature vectors repeat
    if adata.n_obs == 0:
        filled_clusters = adata_raw.obs.loc[adata_raw.obs.kmeans != "initial filter", "kmeans"].nunique()
        raise ValueError(
            f"No barcode is in the selected clusters {clusters_to_select}: k-means filled only "
            f"{filled_clusters} of {clustering_k} clusters."
        )

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

    if not (adata_raw.obs["qclus"] == "passed").any():
        warnings.warn("No barcode passed QClus.")

    return adata_raw


def quickstart_qclus(counts_path: str,
                     fraction_unspliced: Union[pd.Series, pd.DataFrame],
                     tissue: str = 'heart'):
    """
    A function that performs quick clustering based on the specified tissue type.
    The function processes single-cell RNA sequencing data to identify specific
    clusters of interest, accommodating different clustering workflows by tissue type.
    If the tissue is "heart", it runs a standard protocol, while for other tissues,
    custom parameters can be applied. The function raises a ValueError if an invalid
    tissue type is specified.

    :param counts_path: Path to the count matrix used as input data.
    :type counts_path: str
    :param fraction_unspliced: A Series or DataFrame with the fraction of unspliced RNA for each cell.
    :type fraction_unspliced: pd.Series or pd.DataFrame
    :param tissue: The type of tissue to analyze, which determines the clustering workflow. Default is 'heart'.
    :type tissue: str

    :return: An AnnData object that contains the results of the clustering.
    :rtype: AnnData
    """
    if tissue == 'heart':
        # Run the standard "heart" workflow
        adata = run_qclus(counts_path, fraction_unspliced)

    elif tissue == 'other':
        # Run the "other" tissue workflow with custom parameters
        adata = run_qclus(counts_path, fraction_unspliced,
                          clustering_features=['pct_counts_nuclear',
                                               'pct_counts_MT',
                                               'fraction_unspliced'],
                          clustering_k=3,
                          clusters_to_select=["0", "1"],
                          )

    else:
        # Raise an error for any invalid tissue types
        raise ValueError(f"Invalid tissue: {tissue}. Please choose 'heart' or 'other'.")

    return adata
