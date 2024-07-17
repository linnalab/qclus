from qclus.utils import *
from qclus.gene_lists import *
import scanpy as sc

def run_qclus(counts_path, fraction_unspliced,
                    nucl_gene_set=nucl_30,
                    celltype_gene_set_dict=celltype_gene_set_dict,
                    minimum_genes=500, 
                    maximum_genes=6000, 
                    max_mito_perc=40, 
                    clustering_features=['pct_counts_nonCM', 
                                    'pct_counts_nuclear',
                                    'pct_counts_MT', 
                                    'pct_counts_CM_cyto', 
                                    'pct_counts_CM_nucl', 
                                    'fraction_unspliced'], 
                    clustering_k=4, 
                    clusters_to_select=["0", "1", "2"], 
                    scrublet_filter=True,
                    scrublet_expected_rate=0.06, 
                    scrublet_minimum_counts=2, 
                    scrublet_minimum_cells=3, 
                    scrublet_minimum_gene_variability_pctl=85, 
                    scrublet_n_pcs=30, 
                    scrublet_thresh=0.1, 
                    outlier_filter=True, 
                    outlier_unspliced_diff=0.1, 
                    outlier_mito_diff=5):

    sc.settings.verbosity = 0

    # Initialize AnnData object
    adata = sc.read_10x_h5(f"{counts_path}")
    adata.var_names_make_unique()
    adata.obs.index = create_new_index(adata.obs.index)
    adata_raw = adata.copy()

    # Filter adata and adata_raw to keep only cells that have corresponding splicing info
    common_barcodes = adata.obs.index.intersection(fraction_unspliced.index)
    if len(common_barcodes) < len(adata.obs.index):
        print(f"Removing {len(adata.obs.index) - len(common_barcodes)} barcodes without splicing information.")

    adata = adata[common_barcodes]
    adata_raw = adata_raw[common_barcodes]

    # Add fraction_unspliced annotation from the subset .loom file
    adata.obs["fraction_unspliced"] = fraction_unspliced.loc[common_barcodes]

    #add cell type specific annotations from given gene sets
    for entry in celltype_gene_set_dict:
        adata.var[entry] = [True if x in celltype_gene_set_dict[entry] else False for x in adata.var.index]
        sc.pp.calculate_qc_metrics(adata, qc_vars=[entry], percent_top=None, log1p=False, inplace=True)
        sc.tl.score_genes(adata, gene_list = celltype_gene_set_dict[entry], score_name = f"score_{entry}")

    #initial filter
    adata.obs["initial_filter"] = [False if maximum_genes >= x >= minimum_genes and y <= max_mito_perc else True for x,y in zip(adata.obs.n_genes_by_counts, adata.obs.pct_counts_MT)]
    initial_filter_list = adata[adata.obs.initial_filter==True].obs.index.to_list()
    adata = adata[adata.obs.initial_filter==False]
        
    #create nonCM annotations
    adata.obs["pct_counts_nonCM"] = adata.obs[['pct_counts_VEC', 
                                                'pct_counts_PER',  
                                                'pct_counts_SMC',  
                                                'pct_counts_AD',  
                                                'pct_counts_SC',  
                                                'pct_counts_N',  
                                                'pct_counts_EEC',  
                                                'pct_counts_FB',  
                                                'pct_counts_L',  
                                                'pct_counts_MESO',  
                                                'pct_counts_MP']].max(1)
    adata.obs["score_nonCM"] = adata.obs[['score_VEC', 
                                                'score_PER',  
                                                'score_SMC',  
                                                'score_AD',  
                                                'score_SC',  
                                                'score_N',  
                                                'score_EEC',  
                                                'score_FB',  
                                                'score_L',  
                                                'score_MESO',  
                                                'score_MP']].max(1)
    
    #Annotate raw counts with some basic quality metrics

    #calculate scrublet score for each cell
    if scrublet_filter:
        adata.obs['score_scrublet'] = calculate_scrublet(adata, 
                                                         expected_rate=scrublet_expected_rate, 
                                                         minimum_counts=scrublet_minimum_counts, 
                                                         minimum_cells=scrublet_minimum_cells, 
                                                         minimum_gene_variability_pctl=scrublet_minimum_gene_variability_pctl, 
                                                         n_pcs=scrublet_n_pcs, 
                                                         thresh=scrublet_thresh)
    
    #normalizeand logarithmize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    #calculate nuclear score for each cell
    adata.var["nuclear"] = [True if x in nucl_gene_set else False for x in adata.var.index]
    sc.pp.calculate_qc_metrics(adata, qc_vars=["nuclear"], percent_top=None, log1p=False, inplace=True)
    sc.tl.score_genes(adata, gene_list=nucl_gene_set, score_name=f"score_nuclear")

    adata_raw.obs = adata.obs

    cluster_embedding = add_qclus_embedding(adata, clustering_features, random_state=1, n_components=2)
    adata_raw.uns["QClus_umap"] = cluster_embedding  # Store the embedding in the unstructured data

    #perform unsupervised clustering with the created cell statistics
    adata.obs["kmeans"] = do_kmeans(adata.obs.loc[:,clustering_features], k=clustering_k)

    # add cluster results to adata_raw with the same index and missing values set to "initial filter"
    adata_raw.obs["kmeans"] = "initial filter"
    adata_raw.obs.loc[adata.obs.index, "kmeans"] = adata.obs.kmeans

    adata.obs["clustering_filter"] = [False if x in clusters_to_select else True for x in adata.obs.kmeans]

    #add filter annotations
    #clustering filter
    clustering_filter_list = adata[adata.obs.clustering_filter==True].obs.index.to_list()
    adata = adata[adata.obs.clustering_filter==False]
    #outlier filter
    if outlier_filter:
        adata.obs["outlier_filter"] = annotate_outliers(adata.obs[["fraction_unspliced", "pct_counts_MT", "kmeans", "pct_counts_nonCM"]], unspliced_diff=outlier_unspliced_diff, mito_diff=outlier_mito_diff)
        outlier_filter_list = adata[adata.obs.outlier_filter==True].obs.index.to_list()
        adata = adata[adata.obs.outlier_filter==False]
    #doublet filter
    if scrublet_filter:
        adata.obs["scrublet_filter"] = [False if x < scrublet_thresh else True for x in adata.obs.score_scrublet]
        scrublet_filter_list = adata[adata.obs.scrublet_filter==True].obs.index.to_list()
        adata = adata[adata.obs.scrublet_filter==False]

    #annotate raw counts with the results of QClus
    adata_raw.obs["qclus"] = "passed"

    adata_raw.obs.loc[initial_filter_list, "qclus"] = "initial filter"
    adata_raw.obs.loc[clustering_filter_list, "qclus"] = "clustering filter"
    adata_raw.obs.loc[outlier_filter_list, "qclus"] = "outlier filter"
    adata_raw.obs.loc[scrublet_filter_list, "qclus"] = "scrublet filter"

    return adata_raw
