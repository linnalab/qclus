import scanpy as sc

from qclus.utils import *
from qclus.gene_lists import *

def run_qclus(counts_path, loompy_path, 
                    gene_set_dict=celltype_gene_set_dict, 
                    nucl_gene_set=nucl_genes_50, 
                    minimum_genes=500, 
                    maximum_genes=6000, 
                    max_mito_perc=40, 
                    clustering_features=['pct_counts_nonCM', 
                                    'pct_counts_nucl_30', 
                                    'pct_counts_MT', 
                                    'pct_counts_CM_cyto', 
                                    'pct_counts_CM_nucl', 
                                    'fraction_unspliced'], 
                    clustering_k=4, 
                    clusters_to_select=["0", "1", "2"], 
                    scrublet=True,
                    scrublet_expected_rate=0.06, 
                    scrublet_minimum_counts=2, 
                    scrublet_minimum_cells=3, 
                    scrublet_minimum_gene_variability_pctl=85, 
                    scrublet_n_pcs=30, 
                    scrublet_thresh=0.1, 
                    outlier_filter=True, 
                    outlier_unspliced_diff=0.1, 
                    outlier_mito_diff=5):
        
    adata = sc.read_10x_h5(f"{counts_path}")
    adata.var_names_make_unique()

    adata.obs.index = create_new_index(adata.obs.index)
    adata.obs["fraction_unspliced"] = create_fraction_unspliced_metric(loompy_path, adata.obs.index)
    
    for entry in gene_set_dict:
        adata.var[entry] = [True if x in gene_set_dict[entry] else False for x in adata.var.index]
        sc.pp.calculate_qc_metrics(adata, qc_vars=[entry], percent_top=None, log1p=False, inplace=True)
        sc.tl.score_genes(adata, gene_list = gene_set_dict[entry], score_name = f"score_{entry}")
        
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


    adata.obs["initial_filter"] = ["keep" if maximum_genes >= x >= minimum_genes and y <= max_mito_perc else "remove" for x,y in zip(adata.obs.n_genes_by_counts, adata.obs.pct_counts_MT)]
    adata = adata[adata.obs["initial_filter"]=="keep"]
    
    if scrublet:
        adata.obs['score_scrublet'] = calculate_scrublet(adata, 
                                                         expected_rate=scrublet_expected_rate, 
                                                         minimum_counts=scrublet_minimum_counts, 
                                                         minimum_cells=scrublet_minimum_cells, 
                                                         minimum_gene_variability_pctl=scrublet_minimum_gene_variability_pctl, 
                                                         n_pcs=scrublet_n_pcs, 
                                                         thresh=scrublet_thresh)
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    adata.var["nucl_30"] = [True if x in nucl_gene_set[:30] else False for x in adata.var.index]
    sc.pp.calculate_qc_metrics(adata, qc_vars=["nucl_30"], percent_top=None, log1p=False, inplace=True)
    sc.tl.score_genes(adata, gene_list = nucl_gene_set[:30], score_name = "score_nucl_30")
        
    adata.obs["kmeans"] = do_kmeans(adata.obs.loc[:,clustering_features], k=clustering_k)

    #add filter annotations
    adata.obs["clustering"] = ["keep" if x in clusters_to_select else "remove" for x in adata.obs.kmeans]
    if outlier_filter:
        adata.obs["outlier"] = annotate_outliers(adata.obs[["fraction_unspliced", "pct_counts_MT", "kmeans", "pct_counts_nonCM"]], unspliced_diff=outlier_unspliced_diff, mito_diff=outlier_mito_diff)
    if scrublet:
        adata.obs["scrublet"] = ["keep" if x < scrublet_thresh else "remove" for x in adata.obs.score_scrublet]

    return adata