import numpy as np 
import pandas as pd
import scrublet as scr
import loompy
import scanpy as sc
import anndata
from sklearn.cluster import KMeans
from sklearn.preprocessing import MinMaxScaler
sc.settings.verbosity = 0   # verbosity: errors (0), warnings (1), info (2), hints (3)

from qclus.utils import *
from qclus.gene_lists import nucl_genes_50, celltype_gene_set_dict


def run_qclus(counts_path, loompy_path, 
                    gene_set_dict=celltype_gene_set_dict, 
                    nucl_gene_set=nucl_genes_50, 
                    minimum_genes=500, 
                    maximum_genes=6000, 
                    max_mito_perc=40, 
                    clustering=True, 
                    clustering_features=['pct_counts_nonCM', 
                                    'pct_counts_nucl_30', 
                                    'pct_counts_MT', 
                                    'pct_counts_CM_cyto', 
                                    'pct_counts_CM_nucl', 
                                    'fraction_unspliced'], 
                    clustering_k=4, 
                    clusters_to_select=["0", "1", "2"], 
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
    adata.obs["pct_counts_nonCM"] = adata.obs[['pct_counts_VEC', 'pct_counts_PER',  'pct_counts_SMC',  'pct_counts_AD',  'pct_counts_SC',  'pct_counts_N',  'pct_counts_EEC',  'pct_counts_FB',  'pct_counts_L',  'pct_counts_MESO',  'pct_counts_MP']].max(1)
    adata.obs["score_nonCM"] = adata.obs[['score_VEC', 'score_PER',  'score_SMC',  'score_AD',  'score_SC',  'score_N',  'score_EEC',  'score_FB',  'score_L',  'score_MESO',  'score_MP']].max(1)

    adata = adata[adata.obs.n_genes_by_counts >= minimum_genes]
    adata = adata[adata.obs.n_genes_by_counts <= maximum_genes]
    adata = adata[adata.obs.pct_counts_MT <= max_mito_perc]
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    adata.var["nucl_30"] = [True if x in nucl_genes_50[:30] else False for x in adata.var.index]
    sc.pp.calculate_qc_metrics(adata, qc_vars=["nucl_30"], percent_top=None, log1p=False, inplace=True)
    sc.tl.score_genes(adata, gene_list = nucl_genes_50[:30], score_name = "score_nucl_30")
        
    if clustering:
        adata.obs["kmeans"] = do_kmeans(adata.obs.loc[:,clustering_features], k=clustering_k)
    if outlier_filter:
        adata.obs["outlier"] = annotate_outliers(adata.obs[["fraction_unspliced", "pct_counts_MT", "kmeans", "pct_counts_nonCM"]], unspliced_diff=outlier_unspliced_diff, mito_diff=outlier_mito_diff)
            
    if clustering:
        adata = adata[adata.obs.kmeans.isin(clusters_to_select)]
    if outlier_filter:
        adata = adata[adata.obs.outlier=="0"]

    return adata.obs.index.to_list()