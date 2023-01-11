import numpy as np 
import pandas as pd
import scrublet as scr
import loompy
import scanpy as sc
import anndata
from sklearn.cluster import KMeans
from sklearn.preprocessing import MinMaxScaler
sc.settings.verbosity = 0   # verbosity: errors (0), warnings (1), info (2), hints (3)

nucl_genes_50 = ['MALAT1', 'NEAT1', 'FTX', 'FOXP1', 'RBMS3', 'ZBTB20', 'LRMDA', 'PBX1', 'ITPR2', 'AUTS2', 'TTC28', 'BNC2', 'EXOC4', 'RORA', 'PRKG1', 'ARID1B', 'PARD3B', 'GPHN', 'N4BP2L2', 'PKHD1L1', 'EXOC6B', 'FBXL7', 'MED13L', 'TBC1D5', 'IMMP2L', 'SYNE1', 'RERE', 'MBD5', 'EXT1', 'WWOX', 'EPB41L4A', 'PTK2', 'ST6GAL1', 'CHD9', 'PTPRG', 'JMJD1C', 'WSB1', 'SBF2', 'STAG1', 'GMDS', 'ADAMTS9-AS2', 'PDE7B', 'RALGAPA2', 'PRKN', 'ZFAND3', 'PLXDC2', 'NAALADL2', 'ASAP1', 'NFIA', 'MIR99AHG']

### DEFINE GENE SETS
celltype_gene_set_dict = {
"MT" : ['MT-ND1', 'MT-ND2', 'MT-CO1', 'MT-CO2', 'MT-ATP8', 'MT-ATP6', 'MT-CO3', 'MT-ND3', 'MT-ND4L', 'MT-ND4', 'MT-ND5', 'MT-ND6', 'MT-CYB'],
"MT_nucl" : ["POLG", "POLG2", "TWNK", "SSBP1", "PRIMPOL", "DNA2", "MGME1", "RNASEH1", "POLRMT", "TFAM", "TEFM", "TFB2M", "TK2", "DGUOK", "RRM2B", "TYMP", "SLC25A4", "OPA1", "MFN1", "MFN2", "DNM1L", "MFF", "FIS1", "MPV17", "SPG7"], 
"CM_cyto" : ["TTN", "RYR2", "PAM", "TNNT2", "RABGAP1L", "PDLIM5", "MYL7", "MYH6"],
"CM_nucl" : ["RBM20", "TECRL", "MLIP", "CHRM2", "TRDN", "PALLD", "SGCD", "CMYA5", "MYOM2", "TBX5", "ESRRG", "LINC02248", "KCNJ3", "TACC2", "CORIN", "DPY19L2", "WNK2", "MITF", "OBSCN", "FHOD3", "MYLK3", "DAPK2", "NEXN"],
"VEC" : ["VWF", "ERG", "ANO2", "PTPRB", "EGFL7", "PREX2", "ADGRL4", "FLT1", "CYYR1", "GRB10", "PPP1R16B", "DOCK9", "SHANK3", "PECAM1", "PLEKHG1", "EMCN"],
"PER" : ["RGS5", "DACH1", "GUCY1A1", "ABCC9", "BGN", "NOTCH3", "PDGFRB", "FRMD3", "RNF152", "CCDC102B", "NGF"],
"SMC" : ["MYH11", "KCNAB1", "NTRK3", "CHRM3", "ACTA2", "RGS6", "DGKG", "ITGA8", "TBX2", "LMOD1", "SDK1", "GPC6", "ANTXR1", "FLNA", "CLMN", "ATP10A", "MCAM", "TAGLN", "CCDC3"],
"AD" : ["PLIN4", "PLIN1", "PDE3B", "GPAM", "PTPRS", "PPARG", "MLXIPL", "MGST1", "AQP7", "SLC19A3", "FABP4", "TPRG1", "DIRC3", "LPL", "PNPLA2", "LIPE", "ADH1B", "ADIPOQ", "PRKAR2B", "CIDEA", "LINC00278", "PFKFB3", "LINC02237", "LIPE-AS1", "SVEP1"],
"SC" : ["XKR4", "AC016766.1", "SLC35F1", "ZNF536", "NCAM2", "GPM6B", "KIRREL3", "SORCS1", "ST6GALNAC5", "PRKCA", "GINS3", "PMP22", "ALDH1A1", "IL1RAPL2", "DOCK5", "NKAIN3", "COL28A1", "RALGPS2", "PKN2-AS1", "KLHL29", "PTPRZ1"],
"N" : ["CSMD1", "SYT1", "KCNIP4", "CNTNAP2", "DLGAP1", "PTPRD", "LRRTM4", "ATRNL1", "LRP1B", "CTNND2", "KCNQ5", "NRG3", "SNTG1", "GRIA2", "RIMS2", "CSMD3", "XIST", "KAZN", "DPP10", "HS6ST3", "OPCML"],
"EEC" : ["PCDH7", "PCDH15", "LINC02147", "LINC02388", "MYRIP", "GMDS", "ADAMTSL1", "LEPR", "CALCRL", "CGNL1", "HMCN1", "NPR3", "POSTN"],
"FB" : ["DCN", "ABCA8", "ABCA6", "ABCA10", "FBLN1", "COL15A1", "FBN1", "C7"],
"L" : ["SKAP1", "RIPOR2", "CD247", "IKZF1", "BCL11B", "SLFN12L", "ITGAL", "SAMD3", "CARD11", "CDC42SE2", "CCND3"],
"MESO" : ["C3", "SULF1", "AP000561.1", "PRG4", "GPM6A", "CDON", "DPP6", "CCDC80", "EZR", "FOS", "BNC1", "AC245041.2", "PRKD1", "CYSTM1", "TLL1", "WT1"],
"MP" : ["TBXAS1", "SLC9A9", "MRC1", "MS4A6A", "RBM47", "DOCK2", "MCTP1", "SYK", "MSR1", "ATP8B4", "F13A1", "CD74", "MS4A4E", "ADAP2"]
}

nucl_gene_set_dict = {
"nucl_50" : nucl_genes_50[:50],
"nucl_40" : nucl_genes_50[:40],
"nucl_30" : nucl_genes_50[:30],
"nucl_20" : nucl_genes_50[:20],
"nucl_10" : nucl_genes_50[:10],
"nucl_2" : nucl_genes_50[:2],
"nucl_1" : nucl_genes_50[:1]
}

def qclus(outs_path, loompy_path, 
                    gene_set_dict=celltype_gene_set_dict, 
                    nucl_gene_set_dict=nucl_gene_set_dict, 
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
    
    
    adata = sc.read_10x_h5(f"{outs_path}/filtered_feature_bc_matrix.h5")
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
    
    for entry in nucl_gene_set_dict:
        adata.var[entry] = [True if x in nucl_gene_set_dict[entry] else False for x in adata.var.index]
        sc.pp.calculate_qc_metrics(adata, qc_vars=[entry], percent_top=None, log1p=False, inplace=True)
        sc.tl.score_genes(adata, gene_list = nucl_gene_set_dict[entry], score_name = f"score_{entry}")
        

    if clustering:
        adata.obs["kmeans"] = do_kmeans(adata.obs.loc[:,clustering_features], k=clustering_k)
    if outlier_filter:
        adata.obs["outlier"] = annotate_outliers(adata.obs[["fraction_unspliced", "pct_counts_MT", "kmeans", "pct_counts_nonCM"]], unspliced_diff=outlier_unspliced_diff, mito_diff=outlier_mito_diff)
            
    if clustering:
        adata = adata[adata.obs.kmeans.isin(clusters_to_select)]
    if outlier_filter:
        adata = adata[adata.obs.outlier=="0"]

    return adata.obs.index.to_list()

cells_to_keep = qclus("/Users/johannesojanen/dropbox/work/mit", "/Users/johannesojanen/dropbox/work/mit/counts_counts_CAD_739.loom")
cells_to_keep