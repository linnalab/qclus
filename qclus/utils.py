import loompy
import pandas as pd
import scrublet as scr
from sklearn.cluster import KMeans
from sklearn.preprocessing import MinMaxScaler


def create_new_index(index):
    barcodes = [x[0:16] for x in index.tolist()]
    return barcodes

def create_fraction_unspliced_metric(loompy_path, index):
    loompy_con = loompy.connect(loompy_path)
    barcodes = [x.split(':')[1][0:16] for x in loompy_con.ca['CellID'].tolist()]
    spliced_df = pd.Series(loompy_con['spliced'][:,:].sum(0).tolist(), index = barcodes)
    unspliced_df = pd.Series(loompy_con['unspliced'][:,:].sum(0).tolist(), index = barcodes)
    ambiguous_df = pd.Series(loompy_con['ambiguous'][:,:].sum(0).tolist(), index = barcodes)
    fraction_unspliced_df = unspliced_df / (spliced_df + unspliced_df + ambiguous_df)
    barcodes_data = [x for x in index] 
    loompy_con.close()
    return fraction_unspliced_df[barcodes_data].values


def calculate_scrublet(adata, 
                       expected_rate=0.06, 
                       minimum_counts=2, 
                       minimum_cells=3, 
                       minimum_gene_variability_pctl=85, 
                       n_pcs=30, 
                       thresh=0.10):
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=expected_rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=minimum_counts, 
                                                                        min_cells=minimum_cells, 
                                                                        min_gene_variability_pctl=minimum_gene_variability_pctl, 
                                                                        n_prin_comps=n_pcs)
    scrub.call_doublets(threshold=thresh)
    return scrub.doublet_scores_obs_

def do_kmeans(X_full, k):
    X_minmax = MinMaxScaler().fit_transform(X_full)
    X_minmax = pd.DataFrame(X_minmax, index = X_full.index, columns = X_full.columns)
    kmeans_custom = KMeans(n_clusters=k, random_state=0).fit(X_minmax) 
    X_full['kmeans'] = kmeans_custom.labels_.astype(str)

    #Sort kmeans clusters 0,1,2,...,k by decreasing fraction_unspliced mean.
    dflist = []
    clusters = []
    for i in range(k):
        df = X_full[X_full.kmeans==str(i)]
        dflist.append(df)
    sorteddflist = sorted(dflist,key=lambda x:x["fraction_unspliced"]. mean(axis=0))
    sorteddflist.reverse()
    for dframe in sorteddflist:
        clusters.append(str(int(dframe.kmeans.median())))
    return [str(clusters.index(x)) for x in X_full.kmeans]

def annotate_outliers(df, unspliced_diff, mito_diff):
    #find kmeans cluster with highest mean pct_counts_nonCM
    max_non_CM = 0
    max_value = 0
    for clus in set(df.kmeans):
        if (df.loc[df.kmeans == clus, "pct_counts_nonCM"].mean() > max_value):
            max_non_CM = clus
            max_value = df.loc[df.kmeans == clus, "pct_counts_nonCM"].mean()
            
    #tag cells that are in the max_non_cm cluster as 1, others 0
    df['max_non_CM'] = ["1" if x in max_non_CM else "0" for x in df.kmeans]
    
    #annotate cells as outlier ("1") or not outlier ("0")
    return [False if 
            ((fraction_unspliced > (df[df.max_non_CM == "1"].fraction_unspliced.quantile(.25) - unspliced_diff)) 
             & (pct_counts_MT < (df[df.max_non_CM == "1"].pct_counts_MT.quantile(.75) + mito_diff))) 
            else True 
            for fraction_unspliced,pct_counts_MT 
            in zip(df.fraction_unspliced, df.pct_counts_MT)]
