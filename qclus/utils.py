import loompy
import pandas as pd
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
    return ["0" if 
            ((fraction_unspliced > (df[df.max_non_CM == "1"].fraction_unspliced.quantile(.25) - unspliced_diff)) 
             & (pct_counts_MT < (df[df.max_non_CM == "1"].pct_counts_MT.quantile(.75) + mito_diff))) 
            else "1" 
            for fraction_unspliced,pct_counts_MT 
            in zip(df.fraction_unspliced, df.pct_counts_MT)]