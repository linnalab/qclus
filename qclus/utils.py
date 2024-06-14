import loompy
import pandas as pd
import scrublet as scr
from sklearn.cluster import KMeans
from sklearn.preprocessing import MinMaxScaler
import pysam
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import os
from umap import UMAP


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
                                                                        n_prin_comps=n_pcs, 
                                                                        verbose=False)
    scrub.call_doublets(threshold=thresh, verbose=False)
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
        dframe['kmeans'] = pd.to_numeric(dframe['kmeans'], errors='coerce')
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

def add_embedding(adata, features, embedding_key='X_umap', random_state=1, n_components=2):
    # Extract the features from the AnnData object
    X_full = adata.obs.loc[:, features]
    # Scale the features
    X_minmax = MinMaxScaler().fit_transform(X_full)
    X_minmax = pd.DataFrame(X_minmax, index=X_full.index, columns=X_full.columns)
    # Compute the embedding
    umap = UMAP(random_state=random_state, n_components=n_components).fit_transform(X_minmax)
    # Add the embedding to the AnnData object
    adata.obsm[embedding_key] = umap
    return adata

def parse_bam_tags(interval, bam, bc, CB_tag, RE_tag, EXON_tag, INTRON_tag):
    bam_file = pysam.AlignmentFile(bam, "rb")
    tags = {
        'CB': [],
        'RE': []
    }
    
    for read in bam_file.fetch(interval[0], interval[1], interval[2]):
        cb_tag = read.get_tag(CB_tag) if read.has_tag(CB_tag) else None
        re_tag = read.get_tag(RE_tag) if read.has_tag(RE_tag) else None
        if cb_tag and cb_tag in bc:
            tags['CB'].append(cb_tag)
            tags['RE'].append(re_tag)

    if not tags['CB']:  # No reads found
        return None

    tags_df = pd.DataFrame(tags)
    tags_df = tags_df.dropna()  # Remove rows with NaN values

    # Count occurrences
    count_data = tags_df.groupby(['CB', 'RE']).size().unstack(fill_value=0).reindex(columns=[EXON_tag, INTRON_tag], fill_value=0)
    return count_data

def nuclear_fraction_tags(outs=None, bam=None, bam_index=None, barcodes=None, regions=None, tiles=100, cores=None, cell_barcode_tag="CB", region_type_tag="RE", exon_tag="E", intron_tag="N", verbose=True):

    if outs:
        bam = f"{outs}/possorted_genome_bam.bam"
        bam_index = f"{outs}/possorted_genome_bam.bam.bai"
        barcodes_path = f"{outs}/filtered_feature_bc_matrix/barcodes.tsv.gz"
        barcodes = pd.read_csv(barcodes_path, header=None, compression='gzip')[0].tolist() if barcodes is None else barcodes

    if not all([bam, bam_index, barcodes]):
        raise ValueError("Please provide bam, bam_index, and barcodes")

    bam_file = pysam.AlignmentFile(bam, "rb", index_filename=bam_index)
    if cores is None:
        cores = max(1, os.cpu_count() - 1)  # Use all available cores minus one

    if regions is None:
        # Split the genome into roughly equal-sized regions based on the number of tiles specified
        total_length = sum(bam_file.lengths)
        tile_size = total_length // tiles
        regions = []
        for contig, length in zip(bam_file.references, bam_file.lengths):
            for start in range(0, length, tile_size):
                end = min(start + tile_size, length)
                regions.append((contig, start, end))

    results = []
    with ProcessPoolExecutor(max_workers=cores) as executor:
        futures = {executor.submit(parse_bam_tags, region, bam, barcodes, cell_barcode_tag, region_type_tag, exon_tag, intron_tag): region for region in regions}
        
        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing BAM tiles", unit="tile"):
            result = future.result()
            if result is not None:
                results.append(result)
    
    if not results:
        return None
    
    final_df = pd.concat(results).groupby(level=0).sum()
    final_df['unspliced_fraction'] = final_df[intron_tag] / (final_df[intron_tag] + final_df[exon_tag])
    final_df['unspliced_fraction'].fillna(0, inplace=True)
        
    return final_df['unspliced_fraction']