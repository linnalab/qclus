"""Synthetic snRNA-seq data with planted populations, so that no participant data is needed for tests."""
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from qclus.gene_lists import CM_gene_set_dict, MT_genes, celltype_gene_set_dict, nucl_30

N_GENES = 2000
N_PROGRAM_GENES = 200

# Number of barcodes per planted population
DEFAULT_SIZES = {
    "CM": 150,
    "VEC": 150,
    "FB": 150,
    "debris": 60,
    "empty": 40,
    "high_mito": 30,
    "doublet": 30,
}
NUCLEI = ("CM", "VEC", "FB")


def make_barcodes(n: int, offset: int = 0) -> list:
    """Distinct 16-letter barcodes."""
    return ["".join("ACGT"[((i + offset) >> (2 * j)) & 3] for j in range(16)) for i in range(n)]


def gene_names() -> list:
    """Every gene in the built-in gene sets, padded with background genes."""
    listed = list(nucl_30) + list(MT_genes)
    for gene_set in list(CM_gene_set_dict.values()) + list(celltype_gene_set_dict.values()):
        listed += gene_set
    genes = list(dict.fromkeys(listed))
    return genes + [f"BG{i:04d}" for i in range(N_GENES - len(genes))]


def _rates(genes: pd.Index, population: str) -> np.ndarray:
    """Expected count of each gene in one droplet of the given population."""
    rates = np.full(len(genes), 0.5)

    def set_rate(names, value):
        rates[genes.isin(names)] = value

    # Marker genes are close to silent unless the population expresses them
    for gene_set in list(celltype_gene_set_dict.values()) + list(CM_gene_set_dict.values()):
        set_rate(gene_set, 0.05)
    set_rate(nucl_30, 6.0)
    set_rate(MT_genes, 1.5)

    if population in NUCLEI:
        # Each cell type has its own block of background genes, which gives Scrublet structure to find
        background = genes[genes.str.startswith("BG")]
        start = NUCLEI.index(population) * N_PROGRAM_GENES
        set_rate(background[start:start + N_PROGRAM_GENES], 2.5)
    if population == "CM":
        set_rate(CM_gene_set_dict["CM_nucl"], 3.0)
        set_rate(CM_gene_set_dict["CM_cyto"], 1.5)
    elif population == "VEC":
        set_rate(celltype_gene_set_dict["VEC"], 3.0)
    elif population == "FB":
        set_rate(celltype_gene_set_dict["FB"], 6.0)
    elif population == "debris":
        # Cytoplasmic debris: mitochondrial and cytoplasmic transcripts, few nuclear ones
        set_rate(nucl_30, 0.3)
        set_rate(MT_genes, 25.0)
        set_rate(CM_gene_set_dict["CM_cyto"], 15.0)
    elif population == "high_mito":
        set_rate(MT_genes, 120.0)
    elif population == "empty":
        rates *= 0.04
    return rates


def make_dataset(sizes: dict = None, seed: int = 0, suffix: str = "-1"):
    """
    Build a count matrix and matching fraction_unspliced table.

    Returns:
        (AnnData, pd.Series, pd.Series): counts with barcodes carrying `suffix`, the fraction of unspliced
        reads indexed by the 16-letter barcode, and the planted population of each 16-letter barcode.
    """
    sizes = DEFAULT_SIZES if sizes is None else sizes
    rng = np.random.default_rng(seed)
    genes = pd.Index(gene_names())

    counts, fractions, populations = [], [], []
    for population, n in sizes.items():
        if population == "doublet":
            block = rng.poisson(_rates(genes, "CM"), size=(n, len(genes))) + rng.poisson(_rates(genes, "VEC"), size=(n, len(genes)))
        else:
            block = rng.poisson(_rates(genes, population), size=(n, len(genes)))
        counts.append(block)
        if population in NUCLEI or population == "doublet":
            fractions.append(rng.normal(0.8, 0.03, size=n))
        elif population == "debris":
            fractions.append(rng.normal(0.25, 0.03, size=n))
        else:
            fractions.append(rng.normal(0.3, 0.05, size=n))
        populations += [population] * n

    barcodes = make_barcodes(len(populations))
    adata = ad.AnnData(sp.csr_matrix(np.vstack(counts).astype(np.float32)))
    adata.obs_names = [barcode + suffix for barcode in barcodes]
    adata.var_names = genes
    fraction_unspliced = pd.Series(np.clip(np.concatenate(fractions), 0, 1), index=barcodes, name="fraction_unspliced")
    return adata, fraction_unspliced, pd.Series(populations, index=barcodes, name="population")
