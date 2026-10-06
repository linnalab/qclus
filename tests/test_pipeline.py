import numpy as np
import pandas as pd
import pytest

import qclus as qc
import qclus.qclus as pipeline
from conftest import write_dataset
from synthetic import NUCLEI, make_barcodes


def test_planted_populations_land_in_expected_filters(dataset, default_run):
    labels = default_run.obs["qclus"]
    populations = dataset.populations.reindex(default_run.obs.index)

    assert list(default_run.obs.index) == list(dataset.fraction_unspliced.index)

    removed_early = populations.isin(["empty", "high_mito"])
    assert (labels[removed_early] == "initial filter").all()
    assert (labels == "initial filter").sum() == removed_early.sum()

    assert (labels[populations == "debris"] == "clustering filter").all()
    assert (labels == "clustering filter").sum() == (populations == "debris").sum()

    nuclei = populations.isin(NUCLEI)
    assert labels[nuclei].isin(["passed", "scrublet filter"]).all()
    assert (labels[nuclei] == "passed").mean() > 0.85


def test_output_keeps_original_barcodes(dataset, default_run):
    assert list(default_run.obs["original_barcode"]) == list(dataset.adata.obs_names)
    assert all(len(barcode) == 16 for barcode in default_run.obs.index)


def test_fraction_table_as_dataframe_and_with_suffixes(dataset, default_run):
    table = dataset.fraction_unspliced.to_frame()
    table.index = [barcode + "-1" for barcode in table.index]
    result = qc.run_qclus(dataset.counts_path, table)
    assert result.obs["qclus"].equals(default_run.obs["qclus"])


def test_barcodes_without_splicing_information_are_dropped(dataset, default_run):
    kept = dataset.fraction_unspliced.iloc[::2]
    result = qc.run_qclus(dataset.counts_path, kept)
    assert list(result.obs.index) == list(kept.index)
    assert set(result.obs["qclus"]) <= set(default_run.obs["qclus"])


def test_no_common_barcodes(dataset):
    other_barcodes = make_barcodes(50, offset=10**6)
    table = pd.Series(0.5, index=other_barcodes)
    with pytest.raises(ValueError, match="No common barcodes"):
        qc.run_qclus(dataset.counts_path, table)


def test_duplicate_barcodes_in_counts(tmp_path, dataset):
    # Two samples in one file: the same 16 letters with different suffixes
    adata = dataset.adata[:20].copy()
    barcodes = [name[:16] for name in adata.obs_names[:10]]
    adata.obs_names = [barcode + "-1" for barcode in barcodes] + [barcode + "-2" for barcode in barcodes]
    path = tmp_path / "two_samples.h5ad"
    adata.write_h5ad(path)
    with pytest.raises(ValueError, match="not unique after truncation"):
        qc.run_qclus(str(path), dataset.fraction_unspliced)


def test_single_cardiomyocyte_feature(dataset):
    result = qc.run_qclus(
        dataset.counts_path,
        dataset.fraction_unspliced,
        clustering_features=["pct_counts_nuclear", "pct_counts_MT", "pct_counts_CM_cyto", "fraction_unspliced"],
    )
    assert "pct_counts_CM_cyto" in result.obs
    assert "pct_counts_CM_nucl" not in result.obs


def test_fraction_unspliced_must_be_a_clustering_feature(dataset):
    with pytest.raises(ValueError, match="must be one of the clustering_features"):
        qc.run_qclus(
            dataset.counts_path,
            dataset.fraction_unspliced,
            clustering_features=["pct_counts_nuclear", "pct_counts_MT"],
        )


def test_gene_set_with_no_genes_present(dataset):
    with pytest.raises(ValueError, match="gene set 'nuclear'"):
        qc.run_qclus(dataset.counts_path, dataset.fraction_unspliced, nucl_gene_set=["Malat1", "Neat1"])


def test_integer_cluster_labels(dataset, default_run):
    result = qc.run_qclus(dataset.counts_path, dataset.fraction_unspliced, clusters_to_select=[0, 1, 2])
    assert result.obs["qclus"].equals(default_run.obs["qclus"])


@pytest.mark.parametrize("selection", [["0", "4"], [], ["first"]])
def test_invalid_cluster_selection(dataset, selection):
    with pytest.raises(ValueError, match="clusters_to_select must be a non-empty subset"):
        qc.run_qclus(dataset.counts_path, dataset.fraction_unspliced, clusters_to_select=selection)


def test_selected_cluster_left_empty_by_kmeans(tmp_path, dataset, monkeypatch):
    # Identical droplets give identical feature vectors, so k-means fills a single cluster
    adata = dataset.adata[:1].copy()
    adata = adata[[0] * 20].copy()
    adata.obs_names = [name + "-1" for name in dataset.fraction_unspliced.index[:20]]
    path = tmp_path / "identical.h5ad"
    adata.write_h5ad(path)
    fractions = pd.Series(0.8, index=dataset.fraction_unspliced.index[:20])
    monkeypatch.setattr(pipeline, "add_qclus_embedding", lambda adata, *args, **kwargs: np.zeros((adata.n_obs, 2)))

    with pytest.raises(ValueError, match="k-means filled only 1 of 2 clusters"):
        qc.run_qclus(
            str(path),
            fractions,
            clustering_features=["pct_counts_nuclear", "pct_counts_MT", "fraction_unspliced"],
            clustering_k=2,
            clusters_to_select=["1"],
            scrublet_filter=False,
        )


def test_no_barcodes_after_initial_filter(tmp_path):
    sample = write_dataset(tmp_path, sizes={"empty": 20})
    with pytest.raises(ValueError, match="Only 0 of 20 barcodes remain after the initial filter"):
        qc.run_qclus(sample.counts_path, sample.fraction_unspliced)


def test_fewer_barcodes_than_clusters(tmp_path):
    sample = write_dataset(tmp_path, sizes={"CM": 3, "empty": 5})
    with pytest.raises(ValueError, match="k-means with clustering_k=4 needs at least 4"):
        qc.run_qclus(sample.counts_path, sample.fraction_unspliced, scrublet_filter=False)


def test_as_many_barcodes_as_scrublet_components(tmp_path):
    sample = write_dataset(tmp_path, sizes={"CM": 5, "VEC": 5})
    with pytest.raises(ValueError, match="the doublet filter with scrublet_n_pcs=10 needs at least 11"):
        qc.run_qclus(sample.counts_path, sample.fraction_unspliced, scrublet_n_pcs=10)


OTHER_TISSUE = dict(
    clustering_features=["pct_counts_nuclear", "pct_counts_MT", "fraction_unspliced"],
    clustering_k=3,
    clusters_to_select=["0", "1"],
    scrublet_filter=False,
)


def test_three_barcodes_are_too_few_for_the_embedding(tmp_path):
    sample = write_dataset(tmp_path, sizes={"CM": 2, "VEC": 1})
    with pytest.raises(ValueError, match="the UMAP embedding needs at least 4"):
        qc.run_qclus(sample.counts_path, sample.fraction_unspliced, **OTHER_TISSUE)


@pytest.mark.parametrize("n_barcodes", [4, 15])
def test_small_inputs_still_run(tmp_path, n_barcodes):
    sample = write_dataset(tmp_path, sizes={"CM": n_barcodes - 2, "debris": 2})
    result = qc.run_qclus(sample.counts_path, sample.fraction_unspliced, **OTHER_TISSUE)
    assert result.n_obs == n_barcodes
    assert result.uns["QClus_umap"].shape == (n_barcodes, 2)


def test_too_few_genes_for_scrublet(dataset):
    with pytest.raises(ValueError, match="Scrublet failed on .* barcodes with scrublet_n_pcs=30") as error:
        qc.run_qclus(
            dataset.counts_path,
            dataset.fraction_unspliced,
            scrublet_minimum_gene_variability_pctl=99.9,
        )
    assert error.value.__cause__ is not None


def test_no_barcode_passing_is_a_warning(dataset):
    with pytest.warns(UserWarning, match="No barcode passed QClus"):
        result = qc.run_qclus(dataset.counts_path, dataset.fraction_unspliced, scrublet_thresh=0.0)
    assert not (result.obs["qclus"] == "passed").any()
    assert (result.obs["qclus"] == "scrublet filter").any()
