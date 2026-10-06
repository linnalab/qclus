import json
import subprocess
import sys
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import qclus as qc
from qclus.gene_lists import nucl_30

HEART_FEATURES = [
    "pct_counts_nonCM",
    "pct_counts_nuclear",
    "pct_counts_MT",
    "pct_counts_CM_cyto",
    "pct_counts_CM_nucl",
    "fraction_unspliced",
]


def test_default_run_raises_no_view_warnings(dataset):
    with warnings.catch_warnings():
        warnings.simplefilter("error", ad.ImplicitModificationWarning)
        qc.run_qclus(dataset.counts_path, dataset.fraction_unspliced)


@pytest.mark.parametrize("container", [list, tuple])
def test_features_as_list_or_tuple(dataset, default_run, container):
    result = qc.run_qclus(dataset.counts_path, dataset.fraction_unspliced, clustering_features=container(HEART_FEATURES))
    assert result.obs["qclus"].equals(default_run.obs["qclus"])


def test_anndata_input_gives_the_same_labels_and_is_not_modified(dataset, default_run):
    adata = dataset.adata.copy()
    result = qc.run_qclus(adata, dataset.fraction_unspliced)

    assert result.obs["qclus"].equals(default_run.obs["qclus"])
    assert list(adata.obs_names) == list(dataset.adata.obs_names)
    assert list(adata.obs.columns) == list(dataset.adata.obs.columns)
    assert list(adata.var.columns) == list(dataset.adata.var.columns)
    assert (adata.X != dataset.adata.X).nnz == 0
    assert not adata.uns and not adata.obsm and not adata.layers


def test_normalised_input_is_flagged(dataset):
    adata = dataset.adata.copy()
    adata.X = adata.X * 0.5
    with pytest.warns(UserWarning, match="not whole numbers"):
        qc.run_qclus(adata, dataset.fraction_unspliced, scrublet_filter=False, compute_embedding=False)


def test_embedding_is_aligned_with_the_barcodes(default_run):
    aligned, compact = default_run.obsm["QClus_umap"], default_run.uns["QClus_umap"]
    scored = (default_run.obs["qclus"] != "initial filter").to_numpy()

    assert aligned.shape == (default_run.n_obs, 2)
    np.testing.assert_array_equal(aligned[scored], compact)
    assert np.isnan(aligned[~scored]).all()


def test_embedding_can_be_skipped(dataset, default_run):
    result = qc.run_qclus(dataset.counts_path, dataset.fraction_unspliced, compute_embedding=False)
    assert "QClus_umap" not in result.uns
    assert "QClus_umap" not in result.obsm
    assert result.obs["qclus"].equals(default_run.obs["qclus"])


def test_provenance_survives_an_h5ad_round_trip(tmp_path, dataset):
    gene_sets = {"Endo cells": ["VWF", "ERG", "ANO2"], "FB": ("DCN", "ABCA8")}
    result = qc.run_qclus(
        dataset.counts_path,
        dataset.fraction_unspliced,
        celltype_gene_set_dict=gene_sets,
        clustering_k=np.int64(4),
    )
    path = tmp_path / "result.h5ad"
    result.write_h5ad(path)
    record = ad.read_h5ad(path).uns["qclus"]

    assert record["version"] == qc.__version__
    params = json.loads(record["params_json"])
    assert params["clustering_k"] == 4
    assert params["clusters_to_select"] == ["0", "1", "2"]
    assert params["celltype_gene_set_dict"] == {"Endo cells": ["VWF", "ERG", "ANO2"], "FB": ["DCN", "ABCA8"]}
    assert params["nucl_gene_set"] == list(nucl_30)
    assert params["scrublet_approx_neighbors"] is False
    assert dict(record["n_barcodes"]) == result.obs["qclus"].value_counts().to_dict()
    assert record["n_barcodes_without_splicing_information"] == 0
    assert 0 < record["outlier_unspliced_threshold"] < 1
    assert record["outlier_mito_threshold"] > 5


def test_outlier_thresholds_are_missing_when_the_filter_is_off(dataset):
    result = qc.run_qclus(dataset.counts_path, dataset.fraction_unspliced, outlier_filter=False, compute_embedding=False)
    assert np.isnan(result.uns["qclus"]["outlier_unspliced_threshold"])
    assert "outlier filter" not in set(result.obs["qclus"])


def test_quickstart_applies_the_preset_and_accepts_overrides(dataset):
    preset = qc.TISSUE_PRESETS["other"]
    expected = qc.run_qclus(dataset.counts_path, dataset.fraction_unspliced, **preset)
    result = qc.quickstart_qclus(dataset.counts_path, dataset.fraction_unspliced, tissue="other")
    assert result.obs["qclus"].equals(expected.obs["qclus"])

    # Arguments given explicitly win over the preset, and the rest of the preset still applies
    overridden = qc.quickstart_qclus(
        dataset.counts_path, dataset.fraction_unspliced, tissue="other", clusters_to_select=["0"], compute_embedding=False
    )
    params = json.loads(overridden.uns["qclus"]["params_json"])
    assert params["clusters_to_select"] == ["0"]
    assert params["clustering_k"] == preset["clustering_k"]
    assert "QClus_umap" not in overridden.uns


def test_quickstart_rejects_unknown_tissues(dataset):
    with pytest.raises(ValueError, match="Invalid tissue: liver"):
        qc.quickstart_qclus(dataset.counts_path, dataset.fraction_unspliced, tissue="liver")


def test_do_kmeans_leaves_its_input_alone():
    features = pd.DataFrame({"fraction_unspliced": [0.1, 0.2, 0.8, 0.9], "pct_counts_MT": [30.0, 28.0, 1.0, 2.0]})
    before = features.copy()
    qc.utils.do_kmeans(features, k=2)
    pd.testing.assert_frame_equal(features, before)


def test_heavy_packages_are_imported_only_when_used():
    # Compare with what the packages QClus needs at import time load by themselves
    script = (
        "import sys\n"
        "import scanpy, anndata, pandas, numpy, tqdm, sklearn.cluster, sklearn.preprocessing\n"
        "before = set(sys.modules)\n"
        "import qclus\n"
        "print(sorted({'loompy', 'pysam', 'scrublet', 'umap'} & (set(sys.modules) - before)))\n"
    )
    output = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, check=True).stdout
    assert output.strip() == "[]"
