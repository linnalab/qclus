import numpy as np
import pandas as pd
import pytest

import qclus.utils as utils
from synthetic import make_barcodes


def test_kmeans_labels_are_ordered_by_fraction_unspliced():
    rng = np.random.default_rng(0)
    features = pd.DataFrame(
        {
            "fraction_unspliced": np.concatenate([rng.normal(0.2, 0.01, 30), rng.normal(0.8, 0.01, 30), rng.normal(0.5, 0.01, 30)]),
            "pct_counts_MT": np.concatenate([rng.normal(30, 1, 30), rng.normal(1, 0.1, 30), rng.normal(10, 1, 30)]),
        }
    )
    labels = np.array(utils.do_kmeans(features.copy(), k=3))
    assert set(labels[:30]) == {"2"}
    assert set(labels[30:60]) == {"0"}
    assert set(labels[60:]) == {"1"}


def test_kmeans_orders_empty_clusters_last():
    features = pd.DataFrame({"fraction_unspliced": [0.5] * 10, "pct_counts_MT": [1.0] * 10})
    with pytest.warns(Warning):
        labels = utils.do_kmeans(features, k=3)
    assert set(labels) == {"0"}


def test_outlier_thresholds():
    table = pd.DataFrame(
        {
            "fraction_unspliced": [0.80, 0.80, 0.80, 0.80, 0.71, 0.69, 0.75, 0.75],
            "pct_counts_MT": [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 6.9, 7.1],
            "kmeans": ["0", "0", "0", "0", "1", "1", "1", "1"],
        }
    )
    outliers = utils.annotate_outliers(table, unspliced_diff=0.1, mito_diff=5.0)
    # Cluster 0 has the highest mean fraction_unspliced, so it sets the thresholds:
    # fraction_unspliced above 0.80 - 0.1 and pct_counts_MT below 2.0 + 5.0
    assert list(outliers) == [False, False, False, False, False, True, False, True]


def test_count_file_errors(tmp_path):
    with pytest.raises(FileNotFoundError):
        utils.read_count_file(str(tmp_path / "missing.h5ad"))

    unsupported = tmp_path / "counts.txt"
    unsupported.write_text("not a count matrix")
    with pytest.raises(ValueError, match="Unsupported file format"):
        utils.read_count_file(str(unsupported))

    corrupt = tmp_path / "corrupt.h5ad"
    corrupt.write_text("not an h5ad file")
    with pytest.raises(IOError, match="Failed to read counts file") as error:
        utils.read_count_file(corrupt)
    assert error.value.__cause__ is not None


def fractions(n=5, **kwargs):
    return pd.Series(np.linspace(0.1, 0.9, n), index=make_barcodes(n), **kwargs)


def test_fraction_unspliced_accepts_series_and_dataframes():
    series = fractions()
    expected = series.rename("fraction_unspliced")

    pd.testing.assert_series_equal(utils.prepare_fraction_unspliced(series), expected)
    pd.testing.assert_series_equal(utils.prepare_fraction_unspliced(series.to_frame("fraction_unspliced")), expected)
    pd.testing.assert_series_equal(utils.prepare_fraction_unspliced(series.to_frame("anything")), expected)

    two_columns = pd.DataFrame({"other": 1.0, "fraction_unspliced": series})
    pd.testing.assert_series_equal(utils.prepare_fraction_unspliced(two_columns), expected)


def test_fraction_unspliced_barcodes_are_truncated_without_touching_the_input():
    series = fractions()
    with_suffix = series.copy()
    with_suffix.index = [barcode + "-1" for barcode in series.index]

    prepared = utils.prepare_fraction_unspliced(with_suffix)

    assert list(prepared.index) == list(series.index)
    assert with_suffix.index[0].endswith("-1")


@pytest.mark.parametrize(
    "table, error, message",
    [
        (pd.DataFrame({"a": [0.1], "b": [0.2]}, index=make_barcodes(1)), ValueError, "none is named 'fraction_unspliced'"),
        (pd.Series([0.1, np.nan], index=make_barcodes(2)), ValueError, "1 missing values"),
        (pd.Series([0.1, 55.0], index=make_barcodes(2)), ValueError, "must lie between 0 and 1"),
        (pd.Series([0.1, -0.2], index=make_barcodes(2)), ValueError, "must lie between 0 and 1"),
        (pd.Series(["low", "high"], index=make_barcodes(2)), ValueError, "must be numeric"),
        (pd.Series([0.1, 0.2], index=["ACGTACGTACGTACGT-1", "ACGTACGTACGTACGT-2"]), ValueError, "not unique after truncation"),
        ({"ACGT": 0.5}, TypeError, "must be a pandas Series or DataFrame"),
    ],
)
def test_fraction_unspliced_is_validated(table, error, message):
    with pytest.raises(error, match=message):
        utils.prepare_fraction_unspliced(table)


@pytest.mark.parametrize("container", [list, tuple, set, pd.Index])
def test_gene_sets_can_be_any_sequence(dataset, container):
    adata = dataset.adata[:50].copy()
    utils.get_qc_metrics(adata, container(["MT-ND1", "MT-CO1", "NOT-A-GENE"]), "MT")
    assert (adata.obs["pct_counts_MT"] > 0).any()


def test_gene_set_given_as_a_string(dataset):
    with pytest.raises(TypeError, match="not a single string"):
        utils.get_qc_metrics(dataset.adata[:50].copy(), "MT-ND1", "MT")


def test_available_cpus():
    assert utils.available_cpus() >= 1
