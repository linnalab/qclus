"""Compare a run on a real sample against a reference captured from QClus 0.2.0.

The sample is participant data and is not part of the repository, so these tests are skipped
unless QCLUS_TEST_DATA points to a directory holding the sample and its reference, and
QCLUS_TEST_SAMPLE gives the name of the sample:

    <sample>.h5
    <sample>_fraction_unspliced.csv
    qclus_regression/reference_obs.csv.gz
    qclus_regression/reference_extra.csv.gz
    qclus_regression/reference_meta.json

Exact agreement is expected in the environment the reference was captured in. Failure messages
report counts only, never barcodes.
"""
import importlib.metadata
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import qclus as qc
import qclus.utils as utils

DATA = Path(os.environ.get("QCLUS_TEST_DATA", "missing"))
REFERENCE = DATA / "qclus_regression"
SAMPLE = os.environ.get("QCLUS_TEST_SAMPLE", "sample")

pytestmark = pytest.mark.skipif(
    not ((REFERENCE / "reference_obs.csv.gz").exists() and (DATA / f"{SAMPLE}.h5").exists()),
    reason="QCLUS_TEST_DATA and QCLUS_TEST_SAMPLE do not point to the sample and its reference",
)


@pytest.fixture(scope="module")
def reference():
    obs = pd.read_csv(REFERENCE / "reference_obs.csv.gz", index_col=0, dtype={"kmeans": str})
    extra = pd.read_csv(REFERENCE / "reference_extra.csv.gz", index_col=0)
    for column in extra.columns:
        if column.startswith("kmeans"):
            extra[column] = extra[column].astype(str)
    meta = json.loads((REFERENCE / "reference_meta.json").read_text())

    changed = {
        package: (version, importlib.metadata.version(package))
        for package, version in meta["packages"].items()
        if package in {"numpy", "scipy", "scikit-learn", "scanpy", "anndata", "pandas", "umap-learn", "scrublet", "annoy"}
        and importlib.metadata.version(package) != version
    }
    if changed:
        print(f"Not the reference environment (reference, here): {changed}")
    return obs, extra, meta


def run(**kwargs):
    fractions = pd.read_csv(DATA / f"{SAMPLE}_fraction_unspliced.csv", index_col=0)
    return qc.run_qclus(str(DATA / f"{SAMPLE}.h5"), fractions, **kwargs)


def assert_same_labels(result, expected, what):
    differing = int((result.astype(str).to_numpy() != expected.astype(str).to_numpy()).sum())
    assert differing == 0, f"{differing} of {len(expected)} barcodes differ in {what}"


def expected_labels(obs, extra, column):
    labels = pd.Series("initial filter", index=obs.index)
    labels[extra.index] = extra[column]
    return labels


def assert_features_unchanged(result, obs):
    assert list(result.obs.index) == list(obs.index)
    for column in obs.columns:
        if pd.api.types.is_numeric_dtype(obs[column]):
            np.testing.assert_allclose(
                result.obs[column].to_numpy(dtype=float), obs[column].to_numpy(dtype=float),
                rtol=1e-6, atol=1e-9, err_msg=column,
            )


def test_approximate_neighbours_reproduce_the_reference(reference):
    if not utils.annoy_works():
        pytest.skip("annoy returns wrong neighbours in this environment")
    obs, extra, _ = reference
    result = run(scrublet_approx_neighbors=True)

    assert_features_unchanged(result, obs)
    assert_same_labels(result.obs["kmeans"], obs["kmeans"], "k-means labels")
    assert_same_labels(result.obs["qclus"], obs["qclus"], "QClus labels")
    np.testing.assert_allclose(
        result.obs.loc[extra.index, "score_scrublet"], extra["score_scrublet_approx"], rtol=0, atol=1e-9
    )


def test_defaults_change_only_the_doublet_step(reference):
    obs, extra, _ = reference
    result = run()

    assert_features_unchanged(result, obs)
    assert_same_labels(result.obs["kmeans"], obs["kmeans"], "k-means labels")
    np.testing.assert_allclose(
        result.obs.loc[extra.index, "score_scrublet"], extra["score_scrublet_exact"], rtol=0, atol=1e-9
    )
    assert_same_labels(result.obs["qclus"], expected_labels(obs, extra, "qclus_n_init_1_exact"), "QClus labels")

    # The first three filters are untouched: a barcode they removed keeps its label, and no other gets one
    early = ["initial filter", "clustering filter", "outlier filter"]
    before, after = obs["qclus"].where(obs["qclus"].isin(early)), result.obs["qclus"].where(result.obs["qclus"].isin(early))
    assert_same_labels(after.fillna("later"), before.fillna("later"), "the initial, clustering and outlier labels")


def test_ten_restarts_match_the_reference(reference):
    obs, extra, _ = reference
    result = run(kmeans_n_init=10)

    assert_same_labels(result.obs.loc[extra.index, "kmeans"], extra["kmeans_n_init_10"], "k-means labels")
    assert_same_labels(result.obs["qclus"], expected_labels(obs, extra, "qclus_n_init_10_exact"), "QClus labels")
