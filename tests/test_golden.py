"""Exact results on the synthetic sample.

The other tests check that the planted groups land in the expected filters, which leaves room for a
few barcodes to change label unnoticed. These tests compare every barcode with the values stored in
tests/data, so that a change which moves even one barcode fails.

The synthetic sample is generated from a fixed seed. Its fingerprint is checked first: if that test
fails, the sample itself has changed (for example through a new random number generator) and the
other failures in this file say nothing about QClus.

Everything is compared exactly, except the doublet scores. They are computed from nearest neighbours
in a principal component space, and the linear algebra library behind that differs between platforms.
When two neighbours are almost equally near, the choice between them can then go the other way, which
moves one score by one step. A handful of such barcodes is tolerated; a real change in the doublet
step moves many scores.

If results are changed on purpose, regenerate the stored values and review the difference:

    QCLUS_UPDATE_EXPECTED=1 python -m pytest tests/test_golden.py
"""
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

DATA = Path(__file__).parent / "data"
LABELS = DATA / "expected_labels.csv"
FINGERPRINT = DATA / "synthetic_sample.json"
UPDATE = os.environ.get("QCLUS_UPDATE_EXPECTED") == "1"

# The filters that do not involve doublet scores
EARLY_FILTERS = ["initial filter", "clustering filter", "outlier filter"]
# Share of doublet scores that may differ between platforms, and the most a single score may move
MAX_SHARE_OF_SCORES_MOVED = 0.01
MAX_SCORE_MOVE = 0.1


def fingerprint(dataset) -> dict:
    counts = dataset.adata.X.toarray().astype(np.int32)
    fractions = np.round(dataset.fraction_unspliced.to_numpy(), 10)
    return {
        "n_barcodes": int(counts.shape[0]),
        "n_genes": int(counts.shape[1]),
        "total_counts": int(counts.sum()),
        "counts_sha256": hashlib.sha256(counts.tobytes()).hexdigest(),
        "fractions_sha256": hashlib.sha256(fractions.tobytes()).hexdigest(),
    }


def results(default_run) -> pd.DataFrame:
    table = default_run.obs[["kmeans", "qclus", "score_scrublet"]].copy()
    table.index.name = "barcode"
    return table


@pytest.fixture(scope="module")
def expected(dataset, default_run):
    if UPDATE:
        FINGERPRINT.write_text(json.dumps(fingerprint(dataset), indent=2) + "\n")
        results(default_run).to_csv(LABELS)
    return pd.read_csv(LABELS, index_col="barcode", dtype={"kmeans": str, "qclus": str})


def differing(actual: pd.Series, wanted: pd.Series) -> str:
    changed = actual.astype(str).to_numpy() != wanted.astype(str).to_numpy()
    moves = pd.crosstab(wanted[changed].astype(str), actual[changed].astype(str), rownames=["expected"], colnames=["got"])
    return f"{int(changed.sum())} of {len(wanted)} barcodes differ:\n{moves}"


def test_synthetic_sample_is_unchanged(dataset, expected):
    assert fingerprint(dataset) == json.loads(FINGERPRINT.read_text())


def test_kmeans_labels_are_unchanged(default_run, expected):
    actual = results(default_run)
    assert list(actual.index) == list(expected.index)
    assert (actual["kmeans"] == expected["kmeans"]).all(), differing(actual["kmeans"], expected["kmeans"])


def test_first_three_filters_are_unchanged(default_run, expected):
    actual = results(default_run)["qclus"].where(lambda labels: labels.isin(EARLY_FILTERS), "later")
    wanted = expected["qclus"].where(lambda labels: labels.isin(EARLY_FILTERS), "later")
    assert (actual == wanted).all(), differing(actual, wanted)


def score_moved(default_run, expected) -> pd.Series:
    """Per barcode, whether its doublet score differs from the stored one."""
    actual, wanted = results(default_run)["score_scrublet"], expected["score_scrublet"]
    assert (actual.isna() == wanted.isna()).all(), "different barcodes were scored"
    return (actual - wanted).abs().fillna(0) > 1e-6


def test_doublet_scores_are_unchanged(default_run, expected):
    actual, wanted = results(default_run)["score_scrublet"], expected["score_scrublet"]
    moved = score_moved(default_run, expected)
    scored = wanted.notna()

    assert moved[scored].mean() <= MAX_SHARE_OF_SCORES_MOVED, f"{int(moved.sum())} of {int(scored.sum())} scores differ"
    assert (actual - wanted).abs().max() <= MAX_SCORE_MOVE


def test_labels_are_unchanged(default_run, expected):
    # A barcode whose score moved may end up on the other side of the threshold; all others must agree
    comparable = ~score_moved(default_run, expected)
    actual, wanted = results(default_run)["qclus"][comparable], expected["qclus"][comparable]
    assert (actual == wanted).all(), differing(actual, wanted)
