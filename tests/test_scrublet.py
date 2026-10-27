import numpy as np
import pytest

import qclus as qc
import qclus.utils as utils
from synthetic import NUCLEI


def fake_scrublet(scores):
    """A stand-in for scrublet.Scrublet that returns the given score for every barcode."""

    class FakeScrublet:
        def __init__(self, counts, expected_doublet_rate):
            self.n_obs = counts.shape[0]

        def scrub_doublets(self, **kwargs):
            self.doublet_scores_obs_ = np.resize(np.asarray(scores, dtype=float), self.n_obs)
            return self.doublet_scores_obs_, None

    return FakeScrublet


def test_identical_scores_are_an_error(dataset, monkeypatch):
    monkeypatch.setattr("scrublet.Scrublet", fake_scrublet([0.25]))
    with pytest.raises(RuntimeError, match="same doublet score") as error:
        utils.calculate_scrublet(dataset.adata)
    assert "scrublet_filter=False" in str(error.value)
    assert "annoy" not in str(error.value)


def test_scores_that_are_not_finite_are_an_error(dataset, monkeypatch):
    monkeypatch.setattr("scrublet.Scrublet", fake_scrublet([0.1, np.nan, 0.3]))
    with pytest.raises(RuntimeError, match="not finite"):
        utils.calculate_scrublet(dataset.adata)


def test_broken_annoy_is_reported_when_approximate_neighbours_are_requested(dataset, monkeypatch):
    monkeypatch.setattr(utils, "annoy_works", lambda: False)
    with pytest.raises(RuntimeError, match="annoy library returns wrong neighbours") as error:
        utils.calculate_scrublet(dataset.adata, approx_neighbors=True)
    assert "scrublet_approx_neighbors=False" in str(error.value)
    assert "conda-forge" in str(error.value)


def test_exact_neighbours_do_not_depend_on_annoy(dataset, default_run, monkeypatch):
    def fail():
        raise AssertionError("annoy must not be tested when exact neighbours are used")

    monkeypatch.setattr(utils, "annoy_works", fail)
    result = qc.run_qclus(dataset.counts_path, dataset.fraction_unspliced)
    assert result.obs["qclus"].equals(default_run.obs["qclus"])


def test_scores_are_in_the_output(default_run):
    scores = default_run.obs["score_scrublet"]
    removed_early = default_run.obs["qclus"] == "initial filter"
    assert scores[removed_early].isna().all()
    assert np.isfinite(scores[~removed_early]).all()


def test_planted_doublets_score_higher_than_singlets(dataset, default_run):
    scores = default_run.obs["score_scrublet"]
    populations = dataset.populations.reindex(default_run.obs.index)
    assert scores[populations == "doublet"].median() > scores[populations.isin(NUCLEI)].median()


def test_approximate_neighbours(dataset):
    if not utils.annoy_works():
        pytest.skip("annoy returns wrong neighbours in this environment")
    result = qc.run_qclus(dataset.counts_path, dataset.fraction_unspliced, scrublet_approx_neighbors=True)
    scores = result.obs["score_scrublet"].dropna()
    assert scores.nunique() > 1


@pytest.mark.parametrize("n_init", [1, 10])
def test_kmeans_restarts_reach_scikit_learn(dataset, monkeypatch, n_init):
    seen = {}
    original = utils.KMeans

    def spy(*args, **kwargs):
        seen.update(kwargs)
        return original(*args, **kwargs)

    monkeypatch.setattr(utils, "KMeans", spy)
    qc.run_qclus(dataset.counts_path, dataset.fraction_unspliced, kmeans_n_init=n_init)
    assert seen["n_init"] == n_init


def test_kmeans_runs_one_restart_by_default(dataset, monkeypatch):
    seen = {}
    original = utils.KMeans

    def spy(*args, **kwargs):
        seen.update(kwargs)
        return original(*args, **kwargs)

    monkeypatch.setattr(utils, "KMeans", spy)
    qc.run_qclus(dataset.counts_path, dataset.fraction_unspliced)
    assert seen["n_init"] == 1
