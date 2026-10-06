from types import SimpleNamespace

import pytest

import qclus as qc
from synthetic import make_dataset


def write_dataset(directory, sizes=None, seed=0, suffix="-1"):
    """Write a synthetic count matrix to disk and return it with its fraction_unspliced table."""
    adata, fraction_unspliced, populations = make_dataset(sizes=sizes, seed=seed, suffix=suffix)
    counts_path = directory / "counts.h5ad"
    adata.write_h5ad(counts_path)
    return SimpleNamespace(
        adata=adata,
        counts_path=str(counts_path),
        fraction_unspliced=fraction_unspliced,
        populations=populations,
    )


@pytest.fixture(scope="session")
def dataset(tmp_path_factory):
    """The default synthetic sample. Shared by all tests, so tests must not modify it."""
    return write_dataset(tmp_path_factory.mktemp("dataset"))


@pytest.fixture(scope="session")
def default_run(dataset):
    """Output of run_qclus with every argument at its default. Shared, so tests must not modify it."""
    return qc.run_qclus(dataset.counts_path, dataset.fraction_unspliced)
