import os

import loompy
import numpy as np
import pytest

import qclus.utils as utils
from synthetic import make_barcodes


def write_loom(path, cell_ids, spliced, unspliced, ambiguous):
    layers = {"": spliced + unspliced + ambiguous, "spliced": spliced, "unspliced": unspliced, "ambiguous": ambiguous}
    loompy.create(
        str(path),
        layers,
        row_attrs={"Gene": np.array([f"gene{i}" for i in range(spliced.shape[0])])},
        col_attrs={"CellID": np.array(cell_ids)},
    )


@pytest.fixture
def counts():
    rng = np.random.default_rng(0)
    spliced, unspliced, ambiguous = (rng.poisson(rate, size=(40, 7)).astype("uint16") for rate in (2.0, 5.0, 0.5))
    # A cell without any counts
    spliced[:, 3] = unspliced[:, 3] = ambiguous[:, 3] = 0
    return spliced, unspliced, ambiguous


@pytest.mark.parametrize("batch_size", [3, 512])
def test_fractions_from_a_read_only_loom_file(tmp_path, counts, batch_size):
    spliced, unspliced, ambiguous = counts
    barcodes = make_barcodes(7)
    path = tmp_path / "sample.loom"
    write_loom(path, [f"sample:{barcode}x" for barcode in barcodes], spliced, unspliced, ambiguous)
    os.chmod(path, 0o444)

    result = utils.fraction_unspliced_from_loom(str(path), batch_size=batch_size)

    total = (spliced + unspliced + ambiguous).sum(axis=0)
    expected = np.divide(unspliced.sum(axis=0), total, out=np.zeros(7), where=total > 0)
    assert list(result.index) == barcodes
    np.testing.assert_array_equal(result["fraction_unspliced"].to_numpy(), expected)
    assert result["fraction_unspliced"].iloc[3] == 0.0


def test_cell_ids_without_a_sample_prefix(tmp_path, counts):
    path = tmp_path / "sample.loom"
    write_loom(path, make_barcodes(7), *counts)
    with pytest.raises(ValueError, match="'sample:barcode' form"):
        utils.fraction_unspliced_from_loom(str(path))


def test_missing_loom_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        utils.fraction_unspliced_from_loom(str(tmp_path / "missing.loom"))


def test_command_line_output_can_be_read_by_the_pipeline(tmp_path, counts):
    import pandas as pd

    from qclus.cli import main

    barcodes = make_barcodes(7)
    path, output = tmp_path / "sample.loom", tmp_path / "fraction_unspliced.csv"
    write_loom(path, [f"sample:{barcode}x" for barcode in barcodes], *counts)

    assert main(["splicing-from-loom", "--loom", str(path), "-o", str(output)]) == 0

    prepared = utils.prepare_fraction_unspliced(pd.read_csv(output, index_col=0))
    assert list(prepared.index) == barcodes
