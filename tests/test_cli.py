import argparse
import inspect
import json
import subprocess
import sys

import anndata as ad
import pandas as pd
import pytest

import qclus as qc
from qclus.cli import build_parser, main, pipeline_options


@pytest.fixture
def inputs(dataset, tmp_path):
    """Command line arguments for the shared synthetic sample, with the fraction table written to disk."""
    table = tmp_path / "fraction_unspliced.csv"
    dataset.fraction_unspliced.to_frame().to_csv(table)
    return ["run", "--counts", dataset.counts_path, "--fraction-unspliced", str(table)]


def run_options():
    """The options of `qclus run` that correspond to arguments of run_qclus, by destination."""
    run_parser = build_parser()._subparsers._group_actions[0].choices["run"]
    group = next(group for group in run_parser._action_groups if group.title == "pipeline options")
    return {action.dest: action for action in group._group_actions}


def test_options_match_the_run_qclus_signature():
    parameters = dict(inspect.signature(qc.run_qclus).parameters)
    del parameters["counts_path"], parameters["fraction_unspliced"]
    options = run_options()

    assert set(options) == set(parameters)
    for name, parameter in parameters.items():
        option, default = options[name], parameter.default
        assert option.default is argparse.SUPPRESS, name
        if isinstance(default, bool):
            assert isinstance(option, argparse.BooleanOptionalAction), name
        elif isinstance(default, (int, float)):
            assert option.type is type(default), name
        elif isinstance(default, tuple):
            assert option.nargs == "+", name


def test_nothing_is_set_unless_given(inputs):
    args = build_parser().parse_args(inputs + ["-o", "result.h5ad"])
    assert pipeline_options(args) == {}
    assert args.tissue == "heart"


def test_list_and_boolean_options(inputs):
    args = build_parser().parse_args(
        inputs + ["-o", "result.h5ad", "--clusters-to-select", "0", "1", "--clustering-features", "pct_counts_MT",
                  "fraction_unspliced", "--no-scrublet-filter", "--scrublet-approx-neighbors", "--clustering-k", "3"]
    )
    assert pipeline_options(args) == {
        "clusters_to_select": ["0", "1"],
        "clustering_features": ["pct_counts_MT", "fraction_unspliced"],
        "scrublet_filter": False,
        "scrublet_approx_neighbors": True,
        "clustering_k": 3,
    }


def test_help_names_the_presets(capsys):
    with pytest.raises(SystemExit) as exit_info:
        main(["run", "--help"])
    assert exit_info.value.code == 0
    text = " ".join(capsys.readouterr().out.split())
    assert "Number of k-means clusters. Default set by --tissue (heart: 4; other: 3)." in text
    assert "Default set by --tissue (heart: 0 1 2; other: 0 1)." in text
    assert "Default: 500." in text


def test_run_writes_the_same_labels_as_the_function(inputs, default_run, tmp_path, capsys):
    output, table = tmp_path / "result.h5ad", tmp_path / "result.csv"

    assert main(inputs + ["-o", str(output), "--obs-csv", str(table)]) == 0

    written = ad.read_h5ad(output)
    assert list(written.obs["qclus"]) == list(default_run.obs["qclus"])
    assert list(written.obs.index) == list(default_run.obs.index)
    assert "QClus_umap" in written.obsm and "qclus" in written.uns

    obs = pd.read_csv(table, index_col="barcode")
    assert list(obs["qclus"]) == list(default_run.obs["qclus"])
    assert list(obs["original_barcode"]) == list(default_run.obs["original_barcode"])

    printed = dict(line.split("\t") for line in capsys.readouterr().out.strip().splitlines())
    assert int(printed["barcodes"]) == default_run.n_obs
    assert int(printed["passed"]) == (default_run.obs["qclus"] == "passed").sum()


def test_settings_reach_the_pipeline(inputs, dataset, tmp_path):
    table = tmp_path / "result.tsv"
    arguments = ["--minimum-genes", "800", "--no-scrublet-filter", "--no-compute-embedding", "--kmeans-n-init", "3"]
    assert main(inputs + ["--obs-csv", str(table)] + arguments) == 0

    expected = qc.run_qclus(
        dataset.counts_path, dataset.fraction_unspliced,
        minimum_genes=800, scrublet_filter=False, compute_embedding=False, kmeans_n_init=3,
    )
    obs = pd.read_csv(table, sep="\t", index_col="barcode")
    assert list(obs["qclus"]) == list(expected.obs["qclus"])
    assert not (tmp_path / "result.h5ad").exists()


def test_tissue_other_matches_quickstart(inputs, dataset, tmp_path):
    table = tmp_path / "result.csv"
    assert main(inputs + ["--obs-csv", str(table), "--tissue", "other"]) == 0
    expected = qc.quickstart_qclus(dataset.counts_path, dataset.fraction_unspliced, tissue="other")
    assert list(pd.read_csv(table, index_col="barcode")["qclus"]) == list(expected.obs["qclus"])


def test_passed_only(inputs, default_run, tmp_path, capsys):
    output, table = tmp_path / "passed.h5ad", tmp_path / "all.csv"
    assert main(inputs + ["-o", str(output), "--obs-csv", str(table), "--passed-only"]) == 0

    written = ad.read_h5ad(output)
    n_passed = (default_run.obs["qclus"] == "passed").sum()
    assert written.n_obs == n_passed
    assert (written.obs["qclus"] == "passed").all()
    assert written.obsm["QClus_umap"].shape == (n_passed, 2)
    assert "QClus_umap" not in written.uns
    assert len(pd.read_csv(table)) == default_run.n_obs
    assert f"barcodes\t{default_run.n_obs}" in capsys.readouterr().out


def test_gene_set_files(inputs, dataset, tmp_path):
    nuclear, cell_types, table = tmp_path / "nuclear.txt", tmp_path / "cell_types.json", tmp_path / "result.csv"
    nuclear.write_text("MALAT1\nNEAT1\n\nFTX\n")
    cell_types.write_text(json.dumps({"VEC": ["VWF", "ERG"], "FB": ["DCN", "ABCA8"]}))
    arguments = ["--nucl-gene-set", str(nuclear), "--celltype-gene-sets", str(cell_types), "--no-compute-embedding"]
    assert main(inputs + ["--obs-csv", str(table)] + arguments) == 0

    expected = qc.run_qclus(
        dataset.counts_path, dataset.fraction_unspliced, compute_embedding=False,
        nucl_gene_set=["MALAT1", "NEAT1", "FTX"], celltype_gene_set_dict={"VEC": ["VWF", "ERG"], "FB": ["DCN", "ABCA8"]},
    )
    assert list(pd.read_csv(table, index_col="barcode")["qclus"]) == list(expected.obs["qclus"])


def usage_error(arguments, capsys):
    """Run the command line expecting a usage error; return its message."""
    with pytest.raises(SystemExit) as exit_info:
        main(arguments)
    assert exit_info.value.code == 2
    return capsys.readouterr().err


def test_an_output_is_required(inputs, capsys):
    assert "at least one of --output and --obs-csv" in usage_error(inputs, capsys)


def test_existing_output_needs_overwrite(inputs, tmp_path, capsys):
    table = tmp_path / "result.csv"
    table.write_text("keep me")
    assert "already exists; use --overwrite" in usage_error(inputs + ["--obs-csv", str(table)], capsys)
    assert table.read_text() == "keep me"

    assert main(inputs + ["--obs-csv", str(table), "--overwrite", "--no-compute-embedding"]) == 0
    assert "qclus" in pd.read_csv(table).columns


def test_outputs_must_differ(inputs, tmp_path, capsys):
    target = tmp_path / "result.h5ad"
    assert "point to the same file" in usage_error(inputs + ["-o", str(target), "--obs-csv", str(target)], capsys)
    assert not target.exists()


def test_an_input_is_never_overwritten(inputs, dataset, capsys):
    message = usage_error(inputs + ["-o", dataset.counts_path, "--overwrite"], capsys)
    assert "inputs are never overwritten" in message


def test_missing_input_and_missing_directory(inputs, tmp_path, capsys):
    arguments = ["run", "--counts", str(tmp_path / "missing.h5"), "--fraction-unspliced", inputs[4], "-o", str(tmp_path / "r.h5ad")]
    assert "does not exist" in usage_error(arguments, capsys)
    assert "directory" in usage_error(inputs + ["-o", str(tmp_path / "no_such_directory" / "r.h5ad")], capsys)


def test_gene_set_names_with_a_slash_cannot_go_to_h5ad(inputs, tmp_path, capsys):
    cell_types = tmp_path / "cell_types.json"
    cell_types.write_text(json.dumps({"L2/3 IT": ["VWF", "ERG"]}))
    arguments = inputs + ["--celltype-gene-sets", str(cell_types), "-o", str(tmp_path / "result.h5ad")]
    assert "contain '/'" in usage_error(arguments, capsys)
    assert not (tmp_path / "result.h5ad").exists()


def test_errors_while_running(dataset, tmp_path, capsys):
    table = tmp_path / "fraction_unspliced.csv"
    pd.DataFrame({"fraction_unspliced": [0.5]}, index=["TTTTTTTTTTTTTTTT"]).to_csv(table)
    arguments = ["run", "--counts", dataset.counts_path, "--fraction-unspliced", str(table), "--obs-csv", str(tmp_path / "r.csv")]

    assert main(arguments) == 1
    assert capsys.readouterr().err.strip().splitlines()[-1] == "qclus: error: No common barcodes found between counts data and fraction_unspliced."
    assert not (tmp_path / "r.csv").exists()

    with pytest.raises(ValueError, match="No common barcodes"):
        main(arguments + ["--debug"])


def test_version_from_the_module_entry_point():
    result = subprocess.run([sys.executable, "-m", "qclus", "--version"], capture_output=True, text=True)
    assert result.returncode == 0
    assert result.stdout.strip() == f"qclus {qc.__version__}"
