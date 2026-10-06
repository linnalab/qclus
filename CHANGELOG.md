# Changelog


All notable changes to QClus are listed here, newest first.



## 0.3.0 (unreleased)



### Packaging


- `pyproject.toml` replaces `setup.py` and `requirements.txt`.

- The version is defined in one place and available as `qclus.__version__`. Earlier releases reported 0.1.0 whatever the tag.

- JupyterLab, leidenalg and igraph are no longer installed with QClus, which does not use them. Install them for the tutorials with the `tutorials` extra.

- Python 3.11 or newer is required, and the dependencies have tested lower bounds.

- Tests run on GitHub Actions for Python 3.11 to 3.13 on Linux and for Python 3.12 on macOS.
