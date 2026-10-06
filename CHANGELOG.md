# Changelog


All notable changes to QClus are listed here, newest first.



## 0.3.0 (unreleased)


The method is unchanged: the same features, the same four filters in the same order, and the same thresholds and gene sets. One default changes results, and it is listed first.



### Changed results


- **Doublet scores use exact nearest neighbours.** Scrublet's approximate search relies on the `annoy` library, which we have seen return wrong neighbours when pip builds it on an Apple Silicon Mac. Every barcode then got the same score and the doublet filter removed nothing, without a warning. On one heart sample the change alters 249 of the 10,860 doublet calls (2.3%); no other label changes. `scrublet_approx_neighbors=True` restores the previous behaviour where annoy works.

- **Reads on a tile boundary are counted once.** `fraction_unspliced_from_bam` counted a read twice when it spanned two of the tiles it generates. User-supplied `regions` keep their meaning: every read overlapping a region is counted.



### Added


- `kmeans_n_init` sets the number of k-means restarts. The default of 1 is what scikit-learn 1.4 and later already ran, so labels do not change.

- The output has two new columns: `score_scrublet`, the doublet score, and `original_barcode`, the barcode as it appears in the counts file.

- `fraction_unspliced` can be a Series or a DataFrame, and its barcodes may carry a suffix such as `-1`.

- `bam_index_path` is optional in `fraction_unspliced_from_bam`, and `fraction_unspliced_from_loom` takes a `batch_size`.



### Fixed


- `fraction_unspliced_from_bam` left missing values in its result under pandas 3, ignored `bam_index_path` in its worker processes, searched the barcodes as a list for every read, and used every core of a shared node instead of the cores given to the job.

- `fraction_unspliced_from_loom` loaded three dense matrices into memory and needed write access to the file. On a sample of 25,096 barcodes its peak memory falls from 2.7 GB to 0.6 GB, with identical results.

- Requesting only one of the two cardiomyocyte features failed.

- `read_count_file` reported an unsupported file format as an `IOError` instead of a `ValueError`.

- Integer cluster labels in `clusters_to_select` matched nothing and ended in an unrelated error. They are now accepted.



### Now rejected, with an error that says why


- `fraction_unspliced` values that are missing, not numeric or outside 0 to 1, and barcodes that are not unique after truncation to 16 characters.

- `clusters_to_select` entries that k-means cannot produce, and clustering features without `fraction_unspliced`.

- A gene set with none of its genes in the data, as happens with mouse symbols or Ensembl IDs.

- Too few barcodes after the initial filter, and too few genes for Scrublet.

- Doublet scores that are not finite or are identical for every barcode. `scrublet_filter=False` skips the doublet step.



### Packaging


- `pyproject.toml` replaces `setup.py` and `requirements.txt`.

- The version is defined in one place and available as `qclus.__version__`. Earlier releases reported 0.1.0 whatever the tag.

- JupyterLab, leidenalg and igraph are no longer installed with QClus, which does not use them. Install them for the tutorials with the `tutorials` extra.

- Python 3.11 or newer is required, and the dependencies have tested lower bounds.

- Tests run on GitHub Actions for Python 3.11 to 3.13 on Linux and for Python 3.12 on macOS.
