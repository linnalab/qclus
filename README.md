# QClus

This is a novel droplet filtering method targeted to streamline the processing of snRNA-seq samples. The algorithm uses various metrics to cluster nuclei and filter empty and highly contaminated droplets. This approach combined with other filtering steps enables for flexible, automated, and reliable preprocessing of samples with varying number of nuclei, quality, and contamination levels. The robustness of this method has been validated on a large number of heterogeneous datasets, in terms of number of nuclei, overall quality, and contamination.

Additionally, while the method was originally developed for cardiac snRNA-seq data, in our paper we show that given its flexible design it can also be applied to other tissues. The tutorials below cover both cases.

Any and all comments/criticisms/suggestions are enthusiastically received! :-)

## Installation

QClus needs Python 3.11 or newer and runs on Linux and macOS. On Windows, please use WSL.

### With conda

Note: In order to use our environment installation script, you need to have conda ([Anaconda](https://docs.anaconda.com/anaconda/install/)/[Miniconda](https://docs.anaconda.com/miniconda/install/)) installed on your machine. 

1. Clone this repository into a suitable location on your machine or server using the following command:

    ```git clone https://github.com/linnalab/qclus.git```
    
2. In the root directory of the package, run the following command to create an environment named ```qclus``` and install the required packages:

    ``` bash environment.sh```

### With pip

QClus can also be installed straight from this repository:

```pip install git+https://github.com/linnalab/qclus.git```

To also install what the tutorial notebooks use (JupyterLab, matplotlib, leidenalg and igraph):

```pip install "qclus[tutorials] @ git+https://github.com/linnalab/qclus.git"```

pip compiles one dependency (`annoy`) from source, so a C++ compiler has to be available. A release on PyPI is planned.


## Getting Started

You can find tutorials on how to use QClus in the `tutorials` directory. They are written in Jupyter notebooks. 

In order to run QClus, you will need two things:

1. The count matrix of your droplets: the path to an .h5 or .h5ad file, or an AnnData object you have already loaded. It has to hold raw counts. In our benchmarking we begin with the filtered 10X output, but in principle the unfiltered output can also be used.
2. The fraction of unspliced reads for each droplet, e.g. `fraction_unspliced.csv`. They should be in the following format:

|              | fraction_unspliced |
|----------------------|---------------------|
| AAACGCTCAGATACCT     | 0.259897           |
| AAACCCAGTGAATGTA     | 0.774209           |
| AAACGAACATTGACCA     | 0.526074           |
| AAACCCAAGCAGCACA     | 0.832765           |
| AAACGAATCCACAGGC     | 0.794653           |


### How to get the fraction of unspliced reads for each droplet

Note: If you have these values already you can skip this step.

#### Option 1: From FASTQ processing

Many raw data processing methods, such as [alevin-fry](https://alevin-fry.readthedocs.io/en/latest/) and [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md) output splicing information automatically. They calculate spliced, unspliced, and ambiguous counts per cell per gene, from which you can easily calculate the total fraction for each droplet.

#### Option 2: From [Velocyto](https://velocyto.org/) output (.loom file)

If you have previously run Velocyto for your data, we have a tutorial for you! Move to the [splicing_from_loompy.ipynb notebook](https://github.com/linnalab/qclus/blob/main/tutorials/splicing_from_loompy.ipynb), which will show you how to get the unspliced values from the .loom file.

#### Option 3: From 10X BAM files

If you do not have the necessary values from the methods above, we have implemented our own method for calculating unspliced fraction directly from your 10X bam files! Move to the [splicing_from_bam.ipynb notebook](https://github.com/linnalab/qclus/blob/main/tutorials/splicing_from_bam.ipynb), which will show you how to get the unspliced values from the bam files.

### QClus quickstart

From here you can run the `quickstart_qclus` function. For cardiac data, the function runs QClus with the same settings as in the original publication.

```python
import qclus as qc
import pandas as pd

counts_path = 'path/to/your/count_matrix.h5ad'  
fraction_unspliced = pd.read_csv("path/to/your/fraction_unspliced.csv", index_col=0)

adata = qc.quickstart_qclus(counts_path, 
                 fraction_unspliced, 
                 tissue = 'heart')
```

For other tissues, you can use the following function:

```python
adata = qc.quickstart_qclus(counts_path, 
                 fraction_unspliced, 
                 tissue = 'other')
```

This executes QClus without cell type specific metrics. For fine-tuning non-cardiac data, you can check the tutorials below.

The returned object holds the droplets that have splicing information. `adata.obs["qclus"]` is `passed` for droplets that pass, and otherwise names the filter that removed the droplet.

### Running QClus from the command line

Installing QClus also installs a `qclus` command, which processes one sample per call:

```
qclus run \
    --counts filtered_feature_bc_matrix.h5 \
    --fraction-unspliced fraction_unspliced.csv \
    --output sample_qclus.h5ad \
    --obs-csv sample_qclus.csv
```

`--output` writes the annotated AnnData object, and `--obs-csv` writes the per-droplet table on its own, which is all you need to filter the sample in R or Seurat. At least one of the two is required. Use `--tissue other` for non-cardiac data and `--passed-only` to keep only the droplets that pass. `qclus run --help` lists every setting.

The fraction of unspliced reads can be calculated from the command line as well:

```
qclus splicing-from-bam \
    --bam possorted_genome_bam.bam \
    --barcodes filtered_feature_bc_matrix/barcodes.tsv.gz \
    --output fraction_unspliced.csv

qclus splicing-from-loom --loom sample.loom --output fraction_unspliced.csv
```


### Running QClus in a Scanpy workflow

We also provide tutorial notebooks ([heart](https://github.com/linnalab/qclus/blob/main/tutorials/qclus_tutorial_heart.ipynb) and [brain](https://github.com/linnalab/qclus/blob/main/tutorials/qclus_tutorial_brain.ipynb)) which show how to run QClus as part of a simple Scanpy workflow and how to evaluate the results.

## How to cite

We hope you find the QClus package useful for your research! If you do, please remember to cite [our paper in Nucleic Acids Research](https://doi.org/10.1093/nar/gkae1145):

Eloi Schmauch, Johannes Ojanen, Kyriakitsa Galani, Juho Jalkanen, Kristiina Harju, Maija Hollmén, Hannu Kokki, Jarmo Gunn, Jari Halonen, Juha Hartikainen, Tuomas Kiviniemi, Pasi Tavi, Minna U Kaikkonen, Manolis Kellis, Suvi Linna-Kuosmanen, QClus: a droplet filtering algorithm for enhanced snRNA-seq data quality in challenging samples, Nucleic Acids Research, Volume 53, Issue 1, 13 January 2025, gkae1145, https://doi.org/10.1093/nar/gkae1145

You can also find the latest preprint for this method [here](https://www.biorxiv.org/content/10.1101/2022.10.21.513315v2).

## If You Have Issues

If you have any issues with the installation or running the tutorials, please open an issue on this repository. We will do our best to help you out as soon as possible!


