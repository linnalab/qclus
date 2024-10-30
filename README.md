# Robust and reliable preprocessing of snRNA-seq data using QClus

This is a novel nuclei filtering method targeted to streamlien the processing of challenging snRNA-seq samples. Originally developed to process cardiac data, we use metrics such as splicing, mitochondrial gene expression, nuclear gene expression, and non-cardiomyocyte and cardiomyocyte marker gene expression to cluster nuclei and filter empty and highly contaminated droplets. This approach combined with other filtering steps enables for flexible, automated, and reliable cleaning of samples with varying number of nuclei, quality, and contamination levels. The robustness of this method has been validated on a large number of heterogeneous datasets, in terms of number of nuclei, overall quality, and contamination.

Additionally, while the method was originally developed for cardiac snRNA-seq data, in our forthcoming paper we show that given it's felxible design it can also be applied to other tissues. We will provide detailed tutorials for this process.

Any and all comments/criticisms/suggestions are enthusiastically received! :-)


## Table of Contents

- [Preprint](#preprint)
- [Installation](#installation)
- [Getting Started](#getting-started)
- [If You Have Issues](#if-you-have-issues)


## Preprint

You can find the current preprint for this method at the link below. Please cite this preprint if you use QClus in your research.

https://www.biorxiv.org/content/10.1101/2022.10.21.513315v2

Figure 2 from the preprint is shown below:

![Optional Alt Text](figures/FIG2.png)

## Installation

In the future we will be adding QClus to PyPI,  making it available via pip install.

Currently we recommend installing QClus from source. To do so, please follow the instructions below:

Note: In order to use our environment installation script, you need to have conda (Anaconda/miniconda) installed on your machine. 

1. Clone this repository into a suitable location on your machine or server using the following command:

    ```git clone https://github.com/linnalab/qclus.git```
    
2. In the root directory of the package, run the following command to create an environment named ```qclus``` and install the required packages:

    ```./environment.sh```


## Getting Started

You can find tutorials on how to use QClus in the `tutorials` directory. They are written in Jupyter notebooks. 

In order to run QClus, you will need the 10X count matrix of your snRNA-seq data, as well as the unspliced values for each cell.

##### Case 1: You have the unspliced values already

Great! Move directly to the qclus_tutorial.ipynb notebook.

##### Case 2: You don't have the unspliced values, but you have run Velocyto on your data and have the .loom file

We have a tutorial for you! Move to the splicing_from_loompy.ipynb notebook, which will show you how to get the unspliced values from the .loom file.

##### Case 3: You don't have the unspliced values, and you haven't run Velocyto on your data

Not a problem! We have implemented our own method for calculating unspliced fraction directly from your 10X bam files! Move to the splicing_from_bam.ipynb notebook, which will show you how to get the unspliced values from the bam files.

## If You Have Issues

If you have any issues with the installation or running the tutorials, please open an issue on this repository. We will do our best to help you out as soon as possible!