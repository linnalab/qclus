# QClus

This is a novel nuclei filtering method targeted to streamline the processing of snRNA-seq samples. The algorithm uses various metrics to cluster nuclei and filter empty and highly contaminated droplets. This approach combined with other filtering steps enables for flexible, automated, and reliable preprocessing of samples with varying number of nuclei, quality, and contamination levels. The robustness of this method has been validated on a large number of heterogeneous datasets, in terms of number of nuclei, overall quality, and contamination.

Additionally, while the method was originally developed for cardiac snRNA-seq data, in our paper we show that given it's felxible design it can also be applied to other tissues. We will provide detailed tutorials for this process.

Any and all comments/criticisms/suggestions are enthusiastically received! :-)

## Installation

In the future we will be adding QClus to PyPI,  making it available via pip install.

Currently we recommend installing QClus from source. To do so, please follow the instructions below:

Note: In order to use our environment installation script, you need to have conda ([Anaconda](https://docs.anaconda.com/anaconda/install/)/[Miniconda](https://docs.anaconda.com/miniconda/install/)) installed on your machine. 

1. Clone this repository into a suitable location on your machine or server using the following command:

    ```git clone https://github.com/linnalab/qclus.git```
    
2. In the root directory of the package, run the following command to create an environment named ```qclus``` and install the required packages:

    ```./environment.sh```


## Getting Started

You can find tutorials on how to use QClus in the `tutorials` directory. They are written in Jupyter notebooks. 

In order to run QClus, you will need the 10X count matrix of your snRNA-seq data, as well as the unspliced values for each cell.

#### Case 1: You have the unspliced values already

Great! Move directly to the qclus_tutorial.ipynb notebook.

#### Case 2: You don't have the unspliced values, but you have run Velocyto on your data and have the .loom file

We have a tutorial for you! Move to the splicing_from_loompy.ipynb notebook, which will show you how to get the unspliced values from the .loom file.

#### Case 3: You don't have the unspliced values, and you haven't run Velocyto on your data

Not a problem! We have implemented our own method for calculating unspliced fraction directly from your 10X bam files! Move to the splicing_from_bam.ipynb notebook, which will show you how to get the unspliced values from the bam files.

## How to cite

We hope you find the QClus package useful for your research! If you do, please remember to cite [our paper in Nucleic Acids Research](https://doi.org/10.1093/nar/gkae1145):

Eloi Schmauch, Johannes Ojanen, Kyriakitsa Galani, Juho Jalkanen, Kristiina Harju, Maija Hollmén, Hannu Kokki, Jarmo Gunn, Jari Halonen, Juha Hartikainen, Tuomas Kiviniemi, Pasi Tavi, Minna U Kaikkonen, Manolis Kellis, Suvi Linna-Kuosmanen, QClus: a droplet filtering algorithm for enhanced snRNA-seq data quality in challenging samples, Nucleic Acids Research, Volume 53, Issue 1, 13 January 2025, gkae1145, https://doi.org/10.1093/nar/gkae1145

You can also find the latest preprint for this method [here](https://www.biorxiv.org/content/10.1101/2022.10.21.513315v2).

## If You Have Issues

If you have any issues with the installation or running the tutorials, please open an issue on this repository. We will do our best to help you out as soon as possible!


