# QClus: Robust and reliable preprocessing method for human heart snRNA-seq

This is a novel nuclei filtering method targeted to human heart samples. We use specific metrics such as splicing, mitochondrial gene expression, nuclear gene expression, and non-cardiomyocyte and cardiomyocyte marker gene expression to cluster nuclei and filter empty and highly contaminated droplets. This approach combined with other filtering steps enables for flexible, automated, and reliable cleaning of samples with varying number of nuclei, quality, and contamination levels. The robustness of this method has been validated on a large number of heterogeneous datasets, in terms of number of nuclei, overall quality, and contamination. 

You can find the original preprint for this method at the link below:

https://www.biorxiv.org/content/10.1101/2022.10.21.513315v1

Any and all comments/criticism/suggestions enthusiastically received! :-)

## Required packages

- numpy
- pandas
- scrublet    
- loompy
- scanpy
- anndata
- sklearn


## Installation

1. Clone this repository into a suitable location on your machine or server using the following command:

```
    git clone https://github.com/johannesojanen/qclus.git
```
    
2. In the root directory of the package, run the following command to install the required packages:

    ```pip install .```

## Quickstart guide

Below is a quick workflow for QClus. 

For a more thorough tutorial please see the tutorials/tutorial.ipynb notebook.

```python
    #import library
    import qclus

    #Define the path to an .h5 file or faw counts, as well as the path to a .loom file containing the splicing information for those same counts
    counts_path = "filtered_feature_bc_matrix_739.h5"
    loompy_path = "counts_counts_CAD_739.loom"

    #run QClus with default settings
    adata = qclus.run_qclus(counts_path,  
                            loompy_path)
```
