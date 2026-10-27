# The BAM counter against velocyto on ENCODE heart samples


We compared `qclus splicing-from-bam` with velocyto on 75 human heart snRNA-seq experiments from ENCODE4. With the GTF file and the repeat mask, the counter gives velocyto's spliced, unspliced and ambiguous counts for every barcode of every experiment, and QClus gives every barcode the same label with either fraction. Without the GTF file, the fraction is higher than velocyto's and QClus passes more nuclei.



## Data


- **Experiments:** 75 ENCODE4 human heart snRNA-seq experiments (10x Multiome gene expression). The tables give their accessions.

- **Alignment:** Cell Ranger 8.0.1 (`cellranger count --chemistry ARC-v1`) from the FASTQ files, with the reference `refdata-gex-GRCh38-2020-A`.

- **Barcodes:** those Cell Ranger called (`filtered_feature_bc_matrix`), 951 to 41,458 per experiment and 1,306,785 in all.

- **Repeat mask:** RepeatMasker for hg38 from UCSC (`rmsk.txt.gz`), written as GTF, as in velocyto's documentation.



## Commands


velocyto 0.17.17, with the repeats masked:

```bash
velocyto run10x -m hg38_rmsk.gtf --samtools-threads 8 --samtools-memory 2000 \
    <cellranger output folder> refdata-gex-GRCh38-2020-A/genes/genes.gtf
```


The fraction of each barcode is unspliced / (spliced + unspliced + ambiguous) from the loom file, read with `qclus.utils.fraction_unspliced_from_loom`.


The counter with velocyto's rules, from QClus commit 10583f1, whose counting code is that of QClus 0.4.0:

```bash
qclus splicing-from-bam --bam outs/possorted_genome_bam.bam \
    --barcodes outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
    --gtf refdata-gex-GRCh38-2020-A/genes/genes.gtf --repeat-mask hg38_rmsk.gtf \
    --cores 8 --output <experiment>.csv.gz
```


The counter without a GTF file is QClus 0.3.0 with the same command, without `--gtf` and `--repeat-mask`.


QClus 0.3.0 then ran with its default settings on each fraction, and the labels were compared barcode by barcode. velocyto and the counter had 8 cores per job on a Slurm cluster. The scripts that ran the jobs are specific to that cluster and are not included.



## Results


**With the GTF file and the repeat mask**, the counts are identical to velocyto's for every barcode in all 75 experiments: 1,537,735,103 spliced, 2,300,393,927 unspliced and 360,548,586 ambiguous molecules. QClus gives all 1,306,785 barcodes the same label with either fraction, and 436,587 pass with both.


**Without the GTF file**, the fraction is a median of 0.058 higher than velocyto's (per experiment, from −0.015 to 0.120), with a Pearson r of 0.954 to 0.995. QClus passes 491,135 nuclei instead of 436,587, a median of 10% more per experiment, and the Jaccard index of the nuclei passing is 0.70 to 0.97 (median 0.88).


Median compute per experiment:

| Method          | Wall time  | CPU time | Peak memory |
|-----------------|------------|----------|-------------|
| velocyto        | 2 h 12 min | 3.3 h    | 35.6 GB     |
| Counter, GTF    | 11 min     | 68 min   | 4.1 GB      |
| Counter, no GTF | 6 min      | 18 min   | 1.4 GB      |



## Files


- **`counter_with_gtf.csv`:** one row per experiment. It gives the barcodes, velocyto's and the counter's total spliced, unspliced and ambiguous molecules, the percentage of barcodes with all three counts the same, the Pearson r and the mean absolute difference of the fractions, the nuclei passing QClus with each fraction, and the percentage of barcodes with the same label.

- **`counter_without_gtf.csv`:** one row per experiment for the counter without a GTF file. It gives the barcodes, the Pearson r, the median fractions, the median and the mean absolute difference per barcode (counter minus velocyto), the nuclei passing with each fraction, the Jaccard index of the nuclei passing, and the percentage of barcodes with the same label.

- **`compute.csv`:** one row per experiment and method, with the cores, wall time and CPU time in seconds, the peak memory in GB as Slurm measured it, and the number of reads Cell Ranger reports for the experiment.
