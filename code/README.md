# Code

The code used for analysis is contained in the .Rmd documents in this folder. Some documents must be run in a particular order, which is defined by file prefix. Other analyses simply require number-prefixed files to be run first.

### Notebook descriptions

| File                             | Description                                                                                                                                                                                                         |
|----------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1.preprocess.Rmd                 | Downloads reads from SRA, trims adapters and primers, imports into QIIME 2, and runs initial QIIME 2 steps like DADA2 denoising, reference tree insertion with SEPP, taxonomy assignment, and contaminant alignment |
| 2.decontam\_by\_correlation.Rmd  | Runs the decontamination procedure, including *decontam* R package and decontamination by correlation                                                                                                               |
| 3.run\_analysis.Rmd              | Runs QIIME2 diversity analysis and outputs visualization files                                                                                                                                                      |
| 4.emp\_comparison.Rmd            | Combines concrete data with EMP & tapwater datasets and runs the analysis; generates some of the figures & tables in the manuscript                                                                                 |
| indicator\_species\_analysis.Rmd | Runs indicator species analysis on the concrete & precursor data                                                                                                                                                    |
| make\_figures.Rmd                | Generates most figures and tables in the manuscript, some of which were later cleaned up in Illustrator                                                                                                             |

### Notebooks

#### Setup

Running this analysis requires R. The following python tools are also required:

-   [Qiime2](https://docs.qiime2.org/2020.2/install/)
-   [SourceTracker2](https://github.com/biota/sourcetracker2)

#### Running

Several steps are quite computationally intensive and use of an HPC is strongly suggested. In some cases the bash chunks for these steps are formatted for submitting over SSH to a compute cluster. These chunks will have to be modified either with alternate SSH credentials or to run locally. In addition, results will have to be copied back for some subsequent chunks to work.
