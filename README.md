# Bacterial communities in concrete change during weathering

This repository includes scripts to process the amplicon data and generate the figures used in (Kiledal, Keffer, and Maresca 2020). An overview of the analysis is provided [here]().

## Notes to the reader

While the full raw data is not provided in this repository, we provide the BIOM tables, sample information, and diversity tables. The raw sequencing data from which the BIOM table was derived is available at [NCBI SRA PRJNA629592](https://www.ncbi.nlm.nih.gov/sra/PRJNA629592).

### Repository Structure

* code

This directory contains the code (as notebooks) for each step; detailed descriptions are found in that directory's README file.

* data

This directory contains the input data, along with output from some computationally expensive tools. 

* results

Figures and tables from the paper which can be reproduced by running the notebooks in the `code` directory.

### Notebooks

#### Setup

Running this analysis requires R. The following python tools are also required:

* [Qiime2](https://docs.qiime2.org/2020.2/install/)
* [SourceTracker2](https://github.com/biota/sourcetracker2)
* [fastspar](https://github.com/scwatts/fastspar)

#### Running

Several steps are quite computationally expensive and use of an HPC is strongly suggested.