# Bacterial communities in concrete change during weathering

This repository includes scripts to process the amplicon data and generate the figures used in (Kiledal, Keffer, and Maresca 2020).

## Notes to the reader

The raw sequence data is housed at [NCBI SRA PRJNA629592](https://www.ncbi.nlm.nih.gov/sra/PRJNA629592). Scripts include functionality to download this data, and run our entire analysis pipeline.

### Repository Structure

-   code

This directory contains the code (as notebooks) for each step; detailed descriptions are found in that directory's README file.

-   data

This directory is where processed data is stored, and comes populated with items like sample metadata and several other reference files.

-   results

Figures and tables from the paper which can be reproduced by running the notebooks in the `code` directory.

### Notebooks

#### Setup

Running this analysis requires R. The following python tools are also required:

-   [Qiime2](https://docs.qiime2.org/2020.2/install/)
-   [SourceTracker2](https://github.com/biota/sourcetracker2)

#### Running

Several steps are quite computationally expensive and use of an HPC is strongly suggested.
