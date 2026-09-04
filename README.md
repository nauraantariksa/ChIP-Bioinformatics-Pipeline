# ChIP Bioinformatics Pipeline

This GitHub page is a collection of scripts and general code used for the analysis of ChIP datasets

Make sure these dependencies are installed in your system:<br /> 

- anaconda3<br /> 
- python (version 3.9)<br />
- fastp v0.23.4
- FastQC v0.11.9
- Bowtie2 v2.5.1   
- SAMtools v1.22.1 <br />
- Picard v2.25.1 
- BEDTools v2.31.1 <br /> 

In anaconda3, make separate environments for the following Python suites used for the data analyses:<br /> 

- deepTools v3.5.6 (in nfa23 Imperial HPC, in deeptools_env, perform $conda install -c conda-forge -c bioconda deeptools)<br /> 
- MACS2 v2.2.9.1 (in nfa23 Imperial HPC, in macs2_env, perform $conda install -c bioconda macs2)<br /> 

For other analyses, install or clone Github repositories:<br /> 

- R<br />
- HOMER v5.1 <br />

For visualization:

- IGV v2.19.8
