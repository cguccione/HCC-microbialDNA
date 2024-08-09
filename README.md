# Microbial DNA as a Diagnostic Marker for Hepatocellular Carcinoma

## Data Preprocessing
### Host Depletion
Host deplete all samples with 47 Pangenomes from the Human Pangenome Project, T2T-CHM13v2.0, GRCh38.p14 and PhiX ([Host Depletion Repo](https://github.com/cguccione/human_host_depletion))
- *Input*: Raw fastq files with mixed microbial and human reads
- *Output*: Fastq files with just microbial reads

### Microbial Classification
Run SHOGUN + Woltka pipeline using RS210 
- *Input*: Fastq files with just microbial reads
- *Output*: Biom/taxonomy table for all samples including blanks

### Lab-associated Decontamination
Run [SCRuB](https://www.nature.com/articles/s41587-023-01696-w) to remove any lab associated contamination: [1_SCRuB_Decontamination.ipynb](https://github.com/cguccione/HCC-microbialDNA/blob/main/1_SCRuB_Decontamination.ipynb)
- *Input*: Biom/taxonomy table for all samples including blanks
- *Output*: Decontaminated biom/taxonomy table for all biological samples

## Data Analysis
### Alpha & Beta Diversity Analysis
Alpha & Beta diversity analysis using [Qiime2](https://qiime2.org)
- *Input*: Decontaminated biom/taxonomy table for all biological samples with extraneous reads removed
- *Output*: Alpha & Beta diversity differences across disease groups

### Birdman
Performing differential abundance using RPCA, Genome data across disease groups data using [Birdman](https://birdman.readthedocs.io/en/latest/).
- *Input*: Decontaminated biom/taxonomy table for all biological samples with extraneous reads removed
- *Output*: Top microbial features according to Birdman
