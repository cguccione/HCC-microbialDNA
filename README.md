# Microbial DNA as a Diagnostic Marker for Hepatocellular Carcinoma

**Goal**: To determine if algorithms from TCGA are robust on a new cohort of patients with HCC and distinguish liver tumors (HCC from non-HCC tumors)

**Hypothesis**: Tumor-associated microbial DNA can be used to classify liver tumor type

**Approach**: Use Machine Learning analysis using microbiomes as a classifier to determine cancer vs. non-cancer tissue + cancer type

## Data Preprocessing
### Host Depletion
Host deplete all samples with 47 Pangenomes from the Human Pangenome Project, T2T-CHM13v2.0, GRCh38.p14 and PhiX ([Host Depletion Repo](https://github.com/cguccione/human_host_depletion))
- *Input*: Raw fastq files with mixed microbial and human reads
- *Output*: Fastq files with just microbial reads

### Microbial Classification
Run SHOGUN + Woltka pipeline using RS210 and WoLr2 databases on Qiita (Study [15336](https://qiita.ucsd.edu/study/description/15336#), Prep 16181)
- *Input*: Fastq files with just microbial reads
- *Output*: Biom/taxonomy table for all samples including blanks

### Lab-associated Decontamination
Run [SCRuB](https://www.nature.com/articles/s41587-023-01696-w) to remove any lab associated contamination: [1_SCRuB_Decontamination.ipynb](https://github.com/cguccione/HCC-microbialDNA/blob/main/1_SCRuB_Decontamination.ipynb)
- *Input*: Biom/taxonomy table for all samples including blanks
- *Output*: Decontaminated biom/taxonomy table for all biological samples

### Falsely Mapped Read Removal
[Zebra](https://journals.asm.org/doi/full/10.1128/msystems.00758-22) to remove any reads which may have been incorrectly mapped: [2_Zebra_Coverage.ipynb](https://github.com/cguccione/HCC-microbialDNA/blob/main/2_Zebra_Coverage.ipynb)
- *Input*: Decontaminated biom/taxonomy table for all biological sample
- *Output*: Decontaminated biom/taxonomy table for all biological samples with extraneous reads removed

## Data Analysis
### Beta Diversity Analysis
Beta diversity analysis using [Qiime2](https://qiime2.org)
- *Input*: Decontaminated biom/taxonomy table for all biological samples with extraneous reads removed
- *Output*: Beta diversity differences across disease groups

### Feature Selection (& Machine Learning)
Feature selection using WoLr2, RPCA, Genome across disease groups. Machine learning to split between groups was not very successful due to small group sizes but in script as well. 
- *Input*: Decontaminated biom/taxonomy table for all biological samples with extraneous reads removed
- *Output*: Top microbial features according to Chi2 and Random Forest

### Birdman
Performing differential abundance using WoLr2, RPCA, Genome data across disease groups data using [Birdman](https://birdman.readthedocs.io/en/latest/).
- *Input*: Decontaminated biom/taxonomy table for all biological samples with extraneous reads removed
- *Output*: Top microbial features according to Birdman
