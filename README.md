# Code and study data for Dos Santos et al. (2024) 

Last updated: 24th April 2024.

This repository contains step-by-step walkthroughs and R code to reproduce all data collection, processing, analyses and visualisation carried out by the authors. Below is a description of all sub-directories and their contents, as well as a summary of what R scripts were used to create the main and supplementary figures.

### Contents:
- **1_VIRGO:**
  - Directory containing functional and taxonomic lookup tables (among others) for the VIRGO database (also available from https://virgo.igs.umaryland.edu).

- **code:**
  - Directory containing R scripts and tab-delimited text files used for data analysis and visualisation. R script files associated with figures and supplementary figures are indicated after the directory contents.

- **figs_for_paper:**
	- Directory containing all .png files used to make the figures for the above study, as well as descriptions of edits to figures during assembly of multi-panel figures
	- **biorenderFigures:**
		- Directory for the .png files containing figures made with biorender
		- **publicationLicenses:**
			- Directory containing biorender publication licenses in PDF format
	- **R_FigsLondonEurope:**
		- Directory containing .png files created in R pertaining to the exploratory dataset (London/Europe datasest; Macklaim et al. / Deng et al.)
	- **R_FigsPregNonPreg:**
		- Directory containing .png files created in R pertaining to the combined exploratory and validation datasets
	- **R_FigsVirginia:**
		- Directory containing .png files created in R pertaining to the validation dataset (Virginia; MOMS-PI study, Fettweis et al.)
	- **SVGs:**
		- Directory containing .svg files created in Inkscape for final study figures

- **paper:**
	- Directory exported from Overleaf, containing all files used for compilation of the manuscript .tex file as well as the resulting manuscript .pdf file.

- **RData:**
	- Directory containing .Rda /.txt /.tsv /.csv files imported into R for use in the analysis and/or visualisation of study data (<u>NOTE:</u> several .csv and .Rda files have been intentionally omitted from this repository as the underlying data is locked behind an authorised access request through the NCBI Database of Genotypes and Phenotypes).
	- **virginia_output:**
		- Directory containing .txt files containing information on quality control and mapping of raw FASTQ files to various databases during data processing, as well as an R script for summarising these files.

- **Supplement:**
	- Directory containing supplementary information on the changes made to the KEGG orthology terms (i.e. pathway names) and justifications for these changes.

<br>

### Figure code
All R scripts are located in the `code/` directory. All figures were edited in inkscape in some way (details given in `figs_for_paper/multipanel_fig_edits.txt`)

##### <u>Main text figures</u>

- **Figure 1:**
  - Created with BioRender.

- **Figure 2:**
  - A: `lon_eur_species_heatmap.R`
  - B: `lon_eur_species_biplots.R`
  - C: `lon_eur_health_vs_BV.R`
  - D: Created with BioRender

- **Figure 3:**
  - Top: `lon_eur_bv_subgroups.R`
  - Bottom: `lon_eur_bvSubgroups_stackedBars.R`
  
- **Figure 4:**
  - A: `preg_vs_nonpreg.R`
  - B: `preg_vs_nonpreg_gardnerella_biplots.R`

- **Figure 5:**
  - A: `virginia_validating.R`
  - B: `virginia_validating.R`
  - C: `virginia_validating.R`

<br>

##### <u>Supplementary figures</u>
Naming of the supplementary figure files does not follow a numerical order. The file name corresponding to each supplementary figure in the manuscript is indicated below, along with the R script used to produce it.

- **Suppl. Fig. 1:**
  - `supplFig0.png`
  - `scale_vs_noScale_MWplots.R`

- **Suppl. Fig. 2:**
  - `supplFig1.png`
  - `lon_eur_health_vs_BV.R`

- **Suppl. Fig. 3:**
  - `supplFig2.png`
  - `virginia_definingMolecularBV.R`

- **Suppl. Fig. 4:**
  - `supplFig3.png`
  - `virginia_pretermBirth.R` (table created in Inkscape)

- **Suppl. Fig. 5:**
  - `supplFig4.png`
  - `lon_eur_gardnerella_biplots.R` (A) and `virginia_gardnerella_biplots.R` (B)

- **Suppl. Fig. 6:**
  - `supplFig5.png`
  - `virginia_definingMolecularBV.R`

- **Suppl. Fig. 7:**
  - `supplFig6.png`
  - `virginia_validating.R`

- **Suppl. Fig. 8:**
  - `supplFigvirginia_coda4microbiome_fig.png`
  - `virginia_coda4microbiome.R`

- **Suppl. Fig. 9:**
  - `supplFig7.png`
  - `virginia_bvSubgroups_stackedBars.R`

- **Suppl. Fig. 10:**
  - `supplFig8.png`
  - `virginia_validating.R`
