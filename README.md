# Code and study data for Dos Santos et al. (2024) 

Last updated: 16th May 2024

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

- **paper_initial:**
	- Directory exported from Overleaf, containing all files used for compilation of the manuscript .tex file as well as the resulting manuscript .pdf file. This is the INITIAL SUBMISSION to Microbiome.

- **paper_revised:**
	- Directory exported from Overleaf, containing all files used for compilation of the manuscript .tex file as well as the resulting manuscript .pdf file and responses to all reviewer comments. This is the REVISED SUBMISSION to Microbiome.

- **RData:**
	- Directory containing .Rda /.txt /.tsv /.csv files imported into R for use in the analysis and/or visualisation of study data.
	- **virginia_output:**
		- Directory containing .txt files containing information on quality control and mapping of raw FASTQ files to various databases during data processing, as well as an R script for summarising these files.

- **Supplement:**
	- Directory containing supplementary information on the changes made to the KEGG orthology terms (i.e. pathway names) and justifications for these changes, an Excel file containing strain names and accession numbers for all _Gardnerella_ genomes used in the pangenome analysis,  and the old code and input files used to determine _Gardnerella_ marker genes (i.e. corresponding to the initial submission to Microbiome).

<br>

### Figure code
All R scripts are located in the `code/` directory. All figures were edited in Inkscape in some way (details given in `figs_for_paper/multipanel_fig_edits.txt`)

##### <u>Main text figures</u>

- **Figure 1:**
  - Initially created with BioRender and edited in Inkscape.

- **Figure 2:**
  - A: `lon_eur_species_heatmap.R`
  - B: `lon_eur_species_biplots.R`
  - C: `lon_eur_health_vs_BV.R`
  - D: Initially created with BioRender and edited in Inkscape.

- **Figure 3:**
  - Top: `lon_eur_bv_subgroups.R`
  - Bottom: `revisions_stackedBars_redo.R` (revision); `lon_eur_bvSubgroups_stackedBars.R` (initial submission)
  
- **Figure 4:**
  - A: `preg_vs_nonpreg.R`
  - B: `gardnerella_pca_redo` (revision); `preg_vs_nonpreg_gardnerella_biplots.R` (initial submission)

- **Figure 5:**
  - A: `virginia_validating.R`
  - B: `virginia_validating.R`
  - C: `virginia_validating.R`
  
- **Figure 6:**
  - Initially created with BioRender and edited in Inkscape.

<br>

##### <u>Supplementary figures</u>
Naming of the supplementary figure files follows a numerical order. The file name corresponding to each supplementary figure in the manuscript is indicated below, along with the R script used to produce it.

- **Suppl. Fig. 1:**
  -  `revisions_scale_plots.R` (revision); `scale_vs_noScale_MWplots.R` (initial submission)

- **Suppl. Fig. 2:**
  - `lon_eur_health_vs_BV.R`

- **Suppl. Fig. 3:**
  - `revisions_bvSubgroups_noBVAB1.R`

- **Suppl. Fig. 4:**
  - `virginia_definingMolecularBV.R` (initial submission)

- **Suppl. Fig. 5:**
  -  `gardnerella_pca_redo.R` (A \& B, revision); `lon_eur_gardnerella_biplots.R` (A, initial submission) and `virginia_gardnerella_biplots.R` (B, initial submission)

- **Suppl. Fig. 6:**
  - `virginia_definingMolecularBV.R`

- **Suppl. Fig. 7:**
  - `virginia_validating.R`

- **Suppl. Fig. 8:**
  - `virginia_coda4microbiome.R`

- **Suppl. Fig. 9:**
  - `revisions_stackedBars_redo.R` (revision); `virginia_bvSubgroups_stackedBars.R` (initial submission)

- **Suppl. Fig. 10:**
  - `virginia_validating_eggNOG.R`
  
  - **Suppl. Fig. 11:**
  - `virginia_validating.R`
  