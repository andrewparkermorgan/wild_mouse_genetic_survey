# Overview
This repository contains code pertaining to analysis of population structure and inbreeding in wild house mice, as described in this manuscript:

> Morgan AP, Hughes JJ et al. Population structure and inbreeding in wild house mice (Mus musculus) at different geographic scales. bioRxiv. [doi:10.1101/2022.02.17.478179](https://www.biorxiv.org/content/10.1101/2022.02.17.478179)

Large data files (eg. genotype matrix in VCF format) and large intermediate files are not included in this repository. As such the analysis scripts will not be ready to run "out of the box"; file paths will need to be adjusted according to the interested user's local directory structure. But all information needed to verify the conduct of key analyses in the manuscript is available in the code as written.

# Dependencies

For alignment and SNV calling from WGS data:

* bwa
* picard
* samblaster
* samtools
* bcftools
* GATK v.4.1.0.0
* python >= 3.7
* snakemake (any recent version)
See [this repository](https://github.com/IDEELResearch/NGS_Align_QC_Pipelines/tree/master/wgs_pe_improved) for detailed walkthrough of alignment pipeline, and [this repository](https://github.com/IDEELResearch/NGS_Align_QC_Pipelines/tree/master/call_variants_gatk4) for SNV calling.

For post-processing and analyses:

* python >= 3.7
* R v.4.1.1, with these packages plus dependencies therein:
	* tidyverse
	* ggplot2
	* cowplot
	* viridis
	* ggbeeswarm
	* ape
	* BEDMatrix
	* raster
	* rgeos
	* maps
	* [mouser](https://github.com/andrewparkermorgan/mouser) (non-CRAN)
	* [popcorn](https://github.com/andrewparkermorgan/popcorn) (non-CRAN)
* bcftools v.1.9
* plink v.200a3
* akt [v.3beb346](https://github.com/Illumina/akt/tree/3beb3461cf6cc906855c547487dfe3183f52bc5c) (NB: later versions introduced breaking changes to pedigree inference, so version matters)
* TreeMix v.1.12
* ADMIXTURE v.1.3.0

# Files

* Batch calculations, in `bash_scripts`
	* `calc_kinship.sh`: estimate kinship coefficients within taxa, generate lists of putatively unrelated individuals
	* `calc_kinship_subsets.sh`: same as above, but randomize input so that different sets of unrelateds are retained (for robustness checks)
	* `calc_freqs.sh`: calculate allele frequencies within taxa, using unrelateds only
	* `calc_roh.sh`: identify long runs of homozygosity (ROH) using taxon-specific allele frequencies as input
* Python utilities
	* `pyadmix`: wrapper script for ADMIXTURE
* Data management
	* `make_sample_lists.R`: using master sample table and list of putatively unrelated individuals, make text files with lists of sample IDs by taxon, etc. Includes code for creating random subsets of representative individuals to use for robustness checks
	* `make_sample_manifest.R`: make Table S1
* Figures
	* `draw_maps.R`: Figure 1
	* `draw_PCA.R`: Figures 2, 3
	* `draw_admix_plots.R`: Figure 4
	* `analyze_inbreeding.R`: Figure 6
	* `analyze_ROH.R`: Figure 7
	* `analyze_demes.R`: Figure 8
