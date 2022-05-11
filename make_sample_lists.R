setwd("~/Dropbox/pmdvlab/wild_redux/")

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggtree)

devtools::load_all("~/Dropbox/pmdvlab/mouser")
devtools::load_all("~/Dropbox/pmdvlab/popcorn")
source("helpers.R")

## read sample metadata
samples <- readxl::read_xlsx("samples.xlsx")
keeps <- read_tsv("geno/good_samples.txt", col_names = "iid")
samples <- left_join(keeps, samples)
dim(samples) # should be 814

## urelateds
unrel <- read_tsv("sample_lists/unrel.txt", col_names = "iid")

## create representative sample set for PCA
## idea is to equalize geographic representation as much as possible, since some areas (ie Madeira, Maryland) are oversampled
pca_pops <- samples %>% select(iid, country, region, locale, site, race, landmass, taxon)
pca_pops$pop <- with(pca_pops,
					 case_when(country == "PT" & region == "Madeira" ~ "Madeira",
					 		  country == "PT" & region == "Porto Santo" ~ "Porto Santo",
					 		  region == "Gough Island" ~ "Gough",
					 		  region == "Floreana Island" ~ "Floreana",
					 		  region == "SE Farallon Island" ~ "Farallon",
					 		  region == "Orkney Islands" ~ "Orkney",
					 		  region == "Scotland" ~ "Scotland",
					 		  landmass == "island" & country == "NZ" ~ site,
					 		  country == "US" & region == "MD" ~ "Maryland",
					 		  country == "US" & region != "MD" ~ region,
					 		  country == "GR" & region == "Crete" ~ "Crete",
					 		  country == "NZ" & landmass != "island" ~ "New Zealand",
					 		  country == "DE" & locale == "HGL" ~ "Heligoland",
					 		  taxon != "dom" ~ taxon,
					 		  TRUE ~ country_decoder[country]
					 ))

## keep up to 10 individuals per location, arbitrarily
pca_reps <- group_by(pca_pops, taxon, pop) %>%
	subset(iid %in% unrel$iid) %>%
	mutate(rank = seq_along(iid)) %>%
	mutate(keep = case_when(taxon != "dom" ~ TRUE,
							rank <= 10 ~ TRUE,
							TRUE ~ FALSE))

## write downsampled list
ungroup(pca_reps) %>%
	subset(keep) %>%
	select(iid) %>%
	write_tsv("sample_lists/downsample.txt", col_names = FALSE)

## write downsampled domesticus iids
ungroup(pca_reps) %>%
	subset(keep & taxon == "dom") %>%
	select(iid) %>%
	write_tsv(., "sample_lists/downsample_dom.txt", col_names = FALSE)

## now downsample again, but in randomly chosen unrelated sets for robustness checks
for (ii in 0:9) {
	
	## keep up to 10 individuals per location, permuting the order
	downsample <- read_tsv(paste0("kinship/random_subsets/dom_", ii,".txt"), col_names = "iid")
	these_reps <- left_join(downsample, pca_pops) %>%
		group_by(taxon, pop) %>%
		mutate(rank = sample( seq_along(iid)) )%>%
		mutate(keep = case_when(taxon != "dom" ~ TRUE,
								rank <= 10 ~ TRUE,
								TRUE ~ FALSE))
	ungroup(these_reps) %>%
		subset(taxon == "dom" & keep) %>%
		select(iid) %>%
		arrange(iid) %>%
		write_tsv( paste0("kinship/random_subsets/dom_", ii, "_downsample.txt"), col_names = FALSE )
}


## write some plink *.fam files for later use
fam <- read_fam("geno/isec/merged.renamed.fam")
fam$fid <- pca_pops$pop[ match(fam$iid, pca_pops$iid) ]
write_tsv(fam, "geno/isec/merged.renamed.by_pop.fam", col_names = FALSE)
fam$fid <- samples$taxon[ match(fam$iid, samples$iid) ]
write_tsv(fam, "geno/isec/merged.renamed.by_taxon.fam", col_names = FALSE)

