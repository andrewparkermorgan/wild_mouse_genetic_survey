# check_subsampled_admix.R
setwd("~/Dropbox/pmdvlab/wild_redux/")

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

devtools::load_all("~/Dropbox/pmdvlab/mouser")
devtools::load_all("~/Dropbox/pmdvlab/popcorn")
source("./helpers.R")

## read sample metadata
samples <- readxl::read_xlsx("samples.xlsx")
keeps <- read_tsv("geno/good_samples.txt", col_names = "iid")
samples <- left_join(keeps, samples)
dim(samples) # should be 814

## read population assignments
fam <- read_fam("geno/isec/merged_dom.A.fam")
pops <- read_fam("geno/isec/merged.renamed.by_pop.fam") %>%
	select(fid, iid) %>%
	rename(pop = fid)

## read main admixture result
Qa <- read_Q_matrix("admix/merged_dom.A.2.Q", iids = fam$iid)
so <- sort_by_cluster(Qa)
Qa <- tidy(Qa) %>% 
	left_join(pops)
Qa2 <- Qa
so2 <- so

## plot actual result
plot_admixture(Qa, sort_order = so) +
	facet_grid(. ~ pop, space = "free_x", scale = "free_x") +
	scale_fill_brewer(palette = "Paired") +
	guides(fill = "none")

## now loop on subsampled runs
plots <- vector("list", 10)
for (ii in 0:4) {
	
	ff <- paste0("scratch/random_subsets/proj_", ii, ".A.2.Q")
	Q <- read_Q_matrix(ff, iids = fam$iid)
	QQ <- tidy(Q) %>% left_join(pops)
	p0 <- plot_admixture(Qa, sort_order = so) +
		facet_grid(. ~ pop, space = "free_x", scale = "free_x") +
		scale_fill_brewer(palette = "Paired") +
		guides(fill = "none")
	
	plots[[ ii+1 ]] <- p0
	
}

cowplot::plot_grid(plotlist = plots[1:5], ncol = 1)
ggsave("figures/admix_dom_random_subsets_K2.pdf", width = 12, height = 18)
