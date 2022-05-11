setwd("~/Dropbox/pmdvlab/wild_redux/")

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

devtools::load_all("~/Dropbox/pmdvlab/mouser")
devtools::load_all("~/Dropbox/pmdvlab/popcorn")
source("helpers.R")

samples <- readxl::read_xlsx("samples.xlsx")
keeps <- read_tsv("geno/good_samples.txt", col_names = "iid")
samples <- left_join(keeps, samples)
samples$taxon <- factor_mus(samples$taxon)
dim(samples) # should be 814

fnames <- expand_grid(taxon = c("dom","mus","cas","gen"), iter = seq(0,9)) %>%
	mutate(fname = paste0("freq/random_subsets/", taxon, "_", iter, ".txt.gz"))

read_freq <- function(path, ...) {
	read_tsv(path, col_names = c("chr","pos","alleles","freq"), col_types = "cicd")
}

freqs_ss <- fnames %>%
	group_by(taxon, iter) %>%
	do(read_freq(.$fname))

freqs_ss$chr <- factor_chrom(freqs_ss$chr)
freqs_ss$chrtype <- factor_chromtype(freqs_ss$chr)

freqs_summ <- freqs_ss %>%
	group_by(taxon, chrtype, chr, pos) %>%
	mutate(maf = pmin(freq, 1 - freq)) %>%
	summarise(mu = mean(maf), sd = sd(maf), fixed = all(freq == 0 | freq == 1))
	
subset(freqs_summ, taxon == "dom" & !fixed & chrtype == "A") %>%
	ggplot() +
	geom_point(aes(x = mu, y = sd/mu))


## read population labels
pca_reps <- read_fam("geno/isec/merged.renamed.by_pop.fam") %>%
	select(iid, fid) %>%
	rename(pop = fid) %>%
	left_join(samples)

## read continent-level cluster assignments
dom_cl <- read_tsv("sample_lists/dom_subgroups.txt", col_names = c("iid","cluster"), col_types = "cc")

## target VCF file
GOOD_VCF <- "geno/isec/merged.renamed.G90.MAF01.vcf.gz"

## iterate over random unrelated subsets and do PCA
for (ii in 0:9) {
	
	message("iteration ", ii, " ...")
	
	## keep up to 10 individuals per location, permuting the order
	downsample <- read_tsv(paste0("kinship/random_subsets/dom_", ii,".txt"), col_names = "iid")
	these_reps <- left_join(downsample, pca_reps) %>%
		group_by(taxon, pop) %>%
		mutate(rank = sample( seq_along(iid)) )%>%
		mutate(keep = case_when(taxon != "dom" ~ TRUE,
								rank <= 10 ~ TRUE,
								TRUE ~ FALSE))
	subset()
	
	## assign population symbols/colours
	pop_list <- unique(subset(these_reps, taxon == "dom" & keep)$pop)
	POP_COL <- setNames( repeat_colours(pop_list, "Paired"), pop_list )
	POP_COL <- c(POP_COL, "(other)" = "grey92")
	POP_SYM <- setNames( repeat_shapes(pop_list), pop_list )
	POP_SYM <- c(POP_SYM, "(other)" = 19)
	
	## run PCA
	this_pc <- pca_vcf(GOOD_VCF, K = 10,
					   skeleton = subset(these_reps, taxon == "dom" & keep)$iid,
					   samples = subset(pca_reps, taxon == "dom")$iid)
	write_pca_result(this_pc, paste0("pca/dom_subsampled_", ii, ".pca"))
	
}

## now iterate again to make the plots
plots_1 <- vector("list", 10)
plots_2 <- vector("list", 10)
for (ii in 0:9) {
	
	## read PCA result
	this_pc <- read_pca_result(paste0("pca/dom_subsampled_", ii, ".pca.pca"))
	
	## draw plot
	pp <- merge_pca_result(this_pc, dom_cl) %>%
		ggplot() +
		geom_point(aes(x = PC1, y = PC2, fill = cluster), shape = 21,
				   colour = "grey50") +
		scale_fill_manual(values = SUBGROUP_COLS, labels = subgroup_cluster_labeller) +
		pca_axes(this_pc) +
		guides(fill = "none") +
		#ggtitle("PCA on autosomal genotypes (downsampled, domesticus only)") +
		theme_classic2()
	plots_1[[ ii+1 ]] <- pp
	
	pp <- merge_pca_result(this_pc, dom_cl) %>%
		ggplot() +
		geom_point(aes(x = PC2, y = PC3, fill = cluster), shape = 21,
				   colour = "grey50") +
		scale_fill_manual(values = SUBGROUP_COLS, labels = subgroup_cluster_labeller) +
		pca_axes(this_pc, K = 2:3) +
		guides(fill = "none") +
		#ggtitle("PCA on autosomal genotypes (downsampled, domesticus only)") +
		theme_classic2()
	plots_2[[ ii+1 ]] <- pp
	
}

cowplot::plot_grid(plotlist = plots_1, nrow = 5, ncol = 2)
ggsave("figures/pca_dom_subsample_reps_10up.pdf", width = 8, height = 12)

cowplot::plot_grid(plotlist = plots_2, nrow = 5, ncol = 2)

