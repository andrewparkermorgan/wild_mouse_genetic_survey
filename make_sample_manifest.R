# make_sample_mainfest.R
setwd("~/Dropbox/pmdvlab/wild_redux/")

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggtree)

devtools::load_all("~/Dropbox/pmdvlab/mouser")
devtools::load_all("~/Dropbox/pmdvlab/popcorn")
source("helpers.R")
rb_table <- readRDS("rb_table.rds")

## read sample metadata
samples <- readxl::read_xlsx("samples.xlsx")
keeps <- read_tsv("geno/good_samples.txt", col_names = "iid")
samples <- left_join(keeps, samples)
dim(samples) # should be 814

## representative sample table
pca_reps <- read_fam("geno/isec/merged.renamed.by_pop.fam") %>%
	select(iid, fid) %>%
	rename(pop = fid) %>%
	left_join(samples)
downsample <- read_tsv("sample_lists/downsample.txt", col_names = "iid")
pca_reps$keep <- pca_reps$iid %in% downsample$iid

## continental-level subgroups in domesticus
subgroups <- read_tsv("sample_lists/dom_subgroups.txt", col_names = c("iid","cluster"))

# keep important columms
keep_cols <- c("iid","taxon","sex","dipnum","race","country","locale","site","lat","long")
fsamples <- samples[ ,keep_cols ]
# add country-level labels
fsamples <- left_join(fsamples, select(pca_reps,iid,pop,keep)) %>%
	rename(pca_label = pop, is_core = keep)
# add continental-level clusters for domesticus
fsamples <- left_join(fsamples, subgroups) %>%
	mutate(cluster = subgroup_cluster_labeller(cluster))
# what was source of genotype data?
geno_src <- list(
	wgs = read_tsv("geno/raw.vcf.gz.samples", col_names = "iid"),
	mm = read_tsv("geno/mm_gold.mm10.fixed.vcf.gz.samples", col_names = "iid"),
	giga = read_tsv("geno/nz_giga.gold.fixed.vcf.gz.samples", col_names = "iid")
) %>%
	bind_rows(.id = "platform")
excl <- read_tsv("geno/exclusions.txt", col_names = "iid")
id_matcher <- read_tsv("geno/id_matcher.txt", col_names = c("iid","new_iid"))
geno_src$exclude <- geno_src$iid %in% excl$iid
geno_src <- subset(geno_src, !exclude) %>%
	left_join(id_matcher) %>%
	mutate(final_iid = case_when(!is.na(new_iid) ~ new_iid,
								 TRUE ~ iid))
fsamples <- left_join(fsamples, select(geno_src, final_iid, platform), by = c("iid" = "final_iid"))
# cheat (looking ahead to future analyses) and add inbreeding coefs
fsamples <- left_join(fsamples,
					  read_csv("roh/inbreeding_coefs_A.csv", col_names = c("iid","f_A")))
fsamples <- left_join(fsamples,
					  read_csv("roh/inbreeding_coefs_X.csv", col_names = c("iid","f_X")))
fsamples$f_A[ is.na(fsamples$f_A) ] <- 0.0
fsamples$f_X[ is.na(fsamples$f_X) ] <- 0.0
fsamples$f_X[ fsamples$sex != "F" ] <- NA_real_
fsamples$f_A[ !(fsamples$taxon %in% c("dom","mus","cas")) ] <- NA_real_
fsamples$f_X[ !(fsamples$taxon %in% c("dom","mus","cas")) ] <- NA_real_
# write result
arrange(fsamples, taxon, country, locale, iid) %>%
	write_tsv(., "sample_lists/final_sample_table.tsv", col_names = TRUE)
