# confirm_sex.R
setwd("~/Dropbox/pmdvlab/wild_redux/")

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

devtools::load_all("~/Dropbox/pmdvlab/mouser")
devtools::load_all("~/Dropbox/pmdvlab/popcorn")
rb_table <- readRDS("rb_table.rds")
source("./helpers.R")

samples <- readxl::read_xlsx("samples.xlsx")
keeps <- read_tsv("geno/good_samples.txt", col_names = "iid")
samples <- left_join(keeps, samples)
samples$karyo <- with(samples, ifelse(dipnum < 40, "Rb", "standard"))
samples$taxon_old <- samples$taxon
samples$taxon <- factor_mus(samples$taxon)
dim(samples) # should be 814

geno <- BEDMatrix::BEDMatrix("geno/isec/merged.renamed.XY.bed")
fam <- read_fam("geno/isec/merged.renamed.XY.fam")
bim <- read_map("geno/isec/merged.renamed.XY.bim")

geno <- as.matrix(geno)
rownames(geno) <- fam$iid
colnames(geno) <- bim$marker

geno_x <- geno[ ,bim$chr == "chrX" & bim$pos < pseudoautosomal_boundary("mm10")[2], drop = FALSE ]
geno_y <- geno[ ,bim$chr == "chrY", drop = FALSE ]

y_good <- rowSums(geno_y == 0 | geno_y == 2, na.rm = TRUE)
x_het <- rowSums(geno_x == 1, na.rm = TRUE)
xy <- tibble(iid = rownames(geno),
			 x_het = unname(x_het), y_good = unname(y_good))
xy <- left_join(xy, select(samples, iid, sex))

ggplot(xy) +
	geom_point(aes(x = x_het, y = y_good, colour = sex)) +
	scale_x_continuous("heterozygous calls on X") +
	scale_y_continuous("non-missing calls on Y") +
	scale_color_sex() +
	theme_classic2() +
	legend_inside("topright")

ggsave("figures/sex_qc_merged.pdf", width = 3, height = 3)
