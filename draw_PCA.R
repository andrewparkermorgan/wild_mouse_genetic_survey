# draw_PCA.R
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

## path to VCF file, filtered for missingness and MAF
GOOD_VCF <- "geno/isec/merged.renamed.G90.MAF01.vcf.gz"

## representative sample table
pca_reps <- read_fam("geno/isec/merged.renamed.by_pop.fam") %>%
	select(iid, fid) %>%
	rename(pop = fid) %>%
	left_join(samples)
downsample <- read_tsv("sample_lists/downsample.txt", col_names = "iid")
pca_reps$keep <- pca_reps$iid %in% downsample$iid

## perform PCA including all subspecies
## use downsampled invididuals to estimate PCs, then project whole cohort onto those PCs
pc_pruned <- pca_vcf(GOOD_VCF, K = 10,
					 skeleton = subset(pca_reps, keep)$iid,
					 samples = pca_reps$iid)

## add some manual annotations to plot
focus_pops <- merge_pca_result(pc_pruned, samples) %>%
	subset( (country %in% c("TW","IN","AF","KZ","CZ","IR","LB") & taxon %in% c("dom","mus","cas")) ) %>%
	group_by(country, taxon) %>%
	summarise(xcen = mean(PC1), ycen = mean(PC2),
			  pop = country_decoder[ country[1] ])

## taxon label coordinates in PC space
sslabs <- tibble( taxon = c("dom","mus","cas","gen","spretus"),
				  x = c(0, -150, -125, -200, -150),
				  y = c(60, 200, -130,   30,  -65) )

## draw PCA
p1 <- merge_pca_result(pc_pruned, samples) %>%
	subset(taxon != "hybrid") %>%
	mutate(taxon = factor_mus(taxon)) %>%
	ggplot() +
	geom_point(aes(x = PC1, y = PC2, colour = taxon)) +
	ggrepel::geom_text_repel(data = focus_pops,
							 aes(x = xcen, y = ycen, label = pop),
							 nudge_x = 50, nudge_y = -25, force = 3, size = 3,
							 segment.colour = "grey70", colour = "grey40") +
	geom_text(data = sslabs,
			  aes(x = x, y = y, label = taxon_labeller_full(taxon), colour = taxon),
			  fontface = "italic") +
	scale_colour_mus("", label = taxon_labeller_full, guide = FALSE) +
	pca_axes(pc_pruned) +
	theme_classic2() +
	theme(legend.text = element_text(face = "italic"))
ggsave("figures/pca_downsample_by_taxon_labelled.pdf", width = 5, height = 5)

## save this PCA
write_pca_result(pc_pruned, "pca/all_downsample_by_taxon.G90_MAF01")

## estimate Y chromosome tree
# read genotypes from *.bed file into matrix
X <- BEDMatrix::BEDMatrix("geno/isec/merged.renamed.Y", simple_names = TRUE)
X <- X[ rownames(X) %in% subset(samples, sex == "M")$iid, ]
unrel <- read_tsv("sample_lists/downsample.txt", col_names = "iid")$iid
unrel <- c(unrel, "Mus_spicilegus")

# keep only sites and samples without too much missing data
keep_sites <- colSums( !is.na(X) ) > 420
keep_inds <- rowSums( is.na(X) ) < 5

# calculate distance matrix
d <- dist(X[ rownames(X) %in% unrel & keep_inds, keep_sites ])
# infer NJ tree
tree <- ape::bionjs(d)
# root the tree with M. spicilegus as outgroup
tree <- ape::root(tree, "Mus_spicilegus")

## identify some key nodes in the tree, for plot annotations
# MRCA of 'Mediterranean' domesticus
x <- subset(samples, iid %in% unrel) %>%
	subset(taxon == "dom") %>%
	subset(country == "CY" | country == "IR") %>%
	.$iid
n1 <- ape::getMRCA(tree, intersect(x, tree$tip.label))

# MRCA of 'north/west European' domesticus
x <- subset(samples, iid %in% unrel) %>%
	subset(taxon == "dom") %>%
	subset(country == "AU" | country == "DE") %>%
	.$iid
n2 <- ape::getMRCA(tree, intersect(x, tree$tip.label))

# MRCA of Taiwan castaneus and all musculus
x <- subset(samples, iid %in% unrel) %>%
	subset(taxon == "mus") %>%
	#subset(country == "TW") %>%
	.$iid
n3 <- ape::getMRCA(tree, intersect(x, tree$tip.label))

## render tree using `ggtree`
tree_df <- fortify(tree, layout = "fan", branch.length = "none") %>%
	mutate(iid = label) %>%
	left_join(select(samples, iid, sex, taxon, country, dipnum, race))

p2 <- ggtree(tree_df, layout = "fan", branch.length = "none") +
	geom_segment(data = subset(tree_df, isTip),
				 aes(x = x, xend = x+8, y = y, yend = y),
				 lty = "solid", colour = "grey80", size = 0.5) +
	geom_cladelabel(n1, label = "", colour = "grey40",
					geom = "text", barsize = 1.5, offset = 2) +
	geom_cladelabel(n2, label = "", colour = "grey70",
					geom = "text", barsize = 1.5, offset = 2) +
	geom_tippoint(aes(colour = taxon)) +
	geom_tiplab2(aes(label = substr(iid, 1, 6)), size = 1.5,
				 offset = 8, hjust = 0) +
	scale_colour_mus()
print(p2)

## plot PCA and Y chromosome phylogeny combo
cowplot::plot_grid(p1, p2 + guides(colour = FALSE), nrow = 1,
				   labels = "AUTO", label_x = 0.9)
ggsave("figures/pca_downsample_by_taxon_labelled_with_Ytree.pdf", width = 10, height = 5)


## now do PCA in domesticus only for more detailed view of domesticus populations
pc_dom_pruned <- pca_vcf(GOOD_VCF, K = 10,
						 skeleton = subset(pca_reps, taxon == "dom" & keep)$iid,
						 samples = subset(pca_reps, taxon == "dom")$iid)

## assign population labels, color palettes and plotting symbols
pop_list <- unique(subset(pca_reps, taxon == "dom" & keep)$pop)
POP_COL <- setNames( repeat_colours(pop_list, "Paired"), pop_list )
POP_COL <- c(POP_COL, "(other)" = "grey92")
POP_SYM <- setNames( repeat_shapes(pop_list), pop_list )
POP_SYM <- c(POP_SYM, "(other)" = 19)

## assign samples to some continental-level subgroups in PCA space
## NB: this is just done by manual inspection
subgroups <- pc_dom_pruned %>%
	mutate(cluster = case_when(PC1 < -30 | PC2 > 60 ~ "1",
							   PC2 < -30 & PC1 > -30 & PC1 < 40 ~ "2",
							   TRUE ~ "3")) %>%
	select(iid, cluster) %>%
	as_tibble()

## write cluster result for later
write_tsv(subgroups, "sample_lists/dom_subgroups.txt", col_names = FALSE)

## draw PCA plot, colored by subgroup assignment
pp0 <- merge_pca_result(pc_dom_pruned, subgroups) %>%
	ggplot() +
	geom_point(aes(x = PC1, y = PC2, fill = cluster), shape = 21,
			   colour = "grey50") +
	scale_fill_manual(values = SUBGROUP_COLS, labels = subgroup_cluster_labeller) +
	pca_axes(pc_dom_pruned) +
	theme_classic2() +
	theme(legend.position = c(0.9,0.15)) +
	coord_equal()

## make seperate plots per subgroup
plot_pca_with_focus <- function(pcs, labs, ...) {
	
	keep_groups <- sort(unique(labs$focal_group))
	keep_pch <- c(POP_SYM[ keep_groups ])
	keep_col <- c(POP_COL[ keep_groups ])
	
	merge_pca_result(pcs, labs) %>%
		ggplot() +
		geom_point(aes(x = PC1, y = PC2, colour = focal_group, shape = focal_group)) +
		scale_shape_manual("population", values = keep_pch, na.translate=FALSE) +
		scale_colour_manual("population", values = keep_col, na.translate=FALSE) +
		pca_axes(pc_dom_pruned) +
		theme_classic2() +
		coord_equal()
	
}

pp1 <- left_join(subgroups, pca_pops) %>%
	mutate(focal_group = ifelse(cluster == "1", pop, "(other)")) %>%
	plot_pca_with_focus(pc_dom_pruned, .) +
	ggtitle(SUBGROUP_NAMES[1])

pp2 <- left_join(subgroups, pca_pops) %>%
	mutate(focal_group = ifelse(cluster == "2", pop, "(other)")) %>%
	plot_pca_with_focus(pc_dom_pruned, .) +
	ggtitle(SUBGROUP_NAMES[2])

pp3 <- left_join(subgroups, pca_pops) %>%
	mutate(focal_group = ifelse(cluster == "3", pop, "(other)")) %>%
	plot_pca_with_focus(pc_dom_pruned, .) +
	ggtitle(SUBGROUP_NAMES[3])

## show spatial gradient in PC1 by plotting samples on map of Eurasia
p4 <- merge_pca_result(pc_dom_pruned, samples) %>%
	ggplot() +
	geom_polygon(data = map_data("world", xlim = c(-18, 52), ylim = c(30,61)),
				 aes(x = long, y = lat, group = group),
				 fill = "grey92", colour = "white") +
	geom_point(aes(x = long, y = lat, fill = PC1),
			   position = position_jitter(width = 0.5, height = 0.5),
			   colour = "grey50", size = 1.5, shape = 21) +
	scale_fill_viridis_b(option = "viridis", direction = -1) +
	scale_x_continuous("longitude") +
	scale_y_continuous("latitude") +
	theme_bw() +
	coord_equal(xlim = c(-18, 52), ylim = c(30,61))

## multi-panel plot with call-outs for subgroups on PCA plot, plus map
cowplot::plot_grid(cowplot::plot_grid(pp0, p4, nrow = 1, labels = c("A","B")), 
				   cowplot::plot_grid(pp1, pp2, pp3, nrow = 1, labels = c("C","D","E")),
				   nrow = 2, rel_heights = c(1, 1.5))
ggsave("figures/pca_callouts_plus_map_eur.pdf", width = 14, height = 10)