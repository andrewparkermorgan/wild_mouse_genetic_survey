# draw_admix_plots.R
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

## read list of samples kept in this analysis
fam <- read_fam("geno/isec/merged_dom.fam") %>%
	select(iid)

## read population labels
pops <- read_fam("geno/isec/merged.renamed.by_pop.fam") %>%
	select(fid, iid) %>%
	rename(pop = fid)

## function to shorten some population labels to make plot more legible
pop_relabeller <- function(f) {
	
	newpops <- c(
		"Georgia" = "US:GA",
		"Florida" = "US:FL",
		"Virginia" = "US:VA",
		"Maryland" = "US:MD",
		"Pennsylvania" = "US:PA",
		"New Hampshire" = "US:NH",
		"Orkneys" = "Orkney",
		"New Zealand" = "NZ"
	)
	newlabs <- newpops[f]
	nas <- is.na(newlabs)
	newlabs[ nas ] <- f[nas]
	return( unname(newlabs) )
	
}


## order populations by longitude
pop_order <- left_join(fam, pops) %>%
	left_join(samples) %>%
	mutate(pop = reorder(pop, long)) %>%
	with(., levels(pop))

## population 'centroids' to use when drawing map
pop_centroids <- left_join(fam, pops) %>%
	left_join(select(samples, iid, long, lat)) %>%
	group_by(pop) %>%
	summarize(long = median(long, na.rm = TRUE),
			  lat = median(lat, na.rm = TRUE))

## read admixture result for K = 2
Qa <- read_Q_matrix("admix/merged_dom.A.2.Q", iids = fam$iid)
so <- sort_by_cluster(Qa)
Qa <- tidy(Qa) %>% 
	left_join(pops)
Qa2 <- Qa
so2 <- so

## render barplot
p1 <- plot_admixture_labelled(Qa, sort_order = so, pop_order = pop_order, label_size = 2) +
	scale_fill_brewer(palette = "Paired", guide = FALSE) +
	ylab("ancestry proportion") +
	xlab("K = 2") +
	theme(axis.title.x = element_text())

## calculate coordinates for arrows from subpopulation labels to map locations
pop_labels <- attr(p1, "pop_labels")
pop_centroids <- left_join(pop_labels, pop_centroids) %>%
	mutate(relpos = pos/max(last)) %>%
	mutate(newpos = 365*relpos-180)

## load map data
# preloaded from `mapdata` package, good enough for low-res
map_df <- map_data("world")

toplot <- Qa2 %>%
	left_join(samples) %>%
	subset(variable == "pop2")

## plot N Eur-like ancestry component on world map
point_cols <- RColorBrewer::brewer.pal(2,"Paired")[1:2]
p5 <- subset(map_df, region != "Antarctica") %>%
	ggplot() +
	geom_polygon(aes(x = long, y = lat, group = group), fill = "grey92", colour = "white") +
	geom_point(data = toplot,
			   aes(x = long, y = lat, colour = value)) +
	geom_segment(data = pop_centroids,
				 aes(x = newpos, xend = long, y = -Inf, yend = lat),
				 colour = "grey70") +
	scale_colour_gradientn("ancestry\nproportion", colours = c(point_cols[1], point_cols[1], point_cols[2], point_cols[2]),
						   breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
	coord_cartesian(xlim = c(-180,185), expand = FALSE) +
	theme_void() +
	theme(legend.position = c(0.9,0.9),
		  legend.justification = c(1,1))

## re-draw barplot, this time for real
pops2 <- mutate(pops, pop = pop_relabeller(pop))
pop_order2 <- pop_relabeller(pop_order)
Qa <- read_Q_matrix("admix/merged_dom.A.2.Q", iids = fam$iid)
so <- sort_by_cluster(Qa)
Qa <- tidy(Qa) %>% 
	left_join(pops2)

p1 <- plot_admixture_labelled(Qa, sort_order = so, pop_order = pop_order2, ymax = 1.37) +
	scale_fill_brewer(palette = "Paired", guide = FALSE) +
	ylab("ancestry proportion") +
	xlab("K = 2") +
	theme(axis.title.x = element_text())

cowplot::plot_grid(p5, p1, nrow = 2, align = "v", labels = "AUTO")
ggsave("figures/admix_dom_projected_bypop_K2_plus_map.pdf", width = 10, height = 10)
