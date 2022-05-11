# draw_maps.R
# Draw pretty maps for main figures
setwd("~/Dropbox/pmdvlab/wild_redux/")

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

devtools::load_all("~/Dropbox/pmdvlab/mouser")
devtools::load_all("~/Dropbox/pmdvlab/popcorn")
source("helpers.R")
rb_table <- readRDS("rb_table.rds")

## read and prune sample metadata
samples <- readxl::read_xlsx("~/Dropbox/pmdvlab/farallon/samples.xlsx")
keeps <- read_tsv("geno/good_samples.txt", col_names = "iid")
samples <- left_join(keeps, samples)
samples$karyo <- with(samples, ifelse(dipnum < 40, "Rb", "standard"))
dim(samples) # should be 814

## load map data
# preloaded from `mapdata` package, good enough for low-res
map_df <- map_data("world")

# load Portugal shapefiles
shp <- rgdal::readOGR("geo/PRT_adm1.shp")
lat_bounds <- c(32.6, 33.5)
long_bounds <- c(-17.5, -16)
newbb <- matrix(c(long_bounds, lat_bounds), byrow = TRUE, nrow = 2)
shp_df <- clip_polygon(shp, newbb) %>%
	broom::tidy(region = "NAME_1")

p5 <- subset(map_df, region != "Antarctica") %>%
	ggplot() +
	geom_polygon(aes(x = long, y = lat, group = group), fill = "grey92", colour = "white") +
	geom_point(data = subset(samples, taxon != "hybrid") %>% mutate(taxon = factor_mus(taxon)),
			   aes(x = long, y = lat, colour = taxon)) +
	scale_x_continuous(position = "top") +
	#scale_x_continuous(limits = c(-180,180)) +
	scale_colour_manual("", values = mus_colors()[ c("dom",'mus',"cas","gen","spretus","spicilegus") ],
						labels = taxon_labeller_full) +
	#scale_colour_mus("", label = taxon_labeller_full, drop = TRUE) +
	coord_cartesian(xlim = c(-180,180), expand = TRUE) +
	theme_void() +
	theme(legend.title = element_text(face = "bold"),
		  legend.text = element_text(face = "italic"),
		  #legend.background = element_rect(colour = "grey92"),
		  legend.position = c(0.05,0.05),
		  legend.justification = c(0,0))

p6 <- subset(samples, taxon == "dom") %>%
	mutate(taxon = factor_mus(taxon), is_rb = (dipnum != 40)) %>%
	ggplot() +
	geom_polygon(data = map_data("world", xlim = c(-18, 52), ylim = c(30,61)),
				 aes(x = long, y = lat, group = group),
				 colour = "white", fill = "grey92") +
	geom_point(aes(x = long, y = lat, colour = taxon, shape = is_rb),
			   position = position_jitter(width = 0.5, height = 0.5),
			   fill = "white") +
	scale_x_continuous("longitude") +
	scale_y_continuous("latitude") +
	scale_shape_manual("karyotype", values = c(19,21), labels = c("standard","Rb")) +
	scale_colour_mus(guide = FALSE) +
	theme_bw() +
	theme(legend.position = "top") +
	coord_equal(xlim = c(-18, 52), ylim = c(30,61))

recenter <- function(x) {
	ifelse(x < 0, x+360, x)
}
unrecenter <- function(x) {
	ifelse(x >= 180, x-360, x)
}

p7 <- subset(samples, taxon == "dom" & country == "NZ") %>%
	mutate(taxon = factor_mus(taxon), is_rb = (dipnum != 40)) %>%
	ggplot() +
	geom_polygon(data = map_data("world", region = "New Zealand"),
				 aes(x = recenter(long), y = lat, group = group),
				 fill = "grey92") +
	geom_point(aes(x = recenter(long), y = lat, colour = taxon),
			   position = position_jitter(width = 0.5, height = 0.5),
			   fill = "white") +
	scale_x_continuous("longitude", labels = unrecenter) +
	scale_y_continuous("latitude") +
	scale_colour_mus(guide = FALSE) +
	annotate("text", label = "New Zealand", x = 0.1, y = 0.9) +
	theme_bw() +
	coord_equal(ylim = c(-52,-33), xlim = c(165,185))

p8 <- subset(samples, taxon == "dom" & country == "PT") %>%
	mutate(taxon = factor_mus(taxon), is_rb = (dipnum != 40)) %>%
	subset(long < -12.5) %>%
	ggplot() +
	geom_polygon(data = shp_df,
				 aes(x = long, y = lat, group = group),
				 fill = "grey92") +
	geom_point(aes(x = long, y = lat, colour = taxon, shape = is_rb),
			   position = position_jitter(width = 0.005, height = 0.005),
			   fill = "white") +
	scale_x_continuous("longitude", labels = unrecenter) +
	scale_y_continuous("latitude") +
	scale_shape_manual("karyotype", values = c(19,21), labels = c("standard","Rb"),
					   guide = FALSE) +
	scale_colour_mus(guide = FALSE) +
	theme_bw() +
	coord_equal()

print(p5 + coord_equal())
ggsave("figures/world_map_by_taxon.pdf", width = 8, height = 4)

cowplot::plot_grid(p6,
				   cowplot::plot_grid(p8, p7, nrow = 1,
				   				   labels = c("B","C")),
				   nrow = 2, rel_heights = c(2,1),
				   labels = c("A",""))
ggsave("figures/domesticus_eur_with_insets.pdf", width = 8, height = 8)
