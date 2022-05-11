# helpers.R
# constants, helper functions etc for 'wild mice redux' manuscript

# this one needs work
inverse_rank_transform <- function(x) {
	
	r <- rank(x)/length(x)
	print(r)
	z <- qnorm(r)
	return(z)
	
}

logit <- function(p, eps = 1e-4) {
	p[ p < eps ] <- eps
	log(p/(1-p))
}

inverse_logit <- function(p, eps = 1e-4) {
	exp(p)/(1+exp(p))
}

pca_axis_labs <- function(pc, K = 1:2, ...) {
	varprop <- round(attr(pc, "explained")*100, 1)
	paste0("PC", K, " (", varprop[K], "%)")
}

pca_axes <- function(pc, K = 1:2, ...) {
	labs <- pca_axis_labs(pc, K)
	list( scale_x_continuous(labs[1]),
		  scale_y_continuous(labs[2]) )
}

read_plink2_pca <- function(ff, K = 10, ...) {

	pcs <- read_tsv(paste0(ff, ".eigenvec"))
	colnames(pcs) <- c("fid","iid", paste0("PC", 1:K))
	evs <- readLines(paste0(ff, ".eigenval")) %>% as.double()

	class(pcs) <- c("pca_result", class(pcs))
	attr(pcs, "eigvals") <- evs
	attr(pcs, "explained") <- evs/sum(evs)

	return(pcs)

}

country_decoder <-
	c("AU" = "Australia",
	  "BE" = "Belgium",
	  "CH" = "Switzerland",
	  "CY" = "Cyprus",
	  "DE" = "Germany",
	  "DK" = "Denmark",
	  "ES" = "Spain",
	  "FR" = "France",
	  "GR" = "Greece",
	  "IR" = "Iran",
	  "IT" = "Italy",
	  "LB" = "Lebanon",
	  "PT" = "Portugal",
	  "TN" = "Tunisia",
	  "UK" = "United Kingdom",
	  "AF" = "Afghanistan",
	  "KZ" = "Kazakhstan",
	  "TW" = "Taiwan",
	  "PL" = "Poland",
	  "RU" = "Russia",
	  "CZ" = "Czech Rep",
	  "IN" = "India")

SUBGROUP_COLS <- setNames(rev(viridis::viridis(3)), c("1","2","3"))
SUBGROUP_NAMES <- c("1" = "Mediterranean", "2" = "Central Europe", "3" = "Northern Europe")

subgroup_cluster_labeller <- function(f) {

	return( SUBGROUP_NAMES[ as.character(f) ] )

}

recenter <- function(x) {
	ifelse(x < 0, x+360, x)
}

unrecenter <- function(x) {
	ifelse(x >= 180, x-360, x)
}

clip_polygon <- function(shp, bb, ...) {

	if (class(bb) == "matrix")
		b_poly <- as(raster::extent(as.vector(t(bb))), "SpatialPolygons")
	else
		b_poly <- as(raster::extent(bb), "SpatialPolygons")

	rgeos::gIntersection(shp, b_poly, byid = TRUE)

}

find_bounds <- function(pp, buffer = 0.05, ...) {
	
	lat_bounds <- range(pp$lat)
	long_bounds <- range(pp$long)
	
	yr <- diff(lat_bounds)
	xr <- diff(long_bounds)
	
	ya <- buffer*yr/2
	xa <- buffer*xr/2
	
	lat_bounds <- lat_bounds + c(-ya, ya)
	long_bounds <- long_bounds + c(-xa, xa)
	
	rez <- matrix(c(long_bounds, lat_bounds), byrow = TRUE, nrow = 2)
	return(rez)
	
}

theme_legend_lowerright <- function(...) {
	theme(legend.position = c(1,0),
		  legend.justification = c(1,0))
}

# for some reason ggh4x package needs but cannot find `isFALSE()`
# not sure what it is supposed to do, but this is my guess
isFALSE <- function(x) { 
	identical(x,FALSE)
	#if(length(x) > 0)
	#	!as.logical(x)
	#else
	#	TRUE
}