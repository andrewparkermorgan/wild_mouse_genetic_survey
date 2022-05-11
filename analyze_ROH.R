# analyze_ROH.R
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
sampels$taxon_old <- samples$taxon
samples$taxon <- factor_mus(samples$taxon)
dim(samples) # should be 814

add_roh_cM_pos <- function(df, themap, ...) {
	
	df$pos <- df$start
	dfs <- left_join(df, select(themap, chr, cM, pos))
	df$pos <- df$end
	dfe <- left_join(df, select(themap, chr, cM, pos))
	df$start_cM <- dfs$cM
	df$end_cM <- dfe$cM
	df$width_cM <- with(df, end_cM - start_cM)
	return(df)
	
}

rb_race_labeller <- function(x, ...) {
	
	df <- select(samples, race, dipnum) %>% unique()
	rr <- with(df, paste0(race, " (", dipnum, ")"))
	races <- setNames(rr, df$race)
	rez <- races[ x ]
	nas <- is.na(rez)
	rez[ nas ] <- paste0(x[ nas ], " (0)")
	return(rez)
	
}

## read genetic map
chrlen <- chromsizes_cM()[1:20] %>% tibble::enframe("chr","len")
map <- read_map("geno/mm_gold.mm10.bim")
thebim <- read_tsv("geno/isec/merged.renamed.bim", col_names = c("chr","marker","cM","pos","A1","A2"), col_types = c("ccdicc"))
thebim$cM <- map$cM[ match(thebim$marker, map$marker) ]
thebim$chr <- gsub("^MT","M", thebim$chr)
thebim$chr <- factor_chrom(thebim$chr)

## read inferred ROH segments
roh <- bind_rows( read_roh("roh/dom.roh"),
				  read_roh("roh/mus.roh"),
				  read_roh("roh/cas.roh"),
				  read_roh("roh/gen.roh") )
nrow(roh)

## add genetic positions
roh <- add_roh_cM_pos(roh, map)
nrow(roh)

## summarize length of ROH segments
left_join(roh, select(samples, iid, taxon, sex)) %>%
	subset(is_autosome(chr) & width_cM > 2) %>%
	summarise(mu = mean(width/1e6), mu_sites = mean(nsites),
			  med = median(width/1e6), med_sites = median(nsites))

## basic diagnostics: count, total lenght (Mb, cM) of ROH segs per individual
roh2 <- roh %>%
	mutate(chrtype = factor_chromtype(chr)) %>%
	subset(width_cM > 2) %>%
	group_by(iid, chrtype) %>%
	summarize(roh_tot_cM = sum(width_cM),
			  roh_tot_Mb = sum(width/1e6),
			  roh_count = sum(width > 0))
roh2 <- select(samples, iid) %>%
	left_join(roh2) %>%
	replace_na(list(chrtype = "A", roh_tot_cM = 0, roh_tot_Mb = 0, roh_count = 0L))

## get rid of ROH segs on X in males, which are fake
roh2 <- left_join(roh2, samples) %>%
	subset(!(chrtype == "X" & sex == "M")) %>%
	mutate(taxon = factor_mus(taxon))

## estimate F_hat per individual
fhat <- subset(roh2, chrtype == "A") %>%
	group_by(iid, chrtype) %>%
	summarise(fhat = sum(roh_tot_Mb)/sum((chromsizes_mm10()[1:19]/1e6-3))) %>%
	inner_join(samples) %>%
	mutate(fhat_logit = logit(fhat))
select(fhat, iid, fhat) %>%
	write_csv("roh/inbreeding_coefs_A.csv", col_names = FALSE)

fhat_x <- subset(roh2, chrtype == "X") %>%
	group_by(iid, chrtype) %>%
	summarise(fhat = sum(roh_tot_Mb)/sum((chromsizes_mm10()[20]/1e6-3))) %>%
	inner_join(samples) %>%
	mutate(fhat_logit = logit(fhat))
select(fhat_x, iid, fhat) %>%
	write_csv("roh/inbreeding_coefs_X.csv", col_names = FALSE)

## read population labels
pca_pops <- read_tsv("geno/isec/merged.renamed.by_pop.fam", col_names = c("fid","iid","mom","dad","sex","pheno")) %>%
	rename(pop = fid) %>%
	select(iid, pop)

pop_regrouper <- list(
	"South Pacific" = c("Ruapuke","Antipodes","Australia","Pitt","Chatham","Auckland"),
	"North America" = c("Florida","Georgia","Virginia","Maryland","Pennsylvania","New Hampshire")
) %>%
	lapply(., tibble::enframe, value = "pop") %>%
	bind_rows(.id = "newpop")

median_hilo <- function(x, alpha = 0.05) {
	med <- stats::median(x, na.rm = TRUE)
	lo <- stats::quantile(x, alpha)
	hi <- stats::quantile(x, 1-alpha)
	data.frame(y = med, ymin = lo, ymax = hi)
}

fhat_mus <- fhat %>%
	mutate(taxon = factor_mus(taxon)) %>%
	subset(taxon %in% c("dom","mus","cas"))
pwtest <- ggpubr::compare_means(fhat ~ taxon, fhat_mus, method = "wilcox.test")

lm(fhat_logit ~ taxon, fhat_mus) %>% anova()
# Analysis of Variance Table
# 
# Response: fhat_logit
# Df  Sum Sq Mean Sq F value   Pr(>F)
# taxon       2   27.63 13.8163  4.8826 0.007816 **
# Residuals 760 2150.55  2.8297

lm(fhat_logit ~ taxon, fhat_mus) %>% lsmeans::lsmeans(pairwise~taxon)
# $lsmeans
# taxon    lsmean         SE  df  lower.CL  upper.CL
# dom   -2.253560 0.06353446 760 -2.378284 -2.128836
# mus   -2.683042 0.35075544 760 -3.371606 -1.994477
# cas   -3.064510 0.26936182 760 -3.593291 -2.535728
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast   estimate        SE  df t.ratio p.value
# dom - mus 0.4294815 0.3564632 760   1.205  0.4506
# dom - cas 0.8109496 0.2767533 760   2.930  0.0098
# mus - cas 0.3814681 0.4422501 760   0.863  0.6641
# 
# P value adjustment: tukey method for comparing a family of 3 estimates 

## inbreeding by taxon
p0 <- fhat_mus %>%
	ggplot(., aes(x = taxon, y = fhat, fill = taxon, colour = taxon)) +
	ggbeeswarm::geom_quasirandom(shape = 21, colour = "white", varwidth = TRUE) +
	stat_summary(fun.data = median_hilo, position = position_nudge(x = 0.3),
				 size = 0.5, geom = "linerange") +
	stat_summary(fun.data = function(f) median_hilo(f, 0.25), position = position_nudge(x = 0.3),
				 size = 1.25, geom = "linerange") +
	stat_summary(fun.data = median_hilo, position = position_nudge(x = 0.3),
				 size = 3, geom = "point", colour = "white", shape = 21) +
	#stat_summary(fun.data = function(f) ggpubr::mean_sd(f), position = position_nudge(x = 0.3)) +
	scale_y_continuous(expression(hat(italic(F))["ROH"]), limits = c(-0.05,1)) +
	scale_x_discrete("taxon", labels = taxon_labeller) +
	scale_fill_mus(guide = FALSE) + 
	scale_colour_mus(guide = FALSE) +
	theme_classic2() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1),
		  axis.title.x = element_blank())
print(p0)

# compare autos to X
avsx <- bind_rows(fhat, fhat_x) %>%
	mutate(taxon = factor_mus(taxon)) %>%
	select(iid, sex, taxon, fhat, chrtype) %>%
	spread(chrtype, fhat, fill = 0) %>%
	subset(sex == "F")
x <- cor.test(~ A + X, data = avsx, method = "spearman", use = "pairwise.complete.obs")
x$p.value

# g+ xlab( expression(paste("Value is ", sigma,",", R^{2},'=',r2.value)))

p00 <- avsx %>%
	subset(taxon %in% c("dom","mus","cas")) %>%
	ggplot(., aes(x = A, y = X)) +
	geom_smooth(colour = "black", fill = "grey90", size = 0.5, method = "lm") +
	geom_point(aes(fill = taxon), pch = 21, colour = "white", size = 2) +
	annotate("text", x = 0.9, y = 0.9, parse = TRUE,
			 label = as.expression( bquote(atop(rho==.(round(unname(x$estimate),3)),
			 								   italic(p)==.(signif(x$p.value,3)))) )) +
	scale_x_continuous(expression(hat(italic(F))["ROH"]~" autosomes"), limits = c(-0.01,1)) +
	scale_y_continuous(expression(hat(italic(F))["ROH"]~" X chromosome"), limits = c(-0.01,1)) +
	scale_fill_mus(label = taxon_labeller) +
	theme_classic2() +
	theme(legend.position = c(1,0),
		  legend.justification = c(1,0))

toplot <- left_join(fhat, pca_pops) %>%
	subset(taxon == "dom") %>%
	left_join(pop_regrouper) %>%
	mutate(newpop = case_when(is.na(newpop) ~ pop,
							  TRUE ~ newpop))

## inbreeding by population
p1 <- toplot %>%
	ggplot(., aes(x = reorder(newpop, fhat), y = fhat)) +
	ggbeeswarm::geom_quasirandom(shape = 21, fill = "black", colour = "white", groupOnX = TRUE) +
	scale_y_continuous(expression(hat(italic(F))["ROH"]), limits = c(0,1)) +
	scale_x_discrete("location") +
	scale_shape_manual(values = c(21,19), guide = FALSE) +
	theme_classic2() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1),
		  axis.title.x = element_blank(),
		  panel.grid.major.x = element_line(colour = "grey92", size = 0.5))
print(p1)

lm(fhat_logit ~ newpop, toplot) %>% anova()
# Analysis of Variance Table
# 
# Response: fhat_logit
# Df Sum Sq Mean Sq F value    Pr(>F)    
# newpop     24  649.6 27.0668  14.362 < 2.2e-16 ***
# Residuals 676 1274.0  1.8846 

p2 <- toplot %>%
	mutate(is_rb = dipnum < 40) %>%
	ggplot(., aes(x = is_rb, y = fhat)) +
	#ggbeeswarm::geom_quasirandom(shape = 21, fill = "black", colour = "white", groupOnX = TRUE) +
	ggbeeswarm::geom_quasirandom(colour = "grey70", groupOnX = TRUE) +
	stat_summary(fun.data = median_hilo, position = position_nudge(x = 0.3),
				 size = 0.5, geom = "linerange") +
	stat_summary(fun.data = function(f) median_hilo(f, 0.25), position = position_nudge(x = 0.3),
				 size = 1.25, geom = "linerange") +
	stat_summary(fun.data = median_hilo, position = position_nudge(x = 0.3),
				 size = 3, geom = "point", fill = "black", colour = "white", shape = 21) +
	scale_y_continuous(expression(hat(italic(F))["ROH"]), limits = c(0,1)) +
	scale_x_discrete(labels = c("standard","Rb")) +
	theme_classic2() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1),
		  axis.title.x = element_blank(),
		  panel.grid.major.x = element_line(colour = "grey92", size = 0.5))
print(p2)

## combined plot, beeswarm by taxon and A vs X
pp1 <- cowplot::plot_grid(p0, p00, nrow = 1, align = "v", labels = c("A","B"))
## combined plot, beeswarm by pop plus stnd vs Rb
pp2 <- cowplot::plot_grid(p1, p2, rel_widths = c(2,1), labels = c("C","D"),
						  align = "v")

## 4-up composite
cowplot::plot_grid(pp1, pp2, nrow = 2)
ggsave("figures/inbreeding_by_subspecies_4up.pdf", width = 8, height = 6)


## now investigate inbreeding at finer scale: Farallon, Floreana, Gough ,Maryland
keep_sites <- c("HBA" = "Chevy Chase",
				"HBB" = "North Potomac",
				"HBC" = "Laurel",
				"CVL" = "Ctrvle",
				"SEF" = "Farallon",
				"FLO" = "Floreana",
				"GOU" = "Gough")
site_grouper <- c("HBA" = "Maryland",
				  "HBB" = "Maryland",
				  "HBC" = "Maryland",
				  "CVL" = "Maryland",
				  "SEF" = "Islands",
				  "FLO" = "Islands",
				  "GOU" = "Islands")
toplot <- subset(fhat, locale %in% names(keep_sites)) %>%
	ungroup() %>%
	mutate(site2 = keep_sites[locale],
		   group = factor(site_grouper[locale])) %>%
	group_by(country, locale, site2, group)
site_order <- with(toplot, reorder(reorder(site2, fhat), as.numeric(group))) %>% levels()

f_by_site <- toplot %>%
	summarise(fhat_mean = mean(fhat, na.rm = TRUE),
			  fhat_med = median(fhat, na.rm = TRUE),
			  sigma = sd(fhat, na.rm = TRUE),
			  mad = mad(fhat, na.rm = TRUE),
			  q25 = quantile(fhat, 0.25, na.rm = TRUE),
			  q75 = quantile(fhat, 0.75, na.rm = TRUE))

## test for variation across sites
lm(fhat_logit ~ site2, toplot) %>% anova()
# Analysis of Variance Table
# 
# Response: fhat_logit
# Df Sum Sq Mean Sq F value    Pr(>F)    
# site2       6 170.44 28.4072  40.185 < 2.2e-16 ***
# Residuals 154 108.86  0.7069  

## inbreeding by locale, one point per mouse
p1 <- toplot %>%
	ggplot(., aes(x = reorder(site2, fhat), y = fhat)) +
	ggbeeswarm::geom_quasirandom(shape = 21, fill = "black", colour = "white", groupOnX = TRUE) +
	scale_y_continuous(expression(hat(italic(F))["ROH"]), limits = c(0,1)) +
	scale_x_discrete("location") +
	scale_shape_manual(values = c(21,19), guide = FALSE) +
	facet_grid(. ~ group, scale = "free_x", space = "free_x") +
	theme_classic2() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1),
		  axis.title.x = element_blank(),
		  panel.grid.major.x = element_line(colour = "grey92", size = 0.5))

## examine distribution of ROH segment lenghts
cM_bin_labels <- c("2-4 cM","4-8 cM","8-12 cM","12-20 cM",">20 cM")
# sum ROH by length bin, per mouse
roh_binned <- roh %>%
	mutate(chrtype = factor_chromtype(chr)) %>%
	subset(chrtype == "A") %>%
	subset(width_cM > 2) %>%
	mutate(roh_bin = cut(width_cM, c(2,4,8,12,20,Inf), labels = cM_bin_labels)) %>%
	group_by(iid, roh_bin) %>%
	summarise(bin_tot_cM = sum(width_cM)) %>% 
	complete(roh_bin, fill = list(bin_tot_cM = 0)) %>% # keep bins which are empty, make them 0
	mutate(genome_tot_cM = sum(bin_tot_cM),
		   prop_genome = bin_tot_cM/genome_tot_cM)

# summarize binned ROH by location
roh_bin_summ <- roh_binned %>%
	left_join(select(samples, iid, country, locale)) %>%
	subset(locale %in% names(keep_sites)) %>%
	mutate(site2 = keep_sites[locale], group = site_grouper[locale]) %>%
	group_by(country, site2, roh_bin) %>%
	summarise(mean_cM = mean(bin_tot_cM),
			  sd_cM = sd(bin_tot_cM),
			  med_cM = median(bin_tot_cM),
			  mad_cM = mad(bin_tot_cM),
			  llq = quantile(bin_tot_cM, 0.05),
			  lq = quantile(bin_tot_cM, 0.25),
			  uq = quantile(bin_tot_cM, 0.75),
			  uuq = quantile(bin_tot_cM, 0.95))
print(roh_bin_summ)

roh_binned %>%
	left_join(select(samples, iid, country, locale)) %>%
	subset(locale %in% names(keep_sites)) %>%
	mutate(site2 = keep_sites[locale], group = site_grouper[locale]) %>%
	ggplot(., aes(x = site2, y = bin_tot_cM, colour = roh_bin, fill = roh_bin)) +
	stat_summary(fun.data = median_hilo, position = position_dodge(width = 0.5),
				 size = 0.5, geom = "linerange") +
	stat_summary(fun.data = function(f) median_hilo(f, 0.25), position = position_dodge(width = 0.5),
				 size = 1.25, geom = "linerange") +
	stat_summary(fun.data = median_hilo, position = position_dodge(width = 0.5),
				 size = 3, geom = "point", colour = "white", shape = 21) +
	scale_y_continuous("sum ROH length (cM)") +
	scale_colour_brewer("length bin", palette = "Spectral", direction = -1) +
	scale_fill_brewer("length bin", palette = "Spectral", direction = -1) +
	theme_classic2() +
	theme_slanty_x()

# formal test of difference in distribution?
tmp <- roh %>%
	mutate(chrtype = factor_chromtype(chr)) %>%
	subset(chrtype == "A") %>%
	subset(width_cM > 2) %>%
	left_join(samples) %>%
	subset(locale %in% names(keep_sites)) %>%
	mutate(site2 = factor(keep_sites[locale], site_order))

kSamples::ad.test(width_cM ~ site2, data = tmp)
kSamples::ad.test(width_cM ~ site2, data = subset(tmp, country == "US") %>% droplevels())

# obtain matrix of (iid x length bin)
tmp <- left_join(roh_binned, samples) %>% 
	subset(locale %in% names(keep_sites)) %>%
	mutate(site2 = keep_sites[locale], group = site_grouper[locale])
tmp_wide <- select(tmp, iid, roh_bin, bin_tot_cM) %>%
	spread(roh_bin, bin_tot_cM, fill = 0) %>%
	ungroup()
X <- select(tmp_wide, -iid) %>% as.matrix()
rownames(X) <- tmp_wide$iid
D <- vegan::vegdist(X)
meta_df <- tibble(iid = rownames(X)) %>%
	left_join(select(samples, iid, locale)) %>%
	mutate(site2 = keep_sites[locale], group = site_grouper[locale])

# how many individuals in each pop with at least one ROH >20 cM?
subset(tmp) %>% 
	group_by(site2) %>%
	subset(roh_bin == ">20 cM") %>%
	summarise(p_no_long_ROH = mean(bin_tot_cM == 0),
			  n_no_long_ROH = sum(bin_tot_cM == 0),
			  med_cM_long_ROH = median(bin_tot_cM),
			  mean_cM_long_ROH = mean(bin_tot_cM))
subset(tmp) %>%
	group_by(site2)

# test for difference in binned distribution of ROH by location
vegan::adonis2(D ~ site2, data = meta_df)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# vegan::adonis2(formula = D ~ site2, data = meta_df)
# Df SumOfSqs      F Pr(>F)    
# site2      6   16.046 32.072  0.001 ***
# Residual 153   12.758    

## plot binned distribution of ROH seg lengths, one column per mouse
p2 <- left_join(roh_binned, samples) %>% 
	subset(locale %in% names(keep_sites)) %>%
	mutate(site2 = keep_sites[locale], group = site_grouper[locale]) %>%
	mutate(site2 = factor(site2, site_order)) %>%
	ggplot() + 
	geom_bar(aes(x = reorder(iid, genome_tot_cM), y = bin_tot_cM, fill = roh_bin), stat = "identity") + 
	ggh4x::facet_nested(~ group + site2, space = "free_x", scale = "free_x",
						nest_line = element_line(linetype = "solid", color = "grey70")) + 
	scale_y_continuous("sum ROH length (cM)") +
	scale_x_discrete("individuals") +
	scale_fill_brewer("length bin", palette = "Spectral", direction = -1) +
	theme_classic2() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))
print(p2)

## CDF for each mouse with group summary, for supplement
roh %>%
	mutate(chrtype = factor_chromtype(chr)) %>%
	subset(chrtype == "A") %>%
	subset(width_cM > 2) %>%
	left_join(samples) %>%
	subset(locale %in% names(keep_sites)) %>%
	mutate(site2 = factor(keep_sites[locale], site_order)) %>% 
	ggplot() +
	stat_ecdf(aes(x = width_cM, group = iid), alpha = 0.2) +
	stat_ecdf(aes(x = width_cM), colour = "#2B83BA", size = 1) +
	scale_y_continuous("proportion of segments") +
	scale_x_log10("segment length (cM)") +
	facet_grid(site2 ~ .) +
	theme_classic2()
ggsave("figures/roh_cdf_by_indiv_x_locale.pdf", width = 4, height = 8)

p3 <- roh %>%
	mutate(chrtype = factor_chromtype(chr)) %>%
	subset(chrtype == "A") %>%
	subset(width_cM > 2) %>%
	left_join(samples) %>%
	subset(locale %in% names(keep_sites)) %>%
	mutate(site2 = factor(keep_sites[locale], site_order)) %>% 
	ggplot() +
	stat_ecdf(aes(x = width_cM, group = iid), alpha = 0.2) +
	stat_ecdf(aes(x = width_cM), colour = "#2B83BA", size = 1) +
	scale_y_continuous("proportion of segments", breaks = c(0, 0.5, 1.0)) +
	scale_x_log10("segment length (cM)") +
	facet_wrap(~ site2) +
	theme_classic2()

cowplot::plot_grid( cowplot::plot_grid(p1, NULL, p3, ncol = 3, rel_widths = c(1, 0.1, 1), labels = c("A","","B")),
					p2 + theme(axis.text.x = element_blank()), nrow = 2, rel_heights = c(1,1), labels = c("","C") )
ggsave("figures/roh_segment_breakdown_focus_pops.pdf", width = 10, height = 6)