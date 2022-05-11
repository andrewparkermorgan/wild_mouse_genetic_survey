#! /bin/bash

ROOT=$HOME/Dropbox/pmdvlab/wild_redux
WHERE=$ROOT/geno
EXCLUDE=$WHERE/exclusions.txt
RENAMER=$WHERE/id_matcher.txt

# 3-way intersection of Mega, Giga, WGS sets
bcftools isec \
	-Oz \
	-p $WHERE/isec \
	-n =3 \
	mm_gold.mm10.fixed.vcf.gz \
	nz_giga.gold.fixed.vcf.gz \
	raw.vcf.gz

# merge sites in the intersection, exclude dups, keep only biallelic sites, add AF/MAF tags
bcftools merge \
	-Ou \
	$WHERE/isec/0002.vcf.gz \
	$WHERE/isec/0001.vcf.gz \
	$WHERE/isec/0000.vcf.gz | \
bcftools view \
	-Ou \
	--samples-file ^$EXCLUDE | \
bcftools view \
	-Ou \
	-m2 -M2 \
	--type snps | \
bcftools +fill-tags \
	-Oz \
	>$WHERE/isec/merged.vcf.gzi && \
bcftools index --tbi $WHERE/isec/merged.vcf.gz

# rename WGS samples to match array conventions
bcftools reheader \
	-Oz \
	--samples $RENAMER \
	>$WHERE/isec/merged.renamed.vcf.gz && \
bcftools index --tbi $WHERE/isec/merged.renamed.vcf.gz
