#! /bin/bash

ROOT=$HOME/Dropbox/pmdvlab/wild_redux
VCF=$ROOT/geno/isec/merged.renamed.vcf.gz
OUT=$ROOT/scratch

plink2 \
	--vcf $VCF \
	--autosome \
	--make-bed \
	--out $OUT/merged.renamed.A

plink2 \
	--vcf $VCF \
	--keep sample_lists/final_dom.txt \
	--autosome \
	--make-bed \
	--out $OUT/merged_dom.A

plink2 \
	--vcf $VCF \
	--keep sample_lists/unrel_dom.txt \
	--autosome \
	--make-bed \
	--out $OUT/merged_unrel_dom.A

plink2 \
	--vcf $VCF \
	--keep sample_lists/downsample_dom.txt \
	--autosome \
	--make-bed \
	--out $OUT/merged_unrel_dom_downsample.A


