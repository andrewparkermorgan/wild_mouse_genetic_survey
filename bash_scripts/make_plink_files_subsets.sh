#! /bin/bash

ROOT=$HOME/Dropbox/pmdvlab/wild_redux
WHERE=$ROOT/scratch/random_subsets/
VCF=$ROOT/geno/isec/merged.renamed.vcf.gz

for ii in 0 1 2 3 4 5 6 7 8 9;
do

	plink2 \
		--vcf $VCF \
		--keep $ROOT/kinship/random_subsets/dom_${ii}.txt \
		--autosome \
		--make-bed \
		--out $WHERE/dom_${ii}.A

	plink2 \
		--vcf $VCF \
		--keep $ROOT/kinship/random_subsets/dom_${ii}_downsample.txt \
		--autosome \
		--make-bed \
		--out $WHERE/dom_downsample_${ii}.A

done
