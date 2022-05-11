#! /bin/bash

ROOT=$HOME/Dropbox/pmdvlab/wild_redux
VCF=$ROOT/geno/isec/merged.renamed.vcf.gz
FREQS=$ROOT/geno/freq

for taxon in dom mus cas gen spretus;
do
	bcftools roh \
		-S $ROOT/sample_lists/final_${taxon}.txt \
		-M 0.5e-6 -G30 \
		--AF-file $FREQS/freq_${taxon}.txt.gz \
		-R $ROOT/regions/chromosomes_AX.bed \
		-Or $VCF \
	>$ROOT/roh/${taxon}.roh
done
