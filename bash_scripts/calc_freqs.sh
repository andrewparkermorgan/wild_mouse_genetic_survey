#! /bin/bash

ROOT=$HOME/Dropbox/pmdvlab/wild_redux
VCF=$ROOT/geno/isec/merged.renamed.vcf.gz
OUT=$ROOT/geno/freq

for taxon in dom mus cas gen;
do
	echo "calculating allele freqs for taxon '$taxon' ..."
	bcftools view -S $ROOT/sample_lists/unrel_${taxon}.txt $VCF | \
	bcftools +fill-tags | bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' | \
	bgzip >$OUT/freq_${taxon}.txt.gz && tabix -s1 -b2 -e2 $OUT/freq_${taxon}.txt.gz
done
