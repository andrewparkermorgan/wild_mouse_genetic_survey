#! /bin/bash

ROOT=$HOME/Dropbox/pmdvlab/wild_redux
VCF=$ROOT/geno/isec/merged.renamed.vcf.gz
OUT=$ROOT/freq/random_subsets

mkdir -p $OUT

for taxon in dom mus cas gen;
do
	echo "calculating allele freqs for taxon '$taxon' ..."
	for ii in 0 1 2 3 4 5 6 7 8 9;
	do
		echo " ... $ii"
		outfile=$OUT/${taxon}_${ii}.txt.gz
		bcftools view -S $ROOT/kinship/random_subsets/${taxon}_${ii}.txt $VCF | \
		bcftools +fill-tags | bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' | \
		bgzip >$outfile && tabix -s1 -b2 -e2 $outfile
	done
done
