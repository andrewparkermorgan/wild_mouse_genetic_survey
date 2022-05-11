#! /bin/bash

ROOT=$HOME/Dropbox/pmdvlab/wild_redux
VCF=$ROOT/geno/isec/merged.renamed.vcf.gz

rm -rf $ROOT/sample_lists/unrel.txt

for taxon in dom mus cas gen spretus;
do
	echo "kinship for $taxon ..."
	akt kin \
		-M 0 \
		-R $ROOT/regions/chromosomes_A.bed \
		--samples-file $ROOT/sample_lists/final_${taxon}.txt $VCF \
		>$ROOT/kinship/${taxon}.kin
	akt relatives -k 0.1 -i 10 $ROOT/kinship/${taxon}.kin >$ROOT/kinship/${taxon}.rel
	grep "^Unrel" $ROOT/kinship/${taxon}.rel | cut -f2 >$ROOT/kinship/unrel_${taxon}.txt
	#cp $ROOT/kinship/unrel_${taxon}.txt $ROOT/sample_lists/
	#cat $ROOT/kinship/unrel_${taxon}.txt >>$ROOT/sample_lists/unrel.txt
	mkdir -p $ROOT/kinship/random_subsets
	for ii in 0 1 2 3 4 5 6 7 8 9;
	do
		akt relatives -k 0.1 -i 10 <(sort -R $ROOT/kinship/${taxon}.kin)| grep "^Unrel" | cut -f2 >$ROOT/kinship/random_subsets/${taxon}_${ii}.txt
	done
done

cat $ROOT/kinship/unrel_*.txt >$ROOT/kinship/unrel.txt
cp $ROOT/kinship/unrel_*.txt $ROOT/sample_lists/
cp $ROOT/kinship/unrel.txt $ROOT/sample_lists/
cp $ROOT/kinship/unrel.txt $ROOT/kinship/all_unrel.txt
