Genotypes for 814 wild house mice
===

This repository contains genotypes in VCF format for 814 wild *Mus musculus*, *Mus spretus* and *Mus spicilegus*.

## File listing ##
* `merged.renamed.vcf.gz`: genotype matrix in VCF format, 57945 sites x 814 samples
* `final_sample_table.tsv`: sample metadata in tab-separated text format. `NA` = missing value.

## Data dictionary ###
| column    | definition                                                                                                                                                                                               |
|-----------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| iid       | unique individual identifier                                                                                                                                                                             |
| taxon     | taxon label; dom = *M. m. domesticus*, mus = *M. m. musculus*, cas = *M. m. castaneus*, gen = *M. m.* "undefined", hybrid = inter-subspecies hybrid, spretus = *M. spretus*, spicilegus = *M. spicilegus* |
| sex       | chromosomal sex; F = female, M = male                                                                                                                                                                    |
| dipnum    | diploid chromosome number                                                                                                                                                                                |
| race      | chromsomal race, as 4-letter code                                                                                                                                                                        |
| country   | country where sample was collected, as ISO 2-letter code                                                                                                                                                 |
| locale    | region where sample was collected, as 4-letter code                                                                                                                                                      |
| site      | geographic place name where sample was collected                                                                                                                                                         |
| lat       | latitude where sample was collected, approximate                                                                                                                                                         |
| long      | latitude where sample was collected, approximate                                                                                                                                                         |
| pca_label | geographic grouping to which this sample is assigned in PCA plots                                                                                                                                        |
| is_core   | if TRUE, this sample is part of the subset of unrelated, geographically non-redundant samples                                                                                                            |
| cluster   | for *M. m. domesticus* samples, continental-level cluster to which this sample is assigned                                                                                                               |
| platform  | genotyping platform; mm = MegaMuga (78K) array, giga = GigaMUGA (143K array), wgs = whole-genome sequencing                                                                                              |
| f_A       | estimated inbreeding coefficient ($\hat{F}_{ROH}$) on autosomes                                                                                                                                          |
| f_X       | estimated inbreeding coefficient ($\hat{F}_{ROH}$) on X chromosome, for females only                                                                                                                     |
