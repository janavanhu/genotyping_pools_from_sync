# genotyping_pools_from_sync
This script requires a sync file as an input and calls both variant and non-variant sites based on user specified requirements from pooled samples and generates a new sync file for both all and variant only sites. 
For variant sites pseudo-haplotypes are generated based on observed allele frequencies and a genepop outpup is created.

The usage is as follows:

python call_pooled_sync.py input.sync minT maxTall maxTpop nPops minPops minAF minCount

to be used within the folder where the script is stored. The user is required to provide values for the following mandatory inputs:

input.sync - the sync file to be used as an input (provide the full path if the file is not wihtin the same the directory as the python script)
minT - minimum coverage for a site to be considered in any population
maxTall - maximum coverage allowed across all populations
maxTpop - maximum coverage allowed within any populations
nPops - the number of population samples in the input.sync file
minPops - the minimum number of populations a site needs to be covered in to be considered
minAF - the minimum allele frequency a site is required to show in a population to be considered
minCount - the minimum count a variant needs to have across populatons to be considered

three output files will be generated in the same folder as the input file:
a new sync file with all varibale and non-variable sites passing user specified requirements (extension: _All_sites.sync).
a new sync file with all varibale sites passing user specified requirements (extension: _All_variable_sites.sync).
a genepop file containing pseudo-haplotypes based on observed allele frequencies (extension: .genepop).


This script was written in python 3.6.3
