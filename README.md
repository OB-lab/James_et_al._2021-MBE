
This repository contains step-by-step instructions, code and data for:
**James ME, Arenas-Castro H, Groh JS, Allen SL, Engelstaedter J and Ortiz-Barrientos D. (2020) *Highly replicated evolution of parapatric ecotypes* [journal] [volume.page] [doi]**

# DNA extraction and Genotyping-by-Sequencing

See [CTAB_protocol](laboratory/CTAB_protocol.doc) for the CTAB DNA extraction protocol, and [GBS_protocol](laboratory/GBS_protocol.doc) for the Genotyping-by-Sequencing protocol.

# Bioinformatics

## Quality filtering and trimming

We received de-multiplexed forward and reverse sequencing files for each individual from The Beijing Genomics Institute (BGI). BGI also removed forward barcodes and quality filtered the raw reads to remove reads containing Illumina adaptors, low quality reads (>50% of bases <Q10), and reads with >10% Ns.

## Read alignment

We used ```TagCleaner``` to remove reverse barcodes from each individual. ```-tag5``` is the reverse barcode sequence (in this example GTCA). ```-mm5``` is the maximum number of mismatches, in this case 1 (we allowed a mismatch of 2 for barcodes 7 nucleotides in length). ```-trim_within``` is the number of bp from the end of the sequence to search for the barcode, in this case 5 (we used one more than the length of the barcode). 

```
perl tagcleaner.pl -fastq ind1_2.fq -tag5 GTCA -mm5 1 -trim_within 5 -out ind1_trim_2 -log ind1_trim_2.log
```

The reference genome was indexed with ```BWA```.

```
bwa index -a is reference.fasta
```

For each individual, we aligned reads to the reference genome with the ```BWA-MEM``` algorithm. The _1 and _2 files correspond to the forward and reverse read files for each individual.

```
bwa mem -M reference.fasta ind1_1.fastq ind1_trim_2.fastq >ind1.sam 2>ind1_sam.log
```

The *sam* files were converted to *bam* files with ```SamTools```.

```
samtools view -bS ind1.sam > ind1.bam
```

```SamTools``` was used to produce a summary of alignment stats per individual (including the number and % of reads mapped to the reference).

```
samtools flagstat ind1.bam &> ind1_stats.txt
```

```PicardTools``` was used to clean bam files (to set soft-clipping for reads beyond end of the reference, and MAPQ to 0 for unmapped reads). PCR duplicates were not marked for removal. 

```
java -jar picard.jar CleanSam INPUT=ind1.bam OUTPUT=ind1.cleaned.bam 2>ind1.cleaned.bam.log
```

```PicardTools``` was used to add read groups to each individual. ```RGID``` is the read group ID, (i.e. the name of the individual). ```RGLB``` is the read group library, (i.e. the sequencing lane number). ```RGPU``` is the read group platform unit (i.e. run barcode of the lane). ```RGSM``` is the read group sample name (i.e. the name of the individual).

```
java -jar picard.jar AddOrReplaceReadGroups INPUT=ind1.cleaned.bam OUTPUT=ind1.sort.rg.bam SORT_ORDER=coordinate RGID=ind1 RGLB=lib1 RGPL=ILLUMINA RGPU= HC2TFBBXX:2 RGSM=ind1 CREATE_INDEX=True 2>ind1.sort.rg.log
```

## SNP calling

We used ```FreeBayes``` to jointly call SNPs per population. Joint calling uses information across all samples to call a SNP (which can be beneficial for samples with low coverage). It thus assumes all samples are genetically similar. It is typical to joint call on “cohorts” of genetically similar individuals (aka populations or species). SNPs from these cohorts are then combined into one VCF file. 

We jointly called SNPs independently for each of the 23 populations because **1)** our samples were collected within distinct populations that occupy discrete geographic ranges, and **2)** previous data shows that individuals within a population cluster together (are quite genetically similar). 

However, when jointly calling SNPs, you must be careful with how the variant caller outputs the data. Typically, variant-calling software only outputs variant sites. However, if SNPs are being called on separate populations, and then combined, you will be unable to distinguish between monomorphic sites, and those that are missing data. Therefore, it is important to output all variant and invariant sites when you have multiple populations that are independently being jointly called. For instance, the ```--report-monomorphic``` flag within ```FreeBayes``` achieves this.

Due to computational constraints, we excluded contigs with extremely high coverage (removal of a contig if any site was >1,000X for an individual). This was an arbitrary cut-off value. Even if we called SNPs within these high coverage regions, these sites would have been removed in the downstream filtering steps due to their high coverage.

We ran ```FreeBayes``` for each of the 23 populations. The code below shows an example with two populations.

```
./freebayes -f reference.fasta ind1.sort.rg.bam ind2.sort.rg.bam --use-best-n-alleles 4 --report-monomorphic --genotype-qualities  –targets high_coverage_regions.bed >pop1.vcf
```

## Normalisation and merging of SNP cohorts

Because SNPs were jointly called per population, it is crucial to normalise these files before merging them together.

The first step of the normalisation process was to split multiallelic sites into biallelic records for each jointly-called VCF file with ```BCFtools```. This was to ensure each that each VCF file was recording multiallelic sites in the same way (see next step).

```
bcftools norm pop1.vcf -m +any -o pop1_unnorm.vcf
```

Each file was then normalised by re-joining biallelic sites into multiallelic records. This was to ensure that the multiallelic records were consistently recorded for each population.

```
bcftools norm pop1_unnorm.vcf -m +any -o pop1_unnorm_norm.vcf
```

For each population, indels were left-aligned and normalised. This was to ensure complex regions (such as indels or short tandem repeats) were normalised the same way for each population. See: https://genome.sph.umich.edu/wiki/Variant_Normalization

```
bcftools norm pop1_unnorm_norm t.vcf -f reference.fasta -o pop1_unnorm_norm_align.vcf
```

VCF output files from ```FreeBayes``` sometimes contain biallelic block substitutions. For instance, instead of having one SNP per line, blocks of multiple SNPs sometimes exist. As we have a separate VCF file per population, different biallelic blocks may exist for each population. If these are merged without decomposing them, repeated SNPs will exist within the file (a variant may be present within a biallelic block, but also on a separate line). Therefore, we used ```vt``` to decompose biallelic block substitutions into separate SNPs for each population.

```
vt decompose_blocksub -p pop1_unnorm_norm_align.vcf -o pop1_unnorm_norm_align_decomp.vcf
```

Each VCF file was zipped and then indexed.

```
bgzip -c pop1_unnorm_norm_align_decomp.vcf > pop1_unnorm_norm_align_decomp.vcf.gz  ;  tabix -p vcf pop1_unnorm_norm_align_decomp.vcf.gz
```

All VCF files were merged into one joint file. The example below is for two populations

```
bcftools merge pop1_unnorm_norm_align_decomp.vcf.gz pop2_unnorm_norm_align_decomp.vcf.gz -o all_joint.vcf
```

## SNP filtering

We closely followed the ```dDocent``` guidelines for SNP filtering: http://ddocent.com/filtering/. Some of the below is directly from the ```dDocent``` pipeline, so kudos to J. Puritz!

Using ```VCFtools```, we first kept variants genotyped in 50% of individuals that have a minimum quality score of 30, and a minor allele count of 1.

```
vcftools --gzvcf all_joint.vcf.gz --max-missing 0.5 --mac 1 --minQ 30 --recode --recode-INFO-all --stdout | gzip -c > all_joint.Q30mac1.vcf.gz
```

We applied a minimum depth per sample of 3 reads (genotypes with less than 3 reads were recoded as missing data).

```
vcftools --gzvcf all_joint.Q30mac1.vcf.gz --minDP 3 --recode --recode-INFO-all --stdout | gzip -c > all_joint.Q30mac1dp3.vcf.gz
```

To remove low-quality individuals (those with a lot of missing data), we first created a list of the proportion of missing data per individual. The following code creates an output file called *out.imiss*, with the fifth column containing the proportion of missing data. 

```
vcftools --gzvcf all_joint.Q30mac1dp3.vcf.gz --missing-indv
```

We graphed the distribution of missing data per individual in R.  

![Alt text](images/missing_data_all.jpeg?raw=true "Title")

We removed the upper tail of the distribution, discarding all individuals with more than 40% missing data. This was achieved by first creating a list of individuals with >40% missing data. 

```
mawk '$5 > 0.4' out.imiss | cut -f1 > lowDP.indv
```

We used this list in VCFtools to remove the low-quality individuals.

```
vcftools --gzvcf all_joint.Q30mac1dp3.vcf.gz --remove lowDP.indv --recode --recode-INFO-all --out all_joint.Q30mac1dp3ir
```

Sites with high coverage could be multi-copy regions of the genome, so we wanted to remove these potential paralogues. We first examined the distribution of mean read depth across all sites by calculating the mean depth per site and graphing the distribution in R. 

```
vcftools --vcf all_joint.Q30mac1dp3ir.recode.vcf --site-mean-depth --out mean_depth
```

![Alt text](images/mean_read_depth_all.jpeg?raw=true "Title")

In general, the mean read depth per locus should be approximately normally distributed. Within the literature, various approaches have been used to select the maximum mean read depth. For instance, Li (2014) suggests using the equation: **d+3*sqrt(d)**, d=mean depth (which is a value of 63 for our dataset). However, this method has been suggested as too conservative for RADseq data. Others use the 90th quantile (which is 88 for our data), two times the mode (40 for our data, after rounding the mean depths to the nearest 10), or, others just eyeball the mean depth distribution and remove the upper tail (~120 for our data). After examining these multiple approaches, we chose a maximum mean depth of 100. 

```
vcftools --vcf  all_joint.Q30mac1dp3ir.recode.vcf --recode --recode-INFO-all --out all_joint.Q30mac1dp3irMaxDP100 --max-meanDP 100 
```

We then filtered for a minimum mean depth of 10. 

```
vcftools --vcf all_joint.Q30mac1dp3irMaxDP100.recode.vcf --min-meanDP 10 --recode --recode-INFO-all --out all_joint.Q30mac1dp3irMaxDP100MinDP10
```

To ensure that each SNP is sequenced in every population, we filtered for missing data per population. We first created a file containing the individual name in the first column, and the population identifier in the second column: *all_inds_pops.txt*. For each population, the individual names were extracted into a separate file. 

```
mawk '$2 == "D00"' all_inds_pops.txt > D00_keep && mawk '$2 == "D01"' all_inds_pops_mac1.txt > D01_keep && mawk '$2 == "D02"' all_inds_pops_mac1.txt > D02_keep
```

For each of the 23 populations we used VCFtools to calculate the proportion of missing data per variant site. The below lines of code are an example for three populations. They produce output files called *D00.lmiss*, *D01.lmiss* and *D02.lmiss*, the 6th column containing the proportion of missing data for each site. 

```
vcftools --vcf all_joint_FB.Q30mac1dp3irMaxDP100MinDP10.recode.vcf --keep D00_keep --missing-site --out D00
vcftools --vcf all_joint_FB.Q30mac1dp3irMaxDP100MinDP10.recode.vcf --keep D01_keep --missing-site --out D01
vcftools --vcf all_joint_FB.Q30mac1dp3irMaxDP100MinDP10.recode.vcf --keep D02_keep --missing-site --out D02
```

For each population file, we combined all loci that had greater than 50% missing data into one file: *badloci0.5*.

```
cat D00.lmiss D01.lmiss D02.lmiss | mawk '!/CHR/' | mawk '$6 > 0.5' | cut -f1,2 >> badloci0.5
```

We used ```VCFtools``` to remove these loci from the VCF file.

```
vcftools --vcf all_joint_FB.Q30mac1dp3irMaxDP100MinDP10.recode.vcf --exclude-positions badloci0.5 --recode --recode-INFO-all --out all_joint_FB.Q30mac1dp3irMaxDP100MinDP10mpp50
```

We additionally filtered for an overall missing data, removing sites if they had greater than 20% missing data. 

```
vcftools --vcf all_joint_FB.Q30mac1dp3irMaxDP100MinDP10mpp50.recode.vcf --recode --recode-INFO-all --max-missing 0.8 --out all_joint_FB.Q30mac1dp3irMaxDP100MinDP10mpp50md80
```

We then removed indels, by first converting variant calls to separate SNP and indel genotypes with ```vcflib```.

```
vcfallelicprimitives all_joint_FB.Q30mac1dp3irMaxDP100MinDP10mpp50md80.recode.vcf --keep-info --keep-geno > all_joint_FB.Q30mac1dp3irMaxDP100MinDP10mpp50d80prim.vcf
```

Indels were then removed.

```
vcftools --vcf all_joint_FB.Q30mac1dp3irMaxDP100MinDP10mpp50d80prim.vcf --remove-indels --recode --recode-INFO-all --out all_joint_FB.Q30mac1dp3irMaxDP100MinDP10mpp50md80ind
```
		
We filtered for Hardy Weinberg Equilibrium within each population using the ```dDocent``` *filter_hwe_by_pop.pl* custom script, found here: https://github.com/jpuritz/dDocent/raw/master/scripts/filter_hwe_by_pop.pl. 

```
./filter_hwe_by_pop.pl -v all_joint_FB.Q30mac1dp3irMaxDP100MinDP10mpp50md80ind.recode.vcf -p all_inds_pops.txt -o all_joint_FB.Q30mac1dp3irMaxDP100MinDP10mpp50md80indHWE
```

(This file was renamed to: *all_rel_50pp_80md_HWE.vcf*)

We also explored variations of the above SNP filtering parameters, ranging from ‘relaxed’ to ‘intermediate’ to ‘stringent’ parameters. We plotted PCAs to examine whether the clustering of populations was robust to variations in the way SNPs were filtered. We also explored further filtering steps outlined in the ```dDocent``` SNP filtering pipeline (http://ddocent.com/filtering/). We found that populations consistently grouped into the same clusters, suggesting our data is robust to slight variations in SNP quality and SNP numbers. Within our paper (the code presented above), we have used an ‘intermediate’ level of filtering.

However, despite the normalisation process before merging the per-population VCF files (outlined above), a few duplicate SNPs remined within the final VCF file. These were removed. First, the chromosome number and variant position were extracted from: *all_rel_50pp_80md_HWE.vcf*, with a colon separating them. 

```
bcftools query -f '%CHROM:%POS\n' all_rel_50pp_80md_HWE.vcf > all_rel_50pp_80md_HWE_stats.txt
``` 

The duplicate sites were identified. 

```
perl -ne 'print if $seen{$_}++' all_rel_50pp_80md_HWE_stats.txt > dups.txt
```
		
```PLINK``` was used to set the variant IDs from '.' to the chromosome and position (separated by a colon). We also processed VCF half-calls as missing data. Half-calls occur when one allele for an individual at a particular site is called with high confidence, but the other is not, leaving that site with only one allele. To be conservative, we treated the whole SNP for that individual as missing. 
	
```	
./plink2 --vcf all_rel_50pp_80md_HWE.vcf --export vcf --out all_rel_50pp_80md_HWE_ids --vcf-half-call m --allow-extra-chr --set-all-var-ids @:#
```
	
The duplicate SNPs were removed from the VCF file. 

```
./plink2 --vcf all_rel_50pp_80md_HWE_ids.vcf --allow-extra-chr --exclude dups.txt --export vcf --out all_rel_50pp_80md_HWE_ids_removed
```

We then re-filtered for overall missing data of 20%.
	
```	
./plink --vcf all_rel_50pp_80md_HWE_ids_removed.vcf --geno 0.2 --export vcf --allow-extra-chr --out all_rel_50pp_80md_HWE_ids_removed_80md 
```

We filtered for an overall minor allele frequency of 0.05.

```
./plink --vcf all_rel_50pp_80md_HWE_ids_removed.vcf --maf 0.05 --export vcf --allow-extra-chr --out all_rel_50pp_80md_HWE_MAF0.05
```

## Selection of unlinked SNPs

In some datasets, only unlinked SNPs were used. To obtain unlinked SNPs (~one SNP per rad tag), we used ```PLINK``` to retain one SNP per 2000bp. 

```
./plink --vcf all_rel_50pp_80md_HWE_MAF0.05.vcf --make-bed --bp-space 2000 --allow-extra-chr --vcf-half-call m  --export vcf --out all_rel_50pp_80md_HWE_MAF0.05_unlinked
```

## The final datasets we used within our analyses are as follows:

* All populations, MAF 0.05, unlinked SNPS: [all_rel_50pp_80md_HWE_MAF0.05_unlinked.vcf.gz](vcf_files/all_rel_50pp_80md_HWE_MAF0.05_unlinked.vcf.gz)
* All populations, MAF 0.01: [all_rel_50pp_80md_HWE_MAF0.01.vcf.gz](vcf_files/all_rel_50pp_80md_HWE_MAF0.01.vcf.gz)
* Western Australia populations removed, MAF 0.05: [ESC_rel_50pp_80md_HWE_MAF0.05.vcf.gz](vcf_files/ESC_rel_50pp_80md_HWE_MAF0.05.vcf.gz) (Note: these Western Australia populations were removed before filtering, when all populations were merged into one joint file)
* Western Australia populations removed, MAF 0.05, unlinked SNPs: [ESC_rel_50pp_80md_HWE_MAF0.05_unlinked.vcf.gz](vcf_files/ESC_rel_50pp_80md_HWE_MAF0.05_unlinked.vcf.gz)
* Western Australia populations removed, MAC 1: [ESC_rel_50pp_80md_HWE_MAC1.vcf.gz](vcf_files/ESC_rel_50pp_80md_HWE_MAC1.vcf.gz)

# Do populations cluster by geography or by ecotype?

We used ```IQ-TREE``` to generate a maximum likelihood phylogeny with using the dataset: [all_rel_50pp_80md_HWE_MAF0.05_unlinked.vcf.gz](vcf_files/all_rel_50pp_80md_HWE_MAF0.05_unlinked.vcf.gz). We used the polymorphisms-aware phylogenetic model. We first used ```PGDspider``` to convert the VCF file to a fasta file. ```FastaToCounts.py``` (https://github.com/pomo-dev/PoMo/blob/master/scripts/FastaToCounts.py) was then used to convert the fasta file to a counts file. This [IQtree_renamed.counts](IQ-TREE/IQtree_renamed.counts) file summarises the allele frequencies for each population and is the input file to ```IQ-TREE```. 

```
FastaToCounts.py all_rel_50pp_80md_HWE_MAF0.05_unlinked_renamed.fasta.gz IQtree_renamed.counts
```

We then used ```ModelFinder``` to determine the best-fit substitution model for the data. The D09 population from Western Australia was assigned as the outgroup. 

```
iqtree -s IQtree_renamed.counts -m MF -o D09 -pre ModelFinder
```

*TVMe+FQ+P+N9+G4* was the best substitution model, and we increased the virtual population size (*N*) to maximum value of 19, as recommended by Schrempf et al. (2016). We then constructed the phylogeny. Branch support was performed using ```UFboot``` (10,000 replicates) and ```SH-aLRT``` (10,000 replicates).

```
iqtree -s IQtree_renamed.counts -m TVMe+FQ+P+N19+G4 -o D09 -bb 10000 -alrt 10000 -pre run1
```

To assess convergence, we undertook 10 separate runs of above ```IQ-TREE``` code and examined tree topology (which remained unchanged with 10 independent runs). We also ensured that the log-likelihood values were stable at the end of each run. 


We then used ```fastSTRUCTURE``` to explore population structure across all populations. We used the dataset: [all_rel_50pp_80md_HWE_MAF0.05_unlinked.vcf.gz](vcf_files/all_rel_50pp_80md_HWE_MAF0.05_unlinked.vcf.gz). We used ```PGDspider``` as well as manipulations in excel to convert the VCF file into ```fastSTRUCTURE``` format, see [fastSTRUCTURE.str](fastSTRUCTURE/fastSTRUCTURE.str). We ran the simple prior (K=1-30) with 100 independent runs per K-value. ```choosek.py``` (https://github.com/rajanil/fastStructure/blob/master/chooseK.py) was used to choose the most likely number of genetic clusters. Results were summarized and plotted in the R package ```pophelper```, by following the tutorial here: http://www.royfrancis.com/pophelper/articles/index.html 

# Has gene flow shaped patterns of divergence across the system?

We used ```TreeMix``` to explore patterns of gene flow in a phylogenetic context. We used the dataset: [all_rel_50pp_80md_HWE_MAF0.01.vcf.gz](vcf_files/all_rel_50pp_80md_HWE_MAF0.01.vcf.gz). To convert this VCF file into the ```TreeMix``` format, we used ```PLINK``` to make bed, bim and fam files. 

```
./plink --vcf all_rel_50pp_80md_HWE_MAF0.01.vcf --allow-extra-chr --make-bed --out all_rel_50pp_80md_HWE_MAF0.01
```

This was converted to a frequency file.

```
./plink --bfile all_rel_50pp_80md_HWE_MAF0.01 --allow-extra-chr --freq --missing --within all_MAF0.01.clusters.txt
```

[all_MAF0.01.clusters.txt](TreeMix/all_MAF0.01.clusters.txt) specifies the individuals in the first two columns, and the population code in the third column (in the order of the VCF file). The output file from above was zipped. 

```
gzip plink.frq.strat
```

The ```plink2treemix.py``` (https://github.com/ekirving/ctvt/blob/master/plink2treemix.py) script was used to convert the zipped file into ```TreeMix``` format, which was used as the input file: [TreeMix.frq.gz](TreeMix/TreeMix.frq.gz) 

```
python plink2treemix.py plink.frq.strat.gz TreeMix.frq.gz
```

We constructed an initial 25 maximum likelihood trees with no migration, 1,000 bootstrap replicates in blocks of 50 SNPs with D09 as the assigned outgroup. For instance, for one replicate:

```
treemix -i TreeMix.frq.gz -bootstrap 1000 -k 50 -root D09 -m 0 -o TreeMix_m0_1
```

We selected the tree with the highest log-likelihood (which was replicate no. 19) as the input tree for all subsequent analyses. We then tested between 1-25 migration events in blocks of 50 SNPs. For instance, for one migration event:

```
treemix -i TreeMix.frq.gz -g TreeMix_m0_19.vertices.gz TreeMix_m0_19.edges.gz -k 50 -m 1 -o TreeMix_m1
```

To select the number of migration events, we examined the log-likelihoods and cumulative variance explained by each model by extracting the log likelihood and variance from each of the 25 migration events and plotting the distribution in R. We also performed jackknife estimates to obtain the standard error and significance of the weight of each migration event. 

```
treemix -i 1 TreeMix.frq.gz -g TreeMix_m0_19.vertices.gz 1 TreeMix _m0_19.edges.gz -se -m 25 -o TreeMix_jackknife_m25
```

We also more formally tested for genetic admixture, we calculated *f3* in ```TreeMix``` for all triads of populations with jackknifing in blocks of 50 SNPs.

```
threepop -i TreeMix.frq.gz -g TreeMix _m0_19.vertices.gz TreeMix_m0_19.edges.gz -k 50
```

We then tested for isolation by distance (IBD) using migration rates inferred from ```fastsimcoal2``` (see below for details), as well as Slatkin’s M. For ```fastsimcoal2``` we tested for IBD between the Dune and Headland ecotypes at each locality. See [IBD_fsc.txt](IBD/IBD_fsc.txt) for input file. We performed a linear model in R using the average gene flow rate between the Dune and Headland at each locality.

```
IBD_fsc <- read.delim ("path/to/file/IBD_fsc.txt", header=T)
summary(lm(GFavrg~GeodistLOGkm, data=IBD_fsc))
```

For Slatkin’s M, we used the dataset: [ESC_rel_50pp_80md_HWE_MAF0.05.vcf.gz](vcf_files/ESC_rel_50pp_80md_HWE_MAF0.05.vcf.gz). We calculated pairwise F<sub>ST</sub> in ```VCFtools```, for instance, for population D00 and D01:

```
vcftools --vcf ESC_rel_50pp_80md_HWE_MAF0.05.vcf --weir-fst-pop pops/D00.txt --weir-fst-pop pops/D01.txt --out D00-D01
```

Where D00.txt and D01.txt are files specifying the individuals within each population. We calculated Slatkin’s M using the formula: (1 / F<sub>ST</sub> - 1) / 4

Using Slatkin’s M, we tested for IBD within the Dunes and also within the Headlands using Mantel tests in R. See [Dgendist.txt](IBD/Dgendist.txt) for the Dune matrix of pairwise genetic distances (log-scale), [Dgeodist.txt](IBD/Dgeodist.txt) for the Dune matrix of pairwise geographic distances (log-scale), [Hgendist.txt](IBD/Hgendist.txt)for Headland genetic distances, and [Hgeodist.txt](IBD/Hgeodist.txt) for Headland geographic distances. For instance, for the Dunes in R: 

```
library(vegan)
Dgen<- read.delim ("path/to/file/Dgendist.txt", header=F)
Dgeo<- read.delim ("path/to/file/Dgeodist.txt", header=F)
DgenM <- as.dist(Dgen)
DgeoM <- as.dist(Dgeo)
#mantel test
mantel(DgeoM, DgenM, permutations = 9999)
```

We also performed a linear model in R using the average gene flow rate between the Dune and Headland at each locality. See [IBD_sm.txt](IBD/IBD_sm.txt) for input file.

```
IBD_sm <- read.delim ("path/to/file/IBD_sm.txt", header=T)
summary(lm(FstSkMeanLOG~GeodistLOGkm, data=IBD_sm))
```

# Is there gene flow between parapatric populations?

```STRUCTURE``` was used to estimate levels of admixture between ecotypes at each locality. We used the dataset: [ESC_rel_50pp_80md_HWE_MAF0.05_unlinked.vcf.gz](vcf_files/ESC_rel_50pp_80md_HWE_MAF0.05_unlinked.vcf.gz). We then extracted each population pair and removed SNPs with MAF < 0.05 per pair and used ```PGDspider``` to convert each VCF file into ```STRUCTURE``` format. ```STRUCTURE``` was run using the admixture model and the correlated allele frequency model with 10 independent runs for K=1-6 (50,000 burn-in and 200,000 MCMC iterations). See [STRUCTURE](STRUCTURE) for all input files as well as example ```mainparams``` and ```extraparams``` files. Results were summarized and plotted in the R package ```pophelper```, by following the tutorial here: http://www.royfrancis.com/pophelper/articles/index.html


# Gene flow detection with fastsimcoal2

```fastsimcoal2``` (available at <http://cmpg.unibe.ch/software/fastsimcoal2/>) is a continuous-time coalescent simulator of genomic diversity under arbitrarily complex evolutionary scenarios. It can estimate demographic parameters from the site frequency spectrum through a composite likelihood maximisation procedure. 

We used ```fastsimcoal2``` to test 10 demographic models (see below) in 30 population pairs of interest given their ecology, phylogenetic relationships, and occurrence patterns:

  + Dune-Headland pairs: D00-H00, D01-H01, D03-H02, D04-H05, D05-H06, D12-H14, D14-H15, D32-H12.
  + Dune-Dune pairs: D00-D02, D01-D03, D01-D04, D02-D03, D04-D05, D05-D12, D12-D14, D14-D32.
  + Headland-Headland pairs: H00-H02, H01-H04, H01-H05, H02-H04, H03-H07, H03-H14, H05-H06, H06-H07, H12-H12A, H12-H15, H14-H15.
  + Allopatric pairs: D03-D32, D03-H12, H02-H12.

The ten demographic models were:

  + M1: No migration.
  + M2: Bidirectional migration.
  + M3: Population 2 to population 1 migration.
  + M4: Population 1 to population 2 migration.
  + M5: Bidirectional migration after secondary contact. 
  + M6: Population 2 to population 1 migration after secondary contact.
  + M7: Population 1 to population 2 migration after secondary contact.
  + M8: bidirectional migration after population splitting with cessation of gene flow.
  + M9: Ancient bidirectional secondary contact 
  + M10: Population size reduction after bidirectional secondary contact.

To run each model, ```fastsimcoal2``` requires three input files:

  + A site frequency spectrum file.
  + A *template* file.
  + An *estimation* file.

All the ```fastsimcoal2``` input files used in this study are available at [James_et_al._2021-MBE-fastsimcoal2_files.zip](fastsimcoal2/James_et_al._2021-MBE-fastsimcoal2_files.zip)

Please note that the header of the template file for models 3, 4, 6, and 7 makes reference to the migration direction backward in time whereas in the paper both the name of the models and migration rates are presented forward in time. Also, population 1 and population 2 follow the order assigned in the population pair name. For instance, for the population pair D00-H00, D00 is refered here as population 1 and H00 as population 2. Population 1 corresponds to the rows of the SFS file and is labeled as DUNE in the template and estimation files. Population 2 corresponds to the columns of the SFS file and is labeled as HEAD in the template and estimation files. Both sample and population sizes are in haploid numbers. 

The rest of this document uses the population pair D00-H00 as an example to name the files and command arguments.

## Getting the site frequency spectrum (SFS) file

The SFS is a summary of genome-wide data describing the distribution of allele frequencies in a sample of one or more populations. For instance, it tells how many SNPs are derived in only one chromosome, two chromosomes, three chromosomes, and so on in the population. Its shape is influenced by the history of the population: migration, population size changes, substructure, etc. The SFS treats all SNPs in the data set as independent of one another. 

A SFS is referred as *folded* when the information about the ancestral/derived state of the SNP is unavailable. Instead, the minor allele frequency is used as criterion for assigning an ancestral/derived-like state. A SFS is referred as *joint* when it summarises information from two or more populations. 

Since the SFS strictly only considers SNPs without missing data, many SNPs can be lost when all the individuals in a VCF file are used to generate the SFS. To avoid this, the data can be downsampled to the number of chromosomes (haploid samples) that maximises the number of SNPs without missing data. The amount of missing data heavily relies on the quality of the sequence (sequencing depth). Since we were interested in retaining rare alleles, we didn't downsample the VCF files after checking that the number of SNPs didn't drastically drop due to missing data.

To generate the SFS file, we used the functions *vcf2dadi* and *dadi2fsc.2pop* from the R script ```vcf2sfs.r``` (April 13th, 2016 version by Shenglin Liu. Repository available at <https://github.com/shenglin-liu/vcf2sfs>). Along the VCF file, ```vcf2dadi``` requires a tab-delimited population specification file. It contains the sample names (as specified in the VCF file) in the first column and the corresponding population name in the second column. This file can be directly generated from the VCF file using the custom Perl script ```getpopmap.pl```. It works for files containing either two or three populations.

```
perl getpopmap.pl D00_H00.vcf D00_H00.popmap D00 H00
```

The first argument is the name of the input VCF file, the second argument is the name of the output popmap file, and the following two or three arguments are the names of the populations. Please note that the names of the samples in the VCF file should start with the population name followed by a hyphen. *vcf2dadi* and *dadi2fsc.2pop* were run sequentially in R. Note that the SFS file should be named *D00_H00_jointMAFpop1_0.obs* to be called by ```fastsimcoal2```.

```
source("[PATH]/vcf2sfs.r")
vcf2dadi("D00_H00.vcf", "D00_H00.popmap", "D00_H00.dadi", c("D00", "H00"), ploidy=2, n.digit=4, filter.indi=NA, filter.snp=0)
dadi2fsc.2pop("D00_H00.dadi", "D00_H00_jointMAFpop1_0.obs", c("D00", "H00"), fold=T)
```

## Getting the template file

The *template* file specifies the evolutionary model and the parameters that should be estimated. It is divided into ten fixed sections in the following order:

  + Number of populations.
  
  + Population effective sizes: It should contain one line per population. If the effective population size is unknown, a parameter name should be specified instead to be calculated. The population in the first line should correspond to the one in the columns of the SFS file and the population in the second line to the one in the rows of the SFS file. The order of the populations in the SFS file is alphabetically sorted based on the population names. The program indexes the first population as 0 and the second one as 1.
  
  + Sample sizes: It should contain one line per population in the order specified in the previous section.
  
  + Growth rates: It should contain one line per population. A zero growth rate means stationary population size. Negative growth implies population expansion since it is estimated backward in time.
  
  + Number of migration matrixes.
  
  + Migration matrix: It should contain one matrix per migration matrix, as mentioned in the previous section. The program indexes the first migration matrix as 0, the second one as 1, and so on. A cell value in the migration matrix is read as the probability that a haploid individual migrates backward in time from the population specified in the row to the population specified in the column. If the migration rate is unknown, a parameter name should be specified instead in the corresponding cell to be calculated.
  
  + Historical events: This section is the heart of the evolutionary model. The first line specifies the number of historical events, namely gene flow, bottlenecks, admixture events, coalescence of populations, etc. Each historical event should be specified in a different line and be defined by 7 parameters in the following order. If any of those values are unknown, a parameter name should be specified instead to be calculated.
    
    - Number of generations in the past at which the historical event happened.
    - Source population.
    - Sink population.
    - Proportion of migrants moving from source to sink population. It should be 1 if both populations coalesce at this event.
    - New relative size for the sink population backward in time.
    - New growth rate for the sink population backward in time.
    - New migration matrix index to be used further back in time.
  
  + Number of independent loci: The number of independent chromosomes to be simulated and a flag indicating if the different chromosomes have a different (1) or a similar (0) structure. 
  
  + Number of linkage blocks: The number of blocks per chromosome that may differ by the type of markers to be simulated, the recombination rate, or the mutation rate.
  
  + Genetic properties: Data type, number of markers, recombination rate, and mutation rate. FREQ in the data type field specifies to estimate the SFS with the simulations.
  
Every section is introduced by a comment line starting with the characters ```//```. The name of the *template* file must be the root name of the SFS file plus the extention *.tpl*. For instance, the corresponding *template* file name of *D00_H00_jointMAFpop1_0.obs* SFS should be *D00_H00.tpl*.

This is how a template file looks like:

```
// Two populations - M5 Bidirectional secondary contact
2
//Population effective sizes: to be estimated.
HEAD
DUNE
//Samples haploid size: the same one defined in the projections. 
126
124
//Growth rates: negative growth implies population expansion
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//Migration matrix 0: to be estimated.
0 H2D
D2H 0
//Migration matrix 1
0 0
0 0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index
2 historical event
TDIV 0 1 1 RESIZE 0 1
TSEC 0 1 0 1 0 1
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per gen recomb and mut rates
FREQ 1 0 1e-8
```

In this file, *HEAD*, *DUNE*, *H2D*, *D2H*, *TDIV*, *RESIZE*, and *TSEC* are the unknown parameters to be estimated by ```fastsimcoal2```. They should be specified in the corresponding *estimation* file as well.

## Getting the estimation file

The *estimation* file specifies the search range of the parameters defined in the *template* file. The lower range limit is an absolute minimum, whereas the upper range is only used as a maximum for choosing a random initial value for the parameter. It is divided into three fixed sections in the following order:
  
  + Parameters: It lists the main parameters and its respective initial search range. Each parameter can either be uniformly or log-uniformly distributed and its values can be either recorded or omitted in the output files. Some parameters, such as *ANCSIZE*, are specified in this section despite not being invoked in the *template* file. This is because they are later used to calculate other complex parameters that were specified in the *template* file, such as *RESIZE* (see below).
  
  + Rules: It includes a set of conditions to be met among the above simple parameters. For instance, in a secondary contact model, it should be specified that the secondary contact event happened after the divergence of the involved populations. 
  
  + Complex parameters: It defines parameters that are obtained as simple operations between other parameters. For instance, the migration rate *H2D* is calculated as the number of migrants *HtoD* divided by the population size of the source population *HEAD*.
  
The name of the *estimation* file must be the root name of the SFS file plus the extension *.est*. For instance, the corresponding *estimation* file name of the *D00_H00_jointMAFpop1_0.obs* SFS should be *D00_H00.est*.

This is how an *estimation* file for a population pair looks like:

```
// Priors and rules file
// *********************

[PARAMETERS]
//#isInt? #name   #dist.#min  #max
//all N are in number of haploid individuals
1  ANCSIZE     unif     100  100000   output
1  HEAD       unif     100  100000   output
1  DUNE       unif     100  100000   output
0  HtoD       logunif  1e-2 20       hide
0  DtoH       logunif  1e-2 20       hide
1  TDIV        unif     100   20000   output
1  TSEC        unif     100   20000   output

[RULES]
TDIV > TSEC

[COMPLEX PARAMETERS]

0  RESIZE = ANCSIZE/DUNE     hide
0  H2D  = HtoD/HEAD       output
0  D2H  = DtoH/DUNE       output
```

## Running fastsimcoal

Running ```fastsimcoal2``` is easy if the three input files are well formatted. All the files should be in the same directory from where the command is invoked. 

```
fsc26 -t D00_H00.tpl -n100000 -m -e D00_H00.est -M -L60 -c10 -q
```

The ```-t``` and  ```-e``` flags specify the *template* and *estimation* files, respectively, the ```-n``` flag specifies the number of SFS simulations to fit to the data, the ```-M``` flag indicates to perform parameter estimation by maximum composite likelihood from the SFS, the ```-L``` flag specifies the number of optimisation cycles, the ```-c``` flag specifies the number of threads to be used for simulations, and the ```-q``` flag keeps output messages at the minimum level.

Reaching reliable results depends to some extent on how well the parameter space is explored. For this, every model should be independently run between 50 and 100 times and the results of the run with the maximum likelihood value should be picked up for further analysis. The custom Perl script ```mkdir_in.pl``` creates multiple directories (one for each independent run per model) and copies the input files into them before launching the independent runs.

```
perl mkdir_in.pl 10 50
```

The first argument corresponds to the number of models to test and the second argument indicates how many independent runs will be launched per model. Please note that this script should be invoked from a directory containing one folder per model named from 1 up to the total number of models. Each model directory should contain the three input files above mentioned.

Since testing the different models demands running ```fastsimcoal2``` multiple times, this is best done in a server. The custom shell executable script ```runD00_H00.sh``` executes the program 50 times per model.

```
qsub runD00_H00.sh
```

## Summarising the results

```fastsimcoal2``` generates an output directory per run, named as the root name of the input files (*i.e. D00_H00*). Within this directory, *D00_H00.bestlhoods* is particularly useful for summarising the estimated parameter values and *D00_H00.brent_lhoods* for evaluating the run performance over the optimisation cycles. 

The custom Perl script ```extract_ml.pl``` reads the *D00_H00.bestlhoods* file across all the runs per model and extracts their maximum likelihood value in a table, exported in the output file *results_D00_H00_bestlhoods.txt*. Based on those values, it identifies the maximum likelihood run per model and extracts its *\*.bestlhoods* and *\*.brent_lhoods* files in a new directory named ```results_D00_H00```. The name of those files are modified to distinguish them among models.

```
perl extract_ml.pl 10 50 D00_H00
```

The first argument corresponds to the number of models, the second argument corresponds to the number of independent runs per model, and the third argument corresponds to the root name of the files. ```extract_ml.pl``` should be invoked from the location that contains all the models in separate directories named as consecutive numbers from 1 to the maximum number of different models. 

The custom R script ```summaryplots_fsc.R``` graphically summarises the performance of the ```fastsimcoal2``` runs across the different models. It generates a PDF file showing a boxplot of the maximum likelihood values of all runs per model and a scatter plot of the likelihood values along the optimisation cycles of the best run per model. This is how the plots look like:

![Alt text](images/models_likelihoods1.png?raw=true "Title")


![Alt text](images/models_likelihoods2.png?raw=true "Title")

## Selecting the best demographic models

Information theory offers an objective way to calculate the probability of multiple demographic models and rank them accordingly, instead of only relying in the Akaike information criterion (AIC) values to decide which is the best model. This approach, namely Akaike weight, is based on the computation of the Kullback-Leibler information of every model using the AIC, followed by the normalization of these values. 

For a given model **i**, its Akaike weight **w <sub>i</sub>** is equal to **exp(−Δ <sub>i</sub> /2) / ∑ from r=1 to R of exp(−Δ <sub>r</sub> /2)**.


```Δ``` being the difference between the AIC value of a particular model and the AIC value of the best model of the set and **R** being the total number of models in the set.  The Akaike weight values of all the models should sum up 1.0. It is read as the weight of evidence in favour of the model assuming that the actual best model is present in the set of models.

## Estimating the confidence intervals

We calculated the confidence intervals for the parameters of the model that had an Akaike weight greater than 0.50 using parametric bootstrap. This approach simulates DNA sequences, and their corresponding SFS, given the chosen model and the parameter values of its best run. Then, it recalculates the parameter values from the simulated SFS. This process was done 100 times.

For simulating the DNA and SFS of the chosen model, the parameter values of its best run should be specified in a *parameter* (*\*.par*) input file. This file can be generated editing the *maximum likelihood parameter* file (*\*_maxL.par*) output by ```fastsimcoal2```. The last three sections of this file usually look like:

```
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per gen recomb and mut rates
FREQ 1 0 1e-8
```

They should be slightly modified to specify that a DNA sequence, representing a given number of independent loci (*i.e.* 10,000) of a particular length (*i.e.* 100 bp), should be simulated. After this, those sections of the new *parameter* file should look like:

```
//Number of independent loci [chromosome] 
10000 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per gen recomb and mut rates
DNA 100 0 1e-8 OUTEXP
```

Please note that after editing the D00_H00_maxL.par it should be renamed as D00_H00.par and put in a new directory along with the original SFS (*D00_H00_jointMAFpop1_0.obs*) and *template* (*D00_H00.tpl*) files. The *parameter values* file (*D00_H00.pv*) from the best run of the chosen model should be included in this new directory as well.

The custom Perl script ```runboot1.pl``` internally runs ```fastsimcoal2``` to simulate the specified number of SFS files in independent directories and then copies the  *estimation*, *template*, and *parameter values* files to the newly created directories. The *parameter values* file is used to specify the initial parameter values of the new runs. ```runboot1.pl``` should be invoqued from the newly created directory.

```
perl runboot1.pl D00_H00 100
```

The first argument corresponds to the name of the population pair and the second argument to the number of bootstrap replicates. It can take a couple of minutes to run. The custom Bash file ```runboot_D00_H00.sh``` recalculates the parameter values given the simulated SFS. 

```
qsub runboot_D00_H00.sh
```

The custom Perl script ```runboot2.pl``` extracts the parameter values from all the bootstrap runs in a text file named *D00_H00_boot.txt*.

```
perl runboot2.pl D00_H00 100
```

Please note that the order of the columns in *D00_H00_boot.txt* corresponds to the order of the output parameters specified in the *estimation* file. For instance, for the *D00_H00.est* file (model 5, as shown above) it would be: ANCSIZE, HEAD, DUNE, TDIV, TSEC, H2D, D2H.

The 95% confidence intervals of every parameter can be estimated in ```R``` environment using the function ```groupwiseMean``` of the ```rcompanion``` library. The R script ```getIC.R``` shows how to do it when model 5 is the best model.

## Calculating the number of migrants per generation

Population sizes (N) and migration rates (m) are output by ```fastsimcoal2``` in haploid numbers and backward in time, however we are interetested in calculating the number of diploid migrants from one population to other forward in time (2Nm). To do this, consider the migration rate from Headland to Dune. This is the probability of a Headland population sending haploid migrants to a Dune population backward in time. We multiplied it by the haploid population size of the Headland to get the number of haploid migrants from a Headland to a Dune backward in time. Then we halved this value to get the number of diploid migrants from a Headland to a Dune backward in time. This number is exactly as the number of diploid migrants from a Dune to a Headland population forward in time (2Nm).

For instance, for the following parater values of the D00-H00 pair:

D00 pop size (Pop1size): 47,926
H00 pop size (Pop2size): 134,364
Migration backward in time from H00 to D00 (mP2->P1): 3.24E-06
Migration backward in time from D00 to H00 (mP1->P2): 1.18E-05

The number of diploid migrants from D00 to H00 forward in time are:

```
2NmP1->P2 	= (mP2->P1 * Pop2size) / 2
		= (3.24E-06 * 134,364) /2
		= 0.2176 
```

And the number of diploid migrants from H00 to D00 forward in time are:

```
2NmP2->P1 	= (mP1->P2 * Pop1size) / 2
		= (1.18E-05* 47,926) /2
		= 0.2830
```

## Testing alternative models with very low gene flow

As ```fastsimcoal2``` uses simulations to approximate the likelihood values, there is variance in the likelihood estimates. Since we are interested in detecting gene flow between the population pairs, it is possible that our best model and an alternative model with exacly the same parameter values but without or very low gene flow are indistinguishable. To test this, we compared the likelihood distribution of the best model with the likelihood distributions of a set of alternative models with very low gene flow (2Nm = 0.01, that is one migrant every 100 generations).

For instance, when bidirectional migration after secondary contact (model 5) was the best model, there were three alternative models: either Dune to Headland migration rate (alternative model A) or Headland to Dune migration rate (alternative model B), or both (alternative model C), were constrained to be 2Nm = 0.01. To do this, either H2D or D2H values, or both, should be edited in the migration matrix of the *D00_H00_maxL.par* file. Please note that those values should correspond to the equivalent proportion of migrants sent backwards in time, not the actual number of migratns (it can be estimated as 0.01/(haploid size of the source population/2)).

To estimate the likehood distribution of the best model and the three alternative models in ```fastsimcoal2```, we need one *parameter* file per model: *D00_H00_maxL.par* (best model), *D00_H00_A.par*, *D00_H00_B.par*, and *D00_H00_C.par*. For each one of these files, we need a copy of the original SFS file named accordingly: *D00_H00_maxL_jointMAFpop1_0.obs*, *D00_H00_A_jointMAFpop1_0.obs*, *D00_H00_B_jointMAFpop1_0.obs*, and *D00_H00_C_jointMAFpop1_0.obs*. All these files should be in the same folder.

The custom Bash file ```mldist_D00_H00.sh``` iteratively runs ```fastsimcoal2``` to produce a distribution of 100 likelihood values for each model: *D00_H00.maxL.lhoods*, *D00_H00.A.lhoods*, *D00_H00.B.lhoods*, and *D00_H00.C.lhoods*. 

```
qsub mldist_D00_H00.SH
```
The custom R script ```mldist_plot.R``` read the output files and plot the likelihood distribution of each model. If the likelihood distribution of any of the alternative models overlap with the likelihood ditribution of the best model, it means that there is no significant differences between the fit of both models.

![Alt text](images/alternative_models.png?raw=true "Title")
