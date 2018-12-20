# Past Demography

## Infering the demographic history of *Senecio lautus* populations from the site frequency spectrum

This repository contains step-by-step instructions and code for infering demographic parameters of population pairs and triads of *Senecio lautus* using ```fastsimcoal```.

```fastsimcoal``` is... it requires three input files...


### Getting the site frequency spectrum (SFS) file

The SFS is a summary of genome-wide data describing the distribution of allele frequencies in a sample of one or more populations. For instance, it tells how many SNPs are derived in only one chromosome, two chromosomes, three chromosomes, and so on. Its shape is influenced by the history of the population: migration, population size changes, substructure, etc. The SFS treats all SNPs in the data set as independent of one another. Sample size is counted in haploid numbers.

A SFS is refered as *folded* when the information about the ancestral/derived state of the SNP is unavailable. Instead, the minor allele frequency is used as criteria for assigning an ancentral/derived-like state. A SFS is refered as *joint* when it summarises information from two or more populations. 

Since the SFS strictly only consider SNPs without missing data, most of the SNPs can be lost when all the individuals in a VCF file are used to generate the SFS. To avoid this, the data should be downsampled to the number of chromosomes that maximises the number of SNPs without missing data. This heavily relies on the quality of the sequence (coverage). In this regard, sequencing depth outweights sample size in importance.

[```easySFS```](https://github.com/isaacovercast/easySFS) is a Python script that generates a SFS file from a VCF file and a tab-delimited population specification file. The later file contains the sample names in the first column and the corresponding population names in the second column. This file could be directly generated from the VCF file using the custum Perl script ```getpopmap.pl```. It works for files containing either two or three populations.

```
perl getpopmap.pl input.vcf popmap.txt NamePop1 NamePop2 NamePop3
```

Before running the actual program, the number of chromosomes or haploid samples (herein projections) should be picked it up. For this, the program should be run in preview mode: 

```
easySFS.py -i input.vcf -p popmap.txt --preview -a
```
The prompted output looks like this:

```
Pop1
(2, 45.0)   (3, 59.0)   (4, 58.0)   (5, 59.0)   (6, 61.0)   (7, 65.0)   (8, 57.0)   (9, 50.0)   (10, 43.0)

Pop2
(2, 68.0)   (3, 96.0)   (4, 106.0)  (5, 110.0)  (6, 118.0)  (7, 109.0)   (8, 106.0)   (9, 96.0)   (10, 86.0)
```

In this example, projections ```7``` and ```6``` maximise the number of kept SNPs in population 1 and 2, respectively. These numbers should be specified with the ```--proj``` flag in the next step. ```-a``` flag means all the SNPs are considered, otherwise, a single SNP would be randomly sampled per RAD locus. ```-o``` flag specifies the output directory.

```
easySFS.py -i input.vcf -p popmap.txt --proj 7,6 -o output_folder -a
```

```easySFS.py``` generates several SFS files by default in two directories contained in the main output directory. Since ```fastsimcoal``` is picky with the format and naming of the input files, it should only be used the SFS files contained in the fastsimcoal directory. In this case, the ```input_jointMAFpop0_1.obs``` file. Beware this file name should be slightly modified as follow to be read by ```fastsimcoal```: ```input_jointMAFpop1_0.obs```.

### Getting the template file





### Getting the parameter file



