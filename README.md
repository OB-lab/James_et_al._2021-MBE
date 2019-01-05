# Past Demography

## Infering the demographic history of *Senecio lautus* populations

This repository contains step-by-step instructions and code for infering demographic parameters of population pairs and triads of *Senecio lautus* using ```fastsimcoal```.

```fastsimcoal``` is a continuous-time coalescent simulator of genomic diversity under arbitrarily complex evolutionary scenarios. It can estimate demographic parameters from the site frequency spectrum through a composite likelihood maximisation procedure. For working, it requires three input files: **1)** a site frequency spectrum file, **2)** a template file, and **3)** an estimation file.


### Getting the site frequency spectrum (SFS) file

The SFS is a summary of genome-wide data describing the distribution of allele frequencies in a sample of one or more populations. For instance, it tells how many SNPs are derived in only one chromosome, two chromosomes, three chromosomes, and so on. Its shape is influenced by the history of the population: migration, population size changes, substructure, etc. The SFS treats all SNPs in the data set as independent of one another. Sample size is counted in haploid numbers.

A SFS is refered as *folded* when the information about the ancestral/derived state of the SNP is unavailable. Instead, the minor allele frequency is used as criteria for assigning an ancentral/derived-like state. A SFS is refered as *joint* when it summarises information from two or more populations. 

Since the SFS strictly only consider SNPs without missing data, most of the SNPs can be lost when all the individuals in a VCF file are used to generate the SFS. To avoid this, the data should be downsampled to the number of chromosomes that maximises the number of SNPs without missing data. This heavily relies on the quality of the sequence (coverage). In this regard, sequencing depth outweights sample size in importance.

[```easySFS```](https://github.com/isaacovercast/easySFS) is a Python script that generates a SFS file from a VCF file and a tab-delimited population specification file. The later file contains the sample names in the first column and the corresponding population names in the second column. This file could be directly generated from the VCF file using the custum Perl script ```getpopmap.pl```. It works for files containing either two or three populations.

```
perl getpopmap.pl input.vcf popmap.txt NamePop1 NamePop2 NamePop3
```

Before running ```easySFS```, the number of chromosomes or haploid samples (herein projections) should be picked it up. For this, the program should be run in preview mode: 

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

```easySFS``` generates several SFS files by default in two directories contained in the main output directory. Since ```fastsimcoal``` is picky with the format and naming of the input files, it should only be used the SFS files contained in the fastsimcoal directory. Beware all files should be slightly renamed to be read by ```fastsimcoal```. The two numbers at the end of the name should be swapped. Hence, ```input_jointMAFpop0_1.obs``` should be renamed as ```input_jointMAFpop1_0.obs``` .

### Getting the template file

The template file specifies the evolutionary model and the parameters that should be estimated. It is divided in ten fixed sections in the following order:

  + Number of populations.
  
  + Populations effective size: It should contain one line per population, as mentioned in the previous section. If the effective population size is unknown, a parameter name should be specified instead to be calculated (*i.e.* NPOP1, NPOP2). Sizes correspond to the number of genes present in a population, it means twice the number of individuals for a diploid species. The order of the populations all along this file should be the same of the SFS file, which usually is alphabetically sorted. The program indexes the first population as 0, the second one as 1, and so on.
  
  + Samples size: It should contain one line per population, as mentioned in the first section, and follow the populations order specified in the previous section. Sizes should be given in haploid numbers. It corresponds to the number of projections set up in ```easySFS``` when generating the SFS file.
  
  + Growth rates: It should contain one line per population. A zero growth rate means stationary population size. Negative growth implies population expansion since it is estimated backward in time.
  
  + Number of migration matrixes.
  
  + Migration matrix: It should contain one matrix per migration matrix, as mentioned in the previous section. The program indexes the first migration matrix as 0, the second one as 1, and so on. Columns/rows keep the logic from/to. If the migration rate is unknown, a parameter name should be specified instead in the corresponding cell to be calculated (*i.e.* MIG12, MIG21).
  
  + Historical events. This section is the heart of the evolutionary model. The first line specifies the number of historical events, namely gene flow, bottle neck, admixture event, coalescence of population, etc. Every historical event should be specified in a different line. Each historical event is defined by 7 numbers in the following order. If any of those values are unknown, a parameter name should be specified instead to be calculated (*i.e.* TDIV, RESIZE).
    
    - Number of generations in the past at which the historical event happened.
    
    - Source population.
    
    - Sink population.
    
    - Proportion of migrants moving from source to sink population. It should be 1 if both populations coalesce at this event.
    
    - New size for the sink population backward in time.
    
    - New growth rate for the sink population backward in time.
    
    - New migration matrix index to be used further back in time.
  
  + Number of independent loci: The number of independent chromosomes to be simulated and a flag indicating if the different chromosomes have a different (1) or a similar (0) structure. 
  
  + Number of linkage blocks: The number of blocks per chromosome that may differ by the type of markers to be simulated, the recombination rate, or the mutation rate.
  
  + Genetic properties: Data type, number of markers, recombination rate, and mutation rate. FREQ in the data type field specifies to estimate the SFS with the simulations.
  
Every section is introduced by a comment line starting with the characters ```//```. This is how a template file for a population pair looks like:

```
// Two populations - D00_H00 - M5 Bidirectional secondary contact
2
//Population effective sizes: to be estimated
NPOP1
NPOP2
//Samples haploid size: the same one defined in the projections. 
48
54
//Growth rates: negative growth implies population expansion
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//Migration matrix 0: to be estimated
0 MIG21
MIG12 0
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
And this is how a template file for a population triad looks like:

```
// Three populations: D32_H12_H12A - M8 Bidirectional migration all combinations
3
//Population effective sizes: to be estimated.
D32
H12
HA12
//Samples haploid size: the same one defined in the projections. 
64
62
64
//Growth rates: negative growth implies population expansion
0
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//Migration matrix 0
0 MIG21 MIG31
MIG12 0 MIG32
MIG13 MIG23 0
//Migration matrix 1
0 0 0
0 0 0
0 0 0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index
2 historical event
TDIV1 2 0 1 RESIZE1 0 1
TDIV2 1 0 1 RESIZE2 0 1
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per gen recomb and mut rates
FREQ 1 0 1e-8
```


### Getting the parameter file




### Running fsc




### Summarising the results




### Selecting the best demographic model




### Estimating the confidence intervals




### References
fastsimcoal web site...
