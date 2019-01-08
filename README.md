# Past Demography

This repository contains step-by-step instructions and code for infering demographic parameters of population pairs and triads of *Senecio lautus* using ```fastsimcoal``` and ```TreeMix```.

## Infering the demographic history of *Senecio lautus* populations using fastsimcoal

```fastsimcoal``` is a continuous-time coalescent simulator of genomic diversity under arbitrarily complex evolutionary scenarios. It can estimate demographic parameters from the site frequency spectrum through a composite likelihood maximisation procedure. Since this approach requires the a priori formulation of the demographic models to test, we restrict its use to population pairs and triads of interest given their phylogenetic relationships and occurrence patterns.

For working, ```fastsimcoal``` requires three input files: **1)** a site frequency spectrum file, **2)** a template file, and **3)** an estimation file.

### Getting the site frequency spectrum (SFS) file

The SFS is a summary of genome-wide data describing the distribution of allele frequencies in a sample of one or more populations. For instance, it tells how many SNPs are derived in only one chromosome, two chromosomes, three chromosomes, and so on. Its shape is influenced by the history of the population: migration, population size changes, substructure, etc. The SFS treats all SNPs in the data set as independent of one another. Sample size is counted in haploid numbers.

A SFS is refered as *folded* when the information about the ancestral/derived state of the SNP is unavailable. Instead, the minor allele frequency is used as criteria for assigning an ancentral/derived-like state. A SFS is refered as *joint* when it summarises information from two or more populations. 

Since the SFS strictly only consider SNPs without missing data, most of the SNPs can be lost when all the individuals in a VCF file are used to generate the SFS. To avoid this, the data should be downsampled to the number of chromosomes that maximises the number of SNPs without missing data. This heavily relies on the quality of the sequence. In this regard, sequencing depth outweights sample size in importance.

```easySFS``` is a Python script that generates a SFS file from a VCF file and a tab-delimited population specification file. The later file contains the sample names in the first column and the corresponding population names in the second column. This file could be directly generated from the VCF file using the custum Perl script ```getpopmap.pl```. It works for files containing either two or three populations.

```
perl getpopmap.pl input.vcf popmap.txt NamePop1 NamePop2 NamePop3
```

The first argument is the name of the input VCF file, the second argument is the name of the output popmap file, and the following two or three arguments are the names of the populations. Beware the names of the samples in the VCF file should stat with the population name followed by a hyphen.


Before running ```easySFS```, the number of chromosomes or haploid samples, herein projections, should be picked it up. For this, the program should be run in preview mode: 

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

In this example, projections ```7``` and ```6``` maximise the number of kept SNPs in population 1 and 2, respectively. These numbers should be specified with the ```--proj``` flag in the next step. The ```-a``` flag means all the SNPs are considered, otherwise, a single SNP would be randomly sampled per RAD locus. The ```-o``` flag specifies the output directory.

```
easySFS.py -i input.vcf -p popmap.txt --proj 7,6 -o output_directory -a
```

```easySFS``` generates several SFS files by default in two directories contained in the main output directory. Since ```fastsimcoal``` is picky with the format and naming pattern of the input files, it should only be used the SFS files contained in the fastsimcoal directory. Beware all files should be slightly renamed to be read by ```fastsimcoal```. The two numbers at the end of the name should be swapped. Hence, ```input_jointMAFpop0_1.obs``` should be renamed as ```input_jointMAFpop1_0.obs```.

### Getting the template file

The template file specifies the evolutionary model and the parameters that should be estimated. It is divided in ten fixed sections in the following order:

  + Number of populations.
  
  + Populations effective size: It should contain one line per population, as mentioned in the previous section. If the effective population size is unknown, a parameter name should be specified instead to be calculated (*i.e.* NPOP1, NPOP2). Sizes correspond to the number of genes present in a population, it means twice the number of individuals for a diploid species. The order of the populations all along this file should be the same of the populations in the SFS file, which usually is alphabetically sorted based on the population names. The program indexes the first population as 0, the second one as 1, and so on.
  
  + Samples size: It should contain one line per population, as mentioned in the first section, and follow the populations order specified in the previous section. Sizes should be given in haploid numbers. It corresponds to the number of projections set up in ```easySFS``` when generating the SFS file.
  
  + Growth rates: It should contain one line per population. A zero growth rate means stationary population size. Negative growth implies population expansion since it is estimated backward in time.
  
  + Number of migration matrixes.
  
  + Migration matrix: It should contain one matrix per migration matrix, as mentioned in the previous section. The program indexes the first migration matrix as 0, the second one as 1, and so on. Columns/rows keep the logic from/to. If the migration rate is unknown, a parameter name should be specified instead in the corresponding cell to be calculated (*i.e.* MIG12, MIG21).
  
  + Historical events. This section is the heart of the evolutionary model. The first line specifies the number of historical events, namely gene flow, bottle neck, admixture event, coalescence of population, etc. Each historical event should be specified in a different line. Each historical event is defined by 7 numbers in the following order. If any of those values are unknown, a parameter name should be specified instead to be calculated (*i.e.* TDIV, RESIZE).
    
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
  
Every section starts with a comment line starting with the characters ```//```. The name of the template file must be the root name of the SFS file plus the extention ```.tpl```. For instance, the corresponding template file name of ```D00_H00_jointMAFpop1_0.obs``` should be ```D00_H00.tpl```.

For the Dune-Headland population pairs, we considered seven demographic models ranging from no migration, free migration, and migration after secondary contact. Sample template files for the D00-H00 pair are in the directory ```TemplateFiles/Pair```. Beware file names were modified to distinguish among models here. They assume Dune population comes first in the SFS file.

  + Model 1: No migration.
  + Model 2: Bidirectional migration.
  + Model 3: Dune to Headland migration.
  + Model 4: Headland to Dune migration.
  + Model 5: Bidirectional migration after secondary contact.
  + Model 6: Dune to Headland migration after secondary contact.
  + Model 7: Headland to Dune migration after secondary contact.

This is how a template file for a population pair looks like:

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

```NPOP1```, ```NPOP2```, ```MIG21```, ```MIG12```, ```TDIV```, ```RESIZE```, and ```TSEC``` are the unknown parameters to be estimated by ```fastsimcoal```. They all should be as well specified in the corresponding estimation file.

For the population triads, we considered a nested-models approach because the number of possible demographic models considerably increases with more than two populations. Nine initial models, ranging from no migration to bidirectional migration among all populations, are first fitted to data and then different nested variants of the best model are explored, such as unidirectional migration or migration after secondary contact.

  + Model 1: No migration. 
  + Model 2: Bidirectional migration between A and B.
  + Model 3: Bidirectional migration between A and C.
  + Model 4: Bidirectional migration between B and C.
  + Model 5: Bidirectional migration between A and B AND A and C.
  + Model 6: Bidirectional migration between A and B AND B and C.
  + Model 7: Bidirectional migration between A and C AND B and C.
  + Model 8: Bidirectional migration among all populations.
  + Model 9: Bidirectional migration between the A-B ancestor and C.

Beware that, contrary to the population pairs case where the template files of the same model are quite the same among pairs because **1)** the Dune population always comes first in the SFS file and the Headland population comes second and **2)** the phylogenetic relationships among populations are the same regardless their orden, the template files of the same model vary from one triad to another and they should be carefully inspected case by case. Sample template files for the D32-H12-H12A triad are in the directory ```TemplateFiles/Triad```. Beware file names were modified to distinguish among models here.

This is how a template file for a population triad looks like:

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

### Getting the estimation file

The estimation file specifies the search range of the parameters defined in the template file. The lower range limit is an absolute minimum, whereas the upper range is only used as a maximum for choosing a random initial value for the parameter. It is divided in three fixed sections in the following order:
  
  + Parameters: It lists the main parameters and its respective initial search range. Each parameter can either be uniformly or log-uniformly distributed and its values can be either recorded or omitted in the output files. Some parameters, such as ```ANCSIZE```, are specified in this section despite they were not invoked in the template file. This is because they are later used to calculate other complex parameters that indeed were specified in the template file, such us ```RESIZE``` (see below).
  
  + Rules: It includes a set of conditions to be met among the above simple parameters. For instance, in a secondary contact model, it should be specified that the secondary contact event happened after the divergence of the involved populations. 
  
  + Complex parameters: It defines parameters that are obtained as simple operations between other parameters. For instance, the migration rate ```MIG12``` is calculated as the number of migrants ```NM12``` divided by the population size of the sink population ```NPOP2```.
  
The name of the template file must be the root name of the SFS file plus the extention ```.est```. For instance, the corresponding template file name of ```D00_H00_jointMAFpop1_0.obs``` should be ```D00_H00.est```. Sample estimation files for the D00_H00 pair are in the directory ```EstimationFiles/Pair``` and for the D32-H12-H12A triad are in the directory ```EstimationFiles/Triad```.

This is how a template file for a population pair looks like:

```
// Priors and rules file
// *********************

[PARAMETERS]
//#isInt? #name   #dist.#min  #max
//all N are in number of haploid individuals
1  ANCSIZE     unif     100  100000   output
1  NPOP1       unif     100  100000   output
1  NPOP2       unif     100  100000   output
0  NM12       logunif  1e-2 20       hide
0  NM21       logunif  1e-2 20       hide
1  TDIV        unif     100   20000   output
1  TSEC        unif     100   20000   output

[RULES]
TDIV > TSEC

[COMPLEX PARAMETERS]

0  RESIZE = ANCSIZE/NPOP2     hide
0  MIG12  = NM12/NPOP2       output
0  MIG21  = NM21/NPOP1       output
```

And this is how a template file for a population triad looks like:

```
// Priors and rules file
// *********************

[PARAMETERS]
//#isInt? #name   #dist.#min  #max
//all N are in number of haploid individuals
1  ANCSIZE     unif     100  100000   output
1  HSIZE     unif     100  100000   output
1  D32       unif     100  100000   output
1  H12       unif     100  100000   output
1  HA12       unif     100  100000   output
1  TDIV1        unif     100   20000   output
1  TDIV2        unif     100   20000   output
0  NM12       logunif  1e-2 20       hide
0  NM21       logunif  1e-2 20       hide
0  NM13       logunif  1e-2 20       hide
0  NM31       logunif  1e-2 20       hide
0  NM23       logunif  1e-2 20       hide
0  NM32       logunif  1e-2 20       hide

[RULES]

TDIV1 < TDIV2

[COMPLEX PARAMETERS]

0  RESIZE1 = HSIZE/D32     hide
0  RESIZE2 = ANCSIZE/HSIZE     hide
0  MIG12  = NM12/H12       output
0  MIG21  = NM21/D32       output
0  MIG13  = NM13/HA12       output
0  MIG31  = NM31/D32       output
0  MIG23  = NM23/HA12       output
0  MIG32  = NM32/H12       output
```

### Running fsc

Running ```fastsimcoal``` is easy if the three input files are well formatted. All the files should be in the same folder from where the command is invoked. 

```
fsc26 -t D00_H00.tpl -n100000 -m -e D00_H00.est -M -L60 -c10 -q
```

The ```-t``` and  ```-e``` flags specify the template and estimation files, respectively, the ```-n``` flag specifies the number of SFS simulations to fit to the data, the ```-M``` flag indicates to perform parameter estimation by maximum composite likelihood from the SFS, the ```-L``` flag specifies the number of optimisation cycles, the ```-c``` specifies the number of threads to be used for simulations, and the ```-q``` flag keeps output messages at the minimum level.

Reaching realible results depends at some extend on how well the parameters space is explored. For this, every model should be independently run between 50 and 100 times and the results of the run with the best maximum likelihood value should be picked up for further analysis. The custum Perl script ```mkdir_in.pl``` helps to create the required directories and copy the input files into them for launching the independent runs.

```
perl mkdir_in.pl 7 75
```

The first argument corresponds to the number of models to test and the second argument indicates how many independent runs will be launch per model. Running this is best done in parallel in a server. The custum shell executable script ```runPOP1_POP2.sh``` just launch the jobs this way.

```
qsub runPOP1_POP2.sh
```

Both ```mkdir_in.pl``` and ```runPOP1_POP2.sh``` should be invoked from the location that contains all the models in separate directories named as consecutive number from 1 to the maximum number of different models. Each directory should contain the three input files for running ```fastsimcoal```.

### Summarising the results

```fastsimcoal``` generates an output directory, named as the root name of the input files (*i.e.* POP1_POP2), that contains most of the relevant output files. ```POP1_POP2.bestlhoods``` is particularly useful for summarising the estimated parameter values and ```POP1_POP2.brent_lhoods``` for evaluating the run performance. The custum Perl script ```extract_ml.pl``` reads the ```POP1_POP2.bestlhoods``` across all the runs per model and summarises their maximum likelihood value in a table. It also copies both of the above mentioned files of the best run per model in a new directory, named ```results_POP1_POP2```. 

```
perl extract_ml.pl 7 75 D00_H00
```

The first argument corresponds to the number of models, the second argument corresponds to the number of independent runs per model, and the third argument corresponds to the root name of the files. ```extract_ml.pl``` should be invoked from the location that contains all the models in separate directories named as consecutive number from 1 to the maximum number of different models. 

The custum R script ```summaryplots_fsc.R``` graphically summarises the performance of fastsimcoal runs across the different models and populations. It generates a PDF file showing, in the first section, boxplots of the maximum likelihood values of all runs per model per population pair and triad and, in the second section, scatter plots of the likelihood values along the optimisation cycles of the best run per model per population pair and triad. By inspecting these plots, we can get a quick idea about the overall performance of the parameter space exploration.


### Selecting the best demographic model




### Estimating the confidence intervals


## Infering the demographic history of *Senecio lautus* populations using TreeMix



## References
fastsimcoal web site...
