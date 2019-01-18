// Three populations: D32_H12_H12A - M8 Bidirectional migration all combinations
3
//Population effective sizes: to be estimated
D32
H12
HA12
//Samples haploid size: the same one defined in the projections
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
0 MIG12 MIG13
MIG21 0 MIG23
MIG31 MIG32 0
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
