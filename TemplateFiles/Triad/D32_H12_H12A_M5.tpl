// Three populations: D32_H12_H12A - M5 Bidirectional D32-H12A AND H12-H12A migration
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
0 0 MIG31
0 0 MIG32
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