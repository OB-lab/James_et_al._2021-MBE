#Run boostrap for POP1_POP2
#PBS -A UQ-SCI-BiolSci
#PBS -l select=1:ncpus=9:mem=1GB
#PBS -l walltime=1:00:00
#PBS -J 1-100

cd /PATH/POP1_POP2/POP1_POP2_${PBS_ARRAY_INDEX}/; /PATH/fsc26_linux64/fsc26 -t POP1_POP2.tpl -e POP1_POP2.est --initValues POP1_POP2.pv -n100000 -m -L30 -M -c10

