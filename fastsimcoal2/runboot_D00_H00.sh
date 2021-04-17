#Run boostrap for D00-H00
#PBS -A UQ-SCI-BiolSci
#PBS -l select=1:ncpus=9:mem=1GB
#PBS -l walltime=1:00:00
#PBS -J 1-100

cd [PATH]/D00_H00/D00_H00_${PBS_ARRAY_INDEX}/; [PATH]/fsc26 -t D00_H00.tpl -e D00_H00.est --initValues D00_H00.pv -n100000 -m -L30 -M -c10
