#Run fsc26 - D00_H00
#PBS -A UQ-SCI-BiolSci
#PBS -l select=1:ncpus=9:mem=1GB
#PBS -l walltime=5:00:00
#PBS -J 1-50

cd [PATH]/D00_H00/1/run${PBS_ARRAY_INDEX}; [PATH]/fsc26 -t D00_H00.tpl -n100000 -m -e D00_H00.est -M -L40 -c10 -q

cd [PATH]/D00_H00/2/run${PBS_ARRAY_INDEX}; [PATH]/fsc26 -t D00_H00.tpl -n100000 -m -e D00_H00.est -M -L40 -c10 -q

cd [PATH]/D00_H00/3/run${PBS_ARRAY_INDEX}; [PATH]/fsc26 -t D00_H00.tpl -n100000 -m -e D00_H00.est -M -L40 -c10 -q

cd [PATH]/D00_H00/4/run${PBS_ARRAY_INDEX}; [PATH]/fsc26 -t D00_H00.tpl -n100000 -m -e D00_H00.est -M -L40 -c10 -q

cd [PATH]/D00_H00/5/run${PBS_ARRAY_INDEX}; [PATH]/fsc26 -t D00_H00.tpl -n100000 -m -e D00_H00.est -M -L40 -c10 -q

cd [PATH]/D00_H00/6/run${PBS_ARRAY_INDEX}; [PATH]/fsc26 -t D00_H00.tpl -n100000 -m -e D00_H00.est -M -L40 -c10 -q

cd [PATH]/D00_H00/7/run${PBS_ARRAY_INDEX}; [PATH]/fsc26 -t D00_H00.tpl -n100000 -m -e D00_H00.est -M -L40 -c10 -q

cd [PATH]/D00_H00/8/run${PBS_ARRAY_INDEX}; [PATH]/fsc26 -t D00_H00.tpl -n100000 -m -e D00_H00.est -M -L40 -c10 -q

cd [PATH]/D00_H00/9/run${PBS_ARRAY_INDEX}; [PATH]/fsc26 -t D00_H00.tpl -n100000 -m -e D00_H00.est -M -L40 -c10 -q

cd [PATH]/D00_H00/10/run${PBS_ARRAY_INDEX}; [PATH]/fsc26 -t D00_H00.tpl -n100000 -m -e D00_H00.est -M -L40 -c10 -q
