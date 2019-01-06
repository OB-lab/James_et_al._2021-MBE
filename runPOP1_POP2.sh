#Run fsc26 - POP1_POP2
#PBS -A UQ-SCI-BiolSci
#PBS -l select=1:ncpus=9:mem=1GB
#PBS -l walltime=4:00:00
#PBS -J 1-75

cd /PATH/POP1_POP2/1/run${PBS_ARRAY_INDEX}; /PATH/fsc26 -t POP1_POP2.tpl -n100000 -m -e POP1_POP2.est -M -L60 -c10 -q

cd /PATH/POP1_POP2/2/run${PBS_ARRAY_INDEX}; /PATH/fsc26 -t POP1_POP2.tpl -n100000 -m -e POP1_POP2.est -M -L60 -c10 -q

cd /PATH/POP1_POP2/3/run${PBS_ARRAY_INDEX}; /PATH/fsc26 -t POP1_POP2.tpl -n100000 -m -e POP1_POP2.est -M -L60 -c10 -q

cd /PATH/POP1_POP2/4/run${PBS_ARRAY_INDEX}; /PATH/fsc26 -t POP1_POP2.tpl -n100000 -m -e POP1_POP2.est -M -L60 -c10 -q

cd /PATH/POP1_POP2/5/run${PBS_ARRAY_INDEX}; /PATH/fsc26 -t POP1_POP2.tpl -n100000 -m -e POP1_POP2.est -M -L60 -c10 -q

cd /PATH/POP1_POP2/6/run${PBS_ARRAY_INDEX}; /PATH/fsc26 -t POP1_POP2.tpl -n100000 -m -e POP1_POP2.est -M -L60 -c10 -q

cd /PATH/POP1_POP2/7/run${PBS_ARRAY_INDEX}; /PATH/fsc26 -t POP1_POP2.tpl -n100000 -m -e POP1_POP2.est -M -L60 -c10 -q
