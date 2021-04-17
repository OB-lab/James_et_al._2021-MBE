#fsc ML dist - D00-H00
#PBS -A UQ-SCI-BiolSci
#PBS -l select=1:ncpus=9:mem=1GB
#PBS -l walltime=0:30:00
#PBS -J 1-100

mkdir [PATH]/maxL_${PBS_ARRAY_INDEX}; cd [PATH]/maxL_${PBS_ARRAY_INDEX}; cp ../D00_H00_maxL* .; [PATH]/fsc26 -i D00_H00_maxL.par -n100000 -m -q -0; sed -n '2,3p' D00_H00_maxL/D00_H00_maxL.lhoods  >> ../D00_H00.maxL.lhoods; cd ..; rm -rf [PATH]/maxL_${PBS_ARRAY_INDEX}

mkdir [PATH]/A_${PBS_ARRAY_INDEX}; cd [PATH]/A_${PBS_ARRAY_INDEX}; cp ../D00_H00_A* .; [PATH]/fsc26 -i D00_H00_A.par -n100000 -m -q -0; sed -n '2,3p' D00_H00_A/D00_H00_A.lhoods  >> ../D00_H00.A.lhoods; cd ..; rm -rf [PATH]/A_${PBS_ARRAY_INDEX}

mkdir [PATH]/B_${PBS_ARRAY_INDEX}; cd [PATH]/B_${PBS_ARRAY_INDEX}; cp ../D00_H00_B* .; [PATH]/fsc26 -i D00_H00_B.par -n100000 -m -q -0; sed -n '2,3p' D00_H00_B/D00_H00_B.lhoods  >> ../D00_H00.B.lhoods; cd ..; rm -rf [PATH]/B_${PBS_ARRAY_INDEX}

mkdir [PATH]/C_${PBS_ARRAY_INDEX}; cd [PATH]/C_${PBS_ARRAY_INDEX}; cp ../D00_H00_C* .; [PATH]/fsc26 -i D00_H00_C.par -n100000 -m -q -0; sed -n '2,3p' D00_H00_C/D00_H00_C.lhoods  >> ../D00_H00.C.lhoods; cd ..; rm -rf [PATH]/C_${PBS_ARRAY_INDEX}
