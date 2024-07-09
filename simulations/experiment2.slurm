#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --time=0-48:00:00

#explicitly put all parameters for repeatability

# vary nH and nP

#python simulate.py --rB 1.0...2.0 -i 100000 --replicates 100
#python simulate.py --rB 2.5 -i 100000
#...
#python simulate.py --rB 5 ..

steps=10000
replicates=100

# EXPERIMENT 2

nH=50
nP=20
#for nH=20; do
#	for nP in {20}; do
		for rB in {1.0,2.0,3.0,4.0,5.0}; do
#			for pJ in {0.2,0.5,0.8,1}; do
pJ=0.5
				mkdir results-nH$nH-nP$nP-rB$rB-pJ$pJ
				cd results-nH$nH-nP$nP-rB$rB-pJ$pJ
				python ../simulate.py --nH $nH --nP $nP --rB $rB --pJ $pJ --replicates $replicates -i $steps &
				cd ..
#			done
		done
#	done
#done

wait

nH=50
nP=20
#for nH=20; do
#	for nP in {20}; do
		for rB in {1.0,2.0,3.0,4.0,5.0}; do
#			for pJ in {0.2,0.5,0.8,1}; do
pJ=0.5
				cd results-nH$nH-nP$nP-rB$rB-pJ$pJ
				mv results.csv ../results-nH$nH-nP$nP-rB$rB-pJ$pJ.csv
				cd ..
				rm -rf results-nH$nH-nP$nP-rB$rB-pJ$pJ
#			done
		done
#	done
#done


# NEXT expts:
# 2. set pJ=0.5 (or other equally interesting value) and iterate rB=1,2,3,4,5
# 3. set rB=2,pJ=0.5, iterate over nP 10...50
# 4. pJ=0.5, rB=2, nP=20, iterate nH over 20,50,100 or similar

