#!/usr/bin/bash
for i in {1..100} ; do
	../kowhai/Release/kowhai --sim -nH 15 -nP 5 -nR 1 -rB 1.0 -pC 0.5 -pJ 0.5 -rX 1.0 --verbose
	cat ./for-segdup-from-kowhai.txt | xargs ./Release/segdup -n 100000 -Tinit 10 -Tfinal 0.0 -o samples -o trace -d 10 -l 1
done
