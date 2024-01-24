from os import system

segdupDir = "/home/yaoban/segdup/"
kowhaiDir = "/home/yaoban/software/kowhai/"
multrecDir = "/home/yaoban/software/MultRec/Multrec/"

system("rm summary.csv")

for _ in range(2):
    system(kowhaiDir + "kowhai --sim -nH 15 -nP 5 -nR 1 -rB 1.0 -pC 0.5 -pJ 0.5 -rX 1.0 --for-segdup --for-multrec --verbose")
    system("cat ./for-segdup-from-kowhai.txt | xargs " + segdupDir + "segdup -n 100000 -Tinit 10 -Tfinal 0.0 -d 10 -l 1")

f = open("summary.csv")
output = open("results.csv", "w")

output.write("nCospec,nIndividualDups,nAllDupEvents,nJointDups,nXtinc,nHostSwitch,nLineageSort,codivs,dups,losses,cost\n")

for line in f:
    kowhaiOutput = next(f).rstrip()
    next(f)
    segdupOutput = next(f).rstrip()
    output.write(kowhaiOutput + segdupOutput + "\n")

f.close()
output.close()

#system("rm summary.csv")
