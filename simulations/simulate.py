from os import system
from getopt import getopt
import sys

#directories
segdupDir = "/home/yaoban/segdup/"
kowhaiDir = "/home/yaoban/software/kowhai/"
multrecDir = "/home/yaoban/software/MultRec/Multrec/"

#kowhai options
nH = 15
nP = 5
rB = 1.0
pC = 0.5
pJ = 0.5

#segdup/multrec options
d = 10
l = 1

#replicates
replicates = 2

opts, args = getopt(sys.argv[1:], "d:l:r:",["segdup-dir=","kowhai-dir=","multrec-dir=","nH=","nP=","rB=","pC=","pJ=","replicates="])

for opt, arg in opts:
    if opt == "-d":
        d = int(arg)
    elif opt == "-l":
        l = int(arg)
    elif opt == "--segdup-dir":
        segdupDir = arg
        if segdupDir[-1] != "\/":
            segdupDir = segdupDir + "\/"
    elif opt == "--kowhai-dir":
        kowhaiDir = arg
        if kowhaiDir[-1] != "\/":
            kowhaiDir = segdupDir + "\/"
    elif opt == "--multrec-dir":
        multrecDir = arg
        if multrecDir[-1] != "\/":
            multrecDir = multrecDir + "\/"
    elif opt == "--nH":
        nH = int(arg)
    elif opt == "--nP":
        nP = int(arg)
    elif opt == "--rB":
        rB = float(arg)
    elif opt == "--pC":
        pC = float(arg)
    elif opt == "--pJ":
        pJ = float(arg)
    elif opt in ("-r", "--replicates"):
        replicates = int(arg)


system("rm summary.csv")

for _ in range(replicates):
    #run programs
    system(kowhaiDir + "kowhai --sim -nH " + nH + " -nP " + nP + " -nR 1 -rB " + rB + " -pC " + pC + " -pJ " + pJ + " --for-segdup --for-multrec --verbose")
    system("cat ./for-segdup-from-kowhai.txt | xargs " + segdupDir + "segdup -n 100000 -Tinit 10 -Tfinal 0.0 -d " + d + " -l " + l)

    multrecFile = open("for-multrec-from-kowhai.txt")
    multrecInput = multrecFile.readline()
    multrecFile.close()

    system(multrecDir + "Multrec -d " + d + " -l " + l + " " + multrecInput[:-3] + "\" -o multrec-output.txt")

    #parse multrec output

f = open("summary.csv")
output = open("results.csv", "w")

output.write("nCospec,nIndividualDups,nAllDupEvents,nJointDups,nXtinc,nHostSwitch,nLineageSort,codivs,dups,losses,cost,mrDups,mrLosses,mrCost\n")

for line in f:
    kowhaiOutput = next(f).rstrip()
    next(f)
    segdupOutput = next(f).rstrip()
    output.write(kowhaiOutput + segdupOutput + "\n")

f.close()
output.close()

#system("rm summary.csv")
