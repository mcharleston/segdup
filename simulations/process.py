from os import listdir

output = open("allresults.csv", "w")
output.write("nCospec,nIndividualDups,nAllDupEvents,nJointDups,nXtinc,nHostSwitch,nLineageSort,codivs,dups,losses,cost,sdTime,mrDups,mrLosses,mrCost,mrTime,nH,nP,rB,pJ\n")

for filename in listdir("."):
    if filename[:8] == "results-" and filename[-4:] == ".csv":
        params = filename[:-4].split("-")
        nH = int(params[1][2:])
        nP = int(params[2][2:])
        rB = float(params[3][2:])
        pJ = float(params[4][2:])

        #print(f'{nH} {nP} {rB} {pJ}')

        f = open(filename)
        next(f)

        for line in f:
            output.write(f'{line.rstrip()},{nH},{nP},{rB},{pJ}\n')

