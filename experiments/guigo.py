import os

output = open("guigo_temps.csv", "w")
output.write("t,s7,s10\n")

for i in range(1,51):
    t = i/5

    os.system(f"cat guigo_all.txt | xargs ../segdup -n 100000 -o trace -Tinit 20 -d 50 -l 1 -o interval 10 -Tfinal {t} -nfinal 100000 -o final -o samples")

    f = open("segdup-samples.csv")

    counts = [0,0]

    for l in f:
        if "[=]s7" in l:
            counts[0] += int(l.split()[1])

        if "[=]s10" in l:
            counts[1] += int(l.split()[1])

    f.close()

    output.write(f"{t},{','.join([str(c) for c in counts])}\n")

