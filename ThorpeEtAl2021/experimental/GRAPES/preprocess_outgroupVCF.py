#!/usr/bin/python

## python version 3.10.0

genome = []
i = 0
with open("alignment_core_ref/26695_core_aligned.aln", 'rt') as f:
    for l in f:
        if(i % 2 == 1):
            genome.append(l)
        i += 1

outgroup = []
i = 0
with open("outgroupAlignment/Hacinonychis_aligned.aln", 'rt') as f:
    for l in f:
        if(i % 2 == 1):
            outgroup.append(l)
        i += 1

diffSeq = []
for i in range(0,len(genome[0])):
    if((genome[0][i] != outgroup[0][i]) & (genome[0][i] != "-") & (genome[0][i] != ".") & (outgroup[0][i] != "-") & (outgroup[0][i] != ".")):
        diffSeq.append([i+1, genome[0][i], outgroup[0][i]])

with open("outgroup_core_diff.txt", "a") as f:
    f.write("position " + "ref " + "outgroup" + "\n")
    for i in diffSeq:
        [f.write(str(j) + " ") for j in i]
        f.write("\n"