#!/usr/bin/python

## python version 3.10.0

genome = []
i = 0
with open("alignment_core_ref/26695_core_aligned.aln", 'rt') as f:
    for l in f:
        if(i % 2 == 1):
            genome.append(l)
        i += 1

genome = genome[0]
posCore = [i + 1 for i in range(0, len(genome)) if(genome[i] != '-') & (genome[i] != '.')]

coreCDS = []
ref = []
i = 0
with open("../DATA/26695.gff", 'rt') as f:
    for l in f:
        seq = l.rstrip("\n").split("\t")
        if(len(l) == 1667868):
            ref.append(l)
        elif((i > 1) & (len(seq) > 2)):
       	    posCDS = list(range(int(seq[3]), int(seq[4])+1))
            inCore = sum([0 if j in posCore else 1 for j in posCDS])
            if((seq[2] == 'CDS') & (inCore == 0)):
                coreCDS.append(seq)
        print(i)
        i +=1

with open("core_CDS.gff", "a") as f:
    for i in coreCDS:
        [f.write(j + " ") for j in i]
        f.write("\n")

