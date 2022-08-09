#!/usr/bin/python

## look at the frequency depending on ancestry for the mutations present in the population sample

import sys
import csv

# script = sys.argv[0]

pathwd = sys.argv[1]
rec = sys.argv[2]
rep = sys.argv[3]

seqG = list(range(8000, 15000, 100))

for G in seqG:

    mutation = []
    genome = []
    nameFile = pathwd + "/genome2_" + str(G) + ".txt"
    with open(nameFile, 'r') as f:
        for l in f:
            line = l.rstrip('\n').split(" ")
            if (len(line) == 9) & (line[0][0]!='p'):
                subpop = line[6]
                if (subpop == 'p1') & (int(line[7]) > 5000):
                    code = '1'
                if subpop == 'p2':
                    code = '2'
                if (subpop == 'p1') & (int(line[7]) <= 5000):
                    code = '3'
                line.append(code)
                mutation.append(line)
            if line[0][0] == 'p':
    	        genome.append(line[2:]) 

    mutation.sort(key = lambda x: int(x[3]))
    mutation = [mutation[i] + [i] for i in range(0,len(mutation))]
    ancestry = []
    ## extract some information needed thereafter
    id = [i[0] for i in mutation]
    pos = [i[3] for i in mutation]
    ix0 = [i[10] for i in mutation]
    ## extract the index and position of the mutations that are on the same site
    doublePos = [i for i in set(pos) if pos.count(i) > 1]
    doubleIx = [i for i in ix0 if doublePos.count(pos[i]) > 0]
    ## put together the index with the same position
    posIx = [doublePos.index(mutation[i][3]) for i in doubleIx]
    doubleIx1 = [[doubleIx[i-1],doubleIx[i]] for i in range(1,len(doubleIx)) if posIx[i] == posIx[i-1]]
    while len(doubleIx1) != len(doublePos):
        irm = []
        for i in range(1,len(doubleIx1)):
            if doubleIx1[i][0] == doubleIx1[i-1][len(doubleIx1[i-1])-1]:
                doubleIx1[i-1].append(doubleIx1[i][1])
                irm.append(doubleIx1[i])
            else : 
                doubleIx1[i-1].append(-1)
        irm = [doubleIx1.remove(i) for i in irm]
        if len(doubleIx1[len(doubleIx1)-1]) < len(doubleIx1[len(doubleIx1)-2]):
            doubleIx1[len(doubleIx1)-1].append(-1)

    doubleIx1 = [[i for i in r] for r in zip(*doubleIx1)]

    for i in genome:
        ## extract the 'code' for the mutation with a known origine
        ibis = [j for j in i if mutation[id.index(j)][9] != '3']
        ## get the index of these mutations
        ix = [mutation[id.index(j)][10] for j in ibis]
        ## get the index of the mutations that are not in the genome (or with code = 3)
        missingIx = [j for j in range(0,len(id)) if ix.count(j) == 0]
        ## extract the index of the mutation (code != 3) that are on a site with multiple mutations
        iix = [doubleIx1[k].index(j) for j in ix if doubleIx.count(j) > 0 for k in range(0,len(doubleIx1)) if doubleIx1[k].count(j) > 0]
	    ## same but for the mutations not in the genome (or thise with code = 3)
        missingIix = [doubleIx1[k].index(j) for j in missingIx if doubleIx.count(j) > 0 for k in range(0,len(doubleIx1)) if doubleIx1[k].count(j) > 0]
        ## remove the doublon (with set) 
        missingIix = list(set([j for j in missingIix if iix.count(j) == 0]))
        ## keep only the mutation that are on a site with only one mutation
        ## and concatenate it with the mutations that are on a site with multiple mutations (after removing those that are in double
        missingIx = [j for j in missingIx if doubleIx.count(j) == 0] + missingIix
        missingIx.sort()
        cn = []
        for m in missingIx:
            ixp = [l for l in ix if l-m > 0]
            if len(ixp) > 0: 
                dp = int(pos[min(ixp)]) - int(pos[m])
            else:
                dp = 1600000
            ixm = [l for l in ix if m-l > 0]
            if len(ixm) > 0: 
                dm = int(pos[m]) - int(pos[max(ixm)]) 
            else:
                dm = 1600000
            if dp > dm:
                cn.append(max(ixm))
            else:
                cn.append(min(ixp))
        ix = cn + ix
        ix.sort()
        ancestry.append([mutation[ix0.index(j)][9] for j in ix])

    if (len(set([len(i) for i in ancestry])) == 1) & (len(ancestry[0]) == len(set(pos))):
        ancestry = [[i for i in r] for r in zip(*ancestry)]
        ## only need to know how many genomes have ancestry = 1
        ## ancestry = 2 will just be the total number of genomes minus the number of ancestry 1
        ancestry = [i.count('1') for i in ancestry]

    with open("mutationAncestry2Proportion_" + pathwd + ".txt", "a") as f:
        f.write(str(G) + " ")
        [f.write(str(i) + " ") for i in ancestry]
        f.write("\n")

