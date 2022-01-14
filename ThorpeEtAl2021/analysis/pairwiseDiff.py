#!/usr/bin/python

# calculate the pairwise difference 
# of the strains within a population

import sys
import random
import csv

# script = sys.argv[0]

pathwd = sys.argv[1]
rec = sys.argv[2]
rep = sys.argv[3]

ns = 14
nsample = 100
typeMut = [ "m" + str(i+1) for i in range(ns) ]
npop = 2

seqG1 = list(range(100,5000,300)) + list(range(5000,5700,100)) + list(range(5700,15000,100))
seqG2 = list(range(5000,5700,100)) + list(range(5700,15000,100))
seqG3 = list(range(8000,15000,100))


generation = []
population = []
popSize = []



with open(pathwd + "/out.csv") as csvDataFile:
    raw = csv.reader(csvDataFile)
    for row in raw:
        generation.append(row[0])
        population.append(row[1])
        popSize.append(row[2])


for pop in range(npop):
    pop += 1
    
    if pop == 1:
        seqG = seqG1
    elif pop == 2:
        seqG = seqG2
    elif pop == 3:
        seqG = seqG3
        

    seqG.append(15000)
    
    for G in seqG:

        print("G: " + str(G) + " pop: " + str(pop))

        N = popSize[generation.index(str(G)) + pop -1]
        nameFile = pathwd + "/genome" + str(pop) + "_" + str(G) + ".txt"

        with open(nameFile, 'r') as f:
            # get No of columns in each line
            nbCol = [ len(l.split(" ")) for l in f.readlines() ]
   
        nbLine = len(nbCol)
        ix = []
        for i,j in enumerate(nbCol):
            if j == 1:
                ix.append(i)
        

        dm = []
        dg = []
        nl = 0
        with open(nameFile, 'rt') as f:
            for l in f:
                nl += 1
                if(nl > ix[0]) & (nl <= ix[1]):
                    dm.append(l.rstrip('\n').split())
                elif nl > ix[1]:
                    dg.append(l.rstrip('\n').split())
        


        dm = dm[1:]
        dg = dg[1:]
        
        
        mutID = [i[0] for i in dm]
        mutType = [i[2] for i in dm]
        
        # sample = random.sample(range(int(N)), nsample)
        totDiff = []

        for ind0 in range(nsample-1):
            for ind1 in range(ind0+1,nsample):
                x = dg[ind0][2:]
                y = dg[ind1][2:]
                diff = list(set(x) - set(y)) + list(set(y) - set(x))
                ix = [mutID.index(i) for i in diff]
                diffType = [mutType[i] for i in ix]
                totDiff.append([diffType.count(i) for i in typeMut])
 

        with open("pairwiseDiff_" + pathwd + ".txt", 'a') as f:
            for j in range(len(typeMut)):
                diffj = [rec, rep, typeMut[j], pop, G]
                diffj = diffj + [int(i[j]) for i in totDiff]
                [f.write(str(i) + " ") for i in diffj]
                f.write("\n")
