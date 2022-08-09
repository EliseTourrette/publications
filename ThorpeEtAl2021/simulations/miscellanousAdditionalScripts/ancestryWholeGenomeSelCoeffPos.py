#!/usr/bin/python

## look at the frequency depending on ancestry for the mutations present in the population sample

import sys
import csv

# script = sys.argv[0]

pathwd = sys.argv[1]
rec = sys.argv[2]
rep = sys.argv[3]

seqG = list(range(8000, 15100, 100))

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
    info = [i[3] + "_" + i[1] + "_" + i[4] + "_" + i[8] + "_" + i[9] for i in mutation]

    with open("mutationAncestryProportionSelCoeffPos_" + pathwd + ".txt", "a") as f:
        f.write(str(G) + " ")
        [f.write(str(i) + " ") for i in info]
        f.write("\n")

