#!/usr/bin/python

## create a population under free rec
## from the mutation frequencies of one population

import sys
import csv
import random

def prod(x):
    res = 1
    for i in x:
        res *= i
    return res

# script = sys.argv[0]

pathwd = sys.argv[1]
rec = sys.argv[2]
rep = sys.argv[3]

seqG = [5000, 5200, 5500, 5900, 8000, 8300]


nameFile = pathwd + "/distribFit2.txt"
distribFitPop2 = []
with open(nameFile, 'r') as d:
    for l in d:
        tmp = l.rstrip('\n').split(" ")
        if int(tmp[0]) in seqG:
            distribFitPop2.append(tmp[2:])


nameFile = pathwd + "/distribFit1.txt"
distribFitPop1 = []
with open(nameFile, 'r') as d:
    for l in d:
        tmp = l.rstrip('\n').split(" ")
        if int(tmp[0]) in seqG:
            distribFitPop1.append(tmp[2:])

for G in seqG:

    nameFile = pathwd + "/genome1_" + str(G) + ".txt"

    mut = []
    with open(nameFile, 'r') as f:
        for l in f:
            m = l.rstrip('\n').split(" ")
            if (len(m) == 9) & (m[0][0]!='p'):
                mut.append(m)


    N = 100
    Nnew = 10000
    position = list(set([i[3] for i in mut]))
    data = [[] for i in position]
    for i in mut:
        data[position.index(i[3])].extend([1+float(i[4])]*int(Nnew*int(i[8])/N))

    data = [i + [1]*(Nnew-len(i)) for i in data]
    data = [random.sample(i, len(i)) for i in data]
    fitness = []
    for i in list(range(0,Nnew)):
        strain = [j[i] for j in data]
        fitness.append(prod(strain))

    with open("distribFitness_" + pathwd + ".txt", 'a') as f:
        k = [str(rec), str(rep), str(G), "p1", "freeRec"] + fitness
        [f.write(str(i) + " ") for i in k]
        f.write("\n")

    genome = []
    idm = [i[0] for i in mut]
    s = [i[4] for i in mut]
    with open(nameFile, 'r') as f:
        for l in f:
            g = l.rstrip('\n').split(" ")
            if g[0][0]=='p':
                genome.append([1+float(s[idm.index(i)]) for i in g[2:]])


    fitness1 = [prod(i) for i in genome]

    nameFile = pathwd + "/genome2_" + str(G) + ".txt"

    mut = []
    with open(nameFile, 'r') as f:
        for l in f:
            m = l.rstrip('\n').split(" ")
            if (len(m) == 9) & (m[0][0]!='p'):
                mut.append(m)


    N = 100
    Nnew = 10000
    position = list(set([i[3] for i in mut]))
    data = [[] for i in position]
    for i in mut:
        data[position.index(i[3])].extend([1+float(i[4])]*int(Nnew*int(i[8])/N))

    data = [i + [1]*(Nnew-len(i)) for i in data]
    data = [random.sample(i, len(i)) for i in data]
    fitness = []
    for i in list(range(0,Nnew)):
        strain = [j[i] for j in data]
        fitness.append(prod(strain))

    with open("distribFitness_" + pathwd + ".txt", 'a') as f:
        k = [str(rec), str(rep), str(G), "p2", "freeRec"] + fitness
        [f.write(str(i) + " ") for i in k]
        f.write("\n")

    genome = []
    idm = [i[0] for i in mut]
    s = [i[4] for i in mut]
    with open(nameFile, 'r') as f:
        for l in f:
            g = l.rstrip('\n').split(" ")
            if g[0][0]=='p':
                genome.append([1+float(s[idm.index(i)]) for i in g[2:]])


    fitness2 = [prod(i) for i in genome]

    with open("distribFitness_" + pathwd + ".txt", 'a') as f:
        k = [str(rec), str(rep), str(G), "p1", "linkageSample"] + fitness1 ## fitness distribution for the sample population
        [f.write(str(i) + " ") for i in k]
        f.write("\n")
        k = [str(rec), str(rep), str(G), "p1", "linkageTrue"] + distribFitPop1[seqG.index(G)] ## true fitness distribution (for the whole population)
        [f.write(str(i) + " ") for i in k]
        f.write("\n")
        k = [str(rec), str(rep), str(G), "p2", "linkageSample"] + fitness2 ## fitness distribution for the sample population
        [f.write(str(i) + " ") for i in k]
        f.write("\n")
        k = [str(rec), str(rep), str(G), "p2", "linkageTrue"] + distribFitPop2[seqG.index(G)] ## true fitness distribution (for the whole population)
        [f.write(str(i) + " ") for i in k]
        f.write("\n")


    