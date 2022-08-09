#!/usr/bin/python

# calculate the independant variance
# ie fitness variance assuming that the loci are independant

import sys
import csv

# script = sys.argv[0]

pathwd = sys.argv[1]
rec = sys.argv[2]
rep = sys.argv[3]

def prod(x):
    res = 1
    for i in x:
        res *= i
    return res

ns = 14
nsample = 100
typeMut = [ "m" + str(i+1) for i in range(ns) ]
npop = 2

seqG = range(8000, 15500, 500)

generation = []
population = []
popSize = []

with open(pathwd + "/out.csv") as csvDataFile:
    raw = csv.reader(csvDataFile)
    for row in raw:
        generation.append(row[0])
        population.append(row[1])
        popSize.append(row[2])


for G in seqG:
    
    nameFile = pathwd + "/mutations_" + str(G) + ".txt"
    
    mut1 = []
    mut2 = []
    mut3 = []

    with open(nameFile, 'rt') as f:
        for l in f:
            mut = l.rstrip('\n').split(" ")
            if len(mut) == 12:
                if mut[3] == 'p1':
                    mut1.append(mut)
                elif mut[3] == 'p2':
                    mut2.append(mut)
                elif mut[3] == 'p3':
                    mut3.append(mut)
   
    for pop in range(npop):
        
        print("G: " + str(G) + " pop: " + str(pop))

        pop += 1
    
        if pop == 1:
            mut = mut1

        elif pop == 2:
            mut = mut2

        elif pop == 3:
            mut = mut3

        N = int(popSize[generation.index(str(G)) + pop -1])
        
        pos = [int(i[6]) for i in mut]
        pos = list(set(pos))

        fit = [[] for i in range(len(pos))]

        for i in mut:
            s = 1 + float(i[7])
            n = int(i[11])
            fit[pos.index(int(i[6]))] += [s] * n 
  
        for i in fit:
            if len(i) < N:
	            s = 1
	            n = N - len(i)
	            i += [s] * n

        var = []
        ms = []
        for i in fit:
            m = 0
            v = 0    
            for j in list(set(i)):
                n = i.count(j)
                m += n*j
            m = m/len(i)
            for j in list(set(i)):
                n = i.count(j)
                v += n*(j-m)**2
            v = v/len(i)
            var.append(v)
            ms.append(m**2)


        ssum = [i+j for i,j in zip(var,ms)]
    
        varIndpdt = prod(ssum) - prod(ms) 
        print(varIndpdt)

        with open("indpdtVar_" + pathwd + ".txt", "a") as f:
            indpdtVar = [str(rec), str(rep), str(pop), str(G), str(round(varIndpdt,10))]
            [f.write(i + " ") for i in indpdtVar]
            f.write("\n")
    

  
 

