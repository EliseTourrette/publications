#!/usr/bin/python

# calculate the fitness mean, variance, per mutation clas

# calculate the number of differences with the ancestor
# ie the number of mutations per strains
# the mutations 1:7 are the same as the mutations 8:14 (respectively) 
# although they originated in two different populations
# if too slow could remove the mutations 8:14 for the first (should not get them as these are the mutations that arose in pop2)
# thus add the number of mutations when they have the same selection coeff

import sys
import csv

# script = sys.argv[0]

pathwd = sys.argv[1]
rec = sys.argv[2]
rep = sys.argv[3]

ns = 14
typeMut = [ "m" + str(i+1) for i in range(ns) ]
selCoeff = [0, -0.005, -0.002, -0.001, -0.0005, -0.0002, -0.0001]
npop = 2

# diffAncestor = [[[] for i in range(int(ns/2))] for j in range(npop)] ## 1st level = pop, 2nd level = mutation type -> each sub-sub list will n_generation elements

generation = []
population = []
popSize = []

fit = open("fitnessOverall_" + pathwd + ".txt", 'a')

with open(pathwd + "/out.csv") as csvDataFile:
    raw = csv.reader(csvDataFile)
    for row in raw:
        generation.append(row[0])
        population.append(row[1])
        popSize.append(row[2])
        diff = [rec, rep, row[1], row[0], row[3], row[4]]
        [fit.write(str(i) + " ") for i in diff]
        fit.write("\n")



fit.close()

for pop in range(npop):
    pop += 1

    for mut in range(int(ns/2)):
        mut += 1

        nameFile1 = pathwd + "/distribMut" + str(pop) + "_" + str(mut) + ".txt"
        nameFile2 = pathwd + "/distribMut" + str(pop) + "_" + str(int(mut + ns/2)) + ".txt"


        nbMut = []
        gMut = []
        with open(nameFile1, 'rt') as f:
            for l in f:
                m = l.rstrip('\n').split(" ")
                nbMut.append(m[2:])
                gMut.append(m[0])



        nbMut2 = []
        gMut2 = []
        with open(nameFile2, 'rt') as f:
            for l in f:
                m = l.rstrip('\n').split(" ")
                nbMut2.append(m[2:])
                gMut2.append(m[0])



        for k in range(len(nbMut)):
            m = [int(i) + int(j) for i,j in zip(nbMut[k], nbMut2[k])]
            # diffAncestor[pop-1][mut-1].append(sum(m) / len(m)) # mean number of differences for the mutations of type mut and for population pop

            with open("diffAncestor_" + pathwd + ".txt", 'a') as f:
                diff = [rec, rep, typeMut[mut-1], pop, gMut[k], sum(m)/len(m)]
                [f.write(str(i) + " ") for i in diff]
                f.write("\n")


            fit = [(1+selCoeff[mut-1])**i for i in m]
            with open("fitnessPerClass_" + pathwd + ".txt", 'a') as f:
                diff = [rec, rep, typeMut[mut-1], selCoeff[mut-1], pop, gMut[k], sum(fit)/len(fit)]
                [f.write(str(i) + " ") for i in diff]
                f.write("\n")














