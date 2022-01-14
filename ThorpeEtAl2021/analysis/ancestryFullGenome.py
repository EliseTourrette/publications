#!/usr/bin/python

## look at the frequency depending on ancestry for the mutations present in the population sample

import sys
import csv
import random

# script = sys.argv[0]

pathwd = sys.argv[1]
rec = sys.argv[2]
rep = sys.argv[3]
G = sys.argv[4]
cutoff = sys.argv[5]

posMut = []
with open("posMut_" + pathwd + ".txt", 'r') as f:
    for l in f:
        line = l.rstrip('\n')
        posMut.append(int(line))

posMut.sort()

G = int(G)
cutoff = int(cutoff)

mutation = []
genome = []
## can do it with the file fullGenome2 when it's available i.e. for the last generation
nameFile = pathwd + "/genome2_" + str(G) + ".txt"
if G == 15000:
    nameFile = pathwd + "/fullGenome2_" + str(G) + ".txt"

## focus on the mutations that appeared after the split but before the admixture
with open(nameFile, 'r') as f:
    for l in f:
        line = l.rstrip('\n').split(" ")
        if (len(line) == 9) & (line[0][0]!='p'):
            subpop = line[6]
            if (subpop == 'p1') & (int(line[7]) > 5000) & (int(line[7]) < 8000):
                code = '1'
            if (subpop == 'p2') & (int(line[7]) < 8000):
                code = '2'
            else:
                code = '3'
            line.append(code)
            mutation.append(line)
        if line[0][0] == 'p':
            genome.append(line[2:]) 

mutation.sort(key = lambda x: int(x[3]))
mutation = [mutation[i] + [i] for i in range(0,len(mutation))]
ancestry = []
posAll = []
## extract some information needed thereafter
id = [i[0] for i in mutation]
pos = [i[3] for i in mutation]
ix0 = [i[10] for i in mutation]

tmpi = 0
## use a sample of the whole population (if has the whole population i.e at the last generation, else use the given genome
genomeSample = genome
if G == 15000:
    genomeSample = random.sample(genome,1000)

for i in genomeSample:
    ## extract the 'code' for the mutation with a known origine
    ibis = [j for j in i if mutation[id.index(j)][9] != '3']
    ## get the index (inside the ordered list of mutations) of these mutations
    ix = [mutation[id.index(j)][10] for j in ibis]
    ## idem but get their position
    posix = [int(mutation[id.index(j)][3]) for j in ibis]
    ## sites to be done (that are not already in ibis), identified via their position
    # sites = [i for i in posMut if posix.count(int(i)) == 0]
    cn = []
    deltaPos = []
    ## save the delta used and print them for later checks
    for m in posMut:
        m = int(m)
        ixp = [l for l in posix if l-m > 0]
        ## if consider that ALL sites are of unknown ancestry, remove the case of l-m = 0: l-m > 0
        if len(ixp) > 0: 
            dp = min(ixp) - m
        else:
            dp = 1600000
        ixm = [l for l in posix if m-l > 0]
        if len(ixm) > 0: 
            dm = m - max(ixm) 
        else:
            dm = 1600000
        if min([dp,dm]) > cutoff:
            cn.append(-1) 
            deltaPos.append(min([dp,dm]))
            ## append a non existing index if the distance between the site and its closest neighbour with known ancestry is more than the cutoff value 
        elif dp > dm:
            cn.append(ix[posix.index(max(ixm))])
            deltaPos.append(dm)
        else:
            cn.append(ix[posix.index(min(ixp))])
            deltaPos.append(dp)
    ancestry.append([mutation[ix0.index(j)][9] if j != -1 else 0 for j in cn])
    tmpi += 1
    with open("deltaPos_" + pathwd + "_" + str(G) + "_" + str(cutoff) + ".txt", "a") as f:
        [f.write(str(i) + " ") for i in deltaPos]
        f.write("\n")
    with open("timeLeft_" + pathwd + ".txt", "a") as f:
        f.write(str(cutoff) + str(G) + str(tmpi))
        f.write("\n")

ancestry = [[i for i in r] for r in zip(*ancestry)]
## only need to know how many genomes have ancestry = 1
## ancestry = 2 will just be the total number of genomes minus the number of ancestry 1 and ancestry not known
ancestry1 = [i.count('1') for i in ancestry]
ancestry2 = [i.count('2') for i in ancestry]

with open("ancestry_cutoff" + str(cutoff) + "_" + pathwd + ".txt", "a") as f:
    f.write(str(G) + " " + "from_pop1" + " ")
    [f.write(str(i) + " ") for i in ancestry1]
    f.write("\n")
    f.write(str(G) + " " + "from_pop2" + " ")
    [f.write(str(i) + " ") for i in ancestry2]
    f.write("\n")
    f.write(str(G) + " " + "positions" + " ")
    [f.write(str(i) + " ") for i in posMut]
    f.write("\n")


