## this script calculate the ancestry matrix of a population
## look at all sites (the begin and end of the interval looked at is given as parameter)
## for the sites with a unknown origin 
## (i.e. without mutations in the individual or the mutations appeared before the split)
## their origin will be calculated based on the origin of their nearest neighbor with known ancestry
## note that if both neighbor with known ancestry are farther than a cutoff value,
## the site origin will be kept as unknown

## import the necessary libraries
import sys
import csv
import random
# import matplotlib.pyplot as plt
# import numpy as np

## define the parameters of the simulation that we want to analyse
## note that if we want to define them in the command line
## must use sys.argv[nbr], nbr > 0 (0 being the name of the script)
pathwd = sys.argv[1]
G = sys.argv[2]
cutoff = "1000"
beginInterval = "1"
endInterval = "1600000"

## be sure that G and cutoff are integer
G = int(G)
cutoff = int(cutoff)
beginInterval = int(beginInterval)
endInterval = int(endInterval)

## positions of the sites for which we are interested to look at their ancestry
## if want to look at the whole chromosome, 
## put the beginning at 1
## and the end equals to the length of the chromosome
posMut = list(range(beginInterval, endInterval+1, 100))
## save the (ordered) position in a file (one row)
with open(pathwd + "_positions_" + str(G) + "_" + str(cutoff) + ".txt", "a") as f:
    f.write(pathwd + " " + str(cutoff) + " " + str(G) + " " )
    [f.write(str(i) + " ") for i in posMut]
    f.write("\n")


## note that for the last generation, the full genome was saved -> change the nane
nameFile = pathwd + "/genome2_" + str(G) + ".txt"
if G == 15000:
    nameFile = pathwd + "/fullGenome2_" + str(G) + ".txt"

## load the population information
## separate them into two variable
## one containing the information about each mutation
## the other containing information about each individual
## both can be linked using the mutation id (population-specific)
mutation = []
genome = []
with open(nameFile, 'r') as f:
    for l in f:
        line = l.rstrip('\n').split(" ")
        ## get the mutation information
        ## based on the number of elements in the line
        ## and the first element value
        if (len(line) == 9) & (line[0][0]!='p'):
            subpop = line[6]
            code = '0'
            if (subpop == 'p1') & (int(line[7]) > 5000) & (int(line[7]) < 8000):
                ## the mutation appeared in p1 after the split but before admixture
                code = '1'
            #elif (subpop == 'p1') & (int(line[7]) > 5000) & (int(line[7]) >= 8000):
            #    ## the mutation appeared in p1 after admixture
            #    code = '3'
            elif (subpop == 'p2') & (int(line[7]) > 5000) & (int(line[7]) < 8000):
                ## the mutation appeared in p2 after the split but before admixture
                code = '2'
            #elif (subpop == 'p2') & (int(line[7]) > 5000) & (int(line[7]) >= 8000):
            #    ## the mutation appeared in p2 after admixture
            #    code = '4'
            else:
                code = '5' 
            line.append(code)
            mutation.append(line)
        ## get the individuals information
        if line[0][0] == 'p':
            genome.append(line[2:]) 

## order the mutation based on their position
mutation.sort(key = lambda x: int(x[3])) ## sort by position
## and add their index in the mutation information
mutation = [mutation[i] + [i] for i in range(0,len(mutation))]

## create the necessary variable 
ancestry = []
posAll = []
## extract some information needed thereafter
## about the mutation
## namely their population-specific index, their position and their index in the list of mutations 
id = [i[0] for i in mutation]
pos = [i[3] for i in mutation]
ix0 = [i[10] for i in mutation]


## use a sample of the whole population 
## could also use the whole population (for testing, it's faster)
genomeSample = random.sample(genome,100)

## record the number of individual done
tmp = 0
## calculate the ancestry matrix
for i in genomeSample:
    ## keep the mutations (via their index), that are in the focal indivdual, with known origin only
    ## if the mutation appeared before the split (code = 5)
    ## consider that it has an unknown origin as do not know if came from p1 or p2
    ixMutInd = [j for j in i if mutation[id.index(j)][9] != '5']
    ## get the index (inside the ordered list of mutations) of these mutations
    ix = [mutation[id.index(j)][10] for j in ixMutInd]
    ## idem but get their position
    posix = [int(mutation[id.index(j)][3]) for j in ixMutInd]
    ## ixNeighbor will contain the population-specific index of the mutation nearest neighbor with known origin
    ## if there is none, it will be equal to -1
    ## deltaPos will contain the number of bases between the focal position and the position of the nearest neighbor
    ixNeighbor = []
    deltaPos = []
    ## calculate the distance with the nearest neighbor with known ancestry
    ## look on the left and on the right of the focal position
    for m in posMut:
        m = int(m)
        ixp = [l for l in posix if l-m >= 0]
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
            ## if the distance is above a given threshold, consider that the ancestry is unknown
            ixNeighbor.append(-1) 
            deltaPos.append(min([dp,dm]))
        elif dp > dm:
            ixNeighbor.append(ix[posix.index(max(ixm))])
            deltaPos.append(dm)
        elif dp <= dm:
            ## ARBITRARY DECISION
            ## if dp = dm (ie both left and right nearest neighbor are at the same distance)
            ## choose the right neighbor
            ixNeighbor.append(ix[posix.index(min(ixp))])
            deltaPos.append(dp)
        else:
            print(str(dp) + " " + str(dm) + "PROBLEM!!" + str(tmp) + " " + str(m))
    ancestry.append([mutation[ix0.index(j)][9] if j != -1 else 0 for j in ixNeighbor])
    ## save the ancestry in a file
    ## one individual = one row
    with open(pathwd + "_ancestryMatrix_" + str(G) + "_" + str(cutoff) + ".txt", "a") as f:
        [f.write(str(i) + " ") for i in ancestry[tmp]]
        f.write("\n")
    tmp += 1
    ## save the index of the individual done
    ## to keep track of the time left
    with open(pathwd + "_" + str(G) + "_timeLeft.txt", "a") as f:
        f.write(pathwd + " " + str(cutoff) + " " + str(G) + " " + str(tmp))
        f.write("\n")
    

#anc = np.array(ancestry, dtype = float)
#for i in list(range(0,len(anc[0,:])-10000, 10000)):
#    plt.matshow(anc[:,i:(i+9999)])





        
    
    
    
    
    
    
    
    
    
    
    
