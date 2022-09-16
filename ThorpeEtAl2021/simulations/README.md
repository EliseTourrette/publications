REQUIREMENTS;
-> simulations: 
	slim v3.6 (https://messerlab.org/slim/)
-> analysis and plotting: 
	python v3.9.6 (libraries sys, csv, random)
	R v4.1.1 (packages ggplot2, cowplot, reshape2, gridExtra)

run time: 
	with the parameters used in the article (population size, mutation rate, number of generations), one simulation needed less than one day to be run
	the analysis was done quasi immediately, except for the calculation of the pairwise differences (time dependent on the sample size)

###################################################################################################
###################################################################################################
## DESCRIPTION DATA
## AND HOW TO OBTAIN IT
## note that the data given in example is not fit to directly use it in the scripts
## it's just there to give an example of the format that is expected
###################################################################################################
###################################################################################################

OUTPUT OF THE SIMULATION
-> distribMutPOP_MUT.txt
-> out.csv
-> mutations_G.txt
-> genomePOP_G.txt
distribFitPOP.txt
IDfixed.txt
fixed.txt
fullGenomePOP_G.txt

DATA NEEDED (FOR THE PLOTS)
-> the RData objects are the output txt files of the same name
sampleData.RData -> generateData.r
dataAncestryCUTOFF_PATHWD.RData -> ancestryFullGenome.py
fitnessOverall.RData -> fitDiffAncestor.py
genome2_G.txt -> output simu
dNdSancestorSumSelCoeff_PATHWD_g8000.RData -> generateData.r -> lastgPOP_MUT.txt -> distribMutPOP_MUT.txt -> output simu
pairwiseDiff.RData -> pairwiseDiff.py


ANALYSES SCRIPTS  (ARGS, INPUT AND OUTPUT)
fitDiffAncestor.py
	arg = pathwd, rec, rep
	input = out.csv [output simu] // distribMutPOP_MUT.txt  [output simu]
	output = fitnessOverall_PATHWD.txt (// diffAncestor_PATHWD.txt // fitnessPerClass_PATHWD.txt)

ancestryFullGenome.py  
	arg = pathwd, rec, rep, g, cutoff
	input = posMut_PATHWD.txt [generateData.r] // genome2_G.txt (or fullGenome2_G.txt, for the last generation)  [output simu]
	output = ancestry_cutoffCUTOFF_PATHWD.txt (// timeLeft_PATHWD.txt // deltaPos_PATHWD_G_CUTOFF.txt)

pairwiseDiff.py
	arg = pathwd, rec, rep
	input = out.csv  [output simu] // genomePOP_G.txt  [output simu]
	output = pairwiseDiff_PATHWD.txt


###################################################################################################
###################################################################################################
## HOW AND WHAT TO RUN
## note that the name of the files as well as the directory in the different scripts need to be changed !!
###################################################################################################
###################################################################################################

## SIMULATION SCRIPT
# to run the simulation script with the command line
# the parameters are defined in the simulation script
# the folder where the simulation will be saved, named folder_simu_rep, need to be created beforehand (but no need to be in the result folder to run the simulations)
mkdir 202105121_1_3
nohup slim -t -d DELTA=0 -d MU=1e-9 -d folder=202105121 -d simu=1 -d rep=3 bottleneckAdmixture_2popBis.slim > 202105121_1_3/out.out &

## ANALYSES
# run using python3
# these script can be run through the command line
# package: sys, random, csv
nohup python3 pairwiseDiff.py "202105122_1_1" "5000bp" 1 > outA211.out  &
nohup python3 fitDiffAncestor.py "202105122_1_1" "5000bp" 1 > outC211.out  &
nohup python3 ancestryFullGenome.py "20201105_1_1" "0bp" 1 15000 1000 > out1.out &

# dN/dS to the ancestor
# only interested in the generation before admixture
# note that these ttwo lines of codes need to be run for every mutation type
# ie for distribMutPOP_MUT.txt, POP = {1,2} and MUT = {1,...,14}
# output: lastgPOP_MUT.txt
# did for population 1
sed '8000q;d' distribMut1_1.txt > tmp11.txt &
tr -s ' '  '\n'< tmp11.txt > lastg1_1.txt &
# and for population 2 
# here the number of the row is different as this population is created only from generation 5000 in the simulation
sed '3001q;d' distribMut2_1.txt > tmp21.txt &
tr -s ' '  '\n'< tmp21.txt > lastg2_1.txt &

## PLOTS
# done using R
# note that this file were not written to be run in the command line
# but to be run in a R session
# package: ggplot2, cowplot, reshape2, gridExtra
# first need to generate the data as RData object using the script
generateData.r
# then plot them using the script
PLOT.r
# and, for the fig4E
plotFig4E.r

## added the additional scripts for the ancestry proportion as a function of the generations in the folder reviews