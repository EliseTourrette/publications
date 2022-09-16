## python version 3.10.0

## read the alignment
## filter the CDS sequence 
## then go back to R
## to read the ref sequence and this list of CDS
## and to continue the analysis

## read the sequence of the core ref genome
## as we want only the CDS that are present in the core genome
genome = []
i = 0
with open("alignments/alignment_core_ref/26695_core_aligned.aln", 'rt') as f:
    for l in f:
        if(i % 2 == 1):
            genome.append(l)
        i += 1

## read the annotation file
## to get the positions of the CDS
## filter them to only keep the ones in the core genome
## reuse what I did in the 1st point
coreCDS = []
ref = []
i = 0
with open("gff/core_CDS.gff", 'rt') as f:
    for l in f:
        seq = l.rstrip("\n").split(" ")
        if(len(l) == 1667868):
            ref.append(l)
        elif((i > 1)):
            coreCDS.append(seq)
        i +=1

seqCDS = [list(genome[0][(int(i[3])-1):int(i[4])]) if i[6] == '+' else [genome[0][j] for j in range(int(i[4])-1, (int(i[3])-2), -1)] for i in coreCDS]

with open("sequence_CDS_ref.fa", "a") as f:
    [f.write(''.join(i) + "\n") for i in seqCDS]

