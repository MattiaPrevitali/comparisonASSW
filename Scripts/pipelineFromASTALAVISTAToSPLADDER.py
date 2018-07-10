#!/usr/bin/python

'''
TO DO:
    -   cambiare destinazione file log di SPLADDER
'''

import sys
import subprocess
import re

#--------------------------- CONFIGURATION SETTINGS ----------------------------

# We set chr = 'Y' beacuse we tried ASGAL on chromosome Y
chr = "Y"
# We set nreads = 9000 becasue we used RNASeqReadSim with 9000 reads
nreads = 9000
# Directory with RNASeqReadSim command
RNASEQSIMDIR = "../RNASeqSim/"
# Directory with STAR command
STARDIR = "../STAR/STAR-2.6.0a/bin/Linux_x86_64/"
# Directory with SPLADDER command
SPLADDERDIR = "../SPLADDER/python/"

#------------------------------ INPUT ARGUMENTS --------------------------------

#argv[1] = id_transcript 1      --> transcript_1 id
#argv[2] = id_transcript 2      --> transcript_2 id
#argv[3] = "annotation.gtf"     --> GTF file about annotation
#argv[4] = "genome.fa"          --> FASTA file about reference genome
#argv[5] = flag                 --> 'y' o 'n' for STAR indexes already existing
#argv[6] = indice STAR          --> directory with STAR indexes fro allignment

#--------------------------------- FUNCTIONS -----------------------------------

# Function to extract information from GTF file about a id_transcript
def extractGtf(idTranscr,inputGtfFile):
    gtfList = []
    newGtfList = []
    gtfList = inputGtfFile.split("\n")
    flag = True
    i=0
    p = '^'+chr+'\s.+\s.+\s[0-9]+\s[0-9]+\s.\s.\s.\s.*transcript_id\s\"' + idTranscr + '\"'
    while (flag and i<len(gtfList)):
        obj = re.match(p,gtfList[i])
        if obj:
            p = '^'+chr+'\s.+\s.+\s[0-9]+\s[0-9]+\s.\s.\s.\s.*gene_id\s\"([A-Za-z0-9]+)\"'
            obj = re.match(p,gtfList[i])
            if obj:
                idGene = obj.group(1)
                flag = False
            else:
                print ("ERROR - gene ID not found")
        i = i + 1
    if i >= len(gtfList):
        print ("ERROR - transcript ID (" + idTranscr + ") not found")
    else:
        pGene = '^'+chr+'\s.+\sgene\s[0-9]+\s[0-9]+\s.\s.\s.\s.*gene_id\s\"' + idGene + '\".*'
        pTranscr = '^'+chr+'\s.+\s.+\s[0-9]+\s[0-9]+\s.\s.\s.\s.*transcript_id\s\"' + idTranscr + '\".*'
        for i in gtfList:
            objGene = re.search(pGene,i)
            objTranscr = re.search(pTranscr,i)
            if objGene:
                newGtfList.append(objGene.group())
            elif objTranscr:
                newGtfList.append(objTranscr.group())
    return newGtfList

#----------------------------------- MAIN --------------------------------------

# Read the ID of transcript 1
idTranscr_1 = sys.argv[1]
# Read the ID of transcript 2
idTranscr_2 = sys.argv[2]
# Read the name of annotation (GTF) file
annotationFileName = sys.argv[3]
# Read the name of genome (FASTA) file
genomeFileName = sys.argv[4]

# Create folder to output files
directoryName = idTranscr_1+"vs"+idTranscr_2+"_SPLADDER"
pr1 = subprocess.Popen(["mkdir",directoryName])
pr1.wait()
if pr1.returncode != 0:
    print ("mkdir "+directoryName+" error - returncode = "+str(pr1.returncode))
else:
    gtfDirName = directoryName+"/AnnotationOutputFiles"
    pr2 = subprocess.Popen(["mkdir",gtfDirName])
    pr2.wait()
    if pr2.returncode != 0:
        print ("mkdir "+gtfDirName+" error - returncode = "+str(pr2.returncode))
    else:
        sampleDirName = directoryName+"/SampleFiles"
        pr3 = subprocess.Popen(["mkdir",sampleDirName])
        pr3.wait()
        if pr3.returncode != 0:
            print ("mkdir "+sampleDirName+" error - returncode = "+str(pr3.returncode))
        else:
            bamStarDirName = sampleDirName+"/BAM_"+idTranscr_1+"and"+idTranscr_2
            prBamDir = subprocess.Popen(["mkdir",bamStarDirName])
            prBamDir.wait()
            if prBamDir.returncode != 0:
                print ("mkdir "+bamStarDirName+" error - returncode = "+str(prBamDir.returncode))
            else:
                indexDirName = sampleDirName+"/StarIndexes"
                prIndexDir = subprocess.Popen(["mkdir",indexDirName])
                prIndexDir.wait()
                if prIndexDir.returncode != 0:
                    print ("mkdir "+indexDirName+" error - returncode = "+str(prIndexDir.returncode))
                else:
                    spladderDirName = directoryName+"/SpladderOutputFiles"
                    pr4 = subprocess.Popen(["mkdir",spladderDirName])
                    pr4.wait()
                    if pr4.returncode != 0:
                        print ("mkdir "+spladderDirName+" error - returncode = "+str(pr4.returncode))
                    else:
                        # Read the content of annotation (GTF) file
                        with open(sys.argv[3], 'r') as input_file_gtf:
                          inputGtfFile = input_file_gtf.read()

                        # Read the content of genome (FASTA) file
                        with open(sys.argv[4], 'r') as input_file_seq:
                          inputFastaFile = input_file_seq.read()

                        # Find the ID of gene of transcript 1
                        gtfList_1 = []
                        gtfList_1 = extractGtf(idTranscr_1,inputGtfFile)

                        # Find the ID of gene of transcript 2
                        gtfList_2 = []
                        gtfList_2 = extractGtf(idTranscr_2,inputGtfFile)

                        print ("> Creating gtf file about " + idTranscr_1 + " and " + idTranscr_2 + " ...")
                        # Write on file T1andT2.gtf
                        gtfT1andT2 = open(gtfDirName+"/annotation"+idTranscr_1+"and"+idTranscr_2+".gtf","w")
                        for i in gtfList_1:
                            gtfT1andT2.write(i)
                            gtfT1andT2.write("\n")
                        for i in gtfList_2:
                            gtfT1andT2.write(i)
                            if i != gtfList_2[len(gtfList_2)-1]:
                                gtfT1andT2.write("\n")
                        gtfT1andT2.close()

                        # Run RNASeqReadSim to generate sample about T1andT2
                        print ("> Creating sample file about " + idTranscr_1 + " and " + idTranscr_2 + " ...")
                        pr7 = subprocess.Popen([RNASEQSIMDIR+"RNASeqSim",genomeFileName,gtfDirName+"/"+"annotation"+idTranscr_1+"and"+idTranscr_2+".gtf",sampleDirName+"/sample"+idTranscr_1+"and"+idTranscr_2+".fa",str(nreads)])
                        pr7.wait()
                        if pr7.returncode != 0:
                            print ("RNASeqSim "+idTranscr_1+" and "+idTranscr_2+" - error "+str(pr7.returncode))
                        else:
                            if sys.argv[5]!='y' and sys.argv[5]!='Y':
                                # Create index of genome with STAR
                                print ("> Creating genome index files with STAR about " + idTranscr_1 + " and " + idTranscr_2 + " ...")
                                prIndexStar = subprocess.Popen(["../STAR/STAR-2.6.0a/bin/Linux_x86_64/STAR","--runThreadN","1","--runMode","genomeGenerate","--genomeDir",indexDirName,"--genomeFastaFiles",genomeFileName,"--sjdbGTFfile",annotationFileName,"--sjdbOverhang","100"])
                                prIndexStar.wait()
                                if prIndexStar.returncode != 0:
                                    print ("STAR index - error "+str(prIndexStar.returncode))
                            else:
                                print ("> Genome STAR index files already exists.")
                                indexDirName = sys.argv[6]
                            # Generate BAM file T1+T2 with STAR
                            print ("> Creating sorted BAM file about "+idTranscr_1 +" and "+idTranscr_2+" ...")
                            prBamStar = subprocess.Popen([STARDIR+"STAR","--runThreadN","1","--genomeDir",indexDirName,"--readFilesIn",sampleDirName+"/sample"+idTranscr_1+"and"+idTranscr_2+".fa","--outSAMtype","BAM","SortedByCoordinate","--outFileNamePrefix",bamStarDirName+"/"])
                            prBamStar.wait()
                            if prBamStar.returncode != 0:
                                print ("STAR BAM - error "+str(prBamStar.returncode))
                            else:
                                # Create an index (.bai) of BAM file just created
                                print ("> Creating "+idTranscr_1 +" and "+idTranscr_2+".bai file ...")
                                prBai = subprocess.Popen(["samtools","index",bamStarDirName+"/Aligned.sortedByCoord.out.bam"])
                                prBai.wait()
                                if prBai.returncode != 0:
                                    print ("Samtools index - error "+str(prBai.returncode))
                                else:
                                    # Run SPLADDER tool
                                    print ("> Processing SPLADDER to search new splicing events ...")
                                    prSpladder = subprocess.Popen(["python",SPLADDERDIR+"spladder.py","-a",gtfDirName+"/annotation"+idTranscr_1+"and"+idTranscr_2+".gtf","-b",bamStarDirName+"/Aligned.sortedByCoord.out.bam","-o",spladderDirName,"-T","y","--ignore_mismatches","y"])
                                    prSpladder.wait()
                                    if prSpladder.returncode != 0:
                                        print ("SPLADDER error n. "+str(prSpladder.returncode))
                                    # Print a message to communicate the end of the script
                                    print("\nPipeline completed.")
