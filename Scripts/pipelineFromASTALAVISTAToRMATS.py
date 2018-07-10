#!/usr/bin/python

'''
TO DO:
    -   dare la possibilita` di impostare 'paired' o 'single' per rMATS
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
RMATSDIR = "../RMATS/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/"

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
directoryName = idTranscr_1+"vs"+idTranscr_2+"_rMATS"
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
                bamStarDirNameT1 = bamStarDirName+"/BAM_"+idTranscr_1
                prBamDirT1 = subprocess.Popen(["mkdir",bamStarDirNameT1])
                prBamDirT1.wait()
                if prBamDirT1.returncode != 0:
                    print ("mkdir "+bamStarDirNameT1+" error - returncode = "+str(prBamDirT1.returncode))
                else:
                    bamStarDirNameT2 = bamStarDirName+"/BAM_"+idTranscr_2
                    prBamDirT2 = subprocess.Popen(["mkdir",bamStarDirNameT2])
                    prBamDirT2.wait()
                    if prBamDirT2.returncode != 0:
                        print ("mkdir "+bamStarDirNameT2+" error - returncode = "+str(prBamDirT2.returncode))
                    else:
                        indexDirName = sampleDirName+"/StarIndexes"
                        prIndexDir = subprocess.Popen(["mkdir",indexDirName])
                        prIndexDir.wait()
                        if prIndexDir.returncode != 0:
                            print ("mkdir "+indexDirName+" error - returncode = "+str(prIndexDir.returncode))
                        else:
                            rmatsDirName = directoryName+"/rMATSOutputFiles"
                            prDirRMATS = subprocess.Popen(["mkdir",rmatsDirName])
                            prDirRMATS.wait()
                            if prDirRMATS.returncode != 0:
                                print ("mkdir "+rmatsDirName+" error - returncode = "+str(prDirRMATS.returncode))
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

                                gtfFilesName = ["annotation"+idTranscr_1+".gtf","annotation"+idTranscr_2+".gtf","annotation"+idTranscr_1+"and"+idTranscr_2+".gtf"]

                                print ("> Creating gtf file about " + idTranscr_1 + " ...")
                                # Write on file T1.gtf
                                gtfT1 = open(gtfDirName+"/"+gtfFilesName[0],"w")
                                for i in gtfList_1:
                                    gtfT1.write(i)
                                    if i != gtfList_1[len(gtfList_1)-1]:
                                        gtfT1.write("\n")
                                gtfT1.close()

                                print ("> Creating gtf file about " + idTranscr_2 + " ...")
                                # Write on file T2.gtf
                                gtfT2 = open(gtfDirName+"/"+gtfFilesName[1],"w")
                                for i in gtfList_2:
                                    gtfT2.write(i)
                                    if i != gtfList_2[len(gtfList_2)-1]:
                                        gtfT2.write("\n")
                                gtfT2.close()

                                print ("> Creating gtf file about " + idTranscr_1 + " and " + idTranscr_2 + " ...")
                                # Write on file T1andT2.gtf
                                gtfT1andT2 = open(gtfDirName+"/"+gtfFilesName[2],"w")
                                for i in gtfList_1:
                                    gtfT1andT2.write(i)
                                    gtfT1andT2.write("\n")
                                for i in gtfList_2:
                                    gtfT1andT2.write(i)
                                    if i != gtfList_2[len(gtfList_2)-1]:
                                        gtfT1andT2.write("\n")
                                gtfT1andT2.close()

                                sampleFilesName = ["sample"+idTranscr_1+".fa","sample"+idTranscr_2+".fa"]

                                # Run RNASeqReadSim to generate sample about T1
                                print ("> Creating sample file about " + idTranscr_1 + " ...")
                                pr5 = subprocess.Popen([RNASEQSIMDIR+"RNASeqSim",genomeFileName,gtfDirName+"/"+gtfFilesName[0],sampleDirName+"/"+sampleFilesName[0],str(nreads)])
                                pr5.wait()
                                if pr5.returncode != 0:
                                    print ("RNASeqSim "+idTranscr_1+" - error "+str(pr5.returncode))
                                else:
                                    # Run RNASeqReadSim to generate sample about T2
                                    print ("> Creating sample file about " + idTranscr_2 + " ...")
                                    pr6 = subprocess.Popen([RNASEQSIMDIR+"RNASeqSim",genomeFileName,gtfDirName+"/"+gtfFilesName[1],sampleDirName+"/"+sampleFilesName[1],str(nreads)])
                                    pr6.wait()
                                    if pr6.returncode != 0:
                                        print ("RNASeqSim "+idTranscr_2+" - error "+str(pr6.returncode))
                                    else:
                                        if sys.argv[5]!='y' and sys.argv[5]!='Y':
                                            # Create index of genome with STAR
                                            print ("> Creating genome index files with STAR about " + idTranscr_1 + " and " + idTranscr_2 + " ...")
                                            prIndexStar = subprocess.Popen([STARDIR+"STAR","--runThreadN","1","--runMode","genomeGenerate","--genomeDir",indexDirName,"--genomeFastaFiles",genomeFileName,"--sjdbGTFfile",annotationFileName,"--sjdbOverhang","100"])
                                            prIndexStar.wait()
                                            if prIndexStar.returncode != 0:
                                                print ("STAR index - error "+str(prIndexStar.returncode))
                                        else:
                                            print ("> Genome STAR index files already exists.")
                                            indexDirName = sys.argv[6]
                                        # Generate BAM file T1 with STAR
                                        print ("> Creating sorted BAM file about "+idTranscr_1 +" ...")
                                        prBamStarT1 = subprocess.Popen([STARDIR+"STAR","--runThreadN","1","--genomeDir",indexDirName,"--readFilesIn",sampleDirName+"/"+sampleFilesName[0],"--outSAMtype","BAM","SortedByCoordinate","--outFileNamePrefix",bamStarDirNameT1+"/"])
                                        prBamStarT1.wait()
                                        if prBamStarT1.returncode != 0:
                                            print ("STAR BAM T1 - error "+str(prBamStarT1.returncode))
                                        else:
                                            # Generate BAM file T2 with STAR
                                            print ("> Creating sorted BAM file about "+idTranscr_2 +" ...")
                                            prBamStarT2 = subprocess.Popen([STARDIR+"STAR","--runThreadN","1","--genomeDir",indexDirName,"--readFilesIn",sampleDirName+"/"+sampleFilesName[1],"--outSAMtype","BAM","SortedByCoordinate","--outFileNamePrefix",bamStarDirNameT2+"/"])
                                            prBamStarT2.wait()
                                            if prBamStarT2.returncode != 0:
                                                print ("STAR BAM T2 - error "+str(prBamStarT2.returncode))
                                            else:
                                                # Run rMATS tool
                                                print ("> Processing rMATS to search new splicing events ...")
                                                prMATS = subprocess.Popen(["python",RMATSDIR+"rmats.py","--b1",bamStarDirNameT1+"/Aligned.sortedByCoord.out.bam","--b2",bamStarDirNameT1+"/Aligned.sortedByCoord.out.bam","--gtf",gtfDirName+"/"+gtfFilesName[2],"--od",rmatsDirName,"-t","paired","--nthread","1","--readLength","100"])
                                                prMATS.wait()
                                                if prMATS.returncode != 0:
                                                    print ("rMATS error n. "+str(prMATS.returncode))
                                                # Print a message to communicate the end of the script
                                                print("\nPipeline completed.")
