#!/usr/bin/python

'''
TO DO:
    -   ciclo per run di RNASeqSim per evitare ripetizione di codice
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
# Directory with ASGAL command
ASGALDIR = "../ASGAL/galig/"

#------------------------------ INPUT ARGUMENTS --------------------------------

#argv[1] = id_transcript 1      --> transcript_1 id
#argv[2] = id_transcript 2      --> transcriot_2 id
#argv[3] = "annotation.gtf"     --> GTF file about annotation
#argv[4] = "genome.fa"          --> FASTA file about reference genome

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
directoryName = idTranscr_1+"vs"+idTranscr_2+"_ASGAL"
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
            asgalDirName = directoryName+"/AsgalOutputFiles"
            pr4 = subprocess.Popen(["mkdir",asgalDirName])
            pr4.wait()
            if pr4.returncode != 0:
                print ("mkdir "+asgalDirName+" error - returncode = "+str(pr4.returncode))
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

                sampleFilesName = ["sample"+idTranscr_1+".fa","sample"+idTranscr_2+".fa","sample"+idTranscr_1+"and"+idTranscr_2+".fa"]

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
                        # Run RNASeqReadSim to generate sample about T1andT2
                        print ("> Creating sample file about " + idTranscr_1 + " and " + idTranscr_2 + " ...")
                        pr7 = subprocess.Popen([RNASEQSIMDIR+"RNASeqSim",genomeFileName,gtfDirName+"/"+gtfFilesName[2],sampleDirName+"/"+sampleFilesName[2],str(nreads)])
                        pr7.wait()
                        if pr7.returncode != 0:
                            print ("RNASeqSim "+idTranscr_1+" and "+idTranscr_2+" - error "+str(pr7.returncode))
                        else:

                            # Run all possible (9) combinations about gtf and sample files
                            listProcesses = [[1,2,3],[4,5,6],[7,8,9]]
                            for i in range(len(gtfFilesName)):
                                for j in range(len(sampleFilesName)):
                                    #print("\n\n\n\ni = "+str(i)+" ,j = "+str(j)+"\n\n\n\n")
                                    print ("> Processing ASGAL to search new splicing events from "+gtfFilesName[i]+" and "+sampleFilesName[j]+" ...")
                                    listProcesses[i][j] = subprocess.Popen([ASGALDIR+"asgal","-g",genomeFileName,"-a",gtfDirName+"/"+gtfFilesName[i],"-s",sampleDirName+"/"+sampleFilesName[j],"-o",asgalDirName+"/outputASGAL_"+gtfFilesName[i]+"_"+sampleFilesName[j]])
                                    listProcesses[i][j].wait()
                                    if listProcesses[i][j].returncode != 0:
                                        print ("ASGAL error n. "+str(listProcesses[i][j].returncode)+" on "+gtfFilesName[i]+" and "+sampleFilesName[j])
                            # Print a message to communicate the end of the script
                            print("\nPipeline completed.")
