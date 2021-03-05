'''COMP 383 Mini-Project'''
#import modules
import os
from Bio import SeqIO
from Bio import Entrez

#log file
log_outfile = open('miniProject.log','w')

currwd=os.getcwd()

#function to take in data & run fastq command
def srrData():
    infile = open('SRRs.txt','r').readlines()
    for srr in infile:
        srr = srr.strip()
        print(srr)
        #file from SRA
        os.system('wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/' + srr+ '/' + srr + '.1')
        #uncompress data
        #convert to paired-end fastq files
        #split reads into two files
        os.system('fastq-dump -I --split-files ' + srr + '.1')

#function for transcriptomes index to find CDS records
def CDSTranscriptomesIndex():
    fasta_outfile = open('EF999921.fasta','w')
    CDS_outfile = open('EF999921_CDS.fasta','w')
    Entrez.email = 'mchoy1@luc.edu'
    #fasta record
    fastaHandle = Entrez.efetch(db = 'nucleotide', id = 'EF999921', rettype = 'fasta')
    fastaRecords = list(SeqIO.parse(fastaHandle,'fasta'))
    fasta_outfile.write('>' + str(fastaRecords[0].description) + '\n' + str(fastaRecords[0].seq))
    fasta_outfile.close()

    #genbank record for CDS
    genbankHandle = Entrez.efetch(db='nucleotide',id = 'EF999921', rettype = '.gb', retmode='text')
    #CDS count
    count = 0
   #loop through genbank records for CDS
    for record in SeqIO.parse(genbankHandle,'genbank'):
        for i in record.features:
            if i.type == 'CDS':
                count += 1
                #extract sequence
                CDS_outfile.write('>' + str(i.qualifiers['protein_id']).replace('[','').replace(']','').replace("'",'') +
                                 "\n" + str(i.location.extract(record.seq)) + '\n')
    log_outfile.write('The HCMV genome (EF999921) has ' + str(count) + ' CDS.')

def kallisto(infile):
   #run kallisto commands
   os.system('time kallisto index -i transcriptome_index.idx EF999921_CDS.fasta')
   os.system('time kallisto quant -i transcriptome_index.idx -o results/' + infile + ' b 30 -t 2 data/' + infile + '_1.fastq data' + infile + '_2.fast2q')

def bowtie2(infile):
    #running bowtie commands
    os.system('bowtie2-build ./EF999921.fasta EF999921_1')
    os.system('bowtie2 --quiet --no-unal --al-conc bowtie2_' + infile
              + '.fastq -x EF999921_1 -1 ' + infile + '.1_1.fastq -2 '
              + infile + '.1_2.fastq -S bowtie_' + infile + '.sam')
    
'''
#count bowtie2 reads
def bowtie2Reads(infile):
    #after filtering
    SRR1_after = open('bowtie2_' + infile + '.1.fastq').readlines()
    SRR2_after = open('bowtie2_ ' + infile + '.2.fastq').readlines()
    donor = ' '
    if infile == 'SRR5660030':
        donor = 'Donor 1 (2dpi)'
    elif infile == 'SRR5660033':
        donor = 'Donor 1 (6dpi)'
    elif infile == 'SRR5660044':
        donor = 'Donor 3 (2dpi)'
    else:
        donor = 'Donor 3 (6dpi)'
    afterBowtie2 =(len(SRR1_after) + len(SRR2_after))/8
    #before filtering
    SRR1_before = open(infile + '._1.fastq').readlines()
    SRR2_before = open(infile + '._2.fastq').readlines()
    beforeBowtie2 = (len(SRR1_before) + len(SRR2_before))/8

    log_outfile.write(donor + ' had ' + str(beforeBowtie2) +
                      ' read pairs before Bowtie 2 filtering and ' +
                      str(afterBowtie2) + ' read pairs after.')
#input for sleuth function
def sleuthInput(infile):
    sleuthOutput = open('SleuthInput.txt','w')
    #columns for file
    sleuthOutput.write('sample' + '\t' +  'condition' + '\t' + 'path' + '\n')
    for srr in infile:
        path = ('/homes/mchoy1/COMP383_MiniProject/' + srr)
        if int(str[3:])%2 == 0:
           #condition 1 (6dpi) is even
           sleuthOutput.write(str(srr) + '\t' + '2dpi' + '\t' + path + '\n')
        else:
            #condition 2 (6dpi) is odd
            sleuthOutput.write(str(srr) + '\t' + '6dpi' + '\t' + path + '\n')

#differiential expression analysis between 2dpi and 6dpi
def sleuth(infile):
    #run sleuth command
    sleuthCommand = 'sleuthScript.R'
    os.system(sleuthCommand)
    sleuthOutfile = open('sleuth_outfile.txt','r').readlines()
    for line in sleuthOutfile:
        log_outfile.write(line)

#SPAdes assembly
def SPAdes(infile):
    srr1 = infile[0]
    srr2 = infile[1]
    srr3 = infile[2]
    srr4 = infile[3]
    #run spades command
    spades = ('spades -k 55,77,99,127 --only-assembler -t 2 -pe1-1 HCMV_' +
              srr1 + '.1.fastq --pe1-2 HCMV_' + srr1 + '.2fastq --pe2-1' +
              srr2 + '.1.fastq --pe2-2 HCMV_' + srr2 + '.2fastq --pe3-1' +
              srr3 + '.1.fastq --pe3-2 HCMV_' + srr3 + '.2fastq --pe4-1' +
              srr4 + '.1.fastq --pe4-2 HCMV_' + srr4 + '.2fastq -o SPAdesAssembly/')
    os.system(spades)
    log_outfile.write(spades)
'''
srrData()
CDSTranscriptomesIndex()
infile = open('SRRs.txt','r').readlines()
for srr in infile:
    srr = srr.strip()
    kallisto(srr)
    bowtie2(srr)
   # bowtie2Reads(srr)
   # sleuthInput(srr)
   # sleuth(srr)
   # SPAdes(srr)
