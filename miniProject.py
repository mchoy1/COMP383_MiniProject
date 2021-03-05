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
