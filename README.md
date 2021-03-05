# COMP383_MiniProject
In order to run this script, clone this github repository into your workspace:

git clone https://github.com/mchoy1/COMP383_MiniProject.git

Software Required: -Python3 -Linux -os -Biopython -Kallisto -R Studio -Bowtie2 -SPAdes

Main file is the miniProject.py (code pipeline)

Project Instructions:

1) We would like to compare HCMV transcriptomes 2- and 6-days post-infection (dpi). First, retrieve the following transcriptomes from two patient donors from SRA and convert to paired-end fastq files. You can use wget (by constructing the path based on the SRR numbers for each of these samples). Donor 1 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896360 Donor 1 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896363 Donor 3 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896374 Donor 3 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896375

2) We will quantify TPM in each sample using kallisto, but first we need to build a transcriptome index for HCMV (NCBI accession EF999921). Use Biopython to retrieve and generate the appropriate input and then build the index with kallisto. You will need to extract the CDS features from the GenBank format. Write the following to your log file: The HCMV genome (EF999921) has (count # extracted) CDS.

3) Quantify the TPM of each CDS in each transcriptome using kallisto and use these results as input to find differentially expressed genes between the two timepoints (2pi and 6dpi) using the R package sleuth. Write the following details for each significant transcript (FDR < 0.05) to your log file, include a header row, and tab-delimit each item: target_id test_stat pval qval

4) It has been proposed that HCMV disease and pathogenesis may be related to the genetic diversity of the virus (Renzette et al. https://www.ncbi.nlm.nih.gov/pubmed/25154343/). Which publicly available strains are most similar to these patient samples? To compare to other strains, we will assemble these transcriptome reads. We don’t expect assembly to produce the entire genome, but enough to be useful in BLAST. Virus sequencing experiments often include host DNAs. It is difficult to isolate the RNA of just the virus (as it only transcribes during infection of the host cell). Before assembly, let’s make sure our reads map to HCMV. Using Bowtie2, create an index for HCMV (NCBI accession EF999921). Next, save the reads that map to the HCMV index for use in assembly. Write to your log file the number of reads in each transcriptome before and after the Bowtie2 mapping. For instance, if I was looking at the Donor 1 (2dpi) sample, I would write to the log (numbers here are arbitrary): Donor 1 (2dpi) had 230000 read pairs before Bowtie2 filtering and 100000 read pairs after.

5) Using the Bowtie2 output reads, assemble all four transcriptomes together to produce 1 assembly via SPAdes. Write the SPAdes command used to the log file.

6) Write Python code to calculate the number of contigs with a length > 1000 and write the # to the log file as follows (replace # with the correct integer): There are # contigs > 1000 bp in the assembly.

7) Write Python code to calculate the length of the assembly (the total number of bp in all of the contigs > 1000 bp in length) and write this # to the log file as follows (replace # with the correct integer): There are # bp in the assembly.
