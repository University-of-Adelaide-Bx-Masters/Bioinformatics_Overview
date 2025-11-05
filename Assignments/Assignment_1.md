# Assignment 1

This assignment covers all material in Module 1 - Fundamentals of Bioinformatics. 
It is split into two parts with 30 marks for each part. 

### Submission Details
- All bash/shell scripts should be uploaded as separate `.sh` files. 
- All other answers should be provided in a single pdf. 
- Please ensure that this pdf also contains your name, student number, and the assignment name (Assignment 1).

## Part 1 - Bash 


## Part 2 - Next Generation Sequencing
You will write two bash scripts. 
Your bash scripts should run without errors. This means it should run regardless of where it is located in the directory tree. Meaningful comments are strongly advised [2 marks]
### Script 1 - analyse.sh [? marks]

Using the tools and techniques you have learnt in the course so far, write a script called `analyse.sh` that does the following:
1. Downloads the sequencing data contained at the following links (using curl or wget) to the directory 01_rawData/fastq (note that the filename in the Box link is not informative; make sure you rename your files as below()): [3 marks]

SRR5882797_10M_1.fastq.gz https://adelaideuniversity.box.com/shared/static/egl3n16r0ziaxlvbs9074xqd1liktnuz.gz
SRR5882797_10M_2.fastq.gz https://adelaideuniversity.box.com/shared/static/g2ly4kzz1blus5juy426i37zl45o38pu.gz

2. Trim your data for poor quality bases and remove adapters using fastp. Write your output to the directory `02_trimmedData/fastq` [3 marks]
3. Align paired-end reads to the genome index using bwa mem, resulting in a single .bam file in the directory `03_alignedData/bam`. Ensure there are no intermediary .sam files saved. [3 marks]
4. Sort and index your bam file [2 marks]
5. Your script should be clearly commented [2 marks]
6. Your script should run without errors [2 marks]
### Script 2 - summarise.sh [? marks]
This script should export the following information as a tab-delimited file to the parent directory. Include the category as the first column, and the numbers as the second column
- How many reads were mapped in total? [1 marks]
- How many reads were mapped as a pair? [1 marks]
- How many reads were mapped as a "proper" pair? [1 marks]

### Short Answer Questions
1. Describe the 4 main components of a FASTQ read/record. (1 mark)
2. Illumina short reads suffer from a deterioration in quality towards the 3' end. Describe the process which causes this. (1 mark)
3. Illumina short reads may contain portions of adapter sequences at their 3' end. Describe how and why some reads may contain parts of an adapter while others may not. (1 mark)
4. What is the difference between a SAM and a BAM file? (2 marks)
5. Index files are regularly encountered in bioinformatics. For example, `.bai` (and `.csi`) is an index file for BAM files and `.fai` are index files of FASTA files. Describe, in general terms, what index files are and what they facilitate. (1 mark)
6. Why are paired-end Illumina reads usually considered superior to Illumina single-end reads? (1 mark)




7. Oxford Nanopore (ONT) and Pacific Biosciences (PacBio) are both single molecule, long read sequencing technologies. Describe their error profile and how this differs from that of Illumina short reads. (1 mark)



