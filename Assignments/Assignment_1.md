# Assignment 1

This assignment covers all material in Module 1 - Fundamentals of Bioinformatics. 
It is split into two parts with 30 marks for each part. 

### Submission Details

- All bash/shell scripts should be uploaded as separate `.sh` files. 
- All other answers should be provided in a single pdf. An easy way to do this is to complete your assignment in Word (or another word processor) and export/save your final document as a pdf.  Alternatively, you are welcome to use Rmarkdown, LaTeX, or something similar that can produce pdf output. 
- Please ensure that this pdf contains your name, student number, and the assignment name (Assignment 1).

## Part 1 - Bash 



## Part 2 - Next Generation Sequencing

You will write two bash scripts (described below) using tools and techniques that were covered in the practicals [2 marks]. Clear and concise comments should be included [2 marks].
#### Data
Athaliana reference genome:
https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-51/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

Annotation file: 
https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-51/plants/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.51.gff3.gz

SRR5882797_10M_1.fastq.gz https://adelaideuniversity.box.com/shared/static/egl3n16r0ziaxlvbs9074xqd1liktnuz.gz

SRR5882797_10M_2.fastq.gz https://adelaideuniversity.box.com/shared/static/g2ly4kzz1blus5juy426i37zl45o38pu.gz

### Script 1 [5 marks]
Script 1 should do the following:
1. Download the genomic sequence (i.e. fasta file) and annotation (i.e. gff file) of the model plant _Arabidopsis thaliana_ to the subdirectory `refs/Athaliana` from the Ensembl ftp directory using the links provided above.  [_2 marks_]
2. Create a genome index for use with `bwa`[_1 marks_]
3. Identify how many chromosomes are in the genome and write this information to standard output (`stdout`) [2 marks]

### Script 2 [13 marks]

Script 2 should do the following: 
1. Download the sequencing data (Illumina paired-end reads) at the links provided to the directory `0_raw/fastq` using `curl` or `wget` and rename these files to the names provided alongside the links.  [3 marks]
2. Trim your data for poor quality bases and remove adapters using `fastp`. Write your trimmed reads to the directory `1_trim/fastq` [3 marks]
3. Align paired-end reads to the genome using `bwa mem`, resulting in a single BAM file in the directory `2_align/bam`. Ensure there are no intermediary SAM (`.sam`) files saved. [3 marks]
4. Sort and index your BAM file [2 marks]
5. Count the number of reads that aligned as proper pairs and output this number to the terminal. [2 marks]

### Short Answer Questions [8 marks]
1. Describe the 4 main components of a FASTQ read/record. [1 mark]
2. Illumina short reads suffer from a deterioration in quality towards the 3' end. Describe the process which causes this. [1 mark]
3. What does it mean if reads align as a proper pair? [1 mark]
4. Illumina short reads may contain portions of adapter sequences at their 3' end. Describe how and why some reads may contain parts of an adapter while others may not. [1 mark]
5. What is the difference between a SAM and a BAM file? [2 marks]
6. Index files are regularly encountered in bioinformatics. For example, `.bai` (and `.csi`) is an index file for BAM files and `.fai` are index files of FASTA files. Describe, in general terms, what index files are and what they facilitate. [1 mark]
7. Why are paired-end Illumina reads usually considered superior to Illumina single-end reads? [1 mark]



