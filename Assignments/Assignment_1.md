# Assignment 1

This assignment covers all material in Module 1 - Fundamentals of Bioinformatics. 
It is split into two parts with 30 marks for each part. 

### Submission Details

- All bash/shell scripts should be uploaded as separate `.sh` files. 
- All other answers should be provided in a single pdf. An easy way to do this is to complete your assignment in Word (or another word processor) and export/save your final document as a pdf.  Alternatively, you are welcome to use Rmarkdown, LaTeX, or something similar that can produce pdf output. 
- Please ensure that this pdf contains your name, student number, and the assignment name (Assignment 1).

## Part 1 - Bash 

You will need to produce two valid bash scripts for Part 1 of this assignment.

Q1. Write a bash script to: 
- Download the gff3 file for your assigned species (see bottom of page) to your current directory from Ensembl [1 mark]
- Count how many of each feature type there is, sorted in numerical order [4 marks]
- Count how many gene features have a gene_biotype attribute of protein_coding, and how many have a gene_biotype attribute of lncRNA [4 marks]
- Export the results to a file with a name of the form `my_species_gff_features.txt` where you use your assigned species name instead of my_species [1 mark]. NB: If the species name is missing or incorrect, no marks will be awarded for this section. The script must also include code to generate one or more comment lines in the output file/table before the table with the genome-build used, (hint: `grep` your gff to find the genome build info as the header is very large in most cases)
- Include informative comments in your script to explain each step [1 mark]
- The script must also write the code used to generate the summary (counts) data and count gene features to the output file as part of the file header. [3 marks]

Q2. For the file we used in the practicals (Drosophila_melanogaster.BDGP6.ncrna.fa), add to the final practical script provided so that:
- Two separate output files are now created, titled `dmel_ncrna_summary.txt` and `dmel_ncrna_chromosome_summary.txt` [1 mark]
- For the first output file, ensure that:
  - This output contains a meaningful header [1 mark] 
  - This output contains column names [1 marks]
  - This output includes: a) gene id; b) chromosome; c) start; d) stop; e) strand and f) gene_biotype [3 marks]
- For the second output file, ensure:
  - This output contains a meaningful header [1 mark] 
  - Count the number of gene transcripts assigned to each chromosome [4 marks]
  - Print the ID names for ncRNAs found on the Y chromosome [4]
- Include appropriate comments which make the script easier to understand [1 mark]


Note: If identical comments are identified in any submissions, a mark of **zero** will be given for this question for all suspicious submissions. 
We also strongly, **strongly** suggest completing these scripts without the use of generative AI tools. You will be able to write your assignment scripts by repurposing the code and bioinformatics tools we went over in the first three practicals. 

## Part 2 - Next Generation Sequencing

#### Data
Athaliana reference genome:
https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-51/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

Annotation file: 
https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-51/plants/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.51.gff3.gz

SRR5882797_10M_1.fastq.gz https://adelaideuniversity.box.com/shared/static/egl3n16r0ziaxlvbs9074xqd1liktnuz.gz

SRR5882797_10M_2.fastq.gz https://adelaideuniversity.box.com/shared/static/g2ly4kzz1blus5juy426i37zl45o38pu.gz

## Bash scripts

You will write two bash scripts as described below. Both scripts must use the tools and techniques that were covered in the practicals. 
Clear and concise comments should be included [2 marks].

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



