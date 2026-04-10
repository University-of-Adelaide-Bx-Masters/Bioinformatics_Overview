# **Practical Assignment 2 - Pathogens and Metagenomics**

The data you need for this assignment can be found in `~/data/Assignment2`. This includes:
- Illumina data for 7 *Citrobacter freundi* isolates (Part A)
- A reference *C. freundii* genome, `CP011612.fasta` (Part B)
- A file containing several AMR gene sequences as used in the practicals, `AMR_genes.fasta` (Parts B and C)
- Complete long-read assemblies for two of the *C. freindii* isolates, `isolate1_assembly.fasta` and `isolate4_assembly.fasta` (Part C)

**You must use tools and approaches covered during the practicals. All analysis must be performed on your allocated VM in a directory named `~/Assignment2`.**

Your submission should include the following files:
- Bash scripts named `cfre_qc.sh` and `cfre_phylo.sh` for Parts A and B, respectively
- A single pdf file containing answers to all other questions and screenshots as specified below

## Part A - 26 marks
You have been provided with paired-end Illumina whole-genome sequencing (WGS) data for 7 isolates that have been classified as *Citrobacter freundi*. Your task for Part A is to perform quality control on the WGS data to determine which samples are suitable for downstream analysis. The expected genome size for *C. freundii* is 4.5-5.5 Mb.

1. Write a bash script called `cfre_qc.sh` **[14 marks]** that:
    - Creates a suitable directory structure and creates symlinks for required data files
    - Performs read QC and trimming for each sample
    - Performs taxonomic classification of reads for each sample
    - Generates a genome assembly for each sample
    - Summarises assembly statistics across all samples
    - Includes informative comments
    - Is clear, concise and easy to understand

2. Briefly summarise the main findings from your read QC analysis **[5 marks]**. Make sure you include the following information:
    - What is the read length (or range of read lengths) before and after trimming?
    - What is the expected depth of coverage for each sample after trimming? How did you calculate this value?
    - Are there any concerning quality issues that mean one or more samples should be excluded from downstream analysis, or are all samples suitable to proceed?

3. Based on the results of your taxonomic classification, which samples are NOT consistent with the expectation of a *Citrobacter freundii* isolate? Why? **[3 marks]**

4. Based on your assembly statistics, which samples should be excluded from downstream analysis? Why? Please include a screenshot from your QUAST report showing the key information you rely on for making these conclusions. **[4 marks]**


## Part B - 14 marks
For Part B, you should EXCLUDE any samples that failed your quality control checks from part A.

1. Write a bash script called `cfre_phylo.sh` to generate a phylogeny using all samples that passed your quality control checks from part A, using `CP011612.fasta` as a reference genome. **[5 marks]**

2. Visualise your phylogeny using iTOL. Midpoint root your tree and include a screenshot of your midpoint rooted tree. Answer the following questions: **[4 marks]**
    - a. Which isolates are clonally related?
    - b. What is the length of the branch leading to isolate1? 
    - c. Assuming a genome size of 5 Mb, How many SNPs does this correspond to?

3. Using BLAST, determine whether the *bla*<sub>KPC</sub> gene is present in each of your assemblies and answer the following questions. Include screenshot(s) from your VM showing the key information you rely on for making these conclusions. **[5 marks]**
    - a. Which assemblies contain *bla*<sub>KPC</sub>?
    - b. Which allele is present in each case? 
    - c. What is the length of the contig it is found on in each case? 


## Part C - 12 marks
Long-read data has been generated for two of the *C. freundii* isolates. You have been provided with complete genome assemblies generated from the long-read data (`isolate1_assembly.fasta` and `isolate4_assembly.fasta`).

1. How many plasmids are present in each of the two complete assemblies? Please include screenshot(s) from your VM terminal showing how you determined this. **[4 marks]**

2. Using BLAST, determine whether the *bla*<sub>KPC</sub> gene is present in each of the two complete assemblies and answer the following questions. Include screenshot(s) from your VM showing the key information you rely on for making these conclusions. **[8 marks]**
    - a. How many copies of *bla*<sub>KPC</sub> are present in each complete genome assembly?
    - b. Is this the same as the number of BLAST hits from Part B where you queried the corresponding short-read assembly? Explain why or why not.
    - c. How many different plasmids carry *bla*<sub>KPC</sub> in isolate1 and isolate4?


## Part D - 8 marks

Please answer the following questions about the metagenomics practical:

1. What value did you use for the parameter `-r` when running `bracken`. Why did you choose this value? **[2 marks]**

2. What is the second most abundant species in the metagenomic dataset you analysed? Please include a screenshot from your VM terminal showing how you determined this. **[4 marks]**

3. Which AMR genes are present in the metagenomic assembly? **[2 marks]**
