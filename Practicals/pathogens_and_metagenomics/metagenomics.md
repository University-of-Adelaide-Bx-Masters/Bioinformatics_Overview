# Metagenomics
{:.no_toc}

* TOC
{:toc}

# **1. Introduction**
In today's practical, you will be investigating metagenomic Illumina sequencing data from a hospital wastewater sample.

## 1.1 Learning Outcomes
1. Understand the main differences between bacterial isolate WGS data and metagenomic WGS data
2. Understand the higher computational resource requiremens when working with metagenomic data
3. Learn how to run a simple metagenomics workflow

# **2. Setup**

## 2.1 Activate software
For today's practical, you will need to activate the `bioinf` conda environment:

```bash
source activate bioinf
```

## 2.2 Create directory structure
Let's create a new directory for today's practical and create subdirectories that reflect the main steps in our analysis. This will help us stay organised.
```bash
mkdir -p ~/Practical_metagenomics/{0_raw,1_species,2_assembly,3_blast,db}
mkdir -p ~/Practical_metagenomics/0_raw/FastQC
mkdir -p ~/Practical_metagenomics/2_assembly/quast
mkdir -p ~/Practical_metagenomics/3_blast/{dbs,queries,results}
```

## 2.3 Get data
The data for today's practical is located in `~/data/metagenomics`. As in previous practicals, we will use symlinks instead of copying large data files.
```bash
# navigate to working directory
cd ~/Practical_metagenomics
# create symlinks for Illumina read files
ln -s ~/data/metagenomics/ERR11269302_[12].fastq.gz 0_raw/
# create symlink for kraken database
ln -s ~/data/dbs/kraken/std_8g db/
# create symlink for AMR gene query sequences
ln -s ~/data/metagenomics/AMR_genes.fasta 3_blast/queries/
```

# **3. Read Quality Control**

Run FastQC on the raw Illumina reads. Be aware that the fastq files are very large and this step will take ~15 minutes.

Look at the FastQC reports and answer the following questions:
- What is the read length?
- What is the total number of bases (including both forward and reverse reads)? If you were sequencing a 5 Mb *E. coli* isolate genome, what coverage would this correspond to? Why do we need sequencing depth to be so much higher for metagenomic sequencing compared to isolate WGS?
- Are there any notable findings in the FastQC reports in terms of sequence quality?

While trimming would be best-practice, we're not going to do that today due to time and resource constraints.

# **4. Species Composition**
There are a number of bioinformatic tools available for investigating species composition of metagenomic sequencing data. Many of these tools have high memory requirements. For this practical we're going to use `kraken2` with a scaled down database that can run on our VMs, as in previous practicals.

Run `kraken2` and `bracken` on the Illumina reads. Be aware that this will take some time to run (~15 minutes). 

Have a look at the inferred species composition in the bracken output file. Remember that you can use `grep` if you want to search for a particular species of interest.

Questions:
- What is the most abundant species?
- What is the second most abundant species? Hint: you may need to use your bash skills to sort the file according to abundance.
- What is the abundance of *Pseudomonas aeruginosa*? What approximate coverage does this equate to?


# **5. Metagenomic Assembly**
Metagenomic assembly can be very computationally intensive. As your VM only has 16 GB RAM, we've run the assembly for you. For your information, the assembly was generated using SPAdes with the following command, which took ~8 hours to run (**DON'T TRY TO RUN THIS ON YOUR VM**):
```bash
spades --meta -t 72 -m 173 -1 ERR11269302_trim_1.fastq.gz -2 ERR11269302_trim_2.fastq.gz -o ERR11269302_metaspades
```
Questions: 
- What does the --meta option do? 
- How many threads were used?
- What does -m 173 mean?

Run the following command to create a symlink to the metagenomic assembly in your directory for today's practical:
```bash
ln -s ~/data/metagenomics/ERR11269302_metaspades_contigs.fasta 2_assembly/
```

## 5.1 Investigating the metagenomic assembly
Let's start by counting the number of contigs in our metagenomic assembly:
```bash 
grep -c "^>" 2_assembly/ERR11269302_metaspades_contigs.fasta
```

Question:
- How does this compare to the number of contigs in our *E. coli* isolate assembly from the Short-read Genome Assembly practical?

We can also see information about contig length and coverage by extracting the contig headers:
```bash
grep "^>" 2_assembly/ERR11269302_metaspades_contigs.fasta | less
```

Questions:
- How variable is the coverage across different contigs? Is this similar or different to what you saw in the *E. coli* isolate assembly? Why?

Now let's run QUAST to look at a range of assembly metrics. We will include our *E. coli* isolate assembly in the QUAST report so that we can compare the isolate assembly and the metagenomic assembly directly:
```bash
quast -t 2 -o 2_assembly/quast/ERR11269302 2_assembly/ERR11269302_metaspades_contigs.fasta ~/Practical_assembly_short/2_assembly/SRR36298124/contigs.fasta
```

Questions:
- In the metagenomic assembly, what is the total number of contigs, and what is the number of contigs >= 500 bp? What does this tell you about the size of most of the contigs in the assembly?
- Ignoring contigs <500 bp, what is the total length of the metagenomic assembly? How does this compare to the total length of the *E. coli* isolate assembly?
- What is the N50 of the metagenomic assembly? How does this compare to the N50 of the *E. coli* isolate assembly?
- Why are the assembly metrics so different between the metagenomic assembly and the isolate assembly?


# **6. BLAST Analysis**
One of the things an assembly allows us to do is to search for specific sequences of interest.
Let's have a look at whether we can find antimicrobial resistance (AMR) genes in our metagenomic dataset. 
We'll begin by creating a BLAST database of our metagenomic assembly:
```bash
makeblastdb -in 2_assembly/ERR11269302_metaspades_contigs.fasta -dbtype nucl -out 3_blast/dbs/ERR11269302_metaspades
```

Now we can query the database using our set of AMR gene query sequences:
```bash
blastn -db 3_blast/dbs/ERR11269302_metaspades -query 3_blast/queries/AMR_genes.fasta -outfmt "7 qlen slen std" > 3_blast/results/ERR11269302_AMR.blast
```

Questions:
- Which AMR genes are present in the metagenomic assembly?
- How does the coverage of the *sul1* contig compare to the coverage of the *dfrA12* contig? What are the possible reasons for this difference?
