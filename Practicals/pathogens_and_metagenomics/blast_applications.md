# BLAST Applications
{:.no_toc}

* TOC
{:toc}

# **1. Introduction**
In the previous practical you learned the basics of running command line BLAST. In this practical we're going to apply those skills to investigating pathogen genomes.

## 1.1 Learning Outcomes
1. Learn how to use BLAST to determine whether a given query gene is present in a bacterial genome assembly
2. Learn how to use BLAST to characterise sequence rearrangements


# **2. Setup**

## 2.1 Activate software
For today's practical, you will need to activate the `bioinf` conda environment:

```bash
source activate bioinf
```

## 2.2 Create directory structure
Let's create a new directory for today's practical and create subdirectories that reflect the main steps in our analysis. This will help us stay organised.
```bash
mkdir -p ~/Practical_blast/{genomes,dbs,queries,results,contigs}
```

## 2.3 Get data
In today's practical, we will investigate some of the assemblies we generated in the previous topic.
The other data you will need is located in `~/data/blast`. As in previous practicals, we will use symlinks instead of copying large data files.
```bash
# navigate to working directory
cd ~/Practical_blast
# create symlinks for query sequences
ln -s ~/data/blast/AMR_genes.fasta queries/
ln -s ~/data/blast/blaKPC-2.fasta queries/
# create symlinks for complete genomes
ln -s ~/data/blast/CAV*.fasta genomes/
# create symlinks for E. coli assemblies. Note you may need to modify the paths depending on where you saved your assemblies.
ln -s ~/Practical_assembly_short/2_assembly/SRR36298124/contigs.fasta genomes/SRR36298124_contigs.fasta
ln -s ~/Practical_assembly_long/1_assembly/SRR36505805/assembly.fasta genomes/SRR36505805_assembly.fasta
```

If you do not have the assemblies from the previous topic, they are also available in the `~/data/blast` directory. In this case you can replace the last two commands above with the following:
```bash
# create symlinks for E. coli assemblies
ln -s ~/data/blast/SRR36298124_contigs.fasta genomes/
ln -s ~/data/blast/SRR36505805_assembly.fasta genomes/
```


# **3. Searching for AMR genes**
We will use BLAST to determine which antimicrobial resistance (AMR) genes are present in a bacterial genome assembly. To do this, we will query our genomes using the file `queries/AMR_genes.fasta`, which contains just a few select AMR genes that are relevant for our analysis today. Let's begin by taking a look at the contents of this file:
```bash
grep "^>" queries/AMR_genes.fasta
```
You can see that there are only a small number of AMR gene sequences in this file. For a real-world analysis, you would likely include many more AMR gene sequences.


## 3.1 Querying our *E. coli* short-read assembly
We will begin by using BLAST to look for AMR genes in the *E. coli* genome assembly from short-read Illumina data that we generated in the previous topic.

First, we need to create a BLAST database:
```bash
makeblastdb -in genomes/SRR36298124_contigs.fasta -dbtype nucl -out dbs/Ecol_Illm
```

Now we can query our genome:
```bash
blastn -db dbs/Ecol_Illm -query queries/AMR_genes.fasta -outfmt "7 qlen slen std" > results/Ecol_Illm_AMR.blast
```

Let's take a look at the results:
```bash
less results/Ecol_Illm_AMR.blast
```

Questions:
- What does `-outfmt "7 qlen slen std"` mean in the `blastn` command?
- Which AMR genes are present in this genome assembly?
- What size are the corresponding contigs?


## 3.2 Querying our *E. coli* long-read assembly
Next, let's try running the same query, but using the long-read assembly we generated for the same *E. coli* isolate, rather than the short-read assembly.

Create a BLAST database for the ONT assembly and run `blastn` using the same `AMR_genes.fasta` query.

Questions:
- Which AMR genes are present in this genome assembly? Is this the same as for the short-read assembly above?
- How many copies of each AMR gene are present in this genome assembly? Is this the same as for the short-read assembly? Why?


# **4. Investigating plasmid rearrangements using BLAST**
In this section, you will analyse a number of complete genome assemblies generated from long-read sequencing. These are:

| **Isolate Name** | **Species** |
|------------------|-------------|
| CAV1016 | *Klebsiella pneumoniae* |
| CAV1151 | *Phytobacter ursingii* |
| CAV1392 | *Klebsiella pneumoniae* |
| CAV1668 | *Enterobacter hormaechei* |
| CAV1669 | *Enterobacter hormaechei* |
| CAV1311 | *Enterobacter hormaechei* |

All of these isolates carry the *bla*<sub>KPC</sub> gene, and we will focus on the *bla*<sub>KPC</sub> plasmids.

## 4.1 Finding the *bla*<sub>KPC</sub> contigs
First, have a look at the genome assembly for one of the isolates:
```bash
grep "^>" genomes/CAV1016.fasta
```
You can see that this isolate has a complete (circular) chromosome, as well as three plasmids.

Now let's create a BLAST database for this isolate:
```bash
makeblastdb -in genomes/CAV1016.fasta -dbtype nucl -out dbs/CAV1016
```
And query the BLAST database using our `AMR_genes.fasta` file, which contains the *bla*<sub>KPC</sub> gene:
```bash
blastn -db dbs/CAV1016 -query queries/AMR_genes.fasta -outfmt "7 qlen slen std" > results/CAV1016_AMR.blast
```
You can view the results using `cat` or `less`.

Questions:
- How many copies of *bla*<sub>KPC</sub> are present in the assembly?
- Which plasmid contig(s) contain *bla*<sub>KPC</sub>?
- Which *bla*<sub>KPC</sub> allele(s) are present?

Write a script to repeat the above steps on the remaining isolates.

## 4.2 Extracting the *bla*<sub>KPC</sub> contigs
Next, we will extract the contigs containing *bla*<sub>KPC</sub> for further comparison.

If we know which line numbers represent our contigs of interest, we can easily extract these using `sed`. For example:
```bash
grep -n "^>" genomes/CAV1016.fasta
```
gives the following output:
```
1:>CP017934.1 Klebsiella pneumoniae strain CAV1016, complete sequence
76970:>CP017935.1 Klebsiella pneumoniae strain CAV1016 plasmid pCAV1016-90, complete sequence
78253:>CP017936.1 Klebsiella pneumoniae strain CAV1016 plasmid pCAV1016-76, complete sequence
79344:>CP017937.1 Klebsiella pneumoniae strain CAV1016 plasmid pKPC_UVA01, complete sequence
```
From this, we can see that the contig we want to extract (including header) runs from line 79344 until the end of the file. We can find the number of lines in the file using `wc -l`. Therefore, we can extract the desired contig as follows:
```bash
sed -n '79344,79969 p' genomes/CAV1016.fasta > contigs/CAV1016_pKPC_UVA01.fasta
```

Extract the *bla*<sub>KPC</sub> contigs from the other assemblies in a similar manner.

## 4.3 Comparing the *bla*<sub>KPC</sub> contigs
Now that we've extracted the *bla*<sub>KPC</sub> contigs, we can use BLAST to compare them!

First, create a blast database for each of the extracted contigs.

Now run the following BLAST queries. In each case, interpret the output. Specifically:
- Do the two contigs align over their entire length?
- If yes, are there SNP or structural differences between the two contigs?
- Do the two contigs represent identical plasmids in different isolates?
- If the contigs do not represent identical plasmids, find the alignment that includes the *bla*<sub>KPC</sub> gene. How long is this alignment? What does this tell you about shared sequence between the two contigs?

Queries to run:
- CAV1016 plasmid against CAV1151 plasmid
- CAV1016 plasmid against CAV1392 plasmid
- CAV1392 plasmid against CAV1392 chromosome
- CAV1016 plasmid against CAV1668 plasmid
- CAV1016 plasmid against CAV1669 plasmid
- CAV1669 plasmid against CAV1311 plasmid
