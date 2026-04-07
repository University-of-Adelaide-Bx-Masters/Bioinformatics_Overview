# Introduction to BLAST
{:.no_toc}

* TOC
{:toc}

# **1. Introduction**
The original BLAST tool [published](https://doi.org/10.1016/S0022-2836(05)80360-2)/released in 1990 is one of the highest cited papers in the literature. The current version has been reengineered and is called BLAST+. NCBI BLAST is a suite of programs that will find local alignments of query sequences to database entries. You can find out about all things BLAST for command line BLAST [here](https://www.ncbi.nlm.nih.gov/books/NBK279690/) or in the [BLAST Book](https://www.oreilly.com/library/view/blast/0596002998/).  

BLAST can align a variety of sequence types:
- nucleotide sequences to nucleotide sequences (BLASTN)
- protein sequences to protein sequences (BLASTP)
- translated nucleotide sequences to protein sequences (BLASTX)
- protein sequences to translated nucleotide sequences (TBLASTN)
- translated nucleotide sequences to translated nucleotide sequences (TBLASTX)

There are a number of parameters that can be varied to increase sensitivity or speed up a search at the expense of sensitivity. There also numerous [parameters](https://www.ncbi.nlm.nih.gov/books/NBK279684/#appendices.Options_for_the_commandline_a) to modify BLAST behaviour and customise output.

## 1.1 Practical Overview
In this practical, we will use the command line version of BLAST to identify matching nucleotide or protein coding regions in a human genomic sequence and in the SwissProt database. 

## 1.2 Learning Outcomes
1. Learn how to generate a BLAST database and run command line BLAST 
2. Gain familiarity with BLAST output formats
3. Investigate the impact of changing BLAST parameters on output and runtime


# **2. Setup**

## 2.1 Activate software
For today's practical, you will need to activate the `bioinf` conda environment:

```bash
source activate bioinf
```

## 2.2 Create directory structure
Let's create a new directory for today's practical and create subdirectories that reflect the main steps in our analysis. This will help us stay organised.
```bash
mkdir -p ~/Practical_blast_intro/{dbs,queries,results}
```

## 2.3 Get data
The files you need are in `~/data/blast_intro/`. See the `README_BLAST.txt` file for descriptions of the data files.

**Copy** the `hg38.fa.gz`, `hg38_reduced.fa.gz` and `uniprot_sprot.fasta.gz` to `~/Practical_blast_intro/dbs` (do not create symlinks as you will need to decompress the files). Copy the remaining files to `~/Practical_blast_intro/queries`.


# **3. Having a BLAST**

## 3.1 Prepare the BLAST databases

BLAST searches a special database of nucleotide sequences that have been broken into `kmers/words` for faster searching (it uses a [hash table](https://en.wikipedia.org/wiki/Hash_table)). It finds matching words in the database and then extends the matches to create a local alignment. You will need to format the BLAST databases that you will use.

You will need to decompress the `hg 38.fa.gz`, `hg38_reduced.fa.gz` and `uniprot_sprot.fasta.gz` files:

```bash
cd ~/Practical_blast_intro/dbs
unpigz -p 2 uniprot_sprot.fasta.gz
unpigz -p 2 hg38.fa.gz
unpigz -p 2 hg38_reduced.fa.gz
```

- Can you think of a way to decompress these files using one command line instead of three?

Once you have done this, you will need to index/format these files so that BLAST can search them.

The syntax for `makeblastdb` is as follows:

`makeblastdb -in <reference.fa> -dbtype nucl -parse_seqids -out <database_name> -title "Database title"`
 
In this example your input file is a nucleotide sequence file `reference.fa` so you need to specify the database type `-dbtype nucl`. You should use the `-parse-seqids` flag as it preserves the sequence identifiers in the input file. The `-out` flag specifies the name for db files and the `-title` flag is optional. 

You will use the following command.  

```bash
makeblastdb -in ~/Practical_blast_intro/dbs/uniprot_sprot.fasta -dbtype 'prot' -parse_seqids -out ~/Practical_blast_intro/dbs/sprot
```

This will generate ten files that BLAST uses.

You will then need to index/format the human hg38 chromosome sequences so that BLAST can search them.

```bash
makeblastdb -in ~/Practical_blast_intro/dbs/hg38.fa -dbtype 'nucl' -parse_seqids -out ~/Practical_blast_intro/dbs/hg38

makeblastdb -in ~/Practical_blast_intro/dbs/hg38_reduced.fa -dbtype 'nucl' -parse_seqids -out ~/Practical_blast_intro/dbs/hg38_reduced
```
Now you are ready to familiarise yourself with command line BLAST. *Once you make the blast dbs you can delete the original fasta files and recover some disk space*

## 3.2 BLAST it!

For quick BLASTN help you can type:

```bash
blastn -help
```

For BLASTX:

```bash
blastx -help
```

Try these to see what the allowed syntax, flags and parameters are.  You will probably want to pipe the output to `less`.

For detailed documentation for all things BLAST see [here](https://www.ncbi.nlm.nih.gov/books/NBK1762/) and for ***detailed command line options and flags*** see [here](https://www.ncbi.nlm.nih.gov/books/NBK279684/#appendices.Options_for_the_commandline_a)

A typical command line might look something like this:

```bash
blastn -query [file.fasta] -task [megablast] -db [database file] -outfmt [0 through 17] -out [outputfile]
```

- For most analyses, I suggest using outfmt 7 or 6, which gives you a tab delimited file. However you should first familiarise yourself with the standard text output.

Output formatting options shown below. Feel free to experiment with these to see how they differ. Some of these are designed for compatibility with downstream applications that require specific input formats. 

> -outfmt `<number>`  
   alignment view options:  
     0 = Pairwise,  
     1 = Query-anchored showing identities,  
     2 = Query-anchored no identities,  
     3 = Flat query-anchored showing identities,  
     4 = Flat query-anchored no identities,  
     5 = BLAST XML,  
     6 = Tabular,  
     7 = Tabular with comment lines,  
     8 = Seqalign (Text ASN.1),  
     9 = Seqalign (Binary ASN.1),  
    10 = Comma-separated values,  
    11 = BLAST archive (ASN.1),  
    12 = Seqalign (JSON),  
    13 = Multiple-file BLAST JSON,  
    14 = Multiple-file BLAST XML2,  
    15 = Single-file BLAST JSON,  
    16 = Single-file BLAST XML2,  
    17 = Sequence Alignment/Map (SAM),  
    18 = Organism Report  

Note that when using output formats >4 some options for choosing the number and description of hits do not work or are incompatible (see below). 

> -num_descriptions `<Integer, >=0>`  
   Number of database sequences to show one-line descriptions for  
   Not applicable for outfmt > 4  
   Default = `500'  
    * Incompatible with:  max_target_seqs  
>    
> -num_alignments `<Integer, >=0>`  
   Number of database sequences to show alignments for  
   Default = `250'  
    * Incompatible with:  max_target_seqs  
>    
> -max_target_seqs `<Integer, >=1>`  
   Maximum number of aligned sequences to keep   
   (value of 5 or more is recommended)  
   Default = `500'  
    * Incompatible with:  num_descriptions, num_alignments  


### 3.2.1 Alignments to Swissprot proteins

You can try the following 

> A basic BLAST command line that gives the basic text output from BLAST, using `blastp` for a protein query to protein db.

```bash
blastp -query ~/Practical_blast_intro/queries/Q967Q8.fasta -task blastp -db ~/Practical_blast_intro/dbs/sprot
```

> A command line that specifies a set of custom columns to report

```bash
blastp -query ~/Practical_blast_intro/queries/Q967Q8.fasta -task blastp -db ~/Practical_blast_intro/dbs/sprot -outfmt "7 delim=  qaccver qlen sallgi sallacc slen pident length mismatch gapopen qstart qend sstart send evalue bitscore"
```

These will of course just dump everything to `stdout` but you know how to cope with that.

> A command line that writes output to a file and specifies the number of threads to run and limits output to 50 targets. For a `blastx` command that uses a nucleotide query against protein db.

```bash
blastx -query ~/Practical_blast_intro/queries/hg38_gene_query.fasta -task blastx -db ~/Practical_blast_intro/dbs/sprot -num_threads 2 -max_target_seqs 50 -out ~/Practical_blast_intro/results/humgene_blastx_sprot.txt -outfmt 7 
```
Call your output file whatever you like, as long as it makes sense to you. 

Once BLASTX has completed you can look at your output using "head", "less", "more" or "cat" or open it with a text editor. 

# **4. Effects of changing parameters**

## 4.1 You can test the effect of `-word_size` on output and speed for `blastn`:

```bash
time blastn -query ~/Practical_blast_intro/queries/hg38_gene_query.fasta -word_size 28 -db ~/Practical_blast_intro/dbs/hg38_reduced -outfmt 7 -out ~/Practical_blast_intro/results/W28.txt
```

```bash
time blastn -query ~/Practical_blast_intro/queries/hg38_gene_query.fasta -word_size 11 -db ~/Practical_blast_intro/dbs/hg38_reduced -outfmt 7 -out ~/Practical_blast_intro/results/W11.txt
```
Be patient, this will take some time (around 8 minutes). While you are waiting, take some time to think about why this is taking so much longer.

## 4.2 You can test the effect of T `-threshold` (this is only for `blastp`) in conjunction with `-word_size`.

- First set T=21 and `-word_size` = 2, we will reduce the wait time by using `-num_threads`=2 .

```bash
time blastp -query ~/Practical_blast_intro/queries/multi-protein_query.fa  -word_size 2 -threshold 21 -db ~/Practical_blast_intro/dbs/sprot -num_threads 2 -outfmt 7 -max_target_seqs 1 -out ~/Practical_blast_intro/results/W2T21multi.txt
```  
- Now increase `-word_size`to 7

```bash
 time blastp -query ~/Practical_blast_intro/queries/multi-protein_query.fa  -word_size 7 -threshold 21 -db ~/Practical_blast_intro/dbs/sprot -num_threads 2 -outfmt 7 -max_target_seqs 1 -out ~/Practical_blast_intro/results/W7T21multi.txt
```
- keep `-word_size` = 7 and lower T to 11

```bash
 time blastp -query ~/Practical_blast_intro/queries/multi-protein_query.fa  -word_size 7 -threshold 11 -db ~/Practical_blast_intro/dbs/sprot -num_threads 2 -outfmt 7 -max_target_seqs 1 -out ~/Practical_blast_intro/results/W7T11multi.txt
```

- Now use a single protein query instead of multi-protein query 

```bash
time blastp -query ~/Practical_blast_intro/queries/Q967Q8.fasta -num_threads 2 -threshold 21 -db ~/Practical_blast_intro/dbs/sprot -outfmt "7 delim=  qaccver qlen sallgi sallacc slen pident length mismatch gapopen qstart qend sstart send evalue bitscore"
```

```bash
time blastp -query ~/Practical_blast_intro/queries/Q967Q8.fasta -num_threads 2 -threshold 11 -db ~/Practical_blast_intro/dbs/sprot -outfmt "7 delim=  qaccver qlen sallgi sallacc slen pident length mismatch gapopen qstart qend sstart send evalue bitscore"
```
- What happened? Is parameter choice important? Why?

## 4.3 You can play with parameter combinations and with search types

Try other parameters for nucleotide comparisons, such as `reward` and `penalty`. These set the score or penalty for match vs mismatch during alignment. Default values for `blastn` are +2/-3. For increased sensitivity when searching less conserved sequences you can try +4/-5. Cost for gaps can be set with `gapopen` and `gapextend` and default values for `blastn` for these are 5 and 2, while default values for `blastp/blastx` are 11 and 1. 

Try using different scoring [matrices](https://www.nature.com/articles/nbt0804-1035) `-matrix` for protein alignments (BLASTP, BLASTX) such as `blosum62` (default) vs `blosum45` or `blosum90`. 
