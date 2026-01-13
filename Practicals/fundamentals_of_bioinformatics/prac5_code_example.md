
### Your script should now look something ike the script below. 

Take a look through it, noting which commands are in which sections and making sure you understand why they're there and what they do. 


```bash
#!/bin/bash

# Variables

# Load software
source activate bioinf

# Create directories
mkdir --parents ~/Practical_alignment/{ref,0_raw,1_trim,2_align,3_variants}
mkdir -p ~/Practical_alignment/0_raw/FastQC
mkdir -p ~/Practical_alignment/1_trim/{fastp,FastQC}

# Move into the Practical directory
cd ~/Practical_alignment

# Get data
cp  ~/data/intro_ngs/chr2_sub.fa ref/
ln -s ~/data/intro_ngs/*.fq.gz 0_raw/

# Assess raw read quality with fastqc
fastqc -o 0_raw/FastQC -t 2 0_raw/ERR3241917_*.fq.gz

# Trim with fastp
fastp --thread 2 \
-i 0_raw/ERR3241917_1.fq.gz \
-I 0_raw/ERR3241917_2.fq.gz \
-o 1_trim/ERR3241917_1.fq.gz \
-O 1_trim/ERR3241917_2.fq.gz \
--unpaired1 1_trim/ERR3241917_1_orphan.fq.gz \
--unpaired2 1_trim/ERR3241917_2_orphan.fq.gz \
--cut_right \
--cut_window_size 4 \
--cut_mean_quality 20 \
--length_required 75 \
--html 1_trim/fastp/ERR3241917_fastp.html

# Assess read quality after trimming with fastqc
fastqc -o 1_trim/FastQC -t 2 1_trim/FastQC/ERR3241917_*.fq.gz

```

Once you've understood what's going on, go back to the main practical. 
 
