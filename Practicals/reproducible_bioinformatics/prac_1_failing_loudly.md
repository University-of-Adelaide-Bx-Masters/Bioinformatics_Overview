# Prac 1 Failing loudly - handling error messages



#### By Evelyn Collen


## **1. Introduction**

### Why do we need to fail (and loudly)?

If ever you feel a little out of your depth in bioinformatics, and that the commands you run don't seem to always work, I have good news for you! A lot of code works on the basis of 'trial and error' quite literally, since, it is difficult even for the best programmers to Nostradamus their code and know absolutely all possible conditions to account for. Programs that run smoothly are usually thanks to countless developers having dealt with and handled various unexpected bugs previously, and tried to set up their code as robustly as possible.

In diagnostics, every part of the bioinformatics needs to be fully stable so that each time we run a patient sample, we are sure all the correct steps have been applied and we can expect the most "correctly called" genetic result as has been observed during the validation stage. Prior to running patient samples, we validate the test conditions so we know exactly what error rate we can expect from the whole process. Especially in a cancer setting, where we are in an arms-race with finding the best treatment for a cancer vs the cancer's development, we must do this reliably, reproducibly, and quickly so the patient can get their result asap. The best way to achieve this is make sure if any step goes wrong, it is immediately obvious.

In research, reproducibility is also a very important concept which the wonderful Dave will teach you more about. Suffice to say for now, if you contribute some code to research and nobody can run it, and can't understand why, this makes doing science quite difficult. But if you've put some informative error messaging, perhaps they will be able to understand where it went wrong.

In summary, the important thing is take a mindset of embracing failure because it is inevitable - what we need to know is why and when it happens, make our errors loud and obvious, and give a crack at tackling them.


**"If at first you don't succeed, print(error), try, try again" (definitely a word-for-word quote from Edward Hickson 1857)**

### 1.1 Reminder about virtual Machines
By now, you would be used to connecting to your virtual machines (VMs) to complete previous practicals. We are going to be logging in to them and purposefully **crashing some code and "failing" as badly (and loudly!) as possible.**

Just a reminder of info for your [virtual machine](https://en.wikipedia.org/wiki/System_virtual_machine) (VM).
Remember the VM is essentially a program that is running on a server but which behaves as though it is an individual, separate computer within a big computer. 


**Please [go here](../../Course_materials/vm_login_instructions.md) for instructions on connecting to your VM.**

### 1.2 Learning Outcomes

1. Learn about failure and why it's important to handle
2. Write some error messages and make them loud and understandable
3. Relinquishing the fear to try things - just do it, and expect it not to work if you're new at it. It's normal!
4. Staring the problem in the face
5. Learn the basics of the vim text editor
6. Review vcf format

### 1.3 About our dataset

For this practical, we have some human patient vcfs generated from Rhampseq Amplicon sequencing for targeted, 'hot-spot' variant screening. It's called 'hot-spot' as this type of assay is used to cheaply (and quickly!) determine very specific sites in the human genome for screening purposes. 

The vcfs and bams you will be looking at are fully anonymised, and changed a little, but are derived from real life cancer patients who have undergone DPYD pharmacogenetic screening. I will explain more about this in the next prac, but DPYD variants control the amount of DPD enzyme available to metabolise toxic drugs. In some individuals with certain variants, if they are given a 'standard' dose of chemotherapy (usually a fluropyrimidine), this can prove fatal, as they don't have the necessary enzyme to break it down. On the other hand, these chemotherapies are very effective treatments for certain aggressive cancers, and a patient's prognosis may rely on getting the right drug dosage tailored for their personal genetic makeup. 

As always, it is very important to understand your dataset well, to properly do some good bioinformatics to it. Our dataset has been created from rhampseq amplicon data, where reads are amplified across very specific targets and highly duplicated.


### 1.4 Follow the breadcrumb trail of errors 

In this week's (and next week's) practical, we are working on a mini bash pipeline, that runs 3 pipeline scripts for us to help perform fast and reliable DPYD screening. 


The first script, "vcf_validator.py", checks the following cases: 

1. Does the input exist?

2. Is the input a vcf?

3. Is the format of the vcf correct?

The second script, "extract_variants.awk", pulls out some key diagnostic variants for each patient and puts them in a table. It also grabs some key information about each variant that we can use later for some filtering. 

The third script, "generate_DPYD_report.py", generates a nice report that will summarize our screening results for the oncologist. We'll take a look at that in the next prac.

Similar to our simple mini bash pipeline design here, most programs can be distilled into a series of scripts that are chained together and subsetted within each other. If an error occurs in one script, it can propagate through a bunch of other ones that call the script that the error first occured in. Particularly as the scripts get more complex and more intertwined, when an error occurs, you may need to follow the breadcrumbs or 'traceback', through several scripts and functions to pinpoint the offending issue. 

Often when presented with a complex error output, it's a good strategy to work your way through it step-by-step. An in-house pipeline may dump a whole bunch of output full of errors of scripts depending on each other that have broken, and through which you will need to trace through; while a bioinformatics package or program may be neater and only output the error message itself, with little to no traceback. 

In the case of Python scripts, and extending pretty well generally to other scripting languages - often, the traceback will be presented at the top and the actual error itself (causing both its own and other scripts to break) will be near the bottom of the output. This isn't a hard-and-fast rule, but it can help. If you have a particularly wordy and obfuscated traceback or output, try to find exactly the point where the word 'Error' or 'Exception' is. Another helpful hint is to try and look for the name of the particular script you are executing, and what the error message is associated with it. Often, the line number will be given too, so you can open up the script and work out exactly what broke.

When a program or bioinformatics tool breaks, similarly, you will need to carefully check the output, and filter through the error message to find the issue. More often than not, if the script has been written well, it is a simple issue with the inputs/parameters that you are feeding to the script. 

We are going to start by focussing on one script at a time to keep it simple, and then string the 3 scripts together in the next practical. 

### 1.5 Getting ready to look into scripts with vim

First, let's get into debug mode! We're going to be staring at scripts, so let's make sure our text editor (in this case vim) gives us line numbers. You can also use nano if you prefer, but vim can be a bit easier for viewing line numbers.


```bash
vim ~/.vim
```
Press ESC, then "i" key, then type "set number" (note don't type in any of the quotation marks here! They are just to denote what keys to press). Then to exit and save, press ESC, then type ":wq" (which should appear in the bottom left hand corner) then enter. 


load software
```bash
source activate bioinf
```

create all directories and move into project directory
```bash
mkdir -p ~/Practical_Failing_Loudly/{0_scripts,1_vcfs,2_bam,3_reports,4_refs}
cd ~/Practical_Failing_Loudly
```

copy scripts and data and also make symlinks

```bash
cp ~/data/failing_loudly/0_scripts/* 0_scripts/
cp ~/data/failing_loudly/patient_1_dodgy.vcfs 1_vcfs/
ln -s ~/data/failing_loudly/*.vcf 1_vcfs/
ln -s ~/data/failing_loudly/*.bam 2_bam/
ln -s ~/data/failing_loudly/4_refs/* 4_refs/*
```

## **2. Validating input requirements**

### 2.1 Checking if an input is existing

Let's now run just the validator script on one of our vcfs:

```bash
conda activate CHANGEME
python3 ./0_scripts/vcf_validator.py --input_vcf 1_vcfs/patient_1_dody.vcfs
```
You should get something like this:

```
Traceback (most recent call last):
  File "/Users/evelyn/Practicals/failing_loudly/vcf_validator.py", line 67, in <module>
    main(args.input_vcf)
    ~~~~^^^^^^^^^^^^^^^^
  File "/Users/evelyn/Practicals/failing_loudly/vcf_validator.py", line 11, in main
    raise Exception("Unknown Error")
Exception: Unknown Error
```

Oops! There is an error! No worries, let's stare at it. Luckily for us, this is quite a straightforward error message and traceback. You can see here right at the bottom is our error: 'Exception: Unknown Error'. Wow, this is not a helpful error message at all! Who wrote this? We can definitely do better. 

Working our way bottom to top through the traceback, you can see the script that has broken is the one we are running, "vcf_validator.py", at line 11. Going further up, the message then tells us the function the error occured in, function main(), which is called at line 67. 

Now that we've pinpointed our error, let's open the script and look at lines 9-11:

```bash
vim ./0_scripts/vcf_validator.py
```
You should see something like this in your script:

```
9 if not os.path.exists(vcf_file):
10    #Error is output below. Change it to make it more informative!
11    raise Exception("Unknown Error")
```

From the code, we can see it seems the path to the file is missing. We are developing this script, so let's go ahead and improve the error message. While you are still in the vim session, press ESC, then "i" key, then use the arrows to get the cursor to line 9. Change "Unknown Error" to "Error: Could not find any existing vcf file". Press ESC, then ":wq" to save, then enter.

(if you make a mistake in vim and need to exit without saving, just press ESC, then type ":q!", then enter. This will get you out of the file without modifying any content)

Now let's run it again!

```bash
python3 ./0_scripts/vcf_validator.py --input_vcf 1_vcfs/patient_1_dody.vcfs
```
Now your error message should say this: 

```
Traceback (most recent call last):
  File "/Users/evelyn/Practicals/failing_loudly/vcf_validator.py", line 67, in <module>
    main(args.input_vcf)
    ~~~~^^^^^^^^^^^^^^^^
  File "/Users/evelyn/Practicals/failing_loudly/vcf_validator.py", line 11, in main
    raise Exception("Error: Could not find any existing vcf file")
Exception: Error: Could not find any existing vcf file
```
That is a much nicer message and something that is easy to fix. Let's check the names of our vcfs:

```bash
ls 1_vcfs/*vcf*
```
Should result in:
```
patient_1_dodgy.vcfs	
Patient_A.vcf		
Patient_B.vcf		
Patient_C.vcf		
Patient_D.vcf
```

Ah, it looks like I did not spell the vcf correctly - it's patient_1_dod*g*y.vcf, not patient_1_dody.vcf. Let's run it again with the right name:

```bash
python3 ./0_scripts/vcf_validator.py --input_vcf 1_vcfs/patient_1_dodgy.vcfs
```

Now we should have a different error. Yay! A new error! Bioinformaticians often will celebrate when they've been debugging a long time and trying to fix an error - if the message changes, it means we've likely fixed one problem, and can move on to the next one.

Let's go ahead and tackle this new error:

```
Traceback (most recent call last):
  File "/Users/evelyn/failing_loudly/Practicals/./vcf_validator.py", line 67, in <module>
    main(args.input_vcf)
    ~~~~^^^^^^^^^^^^^^^^
  File "/Users/evelyn/failing_loudly/Practicals/failing_loudly/./vcf_validator.py", line 19, in main
    raise ValueError("Unknown Error number 2")
ValueError: Unknown Error number 2
```
Again, another unhelpful error message. Working our way from the bottom, the error is "ValueError: Unknown Error number 2", and again you can see in the traceback this occurs in the script at line 19. 


**Questions:**
1. Based on just this error message, can you work out what line of the script the function main() gets called?
2. We are only running one script here, with one function. How do you think the error message might change if I called this script from another script, let's call it script#2? Would the error message get shorter or longer?
 
<details>
<summary>Answers</summary>
<ul><li>1. Line 67, as can be seen near the top of the traceback</li>
<li>2. The error message would likely get longer because script#2 would also throw an error, which would be added to the overall traceback </li> </ul>
</details>


### 2.2 Checking if the input is a vcf

Most input formats are denoted by their extensions, eg .fastq, .bam, .csv. A good sanity check in your scripts is to check if the extension matches what you expect. Let's open the script again with vim: 

```bash
vim ./0_scripts/vcf_validator.py
```
Our previous error message told us to go to line 19, so let's quickly go there and have a look. You should see this: 

18 if not vcf_file.lower().endswith(".vcf"):
19        raise ValueError("Unknown Error number 2")

Again, this is not very useful information. Press ESC, then "i", then use the arrows to go to line 19 and change "Unknown Error number 2" to something more useful. This time it's your turn to come up with a nice error message that tells us that the vcf extension is not right. Bear in mind when writing a good error message, best practice is to ensure it is providing as much information and context as possible pertaining to the issue, while still being concise and specific, and ideally something actionable. 


When you're done, press ESC, then ":wq" to save your changes, then enter. Run the script again and enjoy a much better articulated error message:

```bash
python3 ./0_scripts/vcf_validator.py --input_vcf 1_vcfs/patient_1_dodgy.vcfs
```
```
Traceback (most recent call last):
  File "/Users/evelyn/failing_loudly/Practicals/./vcf_validator.py", line 67, in <module>
    main(args.input_vcf)
    ~~~~^^^^^^^^^^^^^^^^
  File "/Users/evelyn/failing_loudly/Practicals/failing_loudly/./vcf_validator.py", line 19, in main
    raise ValueError("Error: Input vcf does not have the correct file extension (.vcf)")
ValueError: Error: Input vcf does not have the correct file extension (.vcf)
```
Now that the error reads much nicer, let's name the vcf properly to actually fix the issue:

```bash
mv 1_vcfs/patient_1_dodgy.vcfs 1_vcfs/patient_1_dodgy.vcf
python3 ./vcf_validator.py --input_vcf 1_vcfs/patient_1_dodgy.vcf
```
Woo-hoo! Now onwards to the next error!

```
Traceback (most recent call last):
  File "/Users/evelyn/failing_loudly/Practicals/failing_loudly/./vcf_validator.py", line 79, in <module>
    main(args.input_vcf)
    ~~~~^^^^^^^^^^^^^^^^
  File "/Users/evelyn/failing_loudly/Practicals/failing_loudly/./vcf_validator.py", line 69, in main
    raise Exception("Cannot have data lines preceding header")
```
Ah, this is now interesting. If you look at the bottom of the traceback again, you can see the  Exception causing the new problem. It seems our vcf may have some formatting issues.


## **3. Validating the formatting of the input itself**

### 3.1 Checking if the vcf is correctly formatted

As you know, the VCF format has some key specifications that must be adhered to, and that will crash various bioinformatics tools if not. The 'patient_1_dodgy.vcf' has some key issues that will crash our validator script, so we are going to go ahead and modify this vcf file directly with vim and awk to satisfy the validator. In the real world, you would usually never do this, you would probably just go back to the source of your dodgy vcf and find out what corrupted it in the first instance.

If you need a refresher, here are the vcf file format specifications: 

[VCF file format spec](https://samtools.github.io/hts-specs/VCFv4.2.pdf)


The first error you should notice is "Cannot have data lines preceding header". If you open up the vcf in vim or nano, you may indeed notice an issue, in that one of the lines before the actual main header (starting with #CHROM) is a data line (at line 21, starting with "chr1" and ends with "PASS").

```bash
vim 1_vcfs/patient_1_dodgy.vcf
```

Lines 20-26  look like this:

```
20 ##bcftools_viewVersion=1.10.2+htslib-1.17
 21 chr1    43349345        .       T       A       278     SAMPLE=Patient_A;TYPE=SNV;DP=4197;VD=171;AF=0.0407;BIAS=2:2;REFBIAS=2012:2002;VARBIAS=86:85;PMEAN=20;PSTD=1;QUAL=37.5;QS    TD=1;SBF=1;ODDRATIO=1.00674;MQ=60;SN=27.5;HIAF=0.0399;ADJAF=0;SHIFT3=1;MSI=2;MSILEN=1;NM=1.3;HICNT=165;HICOV=4136;LSEQ=CTGCTGCTGAGGTGGCAGTT;RSEQ=CCTGCACACTACAGGTACCG;GDAMP=1;TL    AMP=1;NCAMP=0;AMPFLAG=0 PASS    GT:DP:VD:AD:AF:RD:ALD:FT        0/1:4197:171:4014,171:0.0407:2012,2002:86,85:PASS
 22 #CHROM
 23 #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Patient_A
 24 #CHROM  POS     ID      REF     ALT     QUAL    INFO    FILTER  FORMAT  Patient_A
 25 ##bcftools_viewCommand=view -O z; Date=Fri May 30 11:19:15 2025
 26 chr1    43349345        .       T       A       278     SAMPLE=Patient_A;TYPE=SNV;DP=4197;VD=171;AF=0.0407;BIAS=2:2;REFBIAS=2012:2002;VARBIAS=86:85;PMEAN=20;PSTD=1;QUAL=37.5;QS    TD=1;SBF=1;ODDRATIO=1.00674;MQ=60;SN=27.5;HIAF=0.0399;ADJAF=0;SHIFT3=1;MSI=2;MSILEN=1;NM=1.3;HICNT=165;HICOV=4136;LSEQ=CTGCTGCTGAGGTGGCAGTT;RSEQ=CCTGCACACTACAGGTACCG;GDAMP=1;TL    AMP=1;NCAMP=0;AMPFLAG=0 PASS    GT:DP:VD:AD:AF:RD:ALD:FT        0/1:4197:171:4014,171:0.0407:2012,2002:86,85:PASS
 ```

You can see line 21 is an erroneous copy of line 26. Let's go ahead and delete line 21. After deleting it, the new line 21 should now be "#CHROM'. 
 

Close and save vim (ESC, then type "wq", then enter), run the validator again, and await the next error message:

```bash
python3 ./vcf_validator.py --input_vcf 1_vcfs/patient_1_dodgy.vcf
```

### 3.2 Fixing the vcf format one error at a time

There are a series of more problems with the vcf. I have given you some instuctions below on how to fix each error, but the errors may not appear in the same order as the order of these solutions, so you will have to go through this list until you find the solution that matches each error. Keep 'fixing' the vcf and matching solutions to the errors, then rerunning the validator script, until no more exceptions are raised. 

**Exception: Header with columns should not be followed by metadata
*Fix*:
The offending line is at line 22 of the vcf, as you can see below.

```
20 ##bcftools_viewVersion=1.10.2+htslib-1.17
 21 #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Patient_A
 22 ##bcftools_viewCommand=view -O z; Date=Fri May 30 11:19:15 2025
```

Metadata should never come after the main header. Move this line up so that it comes just after line 20, so that all the metadata lines go together. After your edit, lines 20-22 should look like this:

```
20 ##bcftools_viewVersion=1.10.2+htslib-1.17
 21 ##bcftools_viewCommand=view -O z; Date=Fri May 30 11:19:15 2025
 22 #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Patient_A
```

**Exception("Too many main headers are present in vcf file")
*Fix*:
We can't have two headers, so go ahead and open the vcf and remove the first header, at line 21. After your edit, your vcf should look like this:

```
20  ##bcftools_viewVersion=1.10.2+htslib-1.17
21 #CHROM  POS     ID      REF     ALT     QUAL    INFO    FILTER  FORMAT  Patient_A
22 ##bcftools_viewCommand=view -O z; Date=Fri May 30 11:19:15 2025
```


*Exception("A column is missing from the vcf main header")
*Fix*:
Use vim or nano to open the vcf file and go to line 21, which just has "CHROM". That's not right as a vcf header need to have at least 8 columns in the main header. Remove this line entirely. After your edits, lines 20-24 should now look like this:
```
20 ##bcftools_viewVersion=1.10.2+htslib-1.17
 21 #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Patient_A
 22 #CHROM  POS     ID      REF     ALT     QUAL    INFO    FILTER  FORMAT  Patient_A
 23 ##bcftools_viewCommand=view -O z; Date=Fri May 30 11:19:15 2025
 24 chr1    43349345        .       T       A       278     SAMPLE=Patient_A;TYPE=SNV;DP=4197;VD=171;AF=0.0407;BIAS=2:2;REFBIAS=2012:2002;VARBIAS=86:85;PMEAN=20;PSTD=1;QUAL=37.5;QS        TD=1;SBF=1;ODDRATIO=1.00674;MQ=60;SN=27.5;HIAF=0.0399;ADJAF=0;SHIFT3=1;MSI=2;MSILEN=1;NM=1.3;HICNT=165;HICOV=4136;LSEQ=CTGCTGCTGAGGTGGCAGTT;RSEQ=CCTGCACACTACAGGTACCG;GDAMP=    1;TL    AMP=1;NCAMP=0;AMPFLAG=0 PASS    GT:DP:VD:AD:AF:RD:ALD:FT        0/1:4197:171:4014,171:0.0407:2012,2002:86,85:PASS
 ```

Exception("The order of columns in the vcf is not correct")
*Fix*: 
The "FILTER" and "INFO" column have gotten swapped somehow. Run this awk command to switch them back:

```bash
awk 'BEGIN{FS=OFS="\t"} /^##/ {print; next} /^#CHROM/ {tmp=$7; $7=$8; $8=tmp; print; next} {tmp=$7; $7=$8; $8=tmp; print}' 1_vcfs/patient_1_dodgy.vcf > 1_vcfs/patient_1_dodgy_fixedcols.vcf && mv 1_vcfs/patient_1_dodgy_fixedcols.vcf 1_vcfs/patient_1_dodgy.vcf
```
Now the order of both the header and the columns of the data section should be switched back in the right order to meet vcf specs. Lines 20-23 should look like this:
```
20 ##bcftools_viewVersion=1.10.2+htslib-1.17
 21 ##bcftools_viewCommand=view -O z; Date=Fri May 30 11:19:15 2025
 22 #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Patient_A
 23 chr1    43349345        .       T       A       278     PASS    SAMPLE=Patient_A;TYPE=SNV;DP=4197;VD=171;AF=0.0407;BIAS=2:2;REFBIAS=2012:2002;VARBIAS=86:85;PMEAN=20;PSTD=1;QUAL    =37.5;QSTD=1;SBF=1;ODDRATIO=1.00674;MQ=60;SN=27.5;HIAF=0.0399;ADJAF=0;SHIFT3=1;MSI=2;MSILEN=1;NM=1.3;HICNT=165;HICOV=4136;LSEQ=CTGCTGCTGAGGTGGCAGTT;RSEQ=CCTGCACACTACAGGTACCG;GD    AMP=1;TLAMP=1;NCAMP=0;AMPFLAG=0 GT:DP:VD:AD:AF:RD:ALD:FT        0/1:4197:171:4014,171:0.0407:2012,2002:86,85:PASS
```

**Question:**
1. Why did we need to use the awk command, instead of just going in with vim or nano and swapping only the column names?

<details>
<summary>Answer</summary>
<ul><li>1. We need to swap the entire column of the data section, not just the header column names. Otherwise, the wrong column will be named incorrectly! Please note also this example is only for teaching purposes. If you ever came across a vcf this messed up, you would to take a break, make a huge cup of coffee, and make a new vcf fresh from the start.</li> </ul>
</details>


## **4. Pull out some QC information about the variants**

In preparation for next week's prac, we are going to run a script that will go through the patient vcfs and pull out information about them. This will be important to ascertain that the variant that has been called isn't just an artefact, and has enough enough depth, the allele frequency we would expect, and passes all filters. All this information is ripe for the plucking straight from the vcf.

Take a look at these four variants commonly screened for in DPYD screening. 

[DPYD variants genomic positions](images_and_refs/DPYD_variants_genome_location.csv)

The name of each variant in the left hand column here is denoted by its HGVS CDNA nomenclature, which describes variants in the context of their effect on the coding region of the transcript. This is commonly used in clinical reports. For example, c.1236G>A is read as cdna, position 1236 of the transcript, sequence change from G to A. 

The right hand column of the file is each diagnostic variant's genomic location. We are going to use this file as a 'key' to pull out info from the vcf, by matching genomic coordinates of these variants to the ones in the vcfs (only if they exist in the vcf) - basically doing a mini annotation ourselves.


Run the following command: 

```bash
cd 0_scripts
awk -f extract_variants.awk ../4_refs/DPYD_variants_genome_location.csv ../1_vcfs/Patient*.vcf > ../3_reports/variant_info.txt
cat ../3_reports/variant_info.txt
cd ..
```

Your output is saved at 3_reports/variant_info.txt and also should be printed to screen. Notice for each of the four diagnostic variants, we have pulled out some information from the vcfs of 3 patients. 


Each of the columns in your output are explained briefly below, and for more information you can search for each one in the vcf specs document:

[VCF file format spec](https://samtools.github.io/hts-specs/VCFv4.2.pdf)

DP : read depth at variant position
FT : sample genotype filter indicating if this genotype was “called” (similar in concept to the FILTER field).
Ref: Reference allele, in this case from hg38 ref
Alt: alternate allele 
GT: Genotype called. 0/0=ref/ref; 0/1=ref/alt; 1/1=alt/alt (assuming we don't have any multi-allelic sites, which we don't here)
FILTER: PASS if this position has passed all filters, i.e., a call is made at this position
AF: allele frequency for each ALT allele in the same order as listed

Note the HGVS nomenclature follows the cdna transcript, and here it is in reverse-complement to the genotypes you see in the vcf because HGVS will always go with the strand direction of the transcript, while the vcf goes in the strand direction of the reference genome (forward strand).

**Questions:**
1. If you list the files in the vcf directory, there is another vcf there, Patient_D.vcf. Why is that vcf missing from variant_info.txt? Why does Patient_B.vcf appear twice?
2. What is the difference between the FILTER column and FT? You may need to refer to the vcf file format document. 
3. Which patient has a homozygous variant?
4. If the HGVS nomenclature and the genotypes fully accord (e.g., both record a substitution of C>T,), which orientation do you expect the transcript to be on (forward or reverse)?
5. Notice all the variant allele frequencies are close to either 0.5 or 1. Why is this, and not something like 0.7?

<details>
<summary>Answers</summary>
<ul><li>1. Patient_D.vcf does not carry a diagnostic variant - that is, the patient has two reference alleles at each of the four sites we tested. </li>
<li>2. FT is the *genotype call* filter, specific for that sample, and indicates the quality of the genotype call. The FILTER column is across all samples in the vcf, and denotes whether the site itself is quality across all present samples (here we only have a single-sample vcf, so more tricky to determine that).</li>
<li>3. Patient A </li>
<li>3. The transcript would be on the foward strand </li>
<li>5. As these are germline variants, they would usually be 0, 1 or 2 copies, inherited equally from mum and dad. So 0 copies is 0% reads are alternate, 1 copy is 50% reads are alternate, 2 copies is 100% reads are alternate. </li> </ul>
</details>


### Bonus tasks if time permitting 

1. The validator script checks a bunch of conditions and raises an error if they are not satisfied. In the first part of the script, if you look at the lines immediately following where these errors are raised, you can see extra 'else - print' statements have been commented out (lines 12 and 20). What do you think these else statements would do? Try removing the comments on the 'elses'/'prints' on lines 12 and 20 and also the series of 'prints' at the bottom of the script, and running the script again on the dodgy vcf: 
```bash
python3 ./vcf_validator.py --input_vcf 1_vcfs/patient_1_dodgy.vcf
```

2. What output did you get? Do you think it is useful to know what checks have been put in place, and that they have passed? Make another directory and write this information to a log. Detailed logs are very useful when debugging as then we can pinpoint what worked and which checks passed before the error. 
```bash
mkdir -p ~/Practical_Failing_Loudly/5_logs
python3 ./vcf_validator.py --input_vcf 1_vcfs/patient_1_dodgy.vcf > 5_logs/vcf_validation.log
```


### Concluding remarks

Hopefully through doing this prac, you can see that through the process of failing loud and proud, we have in the end achieved success! Until next time, when we are going to fail even harder!

