# Reproducibility and Failing loudly - handling error in clinical screening

#### By Evelyn Collen


## **1. Introduction**

### What happens when we can't afford to fail - DPYD pharmacogenomic screening 


In last week's practical, we started to get the hang of handling failing scripts and the concept of failure in general. Here we are going to focus on how we can best handle errors in diagnostic reporting. 

The DPYD gene is responsible for generating the dihydropyrimidine dehydrogenase (DPD) enzyme, which plays a key role in the metabolism of toxic compounds. Deficiency in this enzyme can cause fatal toxicity to fluoropyrimidine chemotherapy treatments (e.g., 5-fluorouracil, capecitabine), which are widely used in the treatment of solid tumours such as colorectal cancer, breast cancer and gastrointestinal cancers. Variants in the DPYD can cause different functionality in the DPD enzyme's function, categorised into zero, decreased, or normal function. 

Getting the variant-tailored dosage just right is absolutely crucial, as the standard dose that is effective for some people can be fatal for certain variant carriers. Severe toxicity occurs in about 10% to 40% of patients, and around 7% of Europeans carry some variant that impaird function. In South Australia, there has unfortunately been a recorded case of patient fatality due to DPD enzyme deficiency and an incorrectly tailored dosage. 


![alt text](dpyd2.jpg)

Figure 1. The relative risk of toxicity in those with versus without the
specified variant when treated with full or individualised dose (taken from [Sonic Genetics](https://www.sonicgenetics.com.au/))
 
 
Depending on the configuration of alleles, and whether those alleles have zero, decreased or normal function, patients will receive a metabolism rating. Have a look at how the rating is worked out in this [DPYD metabolism rating table](4_refs/DPYD_metabolism_rating_and_recommendations.xlsx)


### 1.1 Reminder about virtual Machines

As usual we will be connecting the virtual machines, 

**Please [go here](../../Course_materials/vm_login_instructions.md) for instructions on connecting to your VM.**

## 1.2 Learning Outcomes

1. Understand clinical failure in a clinical bioinformatics and how to minimise it
2. Learn ways to minimise failure in 
4. 


## Streamlining DPYD screening

In this scenario, you are a clinical bioinformatician tasked with finding out whether some cancer patients are carrying particular variants in a DPYD gene, and what phenotype and metabolism rating can be deduced from their allele types. It's absolutely critical we get this right as the oncologist will use this information to determine the best course of treatment that can aggressively handle the cancer, whilst minimising the effect on the patient as much as possible. 
The patient is waiting on this information - so not only we have to get it right, we want do it quickly. Last week we made sure our errors were noisy in our scripting; this week we are going to make our pipeline noisy about the whole process. 

We also want the reporting process to be streamlined and efficient and noisy about any errors, so we can issue the report to the oncologist both promptly and confidentally correctly. If it's very important to get something right, particular with time constraints the best thing to do is second check everything, and make our errors very loud and obvious.  


## About the dataset

rhampseq

**Questions:**
1.  Refering to the DPYD metabolism table, what score would be given to patient who has 1 normal funtion allele and 1 decreased function? Would you classify their phenotype as normal, intermediate or poor?
2. 
 
<details>
<summary>Answers</summary>
<ul><li>1. Line 67, as can be seen near the top of the traceback</li>
<li>2.  </li> </ul>
</details>



### Does the patient sample pass all sanity checks?

The throughput of NGS samples going through clinical laboratories around Australia is really high, and getting higher every year. If you have around 20,000 samples to process each year, and some lab steps are manual, how can we guarantee that no sample has been swapped or contaminated with another? 

One way is to separate the sample into two, right at the beginning when the lab first receives the sample. The first part of the sample goes through the normal testing process, and we generate data for it. The second part goes through a completely independent workflow, where we target just a handful of common SNPs.



There are quite a lot of other sanity checks we can do from the bioinformatics side. ancestry and pedigree



**Questions:**
1.  I haven't mentioned one really common sanity check to interrogate the sample integrity. Can you guess what it is? Hint: what's something obvious you can tell from a patient's genetic makeup? 

<details>
<summary>Answers</summary>
<ul><li>1. Doing a sex check to see that the genetic </li>
<li>2.  </li> </ul>
</details>



match vcfs with gatk fingerprint




**Questions:**
1.  What would happen to the LOD score if the sample swap occured *prior* to the lab receiving the sample? 

<details>
<summary>Answers</summary>
<ul><li>1. Line 67, as can be seen near the top of the traceback</li>
<li>2.  </li> </ul>
</details>


### Is there any contamination? 


High number of SNPs at low vaf 

#filter 

```bash
bcftools filter -i 'DP>10 && QUAL>30' file.vcf
```



## Is the variant correctly called?




Last week, 



No variant caller is perfect, and even the biggest and more robust programs can make mistakes. Take a look at this GATK Haplotype Caller issue page for incorrectly called variants:

use bcftools to check the variant call 

give link to gatk calling bugs (high amplicon reads)


## Running Quality Control and outputting pass or fail 









## What happens if things go wrong in the pipeline?


Incorrect quality threshold??


## 




## Generating the report 
 




#Generates summary statistics
bcftools stats file.vcf


GATK check fingerprint







match patient to report 



## Beware ChatGPT and AI friends







# Bonus tasks if time permitting 
3. Take a look at our bash script. "DPYD_mini_pipeline.sh". In that bash scripts, are the paths to the python and awk scripts absolute or relative? Could this cause issues? Can you change the path to be absolute instead? (hint: to get the path of a script or file, you could run):

```bash
 ls -d "$PWD"/{script_or_file_name} )
```



## Concluding remarks

Hopefully through doing this prac, you will find that even when the stakes for not failing are high, 