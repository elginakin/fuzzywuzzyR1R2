
# This Markdown file is for scratch testing of R1 R2 extraction on unmerged and merged fastq reads for the Bowman Congenital deep sequencing project. 


Using FLASH merge and seqit grep to extract reads from R1 and R2 Files based on sequence homology to reference amplicon sequence. 

TODO: 
1) Merge Reads by R1 and R2 using FLASH
2) Extract forward reads 
    - Extract based off F amplicon (10bp)
        - must be within first 20 bp of sequence
3) Extract reverse reads
    - Extract based off F amplicon (10bp )
    - Use quality cutoffs 
    - Use Length cutoffs (histogram of lengths)
4) Cat extracted merged reads back together.
5) Run QC stats on all generated fastq files

# Find based on F primer grep 100% homology 

Reference amplicon sequence:
 - Lenth:  516bp 
```
>XM_7971521:321-829_Trypanosoma_cruzi_strain_CL_Brener_lathosterol_oxidase_partial_mRNA
ATATTACAATGTTTCTGATTATGGATGGCCGTATCTTTTTCTTAGTATATTGATGTTTTTTATCTTTACG
GACTTTATGGTATATTGGTTTCATCGTGGTTTACATCACCCAACATTATACCGATACCTTCATAAATTAC
ATCATACATACAAATATACAACACCATTTTCATCTCATGCATTTAATCCTTGTGATGGATTTGGTCAAGG
TTCACCATATTACGCATTTATTTTTTTATTTCCTATGCATAATTACCTTTTTGTTATTCTCTTTTTCGCT
GTTAACCTATGGACAATCTCTATTCATGATCAAGTAGACTTTGGAGGTCATTTTGTGAACACAACAGGAC
ACCATACAATTCATCATGTACTTTTTAATTACGACTATGGACAGTACTTTACCGTATGGGACCGTATTGG
CGGAACTTATAAACCCGCACAACAGACGCATCTTTTCCCATTATTTACAAAAGGTGGACGCATTGAAGAG
GTTGAGTCAACGAAGAAGA
```

## F primer used
>XM_7971521321:321-829:1-10
ATATTACAAT

## R Primer used 
>XM_7971521321:321-829:478-488
GGTGGACGCA

## R primer used #reverse of AGTCAACGAA
>XM_7971521321:321-829:501-511
AAGCAACTGA

# Test Scripts 
## sekit extract test script 1 FAILED 

Failed - did not incoperate specific window to look for sequence homology

```sh

cat ./reads/035_0m_MID_1_R1.fastq | seqkit grep -s -i -p -m2 ATATTACAAT -o ./output_seqkit/035_0m_MID_1_R1_F-extract.fastq

seqkit stats ./output_seqkit/035_0m_MID_1_R1_F-extract.fastq ./reads/035_0m_MID_1_R1.fastq -o ./output_seqkit/test_stats.txt

#find based on grep 90% homology 
cat ./reads/035_0m_MID_1_R1.fastq | seqkit grep -P -m1 -s -i -p  ATATTACAAT -o ./output_seqkit/035_0m_MID_1_R1_F-extract-m1.fastq
#compare to original read 
seqkit stats ./output_seqkit/035_0m_MID_1_R1_F-extract-m2.fastq ./reads/035_0m_MID_1_R1.fastq -o ./output_seqkit/test_stats.txt

#stats all reads
seqkit stats ./output_seqkit/*.fastq ./reads/035_0m_MID_1_R1.fastq > ./output_seqkit/seqkit_stats.txt

#find REVERSE reads in R1 file based on grep 80% homology 

cat ./reads/035_0m_MID_1_R1.fastq | seqkit grep -P -m1 -s -i -p  ATATTACAAT -o ./output_seqkit/035_0m_MID_1_R1_F-extract-m1.fastq

#reverse complement
```
## sekit extract test script 2 FAILED 

Failed - did not incoperate specific window to look for sequence homology

```sh Merge and Extract single reads from 035_0m_MID_1_R1 and R2

cd /Users/elgin/3_Monica_and_Bob/Congenital_Chagas_Data_From_Jill/Bowman_congential_fastqs/Elgin/R1R2_Problem

#Merge with Flash with log 
    # Want to have it make a irectory and automatically name files within that directory

flash ./reads/035_0m_MID_1_R1.fastq ./reads/035_0m_MID_1_R2.fastq -o ./output_seqkit | tee 035_0m_MID_1_R1_flash.log | seqkit stats ./output_FLASH/*.fastq -o ./output_FLASH/test_stats.txt

#extract forward and reverse reads from merged file with seq kit (just curious)

cat ./output_FLASH/output_seqkit.extendedFrags.fastq | seqkit grep -P -m1 -s -i -p ATATTACAAT -o ./output_seqkit/FLASH_merged/035_0m_MID_1_R1R2_FLASH_F-extract-m1_.fastq

seqkit stats -a ./output_seqkit/FLASH_merged/035_0m_MID_1_R1R2_FLASH_F-extract-m1_.fastq > ./output_seqkit/FLASH_merged/035_0m_MID_1_R1R2_FLASH_F-extract-m1_stats.txt

```

## seqit extract test script 3 CORRECT

TODO: 

Post flash merge, must re-extract reverse reads using primer homology. Not sure what primers to use for sequence homology approach this. Options: 

- Reverse complement of amplicon primer (~10bp)
- Reverse of amplicon primer 

```sh

cd /Users/elgin/Library/CloudStorage/OneDrive-JohnsHopkins/Rotation_Folders/3_Monica_and_Bob/Congenital_Chagas_Data_From_Jill/Bowman_congential_fastqs/Elgin/Sandbox_attepmpts/01_R1R2_Problem_Sandbox

#extract forward reads from primers within first 20bp
cat ./merged_reads/035_0m_MID_1.extendedFrags.fastq | seqkit grep -s -R 1:20 -i -p -m 1 ATATTACAAT -o ./output_seqkit_R-extract_test/035_0m_MID_1.extendedFrags_F-extract.fastq

#extract reverse reads from primers within first 20bp
cat ./merged_reads/035_0m_MID_1.extendedFrags.fastq | seqkit grep -s -R 1:20 -i -p -m 1 AAGCAACTGA -o ./output_seqkit_R-extract_test/035_0m_MID_1.extendedFrags_R-extract.fastq #primer is reverse of nt 501-511 of nucleotides

TODO:
#for local run testing in wd THIS WORKED 

#Forward Extract
cat 035_0m_MID_1.extendedFrags.fastq | seqkit grep -s -R1:70 -m1 -i -P -p ATATTACAAT 035_0m_MID_1.extendedFrags.fastq -o F_extract.fastq 

#reverse Extract
cat 035_0m_MID_1.extendedFrags.fastq | seqkit grep -s -R1:70 -m1 -i -P -p TCCCATCTTCGTTGACT 035_0m_MID_1.extendedFrags.fastq -o R_extract.fastq

seqkit stats -a *.fastq > F_and_R_extract_stats.txt

Try primers 
AAGCAACTGA #use this 
AATGCGTCCACCT #rev complement of some end of amplicon bp 
TCCCATCTTCGTTGACT #use this #from Jill - Reverse complement of end of read

Flags 

-s 
-p #ignore case 
-P #read only positive strand
-R 1:20 
-i 
-m1 

035_0m_MID_1.extendedFrags.fastq 


-o 035_0m_MID_1.extendedFrags_R-extract.fastq #primer is reverse of nt 501-511 of nucleotides


```
# REDO with FLASH merged reads

Instead of extracting fastq files individually, merge reads to solve directionality issues and extract. 

```sh FLASH merge
#loop for FLASH merging forward and reverse reads for fastq files
#make a variable for date and generate a FLASH 

mkdir correct_orig

for sample in `ls ./*R1.fastq`
do
echo = " merging ${base}_R1.fastq ${base}_R2.fastq" 

dir="/Users/elgin/Library/CloudStorage/OneDrive-JohnsHopkins/Rotation_Folders/3_Monica_and_Bob/Congenital_Chagas_Data_From_Jill/Bowman_congential_fastqs/correct_orig"

base=$(basename $sample "_R1.fastq")

#run flash with output base
flash -o ${base} -d /Users/elgin/Library/CloudStorage/OneDrive-JohnsHopkins/Rotation_Folders/3_Monica_and_Bob/Congenital_Chagas_Data_From_Jill/Bowman_congential_fastqs/correct_orig /${dir}/${base}_R1.fastq /${dir}/${base}_R2.fastq

done
```

```sh Run seqit summary for fastq files

seqkit stats -a *.fastq > Merged_reads_only.stats.txt &

```

# Automate extracting sequences from merged reads
```sh rewrite single extraction event with seqkit

#extract forward read with 90% homology to first 10bp of amplicon
cat ./output_FLASH/output_seqkit.extendedFrags.fastq | seqkit grep -P -m1 -s -i -p ATATTACAAT -o ./output_seqkit/FLASH_merged/035_0m_MID_1_R1R2_FLASH_F-extract-m1_.fastq

```

TODO: 
```sh primitive bash for loop for seqkit extraction 

for sample in `ls ./*.fastq`
do

echo = "extracting forward reads from ${base}.fastq" 

dir="/Users/elgin/Library/CloudStorage/OneDrive-JohnsHopkins/Rotation_Folders/3_Monica_and_Bob/Congenital_Chagas_Data_From_Jill/Bowman_congential_fastqs/Elgin/Sandbox_attepmpts/EXTRACT_loop_test/reads"

outdir="/Users/elgin/Library/CloudStorage/OneDrive-JohnsHopkins/Rotation_Folders/3_Monica_and_Bob/Congenital_Chagas_Data_From_Jill/Bowman_congential_fastqs/Elgin/Sandbox_attepmpts/EXTRACT_loop_test/extract_output_F"

base=$(basename $sample ".fastq") 

#run flash with output base

cat ${base}.fastq | seqkit grep -P -m1 -s -i -p ATATTACAAT -o /${outdir}/${base}_F_extract.fastq

seqkit -o ${base} -d /Users/elgin/Library/CloudStorage/OneDrive-JohnsHopkins/Rotation_Folders/3_Monica_and_Bob/Congenital_Chagas_Data_From_Jill/Bowman_congential_fastqs/correct_orig /${dir}/${base}_R1.fastq /${dir}/${base}_R2.fastq

#run sh /Users/elgin/3_Monica_and_Bob/scripts/R1R2_parse_scripts/seqkit_extract_forward_reads.sh

```

TODO:
```py script to automate extracting true F and R reads from merged reads

"""adopting from Jaclyn's sort reads by primer cutadapt functions"""

import os
import glob
import pandas as pd

def flash_merge_reads(read_1, read_2)

    """ 
    Function to automate FLASH merging of R1 and R2 reads on all sequences within a directory. 

    Locates R1 and R2 files with the extension suffix "_R1 or _R2" and merges them with a name extension of _merge.fastq

    Args 
        read_1: first .fastq file with R1 file extension 
        read_2: second .fastq file with R2 file extension
    """
    # initialize flags/variables
    i = True

    #universal FLASH merge flags
    cut_flags = #for flash
    input_files = seq_file_1 + " " + seq_file_2
    
    #create folder for FLASH output files
    folder_name = "FLASH"
    try: 
        os.mkdir(folder_name)
        print("folder '{}' created ".format(folder_name))
    except FileExistsError:
        print("folder {} already exists".format(folder_name))
    
    #define output files 

    output_merged_reads = "-o" + folder + "_merged.fastq"

    #run FLASH
    #merging reads
    if i: 
        flash_command = "flash" + cutflags + input_files 

        os.system(flash_command + " > /flash_output_" + file_name + ".fasta" 

        #end of run 1 
        i = FALSE

```

# Post Merge, extraction, concatenation, and stats: 
## fastqc on all extracted and merged reads

```sh fastqc then multiqc

fastqc *.fastq.gz -t 8 -o /Users/elgin/Library/CloudStorage/OneDrive-JohnsHopkins/Rotation_Folders/3_Monica_and_Bob/Congenital_Chagas_Data_From_Jill/Bowman_congential_fastqs/Elgin/correct_FLASH-merge/Merged_reads_only/Primer_Extracted_reads/concat_extracted_fastqs/fastqc &

#run multiqc to compile all fastq reports
multiqc .

```

# Final Scripts - Split into 2 shell scripts for execution 

# 01 - FLASH merge of all corrected reads

```sh 01 FLASH MERGE 
#loop for FLASH merging forward and reverse reads for fastq files

#generate directory for flash_merge
mkdir correct_orig

#loop through all fastqs in designated directory and output them to PATH.  
for sample in `ls ./*R1.fastq`
do

echo = "merging ${base}_R1.fastq ${base}_R2.fastq with flash (default parameters)" 

#must change path

#relative path
#dir = "/.correct"

dir="/Users/elgin/Library/CloudStorage/OneDrive-JohnsHopkins/Rotation_Folders/3_Monica_and_Bob/Congenital_Chagas_Data_From_Jill/Bowman_congential_fastqs/correct_orig"
outdir="./FLASH_output"
base=$(basename $sample "_R1.fastq")

#run flash with default parameter with output base
flash -o ${base} -d /${outdir} /${dir}/${base}_R1.fastq /${dir}/${base}_R2.fastq

done
```

# 02 - extraction of all reads based on primer homology

```sh 02 seqkit Extract

# This script is designed to extract reads from one MERGED fastq file (merged via FLASH) 
# based on primer homology of ~90+% identity, reverse complement the reverse reads, and concatenate
# them back into a single file. Final products are merged reads all in positive sense orientation. 

#move script to dir with all fastqs OR change to relative paths
#script extracts forward and reverse reads based on primer homology - seqkit sliding search

#generate directories
mkdir 02_intermediate_fastqs
mkdir 02_concat_extracted_fastqs
mkdir 02_MERGED_F-R-comp_extracted_fastqs
mkdir logs

#list through FLASH output fastqs
for sample in *.fastq
do

base=$(basename $sample ".fastq") 

echo "extracting FORWARD and REVERSE reads from $sample"

# use seqkit to extract sequences based on forward primer sequence.
# Flags indicate scanning of first 70bp on ONLY the positive strand for primer homology (max mismatch of 1) and outputs file with an appropriate tail name

# extract forward reads
seqkit grep -s -R1:70 -m1 -i -P -p ATATTACAAT ${base}.fastq -o ${base}_F_extracted.fastq

# extract reverse reads
seqkit grep -s -R1:70 -m1 -i -P -p TCCCATCTTCGTTGACT ${base}.fastq -o ${base}_R_extracted.fastq

# reverse complement all reverse extracted reads 
seqkit seq -r -p -t dna -v ${base}_R_extracted.fastq -o ${base}_R-comp_extracted.fastq

# cat reads
cat ${base}_F_extracted.fastq ${base}_R-comp_extracted.fastq > ${base}MERGED_F-R-comp_extracted_cat.fastq

# zip fastq files 
gzip ${base}_F_extracted.fastq  
gzip ${base}_R_extracted.fastq
gzip ${base}_R-comp_extracted.fastq
gzip ${base}MERGED_F-R-comp_extracted_cat.fastq

done

echo ""
echo "Summarizing generated fastq Files" #also includes intermediate sample files

seqkit stats -a *.fastq.gz > 02_COMPLETE_extraction_stats.txt #summary file of all fastqs

echo "moving to designated directories"
#move cat files 
mv *MERGED_F-R-comp_extracted_cat.fastq.gz 02_MERGED_F-R-comp_extracted_fastqs

#move intermediate files
mv *F_extracted.fastq.gz *R_extracted.fastq.gz *R-comp_extracted.fastq.gz 02_intermediate_fastqs 
mv *.txt logs


```

# Extraction to seperate R1 R2 files

```sh Extraction to seperate R2 R2 files

#starting from individual files output to 2 files for DADA2 or SeekDeep Analysis

mkdir 02_orig_correct_reads
mkdir 02_intermediate_fastqs
mkdir 02_concat_extracted_fastqs-FandR
mkdir logs

for sample in `ls ./*_R1.fastq`
do

base=$(basename $sample ".fastq") 

echo = "extracting forward and reverse reads from ${base}_R1.fastq ${base}_R2.fastq with seqkit" 

# use seqkit to extract sequences based on forward primer sequence.
# Flags indicate scanning of first 70bp on ONLY the positive strand for primer homology (max mismatch of 1) and outputs file with an appropriate tail name

# Extract from R1 reads
# extract R1 forward reads scanning first 70bp in reads
seqkit grep -s -R1:50 -m1 -i -P -p ATATTACAAT ${base}_R1.fastq -o ${base}_R1-F_extracted.fastq

# extract R1 reverse reads
seqkit grep -s -R1:50 -m1 -i -P -p TCCCATCTTCGTTGACT ${base}_R1.fastq -o ${base}_R1-R_extracted.fastq

# reverse complement all R1 reverse extracted reads 
seqkit seq -r -p -t dna -v ${base}_R1-R_extracted.fastq -o ${base}_R1-R-comp_extracted.fastq




# Extract from R2 reads
# extract R2 forward reads
seqkit grep -s -R1:70 -m1 -i -P -p ATATTACAAT ${base}_R2.fastq -o ${base}_R2-F_extracted.fastq

# extract R2 reverse reads
seqkit grep -s -R1:70 -m1 -i -P -p TCCCATCTTCGTTGACT ${base}_R2.fastq -o ${base}_R2-R_extracted.fastq

# reverse complement all R1 reverse extracted reads 
seqkit seq -r -p -t dna -v ${base}_R2-R_extracted.fastq -o ${base}_R2-R-comp_extracted.fastq

# cat forward reads
cat ${base}_R1-F_extracted.fastq ${base}_R2-F_extracted.fastq > ${base}R1R2-F_extracted.fastq

# cat reverse complement of 

# zip fastq files 
gzip ${base}_F_extracted.fastq  
gzip ${base}_R_extracted.fastq
gzip ${base}_R-comp_extracted.fastq
gzip ${base}MERGED_F-R-comp_extracted_cat.fastq

done

echo ""
echo "Summarizing generated fastq Files" #also includes intermediate sample files

seqkit stats -a *.fastq.gz > 02_COMPLETE_extraction_stats.txt #summary file of all fastqs

echo "moving to designated directories"
#move cat files 
mv *MERGED_F-R-comp_extracted_cat.fastq.gz 02_MERGED_F-R-comp_extracted_fastqs

#move intermediate files
mv *F_extracted.fastq.gz *R_extracted.fastq.gz *R-comp_extracted.fastq.gz 02_intermediate_fastqs 
mv *.txt logs


```


```sh Test retrieval number of reads using seqkit
#rev comp primer from Jill: TCCCATCTTCGTTGACT - taken from MID spreadsheet

# test on a forward read
# test on 4296_1_1m_MID_1_R1.fastq

#extract forward reads
seqkit grep -s -R1:50 -m1 -i -P -p ATATTACAAT 4296_1_1m_MID_1_R1.fastq -o F-extracted_reads.fastq

#using reverse complement 
seqkit grep -s -R1:50 -m1 -i -P -p TCCCATCTTCGTTGACT 4296_1_1m_MID_1_R1.fastq -o rev_rev-comp-primer_extracted_reads_jills.fastq

#using just reverse of last 

seqkit grep -s -R1:50 -m1 -i -P -p AGAAGAAGCAACTGAGTTG 4296_1_1m_MID_1_R1.fastq -o rev_rev-primer_extracted_reads.fastq

#using (+) sense last base pairs of end of amplicon
seqkit grep -s -R1:50 -m1 -i -P -p GTTGAGTCAACGAAGAAGA 4296_1_1m_MID_1_R1.fastq -o rev_pos-sense-endofamplicon_extracted_reads.fastq

#using Jill's rev comp primer TCCCATCTTCGTTGACT

#try rev comp of end of primer: TCTTCTTCGTTGACTCAC 
seqkit grep -s -R1:50 -m1 -i -P -p TCTTCTTCGTTGACTCAC 4296_1_1m_MID_1_R1.fastq -o rev_rev-comp-primer_extracted_reads.fastq

# summary statistics
seqkit stats *.fastq

```
# Result: 
file                                               format  type  num_seqs    sum_len  min_len  avg_len  max_len   Q1     Q2   Q3  sum_gap  N50  Q20(%)  Q30(%)
4296_1_1m_MID_1_R1.fastq                           FASTQ   DNA      9,312  2,508,404      100    269.4      300  252    300  300        0  300   89.64   81.65
F-extracted_reads.fastq                            FASTQ   DNA      1,890    566,646      175    299.8      300  300    300  300        0  300    97.9   94.34
rev_pos-sense-endofamplicon_extracted_reads.fastq  FASTQ   DNA          2        473      173    236.5      300  173  236.5  300        0  300   67.65    53.7
rev_rev-comp-primer_extracted_reads.fastq          FASTQ   DNA      5,495  1,490,106      100    271.2      300  257    300  300        0  300   87.84   78.93
rev_rev-primer_extracted_reads.fastq                                    0          0        0        0        0    0      0    0        0    0       0       0

Proper primer sequences to use 
Scanning first 20 abs pairs results in zero extracted reads...what do the first 20bp look like? are MIDs, adapters, etc still there? 
    - Lowest number of bp possible for scanning of region of homology? 
    - 50 grabs everything 
 
##     Forward Primer: ATATTACAAT (first 10 base pairs, scanning first 70 bp)  - move down to 50 

## Reverse primer 
    - uses reverse complement scanning first 50 bp of read for primer homology
    - 

    


