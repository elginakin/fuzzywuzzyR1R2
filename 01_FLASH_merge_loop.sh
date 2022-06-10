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




