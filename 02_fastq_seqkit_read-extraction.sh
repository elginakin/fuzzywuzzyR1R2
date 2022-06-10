#move script to dir with all fastqs OR change to relative paths
#script extracts forward and reverse reads based on primer homology 

echo "                     ¯\_(ツ)_/¯ we love to see it  "
echo ""
echo "    ___ __"
echo "   (_  ( . ) )__                  '.    \   :   /    .'"
echo "     '(___(_____)      __           '.   \  :  /   .'"
echo "                     /. _\            '.  \ : /  .'"
echo "                .--.|/_/__      -----____   _  _____-----"
echo "_______________''.--o/___  \_______________(_)___________"
echo "       ~        /.'o|_o  '.|  ~                   ~   ~"
echo "  ~            |/    |_|  ~'         ~"
echo "               '  ~  |_|        ~       ~     ~     ~"
echo "      ~    ~          |_|O  ~                       ~"
echo "             ~     ___|_||_____     ~       ~    ~"
echo "   ~    ~      .'':. .|_|A:. ..::''."
echo "             /:.  .:::|_|.\ .:.  :.:\   ~"
echo "  ~         :..:. .:. .::..:  .:  ..:.       ~   ~    ~"
echo "             \.: .:  :. .: ..::    /"
echo "    ~      ~      ~    ~    ~         ~"
echo "               ~           ~    ~   ~  Im an island boi      ~"
echo "        ~         ~            ~   ~                 ~"
echo "   ~                  ~    ~ ~                 ~"
echo "   "
echo ""
echo ""
echo "This script extracts merged forward and reads from a single fastq into 2 seperate files and merges them back into a single zipped gzipped fastq. Sumulative summary statistics are produced for all generated fastqs including intermediate extraction files for the forward and reverse reads"
echo ""
echo "forward primer used: ATATTACAAT"
echo "reverse primer used: TCCCATCTTCGTTGACT"
echo ""

#forward primer extraction

mkdir intermediate_fastqs
mkdir concat_extracted_fastqs

#cycle through FLASH output fastqs
for sample in *.extendedFrags.fastq
do
base=$(basename $sample ".fastq") 
echo "extracting FORWARD and REVERSE reads from $sample"

#use seqkit to extract sequences based on forward primer sequence.
#Flags indicate scanning of first 70bp on ONLY the positive strand for primer homology (max mismatch of 1) and outputs file with an appropriate tail name

#extract forward reads
seqkit grep -s -R1:70 -m1 -i -P -p ATATTACAAT ${base}.fastq -o ${base}_F_extracted.fastq

#extract reverse reads
seqkit grep -s -R1:70 -m1 -i -P -p TCCCATCTTCGTTGACT ${base}.fastq -o ${base}_R_extracted.fastq

#reverse complement all reverse extracted reads 
seqkit seq -r -p ${base}_R_extracted.fastq -o ${base}_R-comp_extracted.fastq

#cat reads
cat ${base}_F_extracted.fastq ${base}_R_extracted.fastq > ${base}_F-R_extracted_cat.fastq

#zip fastq files 
gzip ${base}_F_extracted.fastq  
gzip ${base}_R_extracted.fastq
gzip ${base}_F-R_extracted_cat.fastq

done

echo ""
echo "Summarizing generated fastq Files" #also includes intermediate sample files

seqkit stats -a *.fastq.gz > COMPLETE_extraction_stats.txt

#gzip *.fastq #zip all files post sumamry
mv *extracted.fastq.gz intermediate_fastqs # move zipped files to designated dir 
mv *_extracted_cat.fastq.gz concat_extracted_fastqs # move concatenated extraction files to designated directory


