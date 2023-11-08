# Code for processing raw Fastq data from species-mixture assay
# Put all script files, barcode files, and Raw Fastq files into one same folder.
# Then run bash script to process Raw Fastq files (demultiplex, mapp and Create a fragment file)

# Here we used test_1.fastq.gz and test_2.fastq.gz file as an exampe (10,000 paired end reads)
# test_1.fastq.gz
# test_2.fastq.gz
gunzip test_1.fastq.gz
gunzip test_2.fastq.gz
## extract Round 1 barcode file
fastx_trimmer -Q33 -f 39 -l 46 -i test_1.fastq -o test_R1F.fastq
fastx_trimmer -Q33 -f 24 -l 31 -i test_2.fastq -o test_R1R.fastq
# extract Round 2 and Round3 barcode file
fastx_trimmer -Q33 -f 1 -l 8 -i test_1.fastq -o test_R2BC.fastq
fastx_trimmer -Q33 -f 1 -l 8 -i test_2.fastq -o test_R3BC.fastq
# extract and reformat genome DNA Fastq files from raw Fastq file to Cell Ranger ATAC format
fastx_trimmer -Q33 -f 66 -l 115 -i test_1.fastq -o test_F_1.fastq
fastx_trimmer -Q33 -f 51 -l 100 -i test_2.fastq -o test_F_2.fastq
# Barcode 1bp correction
python3 ss2_py.py -1  key_value1.csv -2 test_R2BC.fastq -o test_R2BC_correction.fastq
python3 ss2_py.py -1  key_value2.csv -2 test_R3BC.fastq -o test_R3BC_correction.fastq
python3 ss2_py.py -1  key_value_round1_1.csv -2 test_R1F.fastq -o test_R1F_correction.fastq
python3 ss2_py.py -1  key_value_round1_2.csv -2 test_R1R.fastq -o test_R1R_correction.fastq
# deleting intermediate files
rm test_R2BC.fastq
rm test_R3BC.fastq
rm test_R1F.fastq
rm test_R1R.fastq
# merge Round1 barcodes (16bp)
paste test_R1F_correction.fastq test_R1R_correction.fastq |sed 's/	//' >./test_R1.fastq
# deleting intermediate files
rm test_R1F_correction.fastq
rm test_R1R_correction.fastq
# merge Round2 and Round3 barcodes (16bp)
paste test_R2BC_correction.fastq test_R3BC_correction.fastq |sed 's/	//' >./test_R2R3BC_correction.fastq
# deleting intermediate files
rm test_R2BC_correction.fastq
rm test_R3BC_correction.fastq
# merge Round1, Round2 and Round3 barcodes (32bp)
paste test_R1.fastq test_R2R3BC_correction.fastq |sed 's/	//' >./test_R1R2R3BC_correction.fastq
# deleting intermediate files
rm test_R1.fastq
rm test_R2R3BC_correction.fastq
# Attach cell barcodes into the read name using sinto (https://timoast.github.io/sinto/scatac.html#create-a-fragment-file)
sinto barcode -b 32 --barcode_fastq test_R1R2R3BC_correction.fastq --read1 test_F_1.fastq --read2 test_F_2.fastq
# deleting intermediate files
rm test_R1R2R3BC_correction.fastq
rm test_F_1.fastq
rm test_F_2.fastq
# rename processed fastq files, which are used for mapping
mv test_F_1.barcoded.fastq test_1.fq
mv test_F_2.barcoded.fastq test_2.fq
# gzip *.fastq
gzip test_1.fastq
gzip test_2.fastq


# Map using bwa
for sample in *_1.fq
do
i=`echo $sample|cut -d '_' -f 1`
bwa mem -t 30 /your_reference_genome_path/genome.fa ${i}_1.fq ${i}_2.fq | samtools view -b - > ${i}.raw.bam
samtools sort -@ 16 ${i}.raw.bam -o ${i}.sort.bam
samtools index ${i}.sort.bam
rm ${i}.raw.bam
done
#
#
# gzip *.fq
gzip test_1.fq
gzip test_2.fq
#

# Create a fragment file using sinto
ls *.sort.bam |while read id;do
sinto fragments -b $id -p 30 -m 10 --use_chrom "[0-9]*" -f $id.fragments.bed  --barcode_regex "[^:]*"
grep -F -f barcode_file_32bp.csv $id.fragments.bed > $id.fragments.filter.bed
sort -k1,1 -k2,2n $id.fragments.filter.bed > $id.fragments.filter.sort.bed
bgzip $id.fragments.filter.sort.bed
tabix -p bed $id.fragments.filter.sort.bed.gz
rm $id.fragments.filter.bed
rm $id.sort.bam.fragments.bed
done


#
# Notes !!!
# Downstream analysis steps can be performed using Signac (https://stuartlab.org/signac/index.html);
# or ArchR (https://github.com/GreenleafLab/ArchR).