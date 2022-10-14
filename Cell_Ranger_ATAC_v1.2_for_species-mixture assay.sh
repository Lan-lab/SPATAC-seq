# Code for processing raw Fastq data from species-mixture assay
# SL1 file as an exampe
# SL1_R1.fastq.gz
# SL1_R2.fastq.gz
gunzip *fastq.gz
# extract Round 1 barcode file
fastx_trimmer -Q33 -f 39 -l 46 -i SL1_R1.fastq -o SL1_R1F.fastq
fastx_trimmer -Q33 -f 24 -l 31 -i SL1_R2.fastq -o SL1_R1R.fastq
# extract Round 2 and Round3 barcode file
fastx_trimmer -Q33 -f 1 -l 8 -i SL1_R1.fastq -o SL1_R2BC.fastq
fastx_trimmer -Q33 -f 1 -l 8 -i SL1_R2.fastq -o SL1_R3BC.fastq
# extract and reformat genome DNA Fastq files from raw Fastq file to Cell Ranger ATAC format
fastx_trimmer -Q33 -f 66 -l 115 -i SL1_R1.fastq -o SL1_S1_L001_R1_001.fastq
fastx_trimmer -Q33 -f 51 -l 100 -i SL1_R2.fastq -o SL1_S1_L001_R3_001.fastq
# Barcode 1bp correction
python3 ./ss2_py.py -1  key_value_round1_1.csv -2 SL1_R1F.fastq -o SL1_R1F_correction.fastq
python3 ./ss2_py.py -1  key_value_round1_2.csv -2 SL1_R1R.fastq -o SL1_R1R_correction.fastq
python3 ./ss2_py.py -1  key_value1.csv -2 SL1_R2BC.fastq -o SL1_R2BC_correction.fastq
python3 ./ss2_py.py -1  key_value2.csv -2 SL1_R3BC.fastq -o SL1_R3BC_correction.fastq
# merge and reformat four barcode Fastq files to one new barcode Fastq file (Cell Ranger ATAC format)
python ./BC_process.py -i SL1_R1F_correction.fastq -i2 SL1_R1R_correction.fastq -i3 SL1_R2BC_correction.fastq -i4 SL1_R3BC_correction.fastq -o SL1_S1_L001_R2_001.fastq
# deleting intermediate files
rm SL1_R1F.fastq
rm SL1_R1R.fastq
rm SL1_R2BC.fastq
rm SL1_R3BC.fastq
rm SL1_R1F_correction.fastq
rm SL1_R1R_correction.fastq
rm SL1_R2BC_correction.fastq
rm SL1_R3BC_correction.fastq
#
# Using cellranger-atac-1.2.0 to process the three reformated Fatsq files
# /fshare2/sky/software/10X_software/ATAC-seq/cellranger-atac-1.2.0/cellranger-atac count --id=test --reference=/fshare2/sky/software/10X_software/ATAC-seq/refdata-cellranger-atac-GRCh38-and-mm10-1.1.0 --fastqs=./ --sample=SL1
/your_software_path/cellranger-atac count --id=test --reference=/your_genome_file_path/refdata-cellranger-atac-GRCh38-and-mm10-1.1.0 --fastqs=./ --sample=SL1