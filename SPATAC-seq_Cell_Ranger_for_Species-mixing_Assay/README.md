## **Mixed-species proof-of-concept experiment**

### Demo **fastq** **data** 

SL1_R1.fastq.gz; SL1_R2.fastq.gz (10,000 paired end reads)

### Modification **in** **Cell** **Ranger** v1.2

Replace the cellranger-atac-cs/1.2.0/lib/python/barcodes/737K-cratac-v1.txt with the new barcodes file in the Modifications_in_Cell_Ranger_ATAC_v1.2;

Replace the cellranger-atac-cs/1.2.0/lib/python/barcodes/<u> <\u>init<u> <\u>.py with the new python file in the Modifications_in_Cell_Ranger_ATAC_v1.2.

### **Data processing**

#### 1. Raw Fastq data

Read 1: contains the Barcode1 (partial), Barcode2 and genome sequences

Read 2: contains the Barcode1 (partial), Barcode3 and genome sequences

#### 2. Reformat raw Fastq file to Cell Ranger ATAC format

Raw read 1 and read 2 > Barcodes file (BC_process.py; SL1_S1_L001_R2_001.fastq.gz)

Raw read 1 > New Read 1: contains the genome sequences (SL1_S1_L001_R1_001.fastq.gz)

Raw read 2 > New Read 2: contains the genome sequences (SL1_S1_L001_R3_001.fastq.gz)

![Image text](https://github.com/Lan-lab/SPATAC-seq/blob/main/SPATAC-seq_Cell_Ranger_for_Species-mixing_Assay/Fastq%20file%20in%2010X%20format.png)

#### 3. Sequence alignment and generation of fragments file by Cell Ranger

The reformated data was processed using Cell Ranger ATAC v1.2 with following references:

wget https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-atac-GRCh38-and-mm10-1.2.0.tar.gz



##### Note: The output files of Cell Ranger are used to evaluate the feasibility of SPATAC-seq.