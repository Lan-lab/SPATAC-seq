## Map the accessible chromatin landscape of zebrafish embryogenesis **at** **single** **cell** **resolution** **by** SPATAC-seq

Note: Put all script files, barcode files, and Raw Fastq files into one same folder. Then run Raw_fastq_demultiplexing_map_fragment.sh to process Raw Fastq files to create a fragement file.

### **Demo** **fastq** **data** 

test_1.fastq.gz; test_2.fastq.gz (10,000 paired end reads)

### Data processing

#### 1. Raw Fastq data

Read 1: contains the Barcode1 (partial), Barcode2 and genome sequences

Read 2: contains the Barcode1 (partial), Barcode3 and genome sequences

#### 2. Extract Barcodes file and genome sequences

Raw read 1 and read 2 > Barcodes file (32bp)

Raw read 1 > New Read 1: contains the genome sequences

Raw read 2 > New Read 2: contains the genome sequences

#### 3. Attach cell barcodes into the read names of New Read 1 and 2 by Sinto

![Image text](https://github.com/Lan-lab/SPATAC-seq/blob/main/Zebrafish_Embryogenesis_single-cell_oPen_chromatin_Atlas_(ZEPA)/fastq_preprocessing/barcodes%20in%20read%20name.png)



#### 4. Map by BWA


Zebrafish reference (built from the zebrafish GRCz11: GCA_000002035.4):

[Download DNA sequence (GRCz11; GCA_000002035.4)](https://ftp.ensembl.org/pub/release-110/fasta/danio_rerio/dna/)


#### 5. Create a fragment file by sinto

![Image text](https://github.com/Lan-lab/SPATAC-seq/blob/main/Zebrafish_Embryogenesis_single-cell_oPen_chromatin_Atlas_(ZEPA)/fastq_preprocessing/bed%20file.png)

The fragment files contain chromosome number (column 1), fragment start and end sites (column 2 and 3), cell barcodes (column 4) and duplicates number (column 5).

#### 6. Downstream analysis by ArchR, Signac or SnapATAC.

ArchR: https://github.com/GreenleafLab/ArchR

Signac: https://stuartlab.org/signac/index.html 

SnapATAC: https://github.com/r3fang/SnapATAC 
#




![Image text](https://github.com/Lan-lab/SPATAC-seq/blob/main/Zebrafish_Embryogenesis_single-cell_oPen_chromatin_Atlas_(ZEPA)/Atlas%20of%20zebrafish%20embryogenesis.png)