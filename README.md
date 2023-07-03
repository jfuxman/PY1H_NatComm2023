# PY1H reproducibility

## Predicting possible TF-TF interactions based on homology.
Possible TF-TF heterodimers were derived from LitBM. Check this [notebook](https://github.com/jfuxman/PY1H_reproducibility/blob/main/Get_TF_pairs.ipynb)

## Bioinformatics analysis of TF-prey sequencing data

### Input files
1. barcodes_fw.fasta and barcodes_rv.fasta: Barcode sequences to perform demultiplexing of Fastq files
2. Fastq Files are deposited at [link](www.google.com)
3. Het array with Ensembl Gene IDs.xlsx

### Scripts
1. trimming.sh: Demultiplexing and trimming of Fastq files
2. BowtieMapping.sh: Indexing and mapping of Fastq files using TF Transcripts as a library
3. get_enst_biomaRt.R: Get all transcript sequences for expected TFs in PY1H experiment
4. fromTxttoFasta.py: Converts to fasta file
5. BowtieMaping.sh: Index library, align FASTQ files, and generates bam files
6. Get_Significant_Pairs.R: Identify significant alignments and TF-TF pairs
 
## Obtaining ChIP-seq data from GTRD

### Input files
1. Final events list 6Feb23.xlsx: 437 TF-TF-cytokine interactions experimentally identified with PY1H (180 cooperative, 257 antagonism)

### Scripts
1. download_gtrd.py: Download MACS2 ChIP peaks from GTRD and metadata information.
2. script_CHIPSEQ.py: Get information on peaks between our TFs of interest and the respective cytokine promoters.


## Identification of binding sites of TF-pairs in cytokine promoters

### Input files
1. Final events list 6Feb23.xlsx: 437 TF-TF pairs experimentally identified with PY1H (180 cooperative, 257 antagonism)
2. promoters.fasta: Fasta file of cytokine promoters coded explained in the previous Excel file.

### Scripts
1. PWMs were downloaded from [CISBP 2.0 database](http://cisbp.ccbr.utoronto.ca/bulk.php) including `pwms_all_motifs` directory and `TF_information.txt` file.
2. TFPM_motifs_identification.R: Script to calculate significant TF motifs identified in cytokine promoters.
3. Promoter_Analysis_merged.Rmd: Script to merge overlapping TF motifs that were previously identified.

## Network Randomization Analysis (Promoters)

### Input files
1. Final events list 6Feb23.xlsx: 437 TF-TF pairs experimentally identified with PY1H (180 cooperative, 257 antagonism)
2. DF_all_motifs.xlsx: Table containing motifs information for each TF-TF-cytokine interaction.

### Scripts
Promoter_Randomization.Rmd: Script to perform randomization and reproduce Supplementary Figure 6 C, D, E 

## Network Randomization Analysis (ChIP peaks)

### Input files
1. Final events list 6Feb23.xlsx: 437 TF-TF pairs experimentally identified with PY1H (180 cooperative, 257 antagonism)
2. all_tfs_chipseq.txt: Table with TF name and Uniprot ID.
3. ChIP-seq.metadata.txt: Metadata file from GTRD database.
4. MACS2_Peaks_TFs_final_16_feb.xlsx:

### Scripts
CHIP_Randomization_antagonism.Rmd: Script to perform randomization and reproduce Figure 2G; and Supplementary Figure 6F, G

## Paralog partner similarity

### Input files

### Scripts
