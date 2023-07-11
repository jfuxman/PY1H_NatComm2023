# PY1H reproducibility

## Predicting possible TF-TF interactions based on homology.
Possible TF-TF heterodimers were derived from LitBM. Check this [notebook](https://github.com/jfuxman/PY1H_reproducibility/blob/main/Get_TF_pairs.ipynb)

## Bioinformatics analysis of TF-prey sequencing data
### Input Files
1. barcodes_fw.fasta and barcodes_rv.fasta: Barcode sequences to perform demultiplexing of Fastq files
2. Fastq Files are deposited at [link](www.google.com)
3. Het array with Ensembl Gene IDs.xlsx
### Scripts
1. trimming.sh: Demultiplexing and trimming of Fastq files.
2. get_enst_biomaRt.R: Get all transcript sequences for expected TFs in PY1H experiment.
3. fromTxttoFasta.py: Converts to fasta file.
4. BowtieMaping.sh: Index library, align FASTQ files, and generates bam files.
5. Get_Significant_Pairs.R: Identify significant alignments and TF-TF pairs.
### Output Files
1. PGADT7_AD2U_results.txt

## Manual revision of detected TF-TF-cytokine interactions
### Input Files
1. PGADT7_AD2U_results.txt
### Output Files
1. Final events list 6Feb23.xlsx: 437 TF-TF-cytokine interactions experimentally identified with PY1H (180 cooperative, 257 antagonism).

## Obtaining ChIP-seq data from GTRD
### Input Files
1. Final events list 6Feb23.xlsx: 437 TF-TF-cytokine interactions experimentally identified with PY1H (180 cooperative, 257 antagonism).
### Scripts
1. download_gtrd.py: Download MACS2 ChIP peaks from GTRD and metadata information.
2. script_CHIPSEQ.py: Get information on peaks between our TFs of interest and the respective cytokine promoters.
### Output Files
1. MACS2_Peaks_TFs_final_16_feb.xlsx

## Identification of binding sites of TF-pairs in cytokine promoters
### Input Files
1. Final events list 6Feb23.xlsx: 437 TF-TF pairs experimentally identified with PY1H (180 cooperative, 257 antagonism)
2. promoters.fasta: Fasta file of cytokine promoters coded explained in the previous Excel file.
### Scripts
1. PWMs were downloaded from [CISBP 2.0 database](http://cisbp.ccbr.utoronto.ca/bulk.php) including `pwms_all_motifs` directory and `TF_information.txt` file.
2. TFPM_motifs_identification.R: Script to calculate significant TF motifs identified in cytokine promoters.
3. Promoter_Analysis_merged.Rmd: Script to merge overlapping TF motifs that were previously identified.
### Output Files
1. motifs_results_merged_motifs.xlsx
2. DF_all_motifs.xlsx

## Network Randomization Analysis (Promoters)
### Input Files
1. Final events list 6Feb23.xlsx: 437 TF-TF pairs experimentally identified with PY1H (180 cooperative, 257 antagonism)
2. DF_all_motifs.xlsx: Table containing motifs information for each TF-TF-cytokine interaction.
### Scripts
Promoter_Randomization.Rmd: Script to perform randomization and reproduce Supplementary Figure 6 C, D, E 
### Output Files
Supplementary Figure 6C, 6D, 6E

## Network Randomization Analysis (ChIP peaks)
### Input files
1. Final events list 6Feb23.xlsx: 437 TF-TF pairs experimentally identified with PY1H (180 cooperative, 257 antagonism)
2. all_tfs_chipseq.txt: Table with TF name and Uniprot ID.
3. ChIP-seq.metadata.txt: Metadata file from GTRD database.
4. MACS2_Peaks_TFs_final_16_feb.xlsx: Peaks information for TF in pY1H
### Scripts
CHIP_Randomization_antagonism.Rmd: Script to perform randomization and reproduce Figure 2G; and Supplementary Figure 6F, G
### Output Files
Figure 2G; Supplementary Figure 6F, 6G

## Network Randomization Analysis (EY1H and PY1H with ChIP data)
### Input files
1. py1h ey1h chip-seq overlap list for Luis.xlsx: 270 TF-promoter from eY1H and 256 TF-promoter from pY1H
2. all_tfs_chipseq.txt: Table with TF name and Uniprot ID.
3. MACS2_Peaks_TFs_final_26_june.xlsx:  Peaks information for TF in eY1H
4. MACS2_Peaks_TFs_final_16_feb.xlsx: Peaks information for TF in pY1H
### Scripts
Reviewers_randomization.Rmd
### Output Files
Figure X.

## Paralog partner similarity

### Input files
1. proteins.fasta: Protein sequences for TFs with Uniprot ID obtained from UNIPROT
2. uniprot.tsv: Metadata table from Uniprot
3. Array TF families.xlsx: TFs from pY1H and TF Family
4. 297 pairs coop+ant events.xlsx:
### Scripts
1.Jaccard_Analysis.Rmd
### Output Files
1. Paralog_analysis_output.xlsx: Jaccard score and alignment score for each pair of TF
2. Figure 3D
  
