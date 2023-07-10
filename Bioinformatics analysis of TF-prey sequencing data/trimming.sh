#!/bin/bash -l

#$ -l h_rt=12:00:00
#$ -l mem_total=125G
#$ -pe omp 16
#$ -N cuttingp_new
#$ -o cuttingp_new.log
#$ -e cuttingp_new.err

module load python3
module load cutadapt

score=118
for i in 8
do
	other=`expr ${score} + ${i}`
	## Without trimming adapters

	## Trimming 'barcoding' primers in HIS
	## We can increase error to 0.2
	cutadapt -e 0.1 \
	--pair-filter=both \
	--action=trim \
	-j 12 \
	-O 10 \
	-g file:barcodes_fw.fasta \
	-A file:barcodes_rv.fasta \
	-o ${i}_{name}_F_trimmed.fastq \
	-p ${i}_{name}_R_trimmed.fastq \
	Jan_data/20221230_JFB${i}_JFB10613_S${other}_R1_001.fastq.gz Jan_data/20221230_JFB${i}_JFB10613_S${other}_R2_001.fastq.gz

	## Trimming again 'Primer seed sequences'

	for k in {A..H}{01..12}
	do
		## Trimming 'AD2U Primer'
		cutadapt -e 0.2 \
		--action=trim \
		--pair-filter=both \
		-j 16 \
		-O 10 \
		-g GGTGGGTCGAATCAA \
		-A TTGATTCGACCCACC \
		-o ${i}_ROI_${k}_AD2U_F_trimmed.fastq \
		-p ${i}_ROI_${k}_AD2U_R_trimmed.fastq \
		--untrimmed-output ${i}_ROI_${k}_AD2U_F_untrimmed.fastq \
		--untrimmed-paired-output ${i}_ROI_${k}_AD2U_R_untrimmed.fastq \
		${i}_${k}_F_trimmed.fastq ${i}_${k}_R_trimmed.fastq

		## Trimming 'PGADT7 Primer'
		cutadapt -e 0.2 \
		--action=trim \
		--pair-filter=both \
		-j 16 \
		-O 20 \
		-g ATCTTTAATACGACTCACTATAGGGCG \
		-A CATATGAGCGTAATCTGGTACGTCGTA \
		-o ${i}_ROI_${k}_PGADT7_F_trimmed.fastq \
		-p ${i}_ROI_${k}_PGADT7_R_trimmed.fastq \
		--untrimmed-output ${i}_ROI_${k}_PGADT7_F_untrimmed.fastq \
		--untrimmed-paired-output ${i}_ROI_${k}_PGADT7_R_untrimmed.fastq \
		${i}_ROI_${k}_AD2U_F_untrimmed.fastq ${i}_ROI_${k}_AD2U_R_untrimmed.fastq
	done

	mkdir ${i}_clean_files_relaxed_prev
	mv ${i}_ROI*_trimmed.fastq ${i}_clean_files_relaxed_prev

	rm *trimmed.fastq
done
