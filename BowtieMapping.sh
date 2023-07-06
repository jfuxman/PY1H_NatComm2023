#!/bin/bash -l

#$ -l h_rt=12:00:00
#$ -N BowtieMapping3
#$ -pe omp 16
#$ -o BowtieMapping3.log
#$ -e BowtieMapping3.err
module load bwa
module load samtools

# Create index
#mkdir TF_index
bwa index -p TF_index_new/TF_transcripts_db Transcript_Seq_new.fasta

mkdir bwa_outputs
for i in {1..8}
do
	for k in {A..H}{01..12}
	do
		bwa mem -t 12 -M TF_index_new/TF_transcripts_db ${i}_clean_files_relaxed_prev/${i}_ROI_${k}_PGADT7_F_trimmed.fastq ${i}_clean_files_relaxed_prev/${i}_ROI_${k}_PGADT7_R_trimmed.fastq > bwa_outputs/${i}_${k}_PGADT7.sam
		samtools view -@ 12 -S -b bwa_outputs/${i}_${k}_PGADT7.sam > bwa_outputs/${i}_${k}_PGADT7.bam
		samtools sort -@ 12 bwa_outputs/${i}_${k}_PGADT7.bam -o bwa_outputs/${i}_${k}_PGADT7.sorted.bam
		samtools index -@ 12 bwa_outputs/${i}_${k}_PGADT7.sorted.bam bwa_outputs/${i}_${k}_PGADT7.sorted.bai 
		echo ${i}_clean_files_relaxed_prev/${i}_ROI_${k}_PGADT7_F_trimmed.fastq
		echo ${i}_clean_files_relaxed_prev/${i}_ROI_${k}_PGADT7_R_trimmed.fastq

		bwa mem -t 12 -M TF_index_new/TF_transcripts_db ${i}_clean_files_relaxed_prev/${i}_ROI_${k}_AD2U_F_trimmed.fastq ${i}_clean_files_relaxed_prev/${i}_ROI_${k}_AD2U_R_trimmed.fastq > bwa_outputs/${i}_${k}_AD2U.sam
		samtools view -@ 12 -S -b bwa_outputs/${i}_${k}_AD2U.sam > bwa_outputs/${i}_${k}_AD2U.bam
		samtools sort -@ 12 bwa_outputs/${i}_${k}_AD2U.bam -o bwa_outputs/${i}_${k}_AD2U.sorted.bam
		samtools index -@ 12 bwa_outputs/${i}_${k}_AD2U.sorted.bam bwa_outputs/${i}_${k}_AD2U.sorted.bai 
		echo ${i}_clean_files_relaxed_prev/${i}_ROI_${k}_AD2U_F_trimmed.fastq
		echo ${i}_clean_files_relaxed_prev/${i}_ROI_${k}_AD2U_R_trimmed.fastq


	done

done
