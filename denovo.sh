#!/usr/bin/env bash
# Used to aanalyze outcome of gene-editing outcome by de novo assembling with NGS data

baseP=$1
prefix=""
cd ${baseP}

bwa index Ref.txt
samtools faidx Ref.txt
tmpdir=denovo

mkdir ${tmpdir}

ls *R1.fq.gz | while read file
do
    file_tmp=${file##${prefix}}
    file_tmp=${file_tmp%.fq.gz}
    bwa mem -R "@RG\tID:1\tLB:library\tPL:Illumina\tSM:${file_tmp%.R1}\tPU:machine" Ref.txt $file > ${tmpdir}/${file_tmp}.bam
    samtools view -F 4 -S -u ${tmpdir}/${file_tmp}.bam > ${tmpdir}/${file_tmp}.accept_hits
    samtools sort -o ${tmpdir}/${file_tmp}.accept_hits.bam ${tmpdir}/${file_tmp}.accept_hits
    samtools index ${tmpdir}/${file_tmp}.accept_hits.bam
    # extract fastq from accept_hits.bam to fastq with bedtools utils
    bamToFastq \
        -i ${tmpdir}/${file_tmp}.accept_hits.bam \
        -fq ${tmpdir}/${file_tmp}.accept_hits.fastq
    # Trinity version 2.1.1
    Trinity \
        --genome_guided_bam ${tmpdir}/${file_tmp}.accept_hits.bam \
        --genome_guided_max_intron 10000  --max_memory 2G --CPU 4 \
        --no_version_check \
        --output ${tmpdir}/${file_tmp}.guide_trinity
    # Trinity genome-guide assembly mode
    align_and_estimate_abundance.pl \
        --transcripts ${tmpdir}/${file_tmp}.guide_trinity/Trinity-GG.fasta \
        --seqType fq \
        --single ${tmpdir}/${file_tmp}.accept_hits.fastq \
        --est_method RSEM \
        --aln_method bowtie \
        --trinity_mode  \
        --prep_reference \
        --output_dir ${tmpdir}/${file_tmp}.abound
done
