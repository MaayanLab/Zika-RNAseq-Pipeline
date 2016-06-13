#!/bin/bash
WORKDIR='/home/maayanlab/Zika'
SRA_BIN='/home/maayanlab/Downloads/sratoolkit.2.5.7-ubuntu64/bin/'
cmd_name='fastq-dump'
featureCounts='/home/maayanlab/Downloads/subread-1.4.6-p2-Linux-x86_64/bin/featureCounts'

GENOME_GTF='/home/maayanlab/Zika/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf'
GENOME_FA='/home/maayanlab/Zika/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa'
STAR_INDEX='/home/maayanlab/Zika/Homo_sapiens/UCSC/hg19/star/STAR_2.4.1c/'

cd $WORKDIR

## create dirs
mkdir fastqs
mkdir paired_fastqs
mkdir fastQC_output
mkdir star_output
mkdir featureCount_output

## dump .sra to .fastq
for sra in $(ls SINGLE/SRR*/*.sra); do
	echo $sra
	$SRA_BIN$cmd_name -O fastqs $sra
done

## Note that paired-end sequencing reads should be dumped with different params
for sra in $(ls PAIRED/SRR*/*.sra); do
	echo $sra
	$SRA_BIN$cmd_name -I --split-files -O paired_fastqs $sra
done

## make star index
# STAR \
#     --runThreadN 8 \
#     --runMode genomeGenerate \
#     --genomeDir $STAR_INDEX \
#     --genomeFastaFiles $GENOME_FA \
#     --sjdbGTFfile $GENOME_GTF \
#     --sjdbOverhang 100

## align and assemble single-end sequencing reads
cd fastqs
for fq in $(ls); do
	basename=$(echo $fq | cut -f1 -d '.')
	echo $basename

	fastqc $fq -o fastQC_output

	STAR \
		--genomeDir $STAR_INDEX \
		--sjdbGTFfile $GENOME_GTF \
		--runThreadN 8 \
		--outSAMstrandField intronMotif \
		--outFilterIntronMotifs RemoveNoncanonical \
		--outFileNamePrefix $WORKDIR/star_output/$basename \
		--readFilesIn $fq \
		--outSAMtype BAM SortedByCoordinate \
		--outReadsUnmapped Fastx \
		--outSAMmode Full

	suffix="Aligned.sortedByCoord.out.bam"
	outname="$basename.count.txt"
	bam="$WORKDIR/star_output/$basename$suffix"
	$featureCounts \
		-T 8 \
		-t exon \
		-g gene_id \
		-a $GENOME_GTF \
		-o $WORKDIR/featureCount_output/$outname \
		$bam
done

## align and assemble paired-end sequencing reads
cd paired_fastqs
for basename in $(ls | cut -f1 -d '_' | sort | uniq); do
	echo $basename
	fq1="_1.fastq"
	fq2="_2.fastq"
	fq1=$basename$fq1
	fq2=$basename$fq2

	fastqc $fq1 -o fastQC_output
	fastqc $fq2 -o fastQC_output

	STAR \
		--genomeDir $STAR_INDEX \
		--sjdbGTFfile $GENOME_GTF \
		--runThreadN 8 \
		--outSAMstrandField intronMotif \
		--outFilterIntronMotifs RemoveNoncanonical \
		--outFileNamePrefix $WORKDIR/star_output/$basename \
		--readFilesIn $fq1 $fq2 \
		--outSAMtype BAM SortedByCoordinate \
		--outReadsUnmapped Fastx \
		--outSAMmode Full

	suffix="Aligned.sortedByCoord.out.bam"
	outname="$basename.count.txt"
	bam="$WORKDIR/star_output/$basename$suffix"
	$featureCounts \
		-T 8 \
		-t exon \
		-g gene_id \
		-a $GENOME_GTF \
		-o $WORKDIR/featureCount_output/$outname \
		$bam
done
