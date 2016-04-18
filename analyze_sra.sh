#!/bin/bash
WORKDIR='/home/maayanlab/Zika'
SRA_BIN='/home/maayanlab/Downloads/sratoolkit.2.5.7-ubuntu64/bin/'
cmd_name='fastq-dump'
featureCounts='/home/maayanlab/Downloads/subread-1.4.6-p2-Linux-x86_64/bin/featureCounts'

GENOME_GTF='/home/maayanlab/Zika/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf'
STAR_INDEX='/home/maayanlab/Zika/Homo_sapiens/UCSC/hg19/star/STAR_2.4.1c/'

cd $WORKDIR

## create dirs
# mkdir fastqs
# mkdir star_output
# mkdir featureCount_output
# ## dump .sra to .fastq
# for sra in $(ls SRR*/*.sra); do
# 	echo $sra
# 	$SRA_BIN$cmd_name -O fastqs $sra
# done

## make star index
# STAR \
#     --runThreadN 8 \
#     --runMode genomeGenerate \
#     --genomeDir /home/maayanlab/Zika/Homo_sapiens/UCSC/hg19/star/STAR_2.4.1c/ \
#     --genomeFastaFiles /home/maayanlab/Zika/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa \
#     --sjdbGTFfile /home/maayanlab/Zika/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf \
#     --sjdbOverhang 100

cd fastqs
for fq in $(ls); do
	basename=$(echo $fq | cut -f1 -d '.')
	echo $basename
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
