#!/bin/bash
WORKDIR='data/'

## paths relative to WORKDIR
GENOME='../genomes/Homo_sapiens/UCSC/hg19'
GENOME_GTF="$GENOME/Annotation/Genes/genes.gtf"
GENOME_FA="$GENOME/Sequence/WholeGenomeFasta/genome.fa"
STAR_INDEX="$GENOME/star/STAR_2.4.1c/"

cd $WORKDIR

## create dirs
mkdir fastqs
mkdir paired_fastqs
mkdir fastQC_output
mkdir star_output
mkdir featureCount_output

## dump .sra to .fastq
echo "Dumping .sra to .fastq"
for sra in $(ls SINGLE/SRR*/*.sra); do
	echo $sra
	fastq-dump -O fastqs $sra
done

## Note that paired-end sequencing reads should be dumped with different params
for sra in $(ls PAIRED/SRR*/*.sra); do
	echo $sra
	fastq-dump -I --split-files -O paired_fastqs $sra
done

## make star index if not exists
if [ ! -d $STAR_INDEX ]; then
	echo "STAR index does not exist, building STAR index"
	STAR \
		--runThreadN 8 \
		--runMode genomeGenerate \
		--genomeDir $STAR_INDEX \
		--genomeFastaFiles $GENOME_FA \
		--sjdbGTFfile $GENOME_GTF \
		--sjdbOverhang 100
fi

## align and assemble single-end sequencing reads
echo "Started to align reads to the genome and assemble transcriptome"
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
	featureCounts \
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
	featureCounts \
		-T 8 \
		-t exon \
		-g gene_id \
		-a $GENOME_GTF \
		-o $WORKDIR/featureCount_output/$outname \
		$bam
done
