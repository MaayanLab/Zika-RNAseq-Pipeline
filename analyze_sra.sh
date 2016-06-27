#!/bin/bash
# This script start the pipeline from .sra files assuming 
# .sra files are located in SRR*/SRR*.sra relative to WORKDIR.

########## Step 0: check command line args and make sure files exist ##########
## Initialize variables with default values:
WORKDIR="data/"
GENOME="$HOME/genomes/Homo_sapiens/UCSC/hg19"
## Parse command line args
while [[ $# -gt 0 ]]; do
	key="$1"
	case $key in
		-g|--genome)
		GENOME="$2"
		shift # past argument
		;;
		-w|--workdir)
		WORKDIR="$2"
		shift # past argument
		;;
		-h|--help)
		echo "Usage: ./analyze_sra.sh -g <GENOME> -w <WORKDIR>"
		exit
		;;
		*)
		# unknown option
		echo "Unknown option: $key, exiting."
		echo "Usage: ./analyze_sra.sh -g <GENOME> -w <WORKDIR>"
		exit
		;;
	esac
	shift # past argument or value  
done


## Detect number of CPUs and use min(N_CPUS, 8) for jobs
N_CPUS=$(nproc)
N_CPUS=$(($N_CPUS>8?8:$N_CPUS))

## Check $WORKDIR
if [[ ! -d $WORKDIR ]]; then
	echo "Could not find working directory: $WORKDIR, exiting. Please make sure the working directory exists"
	exit 1
else
	## Convert to absolute paths
	GENOME=$(readlink -e $GENOME)
	WORKDIR=$(readlink -e $WORKDIR)
	echo "GENOME=$GENOME"
	echo "WORKDIR=$WORKDIR"
fi

## Check $GENOME
if [[ ! -d $GENOME ]]; then
	echo "Could not find reference genome: $GENOME, exiting. Please make sure the working directory exists"
	exit 1
else
	GENOME_GTF="$GENOME/Annotation/Genes/genes.gtf"
	GENOME_FA="$GENOME/Sequence/WholeGenomeFasta/genome.fa"
	if [[ ! -f $GENOME_GTF ]]; then
		echo "$GENOME_GTF not found, exiting"
		exit 1
	fi
	if [[ ! -f $GENOME_FA ]]; then
		echo "$GENOME_FA not found, exiting"
		exit 1
	fi
	STAR_INDEX="$GENOME/star/STAR_2.4.1c/"
fi
## Make STAR index if not exists
if [ ! -d $STAR_INDEX ]; then
	echo "STAR index does not exist, building STAR index"
	STAR \
		--runThreadN $N_CPUS \
		--runMode genomeGenerate \
		--genomeDir $STAR_INDEX \
		--genomeFastaFiles $GENOME_FA \
		--sjdbGTFfile $GENOME_GTF \
		--sjdbOverhang 100
fi


is_paired() {
	# function to examine whether a .sra file is paired-end sequencing 
	# ref: https://www.biostars.org/p/139422/
	local SRA="$1"
	local x=$(
		fastq-dump -I -X 1 -Z --split-spot "$SRA" 2>/dev/null \
		| awk '{if(NR % 2 == 1) print substr($1,length($1),1)}' \
		| uniq \
		| wc -l
	)
	# echo $SRA $x
	# $x should be 2 if paired-end, 1 if single-end
	if [ $x == 2 ]; then
		return 0 # true
	else
		return 1 # false
	fi
}

dump_sra() {
	# function to dump sra to fastq file with auto-detecting 
	# whether reads are from single-end or paired-end
	local sra="$1"
	if is_paired $sra; then
		echo "$sra is detected as paired-end sequencing reads"
		# Note that paired-end sequencing reads should be dumped into two fastq files
		fastq-dump --gzip -I --split-files -O paired_fastqs "$sra"
	else
		echo "$sra is detected as single-end sequencing reads"
		fastq-dump --gzip -O fastqs "$sra"
	fi
}

cd $WORKDIR

## create dirs if not exists
mkdir -p fastqs
mkdir -p paired_fastqs
mkdir -p fastQC_output
mkdir -p star_output
mkdir -p featureCount_output


########## Step 1: .sra -> .fastq ##########
## Dump .sra to .fastq
echo "Dumping .sra files to .fastq.gz files"
## Dump in parallel using xargs
export -f is_paired
export -f dump_sra
find SRR*/*.sra | xargs --max-args=1 --max-procs=$N_CPUS -I {} bash -c 'dump_sra "$@"' _ {}


########## Step 2: QC, align and assemble sequencing reads ##########
## Align and assemble single-end sequencing reads
echo "Started to align reads to the genome and assemble transcriptome"
cd fastqs
for fq in $(ls); do
	basename=$(echo $fq | cut -f1 -d '.')
	echo "Performing FastQC for $basename"
	fastqc $fq -o ../fastQC_output

	echo "Aligning reads from $basename to the reference genome"
	STAR \
		--genomeDir $STAR_INDEX \
		--sjdbGTFfile $GENOME_GTF \
		--runThreadN $N_CPUS \
		--outSAMstrandField intronMotif \
		--outFilterIntronMotifs RemoveNoncanonical \
		--outFileNamePrefix ../star_output/$basename \
		--readFilesIn $fq \
		--readFilesCommand zcat \
		--outSAMtype BAM Unsorted \
		--outReadsUnmapped Fastx \
		--outSAMmode Full

	suffix="Aligned.out.bam"
	outname="$basename.count.txt"
	bam="../star_output/$basename$suffix"
	featureCounts \
		-T $N_CPUS \
		-t exon \
		-g gene_id \
		-a $GENOME_GTF \
		-o ../featureCount_output/$outname \
		$bam
done

## Align and assemble paired-end sequencing reads
cd ../paired_fastqs
for basename in $(ls | cut -f1 -d '_' | sort | uniq); do
	echo $basename
	fq1="_1.fastq.gz"
	fq2="_2.fastq.gz"
	fq1=$basename$fq1
	fq2=$basename$fq2
	echo "Performing FastQC for $basename"
	fastqc $fq1 -o ../fastQC_output
	fastqc $fq2 -o ../fastQC_output
	echo "Aligning reads from $basename to the reference genome"
	STAR \
		--genomeDir $STAR_INDEX \
		--sjdbGTFfile $GENOME_GTF \
		--runThreadN $N_CPUS \
		--outSAMstrandField intronMotif \
		--outFilterIntronMotifs RemoveNoncanonical \
		--outFileNamePrefix ../star_output/$basename \
		--readFilesIn $fq1 $fq2 \
		--readFilesCommand zcat \
		--outSAMtype BAM Unsorted \
		--outReadsUnmapped Fastx \
		--outSAMmode Full

	suffix="Aligned.out.bam"
	outname="$basename.count.txt"
	bam="../star_output/$basename$suffix"
	featureCounts \
		-T $N_CPUS \
		-t exon \
		-g gene_id \
		-a $GENOME_GTF \
		-o ../featureCount_output/$outname \
		$bam
done
