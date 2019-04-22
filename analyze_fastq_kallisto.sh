#!/bin/bash
# This script start the pipeline (with kallisto instead of STAR+featureCount) from .fastq.gz files assuming 
# all single-end.fastq.gz files are located in fastqs/ relative to $WORKDIR and
# all paired-end.fastq.gz files are located in paired_fastqs/ relative to $WORKDIR.


########## Step 0: check command line args and make sure files exist ##########
## Initialize variables with default values:
WORKDIR="data/"
GENOME="$HOME/genomes/Homo_sapiens/UCSC/hg19"
skip_qc=false
n_cpus=8

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
		-t|--cpus)
		n_cpus="$2"
		shift # past argument
		;;        
		--skip-qc)
		skip_qc=true
		shift # past argument
		;;		
		-h|--help)
		echo "Usage: ./analyze_sra.sh -g <GENOME> -w <WORKDIR>"
		echo "Options:"
		echo "  --skip-qc: Skip fastQC steps for the fastq files"
		echo "  -t, --cpus: Number of CPUs to use, default to 8"
		exit
		;;
		*)
		# unknown option
		echo "Unknown option: $key, exiting."
		echo "Usage: ./analyze_sra.sh -g <GENOME> -w <WORKDIR>"
		echo "Options:"
		echo "  --skip-qc: Skip fastQC steps for the fastq files"
		echo "  -t, --cpus: Number of CPUs to use, default to 8"
		exit
		;;
	esac
	shift # past argument or value  
done


## Detect number of CPUs and use min(N_CPUS, n_cpus) for jobs
N_CPUS=$(nproc)
N_CPUS=$(($N_CPUS>n_cpus?n_cpus:$N_CPUS))
echo "Number of CPUs to be used: $N_CPUS"

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
	GENOME_GTF="$GENOME/Homo_sapiens.GRCh38.92.gtf.gz"
	GENOME_FA="$GENOME/Homo_sapiens.GRCh38.cdna.all.fa.gz"
	if [[ ! -f $GENOME_GTF ]]; then
		echo "$GENOME_GTF not found, exiting"
		exit 1
	fi
	if [[ ! -f $GENOME_FA ]]; then
		echo "$GENOME_FA not found, exiting"
		exit 1
	fi
	KALLISTO_INDEX="$GENOME/Homo_sapiens.GRCh38.cdna.all.idx"
fi
## Make kallisto index if not exists
if [ ! -f $KALLISTO_INDEX ]; then
	echo "kallisto index does not exist, building kallisto index"
	kallisto index \
		-i $KALLISTO_INDEX \
		$GENOME_FA
fi


cd $WORKDIR

## create dirs if not exists
mkdir -p fastQC_output
mkdir -p kallisto_output

########## Step 2: QC, align and assemble sequencing reads ##########
echo "Started to align reads to the genome and assemble transcriptome"
if [[ -d fastqs ]]; then
	## Align and assemble single-end sequencing reads
	cd fastqs
	for fq in $(ls); do
		basename=$(echo $fq | cut -f1 -d '.')
		if [ "$skip_qc" = false ]; then
			echo "Performing FastQC for $basename"
			fastqc $fq -o ../fastQC_output
		fi

		echo "Aligning reads from $basename to the reference genome"
        kallisto quant \
            -t $N_CPUS \
            -i $KALLISTO_INDEX \
            -o ../kallisto_output/$basename \
            $fq

	done
	cd ..
fi	

if [[ -d paired_fastqs ]]; then
	## Align and assemble paired-end sequencing reads
	cd paired_fastqs
	for basename in $(ls | cut -f1 -d '_' | sort | uniq); do
		echo $basename
		fq1="_1.fastq"
		fq2="_2.fastq"
		fq1=$basename$fq1
		fq2=$basename$fq2
		if [ "$skip_qc" = false ]; then
			echo "Performing FastQC for $basename"
			fastqc $fq1 -o ../fastQC_output
			fastqc $fq2 -o ../fastQC_output
		fi
		echo "Aligning reads from $basename to the reference genome"
        kallisto quant \
            -t $N_CPUS \
            -i $KALLISTO_INDEX \
            -o ../kallisto_output/$basename \
            $fq1 $fq2

	done
	cd ..
fi
