#!/bin/bash
#SBATCH --job-name=fastq-fastqc
#SBATCH --time=01:00:00
#SBATCH --partition=small
#SBATCH --mem=25GB
#SBATCH --cpus-per-task=8

PROGNAME="$( basename $0 )"

# Usage
function usage() {
        cat <<EOS >&2
Usage: ${PROGNAME} [-w <path (required)>]
Options:
  -w, --working           /PATH/to/the working directory. Required!
  -h, --help              Show usage.
  -v, --version           Show version.
EOS
	exit 1
}

function version() {
          cat <<EOS >&2
STRT2-NextSeq-automated-pipeline_fastq-fastQC ver2022
EOS
  exit 1
}

# Parameter settings
PARAM=()
for opt in "$@"; do
    case "${opt}" in
	'-w' | '--working' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
		echo "${PROGNAME}: option requires an argument -- $( echo $1 | sed 's/^-*//' )" 1>&2
		exit 1
	    fi
	    WorkingDir=true
	    WorkingDir_PATH="$2"
	    shift 2
	    ;;
	'-h' | '--help' )
	    usage
	    ;;
	'-v' | '--version' )
	    version
	    ;;
	'--' | '-' )
	    shift
	    PARAM+=( "$@" )
	    break
	    ;;
	-* )
	    echo "${PROGNAME}: illegal option -- '$( echo $1 | sed 's/^-*//' )'" 1>&2
	    exit 1
	    ;;
    esac
done

if [[ -n "${PARAM[@]}" ]]; then
    usage
fi

[ "${WorkingDir}" != "true" ] && usage

module load biokit
module load tykky
export PATH="$WorkingDir_PATH/STRT2-env/bin:$PATH"

# Prepare fastq files
mkdir out/fastq
for file in out/Output_bam/*.output.bam; do name=$(basename $file .output.bam); samtools fastq -F 256 -F 1024 $file | pigz -c > out/fastq/$name.fq.gz; done

# Run fastQC
mkdir out/fastQC
fastqc -t 24 --nogroup -o out/fastQC out/fastq/*.fq.gz

# Run multiQC
OUTPUT_NAME=$(basename out/fastQC/*_1_fastqc.html _1_fastqc.html)
multiqc out/fastQC/ -n out/${OUTPUT_NAME}_MultiQC_report.html 
