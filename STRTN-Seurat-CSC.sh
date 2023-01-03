#!/bin/bash
#SBATCH --job-name=scRNA-seq-analysis
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --mem=20GB

PROGNAME="$( basename $0 )"

# Usage
function usage() {
    cat <<EOS >&2
Usage: ${PROGNAME} [-w <working (required)>]
Options:
  -w, --working		  /PATH/to/the working directory. Required!
  -h, --help              Show usage.
  -v, --version           Show version.
EOS
    exit 1
}

function version() {
      cat <<EOS >&2
STRTN NextSeq-automated-pipeline_scRNAseqAnalysis ver2022
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
    esac
done

if [[ -n "${PARAM[@]}" ]]; then
    usage
fi

[ "${WorkingDir}" != "true" ] && usage

# Load the required module
module load r-env

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temporary folder path
echo "TMPDIR=${WorkingDir_PATH}" >> ~/.Renviron

# Run the R script
cp src/Example-BarcodesStages.txt out/
cd out
srun apptainer_wrapper exec Rscript --no-save ${WorkingDir_PATH}/bin/STRTN-scRNA-seq-analysis.R

