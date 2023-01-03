#!/bin/bash

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
STRTN NextSeq-automated-pipeline_scRNAseqAnalysis ver2023
EOS
  exit 1
}

# Parameter settings
PARAM=()
for opt in "$@"; do
    case "${opt}" in
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


# Run the R script
cp src/Example-BarcodesStages.txt out/
cd out
R CMD BATCH --slave --vanilla  ../bin/STRTN-scRNA-seq-analysis.R STRTN-scRNA-seq-analysis.R.log

