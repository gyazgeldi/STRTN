#!/bin/bash
#SBATCH --job-name=uploadingDatatoAllas
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --mem=8GB

PROGNAME="$( basename $0 )"

# Usage
function usage() {
        cat <<EOS >&2
Usage: ${PROGNAME} [-w <path (required)>]
Options: 
-w, --working	/PATH/to/the working directory. Required! 
-n, --name	The hub in Allas. (default: STRT2-hub)	
EOS
	exit 1
}

function version() {
          cat <<EOS >&2
STRT-N NextSeq-automated-pipeline ver2022
EOS
  exit 1
}

# Default parameter.
NAME_HUB=STRT2-hub

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
	'-n' | '--name' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
		echo "${PROGNAME}: option requires an argument -- $( echo $1 | sed 's/^-*//' )" 1>&2
		exit 1
	    fi
	    name_hub=true
	    NAME_HUB="$2"
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
[ "${name_hub}" != "true" ] && usage

# Before running this script, please make sure that you accessed allas and are able to use swift protocol. Please follow below commands.
# module load allas
# allas-conf <project>

# Create new container for visualization files
swift post $NAME_HUB

# Upload .bam, .bai and .bw files
for file in *.output.bam
do
    name=$(basename $file .output.bam)
    swift upload $NAME_HUB ${WorkingDir_PATH}/$file
    swift upload $NAME_HUB ${WorkingDir_PATH}/$file.bai
    swift upload $NAME_HUB ${WorkingDir_PATH}/${name}_plus.bw
    swift upload $NAME_HUB ${WorkingDir_PATH}/${name}_minus.bw
    swift post $NAME_HUB --read-acl ".r:*"
done

# Upload coding_5´-end file and hub txt file
swift upload $NAME_HUB ${WorkingDir_PATH}/coding_5end.bb
swift upload $NAME_HUB ${WorkingDir_PATH}/hub.txt

# Make accessible and create sharable link
swift post $NAME_HUB --read-acl ".r:*"
account=$(swift stat | grep "Account: " | awk '{print $2}')
link=https://a3s.fi/swift/v1/$account/$NAME_HUB${WorkingDir_PATH}/hub.txt
echo -e "The sharable link is""\t"$link

rm *.output.bam && rm *.output.bam.bai 
rm *.bw && rm coding_5end.bb
