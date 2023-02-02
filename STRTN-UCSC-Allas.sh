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
-e, --email     Email address for hub file.
-n, --name	Hub name to store in Allas or hosting service. (default: STRTN-hub)
-h, --help      Show usage.
-v, --version   Show version.
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
NAME_HUB=STRTN-hub

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
	'-e' | '--email' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
		echo "${PROGNAME}: option requires an argument -- $( echo $1 | sed 's/^-*//' )" 1>&2
		exit 1
	    fi
	    email=true
	    EMAIL="$2"
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

# Creation the tracks for BAM files
shortLabel_hub="STRT-N data sequencing library tracks"
longLabel_hub="STRT-N data sequencing library tracks"
useOneFile=on
echo -e "hub "${NAME_HUB}"\n""shortLabel "${shortLabel_hub}"\n""longLabel "${longLabel_hub}"\n""useOneFile "${useOneFile}"\n""email "${EMAIL}"\n""\n""genome "${GENOME_VALUE}"\n" >> hub.txt

readarray -t bam_array < <(printf '%s\n' *.output.bam | sort -V)

# First barcode will be shown as example. Other barcodes will be hidden so one can set up configurations as desired.
file=${bam_array[0]}
name=$(basename $file .output.bam)
longLabel=${name}", read alignments tracks"
type=bam
visibility=squish
autoScale=on
alwaysZero=on
graphType=bar
windowingFunction=mean
bamColorMode=strand
colorByStrand="0,0,255 255,0,0"
bigDataUrl=$file
echo -e "track "${name}"\n""shortLabel "${name}"\n""longLabel "$longLabel"\n""type "${type}"\n""bigDataUrl "${bigDataUrl}"\n""visibility "${visibility}"\n""autoScale "${autoScale}"\n""alwaysZero "${alwaysZero}"\n""graphType "${graphType}"\n""windowingFunction "${windowingFunction}"\n""bamColorMode "${bamColorMode}"\n""colorByStrand "\"${colorByStrand}\""\n" >> hub.txt

for file in "${bam_array[@]:1:47}"
do
    name=$(basename $file .output.bam)
    longLabel=${name}", read alignments tracks"
    type=bam
    visibility=hide
    autoScale=on
    alwaysZero=on
    graphType=bar
    windowingFunction=mean
    bamColorMode=strand
    colorByStrand="0,0,255 255,0,0"
    bigDataUrl=${file}
    echo -e "track "${name}"\n""shortLabel "${name}"\n""longLabel "${longLabel}"\n""type "${type}"\n""bigDataUrl "${bigDataUrl}"\n""visibility "${visibility}"\n""autoScale "${autoScale}"\n""alwaysZero "${alwaysZero}"\n""graphType "${graphType}"\n""windowingFunction "${windowingFunction}"\n""bamColorMode "${bamColorMode}"\n""colorByStrand "\"${colorByStrand}\""\n" >> hub.txt
done

# Creation the tracks for bigWig files
# BigWig files
# First barcode will be shown as example.
file=${bam_array[0]}
name=$(basename $file .output.bam)
longLabel=${name}", forward strand histogram tracks"
type=bigWig
visibility=full
autoScale=on
alwaysZero=on
graphType=bar
windowingFunction=mean
color=0,0,100
bigDataUrl=${name}_plus.bw
echo -e "track "${name}_plus"\n""shortLabel "${name}"\n""longLabel "${longLabel}"\n""type "${type}"\n""bigDataUrl "${bigDataUrl}"\n""visibility "${visibility}"\n""autoScale "${autoScale}"\n""alwaysZero "${alwaysZero}"\n""graphType "${graphType}"\n""windowingFunction "${windowingFunction}"\n""color "${color}"\n" >> hub.txt

trackname=${name}", reverse strand"
longLabel=${name}", reverse strand histogram tracks"
color=100,0,0
bigDataUrl=${name}_minus.bw
echo -e "track "${name}_minus"\n""shortLabel "${name}"\n""longLabel "${longLabel}"\n""type "${type}"\n""bigDataUrl "${bigDataUrl}"\n""visibility "${visibility}"\n""autoScale "${autoScale}"\n""alwaysZero "${alwaysZero}"\n""graphType "${graphType}"\n""windowingFunction "${windowingFunction}"\n""color "${color}"\n" >> hub.txt

# Other samples
for file in "${bam_array[@]:1:47}"
do
    name=$(basename $file .output.bam)
    trackname=${name}", forward strand"
    longLabel=${name}", forward strand histogram tracks"
    type=bigWig
    visibility=hide
    autoScale=on
    alwaysZero=on
    graphType=bar
    windowingFunction=mean
    color=0,0,100
    bigDataUrl=${name}_plus.bw
    echo -e "track "${name}_plus"\n""shortLabel "${name}"\n""longLabel "${longLabel}"\n""type "${type}"\n""bigDataUrl "${bigDataUrl}"\n""visibility "${visibility}"\n""autoScale "${autoScale}"\n""alwaysZero "${alwaysZero}"\n""graphType "${graphType}"\n""windowingFunction "${windowingFunction}"\n""color "${color}"\n" >> hub.txt
    trackname=${name}", reverse strand"
    longLabel=${name}", reverse strand histogram tracks"
    color=100,0,0
    bigDataUrl=${name}_minus.bw
    echo -e "track "${name}_minus"\n""shortLabel "${name}"\n""longLabel "${longLabel}"\n""type "${type}"\n""bigDataUrl "${bigDataUrl}"\n""visibility "${visibility}"\n""autoScale "${autoScale}"\n""alwaysZero "${alwaysZero}"\n""graphType "${graphType}"\n""windowingFunction "${windowingFunction}"\n""color "${color}"\n" >> hub.txt
done

# Creation the track for bigBed file
trackname="coding_5´-end"
longLabel="coding 5'-end track"
type=bigBed
visibility=full
autoScale=on
colorByStrand="0,0,100 100,0,0"
bigDataUrl=coding_5end.bb
echo -e "track "${trackname}"\n""shortLabel "${trackname}"\n""longLabel "${longLabel}"\n""type "${type}"\n""bigDataUrl "${bigDataUrl}"\n""visibility "${visibility}"\n""autoScale "${autoScale}"\n""colorByStrand "\"${colorByStrand}\""\n" >> hub.txt

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

