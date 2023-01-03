#!/bin/bash
#SBATCH --job-name=MainPipeline
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --mem=200GB
#SBATCH --cpus-per-task=8

PROGNAME="$( basename $0 )"

# Usage
function usage() {
    cat <<EOS >&2
Usage: ${PROGNAME} [-o <output>] [-g <genome (required)>] [-a <annotation>] [-b <path (required)>] [-i <path (required)>] [-w <path (required)>]
Options:
  -o, --out               Output file name. (default: OUTPUT)
  -g, --genome            Genome (hg19/hg38/mm9/mm10/mm39/canFam3/canFam6/bosTau9). Required!
  -a, --annotation        Gene annotation (ref{RefSeq}/ens{Ensembl}/kg{UCSC KnownGenes}/wgEncodeGencodeBasis*) for QC and counting. Default : ref. NOTE: no Ensembl for hg38 & mm10 & mm39 & canFam6, no KnownGenes for canFam3 & canFam6 & bosTau9, no Gencode for mm9 & canFam3 & canFam6 & bosTau9.
  -b, --basecalls         /PATH/to/the Illumina basecalls directory. Required!
  -i, --index             /PATH/to/the directory and basename of the HISAT2 index. Fasta file has to be 'basename.fasta'. Required!
  -w, --working           /PATH/to/the working directory. Required!
  -e, --email             Email address for hub file.
  -n, --name		  Hub name to store in Allas.
  -c, --center            The name of the sequencing center that produced the reads. (default: CENTER)
  -r, --run               The barcode of the run. Prefixed to read names. (default: RUNBARCODE)
  -s, --structure         Read structure (default: 8M3S75T6B)
  -d, --dta               Downstream-transcriptome-assembly for HISAT2, which is useful for TFE-based analysis but leads to fewer alignments with short-anchors.
  -h, --help              Show usage.
  -v, --version           Show version.
EOS
    exit 1
}

function version() {
      cat <<EOS >&2
STRT-N NextSeq-automated-pipeline ver2022
EOS
  exit 1
}

# Default parameters
OUTPUT_NAME=OUTPUT
run_VALUE=RUNBARCODE
center_VALUE=FUGU
READ_STRUCTURE=8M3S75T6B
IF_DTA=true
NAME_HUB=STRT2-hub

# Parameter settings
PARAM=()
for opt in "$@"; do
    case "${opt}" in
	'-o' | '--out' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
		echo "${PROGNAME}: option requires an argument -- $( echo $1 | sed 's/^-*//' )" 1>&2
		exit 1
	    fi
	    OUTNAME=true
	    OUTPUT_NAME="$2"
	    shift 2
	    ;;
	'-g' | '--genome' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
		echo "${PROGNAME}: option requires an argument -- $( echo $1 | sed 's/^-*//' )" 1>&2
		exit 1
	    elif [[ "$2" =  "hg19" ]]; then
		GENOME=true
		GENOME_VALUE="hg19"
		shift 2
	    elif [[ "$2" =  "hg38" ]]; then
		GENOME=true
		GENOME_VALUE="hg38"
		shift 2
	    elif [[ "$2" =  "mm9" ]]; then
		GENOME=true
		GENOME_VALUE="mm9"
		shift 2
	    elif [[ "$2" =  "mm10" ]]; then
		GENOME=true
		GENOME_VALUE="mm10"
		shift 2
	    elif [[ "$2" = "mm39" ]]; then
		GENOME=true
		GENOME_VALUE="mm39"
		shift 2
	    elif  [[ "$2" =  "canFam3" ]]; then
		GENOME=true
		GENOME_VALUE="canFam3"
		shift 2
	    elif  [[ "$2" =  "canFam6" ]]; then
		GENOME=true
		GENOME_VALUE="canFam6"
		shift 2
	    elif  [[ "$2" =  "bosTau9" ]]; then
		GENOME=true
		GENOME_VALUE="bosTau9"
		shift 2
	    else
		usage
		exit 1
	    fi
	    ;;
	'-a' | '--annotation' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
		echo "${PROGNAME}: option requires an argument -- $( echo $1 | sed 's/^-*//' )" 1>&2
		exit 1
	    elif [[ "$2" =  "ref" ]]; then
		ANNO=true
		ANNO_VALUE="ref"
		shift 2
	    elif [[ "$2" =  "kg" ]]; then
		ANNO=true
		ANNO_VALUE="kg"
		shift 2
	    elif [[ "$2" =  "ens" ]]; then
		ANNO=true
		ANNO_VALUE="ens"
		shift 2
	    elif [[ "$2" =  wgEncodeGencodeBasic* ]]; then
		ANNO=true
		ANNO_VALUE=$2
		shift 2
	    else
		usage
		exit 1
	    fi
	    ;;
	'-b' | '--basecalls' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
		echo "${PROGNAME}: option requires an argument -- $( echo $1 | sed 's/^-*//' )" 1>&2
		exit 1
	    fi
	    BaseCallsDir=true
	    BaseCallsDir_PATH="$2"
	    shift 2
	    ;;
	'-i' | '--index' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
		echo "${PROGNAME}: option requires an argument -- $( echo $1 | sed 's/^-*//' )" 1>&2
		exit 1
	    fi
	    Index=true
	    Index_PATH="$2"
	    shift 2
	    ;;
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
	'-c' | '--center' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
		echo "${PROGNAME}: option requires an argument -- $( echo $1 | sed 's/^-*//' )" 1>&2
		exit 1
	    fi
	    center_VALUE="$2"
	    shift 2
	    ;;
	'-r' | '--run' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
		echo "${PROGNAME}: option requires an argument -- $( echo $1 | sed 's/^-*//' )" 1>&2
		exit 1
	    fi
	    run_VALUE="$2"
	    shift 2
	    ;;
	'-s' | '--structure' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
		echo "${PROGNAME}: option requires an argument -- $( echo $1 | sed 's/^-*//' )" 1>&2
		exit 1
	    fi
	    READ_STRUCTURE="$2"
	    shift 2
	    ;;
	'-d' | '--dta' )
	    IF_DTA=true; shift
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

[ "${GENOME}" != "true" ] && usage
[ "${BaseCallsDir}" != "true" ] && usage
[ "${Index}" != "true" ] && usage
[ "${ANNO}" != "true" ] && ANNO_VALUE=ref
[ "${WorkingDir}" != "true" ] && usage

# Assign the base PATH for running the R scripts
r_PATH=$PATH

# Loading required tools
module load tykky
export PATH="${WorkingDir_PATH}/STRT2-env/bin:$PATH"

# Make temporary and output directory
mkdir tmp
mkdir out

# Preparation for barcodes
ALL_LINES=`cat src/barcode.txt | wc -l`
NLINES=`expr $ALL_LINES \- 1`

for i in `seq 1 $NLINES`
do
    echo -e ${OUTPUT_NAME}_${i}_Lane1.bam"\t"${OUTPUT_NAME}_${i}_Lane1"\t"${OUTPUT_NAME}_${i}_Lane1 >> tmp/out
done
paste tmp/out <(awk 'NR>1{print $1}' src/barcode.txt) | cut -f 1-4 > tmp/out2 && rm tmp/out
echo -e ${OUTPUT_NAME}_non-indexed_Lane1.bam"\t"${OUTPUT_NAME}_non-indexed_Lane1"\t" ${OUTPUT_NAME}_non-indexed_Lane1"\t"N >> tmp/out2
echo -e OUTPUT"\t"SAMPLE_ALIAS"\t"LIBRARY_NAME"\t"BARCODE_1  | cat - tmp/out2 > library.param.lane1 && rm tmp/out2

# Number of lanes
nlanes=`ls -l ${BaseCallsDir_PATH} | grep ^d | wc -l`
for i in `seq 2 $nlanes`
do
    sed -e "s/Lane1/Lane${i}/g" library.param.lane1 > library.param.lane${i}
done

# Convert BCL files to BAM files
for i in `seq 1 $nlanes`
do
    java -Xmx5g -Djava.io.tmpdir=tmp -jar /appl/soft/bio/picard/picard-tools-2.27.4/picard.jar ExtractIlluminaBarcodes \
	 BASECALLS_DIR=${BaseCallsDir_PATH}/ \
	 LANE=${i} \
	 READ_STRUCTURE=${READ_STRUCTURE} \
	 BARCODE_FILE=src/barcode.txt  \
	 METRICS_FILE=metrics_output_lane${i}.txt ;
    java -Xmx5g -Djava.io.tmpdir=tmp -jar /appl/soft/bio/picard/picard-tools-2.27.4/picard.jar IlluminaBasecallsToSam \
	 BASECALLS_DIR=${BaseCallsDir_PATH}/ \
	 LANE=${i} \
	 READ_STRUCTURE=${READ_STRUCTURE} \
	 RUN_BARCODE=${run_VALUE} \
	 IGNORE_UNEXPECTED_BARCODES=true \
	 LIBRARY_PARAMS=library.param.lane${i} \
	 SEQUENCING_CENTER=${center_VALUE} \
	 INCLUDE_NON_PF_READS=false
done

rm library.param.lane*
mkdir out/ExtractIlluminaBarcodes_Metrics && mv metrics_output_lane*.txt out/ExtractIlluminaBarcodes_Metrics

# Make the fasta reference / sequence dictionary if they do not exist.
if [[ ! -e ${Index_PATH}.fasta ]]; then
    hisat2-inspect ${Index_PATH} > ${Index_PATH}.fasta
    java -Xmx5g -Djava.io.tmpdir=tmp -jar /appl/soft/bio/picard/picard-tools-2.27.4/picard.jar CreateSequenceDictionary R=${Index_PATH}.fasta O=${Index_PATH}.dict
fi
if [[ ! -e ${Index_PATH}.dict ]]; then
    java -Xmx5g -Djava.io.tmpdir=tmp -jar /appl/soft/bio/picard/picard-tools-2.27.4/picard.jar CreateSequenceDictionary R=${Index_PATH}.fasta O=${Index_PATH}.dict
fi

# Mapping by HISAT2 and merging with the original unaligned BAM files to generate UMI-annotated BAM files
mkdir tmp/UMI
mkdir out/HISAT2_Metrics

if [ $IF_DTA = true ]; then
    for file in *.bam
    do
	name=$(basename $file .bam)
	echo $name >> out/HISAT2_Metrics/Alignment-summary.txt
	java -Xmx5g -Djava.io.tmpdir=tmp -jar /appl/soft/bio/picard/picard-tools-2.27.4/picard.jar SortSam \
	       I=$file \
	       O=tmp/.unmapped.sorted.bam \
	       SORT_ORDER=queryname;
	java -Xmx5g -Djava.io.tmpdir=tmp -jar /appl/soft/bio/picard/picard-tools-2.27.4/picard.jar SamToFastq \
	       I=tmp/.unmapped.sorted.bam\
	       F=tmp/.tmp.fastq \
	    | hisat2 -p 8 --dta -x ${Index_PATH} \
		     -U tmp/.tmp.fastq -S /dev/stdout \
		     2>> out/HISAT2_Metrics/Alignment-summary.txt  \
	    | java -Xmx5g -Djava.io.tmpdir=tmp -jar /appl/soft/bio/picard/picard-tools-2.27.4/picard.jar SortSam \
		     I=/dev/stdin \
		     O=tmp/.mapped.sorted.sam \
		     SORT_ORDER=queryname;
	java -Xmx5g -Djava.io.tmpdir=tmp -jar /appl/soft/bio/picard/picard-tools-2.27.4/picard.jar MergeBamAlignment \
	       ATTRIBUTES_TO_RETAIN=XS \
	       UNMAPPED=tmp/.unmapped.sorted.bam  \
	       ALIGNED=tmp/.mapped.sorted.sam \
	       O=tmp/UMI/$name.umi.bam \
	       R=${Index_PATH}.fasta
    done
else
    for file in *.bam
    do
	name=$(basename $file .bam)
	echo $name >> out/HISAT2_Metrics/Alignment-summary.txt
	java -Xmx5g -Djava.io.tmpdir=tmp -jar /appl/soft/bio/picard/picard-tools-2.27.4/picard.jar SortSam \
	       I=$file \
	       O=tmp/.unmapped.sorted.bam \
	       SORT_ORDER=queryname;
	java -Xmx5g -Djava.io.tmpdir=tmp -jar /appl/soft/bio/picard/picard-tools-2.27.4/picard.jar SamToFastq \
	       I=tmp/.unmapped.sorted.bam \
	       F=tmp/.tmp.fastq \
	    | hisat2 -p 8 -x ${Index_PATH} \
		     -U tmp/.tmp.fastq -S /dev/stdout \
		     2>> out/HISAT2_Metrics/Alignment-summary.txt  \
	    | java -Xmx5g -Djava.io.tmpdir=tmp -jar /appl/soft/bio/picard/picard-tools-2.27.4/picard.jar SortSam \
		     I=/dev/stdin \
		     O=tmp/.mapped.sorted.sam \
		     SORT_ORDER=queryname;
	java -Xmx5g -Djava.io.tmpdir=tmp -jar /appl/soft/bio/picard/picard-tools-2.27.4/picard.jar MergeBamAlignment \
	       ATTRIBUTES_TO_RETAIN=XS \
	       UNMAPPED=tmp/.unmapped.sorted.bam  \
	       ALIGNED=tmp/.mapped.sorted.sam \
	       O=tmp/UMI/$name.umi.bam \
	       R=${Index_PATH}.fasta
    done
fi

rm tmp/.unmapped.sorted.bam
rm tmp/.mapped.sorted.sam
mkdir tmp/merged
mkdir tmp/Unaligned_bam
mv *.bam tmp/Unaligned_bam

# Merging all lanes
for i in `seq 1 $NLINES`
do
    java -Xmx5g -Djava.io.tmpdir=tmp -jar /appl/soft/bio/picard/picard-tools-2.27.4/picard.jar MergeSamFiles \
	 $(printf "I=%s " tmp/UMI/${OUTPUT_NAME}_${i}_Lane*.umi.bam) \
	 O=/dev/stdout |
	java -Xmx5g -Djava.io.tmpdir=tmp -jar /appl/soft/bio/picard/picard-tools-2.27.4/picard.jar AddOrReplaceReadGroups \
	     I=/dev/stdin \
	     O=tmp/merged/${OUTPUT_NAME}_${i}.merged.bam \
	     RGLB=${OUTPUT_NAME}_${i} RGPL=NextSeq RGPU=${i} RGSM=${i}
done

rm -rf tmp/UMI

# Mark potential PCR duplicates
mkdir out/MarkDuplicates_Metrics
for i in `seq 1 $NLINES`
do
    java -Xmx5g -Djava.io.tmpdir=tmp -jar /appl/soft/bio/picard/picard-tools-2.27.4/picard.jar MarkDuplicates \
	 INPUT=tmp/merged/${OUTPUT_NAME}_${i}.merged.bam \
	 OUTPUT=out/${OUTPUT_NAME}_${i}.output.bam \
	 METRICS_FILE=out/MarkDuplicates_Metrics/${OUTPUT_NAME}_${i}.metrics.txt \
	 BARCODE_TAG=RX
done

rm -rf tmp/merged

# Preparation for annotation and QC
if [[ ${GENOME_VALUE} = "hg38" ]] && [[ ${ANNO_VALUE} =  "ens" ]]; then
    echo "No Ensembl gene annotations!! Please use RefSeq, KnownGenes, or Gencode for hg38"
    exit 1
elif [[ ${GENOME_VALUE} = "mm10" ]] && [[ ${ANNO_VALUE} =  "ens" ]]; then
    echo "No Ensembl gene annotations!! Please use RefSeq or KnownGenes, or Gencode for mm10"
    exit 1
elif [[ ${GENOME_VALUE} = "mm39" ]] && [[ ${ANNO_VALUE} =  "ens" ]]; then
    echo "No Ensembl gene annotations!! Please use RefSeq or KnownGenes, or Gencode for canFam6"
    exit 1
elif [[ ${GENOME_VALUE} = "canFam6" ]] && [[ ${ANNO_VALUE} =  "ens" ]]; then
    echo "No Ensembl gene annotations!! Please use RefSeq or KnownGenes, or Gencode for canFam6"
    exit 1   
elif [[ ${GENOME_VALUE} = "canFam3" ]] && [[ ${ANNO_VALUE} =  "kg" ]]; then
    echo "No KnownGenes annotations!! Please use RefSeq or Ensembl for canFam3"
    exit 1
elif [[ ${GENOME_VALUE} = "canFam6" ]] && [[ ${ANNO_VALUE} =  "kg" ]]; then
    echo "No KnownGenes annotations!! Please use RefSeq or Ensembl for canFam6"
    exit 1
elif [[ ${GENOME_VALUE} = "bosTau9" ]] && [[ ${ANNO_VALUE} =  "kg" ]]; then
    echo "No KnownGenes annotations!! Please use RefSeq or Ensembl for bosTau9"
    exit 1
elif [[ ${GENOME_VALUE} = "canFam3" ]] && [[ ${ANNO_VALUE} =  wgEncodeGencodeBasic* ]]; then
    echo "No Gencode annotations!! Please use RefSeq or Ensembl for canFam3"
    exit 1
elif [[ ${GENOME_VALUE} = "canFam6" ]] && [[ ${ANNO_VALUE} =  wgEncodeGencodeBasic* ]]; then
    echo "No Gencode annotations!! Please use RefSeq or Ensembl for canFam6"
    exit 1
elif [[ ${GENOME_VALUE} = "mm9" ]] && [[ ${ANNO_VALUE} =  wgEncodeGencodeBasic* ]]; then
    echo "No Gencode annotations!! Please use RefSeq, KnownGenes, or Ensembl for mm9"
    exit 1
elif [[ ${GENOME_VALUE} = "bosTau9" ]] && [[ ${ANNO_VALUE} =  wgEncodeGencodeBasic* ]]; then
    echo "No Gencode annotations!! Please use RefSeq, KnownGenes, or Ensembl for bosTau9"
    exit 1
elif [[ ${ANNO_VALUE} =  "ens" ]]; then
    echo "Downloading the Ensembl annotation data..."
    curl -o src/ensGene.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/${GENOME_VALUE}/database/ensGene.txt.gz
    curl -o src/ensemblToGeneName.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/${GENOME_VALUE}/database/ensemblToGeneName.txt.gz
    gunzip src/ensGene.txt.gz
    gunzip src/ensemblToGeneName.txt.gz
    join -1 1 -2 2 -t $'\t' <(sort -k 1,1 src/ensemblToGeneName.txt) <(sort -k 2,2 src/ensGene.txt) > src/common.txt
    join -1 1 -2 2 -t $'\t' -v 2 <(sort -k 1,1 src/ensemblToGeneName.txt) <(sort -k 2,2 src/ensGene.txt) | awk 'BEGIN{OFS="\t"}{print $2,$13,$1,$1=$2="",$0}' | cut -f 1-3,7- > src/no-genename.txt
    rm src/ensGene.txt && rm src/ensemblToGeneName.txt
    cat src/common.txt src/no-genename.txt > src/ens-genes.txt
    rm src/common.txt && rm src/no-genename.txt
    ruby bin/ENSEMBL-extract.rb
    shift 2
elif [[ ${ANNO_VALUE} =  "kg" ]]; then
    echo "Downloading the UCSC KnownGenes annotation data..."
    curl -o src/knownGene.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/${GENOME_VALUE}/database/knownGene.txt.gz
    curl -o src/kgXref.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/${GENOME_VALUE}/database/kgXref.txt.gz
    gunzip src/knownGene.txt.gz
    gunzip src/kgXref.txt.gz
    join  -1 1 -2 1 -t $'\t' <(sort -k 1,1 src/kgXref.txt | cut -f 1-5) <(sort -k 1,1 src/knownGene.txt) > src/knowngene-names.txt
    rm src/knownGene.txt && rm src/kgXref.txt
    ruby bin/KnownGenes-extract.rb
    shift 2
elif [[ ${ANNO_VALUE} =  "ref" ]]; then
    echo "Downloading the NCBI RefSeq annotation data..."
    curl -o src/refGene.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/${GENOME_VALUE}/database/refGene.txt.gz
    gunzip src/refGene.txt.gz
    ruby bin/RefSeq-extract.rb
    shift 2
elif [[ ${ANNO_VALUE} =  wgEncodeGencodeBasic* ]]; then
    echo "Downloading the Gencode annotation data..."
    curl -o src/Gencode.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/${GENOME_VALUE}/database/${ANNO_VALUE}.txt.gz
    gunzip src/Gencode.txt.gz
    ruby bin/GENCODE-extract.rb
    shift 2
else
    echo "Something is wrong with the annotation data file."
    exit 1
fi

echo "Downloading the chromosome size data..."
curl -o src/${GENOME_VALUE}.chrom.sizes http://hgdownload.soe.ucsc.edu/goldenPath/${GENOME_VALUE}/bigZips/${GENOME_VALUE}.chrom.sizes
cat src/${GENOME_VALUE}.chrom.sizes | awk '{print $1"\t"1"\t"$2}' | sortBed -i > src/chrom.size.bed
cat src/proxup.bed | grep -v _alt  | grep -v _hap | grep -v _fix | grep -v _random | grep -v ^chrUn | sortBed -i stdin | intersectBed -a stdin -b src/chrom.size.bed > src/proxup_trimmed.bed
cat src/5utr.bed src/proxup_trimmed.bed | grep -v _alt | grep -v _hap | grep -v _fix | grep -v _random | grep -v ^chrUn | sortBed -i stdin | mergeBed -s -o distinct,distinct,distinct -c 4,5,6 -i - | grep -v , > src/coding_5end.bed
cat src/exon.bed src/proxup_trimmed.bed | grep -v _alt | grep -v _hap | grep -v _fix | grep -v _random | grep -v ^chrUn | sortBed -i stdin | mergeBed -s -o distinct,distinct,distinct -c 4,5,6 -i - > src/coding.bed
cat src/ERCC.bed src/coding_5end.bed | awk '{print $4 "\t" $1 "\t" $2+1 "\t" $3 "\t" $6}' > src/5end-regions.saf

rm src/${GENOME_VALUE}.chrom.sizes
rm src/5utr.bed
rm src/exon.bed
rm src/proxup.bed
rm src/proxup_trimmed.bed

# Quality check
cd out
echo -e Barcode"\t"Qualified_reads"\t"Total_reads"\t"Redundancy"\t"Mapped_reads"\t"Mapped_rate\
     "\t"Spikein_reads"\t"Spikein-5end_reads"\t"Spikein-5end_rate"\t"Coding_reads"\t"Coding-5end_reads"\t"Coding-5end_rate > ${OUTPUT_NAME}-QC.txt

for file in *.output.bam
do
    name=$(basename $file .output.bam)
    samtools index $file
    QR=$(samtools view -F 256 $file | wc -l)
    Total=$(samtools view -F 256 -F 1024 $file | wc -l)
    Redundancy=$(echo "scale=2;$QR/$Total" | bc)
    Map=$(samtools view -F 256 -F 1024 -F 4 $file | wc -l)
    Rate=$(echo "scale=1;$Map*100/$Total" | bc)
    Spike=$(samtools view -F 256 -F 1024 -F 4 $file |grep -e ERCC -e NIST| wc -l)
    spikein_5end_reads=$(samtools view -u -F 256 -F 1024 -F 4 $file | intersectBed -abam stdin -wa -bed -b ../src/ERCC.bed | cut -f 4 | sort -u | wc -l)
    spikein_5end_rate=$(echo "scale=1;$spikein_5end_reads*100/$Spike" | bc)
    coding_reads=$(samtools view -u -F 256 -F 1024 -F 4 $file | intersectBed -abam stdin -wa -bed -b ../src/coding.bed | cut -f 4 | sort -u | wc -l)
    coding_5end_reads=$(samtools view -u -F 256 -F 1024 -F 4 $file | intersectBed -abam stdin -wa -bed -b ../src/coding_5end.bed | cut -f 4 | sort -u | wc -l)
    coding_5end_rate=$(echo "scale=1;$coding_5end_reads*100/$coding_reads" | bc)
    echo -e $name"\t"$QR"\t"$Total"\t"$Redundancy"\t"$Map"\t"$Rate"\t"$Spike"\t"$spikein_5end_reads"\t"$spikein_5end_rate"\t"$coding_reads"\t"$coding_5end_reads"\t"$coding_5end_rate >> ${OUTPUT_NAME}-QC.txt
done

# Counting by featureCounts
featureCounts -T 8 -s 1 --largestOverlap --ignoreDup --primary -a ../src/5end-regions.saf -F SAF -o ${OUTPUT_NAME}_byGene-counts.txt *.bam

mkdir Output_bai && mv *.bam.bai Output_bai
mkdir Output_bam && mv *.bam Output_bam

# Quality check plotting
# Change the path to run R script
export "PATH=$r_PATH"
module load r-env

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temporary folder path
echo "TMPDIR=${WorkingDir_PATH}" >> ~/.Renviron

# Run the R script
srun apptainer_wrapper exec Rscript --no-save ${WorkingDir_PATH}/bin/QC-plot.R

cd ..

# Change the path to run script with conda packages
module load biokit
module load tykky
export PATH="${WorkingDir_PATH}/STRT2-env/bin:$PATH"

# Removing reads that have flags 
for file in ${WorkingDir_PATH}/out/Output_bam/*.output.bam
do
    name=$(basename $file .output.bam)
    samtools view -u -F 256 -F 1024 -F 4 $file > ${name}.output.bam
    samtools index ${name}.output.bam
done

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
alwaysZero=on
graphType=bar
windowingFunction=mean
bamColorMode=strand
colorByStrand="0,0,255 255,0,0"
bigDataUrl=$file
echo -e "track "${name}"\n""shortLabel "${name}"\n""longLabel "$longLabel"\n""type "${type}"\n""bigDataUrl "${bigDataUrl}"\n""visibility "${visibility}"\n""alwaysZero "${alwaysZero}"\n""graphType "${graphType}"\n""windowingFunction "${windowingFunction}"\n""bamColorMode "${bamColorMode}"\n""colorByStrand "\"${colorByStrand}\""\n" >> hub.txt

for file in "${bam_array[@]:1:47}"
do
    name=$(basename $file .output.bam)
    longLabel=${name}", read alignments tracks"
    type=bam
    visibility=hide
    alwaysZero=on
    graphType=bar
    windowingFunction=mean
    bamColorMode=strand
    colorByStrand="0,0,255 255,0,0"
    bigDataUrl=${file}
    echo -e "track "${name}"\n""shortLabel "${name}"\n""longLabel "${longLabel}"\n""type "${type}"\n""bigDataUrl "${bigDataUrl}"\n""visibility "${visibility}"\n""alwaysZero "${alwaysZero}"\n""graphType "${graphType}"\n""windowingFunction "${windowingFunction}"\n""bamColorMode "${bamColorMode}"\n""colorByStrand "\"${colorByStrand}\""\n" >> hub.txt
done

# Create spike-ins size file and chromosomal size file
wget https://www-s.nist.gov/srmors/certificates/documents/SRM2374_putative_T7_products_NoPolyA_v2.FASTA
seqtk seq -L 5 SRM2374_putative_T7_products_NoPolyA_v2.FASTA > spikeins_size.fasta
samtools faidx spikeins_size.fasta
cut -f1-2 spikeins_size.fasta.fai > spikeins_sizes

wget https://hgdownload.soe.ucsc.edu/goldenPath/${GENOME_VALUE}/bigZips/${GENOME_VALUE}.chrom.sizes
cat spikeins_sizes >> ${GENOME_VALUE}.chrom.sizes
mv ${GENOME_VALUE}.chrom.sizes ${GENOME_VALUE}.chrom.sizes_with_spike_ins
rm SRM2374_putative_T7_products_NoPolyA_v2.FASTA
rm spikeins_*

cd out

# Change the path to run R script
export "PATH=$r_PATH"

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temporary folder path
echo "TMPDIR=${WorkingDir_PATH}" >> ~/.Renviron

# Run the R script to find scale factor for each sample
srun apptainer_wrapper exec Rscript --no-save ${WorkingDir_PATH}/bin/calculationScaleFactor.R

cd ..

# Change the path to run script with conda packages
module load biokit
module load tykky
export PATH="${WorkingDir_PATH}/STRT2-env/bin:$PATH"

# Create bedGraph files and bigWig files for forward and reverse strand for each sample
for file in *.output.bam
do
    name=$(basename $file .output.bam)
    scale_factor=$(grep $file ${WorkingDir_PATH}/out/scale_factor_df.txt | awk '{print $2}')
    bedtools genomecov -5 -strand + -ibam $file -bg -scale $scale_factor | sort -k1,1 -k2,2n > ${WorkingDir_PATH}/${name}_plus.bedgraph
    bedtools genomecov -5 -strand - -ibam $file -bg -scale $scale_factor | sort -k1,1 -k2,2n > ${WorkingDir_PATH}/${name}_minus.bedgraph
    bedGraphToBigWig ${WorkingDir_PATH}/${name}_plus.bedgraph ${GENOME_VALUE}.chrom.sizes_with_spike_ins ${WorkingDir_PATH}/${name}_plus.bw
    bedGraphToBigWig ${WorkingDir_PATH}/${name}_minus.bedgraph ${GENOME_VALUE}.chrom.sizes_with_spike_ins ${WorkingDir_PATH}/${name}_minus.bw
done

# Creation the tracks for bigWig files
# BigWig files
# First barcode will be shown as example.
file=${bam_array[0]}
name=$(basename $file .output.bam)
longLabel=${name}", forward strand histogram tracks"
type=bigWig
visibility=full
alwaysZero=on
graphType=bar
windowingFunction=mean
color=0,0,100
bigDataUrl=${name}_plus.bw
echo -e "track "${name}_plus"\n""shortLabel "${name}"\n""longLabel "${longLabel}"\n""type "${type}"\n""bigDataUrl "${bigDataUrl}"\n""visibility "${visibility}"\n""alwaysZero "${alwaysZero}"\n""graphType "${graphType}"\n""windowingFunction "${windowingFunction}"\n""color "${color}"\n" >> hub.txt

trackname=${name}", reverse strand"
longLabel=${name}", reverse strand histogram tracks"
color=100,0,0
bigDataUrl=${name}_minus.bw
echo -e "track "${name}_minus"\n""shortLabel "${name}"\n""longLabel "${longLabel}"\n""type "${type}"\n""bigDataUrl "${bigDataUrl}"\n""visibility "${visibility}"\\n""alwaysZero "${alwaysZero}"\n""graphType "${graphType}"\n""windowingFunction "${windowingFunction}"\n""color "${color}"\n" >> hub.txt

# Other samples
for file in "${bam_array[@]:1:47}"
do
    name=$(basename $file .output.bam)
    trackname=${name}", forward strand"
    longLabel=${name}", forward strand histogram tracks"
    type=bigWig
    visibility=hide
    alwaysZero=on
    graphType=bar
    windowingFunction=mean
    color=0,0,100
    bigDataUrl=${name}_plus.bw
    echo -e "track "${name}_plus"\n""shortLabel "${name}"\n""longLabel "${longLabel}"\n""type "${type}"\n""bigDataUrl "${bigDataUrl}"\n""visibility "${visibility}"\n""alwaysZero "${alwaysZero}"\n""graphType "${graphType}"\n""windowingFunction "${windowingFunction}"\n""color "${color}"\n" >> hub.txt
    trackname=${name}", reverse strand"
    longLabel=${name}", reverse strand histogram tracks"
    color=100,0,0
    bigDataUrl=${name}_minus.bw
    echo -e "track "${name}_minus"\n""shortLabel "${name}"\n""longLabel "${longLabel}"\n""type "${type}"\n""bigDataUrl "${bigDataUrl}"\n""visibility "${visibility}"\n""alwaysZero "${alwaysZero}"\n""graphType "${graphType}"\n""windowingFunction "${windowingFunction}"\n""color "${color}"\n" >> hub.txt
done

# Create coding_5´end bigBed file
cp ${WorkingDir_PATH}/src/coding_5end.bed ${WorkingDir_PATH}/
sort -k1,1 -k2,2n coding_5end.bed > sorted_coding_5end.bed
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
chmod a+x bedToBigBed
wget https://hgdownload.soe.ucsc.edu/goldenPath/${GENOME_VALUE}/bigZips/${GENOME_VALUE}.chrom.sizes
cat ${GENOME_VALUE}.chrom.sizes > ${GENOME_VALUE}.chrom.sizes_without_spike_ins
./bedToBigBed sorted_coding_5end.bed ${GENOME_VALUE}.chrom.sizes_without_spike_ins coding_5end.bb

# Creation the track for bigBed file
trackname="coding_5´-end"
longLabel="coding 5'-end track"
type=bigBed
visibility=full
colorByStrand="0,0,100 100,0,0"
bigDataUrl=coding_5end.bb
echo -e "track "${trackname}"\n""shortLabel "${trackname}"\n""longLabel "${longLabel}"\n""type "${type}"\n""bigDataUrl "${bigDataUrl}"\n""visibility "${visibility}"\n""colorByStrand "\"${colorByStrand}\""\n" >> hub.txt

rm coding_5end.bed
rm .bedgraph*
rm ${GENOME_VALUE}.chrom.sizes
rm ${GENOME_VALUE}.chrom.sizes_with_spike_ins
rm out/scale_factor_df.txt
rm sorted_coding_5end.bed
rm ${GENOME_VALUE}.chrom.sizes_without_spike_ins
rm bedToBigBed