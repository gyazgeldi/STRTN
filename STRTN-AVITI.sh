#!/bin/bash

PROGNAME="$( basename $0 )"

# Usage
function usage() {
    cat <<EOS >&2
Usage: ${PROGNAME} [-o <output>] [-g <genome (required)>] [-a <annotation>] [-R <R1 FASTQ (required)>] [-B <I1 FASTQ (required)>] [-i <path (required)>]

Options:
  -o, --out               Output file name. (default: OUTPUT)
  -g, --genome            Genome (hg19/hg38/mm9/mm10/mm39/canFam3/canFam4/canFam6/bosTau9). Required!
  -a, --annotation        Gene annotation (ref{RefSeq}/ens{Ensembl}/kg{UCSC KnownGenes}/wgEncodeGencodeBasis*) for QC and counting. Default : ref. NOTE: no Ensembl for hg38 & mm10 & mm39 & canFam4 & canFam6, no KnownGenes for canFam3 & canFam4 & canFam6 & bosTau9, no Gencode for mm9 & canFam3 & canFam4 & canFam6 & bosTau9.
  -R, --read-fastq        Path to the R1 FASTQ file. Required!
  -B, --barcode-fastq	  Path to the I1 barcode FASTQ file. Required!
  -i, --index             /PATH/to/the directory and basename of the HISAT2 index. Fasta file has to be 'basename.fasta'. Required!
  -c, --center            The name of the sequencing center that produced the reads. (default: CENTER)
  -r, --run               The barcode of the run. Prefixed to read names. (default: RUNBARCODE)
  -s, --structure         Read structure (default: "8M3S74T 6B")
  -d, --dta               Downstream-transcriptome-assembly for HISAT2, which is useful for TFE-based analysis but leads to fewer alignments with short-anchors.
  -h, --help              Show usage.
  -v, --version           Show version.
EOS
    exit 1
}

function version() {
      cat <<EOS >&2
STRT-N AVITI-automated-pipeline ver2026
EOS
  exit 1
}

# Default parameters
OUTPUT_NAME=OUTPUT
run_VALUE=RUNBARCODE
center_VALUE=FUGU
READ_STRUCTURE="8M3S74T 6B"
IF_DTA=true


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
	    elif  [[ "$2" =  "canFam4" ]]; then
		GENOME=true
		GENOME_VALUE="canFam4"
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
	'-R' | '--read-fastq' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
		echo "${PROGNAME}: option requires an argument -- $( echo $1 | sed 's/^-*//' )" 1>&2
		exit 1
	    fi
	    R1Fastq=true
	    R1Fastq_PATH="$2"
	    shift 2
	    ;;
	'-B' | '--barcode-fastq' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
		echo "${PROGNAME}: option requires an argument -- $( echo $1 | sed 's/^-*//' )" 1>&2
		exit 1
	    fi
	    I1Fastq=true
	    I1Fastq_PATH="$2"
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
[ "${R1Fastq}" != "true" ] && usage
[ "${I1Fastq}" != "true" ] && usage
[ "${Index}" != "true" ] && usage
[ "${ANNO}" != "true" ] && ANNO_VALUE=ref

# Make temporary and output directory
mkdir -p tmp
mkdir -p out

# Converting AVITI FASTQ to unmapped BAM file of entire one library
fgbio FastqToBam --input ${R1Fastq_PATH} ${I1Fastq_PATH} --read-structures ${READ_STRUCTURE} --output ${OUTPUT_NAME}_fastqtobam.bam --sample ${OUTPUT_NAME} --library ${OUTPUT_NAME}

# Demultiplexing 48 samples based on barcode.txt
tail -n +2 src/barcode.txt | while read barcode number; do samtools view -b -d BC:${barcode} ${OUTPUT_NAME}_fastqtobam.bam -o ${OUTPUT_NAME}_${number}.bam; done

# Mapping by HISAT2 and merging with original unmapped BAM
mkdir -p tmp/UMI
mkdir -p out/HISAT2_Metrics

for i in $(seq 1 48)
do
    file=${OUTPUT_NAME}_${i}.bam
    name=${OUTPUT_NAME}_${i}
    echo $name >> out/HISAT2_Metrics/Alignment-summary.txt

    picard SortSam I=$file O=tmp/${name}.unmapped.sorted.bam SORT_ORDER=queryname
    picard SamToFastq I=tmp/${name}.unmapped.sorted.bam F=tmp/${name}.fastq
    hisat2 -p 8 --dta -x ${Index_PATH} -U tmp/${name}.fastq -S tmp/${name}.mapped.sam 2>> out/HISAT2_Metrics/Alignment-summary.txt
    picard SortSam I=tmp/${name}.mapped.sam O=tmp/${name}.mapped.sorted.sam SORT_ORDER=queryname
    picard MergeBamAlignment UNMAPPED=tmp/${name}.unmapped.sorted.bam ALIGNED=tmp/${name}.mapped.sorted.sam O=tmp/UMI/${name}.umi.bam R=${Index_PATH}.fasta ATTRIBUTES_TO_RETAIN=XS
done

rm -f tmp/*.unmapped.sorted.bam
rm -f tmp/*.mapped.sorted.sam
rm -f tmp/*.mapped.sam
rm -f tmp/*.fastq

mkdir -p tmp/Unaligned_bam
mv ${OUTPUT_NAME}_fastqtobam.bam tmp/Unaligned_bam/
mv ${OUTPUT_NAME}_[0-9]*.bam tmp/Unaligned_bam/

# Mark potential PCR duplicates
mkdir -p out/MarkDuplicates_Metrics
for i in $(seq 1 48)
do
    picard MarkDuplicates INPUT=tmp/UMI/${OUTPUT_NAME}_${i}.umi.bam OUTPUT=out/${OUTPUT_NAME}_${i}.output.bam METRICS_FILE=out/MarkDuplicates_Metrics/${OUTPUT_NAME}_${i}.metrics.txt BARCODE_TAG=RX
done

# Preparation for annotation, extraction the information on genomic regions from annotation files and quality check
if [[ ${GENOME_VALUE} = "hg38" ]] && [[ ${ANNO_VALUE} =  "ens" ]]; then
    echo "No Ensembl gene annotations!! Please use RefSeq, KnownGenes, or Gencode for hg38"
    exit 1
elif [[ ${GENOME_VALUE} = "mm10" ]] && [[ ${ANNO_VALUE} =  "ens" ]]; then
    echo "No Ensembl gene annotations!! Please use RefSeq, KnownGenes, or Gencode for mm10"
    exit 1
elif [[ ${GENOME_VALUE} = "mm39" ]] && [[ ${ANNO_VALUE} =  "ens" ]]; then
    echo "No Ensembl gene annotations!! Please use RefSeq or KnownGenes, or Gencode for mm39"
    exit 1
elif [[ ${GENOME_VALUE} = "canFam4" ]] && [[ ${ANNO_VALUE} =  "ens" ]]; then
    echo "No Ensembl gene annotations!! Please use RefSeq for canFam4"
    exit 1  
elif [[ ${GENOME_VALUE} = "canFam6" ]] && [[ ${ANNO_VALUE} =  "ens" ]]; then
    echo "No Ensembl gene annotations!! Please use RefSeq for canFam6"
    exit 1   
elif [[ ${GENOME_VALUE} = "canFam3" ]] && [[ ${ANNO_VALUE} =  "kg" ]]; then
    echo "No KnownGenes annotations!! Please use RefSeq or Ensembl for canFam3"
    exit 1
elif [[ ${GENOME_VALUE} = "canFam4" ]] && [[ ${ANNO_VALUE} =  "kg" ]]; then
    echo "No KnownGenes annotations!! Please use RefSeq for canFam4"
    exit 1
elif [[ ${GENOME_VALUE} = "canFam6" ]] && [[ ${ANNO_VALUE} =  "kg" ]]; then
    echo "No KnownGenes annotations!! Please use RefSeq for canFam6"
    exit 1
elif [[ ${GENOME_VALUE} = "bosTau9" ]] && [[ ${ANNO_VALUE} =  "kg" ]]; then
    echo "No KnownGenes annotations!! Please use RefSeq or Ensembl for bosTau9"
    exit 1
elif [[ ${GENOME_VALUE} = "canFam3" ]] && [[ ${ANNO_VALUE} =  wgEncodeGencodeBasic* ]]; then
    echo "No Gencode annotations!! Please use RefSeq or Ensembl for canFam3"
    exit 1
elif [[ ${GENOME_VALUE} = "canFam4" ]] && [[ ${ANNO_VALUE} =  wgEncodeGencodeBasic* ]]; then
    echo "No Gencode annotations!! Please use RefSeq for canFam4"
    exit 1
elif [[ ${GENOME_VALUE} = "canFam6" ]] && [[ ${ANNO_VALUE} =  wgEncodeGencodeBasic* ]]; then
    echo "No Gencode annotations!! Please use RefSeq for canFam6"
    exit 1
elif [[ ${GENOME_VALUE} = "mm9" ]] && [[ ${ANNO_VALUE} =  wgEncodeGencodeBasic* ]]; then
    echo "No Gencode annotations!! Please use RefSeq, KnownGenes, or Ensembl for mm9"
    exit 1
elif [[ ${GENOME_VALUE} = "bosTau9" ]] && [[ ${ANNO_VALUE} =  wgEncodeGencodeBasic* ]]; then
    echo "No Gencode annotations!! Please use RefSeq or Ensembl for bosTau9"
    exit 1
elif [[ ${ANNO_VALUE} =  "ens" ]]; then
    echo "Downloading the Ensembl annotation data..."
    curl -o src/ensGene.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/${GENOME_VALUE}/database/ensGene.txt.gz
    curl -o src/ensemblToGeneName.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/${GENOME_VALUE}/database/ensemblToGeneName.txt.gz
    gunzip src/ensGene.txt.gz
    gunzip src/ensemblToGeneName.txt.gz
    join -1 1 -2 2 -t $'\t' <(sort -k 1,1 src/ensemblToGeneName.txt) <(sort -k 2,2 src/ensGene.txt) > src/common.txt
    join -1 1 -2 2 -t $'\t' -v 2 <(sort -k 1,1 src/ensemblToGeneName.txt) <(sort -k 2,2 src/ensGene.txt) | awk 'BEGIN{OFS="\t"}{print $1,$13,$2,$1=$2="",$0}' | cut -f 1-3,7- > src/no-genename.txt
    rm src/ensGene.txt && rm src/ensemblToGeneName.txt
    cat src/common.txt src/no-genename.txt > src/ens-genes.txt
    rm src/common.txt && rm src/no-genename.txt
    ruby bin/ENSEMBL-extract.rb
elif [[ ${ANNO_VALUE} =  "kg" ]]; then
    echo "Downloading the UCSC KnownGenes annotation data..."
    curl -o src/knownGene.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/${GENOME_VALUE}/database/knownGene.txt.gz
    curl -o src/kgXref.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/${GENOME_VALUE}/database/kgXref.txt.gz
    gunzip src/knownGene.txt.gz
    gunzip src/kgXref.txt.gz
    join  -1 1 -2 1 -t $'\t' <(sort -k 1,1 src/kgXref.txt | cut -f 1-5) <(sort -k 1,1 src/knownGene.txt) > src/knowngene-names.txt
    rm src/knownGene.txt && rm src/kgXref.txt
    ruby bin/KnownGenes-extract.rb
elif [[ ${ANNO_VALUE} =  "ref" ]]; then
    echo "Downloading the NCBI RefSeq annotation data..."
    curl -o src/refGene.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/${GENOME_VALUE}/database/refGene.txt.gz
    gunzip src/refGene.txt.gz
    ruby bin/RefSeq-extract.rb
elif [[ ${ANNO_VALUE} =  wgEncodeGencodeBasic* ]]; then
    echo "Downloading the Gencode annotation data..."
    curl -o src/Gencode.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/${GENOME_VALUE}/database/${ANNO_VALUE}.txt.gz
    gunzip src/Gencode.txt.gz
    ruby bin/GENCODE-extract.rb
else
    echo "Something is wrong with the annotation data file."
    exit 1
fi

# Chromosome size data and annotation regions
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

# Counting reads aligned to 5'end of genes
featureCounts -T 8 -s 1 --largestOverlap --ignoreDup --primary -a ../src/5end-regions.saf -F SAF -o ${OUTPUT_NAME}_byGene-counts.txt *.bam

mkdir -p Output_bai && mv *.bam.bai Output_bai
mkdir -p Output_bam && mv *.bam Output_bam

# Quality check plotting
# Run the R script
R CMD BATCH --slave --vanilla  ../bin/QC-plot.R QC-plot.R.log 
cd ..
