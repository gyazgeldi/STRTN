#!/bin/bash

PROGNAME="$( basename $0 )"

# Usage
function usage() {
      cat << EOS >&2
Usage: ${PROGNAME} 
Options:
  -c, --coverage          Minimum read coverage allowed for the predicted transcripts. (default: 5)
  -l, --length            Minimum length allowed for the predicted transcripts. (default: 75)
  -h, --help              Show usage.
  -v, --version           Show version.
EOS
  exit 1
}

function version() {
  cat << EOS >&2
STRT-N NextSeq-automated-pipeline_TFE-based ver2022
EOS
  exit 1
}

# Default parameters
cover_VALUE=5
len_VALUE=75

# Parameter settings
PARAM=()
for opt in "$@"; do
    case "${opt}" in
    '-c' | '--coverage' )
            if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
                echo "${PROGNAME}: option requires an argument -- $( echo $1 | sed 's/^-*//' )" 1>&2
                exit 1
            fi
            cover_VALUE="$2"
            shift 2
            ;;
    '-l' | '--length' )
            if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
                echo "${PROGNAME}: option requires an argument -- $( echo $1 | sed 's/^-*//' )" 1>&2
                exit 1
            fi
            len_VALUE="$2"
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

# Make temporary and output directories
mkdir byTFE_tmp
mkdir byTFE_out
mkdir byTFE_tmp/class

OUTPUT_NAME=$(basename out/Output_bam/*_1.output.bam _1.output.bam)

# Sample classification
while read row; do
    column1=`echo ${row} | cut -d ' ' -f 1`
    column2=`echo ${row} | cut -d ' ' -f 2`
  if [[ $column2 != "NA" ]]; then
    mkdir -p byTFE_tmp/class/${column2}
cp out/Output_bam/${OUTPUT_NAME}_${column1}.output.bam byTFE_tmp/class/${column2}
  else
    :
  fi
done < src/TFEclass.txt

classes=`find byTFE_tmp/class/* -type d` 
for class in $classes;
do
CLASS_NAME=$(basename $class byTFE_tmp/class/)
  # Merge all BAMs, remove duplicated, non-primary, unmapped reads, and sort
  samtools merge -@ 8 - $class/*.bam | samtools view -@ 8 -b -F 256 -F 1024 -F 4 - | samtools sort -@ 8 -o $class/merged.bam

  # Assembly with Stringtie 
  stringtie $class/merged.bam -o $class/stringtie.gtf -p 8 -m ${len_VALUE} --fr -l ${OUTPUT_NAME}.${CLASS_NAME} -c ${cover_VALUE}

  # Extract 1st-exon
  cat $class/stringtie.gtf | awk '{if($7=="+"||$7=="."){print $0}}'| grep 'exon_number "1"' \
      | awk 'OFS =  "\t" {print $1,$4-1,$5,$12,"0",$7}' | sed -e 's/"//g'| sed -e 's/;//g' > $class/firstExons-fwd.bed
  cat $class/stringtie.gtf | awk 'BEGIN{OFS="\t"}{if($7=="-" && $3=="exon"){print $1,$4-1,$5,$12,"0",$7}}' \
      | sed -e 's/"//g'| sed -e 's/;//g' | sort -k 4,4 -k 1,1 -k 2,2n | bedtools groupby -i stdin -g 4  -c 1,2,3,4,5,6 -o last \
      | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5,$6,$7}' > $class/firstExons-rev.bed
  cat $class/firstExons-fwd.bed  $class/firstExons-rev.bed | sortBed -i stdin > $class/firstExons.bed
  rm $class/firstExons-fwd.bed && rm $class/firstExons-rev.bed

  # Fiveprimes for peak detection
  mkdir $class/bedGraph
  for file in $class/*.output.bam
  do
  name=$(basename $file .output.bam)
  Spike=$(samtools view -F 256 -F 1024 -F 4 $file |grep -e ERCC -e NIST| wc -l)
  samtools view -b -F 256 -F 1024 -F 4 $file | bamToBed -i stdin\
  | gawk 'BEGIN{ FS="\t"; OFS=" " }{if($6=="+"){print $1,$2,$2+1,".",0,"+"}else{print $1,$3-1,$3,".",0,"-"}}'\
  | sort -k 1,1 -k 2,2n\
  | uniq -c\
  | gawk 'BEGIN{ FS=" "; OFS="\t" }{print $2,$3,$4,$5,$1/'$Spike',$7}'\
  | pigz -c > $class/bedGraph/$name.bedGraph.gz
  done
  gunzip -c $class/bedGraph/*.bedGraph.gz | sort -k 1,1 -k 2,2n | mergeBed -s -c 4,5,6 -o distinct,sum,distinct -d -1  > $class/fivePrimes.bed
done
	     
# TFE annotation
cat byTFE_tmp/class/*/firstExons.bed | sort -k 1,1 -k 2,2n |awk '{if($6=="+"){print $0}}' | grep -e ERCC -e NIST \
| mergeBed -s -c 6 -o distinct | bedtools groupby -i stdin -g 1  -c 1,2,3,4 -o first \
| awk 'BEGIN{OFS="\t"}{print $2,$3,$4,"RNA_SPIKE_"$2,0,$5}' > byTFE_tmp/${OUTPUT_NAME}_spike-firstExons_class.bed 
cat byTFE_tmp/class/*/firstExons.bed | sort -k 1,1 -k 2,2n | mergeBed -s -c 6 -o distinct \
| grep -v ERCC| grep -v NIST | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"TFE"NR,0,$4}' > byTFE_tmp/${OUTPUT_NAME}_nonspike-firstExons_class.bed 
cat byTFE_tmp/${OUTPUT_NAME}_spike-firstExons_class.bed  byTFE_tmp/${OUTPUT_NAME}_nonspike-firstExons_class.bed > byTFE_tmp/${OUTPUT_NAME}_TFE-regions.bed
rm byTFE_tmp/${OUTPUT_NAME}_spike-firstExons_class.bed && rm byTFE_tmp/${OUTPUT_NAME}_nonspike-firstExons_class.bed

# Counting
awk '{print $4 "\t" $1 "\t" $2+1 "\t" $3 "\t" $6}' byTFE_tmp/${OUTPUT_NAME}_TFE-regions.bed > byTFE_tmp/${OUTPUT_NAME}_TFE-regions.saf
featureCounts -T 8 -s 1 --largestOverlap --ignoreDup --primary -a byTFE_tmp/${OUTPUT_NAME}_TFE-regions.saf -F SAF -o byTFE_tmp/${OUTPUT_NAME}_byTFE-counts.txt byTFE_tmp/class/*/*.output.bam

# Peaks
cat byTFE_tmp/class/*/fivePrimes.bed | sort -k 1,1 -k 2,2n | mergeBed -s -c 4,5,6 -o distinct,sum,distinct -d -1  > byTFE_tmp/${OUTPUT_NAME}_fivePrimes.bed
intersectBed -wa -wb -s -a byTFE_tmp/${OUTPUT_NAME}_TFE-regions.bed -b byTFE_tmp/${OUTPUT_NAME}_fivePrimes.bed \
  | cut -f 4,7,8,9,11,12 \
  | gawk 'BEGIN{ FS="\t"; OFS="\t" }{p=$6=="+"?$3:-$4;print $2,$3,$4,$1,$5,$6,p,$1}' \
  | sort -k 8,8 -k 5,5gr -k 7,7g \
  | uniq -f 7 \
  | cut -f 1-6 \
  | sort -k 1,1 -k 2,2n > byTFE_out/${OUTPUT_NAME}_peaks.bed

# Annotation of peaks
mkdir src/anno
if test -f src/ens-genes.txt && test ! -f src/knowngene-names.txt && test ! -f src/refGene.txt && test ! -f src/Gencode.txt; then
  echo "Annotation with Ensembl"
  ruby bin/ensGene_annotation.rb
  shift 2
elif test ! -f src/ens-genes.txt && test -f src/knowngene-names.txt && test ! -f src/refGene.txt && test ! -f src/Gencode.txt; then
  echo "Annotation with UCSC KnownGenes"
  ruby bin/knownGene_annotation.rb
  shift 2
elif test ! -f src/ens-genes.txt && test ! -f src/knowngene-names.txt && test -f src/refGene.txt && test ! -f src/Gencode.txt; then
  echo "Annotation with NCBI RefSeq"
  ruby bin/refGene_annotation.rb
  shift 2
elif test ! -f src/ens-genes.txt && test ! -f src/knowngene-names.txt && test ! -f src/refGene.txt && test -f src/Gencode.txt; then
  echo "Annotation with GENCODE"
  ruby bin/GENCODE_annotation.rb
  shift 2
else
  echo "Something is wrong with the annotation data file."
  exit 1
fi

intersectBed -a src/anno/Coding-up.bed -b src/chrom.size.bed > src/anno/Coding-up_trimmed.bed
intersectBed -a src/anno/NC-up.bed -b src/chrom.size.bed > src/anno/NC-up_trimmed.bed

# Coding 5'UTR
intersectBed -s -wa -wb -a byTFE_out/${OUTPUT_NAME}_peaks.bed -b src/anno/Coding-5UTR.bed | awk -F "\t" '{print($4,$10)}' \
|awk -F "|" '{if(a[$1])a[$1]=a[$1]";"$2; else a[$1]=$2;}END{for (i in a)print i, a[i];}' OFS="\t" \
|awk -F " " '{print $1"\t"$2","$3}' \
|awk -F "\t" '{if(a[$1])a[$1]=a[$1]":"$2; else a[$1]=$2;}END{for (i in a)print i, a[i];}' OFS="\t" \
|awk -v 'OFS=\t' '{print $1,$2,"Coding_5UTR"}' | sort -k 1,1 > src/anno/peaks_class1.txt

# Coding upstream
intersectBed -s -wa -wb -a byTFE_out/${OUTPUT_NAME}_peaks.bed -b src/anno/Coding-5UTR.bed -v > src/anno/peaks_nonClass1.bed
intersectBed -s -wa -wb -a src/anno/peaks_nonClass1.bed -b src/anno/Coding-up_trimmed.bed | awk -F "\t" '{print($4,$10)}' \
|awk -F "|" '{if(a[$1])a[$1]=a[$1]";"$2; else a[$1]=$2;}END{for (i in a)print i, a[i];}' OFS="\t" \
|awk -F " " '{print $1"\t"$2","$3}' \
|awk -F "\t" '{if(a[$1])a[$1]=a[$1]":"$2; else a[$1]=$2;}END{for (i in a)print i, a[i];}' OFS="\t" \
|awk -v 'OFS=\t' '{print $1,$2,"Coding_upstream"}' | sort -k 1,1 > src/anno/peaks_class2.txt

# Coding CDS
intersectBed -s -wa -wb -a src/anno/peaks_nonClass1.bed -b src/anno/Coding-up_trimmed.bed -v > src/anno/peaks_nonClass2.bed
intersectBed -s -wa -wb -a src/anno/peaks_nonClass2.bed -b src/anno/Coding-CDS.bed | awk -F "\t" '{print($4,$10)}' \
|awk -F "|" '{if(a[$1])a[$1]=a[$1]";"$2; else a[$1]=$2;}END{for (i in a)print i, a[i];}' OFS="\t" \
|awk -F " " '{print $1"\t"$2","$3}' \
|awk -F "\t" '{if(a[$1])a[$1]=a[$1]":"$2; else a[$1]=$2;}END{for (i in a)print i, a[i];}' OFS="\t" \
|awk -v 'OFS=\t' '{print $1,$2,"Coding_CDS"}' | sort -k 1,1 > src/anno/peaks_class3.txt

# Coding 3'UTR 
intersectBed -s -wa -wb -a src/anno/peaks_nonClass2.bed -b src/anno/Coding-CDS.bed -v > src/anno/peaks_nonClass3.bed
intersectBed -s -wa -wb -a src/anno/peaks_nonClass3.bed -b src/anno/Coding-3UTR.bed | awk -F "\t" '{print($4,$10)}' \
|awk -F "|" '{if(a[$1])a[$1]=a[$1]";"$2; else a[$1]=$2;}END{for (i in a)print i, a[i];}' OFS="\t" \
|awk -F " " '{print $1"\t"$2","$3}' \
|awk -F "\t" '{if(a[$1])a[$1]=a[$1]":"$2; else a[$1]=$2;}END{for (i in a)print i, a[i];}' OFS="\t" \
|awk -v 'OFS=\t' '{print $1,$2,"Coding_3UTR"}' | sort -k 1,1 > src/anno/peaks_class4.txt

# Noncoding 1st-exon
intersectBed -s -wa -wb -a src/anno/peaks_nonClass3.bed -b src/anno/Coding-3UTR.bed -v > src/anno/peaks_nonClass4.bed
intersectBed -s -wa -wb -a src/anno/peaks_nonClass4.bed -b src/anno/NC-1stexon.bed | awk -F "\t" '{print($4,$10)}' \
|awk -F "|" '{if(a[$1])a[$1]=a[$1]";"$2; else a[$1]=$2;}END{for (i in a)print i, a[i];}' OFS="\t" \
|awk -F " " '{print $1"\t"$2","$3}' \
|awk -F "\t" '{if(a[$1])a[$1]=a[$1]":"$2; else a[$1]=$2;}END{for (i in a)print i, a[i];}' OFS="\t" \
|awk -v 'OFS=\t' '{print $1,$2,"Noncoding_1st-exon"}' | sort -k 1,1 > src/anno/peaks_class5.txt

# Noncoding upstream
intersectBed -s -wa -wb -a src/anno/peaks_nonClass4.bed -b src/anno/NC-1stexon.bed -v > src/anno/peaks_nonClass5.bed
intersectBed -s -wa -wb -a src/anno/peaks_nonClass5.bed -b src/anno/NC-up_trimmed.bed | awk -F "\t" '{print($4,$10)}' \
|awk -F "|" '{if(a[$1])a[$1]=a[$1]";"$2; else a[$1]=$2;}END{for (i in a)print i, a[i];}' OFS="\t" \
|awk -F " " '{print $1"\t"$2","$3}' \
|awk -F "\t" '{if(a[$1])a[$1]=a[$1]":"$2; else a[$1]=$2;}END{for (i in a)print i, a[i];}' OFS="\t" \
|awk -v 'OFS=\t' '{print $1,$2,"Noncoding_upstream"}' | sort -k 1,1 > src/anno/peaks_class6.txt

# Noncoding other exon
intersectBed -s -wa -wb -a src/anno/peaks_nonClass5.bed -b src/anno/NC-up_trimmed.bed -v > src/anno/peaks_nonClass6.bed
intersectBed -s -wa -wb -a src/anno/peaks_nonClass6.bed -b src/anno/NC-exon.bed | awk -F "\t" '{print($4,$10)}' \
|awk -F "|" '{if(a[$1])a[$1]=a[$1]";"$2; else a[$1]=$2;}END{for (i in a)print i, a[i];}' OFS="\t" \
|awk -F " " '{print $1"\t"$2","$3}' \
|awk -F "\t" '{if(a[$1])a[$1]=a[$1]":"$2; else a[$1]=$2;}END{for (i in a)print i, a[i];}' OFS="\t" \
|awk -v 'OFS=\t' '{print $1,$2,"Noncoding_other-exon"}' | sort -k 1,1 > src/anno/peaks_class7.txt

# Intron
intersectBed -s -wa -wb -a src/anno/peaks_nonClass6.bed -b src/anno/NC-exon.bed -v > src/anno/peaks_nonClass7.bed
intersectBed -s -wa -wb -a src/anno/peaks_nonClass7.bed -b src/anno/Intron.bed | awk -F "\t" '{print($4,$10)}' \
|awk -F "|" '{if(a[$1])a[$1]=a[$1]";"$2; else a[$1]=$2;}END{for (i in a)print i, a[i];}' OFS="\t" \
|awk -F " " '{print $1"\t"$2","$3}' \
|awk -F "\t" '{if(a[$1])a[$1]=a[$1]":"$2; else a[$1]=$2;}END{for (i in a)print i, a[i];}' OFS="\t" \
|awk -v 'OFS=\t' '{print $1,$2,"Intron"}' | sort -k 1,1 > src/anno/peaks_class8.txt

# Unannotated
intersectBed -s -wa -wb -a src/anno/peaks_nonClass7.bed -b src/anno/Intron.bed -v > src/anno/peaks_nonClass8.bed
cat src/anno/peaks_nonClass8.bed | awk -v 'OFS=\t' '{print($4,$1":"$3";"$6,"Unannotated")}' > src/anno/peaks_class9.txt

for i in {1..9}; do
cat src/anno/peaks_class${i}.txt 
done | sort -k 1,1 > byTFE_out/${OUTPUT_NAME}_annotation.txt

# Peak annotation
join -1 4 -2 4 -t "$(printf '\011')" <(sort -k 4,4 byTFE_tmp/${OUTPUT_NAME}_TFE-regions.bed) <(sort -k 4,4 byTFE_out/${OUTPUT_NAME}_peaks.bed)  \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$9,$6}' > byTFE_tmp/${OUTPUT_NAME}_TFE-region-peak.txt
join -1 1 -2 1 -t "$(printf '\011')" <(sort -k 1,1 byTFE_out/${OUTPUT_NAME}_annotation.txt) <(sort -k 1,1 byTFE_tmp/${OUTPUT_NAME}_TFE-region-peak.txt)  \
> byTFE_tmp/${OUTPUT_NAME}_TFE-region-peak-anno.txt
join -1 1 -2 1 -t "$(printf '\011')" <(echo -e "Geneid""\t""Gene""\t""Annotation""\t""Chr""\t""Start""\t""End""\t""Peak""\t""Strand" \
| cat - <(sort -k 1,1 byTFE_tmp/${OUTPUT_NAME}_TFE-region-peak-anno.txt)) <(cat byTFE_tmp/${OUTPUT_NAME}_byTFE-counts.txt | sed -e '1d'  \
| awk 'NR<2{print $0;next}{print $0| "sort -k 1,1"}') | cut -f-8,13- | awk 'NR<2{print $0;next}{print $0| "sort -k4,4 -k5,5n -k8,8"}'  \
| sed -e "1 s/Geneid/TFE/g" | sed -e "1 s/byTFE_tmp\/class\///g" | sed -e "1 s/.output.bam//g" | sed -e "1 s/\//\|/g"\
> byTFE_out/${OUTPUT_NAME}_byTFE-counts_annotation.txt

rm byTFE_tmp/${OUTPUT_NAME}_TFE-region-peak.txt &&rm byTFE_tmp/${OUTPUT_NAME}_TFE-region-peak-anno.txt
