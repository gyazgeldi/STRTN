# TFE-based analysis

Transcript far 5'-end (TFE) analysis is developed in [Töhönen et al. 2015](https://doi.org/10.1038/ncomms9207) and is modified in [Ezer et al. 2021](https://doi.org/10.1016/j.xpro.2021.100995). TFEs are defined as the first exon (5'-end region) of assembled STRT reads (transcripts) mapped to genome. Here is detail of [TFE analysis pipeline](https://github.com/my0916/STRT2) (by Yoshihara M.) and visualization in UCSC.

## Dependencies
- [StringTie](https://ccb.jhu.edu/software/stringtie/)
- [SAMtools](http://samtools.sourceforge.net/)
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
- [Subread](http://subread.sourceforge.net/)

## Requirements
- Source files (in `src` directory)
  - `TFEclass.txt` : Sample classification with barcode name (1-48) which is used for the transcript assembly. Please set __`NA`__ for those samples that are not used for further analysis (e.g. negative controls or outlier samples). In the default settings, all samples are classfied as the same class.
   #### Example
     |     |     |
     | :-: | :-: |
     | 1 | classA | 
     | 2 | classB | 
     | 3 | classC | 
     | 4 | classA | 
     | 5 | classB | 
     | 6 | classC | 
     | 7 | classA | 
     | 8 | classB | 
     | 9 | NA | 
    
  In this case, transcript assemblies are performed using i) samples 1, 4, 7 for classA, ii) 2, 5, 8 for classB, and iii) 3, 6 for classC, respectively. Then first exons are collected from all these 3 classes and named as TFEs. Sample 9 is not used for the analysis. 
  - `refGene.txt`, `knowngene-names.txt`,  `ens-genes.txt`, or `Gencode.txt`  : Prepared within the STRTN NextSeq-pipeline (STRTN.sh or STRTN-CSC.sh), which is used for the annotation of TFE peaks.
  
## Example usage
```
./STRTN-TFE.sh  
```
For CSC:
```
sbatch -A /scratch/<project_id>/ -p core -n 8 -t 24:00:00 ./STRTN-TFE-CSC.sh -w /scratch/<project_id>
```

## Parameters
- __Optional__

   | Name&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|Default value|Description|
   | :--- | :--- | :--- |
   | `-c, --coverage` | 5 | Minimum read coverage allowed for the predicted transcripts used in StringTie.|
   | `-l, --length` | 75 | Minimum length allowed for the predicted transcripts used in StringTie.|
   | `-h, --help`| | Show usage.|
   | `-v, --version`| | Show version.|

## Outputs
Outputs are provided in `byTFE_out` directory.

- __`OUTPUT`\_byTFE-counts_annotation.txt__ <br>
Read count table output from featureCounts with genomic annotations. Note that the order of samples (columns) are based on `TFEclass.txt` (different from `OUTPUT`\_byGene-counts.txt).

- __`OUTPUT`\_annotation.txt__ <br>
TFE annotation output containing gene name and transcript ID or chromosome position, strand and annotation detail.

- __`OUTPUT`\_peaks.bed__ <br>
TFE peaks output containing chromosome, chromosome positions of peak, TFE name, score and strand.

## Visualization in UCSC
To visualize in UCSC, please follow following steps:

1. Merge with TFE region and TFE peak files based on TFE name column.
```
join -1 4 -2 4 -t "$(printf '\011')" <(sort -k 4,4 byTFE_tmp/${OUTPUT_NAME}_TFE-regions.bed) <(sort -k 4,4 byTFE_out/${OUTPUT_NAME}_peaks.bed) | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1,$5,$6,$8,$9,$10}' > ${OUTPUT_NAME}_TFE-region-peak.txt 
```
2. Remove rows belong to spike-ins.
3. Upload to `OUTPUT`\_TFE-region-peak.txt to UCSC.
