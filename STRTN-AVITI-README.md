# STRT-N AVITI pipeline 

This is an AVITI-compatible version of STRTN preprocessing pipeline for use on CSC Roihu environment. Unlike the original NextSeq pipeline, this version starts directly from AVITI FASTQ files. The FASTQ files are converted to unmapped BAM file (one unmapped BAM file per library) using 'fgbio' and demultiplexed based on the barcode tag ('BC:'). The pipeline assumes 48 samples per one STRT-N library. The remaining preprocessing steps follow the original STRT-N pipeline. Therefore, 'fgbio' was added to 'STRTN-env-AVITI.yml'.  

## Files:
- "STRTN-AVITI-Roihu.sh" (main AVITI preprocessing pipeline)
- "STRTN-env-AVITI.yml" (software environment for the AVITI preprocessing pipeline)
- "src/barcode.txt" (sample barcodes for demultiplexing) 

## Input files:
The pipeline requires two FASTQ files from one AVITI STRT-N library:
- R1 FASTQ: sequencing reads
- I1 FASTQ: sample barcode reads

Default read structures: 
R1: `8M3S74T`
I1: `6B`

## Environment setup on Roihu
Create the environment using Tykky:
```bash
module load tykky
conda-containerize new --mamba \
    --prefix /projappl/<PROJECT>/STRTN-env-AVITI \
    STRTN-env-AVITI.yml
```

## Reference genome
A HISAT2 index built with the corresponding reference genome is required. Please use 'STRTN-Indexes-Dictionary-CSC.sh' as an example.

The required argument '-i' should contain the path to the HISAT2 index and its basename, without .fasta extension (e.g. /path/to/dog_index/dog_reference).

## Annotation
Please check the main STRTN README table for available annotation options.

## Usage
```bash
STRTN-AVITI-Roihu.sh \
    -o {OUTPUT_NAME} \
    -g {GENOME_VALUE} \
    -a {ANNO_VALUE} \
    -R {R1Fastq_PATH} \
    -B {I1Fastq_PATH} \
    -i {Index_PATH} \
    -w {WorkingDir_PATH} \
    -s "{READ_STRUCTURE}"
```

## Output files
The main output files are generated in the 'out/' directory, including:
- sample BAM files
- HISAT2 alignment metrics
- PCR duplicate marking metrics
- QC summary
- gene-level count matrix
- QC plots
