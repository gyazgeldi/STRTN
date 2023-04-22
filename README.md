# STRTN-NextSeq Analysis Pipeline

This is a pipeline for the analysis of STRT-N RNA-sequencing outputs from NextSeq. There are scripts that can be used both on a laptop and on the CSC (IT Center for Science) platform. There are mainly five scripts: STRTN.sh (main pipeline/gene-based analysis), STRTN-TFE.sh (TFE-based analysis), fastq-fastQC.sh, STRTN-UCSC-Allas.sh and STRTN-Seurat.sh. There are also codes that these scripts can run on the CSC platform. This pipeline is based on STRT2 pipeline that is developed by Masahito Yoshihara ([Ezer et al., 2021](https://www.sciencedirect.com/science/article/pii/S2666166721007012?via%3Dihub)). This is a fork study, therefore you see the scripts of the previous study that start with STRT2 and are able to run in UPPMAX platform as well.

NOTE: These requirements are required for this pipeline. Sequence data processing, visualization on UCSC and visualization using Seurat require about 150GB, 15GB, 50MB and 150MB of memory depending on genome size and raw data size. The installation of conda packages and pipeline, and running the STRT-N pipeline should be run on Linux terminal.


## Installation
```
git clone https://github.com/gyazgeldi/STRTN.git
```
## Dependencies
You can see the dependencies for each analysis.

For `STRTN.sh` ([The main pipeline with visualization](https://github.com/gyazgeldi/STRTN/blob/master/Visualization-in-UCSC-README.md))
- [Picard](https://broadinstitute.github.io/picard/)
- [HISAT2](https://daehwankimlab.github.io/hisat2/)
- [SAMtools](http://samtools.sourceforge.net/)
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
- [Subread](http://subread.sourceforge.net/)
- [Seqtk](https://github.com/lh3/seqtk)
- [UCSC-bedGraphToBigWig](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/ucsc-bedgraphtobigwig)

For `STRTN-Seurat.sh` ([The visualization of results using Seurat package](https://github.com/gyazgeldi/STRTN/blob/master/STRTN-Seurat-README.md))
- [stringr](https://stringr.tidyverse.org/)
- [dplyr](https://dplyr.tidyverse.org/)
- [ggplot2](https://ggplot2.tidyverse.org/)
- [cowplot](https://www.rdocumentation.org/packages/cowplot/versions/1.1.1)
- [ggbeeswarm](https://github.com/eclarke/ggbeeswarm)
- [forcats](https://forcats.tidyverse.org/)
- [Seurat](https://satijalab.org/seurat/)

For `fastq-fastQC.sh` 
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC](https://multiqc.info/)

For `STRT-TFE.sh` ([TFE-based analysis](https://github.com/my0916/STRT2/blob/master/TFE-README.md))
- [StringTie](https://ccb.jhu.edu/software/stringtie/)

The conda environment is provided as `STRTN-env.yml`. The environment can be created with the followings:
```
conda env create -n STRTN-test -f STRTN-env.yml
conda activate STRTN-test
```

For [CSC](https://www.csc.fi/), these software are available through the `module` command in the scripts (`STRTN-CSC.sh`, `STRTN-Seurat.sh`, `STRTN-UCSC-Allas.sh`, `STRTN-TFE-CSC.sh`, and `fastq-fastqc-CSC.sh` as well as `STRTN-Indexes-Dictionary-CSC.sh`).

## Requirements
These required files are for the main pipeline. For other optional analyses, please visit [additional documentation](https://github.com/gyazgeldi/STRTN#additional-documentation).
- Illumina BaseCalls files (.bcl). The number of lanes is determined based on the number of directories in the basecalls directory. Here is an example of 4 lanes: 
```
  ├── L001
  ├── L002
  ├── L003
  └── L004
```
- HISAT2 index built with a reference genome, (ribosomal DNA), and ERCC spike-ins 
  - See also [How to build HISAT2 index](https://github.com/gyazgeldi/STRTN#how-to-build-hisat2-index-in-csc).
  - The HISAT2 index directory should include the followings:
```
    ├── [basename].1.ht2
    ├── [basename].2.ht2
    ├── [basename].3.ht2
    ├── [basename].4.ht2
    ├── [basename].5.ht2
    ├── [basename].6.ht2
    ├── [basename].7.ht2
    ├── [basename].8.ht2
    ├── [basename].fasta
    └── [basename].dict
```
- Source files (in `src` directory)
  - `barcode.txt` : Barcode sequence with barcode name (1–48). __Please modify if you used different (number of) barcodes.__
  - `ERCC.bed` : 5'-end 50 nt region of ERCC spike-ins ([SRM2374](https://www-s.nist.gov/srmors/view_detail.cfm?srm=2374)) for annotation and quality check.
  - `Example-BarcodesStages` : Sample explanation for data reduction and visualization using PCA, UMAP and violin plots.

## Usage:
These usage are for the main pipeline. For other optional analyses, please visit [additional documentation](https://github.com/gyazgeldi/STRTN#additional-documentation).

For general users:
```
./STRTN.sh -o {OUTPUT_NAME} -g {GENOME_VALUE} -a {ANNO_VALUE} -b {BaseCallsDir_PATH} -i {Index_PATH} -w {WorkingDir_PATH} -c {center_VALUE} -r {run_VALUE} -s {READ_STRUCTURE}    
```
For CSC users:
```
sbatch -A project_2005262 ./STRTN-CSC.sh -o {OUTPUT_NAME} -g {GENOME_VALUE} -a {ANNO_VALUE} -b {BaseCallsDir_PATH} -i {Index_PATH} -w {WorkingDir_PATH} -c {center_VALUE} -r {run_VALUE} -s {READ_STRUCTURE}   
```

## Example usage
These example usage are for the main pipeline. For other optional analyses, please visit [additional documentation](https://github.com/gyazgeldi/STRTN#additional-documentation).

For general users:
```
./STRTN.sh -o STRTN_MOUSE_LIB -g mm39 -a wgEncodeGencodeBasicVM31 -b /mnt/c/Users/gamyaz/STRTN-Pipeline/Data/Intensities/BaseCalls -i /mnt/c/Users/gamyaz/STRTN-Pipeline/mouse_index/mouse_reference -w /mnt/c/Users/gamyaz/STRTN-Pipeline -p /mnt/c/Users/gamyaz/Downloads/ENTER/pkgs/picard-2.27.4-hdfd78af_0/share/picard-2.27.4-0 -c FUGU -r RUNBARCODE -s 8M3S75T6B
```
For CSC users:
```
sbatch -A project_2005262 ./STRTN-CSC.sh -o STRTN_MOUSE_LIB -g mm39 -a wgEncodeGencodeBasicVM31 -b /scratch/project_2005262/Data/Intensities/BaseCalls -i /scratch/project_2005262/mouse_index/mouse_reference -w /scratch/project_2005262 -c FUGU -r RUNBARCODE -s 8M3S75T6B
```

## Parameters
These parameters are for the main pipeline. For other optional analyses, please visit [additional documentation](https://github.com/gyazgeldi/STRTN#additional-documentation).

- __Mandatory__

   | Name | Description |
   | :--- | :--- |
   | `-g, --genome` | Reference genome. Choose one hg19/hg38/mm9/mm10/mm39/canFam3/canFam6/bosTau9. |
   | `-b, --basecalls` | /PATH/to/the Illumina basecalls directory.|
   | `-i, --index` | /PATH/to/the directory and basename of the HISAT2 index for the reference genome. |
   | `-w, --working` | /PATH/to/the working directory. | 
   | `-p, --picardhome` | /PATH/to/the picard.jar. | 

- __Optional__

   | Name&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|Default value|Description|
   | :--- | :--- | :--- |
   | `-o, --out` | OUTPUT | Output file name.|
   | `-a, --annotation` | ref | Gene annotation for QC and counting. <br> Choose from `ref`(RefSeq)/`ens`(Ensembl)/`kg`(UCSC KnownGenes), or directly input the Gencode annotation file name (eg. `wgEncodeGencodeBasicVM31`) for Gencode. <br>Note that some annotations are unavailable in some cases. Please find the details below.
   | `-c, --center ` | CENTER | The name of the sequencing center that produced the reads.<br>Required for the the Picard IlluminaBasecallsToSam program.|
   | `-r, --run` | RUNBARCODE | The barcode of the run. Prefixed to read names.<br>Required for the the Picard IlluminaBasecallsToSam program.|
   | `-s, --structure` | 8M3S75T6B | Read structure.<br>Required for the the Picard IlluminaBasecallsToSam program.<br>Details are described [here](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.4.0/picard_illumina_IlluminaBasecallsToSam.php#--READ_STRUCTURE).|
   | `-d, --dta` | | Add `-d, --dta` (downstream-transcriptome-assembly) if you plan to perform [TFE-based analysis](https://github.com/my0916/STRT2/blob/master/TFE-README.md).<br>Please note that this leads to fewer alignments with short-anchors.|
   | `-h, --help`| | Show usage.|
   | `-v, --version`| | Show version.|
   
   - `-a, --annotation` availability as of Nov 2022:

    | | RefSeq (ref) | Ensembl (ens) | KnownGenes (kg) | Gencode |
    | :---: | :---: | :---: | :---: | :---: |
    | hg19 (human) | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |
    | hg38 (human) | :heavy_check_mark: | NA | :heavy_check_mark: | :heavy_check_mark: |
    | mm9 (mouse) | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | NA |
    | mm10 (mouse) | :heavy_check_mark: | NA | :heavy_check_mark: | :heavy_check_mark: |
    | mm39 (mouse) | :heavy_check_mark: | NA | :heavy_check_mark: | :heavy_check_mark: |
    | canFam3 (dog) | :heavy_check_mark: | :heavy_check_mark: | NA | NA |
    | canFam6 (dog) | :heavy_check_mark: | NA | NA | NA |
    | bosTau9 (bovine) | :heavy_check_mark: | :heavy_check_mark: | NA | NA |
 
## Outputs
Outputs are provided in `out` directory. These outputs are obtained from the main pipeline. For other optional analyses, please visit [additional documentation](https://github.com/gyazgeldi/STRTN#additional-documentation). Unaligned BAM files generated with Picard IlluminaBasecallsToSam program are found in `tmp/Unaligned_bam`.

- __`OUTPUT`-QC.txt__ <br>
Quality check report for all samples.

   | Column |Value|
   | ------------- | ------------- |
   |`Barcode` | Sample name. `OUTPUT` with numbers|
   |`Qualified_reads` | Primary aligned read count|	
   |`Total_reads`| Read count without redundant (duplicate) reads|
   |`Redundancy` | Qualified reads / Total reads| 
   |`Mapped_reads` | Mapped read count (Total reads without unmapped reads)|
   |`Mapped_rate` | Mapped reads / Total reads|  
   |`Spikein_reads` | Read count mapped to ERCC spike-ins|
   |`Spikein-5end_reads` | Read count mapped to the 5'-end 50 nt region of ERCC spike-ins|
   |`Spikein-5end_rate` | Spikein-5end reads / Spikein reads|
   |`Coding_reads` | Read count aligned within any exon or the 500 bp upstream of coding genes|
   |`Coding-5end_reads` | Read count aligned the 5′-UTR or 500 bp upstream of coding genes| 
   |`Coding-5end_rate` | Coding-5end reads / Coding reads|

- __`OUTPUT`-QC-plots.pdf__ <br>
Quality check report by boxplots.
`Mapped_reads`, `Mapped_rate`, `Spikein_reads`, `Mapped / Spikein`, `Spikein-5end_rate`, and `Coding-5end_rate` are shown for all samples. Barcode numbers of outlier samples are marked with red characters.<br> Please consider these outlier samples for the further downstream analysis.

- __`OUTPUT`\_byGene-counts.txt__ <br>
Read count table output from featureCounts.
https://bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf

- __`OUTPUT`\_byGene-counts.txt.summary__ <br>
Filtering summary from featureCounts.
https://bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf

- __Output_bam__ <br>
Resulting BAM files including unmapped, non-primary aligned, and duplicated (marked) reads.

- __Output_bai__ <br>
Index files (.bai) of the resulting BAM files in the `Output_bam` directory.

- __ExtractIlluminaBarcodes_Metrics__ <br>
Metrics file produced by the Picard ExtractIlluminaBarcodes program. The number of matches/mismatches between the barcode reads and the actual barcodes is shown per lane.
https://gatk.broadinstitute.org/hc/en-us/articles/360037426491-ExtractIlluminaBarcodes-Picard-

- __HISAT2_Metrics__ <br>
Alignment summary of samples from each lane produced by the HISAT2 program. 
https://daehwankimlab.github.io/hisat2/manual/

- __MarkDuplicates_Metrics__ <br>
Metrics file indicating the numbers of duplicates produced by the Picard MarkDuplicates program.
https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-

- __`OUTPUT`.output.bam__ <br>
BAM files containing reads except for duplicate and non-primary reads.

- __`OUTPUT`.minus.bw__ and __`OUTPUT`.plus.bw__ <br>
BigWig files for each strands of each sample.

- __coding_5end.bb__ <br>
BigBed file for coding-5'end annotation file.

## How to build HISAT2 index in CSC
Here is the case for house mouse genome (mm39). The genome indexing step requires big memory and it might not be possible to carry out it on a laptop. Indexes and dictionary was prepared in CSC, see commands at STRTN-Indexes-Dictionary-CSC.sh. The built indexes can be accessed in https://doi.org/10.5281/zenodo.7457660. 
### 1. Create conda environment folder file to install the required packages, install and add the bin directory to the path.
```
mkdir STRTN-env
conda-containerize new --prefix STRTN-env STRTN-env.yml
export PATH="<install_dir>/STRTN-env/bin:$PATH"
```
### 2. Load the required module.
```
module load tykky
export PATH="<install_dir>/STRTN-env/bin:$PATH"

module load r-env
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi
echo "TMPDIR=${WorkingDir_PATH}" >> ~/.Renviron
```
### 3. Obtain the genome sequences of reference and ERCC spike-ins. 
You may add the ribosomal DNA repetitive unit for human (U13369) and mouse (BK000964).
```
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
unpigz -c mm39.fa.gz | ruby -ne '$ok = $_ !~ /^>chrUn_/ if $_ =~ /^>/; puts $_ if $ok' > mouse_reference.fasta
wget https://www-s.nist.gov/srmors/certificates/documents/SRM2374_putative_T7_products_NoPolyA_v2.FASTA
cat SRM2374_putative_T7_products_NoPolyA_v2.FASTA >> mouse_reference.fasta
```
### 4. Extract splice sites and exons from a GTF file. Here we used wgEncodeGencodeBasicVM31 as the annotation file.
```
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/wgEncodeGencodeBasicVM31.txt.gz
unpigz -c wgEncodeGencodeBasicVM31.txt.gz | hisat2_extract_splice_sites.py - | grep -v ^chrUn > splice_sites.txt
unpigz -c wgEncodeGencodeBasicVM31.txt.gz | hisat2_extract_exons.py - | grep -v ^chrUn > exons.txt
```
You may additionally perform `hisat2_extract_snps_haplotypes_UCSC.py` to extract SNPs and haplotypes from a dbSNP file for human and mouse.
### 5. Build the HISAT2 index.
```
hisat2-build mouse_reference.fasta --ss splice_sites.txt --exon exons.txt mouse_index/mouse_reference
```
This outputs a set of files with suffixes. Here, `mouse_reference.1.ht2`, `mouse_reference.2.ht2`, ..., `mouse_reference.8.ht2` are generated.<br>In this case, `mouse_reference` is the basename used for `-i, --index`.
### 6. Create the sequence dictionary for the reference and Spike-in sequences.
```
picard CreateSequenceDictionary R=mouse_reference.fasta O=mouse_reference.dict
```
This is required for the Picard MergeBamAlignment program. Note that the original FASTA file (`mouse_reference.fasta` here) is also required.
### 7. Put the genome indexes, genome fasta file, sequence dictionary to same folder.
```
mv mouse_reference.dict mouse_reference
mv mouse_reference.fasta mouse_reference
```
## Additional documentation
- [STRTN-TFE.sh](https://github.com/gyazgeldi/STRTN/blob/master/STRTN-TFE.sh): This script finds only reads within 5'-UTR or proximal upstream of protein-coding and non-coding genes. The details are in [READ.ME file](https://github.com/gyazgeldi/STRTN/blob/master/STRTN-TFE-README.md).
- [STRTN-UCSC-Allas.sh](https://github.com/gyazgeldi/STRTN/blob/master/STRTN-UCSC-Allas.sh): This script uploads BAM, BigWig, BED files to UCSC-Allas storage service and creates an accessible link to visualize in UCSC genome browser tool. The details are in [READ.ME file](https://github.com/gyazgeldi/STRTN/blob/master/Visualization-in-UCSC-README.md).
- [STRTN-Seurat.sh](https://github.com/gyazgeldi/STRTN/blob/master/STRTN-Seurat.sh): This script performs scRNA-seq analysis using R-Seurat package and creates PCA, UMAP, violin plots. The details are in [READ.ME file](https://github.com/gyazgeldi/STRTN/blob/master/STRTN-Seurat-README.md).
- [fastq-fastQC.sh](https://github.com/gyazgeldi/STRTN/blob/master/fastq-fastQC.sh): After running the [main pipeline](https://github.com/gyazgeldi/STRTN/blob/master/STRTN.sh), you can generate fastq files for each sample from the output BAM files in the `fastq` directory. These fastq files (without duplicated reads) can be submitted to public sequence databases.<br>
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) files are also generated for each fastq file in the `fastqc` directory. Based on the FastQC results, [MultiQC](https://multiqc.info/) report (__MultiQC_report.html__) is generated.

## Flowchart
![image](https://user-images.githubusercontent.com/101990822/233770065-b217d1c8-f608-47d8-8b5a-619da678a430.png)

