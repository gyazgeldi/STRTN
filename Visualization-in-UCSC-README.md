# Visualization in UCSC

STRT-N library results is viewed in UCSC-tool (Lee et al., 2022) by the main pipeline STRTN.sh. 

## Dependencies
- [Samtools](https://www.htslib.org/)
- [Seqtk](https://github.com/lh3/seqtk)
- [Bedtools](https://bedtools.readthedocs.io/en/latest/)
- [BedGraphToBigWig](https://github.com/ucscGenomeBrowser/kent/releases/tag/v377_base)

## Requirements
- `*.output.bam` : BAM files for each sample (in `out` directory).
- `coding_5end.bed`: Coding 5'end BED file as an annotation file to view 5'ends (in `src` directory). 

## Example usage
For general users:
```
./STRTN.sh -o STRTN_MOUSE_LIB -g mm39 -a wgEncodeGencodeBasicVM30 -b /mnt/c/Users/gamyaz/STRTN-Pipeline/Data/Intensities/BaseCalls -i /mnt/c/Users/gamyaz/STRTN-Pipeline/mouse_index/mouse_reference -w /mnt/c/Users/gamyaz/STRTN-Pipeline -p /mnt/c/Users/gamyaz/Downloads/ENTER/pkgs/picard-2.27.4-hdfd78af_0/share/picard-2.27.4-0 -e gamze.yezgeldi@helsinki.fi -n STRTN-hub-mouse -c FUGU -r RUNBARCODE -s 8M3S75T6B
```
For CSC users:
```
sbatch -A project_2005262 ./STRTN-CSC.sh -o STRTN_MOUSE_LIB -g mm39 -a wgEncodeGencodeBasicVM30 -b /scratch/project_2005262/Data/Intensities/BaseCalls -i /scratch/project_2005262/mouse_index/mouse_reference -w /scratch/project_2005262 -e gamze.yezgeldi@helsinki.fi -n STRTN-hub-mouse -c FUGU -r RUNBARCODE -s 8M3S75T6B
```

## Parameters
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
   | `-a, --annotation` | ref | Gene annotation for QC and counting. <br> Choose from `ref`(RefSeq)/`ens`(Ensembl)/`kg`(UCSC KnownGenes), or directly input the Gencode annotation file name (eg. `wgEncodeGencodeBasicVM30`) for Gencode. <br>Note that some annotations are unavailable in some cases. Please find the details below.
   | `-c, --center ` | CENTER | The name of the sequencing center that produced the reads.<br>Required for the the Picard IlluminaBasecallsToSam program.|
   | `-r, --run` | RUNBARCODE | The barcode of the run. Prefixed to read names.<br>Required for the the Picard IlluminaBasecallsToSam program.|
   | `-s, --structure` | 8M3S74T6B | Read structure.<br>Required for the the Picard IlluminaBasecallsToSam program.<br>Details are described [here](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.4.0/picard_illumina_IlluminaBasecallsToSam.php#--READ_STRUCTURE).|
   | `-e, --email` | EMAIL | Email address for hub file to share details.|
   | `-n, --name` | NAME | Hub name to store in CSC-Allas or hosting.|
   | `-d, --dta` | | Add `-d, --dta` (downstream-transcriptome-assembly) if you plan to perform [TFE-based analysis](https://github.com/my0916/STRT2/blob/master/TFE-README.md).<br>Please note that this leads to fewer alignments with short-anchors.|
   | `-h, --help`| | Show usage.|
   | `-v, --version`| | Show version.|

## Outputs
Outputs are provided in {WorkingPATH} directory.

- __hub.txt__ <br>
Parameters for each tracks. 

- __link__ <br>
The link indicating hub.txt location that can be copied directly to Track Hub in the UCSC genome browser

## Hosting of data
Track hub file and results files must be located in web-accessible locations (https://genome.ucsc.edu/goldenPath/help/hgTrackHubHelp#Hosting). 
Here is an example to store files in CSC-Allas container by STRTN-UCSC-Allas.sh. To store files, please run following commands.

### 1. Load required module.
Please make sure that you accessed allas and are able to use swift protocol. 
```
module load allas
allas-conf <project>
```
### 2. Create new container for visualization files.
```
swift post $NAME_HUB
```
### 3. Upload files for example .bam, .bai and .bw files.
```
for file in *.output.bam
do
    name=$(basename $file .output.bam)
    swift upload $NAME_HUB ${WorkingDir_PATH}/$file
    swift upload $NAME_HUB ${WorkingDir_PATH}/$file.bai
    swift upload $NAME_HUB ${WorkingDir_PATH}/${name}_plus.bw
    swift upload $NAME_HUB ${WorkingDir_PATH}/${name}_minus.bw
    swift post $NAME_HUB --read-acl ".r:*"
done
```
### 4. Upload other annotation file for example coding_5´-end.txt and hub.txt file.
```
swift upload $NAME_HUB ${WorkingDir_PATH}/coding_5end.bb
swift upload $NAME_HUB ${WorkingDir_PATH}/hub.txt
```
### 5. Make accessible and create sharable link. Remove unneccessary files.
```
swift post $NAME_HUB --read-acl ".r:*"
account=$(swift stat | grep "Account: " | awk '{print $2}')
link=https://a3s.fi/swift/v1/$account/$NAME_HUB${WorkingDir_PATH}/hub.txt
echo -e "The sharable link is""\t"$link

rm *.output.bam && rm *.output.bam.bai 
rm *.bw && rm coding_5end.bb
```
