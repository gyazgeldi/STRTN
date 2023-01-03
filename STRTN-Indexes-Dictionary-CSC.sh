#!/bin/bash
#SBATCH --job-name=mouse-indexes-dic
#SBATCH --account=project_2005262
#SBATCH --time=05:00:00
#SBATCH --partition=small
#SBATCH --mem=200GB
#SBATCH --cpus-per-task=8

module load tykky
export PATH="/scratch/project_2005262/STRT2-env/bin:$PATH"

# Creating a directory for indexes
mkdir mouse_index

# Obtain the genome sequences of reference and ERCC spike-ins
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
unpigz -c mm39.fa.gz | ruby -ne '$ok = $_ !~ /^>chrUn_/ if $_ =~ /^>/; puts $_ if $ok' > mouse_reference.fasta

wget https://www-s.nist.gov/srmors/certificates/documents/SRM2374_putative_T7_products_NoPolyA_v2.FASTA
cat SRM2374_putative_T7_products_NoPolyA_v2.FASTA >> mouse_reference.fasta

# Downloading annotation file
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/gencode.vM31.annotation.gtf.gz

# Extract splice sites and exons from annotation .gtf file
unpigz -c gencode.vM31.annotation.gtf.gz | hisat2_extract_splice_sites.py - | grep -v ^chrUn > splice_sites.txt
unpigz -c gencode.vM31.annotation.gtf.gz | hisat2_extract_exons.py - | grep -v ^chrUn > exons.txt

# Build the genome HISAT2 indexes
hisat2-build mouse_reference.fasta --ss splice_sites.txt --exon exons.txt mouse_index/mouse_reference

# Creating dictionary
picard CreateSequenceDictionary R=mouse_reference.fasta O=mouse_reference.dict

# Should be in the same folder
mv mouse_reference.dict mouse_index
mv mouse_reference.fasta mouse_index

# Removing files so it doesn't take up space
rm mm39.fa
rm mm39.fa.gz
rm gencode.vM31.annotation.gtf
rm gencode.vM31.annotation.gtf.gz
rm SRM2374_putative_T7_products_NoPolyA_v2.FASTA
rm exons.txt
rm splice_sites.txt
