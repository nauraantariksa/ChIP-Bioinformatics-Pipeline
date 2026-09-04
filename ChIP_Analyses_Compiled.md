# Table of Contents
- [Data analyses](#data-analyses)
  - [List of tools/software](#list-of-toolssoftware)
  - [Adaptor Trimming](#adaptor-trimming)
  - [Genome Download and Indexing](#genome-download-and-indexing)
  - [Alignment of BG4-ChIP data](#alignment-of-bg4-chip-data)
  - [BAM QC](#bam-qc)
  - [Peak Calling with MACS2](#peak-calling-with-macs2)
  - [Peak analyses (fold change, enrichment)](#peak-analyses-fold-change-enrichment)
  - [Differential BG4 binding analyses (CSB+ vs CSB-)](#differential-bg4-binding-analyses-csb-vs-csb-)

# Data analyses
This GitHub page is a collection of scripts and general code used for the analysis of ChIP datasets. Differential analyses were given with a comparison between CS1AN CSB+ and CSB- cell lines as an example.

## List of tools/software
Make sure these dependencies are installed in your system:
- anaconda3
- python (version 3.9)
- fastp v0.23.4
- FastQC v0.11.9
- Bowtie2 v2.5.1
- SAMtools v1.22.1
- Picard v2.25.1
- BEDTools v2.31.1

In anaconda3, make separate environments for the following Python suites used for the data analyses:
- deepTools v3.5.6
- MACS2 v2.2.9.1

For other analyses, install or clone Github repositories:
- R
- HOMER v5.1

For visualization:
- IGV v2.19.8

## Adaptor Trimming
Raw fastq files were trimmed to remove adaptor reads before alignment.
```bash
echo "FILENAME_1_R1 FILENAME_1_R2 FILENAME_2_R1 FILENAME_2_R2" > FILENAME.txt

#PBS -lwalltime=08:00:00
#PBS -lselect=1:ncpus=250:mem=900g

eval "$(/rds/general/user/nfa23/home/anaconda3/bin/conda shell.bash hook)"

for name in $(cat $PBS_O_WORKDIR/FILENAME.txt)
do
base_name=$(echo $name | sed 's/_R[12]$//')
/rds/general/user/nfa23/home/anaconda3/bin/fastp -i $PBS_O_WORKDIR/${base_name}_R1.fq.gz -I $PBS_O_WORKDIR/${base_name}_R2.fq.gz -o trimmed_${base_name}_R1.fq.gz -O trimmed_${base_name}_R2.fq.gz --detect_adapter_for_pe -l 20 -j ${base_name}.fastp.json -h ${base_name}.fastp.html
done

cp * $PBS_O_WORKDIR
```

To check the quality of trimming, quality control was performed with FastQC
```bash
#PBS -lwalltime=08:00:00
#PBS -lselect=1:ncpus=1:mem=900gb

module load fastqc

for file in $PBS_O_WORKDIR/*.fastq
do
fastqc $file -d $TMPDIR -o .
done

cp * $PBS_O_WORKDIR
```

## Genome Download and Indexing
All files were aligned to the human reference genome hg38 (UCSC). Here, the link to the genome of reference is given.
```bash
#PBS -lwalltime=08:00:00
#PBS -lselect=1:ncpus=1:mem=900gb

wget --no-verbose --directory-prefix hg38_genome "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
gunzip hg38.fa.gz

cp * $PBS_O_WORKDIR
```
The genome was indexed using Bowtie2
```bash
#PBS -lwalltime=08:00:00
#PBS -lselect=1:ncpus=1:mem=900gb

eval "$(/rds/general/user/nfa23/home/anaconda3/bin/conda shell.bash hook)"
module load Bowtie2/2.5.1-GCC-12.3.0

bowtie2-build $PBS_O_WORKDIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa hg38

cp * $PBS_O_WORKDIR 
```

## Alignment of BG4-ChIP data
ChIP data files were aligned using Bowtie2, generating SAM files. The aligned files were then filtered (for unmapped reads, non-primary reads, supplementary alignment, and low-quality read mapping < 10), to generate BAM files.
```bash
#PBS -lwalltime=08:00:00
#PBS -lselect=1:ncpus=250:mem=900gb

eval "$(/rds/general/user/nfa23/home/anaconda3/bin/conda shell.bash hook)"
module load Bowtie2/2.5.1-GCC-12.3.0

read_1="path/to/read_1/trimmed/file.fq"
read_2="path/to/read_2/trimmed/file.fq"

bowtie2 -x /rds/general/project/diantonio/live/tools/genomes/hg38_ucsc/hg38 -1 $read_1 -2 $read_2 -S FILE_NAME_Bio1.sam 2> summary_FILE_NAME_Bio1.txt

cp * $PBS_O_WORKDIR
```
```bash
#PBS -lwalltime=08:00:00
#PBS -lselect=1:ncpus=250:mem=900gb

module load SAMtools/1.18-GCC-12.3.0
module load BEDTools/2.31.0-GCC-12.3.0
module load tools/prod
module load picard/2.25.1-Java-11

sam_dir=/PATH/TO/SAM/FILES
tmp_dir=$EPHEMERAL/tmp_filtering
mkdir -p $tmp_dir

for file in $sam_dir/*.sam
do
name=`basename $file .sam`
samtools view -F 2308 -b -q 10 $sam_dir/$name.sam > $tmp_dir/$name.bam
samtools sort -o $tmp_dir/$name.sort.bam -O 'bam' -T $tmp_dir/temp_$name $tmp_dir/$name.bam
java -jar $EBROOTPICARD/picard.jar MarkDuplicates REMOVE_DUPLICATES=TRUE I=$tmp_dir/$name.sort.bam O=$tmp_dir/$name.dedup.sort.bam M=$tmp_dir/$name.metrics.txt TMP_DIR=$tmp_dir
bedtools intersect -v -abam $tmp_dir/$name.dedup.sort.bam -b '/rds/general/project/diantonio/live/tools/blacklisted_regions/hg38-blacklist.v2.bed' > $tmp_dir/$name.bl.dedup.sort.bam
samtools index $tmp_dir/$name.bl.dedup.sort.bam
rm $tmp_dir/$name.bam $tmp_dir/$name.sort.bam $tmp_dir/$name.dedup.sort.bam
done

cp $tmp_dir/*.bl.dedup.sort.bam* $tmp_dir/*.metrics.txt "$PBS_O_WORKDIR"
```

## BAM QC
To ensure the quality of the ChIP data, several quality control processes were performed at the BAM level. This included correlation, fingerprint plot, and strand cross-correlation analyses. Most of these analyses were performed using deepTools.

### Correlation
BAM files were converted into BIGWIG files, normalized by reads per-kilobase million (RPKM)
```bash
#PBS -lwalltime=08:00:00
#PBS -lselect=1:ncpus=128:mem=900gb

eval "$(/rds/general/user/nfa23/home/anaconda3/bin/conda shell.bash hook)"
bam_dir='/rds/general/project/diantonio/live/Naura/ChIP/CS1AN_CSB_Novogene/Deduplication'
for file in $bam_dir/*.dedup2.bl.dedup.sort.bam
do
name=`basename $file .bl.dedup.sort.bam`
bamCoverage -b $file -o $name.bs5.RPKM.bw --binSize 5 --numberOfProcessors max --normalizeUsing RPKM
done
```
```bash
#PBS -lwalltime=72:00:00
#PBS -lselect=1:ncpus=128:mem=900gb

eval "$(~/anaconda3/bin/conda shell.bash hook)"
conda activate deeptools_env
cd $PBS_O_WORKDIR
multiBigwigSummary bins --bwfiles FILE_1.bw FILE_2.bw FILE_3.bw --binSize 5 --numberOfProcessors max -o FILE_NAME_Matrix.npz
plotCorrelation -in FILE_NAME_Matrix.npz --corMethod spearman --whatToPlot heatmap --colorMap RdYlBu --plotNumbers --removeOutliers -o FILE_NAME_Matrix.png
```

### Fingerprint plot
```bash
#PBS -lwalltime=72:00:00
#PBS -lselect=1:ncpus=128:mem=900gb

eval "$(~/anaconda3/bin/conda shell.bash hook)"
conda activate deeptools_env
plotFingerprint -b FILE_1.bam FILE_2.bam FILE_3.bam --labels FILE_1 FILE_2 FILE_3 --binSize 5 --minMappingQuality 30 --numberOfProcessors 128 -T "Fingerprints of FILE_NAME" --plotFile Fingerprints_FILE_NAME.png
```

### Strand cross-correlation


## Peak Calling with MACS2

## Peak analyses (fold change, enrichment)

## Differential BG4 binding analyses (CSB+ vs CSB-)
