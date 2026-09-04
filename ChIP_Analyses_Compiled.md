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
ChIP peaks were identified using MACS2, using several flags to ensure reproducibility and quality of peaks across replicates.
```bash
eval "$(/rds/general/user/nfa23/home/anaconda3/bin/conda shell.bash hook)"
conda activate macs2_env
macs2 callpeak  -t /path/to/ChIP/BAM/File -c /path/to/ChIP/Input/BAM/File -f BAMPE --keep-dup all -p 0.0001 -n FILE_NAME
```
## Peak analyses (fold change, enrichment)
There are several ways to determine the consensus peaks within a cell line. Here, intersecting 3x biological replicates using multiinter tool of BEDTools is shown. Upon determining consensus peaks, they can be annotated using HOMER, further analyzed using MEME suite, etc.

```bash
module load BEDTools/2.31.1-GCC-14.3.0
bedtools multiinter -i FILE_NAME_Bio1.bed FILE_NAME_Bio2.bed FILE_NAME_Bio2.bed | awk '$4 >= 2' OFS="\t" > FILE_NAME_Multiinter.bed
awk 'BEGIN {OFS="\t"} {print $1, $2, $3}' FILE_NAME_Multiinter.bed > FILE_NAME_Multiinter_Min.bed #the multiinter bed file includes many columns that are not compatible with IGV. This command prints the first 3 columns of a BED file to ensure IGV can open the file.
bedtools merge -i FILE_NAME_Multiinter_Min.bed > FILE_NAME_Multiinter_Min_Merged.bed #multiinter causes peak "fragmentation". This command merges fragmented peaks.

annotatePeaks.pl  FILE_NAME_Multiinter_Min_Merged.bed hg38 > Annotated_FILE_NAME_Multiinter_Min_Merged.bed 2> Report_FILE_NAME_Multiinter_Min_Merged.bed.bed
```
Peak annotation gives the number of peaks in certain genomics regions. This can be plotted to show the distribution of peaks across the genome. To understand if the significance of those peaks in certain genomic loci, fold enrichment analysis was performed.
```bash
bedtools shuffle -excl /rds/general/user/nfa23/projects/diantonio/live/tools/blacklisted_regions/hg38-blacklist.v2.bed -incl /rds/general/user/nfa23/projects/diantonio/live/Naura/CnT/CS1AN_WT/PeakFiles_SEACR/UCSC/g4seq_marsico_hg38.bed -maxTries 1000 -noOverlapping -i FILE_NAME_Multiinter_Min_Merged.bed  -g hg38.genome > FILE_NAME_SharedPeaks_Shuffled_1.bed 2> FILE_NAME_SharedPeaks_Shuffled_1_error.log
```

## Differential BG4 binding analyses (CSB+ vs CSB-)
ChIP differential analysis was adapted from Robert Hansel-Hertsch (2016) and requires cloning of the GitHub repository https://github.com/sblab-bioinformatics/dna-secondary-struct-chrom-lands/tree/master. Here, the differential analysis between CS1AN cells with and without CSB is shown. Note that consensus peak determination and merging were performed differently here.
```bash
~/dna-secondary-struct-chrom-lands/utils/mergePeaks.sh *_peaks.narrowPeak | cut -f1-3 > CS1AN_WTvCSB_BG4_Union.bed
for bam in CS1AN_*.dedup2.bl.dedup.sort.bam
do   
coverageBed -sorted -g genome.full.txt -a CS1AN_WTvCSB_Union.main.bed -b "$bam" > "${bam%.bam}.union.bed"
done

echo 'chrom start end nreads nonzero len frac sample_id' | tr ' ' '\t' > union.counts.bed
~/dna-secondary-struct-chrom-lands/utils/tableCat.py -i CS1AN_*.union.bed -r '.dedup2.bl.dedup.sort.union.bed' >> union.counts.bed
for b in CS1AN_*.dedup2.bl.dedup.sort.bam; do
  echo -e "$(basename $b .dedup2.bl.dedup.sort.bam)\t$(samtools view -c -f 0x2 -F 0x400 $b)"
done > libsizes.txt
```
Upon obtaining union.counts.bed and libsizes.txt files, the data was analysed locally using an R-program below, as adapted from Robert Hansel-Hertsch (2016) with the aid of Claude Sonnet 5. The program below is specifically for CSB+/- project.
```bash
library(data.table)
library(edgeR)

setwd(".")   # directory containing the *.union.bed files

# --- 1. read the tableCat.py long-format table into a count matrix ----------

cnt <- fread("/Users/naura/Downloads/union.counts.bed")
cnt[, locus := paste(chrom, start, end, sep = "_")]

cntct  <- dcast.data.table(cnt, locus ~ sample_id, value.var = "nreads")
counts <- as.matrix(cntct[, -1])
rownames(counts) <- cntct$locus

# dcast sorts columns alphabetically - keep the locus order for later use
loci <- cntct$locus

message(sprintf("Loaded %d regions x %d libraries", nrow(counts), ncol(counts)))
print(colnames(counts))


# --- 2. separate IP from input ---------------------------------------------

is_input <- grepl("Input", colnames(counts))
ip  <- counts[, !is_input, drop = FALSE]
inp <- counts[,  is_input, drop = FALSE]


# --- 3. sum technical replicates within each biological replicate -----------
# Tech1/Tech2 are independent IPs from the same chromatin: they are NOT independent biological observations. Summing keeps integer counts and the negative-binomial mean-variance relationship intact.

bioreps <- c("CS1AN_CSB_BG4_Bio1", "CS1AN_CSB_BG4_Bio3", "CS1AN_CSB_BG4_Bio4",
             "CS1AN_WT_BG4_Bio1",  "CS1AN_WT_BG4_Bio2",  "CS1AN_WT_BG4_Bio3")

ip_bio <- sapply(bioreps, function(b) {
  cols <- grepl(b, colnames(ip), fixed = TRUE)
  stopifnot(sum(cols) > 0)
  rowSums(ip[, cols, drop = FALSE])
})
message("\nTechnical replicates summed. Biological replicates:")
print(colSums(ip_bio))


# --- 4. library sizes = TOTAL mapped reads per BAM --------------------------

libsize_file <- "/Users/naura/Downloads/libsizes.txt"
if (file.exists(libsize_file)) {
  ls_tab <- fread(libsize_file, col.names = c("sample", "total"))
  # sum the technical replicates' library sizes to match the summed counts
  libSize <- sapply(bioreps, function(b)
    sum(ls_tab$total[grepl(b, ls_tab$sample, fixed = TRUE)]))
  stopifnot(all(libSize > 0))
} else {
  warning("libsizes.txt not found - falling back to column sums. ",
          "This will absorb the FRiP difference. Generate the file.")
  libSize <- colSums(ip_bio)
}
print(libSize)


# --- 5. edgeR ---------------------------------------------------------------

group <- factor(c(rep("CSBplus", 3), rep("CSBminus", 3)),
                levels = c("CSBminus", "CSBplus"))   # logFC = CSB+ / CSB-
design <- model.matrix(~ group)

run_edger <- function(norm_method) {
  y <- DGEList(counts = ip_bio, group = group)
  y$samples$lib.size <- libSize
  
  keep <- filterByExpr(y, design)
  message(sprintf("\n[%s] regions retained after filtering: %d of %d",
                  norm_method, sum(keep), length(keep)))
  y <- y[keep, , keep.lib.sizes = TRUE]     # TRUE: keep the BAM-derived sizes
  
  y <- calcNormFactors(y, method = norm_method)
  message(sprintf("[%s] normalisation factors:", norm_method))
  print(round(y$samples$norm.factors, 3))
  
  y   <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)                # QL: better for small n than exactTest
  res <- glmQLFTest(fit, coef = 2)
  
  dt <- data.table(topTags(res, n = Inf)$table, keep.rownames = "locus")
  message(sprintf("[%s] FDR < 0.05: %d up, %d down",
                  norm_method, nrow(dt[FDR < 0.05 & logFC > 0]),
                  nrow(dt[FDR < 0.05 & logFC < 0])))
  list(y = y, dt = dt)
}

res_none <- run_edger("none")
res_tmm  <- run_edger("TMM")


# --- 6. diagnostics ---------------------------------------------------------

pdf("bg4_diagnostics.pdf", width = 10, height = 5)
par(mfrow = c(1, 2), las = 1, bty = "l")
plotMDS(res_none$y, labels = colnames(ip_bio), main = "MDS (lib.size = total)")
plotBCV(res_none$y, main = "Biological CV")
dev.off()

# MA plots, one per normalisation
pal <- colorRampPalette(c("white", "lightblue", "yellow", "red"), space = "Lab")
pdf("bg4_maplots.pdf", width = 10, height = 5, pointsize = 10)
par(mfrow = c(1, 2), las = 1, mgp = c(1.75, 0.5, 0), bty = "l",
    mar = c(3, 3, 3, 0.5))
for (nm in c("none", "TMM")) {
  dt <- if (nm == "none") res_none$dt else res_tmm$dt
  smoothScatter(dt$logCPM, dt$logFC, colramp = pal, col = "blue",
                xlab = "logCPM", ylab = "logFC",
                main = sprintf("BG4 [CSB+ - CSB-], norm = %s", nm))
  lines(loess.smooth(dt$logCPM, dt$logFC, span = 0.1), lwd = 2, col = "grey60")
  abline(h = 0, col = "grey30")
  points(dt$logCPM, dt$logFC,
         col = ifelse(dt$FDR < 0.05, "#FF000080", "transparent"),
         cex = 0.5, pch = ".")
  mtext(side = 3, line = -1.2, adj = 1,
        text = sprintf("FDR<0.05 up: %s", nrow(dt[FDR < 0.05 & logFC > 0])))
  mtext(side = 1, line = -1.2, adj = 1,
        text = sprintf("FDR<0.05 down: %s", nrow(dt[FDR < 0.05 & logFC < 0])))
  grid(col = "grey50")
}
dev.off()


# --- 7. concordance between the two normalisations --------------------------
# The intersection is the robust set. The 'none'-only set is the global
# component and carries the no-spike-in caveat.

sig_none <- res_none$dt[FDR < 0.05, locus]
sig_tmm  <- res_tmm$dt[FDR < 0.05, locus]
message(sprintf("\nSignificant: none = %d, TMM = %d, both = %d",
                length(sig_none), length(sig_tmm),
                length(intersect(sig_none, sig_tmm))))


# --- 8. write out, with coordinates restored --------------------------------

annotate <- function(dt) {
  parts <- tstrsplit(dt$locus, "_", fixed = TRUE)
  dt[, `:=`(chrom = parts[[1]],
            start = as.integer(parts[[2]]),
            end   = as.integer(parts[[3]]))]
  setcolorder(dt, c("chrom", "start", "end"))[]
}

fwrite(annotate(res_none$dt), "bg4_diff_normNone.txt", sep = "\t")
fwrite(annotate(res_tmm$dt),  "bg4_diff_normTMM.txt",  sep = "\t")

# BED of the robust set, for downstream motif / T-run analysis
robust <- res_tmm$dt[locus %in% sig_none & FDR < 0.05]
fwrite(robust[logFC > 0, .(chrom, start, end)],
       "diff_CSBplus_up.bed", sep = "\t", col.names = FALSE)
fwrite(robust[logFC < 0, .(chrom, start, end)],
       "diff_CSBminus_up.bed", sep = "\t", col.names = FALSE)


# --- 9. input check (copy number control) -----------------------------------
# CS1AN parental and the CSB-reinstated line may differ in copy number.
# A differential region should NOT show a matching shift in the inputs.

if (ncol(inp) > 0) {
  inp_cpm <- t(t(inp) / colSums(inp)) * 1e6
  plus_i  <- rowMeans(inp_cpm[, grepl("_CSB_", colnames(inp_cpm)), drop = FALSE])
  minus_i <- rowMeans(inp_cpm[, grepl("_WT_",  colnames(inp_cpm)), drop = FALSE])
  input_lfc <- data.table(locus = loci,
                          input_logFC = log2((plus_i + 1) / (minus_i + 1)))
  chk <- merge(res_tmm$dt[, .(locus, logFC, FDR)], input_lfc, by = "locus")
  message(sprintf(
    "\nCorrelation of IP logFC with input logFC at significant regions: %.3f",
    cor(chk[FDR < 0.05, logFC], chk[FDR < 0.05, input_logFC],
        use = "complete.obs")))
  message("A strong positive correlation implicates copy-number differences.")
  fwrite(chk, "input_copynumber_check.txt", sep = "\t")
}

message("\nDone.")

#File outputs:
#bg4_diff_normNone.txt is all regions, without normalization
#bg4_diff_normTMM.txt	is all regions, with TMM normalization
#diff_CSBplus_up.bed is regions UP in CSBplus without normalization
#diff_CSBminus_up.bed regions UP in CSBminus without normalization
#input_copynumber_check.txt is the logFC of input files
#bg4_diagnostics.pdf is a clustering file using MDS (similar to PCA)
#bg4_maplots.pdf is a logFC vs logCPM graph for each normalization

#Checking "stable copy number" in input, namely regions where the input is unchanged between CSB+/CSB-

library(data.table)
chk <- fread("/Users/naura/Downloads/input_copynumber_check.txt")
sig <- chk[FDR < 0.05]
nrow(sig)
nrow(sig[abs(input_logFC) < 0.5])

#Downstream analyses were performed only in regions where the input is "stable", i.e. logFC of input > 0.05

# Extract input-matched significant BG4 regions.
#
# Starting point: bg4_diff_normTMM.txt (differential test) and
# input_copynumber_check.txt (input logFC per region).
#
# A region is kept if:
#   FDR < 0.05                  edgeR called it differential
#   |input_logFC| < 0.5         input is comparable between the two cell
#                               lines, so the IP difference cannot be
#                               explained by DNA abundance / copy number
#
# NAMING: "CSB" files are CSB+ (reinstated); "WT" files are CSB- (parental,
# CSB-deficient). Positive logFC = higher in CSB+.
# ---------------------------------------------------------------------------

library(data.table)

FDR_CUT   <- 0.05
INPUT_CUT <- 0.5

diff <- fread("/Users/naura/Downloads/bg4_diff_normTMM.txt")
chk  <- fread("/Users/naura/Downloads/input_copynumber_check.txt")

# input_logFC is the only column needed from the check file
d <- merge(diff, chk[, .(locus, input_logFC)], by = "locus")

message(sprintf("regions tested:            %d", nrow(d)))
message(sprintf("significant (FDR < %.2f):   %d", FDR_CUT, sum(d$FDR < FDR_CUT)))

keep <- d[FDR < FDR_CUT & abs(input_logFC) < INPUT_CUT]
message(sprintf("input-matched significant: %d", nrow(keep)))

up   <- keep[logFC > 0][order(-logFC)]   # higher in CSB+
down <- keep[logFC < 0][order(logFC)]    # higher in CSB-
message(sprintf("  higher in CSB+: %d", nrow(up)))
message(sprintf("  higher in CSB-: %d", nrow(down)))

# --- BED for downstream tools (HOMER, bedtools, MEME) ----------------------
# 6 columns: chrom, start, end, name, score, strand.
# score is |logFC| scaled to BED's 0-1000 range; strand is "." since BG4
# peaks are unstranded.

to_bed <- function(dt) {
  dt[, .(chrom, start, end,
         name  = locus,
         score = pmin(round(abs(logFC) * 200), 1000),
         strand = ".")][order(chrom, start)]
}

fwrite(to_bed(up),   "BG4_CSBplus_up_inputmatched.bed",
       sep = "\t", col.names = FALSE)
fwrite(to_bed(down), "BG4_CSBminus_up_inputmatched.bed",
       sep = "\t", col.names = FALSE)
fwrite(to_bed(keep), "BG4_diff_all_inputmatched.bed",
       sep = "\t", col.names = FALSE)

# --- background set for motif analysis -------------------------------------
# Regions tested but NOT differential, with matched input. These are real
# BG4 peaks that did not change, so using them as HOMER background asks
# "what distinguishes the CSB-dependent G4s from G4s in general" rather
# than the trivial "are these G4s G-rich".

bg <- d[FDR >= 0.5 & abs(input_logFC) < INPUT_CUT]
message(sprintf("non-differential background: %d", nrow(bg)))
fwrite(to_bed(bg), "BG4_nondiff_background.bed", sep = "\t", col.names = FALSE)

# --- full tables with all statistics, for the record -----------------------

fwrite(up,   "BG4_CSBplus_up_inputmatched.txt",  sep = "\t")
fwrite(down, "BG4_CSBminus_up_inputmatched.txt", sep = "\t")

message("\nWritten:")
message("  BG4_CSBplus_up_inputmatched.bed   (higher in CSB+)")
message("  BG4_CSBminus_up_inputmatched.bed  (higher in CSB-)")
message("  BG4_diff_all_inputmatched.bed     (both directions)")
message("  BG4_nondiff_background.bed        (background for motif analysis)")
message("  ... plus .txt versions with logFC/FDR retained")
```
