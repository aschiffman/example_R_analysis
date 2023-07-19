# An example R script for:
# Loading peaks from a bed file
# Counting the # of reads in those peaks (makes summarizedexperiment object)
# Author: Allison Schiffman

# import libraries (not all are needed)
library(pheatmap)
library(cluster)
library(edgeR)
library(csaw)
library(RColorBrewer)
library(colorspace)
library(assertthat)
library(ggforce)
library(ggsignif)
library(ggpubr)
library(ggcorrplot)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(tidyverse)

# set up ----
# Desktop paths
directory="~/R_analysis/my_experiment/" # main folder
function_directory = "~/R_analysis/" # folder with de_novo_enhancers_util.R file

# Sets directories for saving
results_directory = str_c(directory, "results/", Sys.Date(),"/")
qc_dir = str_c(directory, "qc/", Sys.Date(),"/")
dir.create(qc_dir, recursive = T, showWarnings = F)
dir.create(results_directory, recursive = T, showWarnings = F)

# Function for defining a GRanges object from Macs3 peaks
get_peaks = function(sample_names, directory, bed_paths, skip=25){
  peak_list = vector(mode = "list", length=length(sample_names))
  for (i in c(1:length(sample_names))) {
    sn = sample_names[i]
    
    pks = read_tsv(str_c(directory, bed_paths[i]), skip=skip)
    
    rgs = GRanges(seqnames = Rle(pks$chr),
                  ranges = IRanges(pks$start, end = pks$end),
                  fold_enrichment = pks$fold_enrichment)
    
    peak_list[[i]] = rgs
  }
  
  peak_list = GRangesList(peak_list)
  
  peaks = GenomicRanges::reduce(unlist(peak_list))
  
  return(list(peaks,peak_list))
}


ranged_expmt_to_l2rpkm = function(ranged_expmt, library_size){
  dge = asDGEList(ranged_expmt)
  colnames(dge) = colnames(ranged_expmt)
  rownames(dge) = rownames(ranged_expmt)
  dge$genes$length=ranged_expmt@rowRanges@ranges@width
  dge$samples$lib.size = library_size
  l2rpkm = as_tibble(rpkm(dge, log = T), rownames=NA)
  return(l2rpkm)
}

# Make data frame of sample names and file paths to bed files (from MACS3) and library size
metadata = bind_cols(samples = c("WT_unstim","KO_unstim","WT_8h","KO_8h"),
                     peak_path = c("Rep1/chip-seq/peaks/chip_WT_unstim_peaks.xls",
                                   "Rep1/chip-seq/peaks/chip_KO_unstim_peaks.xls",
                                   "Rep1/chip-seq/peaks/chip_WT_8h_peaks.xls",
                                   "Rep1/chip-seq/peaks/chip_KO_8h_peaks.xls"),
                     path= c("Rep1/chip-seq/bams/chip_WT_unstim.bam",
                             "Rep1/chip-seq/bams/chip_KO_unstim.bam",
                             "Rep1/chip-seq/bams/chip_WT_8h.bam",
                             "Rep1/chip-seq/bams/chip_KO_8h.bam"),
                     library_size=c(14444444, 13333333, 15555555, 12345678))
sample_names = metadata$sample # make a list of all sample names

##### Call peaks in samples----
tmp = get_peaks(sample_names, directory, metadata$peak_path) # Uses my get_peaks function, inputs: sample_names, directory, list of full paths to bed files
stim_peaks = tmp[[1]]
# count reads in peaks to granges object
peaks = regionCounts(str_c(directory, metadata$path), stim_peaks) # uses csaw function
# save peak counts
save(peaks, file=str_c(directory, "all_chip_peaks.RData"))

# Count RPKM (Reads per Kilo base of transcript per million dedupped reads)
peaks_rpkm = ranged_expmt_to_l2rpkm(peaks, metadata$library_size) # my function to count RPKM in each peak