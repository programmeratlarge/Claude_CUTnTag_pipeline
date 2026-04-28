#!/usr/bin/env Rscript
# annotate_peaks.R
# Peak annotation with ChIPseeker.
# Inputs : peaks (narrowPeak/broadPeak/BED), GTF, TSS window
# Outputs: annotation TSV, pie chart, feature distribution barplot, MultiQC TSV
suppressPackageStartupMessages({
    library(optparse)
    library(GenomicFeatures)
    library(ChIPseeker)
    library(rtracklayer)
    library(ggplot2)
})

opt_list <- list(
    make_option("--peaks", type="character", help="Peaks file (narrowPeak/broadPeak/BED)"),
    make_option("--gtf",   type="character", help="Annotation GTF/GFF"),
    make_option("--tss_window", type="integer", default=3000),
    make_option("--label", type="character", help="Identifier for output filenames"),
    make_option("--out_anno", type="character"),
    make_option("--out_pie",  type="character"),
    make_option("--out_feat", type="character"),
    make_option("--out_mqc",  type="character")
)
opt <- parse_args(OptionParser(option_list=opt_list))

# Build TxDb from the GTF on the fly (works for any GTF; no organism assumption)
txdb <- suppressWarnings(makeTxDbFromGFF(opt$gtf, format="auto"))

# Read peaks
peaks <- tryCatch(
    readPeakFile(opt$peaks),
    error=function(e) {
        message("WARN: could not read peaks: ", conditionMessage(e))
        GenomicRanges::GRanges()
    }
)

if (length(peaks) == 0) {
    # Empty peak set: write empty outputs and a placeholder MQC
    write.table(data.frame(), file=opt$out_anno, sep="\t", quote=FALSE, row.names=FALSE)
    cat(sprintf(
"# id: peak_annotation_%s
# section_name: 'Peak annotation: %s'
# plot_type: 'bargraph'
Category\tCount
NoPeaks\t0
", opt$label, opt$label), file=opt$out_mqc)
    quit(status=0)
}

annoData <- annotatePeak(
    peaks,
    TxDb = txdb,
    tssRegion = c(-opt$tss_window, opt$tss_window),
    level = "gene",
    verbose = FALSE
)

# Annotation TSV
df <- as.data.frame(annoData)
write.table(df, file=opt$out_anno, sep="\t", quote=FALSE, row.names=FALSE)

# Pie chart
pdf(opt$out_pie, width=7, height=7)
plotAnnoPie(annoData)
dev.off()

# Feature distribution
pdf(opt$out_feat, width=8, height=4)
plotAnnoBar(annoData)
dev.off()

# MultiQC custom content: feature category counts
af <- annoData@annoStat
af$Feature <- as.character(af$Feature)
mqc_lines <- c(
    sprintf("# id: peak_annotation_%s", opt$label),
    sprintf("# section_name: 'Peak annotation: %s'", opt$label),
    "# description: 'Distribution of peaks across genomic features.'",
    "# plot_type: 'bargraph'",
    sprintf("# pconfig:"),
    sprintf("#     id: 'peak_anno_%s_bg'", opt$label),
    sprintf("#     title: 'Peak annotation distribution: %s'", opt$label),
    sprintf("#     ylab: 'Fraction (%%)'"),
    "Category\tFrequency"
)
for (i in seq_len(nrow(af))) {
    mqc_lines <- c(mqc_lines, sprintf("%s\t%.4f", af$Feature[i], af$Frequency[i]))
}
writeLines(mqc_lines, opt$out_mqc)
