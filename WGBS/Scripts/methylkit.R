# nolint start
suppressMessages(library(methylKit))
suppressMessages(library("genomation"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

run_methylkit <- function(
    dir,
    run,
    regions,
    file.list,
    sample.id,
    treatement,
    assembly,
    pipeline,
    output,
    mincov = 5,
    min_methylation_diff = 10,
    min_pvalue = 0.05,
    min_DMR_length = 200) {
  tryCatch(
    {
      message("--- Reading the methylation files ---")
      myobj <- methRead(file.list,
        sample.id = sample.id,
        assembly = assembly,
        treatment = treatement,
        context = "CpG",
        mincov = mincov,
        pipeline = pipeline
      )
      message("--- Done ---")

      message("--- Reading the regions ---")
      Rmsk_TEs <- read.table(regions)
      Rmsk_TEs <- Rmsk_TEs %>%
        dplyr::select(V1, V2, V3) %>%
        mutate(strand = "*")
      colnames(Rmsk_TEs) <- c("seqnames", "start", "end", "strand")
      Rmsk_TEs$strand <- gsub("C", "*", Rmsk_TEs$strand)
      g_Rmsk_TEs <- GenomicRanges::makeGRangesFromDataFrame(Rmsk_TEs, keep.extra.columns = FALSE)
      message("--- Done ---")

      message("--- Calculating DMRs ---")
      elements <- regionCounts(myobj, g_Rmsk_TEs)
      meth_DMR <- methylKit::unite(elements, destrand = FALSE)
      myDiff <- calculateDiffMeth(meth_DMR)
      message("--- Done ---")

      message("--- Saving results ---")
      message("--- Step 1 ---")
      myDiff$diffexpressed <- "Not_significant"
      myDiff$diffexpressed[myDiff$meth.diff > min_methylation_diff & myDiff$pvalue < min_pvalue] <- "Hypermethylated"
      myDiff$diffexpressed[myDiff$meth.diff < -min_methylation_diff & myDiff$pvalue < min_pvalue] <- "Hypomethylated"

      message("--- Step 2 ---")
      DMRs <- getData(myDiff) %>%
        subset(diffexpressed != "Not_significant")
      DMRs$length <- DMRs$end - DMRs$start
      DMRs <- DMRs %>%
        dplyr::filter(length >= min_DMR_length)

      message("--- Step 3 ---")
      save(myDiff, DMRs, file = paste0(dir, "/METHYLATION/DMR/methylkit_analysis.RData"))
      
      myDiff_header = colnames(myDiff)
      write.table(myDiff, paste0(dir, "/METHYLATION/DMR/ALL_methylkit_regions.bed"), row.names = FALSE, sep = "\t", col.names = FALSE, quote = FALSE)
      write.table(myDiff_header, paste0(dir, "/METHYLATION/DMR/ALL_methylkit_regions_header.txt"), row.names = FALSE, col.names = FALSE)
      write.csv(DMRs, paste0(dir, "/METHYLATION/DMR/methylkit_DMRs_pvalue", min_pvalue, ".csv"), row.names = FALSE)

      message("--- Step 4 ---")
      df_summary <- DMRs %>%
        group_by(diffexpressed) %>%
        summarise(count = n(), .groups = "drop") %>%
        mutate(percentage = count / sum(count) * 100) %>%
        mutate(label = paste0(count, " (", round(percentage, 0), "%)")) %>%
        arrange(desc(diffexpressed)) %>%
        mutate(cumulative_count = cumsum(count))

      total_DMRs <- sum(df_summary$count)

      perc <- ggplot(df_summary, aes(x = "", y = count, fill = diffexpressed)) +
        geom_bar(position = "stack", stat = "identity", width = 0.3) +
        geom_text(aes(label = count), position = position_stack(vjust = 0.5)) +
        theme_minimal(base_size = 15) +
        scale_fill_manual(values = c("Hypermethylated" = "#bb0c00", "Hypomethylated" = "#00AFBB")) +
        labs(fill = "", x = paste0("Total DMRs (pval <", min_pvalue, ") : ", total_DMRs), y = "Number of Regions")
      ggsave(filename = paste0(dir, "/METHYLATION/DMR/DMR_summary.png"), plot = perc, width = 6, height = 4)

      message("--- Step 5 ---")
      volcano <- ggplot(data = myDiff, aes(x = meth.diff, y = -log10(pvalue), col = diffexpressed)) +
        geom_vline(xintercept = c(-min_methylation_diff, min_methylation_diff), col = "gray", linetype = "dashed") +
        geom_hline(yintercept = -log10(min_pvalue), col = "gray", linetype = "dashed") +
        scale_color_manual(values = c("Hypermethylated" = "#bb0c00", "Not_significant" = "grey", "Hypomethylated" = "#00AFBB")) +
        labs(color = "Expression", x = expression("log"[2] * "FC"), y = expression("-log"[10] * "p-value")) +
        geom_point(size = 2) +
        theme_minimal()
      ggsave(filename = paste0(dir, "/METHYLATION/DMR/VolcanoPlot.png"), plot = volcano, width = 8, height = 8)
      message("--- Done ---")
    },
    error = function(e) {
      error_found <<- TRUE
      if (exists("myDiff")){
        save(myDiff, file = paste0(dir, "/METHYLATION/DMR/methylkit_analysis.RData"))
        
      }
      writeLines(as.character(e), con = paste0(dir, "/LOGS/METHYLATION_LOGS/DMR_LOGS/", run, "_error.txt"))
    }
  )
}

design_matrix <- read.csv(snakemake@params$design_matrix)
file.list <- as.list(snakemake@input$calls)
unique_conditions <- unique(design_matrix$Condition)
wt_index <- grep("wt", unique_conditions, ignore.case = TRUE)
if (length(wt_index) == 1) {
  unique_conditions <- c(unique_conditions[-wt_index], unique_conditions[wt_index])
}
condition_1_samples <- as.list(design_matrix$SampleID[design_matrix$Condition == unique_conditions[1]])
condition_2_samples <- as.list(design_matrix$SampleID[design_matrix$Condition == unique_conditions[2]]) 
sample.id <- as.list(c(condition_1_samples, condition_2_samples))
treatement <- c(rep(1, length(condition_1_samples)), rep(0, length(condition_2_samples)))


run_methylkit(
  dir = snakemake@params$run_dir,
  run = snakemake@params$run,
  regions = snakemake@params$regions,
  file.list = file.list,
  sample.id = sample.id,
  treatement = treatement,
  assembly = snakemake@params$genome,
  pipeline = snakemake@params$pp,
  output = snakemake@output,
  mincov = snakemake@params$mincov,
  min_methylation_diff = snakemake@params$min_methylation_diff,
  min_pvalue = snakemake@params$pvalue,
  min_DMR_length = snakemake@params$min_DMR_length
)


# nolint end