# nolint start
library(dplyr)
library(ggplot2)
library(tidyr)

make_plots <- function(dir, regions, pvalue_cutoff = 0.05, meth_diff_cutoff = 0.1, DMR_min_length = 200) {
    tryCatch(
        {
            message("--- Reading the annotation file ---")
            overlaps <- read.table(regions)
            colnames(overlaps) <- c("seqnames", "start", "end", "strand", "pvalue", "qvalue", "meth.diff", "diffexpressed", "annot.seqnames", "annot.Start", "annot.End", "annot.strand", "annot.Rep_Name", "annot.Family.Class", "X", "Y", "Z", "annot.length")
            overlaps <- overlaps %>%
                separate(annot.Family.Class, into = c("Class", "Family"), sep = "/", remove = FALSE) %>%
                subset(annot.Family.Class != ".") %>%
                mutate(region_length = end - start) %>%
                subset(pvalue < pvalue_cutoff) %>%
                subset(meth.diff > meth_diff_cutoff | meth.diff < -meth_diff_cutoff) %>%
                subset(region_length >= DMR_min_length)
            message("--- Done ---")

            message("--- Plotting the annotation ---")
            num <- as.data.frame(table(overlaps$Class)) %>%
                mutate(percentage = Freq / sum(Freq) * 100) %>%
                mutate(label = paste0(round(percentage, 0), "% ", Var1, " (", Freq, ")")) %>%
                mutate(cumulative_count = cumsum(Freq))

            sum <- ggplot(num, aes(x = "", y = Freq, fill = Var1)) +
                geom_bar(stat = "identity", width = 0.4, color = "black", size = 0.2) +
                theme_void() +
                labs(fill = "", x = "", y = "") +
                scale_fill_discrete(labels = num$label)
            ggsave(filename = paste0(dir, "/METHYLATION/DMR/TE_summary.png"), plot = sum, width = 4, height = 4)

            Down_reg_TE <- overlaps %>%
                dplyr::filter(diffexpressed == "Hypomethylated") %>%
                dplyr::select(seqnames, start, end, diffexpressed, annot.Family.Class) %>%
                arrange(annot.Family.Class)

            Up_reg_TE <- overlaps %>%
                dplyr::filter(diffexpressed == "Hypermethylated") %>%
                dplyr::select(seqnames, start, end, diffexpressed, annot.Family.Class) %>%
                arrange(annot.Family.Class)

            Down_reg_TE_freq <- as.data.frame(table(Down_reg_TE$annot.Family.Class)) %>%
                mutate(SUM = sum(Freq), Percent = (Freq / SUM) * 100) %>%
                arrange(Var1)

            Up_reg_TE_freq <- as.data.frame(table(Up_reg_TE$annot.Family.Class)) %>%
                mutate(SUM = sum(Freq), Percent = (Freq / SUM) * 100) %>%
                arrange(Var1)

            combined_data <- rbind(
                Up_reg_TE_freq %>% mutate(Regulation = "Hypermethylated"),
                Down_reg_TE_freq %>% mutate(Regulation = "Hypomethylated")
            )

            combined_data <- combined_data %>%
                separate(Var1, into = c("Class", "Family"), sep = "/", remove = FALSE)

            sum_data <- combined_data %>%
                group_by(Class, Regulation) %>%
                summarise(Sum = sum(Freq))

            df2 <- left_join(combined_data, sum_data, by = c("Class", "Regulation"))

            ove <- ggplot(df2, aes(x = Class, y = Freq, fill = Var1)) +
                geom_bar(position = "stack", stat = "identity", width = 0.8, color = "black", size = 0.1) +
                theme_minimal(base_size = 12) +
                geom_text(aes(label = ifelse(Freq / Sum > 0.1, Freq, "")), position = position_stack(vjust = 0.5), colour = "white") +
                geom_text(aes(y = Sum, label = Sum), vjust = -0.5, colour = "black") +
                labs(fill = "TEs", x = "TEs", y = "") +
                facet_wrap(~Regulation) +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                scale_fill_manual(values = c("DNA/hAT-Charlie" = "#F8766D", "LINE/L1" = "#7CAE00", "LINE/L2" = "#0BB500", "LTR/ERV1" = "#00BE67", "LTR/ERVK" = "#00C19A", "LTR/ERVL" = "#00B8E7", "LTR/ERVL-MaLR" = "#00A9FF", "SINE/Alu" = "#8494FF", "SINE/B4" = "#ED68ED", "SINE/B2" = "#C77CFF"))
            ggsave(filename = paste0(dir, "/METHYLATION/DMR/TE_summary_2.png"), plot = ove, width = 8, height = 10)
            message("--- Done ---")
        },
        error = function(e) {
            error_found <<- TRUE
            writeLines(as.character(e), con = paste0(dir, "/LOGS/METHYLATION_LOGS/DMR_LOGS/", "Plot_annotation_error.txt"))
        }
    )
}

dir <- snakemake@params$run_dir
regions <- snakemake@input$regions
pvalue_cutoff <- snakemake@params$pvalue_cutoff
meth_diff_cutoff <- snakemake@params$meth_diff_cutoff
DMR_min_length <- snakemake@params$min_DMR_length

make_plots(dir, regions, pvalue_cutoff, meth_diff_cutoff, DMR_min_length)

# nolint end
