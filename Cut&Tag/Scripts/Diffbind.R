library(DiffBind)
library(Rsamtools)
library(dplyr)
library(parallel)
library(ggplot2)

run_diffbind <- function(
    des_mat,
    mapq,
    output,
    condtion_name,
    mode = "cutntag",
    pvalue = 0.05,
    flip = TRUE) {
    samples <- read.csv(des_mat, header = TRUE)
    print(packageVersion("DiffBind"))
    print(samples)
    OBJ <- dba(sampleSheet = samples)
    tryCatch(
        {
            if (mode == "cutntag") {
                OBJ$config$yieldSize <- 15000000
                OBJ$config$intersectMode <- "IntersectionNotEmpty"
                OBJ$config$singleEnd <- FALSE
                OBJ$config$fragments <- FALSE
                OBJ$config$mapQCth <- mapq
                OBJ$config$minQCth <- mapq
                OBJ$config$inter.feature <- FALSE
                OBJ$config$fragmentSize <- 0
                OBJ$config$minOverlap <- 0
                OBJ$minOverlap <- 0
                OBJ$config$scanbamparam <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE), mapqFilter = 0)
                up <- 0.31
                down <- -0.5
                OBJ <- dba.count(DBA = OBJ, minOverlap = 0, bSubControl = FALSE, bRemoveDuplicates = TRUE, bUseSummarizeOverlaps = TRUE, bParallel = FALSE, score = DBA_SCORE_READS, summits = FALSE)
                OBJ <- dba.contrast(OBJ, minMembers=2, categories=DBA_CONDITION)
                OBJ <- dba.analyze(OBJ, method=DBA_DESEQ2)
                png(paste0(condtion_name, "_PCA.png"), width = 500, height = 350)
                dba.plotPCA(OBJ, attributes = DBA_CONDITION, label = DBA_ID)
                dev.off()
                png(paste0(condtion_name, "_Volcano.png"), width = 800, height = 600)
                dba.plotVolcano(OBJ, method = DBA_DESEQ2, bUsePval = TRUE, th = pvalue, bFlip = flip)
                dev.off()
                res_all <- dba.report(OBJ, method = DBA_DESEQ2, bUsePval = TRUE, bCalledDetail = TRUE, bCounts = TRUE, bNormalized = FALSE, bFlip = flip, th = pvalue)
                res_all_df <- as.data.frame(res_all)
                data_all <- res_all_df
                data_all$status <- "Not_significant"
                data_all$status[data_all$Fold > up & data_all$p.value <= pvalue] <- "Upregulated"
                data_all$status[data_all$Fold < down & data_all$p.value <= pvalue] <- "Downregulated"

                out_all <- data_all %>% dplyr::select(seqnames, start, end, strand, Fold, p.value, status)
                out_all_header <- c("seqnames", "start", "end", "strand", "Fold", "p.value", "status")
                write.table(out_all, paste0(condtion_name, "_diffbind.bed"), row.names = FALSE, sep = "\t", col.names = FALSE, quote = FALSE)
                write.table(out_all_header, paste0(condtion_name, "_diffbind_header.txt"), row.names = FALSE, col.names = FALSE)

                save(OBJ, data_all, file = output)
            } else {
                OBJ$config$singleEnd <- FALSE
                OBJ$config$mapQCth <- mapq
                OBJ$config$minQCth <- mapq
                OBJ$config$scanbamparam <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE), mapqFilter = 0)
                up <- 0
                down <- 0
                OBJ <- dba.count(DBA = OBJ, minOverlap = 0, bSubControl = FALSE, bRemoveDuplicates = TRUE, bUseSummarizeOverlaps = TRUE, bParallel = FALSE, score = DBA_SCORE_READS, summits = 75)
                OBJ <- dba.contrast(OBJ, minMembers=2, categories=DBA_CONDITION)
                OBJ <- dba.analyze(OBJ, method=DBA_DESEQ2)
                png(paste0(condtion_name, "_PCA.png"), width = 500, height = 350)
                dba.plotPCA(OBJ, attributes = DBA_CONDITION, label = DBA_ID)
                dev.off()
                png(paste0(condtion_name, "_Volcano.png"), width = 800, height = 600)
                dba.plotVolcano(OBJ, method = DBA_DESEQ2, bUsePval = TRUE, th = pvalue, bFlip = flip)
                dev.off()
                res_all <- dba.report(OBJ, method = DBA_DESEQ2, bUsePval = TRUE, bCalledDetail = TRUE, bCounts = TRUE, bNormalized = FALSE, bFlip = filp, th = pvalue)
                res_all_df <- as.data.frame(res_all)
                data_all <- res_all_df
                data_all$status <- "Not_significant"
                data_all$status[data_all$Fold > up & data_all$p.value <= pvalue] <- "Upregulated"
                data_all$status[data_all$Fold < down & data_all$p.value <= pvalue] <- "Downregulated"

                out_all <- data_all %>% dplyr::select(seqnames, start, end, strand, Fold, p.value, status)
                out_all_header <- c("seqnames", "start", "end", "strand", "Fold", "p.value", "status")
                write.table(out_all, paste0(condtion_name, "_diffbind.bed"), row.names = FALSE, sep = "\t", col.names = FALSE, quote = FALSE)
                write.table(out_all_header, paste0(condtion_name, "_diffbind_header.txt"), row.names = FALSE, col.names = FALSE)

                save(OBJ, data_all, file = output)
            }
        },
        error = function(e) {
            error_found <<- TRUE
            save(OBJ, file = output)
            writeLines(as.character(e), con = paste0(condtion_name, "_error.txt"))
        }
    )
}

matrices <- snakemake@input[["matrices"]]

for (mat in matrices) {
    condtion_name <- tools::file_path_sans_ext(mat)
    basename <- basename(condtion_name)
    match <- regmatches(basename, regexpr("(mapq)(\\d+)", basename))
    mapq <- as.integer(gsub("mapq", "", match))
    run_diffbind(
        des_mat = mat,
        mapq = mapq,
        condtion_name = condtion_name,
        output = paste0(condtion_name, ".Rdata"),
        mode = snakemake@params[["mode"]],
        pvalue = snakemake@params[["pvalue"]],
        flip = snakemake@params[["flip"]]
    )
}
