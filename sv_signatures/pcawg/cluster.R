# Use upstream scoring and footprint functions; expose reviewed cutoffs only.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 9)
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(doParallel))
source(args[5])
chr_sizes <- get_chr_sizes(args[2])
genome_size <- sum(chr_sizes)
centromere_telomere_coords <- get_centromere_telomere_coords(args[3])
registerDoParallel(cores = as.integer(args[6]))
# Assert the upstream syntax before altering the two documented defaults.
stopifnot(as.character(body(clust_rgs)[[2]][[2]]) == 'cutoff_1',
          body(clust_rgs)[[2]][[3]] == 0.01,
          as.character(body(clust_rgs)[[3]][[2]]) == 'cutoff_2',
          body(clust_rgs)[[3]][[3]] == 0.05)
body(clust_rgs)[[2]] <- substitute(cutoff_1 <- VALUE, list(VALUE = as.numeric(args[7])))
body(clust_rgs)[[3]] <- substitute(cutoff_2 <- VALUE, list(VALUE = as.numeric(args[8])))
formals(get_footprints)$cutoff <- as.numeric(args[9])
if (file.info(args[1])$size == 0) {
    file.create(paste0(args[4], '.sv_clusters_and_footprints.tsv'))
    file.create(paste0(args[4], '.sv_distance_pvals'))
    quit(status = 0)
}
d <- read.table(args[1], header = FALSE, sep = '\t', stringsAsFactors = FALSE,
                comment.char = '', quote = '', colClasses = c('character', NA,
                NA, 'character', NA, NA, 'character', NA,
                'character', 'character'))
stopifnot(ncol(d) == 10, all(d[,1] %in% names(chr_sizes)), all(d[,4] %in% names(chr_sizes)))
res <- compute_clusters_and_footprints(d)
last <- vapply(seq_len(max(res$clust_and_fps$clust)), function(idx) {
    members <- res$clust_and_fps$clust == idx
    if (sum(members) == 1) return(NA_real_)
    values <- res$m[members, members]
    values <- values[!is.na(values) & values < res$clustering_fdr_cutoff]
    if (length(values)) max(values) else 0
}, numeric(1))
# Preserve upstream values, but do not label this comparison-count-scaled
# clustering statistic as an adjusted FDR or an ordinary raw p-value.
out <- data.frame(res$clust_and_fps, last_agglomeration_value = last[res$clust_and_fps$clust])
write.table(out, paste0(args[4], '.sv_clusters_and_footprints.tsv'), sep = '\t',
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(res$m, paste0(args[4], '.sv_distance_pvals'), sep = '\t',
            row.names = FALSE, col.names = FALSE, quote = FALSE)
stopImplicitCluster()
