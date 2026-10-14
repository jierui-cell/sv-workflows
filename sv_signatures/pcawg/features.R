# Build the 44 core features from the pinned panprostate code, before its NMF step.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(forcats))
inputs <- read.delim(args[1], check.names = FALSE, stringsAsFactors = FALSE)
out <- normalizePath(args[3])
read_input <- function(path, columns) {
    x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, colClasses = 'character')
    absent <- setdiff(columns, colnames(x))
    if (length(absent)) stop(paste(path, 'missing columns:', paste(absent, collapse = ', ')))
    if (anyDuplicated(x$id)) stop(paste('Duplicate annotation id in', path))
    x
}
classes <- repl <- fragile <- cluster_maps <- list()
for (i in seq_len(nrow(inputs))) {
    sample <- inputs$pair_id[i]
    d <- read_input(inputs$classified_sv[i], c(paste0('V', 1:18), 'id', 'tumour'))
    r <- read_input(inputs$replication_timing[i], c('id', 'name_regionA1', 'name_regionA2', 'dist_regionA1', 'dist_regionA2'))
    f <- read_input(inputs$fragile_sites[i], c('id', 'name_regionA1', 'name_regionA2', 'dist_regionA1', 'dist_regionA2'))
    d$V11 <- as.numeric(d$V11); d$V12 <- as.numeric(d$V12)
    f$dist_regionA1 <- as.numeric(f$dist_regionA1); f$dist_regionA2 <- as.numeric(f$dist_regionA2)
    r$dist_regionA1 <- as.numeric(r$dist_regionA1); r$dist_regionA2 <- as.numeric(r$dist_regionA2)
    if (!identical(d$V7, d$id)) stop('Classification id must equal its V7 junction identifier')
    if (any(d$tumour != sample) || !setequal(d$id, r$id) || !setequal(d$id, f$id))
        stop(paste('Classification sample or annotation id mismatch:', sample))
    if (any(!is.finite(d$V11)) || any(!is.finite(d$V12)) || any(d$V11 < 1 | d$V12 < 1)) stop('Invalid refined coordinates')
    if (anyNA(r[,c('name_regionA1','name_regionA2')]) ||
            any(!as.matrix(r[,c('name_regionA1','name_regionA2')]) %in% c('early', 'mid', 'late')))
        stop('Replication timing must be measured and labelled early/mid/late for both endpoints')
    if (anyNA(f[,2:5])) stop('Fragile-site annotations must be supplied, including explicit non-site labels')
    # Avoid cross-patient cluster collisions in upstream distinct(V13) calls.
    if (nrow(d)) {
        if (anyNA(d$V13) || any(!nzchar(d$V13))) stop('Missing classification cluster ID')
        original_ids <- unique(d$V13)
        local_ids <- paste(sample, seq_along(original_ids), sep = ':')
        cluster_maps[[i]] <- data.frame(sample = sample, original_cluster_id = original_ids,
                                        local_cluster_id = local_ids)
        d$V13 <- local_ids[match(d$V13, original_ids)]
        d$id <- paste(sample, d$id, sep = ':')
        d$V7 <- d$id
        r$id <- paste(sample, r$id, sep = ':')
        f$id <- paste(sample, f$id, sep = ':')
    }
    d <- d[, c(paste0('V', 1:18), 'id', 'tumour')]
    classes[[i]] <- d; repl[[i]] <- r; fragile[[i]] <- f
}
d <- bind_rows(classes)
write.table(bind_rows(cluster_maps), file.path(out, 'classification_cluster_ids.tsv'),
            sep = '\t', row.names = FALSE, quote = FALSE)
if (!nrow(d)) stop('No classified events; signature extraction is undefined')
base_classes <- sub(':.*$', '', d$V18)
supported <- c('del', 'td', 'inter', 'fb_then_fb', 'fb', 'inv', 'direct_inversion',
               'inverted_duplication', 'inversion_gain_loss', 'dup_trp_dup',
               'foldback_transloc', 'inv_shard_transloc', 'TD_after_unbal_transloc',
               'chromoplexy_cycle_plus_templated_insertions', 'shard_cycle', 'shard_chain',
               'chromoplexy_cycle', 'chromoplexy_chain')
unknown <- setdiff(unique(base_classes), supported)
if (length(unknown)) stop(paste('Unmodelled cSV classes require review:', paste(unknown, collapse = ', ')))
work <- tempfile('features-')
dir.create(work)
complex_sv_path <- file.path(work, 'classes.tsv')
replicating_timing_path <- file.path(work, 'replication.tsv')
fragileSite_path <- file.path(work, 'fragile.tsv')
write.table(d, complex_sv_path, sep = '\t', row.names = FALSE, quote = FALSE)
write.table(bind_rows(repl), replicating_timing_path, sep = '\t', row.names = FALSE, quote = FALSE)
write.table(bind_rows(fragile), fragileSite_path, sep = '\t', row.names = FALSE, quote = FALSE)
useARBS <- FALSE; useCT <- FALSE; uniquePatients <- FALSE
output_path <- out
code <- readLines(args[2], warn = FALSE)
find_one <- function(pattern) {
    idx <- grep(pattern, code, fixed = TRUE)
    if (length(idx) != 1) stop(paste('Unreviewed upstream feature code marker:', pattern))
    idx
}
start_fn <- find_one('getAvgTempInsertionSize <- function(classData){')
end_fn <- find_one('# load data') - 2L
start_features <- find_one('complex_sv_class <- as_tibble(read.table(complex_sv_path,')
end_features <- find_one('#### nmf') - 1L
helper_code <- paste(code[start_fn:end_fn], collapse = '\n')
# Match an entire generated junction ID: junction 1 must not hide 10, 11, etc.
helper_code <- gsub('paste0("^",i)', 'paste0("^",i,"_")', helper_code, fixed = TRUE)
eval(parse(text = helper_code))
feature_code <- paste(code[start_features:end_features], collapse = '\n')
feature_code <- sub(code[start_features], 'complex_sv_class <- as_tibble(d)', feature_code, fixed = TRUE)
# Exclude exactly reciprocal two-junction cycles, not cycles numbered 20-29.
feature_code <- gsub('^chromoplexy_cycle:(?!2)', '^chromoplexy_cycle:(?!2$)', feature_code, fixed = TRUE)
# Upstream 1:0 loops fail when a cohort has no templated-insertion class.
for (item in c('nrow(tmpIns)', 'length(tmpIns_clust_ID)', 'length(size_tmpIns)')) {
    replacement <- if (startsWith(item, 'nrow')) paste0('seq_len(', item, ')') else
        paste0('seq_along(', sub('^length\\((.*)\\)$', '\\1', item), ')')
    feature_code <- gsub(paste0('1:', item), replacement, feature_code, fixed = TRUE)
}
eval(parse(text = feature_code))
counts <- read.table(file.path(out, 'complex_sv_counts.txt'), header = TRUE, sep = '\t',
                     row.names = 1, check.names = FALSE)
stopifnot(ncol(counts) == 44, all(is.finite(as.matrix(counts))), all(counts >= 0))
all_samples <- inputs$pair_id
complete <- matrix(0, nrow = ncol(counts), ncol = length(all_samples),
                   dimnames = list(colnames(counts), all_samples))
complete[, rownames(counts)] <- t(as.matrix(counts))
write.table(data.frame(Feature = rownames(complete), complete, check.names = FALSE),
            file.path(out, 'PCAWG.features.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
jsonlite::write_json(list(profile = 'panprostate_2026_core', channels = 44,
    samples = length(all_samples), classified_junctions = nrow(d),
    feature_total = sum(complete), use_ARBS = FALSE, use_chromothripsis = FALSE,
    note = 'Counts include junctions, clusters, or template lengths as defined upstream; replication timing and fragile-site channels can overlap.'),
    file.path(out, 'feature_qc.json'), auto_unbox = TRUE, pretty = TRUE)
unlink(work, recursive = TRUE)
