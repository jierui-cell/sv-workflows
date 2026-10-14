# Seeded cohort-specific KL-NMF; input profiles are recorded in the output.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(10,12))
profile <- if (length(args) == 12) args[11] else 'panprostate_2026_core'
prefix <- if (length(args) == 12) args[12] else 'PCAWG.NMF'
stopifnot(prefix %in% c('PCAWG.NMF','ClusterSV.NMF'))
suppressPackageStartupMessages(library(NMF))
suppressPackageStartupMessages(library(doParallel))
out <- args[2]
rank_min <- as.integer(args[3]); rank_max <- as.integer(args[4]); fixed <- as.integer(args[5])
nrun <- as.integer(args[6]); seed <- as.integer(args[7]); cpus <- as.integer(args[8])
epsilon <- as.numeric(args[9]); max_iter <- as.integer(args[10])
stopifnot(is.finite(epsilon), epsilon > 0, max_iter > 0, cpus > 0, nrun > 0,
          rank_min >= 2, rank_max >= rank_min, fixed == 0 || fixed >= 2)
script_file <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE))
source(file.path(dirname(normalizePath(script_file)), 'stable_nmf.R'))
algorithm <- stable_brunet(epsilon)
tab <- read.delim(args[1], check.names = FALSE, row.names = 1)
a_full <- as.matrix(tab)
storage.mode(a_full) <- 'double'
stopifnot(!anyDuplicated(rownames(a_full)), !anyDuplicated(colnames(a_full)),
          all(is.finite(a_full)), all(a_full >= 0), all(a_full == floor(a_full)))
a <- a_full[rowSums(a_full) > 0, colSums(a_full) > 0, drop = FALSE]
if (!nrow(a) || !ncol(a)) stop('No nonzero samples/features for NMF')
if (min(dim(a)) < 2) stop('Too few nonzero samples/features for signature discovery')
registerDoParallel(cores = cpus)
parallel_options <- if (cpus == 1) '-p' else paste0('p', cpus)
set.seed(seed)
if (fixed > 0) {
    if (fixed > min(dim(a))) stop('Fixed rank exceeds data dimensions')
    rank <- fixed
    selection <- 'user_fixed_rank'
} else {
    upper <- min(rank_max, min(dim(a)))
    if (upper < rank_min) stop('Insufficient samples/features for requested rank_min; reduce rank_min or add samples')
    ranks <- seq.int(rank_min, upper)
    estimate <- nmfEstimateRank(a, ranks, method = algorithm, nrun = nrun, seed = seed,
                               .options = parallel_options, eps = epsilon, maxIter = max_iter)
    measures <- estimate$measures
    write.table(measures, file.path(out, 'NMF.rank_survey.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
    values <- measures$cophenetic
    changes <- diff(values)
    candidates <- which(is.finite(changes) & changes < 0)
    if (length(candidates)) {
        best <- candidates[changes[candidates] == min(changes[candidates])]
        rank <- min(measures$rank[best])
        selection <- 'rank_before_steepest_adjacent_cophenetic_drop'
    } else {
        usable <- which(is.finite(values))
        rank <- if (length(usable)) min(measures$rank[usable[values[usable] == max(values[usable])]]) else min(ranks)
        selection <- 'smallest_rank_at_maximum_cophenetic_or_smallest_admissible_rank'
    }
}
fit <- nmf(a, rank = rank, method = algorithm, nrun = nrun, seed = seed,
           .options = parallel_options, eps = epsilon, maxIter = max_iter)
w <- basis(fit); h <- coef(fit)
scale <- colSums(w)
if (any(!is.finite(scale) | scale <= 0)) stop('Degenerate NMF component')
w <- sweep(w, 2, scale, '/')
h <- sweep(h, 1, scale, '*')
names <- paste0('cSV', seq_len(rank))
W <- matrix(0, nrow(a_full), rank, dimnames = list(rownames(a_full), names))
H <- matrix(0, rank, ncol(a_full), dimnames = list(names, colnames(a_full)))
W[rownames(w),] <- w; H[,colnames(h)] <- h
write.table(data.frame(Feature = rownames(W), W, check.names = FALSE),
            file.path(out, paste0(prefix,'.signatures.tsv')), sep = '\t', row.names = FALSE, quote = FALSE)
attribution <- sweep(W, 2, rowSums(H), '*')
attribution_totals <- rowSums(attribution)
attribution[attribution_totals > 0,] <- attribution[attribution_totals > 0,,drop = FALSE] / attribution_totals[attribution_totals > 0]
write.table(data.frame(Feature = rownames(attribution), attribution, check.names = FALSE),
            file.path(out, paste0(prefix,'.feature_attribution.tsv')), sep = '\t', row.names = FALSE, quote = FALSE)
write.table(data.frame(Samples = colnames(H), t(H), check.names = FALSE),
            file.path(out, paste0(prefix,'.activities.tsv')), sep = '\t', row.names = FALSE, quote = FALSE)
fractions <- t(H)
totals <- rowSums(fractions)
fractions[totals > 0,] <- fractions[totals > 0,,drop = FALSE] / totals[totals > 0]
write.table(data.frame(Samples = rownames(fractions), fractions, check.names = FALSE),
            file.path(out, paste0(prefix,'.fractions.tsv')), sep = '\t', row.names = FALSE, quote = FALSE)
reconstructed <- W %*% H
cosine <- vapply(seq_len(ncol(a_full)), function(i) {
    denom <- sqrt(sum(a_full[,i]^2) * sum(reconstructed[,i]^2))
    if (denom == 0) NA_real_ else sum(a_full[,i] * reconstructed[,i]) / denom
}, numeric(1))
qc <- data.frame(Samples = colnames(a_full), input_feature_total = colSums(a_full),
                 reconstructed_total = colSums(reconstructed),
                 absolute_residual = colSums(abs(a_full - reconstructed)), cosine = cosine)
write.table(qc, file.path(out, 'NMF.reconstruction_qc.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
saveRDS(fit, file.path(out, 'NMF.fit.rds'))
jsonlite::write_json(list(status = 'completed', method = 'brunet_KL_NMF_with_per_iteration_factor_floor',
    feature_profile = profile, output_prefix = prefix, rank_selection = selection,
    rank = rank, nrun = nrun, seed = seed, cpus = cpus,
    factor_floor = epsilon, factor_floor_interval = 1, max_iterations = max_iter,
    copied_factor_updates = TRUE,
    stop_connectivity_checks = 40, connectivity_check_interval = 10,
    samples = ncol(a_full), positive_samples = ncol(a), channels = nrow(a_full),
    note = 'Cohort-specific de novo signatures; cSV labels are not COSMIC SV1-SV12 or published PCAWG signature IDs.'),
    file.path(out, 'signature_status.json'), auto_unbox = TRUE, pretty = TRUE)
stopImplicitCluster()
