args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1], mustWork = TRUE)
project <- file.path(root, 'environments/clustersv-r43')
bootstrap <- file.path(root, 'environments/clustersv-bootstrap')
dir.create(bootstrap, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(bootstrap, .libPaths()))
options(repos = c(CRAN = 'https://cloud.r-project.org'))
if (!requireNamespace('renv', quietly = TRUE)) install.packages('renv', lib = bootstrap)
lock <- file.path(root, 'sv_signatures/pcawg/renv.lock')
dir.create(project, recursive = TRUE, showWarnings = FALSE)
renv::init(project = project, bare = TRUE, restart = FALSE)
renv::settings$bioconductor.version('3.17', project = project)
if (file.exists(lock)) {
    renv::restore(project = project, lockfile = lock, prompt = FALSE)
} else {
    renv::install(c('BiocManager', 'bioc::Biobase', 'gtools@3.9.5',
                   'doParallel@1.0.17', 'R.utils@2.13.0', 'jsonlite@2.0.0',
                   'NMF@0.28', 'tibble', 'dplyr', 'tidyr', 'stringr', 'forcats'),
                  project = project)
    renv::snapshot(project = project, lockfile = lock, type = 'all', prompt = FALSE)
}
lib <- renv::paths$library(project = project)
writeLines(lib, file.path(project, 'library-path.txt'))
.libPaths(c(lib, .libPaths()))
pkgs <- c('gtools', 'doParallel', 'R.utils', 'NMF', 'jsonlite', 'tibble',
          'dplyr', 'tidyr', 'stringr', 'forcats')
stopifnot(all(vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)))
write.table(installed.packages(lib.loc = lib)[, c('Package', 'Version')],
            file.path(project, 'packages.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
cat('ClusterSV/NMF environment ready: ', R.version.string, '\n')
