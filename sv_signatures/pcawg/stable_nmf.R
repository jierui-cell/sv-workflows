# NMF 0.28 floors factors only every tenth iteration. Sparse counts can
# underflow before then; preserve its official KL update and floor every step.
stable_brunet <- function(epsilon = .Machine$double.eps) {
    stopifnot(length(epsilon) == 1, is.finite(epsilon), epsilon > 0)
    algorithm <- NMF::nmfAlgorithm('brunet')
    original_update <- algorithm@Update
    algorithm@Update <- function(i, v, x, copy = TRUE, eps = epsilon, ...) {
        basis(x) <- pmax(basis(x), eps)
        coef(x) <- pmax(coef(x), eps)
        x <- original_update(i, v, x, copy = TRUE, eps = eps, ...)
        if (any(!is.finite(basis(x))) || any(!is.finite(coef(x))))
            stop(paste('Nonfinite KL-NMF factors at iteration', i))
        basis(x) <- pmax(basis(x), eps)
        coef(x) <- pmax(coef(x), eps)
        x
    }
    algorithm
}
