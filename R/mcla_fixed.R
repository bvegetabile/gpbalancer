mcla_fixed <- function(y, covmat, tol=1e-3, max_iters =20, verbose=F){
    n_obs <- length(y)
    found_classes <- as.vector(unique(y))
    found_classes <- found_classes[order(found_classes)]
    n_classes <- length(found_classes)
    targets <- matrix(NA, nrow=n_obs, ncol=n_classes)
    for(c in 1:n_classes){
        targets[,c] <- ifelse(y == found_classes[c], 1, 0)
    }

    res <- gp_mcla(covmat, targets, n_classes, tol=tol, max_iters=max_iters, verbose=verbose)
    return(res)
}
