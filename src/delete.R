asynch.samples = function(pheno, probs, expr, covar=NULL) {
    
    samples = intersect(rownames(pheno), rownames(probs))
    samples = intersect(samples, rownames(expr))
    if(!is.null(covar)) {
        samples = intersect(samples, rownames(covar))
    }
    
    if(length(samples) == 0) {
        stop("There are no samples in common. Please check the rownames on all objects.")
    }
    
    samples = sort(samples)
    pheno = pheno[samples,,drop = FALSE]
    probs = probs[samples,,drop = FALSE]
    expr  = expr[samples,,drop = FALSE]
    
    if(!missing(covar)) {
        covar = covar[samples,,drop = FALSE]
    }
    
    message(paste("Scanning with", nrow(pheno), "samples."))
    
    return(list(pheno = pheno, probs = probs, expr = expr, covar = covar))
    
} # synch.samples()




amediation.scan <- function(target, 
                           mediator, 
                           annotation, 
                           qtl.geno, 
                           covar=NULL,
                           method=c("double-lod-diff", "ignore", "lod-diff", 
                                    "lod-ratio"), 
                           verbose=TRUE) {
    
    # calculates log10-Likelihood of linear model y ~ 1 + X
    LL <- function(y, X) {
        -length(y)/2*log10(sum(qr.resid(qr(cbind(X,1)),y)^2))
    }
    
    # Synch sample IDs.
    tmp = synch.samples(pheno = target, probs = qtl.geno, expr = mediator, covar = covar)
    target   = tmp$pheno
    qtl.geno = tmp$probs
    mediator = tmp$expr
    covar    = tmp$covar
    rm(tmp)
    
    # check input
    stopifnot(NROW(target) == NROW(mediator))
    stopifnot(NROW(annotation) == NCOL(mediator))
    stopifnot(NROW(qtl.geno) == NROW(target))
    stopifnot(!any(is.na(qtl.geno)))
    stopifnot(all(is.numeric(target[,1])))
    stopifnot(all(is.numeric(mediator)))
    stopifnot(all(is.numeric(qtl.geno)))
    stopifnot(c("CHR", "MIDDLE_POINT") %in% toupper(names(annotation)))
    method = match.arg(method)
    
    if (!is.null(covar)) {
        stopifnot(NROW(target) == NROW(covar))
        stopifnot(!any(is.na(covar)))
        stopifnot(all(is.numeric(covar)))
    }
    
    
    # data preparation
    mediator <- cbind(mediator) # to ensure 'mediator' is a matrix
    N <- ncol(mediator) # number of points to scan
    if (is.null(covar)) covar <- cbind(rep(1, length(target))) # if no covariates, use just intercept
    LOD <- rep(NA, N) # prepare output
    
    if (method == "double-lod-diff") {
        no.na <- !is.na(target)
        LOD0 <- LL(target[no.na], cbind(covar, qtl.geno)[no.na,]) - LL(target[no.na], covar[no.na,])
    }
    
    # for-loop comparing M0: target~covar+mediator[,i] vs M1: target~covar+mediator[,i]+qtl.geno
    for (i in 1:N) {
        
        if (verbose & i %% 1000 == 0) print(i)
        
        no.na <- !is.na(target) & !is.na(mediator[,i])
        loglik0 <- LL(target[no.na], cbind(covar[no.na,], mediator[no.na,i]))
        loglik1 <- LL(target[no.na], cbind(covar[no.na,], mediator[no.na,i], qtl.geno[no.na,]))
        
        if (method == "ignore" | (method == "double-lod-diff" & all(no.na))) {
            # "double-lod-diff" for no missing observation is identical to "ignore"
            LOD[i] <- loglik1 - loglik0
        } else {
            loglik2 <- LL(target[no.na], covar[no.na,])
            loglik3 <- LL(target[no.na], cbind(covar[no.na,], qtl.geno[no.na,]))
            
            if (method == "lod-diff") {
                LOD[i] <- loglik3 - loglik2 - (loglik1-loglik0)
            } else if (method == "double-lod-diff") {
                LOD[i] <- LOD0 - (loglik3 - loglik2 - (loglik1-loglik0))
            } else if (method == "lod-ratio") {
                LOD[i] <- (10^loglik1-10^loglik0) / (10^loglik3 - 10^loglik2)
            }
        }
    }
    
    output <- annotation
    output$LOD <- LOD
    class(output) <- c("mediation", "data.frame")
    return(output)
}