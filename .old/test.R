grab_residual_mat <- function(variable_mat, 
                              adjust_mat, 
                              varOfInterest, 
                              adjOfInterest,
                              use_qr = TRUE) {
    samples <- intersect(rownames(variable_mat), rownames(adjust_mat))
    
    full_dat <- cbind(variable_mat[samples,], adjust_mat[samples,]) %>% as.data.frame()
    
    if(use_qr) {
        #browser()
        # Fast, but can't handle NAs in y
        X.0 <- model.matrix(formula(paste("~ + 1 +", paste(adjOfInterest, collapse = " + "))), data = full_dat)
        qr.0 <- qr(X.0[samples,])
        
        resid_mat <- sapply(1:length(varOfInterest), function(i) {
            
            
            d <- full_dat[samples, varOfInterest[i]]
            return(qr.resid(qr.0, d))
            
            
            # remove NA's
            #qr.0$qr <- qr.0$qr[!is.na(d), , drop = FALSE]
            #s <- samples[!is.na(d)]
            #d <- d[!is.na(d)]
            
            # construct a matrix with and set the names
            #matrix_residulas <- matrix(qr.resid(qr.0, d))
            #rownames(matrix_residulas) <- s
            
            # construct a matrix to pass back, so it's always the same 
            # sample size
            #ret <- matrix(NA, nrow = length(samples))
            #rownames(ret) <- samples
            
            # set the data
            #ret[rownames(ret) %in% rownames(matrix_residulas), ] <- 
            #    matrix_residulas[rownames(matrix_residulas) %in% rownames(ret), ]
            
            #return(ret)
        }, simplify = TRUE)
        rownames(resid_mat) <- samples
    }
    else{
        ## Way too slow, use QR trick
        resid_mat <- sapply(1:length(varOfInterest), function(i) {
            fit_formula <- paste(varOfInterest[i], "~", paste(adjOfInterest, collapse = " + "))
            fit <- lm(formula(fit_formula), data = full_dat)
            return(fit$residuals[samples])
        }, simplify = TRUE)
    }
    
    colnames(resid_mat) <- varOfInterest
    resid_mat
}



