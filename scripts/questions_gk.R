# grab the correlations between 'G33_ins_secrete_gm' in dataset.exvivo$data
# and all the data in dataset.islet.rnaseq$data$rz
# make sure we are only using the same samples
samples <- intersect(rownames(dataset.exvivo$data), rownames(dataset.islet.rnaseq$data$rz))
# 378 samples

corValues <- cor(dataset.exvivo$data[samples, 'G33_ins_secrete_gm'], 
                 dataset.islet.rnaseq$data$rz[samples, ], use = "pair")

# reorder
corValues <- corValues[1, order(abs(corValues), decreasing = TRUE)]
# corValues is now a list of genes and correlations
head(corValues)

#ENSMUSG00000040435 ENSMUSG00000056749 ENSMUSG00000021108 ENSMUSG00000026628 ENSMUSG00000022718 
#0.6157800          0.6044411          0.5512109          0.5471041         -0.5349566 
#ENSMUSG00000002289 
#0.5334357

# now let's get all the x and y variables between 'G33_ins_secrete_gm' in dataset.exvivo$data
# and 'ENSMUSG00000026628' in dataset.islet.rnaseq$data$rz so we can plot
plotData <- tibble(x = dataset.exvivo$data[samples, 'G33_ins_secrete_gm'], 
                   y = dataset.islet.rnaseq$data$rz[samples, 'ENSMUSG00000026628'])

head(plotData)
# A tibble: 6 x 2
#x      y
#<dbl>  <dbl>
#  1 -1.45  -0.552
#2 -2.50  -0.738
#3 -0.651  0.334
#4 -1.89  -0.341
#5 -1.40  -0.567
#6 -1.04  -0.259

# simple look at the data
ggplot(plotData, aes(x=x, y=y)) + geom_point()


###
# NEW 
###

## Function for partial correlations
# x & y are the two variables of interest
# z is one or more covariates to regress out
grab_partial_correlations <- function(x, y, z_vector, data) {
  fit1_formula <- paste(y, "~", paste(z_vector, collapse = " + "))
  fit1 <- lm(formula(fit1_formula), data = data)
  
  fit2_formula <- paste(x, "~", paste(z_vector, collapse = " + "))
  fit2 <- lm(formula(fit2_formula), data = data)
  
  cor(fit1$residuals, fit2$residuals, use = "complete.obs")
}

# QUESTION 1: how do I get values equivalent to corValues above, but with a covariate adjustment?
# I would like a function where I can pass in 2 datasets, a variable of interest, and a covariate (or several)



###################################################################################################
#######
####### Hmmmmm.... This is a bit tricky because of how things are split across multiple datasets
#######
####### I'm thinking we should break it up a bit more
####### - make a function to pull residuals
####### - then basically do that for data matrices of
####### - - variable of interest
####### - - Comparison variables
####### 
####### For the residual function, I have a fast option that can't handle NAs (QR decomposition)
####### Taking residuals for a big matrix could be kind of slow
####### We could improve the QR decomposition method so it's more flexible
####### 
####### After making these matrices, things follow pretty nicely with the old functions
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
      return(qr.resid(qr.0, full_dat[samples, varOfInterest[i]]))
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
## Works for multiple variables
G33_ins_secrete_gm_resid <- grab_residual_mat(variable_mat = dataset.exvivo$data, 
                              adjust_mat = dataset.exvivo$covar.matrix, 
                              varOfInterest = c("G33_ins_secrete_gm"), 
                              adjOfInterest = "sexM")
rna_resid <- grab_residual_mat(variable_mat = dataset.islet.rnaseq$data$rz, 
                               adjust_mat = dataset.exvivo$covar.matrix, 
                               varOfInterest = colnames(dataset.islet.rnaseq$data$rz), 
                               adjOfInterest = "sexM")

corValuesAjd <- cor(G33_ins_secrete_gm_resid[samples,], 
                    rna_resid[samples,], 
                    use = "pair")
# reorder
corValuesAjd <- corValuesAjd[1, order(abs(corValuesAjd), decreasing = TRUE)]
# corValues is now a list of genes and correlations
head(corValuesAjd)

plotData <- tibble(x = G33_ins_secrete_gm_resid[samples,], 
                   y = rna_resid[samples, "ENSMUSG00000026628"])
head(plotData)
# x      y
# <dbl>  <dbl>
#   1 -0.342  -0.576
# 2 -1.39   -0.762
# 3  0.453   0.310
# 4 -0.787  -0.365
# 5 -0.298  -0.591
# 6  0.0652 -0.283
ggplot(plotData, aes(x=x, y=y)) + geom_point()


###################################################################################################
#######
####### Seems like we have potentially three matrices of interest
####### Two "phenotypes" and some covariate
####### Relaxing that should be simple in principle, it just makes managing all the pieces of data annoying
####### This current code works well for one main variable of interest (G33_ins_secrete_gm)
####### - compared to a matrix of other variables
####### - adjusted for a subset of covariates from another matrix
#######
####### Answering the following question is straightforward because we've basically just made adjusted matrices





# now corValuesAdjusted would be just like corValues above so I can make the list
# to display to users
#

# Let's say in our dataset, we have sex as a covariate and we want to adjust
# The only place we have sex data is in dataset.islet.rnaseq$annot.samples and
# dataset.islet.rnaseq$covar.matrix
#

# Would I create a vector/list called sexCovar and model something like
# sexCovar <- dataset.islet.rnaseq$covar.matrix[,'sexM']
# d <- dataset1[, varOfInterest]
# m <- model.matrix.lm(d~sexCovar, data = dataset1, na.action=na.pass)
# but now what do I do with dataset2 data?  loop through each and ((LOST))?
#


# QUESTION 2: how do I get values equivalent to plotData above, but with a covariate adjustment?
# I would like a function where I can pass in 2 datasets, 2 variables of interest, and a covariate (or several)

grab_plot_data <- function(dataset1, dataset2, varOfInterest1, varOfInterest2, covariate) {}
#

plotDataAdjusted <- grab_plot_data(dataset.exvivo$data, 
                                   dataset.islet.rnaseq$data$rz,
                                   'G33_ins_secrete_gm',
                                   'ENSMUSG00000026628',
                                   c('sex'))

# now plotDataAdjusted would be just like plotData above so I can plot it


