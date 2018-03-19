library(data.table)
library(dplyr)
library(intermediate)
library(jsonlite)
library(qtl2)
library(qtl2db)
library(RSQLite)
library(pryr)

source('/Users/mvincent/work/qtlviewer/qtlapi/qtlapi.R')
load('/Users/mvincent/work/qtlviewer/data/Attie_DO378_eQTL_viewer_v4.RData')

ds <- 'dataset.islet.rnaseq'
id <- 'ENSMUSG00000006732'
chrom <- '3'
blup <- FALSE

for (i in 1:10) {
    for (x in 1:10) {
        ptm <- proc.time()
        expr_results <- GetFoundercoefs('dataset.islet.rnaseq', 'ENSMUSG00000006732', chrom=chrom, blup=blup, nCores=i)
        elapsed <- proc.time() - ptm
        print(paste0(i, '||', x, '||', elapsed["elapsed"], '||', blup))
    }
}

blup <- TRUE

for (i in 1:10) {
    for (x in 1:10) {
        ptm <- proc.time()
        expr_results <- GetFoundercoefs('dataset.islet.rnaseq', 'ENSMUSG00000006732', chrom=chrom, blup=blup, nCores=i)
        elapsed <- proc.time() - ptm
        print(paste0(i, '||', x, '||', elapsed["elapsed"], '||', blup))
    }
}




