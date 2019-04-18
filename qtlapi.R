# #############################################################################
#
# Load the libraries will we need
#
# #############################################################################

library(data.table)
library(tibble)
library(dplyr)
library(intermediate)
library(jsonlite)
library(qtl2)
library(RSQLite)
library(pryr)
library(gtools)


# #############################################################################
#
# Load the main data file.
#
# This was intended to be used in a controlled Docker environment and has
# some assumptions as to where the RData file resides and how it should be
# named.  For SNP association mapping, there also needs to be a SQLite 
# database containing SNP information.
#
# Please see: https://github.com/churchill-lab/qtlapi/
#
# There are 2 methods in which we can utilize this script.
#
#     1. Pre load the RData file into the environment and set db.file to
#        the location of where the SNP database resides on disk.
#            
#                     - OR -
#
#     2. Place a file with ".RData" extension in /app/qtlapi/data as well as
#        a file with a ".sqlite" extension for the SNP database.
# 
# #############################################################################


# set a variable called debugMode to TRUE for step 1 above.

if (exists("debug.mode")) {
    debug.mode <- get("debug.mode")
} else {
    debug.mode <- FALSE
}

if (debug.mode) {
    print('DEBUG MODE: Make sure the data is loaded and db.file is defined')
} else {
    print("Finding the data file...")
    data.files <- list.files("/app/qtlapi/data", "*.RData", full.names=TRUE)

    if (length(data.files) == 0) {
        stop("There needs to be an .RData file in /app/qtlapi/data")
    } else if (length(data.files) > 1) {
        stop("There needs to be only 1 .RData file in /app/qtlapi/data")
    }

    print(paste0("Loading the data file: ", data.files))

    load(data.files, .GlobalEnv)

    db.file <- list.files("/app/qtlapi/data", "*.sqlite", full.names=TRUE)

    if (length(db.file) == 0) {
        stop("There needs to be an .sqlite file in /app/qtlapi/data")
    } else if (length(db.file) > 1) {
        stop("There needs to be only 1 .sqlite file in /app/qtlapi/data")
    }

    print(paste0("Using SNP db file: ", db.file))
}


# #############################################################################
#
# Utility functions
#
# #############################################################################


#' Opposite of %in%  
`%not in%` <- function (x, table) match(x, table, nomatch = 0L) == 0L

#' Check value to see if it could be a boolean.  
#' 
#' Acceptable TRUE boolean values are TRUE, 1, "T", "TRUE", "YES", "Y", "1"
#'
#' @param value the value to check
#' 
#' @return TRUE if it is a boolean value and it's value is "TRUE"
#' 
toBoolean <- function(value) {
    if (is.logical(value)) {
        return(value)
    } else if (is.character(value)) {
        return(toupper(value) %in% c("T", "TRUE", "YES", "Y", "1"))
    } else if (is.numeric(value)) {
        return(value == 1)
    }

    FALSE
}


#' Check a value for NULL and return default if the value is NULL.
#'
#' @param value the value to check
#' @param default the default value to return
#' 
#' @return either value if it is not NULL or default
#' 
nvl <- function(value, default) {
    if (is.null(value)) {
        return(default)
    }

    value
}


#' Convert value to numeric if possible and return it, otherwise return default.
#'
#' @param value the value to check
#' @param default the default value to return
#' 
#' @return either value if it is not NULL and numeric or default
#' 
nvlInteger <- function(value, default) {
    tryCatch({
        n <- as.numeric(value)
        if ((n %% 1) == 0) {
            return(n)
        } else {
            return(default)
        }
    },
    error = function(cond) {
        return(default)
    },
    warning = function(cond) {
        return(default)
    },
    finally={
        # nothing
    })
}


#' Check dataset to see if the datatype value is "phenotype".
#'
#' @param dataset a list containg all values for a data set
#' 
#' @return TRUE if the datatype is phenotype, FALSE otherwise
#' 
isPheno <- function(dataset) {
    if (startsWith(tolower(dataset$datatype), "pheno")) {
        return(TRUE)
    }

    FALSE
}


# #############################################################################
#
# Data methods
#
# #############################################################################


#' Get the dataset by id (a string).
#' 
#' @param id the dataset id as a string, either 'dataset.name' or just 'name'
#' 
#' @return the dataset object
#' 
GetDataSet <- function(id) {
    if (exists(id)) {
        return(get(id))
    } else {
        expanded <-paste0('dataset.', id) 
        if (exists(expanded)) {
            return(get(expanded))
        }
    }

    NULL
}

#' Get the dataset by id (a string).
#' 
#' @param id the dataset id as a string
#' @param data the subset of the dataset's data list by id string
#' 
#' @return the data object
#' 
GetData <- function(id, data=NULL) {
    ds <- GetDataSet(id)
    
    if (gtools::invalid(ds)) {
        stop(paste0("dataset not found '", id, '"'))
    }
    
    ret <- NULL
    
    if (is.matrix(ds$data) && (!is.null(data))) {
        stop(paste0("specified data '", data, "' not found in '", id, "'"))
    }
    
    # rz, norm, raw, log, transformed
    if (!is.null(data)) {
        if (is.matrix(ds$data)) {
            stop(paste0("only 1 data element in '", id, "'"))
        }
        
        ret <- ds$data[[data]]
    } else {
        if (is.matrix(ds$data)) {
            ret <- ds$data
        } else {
            if (!is.null(ds$data$rz)) {
                ret <- ds$data$rz
            } else if (!is.null(ds$data$norm)) {
                ret <- ds$data$norm
            } else if (!is.null(ds$data$raw)) {
                ret <- ds$data$raw
            } else if (!is.null(ds$data$log)) {
                ret <- ds$data$log
            } else if (!is.null(ds$data$transformed)) {
                ret <- ds$data$transformed
            }
        }
    }
    
    ret
}

#' Get all "dataset.*" information.
#' 
#' list(dataSets =, ensemblVersion =)
#'
#' @return a list of all the dataset objects, along with the ensembl.version 
#' and number of markers
#' 
GetDatasetInfo <- function() {
    datasets <- grep("^dataset*", apropos("dataset\\."), value = TRUE)
    ret <- c()

    for (d in datasets) {
        ds <- get(d)

        annotations <- list()

        if (tolower(ds$datatype) == "mrna") {
            annotations <- list(ids = ds$annot.mrna$gene.id)
        } else if(tolower(ds$datatype) == "protein") {
            annotations <- list(ids = tibble(protein.id = ds$annot.protein$protein.id, 
                                             gene.id    = ds$annot.protein$gene.id))
        } else if(isPheno(ds)) {
            annotations <- as_tibble(ds$annot.pheno[which(ds$annot.pheno$omit == FALSE), ])
        }

        temp <- list(id              = d, 
                     annotations     = annotations, 
                     display.name    = nvl(ds$display.name, d), 
                     datatype        = ds$datatype, 
                     covar.info      = ds$covar.info)

        ret <- c(ret, list(temp))
    }

    list(datasets        = ret, 
         ensembl.version = ensembl.version)
}


#' Perform the LOD scan.
#' 
#' @param dataset the dataset identifier
#' @param id the identifier
#' @param interactive.covar the interactive covariate
#' @param cores number of cores to use (0=ALL)
#' 
#' @return a tibble with the following columns: id, chr, pos, lod
#' 
PerformLODScan <- function(dataset, id, interactive.covar = NULL, cores = 0) {
    # get the dataset
    ds = GetDataSet(dataset)

    if (gtools::invalid(ds)) {
        stop(paste0("dataset not found '", dataset, "'"))
    }

    # get the data
    data <- GetData(dataset)
    
    if (gtools::invalid(data)) {
        stop(paste0("data not found for '", dataset, "'"))
    }
    
    # get the index
    idx <- which(colnames(data) == id)

    if (gtools::invalid(idx)) {
        stop(paste0("id not found '", id, "'"))
    }

    # make sure num.cores is appropriate  
    num.cores = nvlInteger(cores, 0)
 
    # set the interactive.covariates, to be used in scan1
    # as scan1(intcovar=interactive.covariate)
    interactive.covariate <- NULL

    if (isPheno(ds)) {
        # get the use.covar variable from the annotations
        i <- which(ds$annot.pheno$data.name == id)
        
        if (gtools::invalid(i)) {
            stop(paste0("id not found '", id, "' in annot.pheno"))
        }
        
        covar.str <- strsplit(ds$annot.pheno$use.covar[i], ":")[[1]]
        covar.str <- paste0("~", paste0(covar.str, collapse = "+"))
        
        # find the mouse.id column and use that to set the rownames
        # this has to be done for the checks in scan1
        temp.samples <- as.data.frame(ds$annot.samples)
        id.column <- match('mouse.id', tolower(colnames(ds$annot.samples)))
        rownames(temp.samples) <- ds$annot.samples[[id.column]]
        covar <- model.matrix(as.formula(covar.str), data = temp.samples)[, -1, drop = FALSE]
    } else {
        covar <- NULL
        
        if (!gtools::invalid(ds$covar.matrix)) {
            covar <- ds$covar.matrix    
        }
    }
    
    if (!gtools::invalid(interactive.covar)) {
        if (interactive.covar %not in% ds$covar.info$sample.column) {
            stop(sprintf("covar: %s not found in %s$covar.info", interactive.covar, dataset))
        }

        n <- ds$covar.info[ds$covar.info$sample.column == interactive.covar, ]
        interactive.covariate <- ds$covar.matrix[, n$covar.column, drop = FALSE]
    }
    
    # perform the scan using QTL2, 
    # - addcovar should always be ALL covars
    # - intcovar should be just the interactive covariate column
    temp <- scan1(genoprobs = genoprobs,
                  kinship   = K,
                  pheno     = data[, idx, drop = FALSE], 
                  addcovar  = covar, 
                  intcovar  = interactive.covariate,
                  cores     = num.cores,
                  reml      = TRUE)

    # construct a 2 dimensional array of data with id, chr, pos, lod as columns
    tempDT <- as_tibble(data.table(id  = markers$marker.id,
                                   chr = markers$chr,
                                   pos = markers$pos,
                                   temp))

    # setting colnames to NULL removes the names in the JSON and return an array
    colnames(tempDT)[4] <- "lod"
    
    tempDT
}


#' Perform the LOD scan for each "value" of the interactive covariate.
#' 
#' @param dataset the dataset identifier
#' @param id the identifier
#' @param interactive.covar the interactive covariate, sample.column in covar.info
#' @param chrom the chromosome
#' @param cores number of cores to use (0=ALL)
#' 
#' @return a tibble with the following columns: id, chr, pos, lod
#' 
PerformLODScanBySample <- function(dataset, id, chrom, interactive.covar,
                                   cores = 0) {
    # get the dataset
    ds = GetDataSet(dataset)

    if (gtools::invalid(ds)) {
        stop(paste0("dataset not found '", dataset, "'"))
    }

    # get the data
    data <- GetData(dataset)
    
    if (gtools::invalid(data)) {
        stop(paste0("data not found for '", dataset, "'"))
    }

    # get the index
    idx <- which(colnames(data) == id)
    
    if (gtools::invalid(idx)) {
        stop(paste0("id not found '", id, "'"))
    }
    
    # make sure the chromosome data exists
    if (gtools::invalid(K[[chrom]])) {
        stop(paste0("chrom not found: ", chrom))
    }
    
    # make sure nCores is appropriate  
    num.cores = nvlInteger(cores, 0)
 
    # create a list with each element is indexed by the covar and the 
    # values are the corresponding unique values in the samples table
    covar.categories <- list()
    
    for (i in 1:nrow(ds$covar.info)) {
        inf <- ds$covar.info[i,]
        stopifnot(!is.null(ds$annot.samples[[inf$sample.column]]))
        
        if (is.factor(ds$annot.samples[[inf$sample.column]])) {
            covar.categories[[inf$sample.column]] <- 
                mixedsort(levels(ds$annot.samples[[inf$sample.column]]))
        } else {
            covar.categories[[inf$sample.column]] <- 
                mixedsort(unique(ds$annot.samples[[inf$sample.column]]))
        }
    }

    if (interactive.covar %not in% ds$covar.info$sample.column) {
        stop(sprintf("covar: %s not found in %s$covar.info", interactive.covar, dataset))
    }

    # convert samples to data.frame
    samples.df <- as.data.frame(ds$annot.samples)

    # extract the index of the column that contains some variation of 'mouse.id'
    # NOTE: this is necessary because rownames are deprecated in tibbles
    id.column <- match('mouse.id', tolower(colnames(samples.df)))
    
    # retrieve the covar.info row for the interactive.covar
    n <- ds$covar.info[ds$covar.info$sample.column == interactive.covar, ]
    
    # retrieve all the samples and the value for interactive.covar
    sampleValues <- samples.df[, n$sample.column, drop = FALSE]

    rownames(samples.df) <- ds$annot.samples[[id.column, ]]
    
    ret <- list()

    # loop through the unique values for the interactive.covar
    for (x in covar.categories[[interactive.covar]]) {
        # samples.idxs and samples.names will contain ONLY the samples that match x
        samples.idxs <- which(samples.df[[n$sample.column]] == x)
        samples.names <- samples.df[samples.df[[n$sample.column]] == x, ][[id.column]]
        
        # set covariates
        if (isPheno(ds)) {
            # get the use.covar variable from the annotations
            i <- which(ds$annot.pheno$data.name == id)
            
            if (gtools::invalid(i)) {
                stop(paste0("id not found '", i, "'"))
            }
            
            covar.str <- strsplit(ds$annot.pheno$use.covar[i], ":")[[1]]
            covar.str <- paste0("~", paste0(covar.str, collapse = "+"))
            
            # subset just the data we need
            covar <- model.matrix.lm(as.formula(covar.str), data = samples.df, na.action=na.pass)
            covar <- covar[samples.names, -1, drop=FALSE]
        } else {
            # subset the covars
            covar <- ds$covar.matrix[samples.names, , drop = FALSE]
        }
        
        # exclude covar columns that contain it's name
        covar <- covar[, -which(grepl(n$covar.column, colnames(covar)))]

        temp <- scan1(genoprobs = genoprobs[samples.names, chrom],
                      kinship   = K[[chrom]][samples.names, samples.idxs],
                      pheno     = data[samples.names, idx, drop = FALSE],
                      addcovar  = covar, 
                      cores     = num.cores,
                      reml      = TRUE)

        temp.markers <- markers[which(markers$chr == chrom), ]

        # construct a 2 dimensional array of data
        tempDT <- as_tibble(data.table(id  = temp.markers$marker.id, 
                                       chr = temp.markers$chr, 
                                       pos = temp.markers$pos, 
                                       temp))

        # setting colnames to NULL removes the names in the JSON and return an array
        colnames(tempDT)[4] <- "lod"
        
        ret[[toString(x)]] <- tempDT
    }
    
    ret
}


#' Get the founder coefficients.
#' 
#' @param dataset the dataset identifier
#' @param id the identifier
#' @param chrom the chromosome
#' @param interactive.covar the interactive covariate, sample.column in covar.info
#' @param blup whether or not to perform BLUP
#' @param center whether or not to center the data
#' @param cores number of cores to use (0=ALL)
#' 
#' @return a data.table with the following columns: id, chr, pos, and A-H
#' 
PerformFounderCoefsScan <- function(dataset, id, chrom, interactive.covar, 
                                    blup = FALSE, center = TRUE, cores = 0) {
    
    # get the dataset
    ds = GetDataSet(dataset)
    
    if (gtools::invalid(ds)) {
        stop(paste0("dataset not found '", dataset, "'"))
    }

    # get the data
    data <- GetData(dataset)
    
    if (gtools::invalid(data)) {
        stop(paste0("data not found for '", dataset, "'"))
    }

    # get the index
    idx <- which(colnames(data) == id)
    
    if (gtools::invalid(idx)) {
        stop(paste0("id not found '", id, "'"))
    }
    
    # make sure the chromosome data exists
    if (gtools::invalid(K[[chrom]])) {
        stop(paste0("chrom not found: ", chrom))
    }

    # make sure ncores is appropriate  
    num.cores = nvlInteger(cores, 0)

    # set covariates
    if (isPheno(ds)) {
        # get the use.covar variable from the annotations
        i <- which(ds$annot.pheno$data.name == id)
        
        if (gtools::invalid(i)) {
            stop(paste0("id not found '", i, "'"))
        }
        
        covar.str <- strsplit(ds$annot.pheno$use.covar[i], ":")[[1]]
        covar.str <- paste0("~", paste0(covar.str, collapse = "+"))

        # find the mouse.id column and use that to set the rownames
        # this has to be done for the checks in scan1
        temp.samples <- as.data.frame(ds$annot.samples)
        id.column <- match('mouse.id', tolower(colnames(ds$annot.samples)))
        rownames(temp.samples) <- ds$annot.samples[[id.column]]
        covar <- model.matrix.lm(as.formula(covar.str), data = temp.samples, na.action=na.pass)[, -1, drop = FALSE]
    } else {
        covar <- NULL
        
        if (!gtools::invalid(ds$covar.matrix)) {
            covar <- ds$covar.matrix    
        }
    }

    ret <- list()

    if (gtools::invalid(interactive.covar)) {
        if (toBoolean(blup)) {
            temp <- scan1blup(genoprobs = genoprobs[, chrom],
                              pheno     = data[, idx, drop = FALSE],
                              kinship   = K[[chrom]],
                              addcovar  = covar,
                              cores     = num.cores)
        } else {
            temp <- scan1coef(genoprobs = genoprobs[, chrom],
                              pheno     = data[, idx, drop = FALSE],
                              kinship   = K[[chrom]],
                              addcovar  = covar)
        }

        if (toBoolean(center)) {
            temp[, LETTERS[1:8]] <- 
                temp[, LETTERS[1:8]] - rowMeans(temp[, LETTERS[1:8]], na.rm = TRUE)
        }

        ret[['additive']] <- as_tibble(data.table(id = names(map[[chrom]]), 
                                                  chr = chrom, 
                                                  pos = map[[chrom]], 
                                                  temp[,LETTERS[1:8]]))
    } else {
        if (interactive.covar %not in% ds$covar.info$sample.column) {
            stop(sprintf("covar: %s not found in %s$covar.info", interactive.covar, dataset))
        }

        # create a list with each element is indexed by the covar and the 
        # values are the corresponding unique values in the samples table
        covar.categories <- list()
        
        for (i in 1:nrow(ds$covar.info)) {
            inf <- ds$covar.info[i,]
            stopifnot(!is.null(ds$annot.samples[[inf$sample.column]]))
            
            if (is.factor(ds$annot.samples[[inf$sample.column]])) {
                covar.categories[[inf$sample.column]] <- 
                    mixedsort(levels(ds$annot.samples[[inf$sample.column]]))
            } else {
                covar.categories[[inf$sample.column]] <- 
                    mixedsort(unique(ds$annot.samples[[inf$sample.column]]))
            }
        }
        
        # retrieve all the samples and the value for interactive.covar
        n <- ds$covar.info[ds$covar.info$sample.column == interactive.covar, ]

        # loop through interactive.covar's unique values
        for (x in covar.categories[[interactive.covar]]) {
            # samplesIdx will contain ONLY the samples that match x

            temp.samples <- as.data.frame(ds$annot.samples)
            id.column <- match('mouse.id', tolower(colnames(ds$annot.samples)))
            rownames(temp.samples) <- ds$annot.samples[[id.column]]
            samples.idxs <- temp.samples[temp.samples[[n$sample.column]] == x, ][[id.column]]

            # subset the covars
            #covar.subset <- covar[samples.idxs, , drop = FALSE]
            covar.subset <- covar#[samples.idxs, , drop = FALSE]

            # exclude covar columns that contain it's name
            covar.subset <- covar.subset[, -which(grepl(n$covar.column, colnames(covar.subset)))]

            if (toBoolean(blup)) {
                temp <- scan1blup(genoprobs = genoprobs[samples.idxs, chrom],
                                  pheno     = data[samples.idxs, idx, drop = FALSE],
                                  kinship   = K[[chrom]][samples.idxs, samples.idxs],
                                  addcovar  = covar.subset[samples.idxs, , drop = FALSE],
                                  cores     = numCores)
            } else {
                temp <- scan1coef(genoprobs = genoprobs[samples.idxs, chrom],
                                  pheno     = data[samples.idxs, idx, drop = FALSE],
                                  kinship   = K[[chrom]][samples.idxs, samples.idxs],
                                  addcovar  = covar.subset[samples.idxs, , drop = FALSE])
            }

            if (toBoolean(center)) {
                temp[, LETTERS[1:8]] <- 
                    temp[, LETTERS[1:8]] - rowMeans(temp[, LETTERS[1:8]], na.rm = TRUE)
            }

            ret[[toString(x)]] <- as_tibble(data.table(id = names(map[[chrom]]), 
                                                       chr = chrom, 
                                                       pos = map[[chrom]], 
                                                       temp[,LETTERS[1:8]]))
        }
    }

    ret
}


#' Get the founder coefficients.
#' 
#' @param dataset the dataset identifier
#' @param id the identifier
#' 
#' @return a list with elements of data.frame and list of the datatypes
#' 
GetExpression <- function(dataset, id) {
    # get the dataset
    ds = GetDataSet(dataset)

    if (gtools::invalid(ds)) {
        stop(paste0("dataset not found '", dataset, "'"))
    }

    # get the data
    data <- GetData(dataset)
    
    if (gtools::invalid(data)) {
        stop(paste0("data not found for '", dataset, "'"))
    }
    
    # get the index
    idx <- which(colnames(data) == id)
    
    if (gtools::invalid(idx)) {
        stop(paste0("id not found '", id, "'"))
    }
    
    if (gtools::invalid(ds$covar.info)) {
        t <- NULL
    } else {
        t <- list()
        for (i in ds$covar.info$sample.column) {
            stopifnot(!is.null(ds$annot.samples[[i]]))
            if (is.factor(ds$annot.samples[[i]])) {
                t[[i]] <- mixedsort(levels(ds$annot.samples[[i]]))
            } else {
                t[[i]] <- mixedsort(unique(ds$annot.samples[[i]]))
            }
        }
    }
    
    # make sure 'mouse.id' is lowercase when passed back
    id.column <- match('mouse.id', tolower(colnames(ds$annot.samples)))
    samples <- ds$annot.samples
    colnames(samples)[id.column] <- 'mouse.id'

    output <- cbind(samples, expression=data[, idx])

    # eliminate the _row column down line for JSON
    rownames(output) <- NULL

    list(data = output, 
         datatypes = t)
}    


#' Get the mediation data.
#' 
#' @param dataset the dataset identifier
#' @param id the identifier
#' @param marker.id the marker identifier
#' @param dataset.mediate the dataset to mediate against
#' 
#' @return a data.frame with the following columns depending on datatype: 
#'         mRNA = gene_id, symbol, chr, pos, LOD
#'         protein = protein_id, gene_id, symbol, chr, pos, LOD
#'         phenotype = 
#'         
PerformMediation <- function(dataset, id, marker.id, dataset.mediate = NULL) {
    # get the dataset
    ds = GetDataSet(dataset)

    if (gtools::invalid(ds)) {
        stop(paste0("dataset not found '", dataset, "'"))
    }

    # get the data
    data <- GetData(dataset)
    
    if (gtools::invalid(data)) {
        stop(paste0("data not found for '", dataset, "'"))
    }
    
    # get the index
    idx <- which(colnames(data) == id)
    
    if (gtools::invalid(idx)) {
        stop(paste0("id not found '", id, "'"))
    }
    
    # get the dataset we are mediating against
    dataset.mediate <- nvl(dataset.mediate, dataset)
    ds.mediate = GetDataSet(dataset.mediate)
    
    if (gtools::invalid(ds.mediate)) {
        stop(paste0("dataset not found '", dataset.mediate, "'"))
    }
    
    # get the data
    data.mediate <- GetData(dataset.mediate)
    
    if (gtools::invalid(data.mediate)) {
        stop(paste0("data not found for '", dataset.mediate, "'"))
    }
    
    # get the annotations
    if (tolower(ds.mediate$datatype) == "mrna") {
        annot <- ds.mediate$annot.mrna[, c("gene.id", "symbol", "chr")]
        # must set middle_point for mediation
        annot$middle_point <- (ds.mediate$annot.mrna$start + ds.mediate$annot.mrna$end) / 2.0
    } else if (tolower(ds.mediate$datatype) == "protein") {
        annot <- ds.mediate$annot.protein[, c("protein.id", "gene.id", "symbol", "chr")]
        # must set middle_point for mediation
        annot$middle_point <- (ds.mediate$annot.protein$start + ds.mediate$annot.protein$end) / 2.0
    } else if (isPheno(ds.mediate)) {
        stop("invalid datatype to mediate against")
    } else {
        stop("invalid datatype")
    }
    
    # get the marker index 
    mrkx <- which(markers$marker.id == marker.id)

    if (gtools::invalid(mrkx)) {
        stop(paste0("marker.id not found '", marker.id, "'"))
    }

    chr.tmp = as.character(markers[mrkx, 'chr'])
    
    if (!gtools::invalid(ds.mediate$covar.matrix)) {
        # perform the mediation
        toReturn <- mediation.scan(target     = data[,idx, drop = FALSE],
                                   mediator   = data.mediate,
                                   annotation = annot,
                                   covar      = ds.mediate$covar.matrix,
                                   qtl.geno   = genoprobs[[chr.tmp]][rownames(data.mediate), , marker.id],
                                   verbose    = FALSE)
    } else {
        toReturn <- mediation.scan(target     = data[,idx, drop = FALSE],
                                   mediator   = data.mediate,
                                   annotation = annot,
                                   qtl.geno   = genoprobs[[chr.tmp]][rownames(data.mediate), , marker.id],
                                   verbose    = FALSE)
    }
    
    rownames(toReturn) <- NULL

    toReturn
}


#' Get the SNP association mapping
#' 
#' @param dataset the dataset identifier
#' @param id the identifier
#' @param chrom the chromosome
#' @param location location on chromosome in base pairs
#' @param window.size the size of the window to scan before and after location
#' @param cores number of cores to use (0=ALL)
#' 
#' @return a data.frame with the following columns: 
#' snp, chr, pos, alleles, sdp, ensembl_gene, csq, index, interval, on_map, lod      
#' 
PerformSNPAssocMapping <- function(dataset, id, chrom, location, window.size=500000,
                                   cores=0) {
    # get the dataset
    ds = GetDataSet(dataset)

    if (gtools::invalid(ds)) {
        stop(paste0("dataset not found '", dataset, "'"))
    }

    # get the data
    data <- GetData(dataset)
    
    if (gtools::invalid(data)) {
        stop(paste0("data not found for '", dataset, '"'))
    }
    
    # get the index
    idx <- which(colnames(data) == id)
    
    if (gtools::invalid(idx)) {
        stop(paste0("id not found '", id, "'"))
    }
    
    sel.chr <- chrom

    loc <- nvlInteger(location, -1)
    if (loc == -1) {
        stop(paste0("location is invalid: ", location))
    }

    pos <- loc / 1000000.0

    tmp <- nvlInteger(window.size, -1)
    if (tmp == -1) {
        stop(paste0("window.size is invalid: ", tmp))
    }
    window.size <- tmp

    window.length <- window.size / 1000000.0
    window.range <- pos + c(-1,1) * window.length

    # make sure ncores is appropriate  
    num.cores = nvlInteger(nCores, 0)
    
    have.SNPS <- FALSE
    tries <- 1

    # extract SNPs from the database
    while (!have.SNPS) {
        myDB <- src_sqlite(db.file, create = FALSE)
        window.snps <- tbl(myDB, sql("SELECT * FROM snps")) %>%
            filter(chr     == sel.chr, 
                   pos_Mbp >= window.range[1], 
                   pos_Mbp <= window.range[2]) %>%
            arrange(pos_Mbp) %>%
            collect(n = Inf)
        
        if (NROW(window.snps) > 0) {
            have.SNPS <- TRUE
        } else {
            window.size <- window.size + nvlInteger(window.size, -1)
            window.length <- window.size / 1000000.0
            window.range <- pos + c(-1,1) * window.length
            tries <- tries + 1

            if (tries > 10) {
                stop(paste0("Cannot find snps in region: ", window.range))
            }
        }
    }

    colnames(window.snps)[c(1,3)] = c("snp", "pos")
    window.snps = index_snps(map = map, window.snps)

    # set covariates
    if (isPheno(ds)) {
        # get the use.covar variable from the annotations
        i <- which(ds$annot.pheno$data.name == id)
        
        if (gtools::invalid(i)) {
            stop(paste0("id not found '", i, "'"))
        }
        
        covar.str <- strsplit(ds$annot.pheno$use.covar[i], ":")[[1]]
        covar.str <- paste0("~", paste0(covar.str, collapse = "+"))
        
        # find the mouse.id column and use that to set the rownames
        # this has to be done for the checks in scan1
        temp.samples <- as.data.frame(ds$annot.samples)
        id.column <- match('mouse.id', tolower(colnames(ds$annot.samples)))
        rownames(temp.samples) <- ds$annot.samples[[id.column]]
        covar <- model.matrix(as.formula(covar.str), data = temp.samples)[, -1, drop = FALSE]
        
    } else {
        covar <- NULL
        
        if (!gtools::invalid(ds$covar.matrix)) {
            covar <- ds$covar.matrix    
        }
    }

    # convert allele probs to SNP probs
    snppr <- genoprob_to_snpprob(genoprobs, window.snps)

    # finally the scan
    if (isPheno(ds)) {
        out.snps <- scan1(pheno     = data[, idx, drop = FALSE], 
                          kinship   = K[[sel.chr]], 
                          genoprobs = snppr, 
                          addcovar  = covar, 
                          cores     = num.cores)
    } else {
        out.snps <- scan1(pheno     = data[, idx, drop = FALSE], 
                          kinship   = K[[sel.chr]], 
                          genoprobs = snppr, 
                          addcovar  = covar, 
                          cores     = num.cores)
    }
    
    map.tmp <- qtl2:::snpinfo_to_map(window.snps)
    tmp <- qtl2:::expand_snp_results(out.snps, map.tmp, window.snps)

    ret <- window.snps
    ret$lod <- tmp$lod[, 1]
    
    ret
}


#' Get the LOD peaks
#' 
#' @param dataset the dataset identifier
#' @param interactive.covariate the interactive covarariate, if NULL than 
#' the additive is used
#' 
#' @return a data.frame with the following columns: marker, chr, pos
#' and depending upons dataset$dataType the following columns:
#' mRNA = gene_id, symbol, gene_chrom, middle, lod
#' protein = protein_id, gene_id, symbol, gene_chrom, middle, lod
#' phenotype = data_name, short_name, description, lod
#' 
GetLODPeaks <- function(dataset, interactive.covariate = NULL) {
    # get the dataset
    ds = GetDataSet(dataset)

    if (gtools::invalid(ds)) {
        stop(paste0("dataset not found '", dataset, "'"))
    }

    peaks <- NULL
    
    if (gtools::invalid(interactive.covariate)) {
        peaks <- ds$lod.peaks$additive
    } else {
        # find the covar and get the name of the lod peaks
        if (interactive.covariate %in% ds$covar.info$sample.column) {
            n <- ds$covar.info[ds$covar.info$sample.column == interactive.covariate, ]
            peaks <- ds$lod.peaks[[n$lod.peaks]]  
        } else {
            stop(sprintf("covar: %s not found in %s$covar.info", interactive.covariate, dataset))
        }
    }

    if (gtools::invalid(peaks)) {
        stop(sprintf("no peaks found for covar %s in %s", interactive.covariate, dataset))
    }

    if (tolower(ds$datatype) == "mrna") {
        temp <- merge(x    = ds$annot.mrna[, c("gene.id", "symbol", "chr", "start", "end")], 
                      y    = peaks[, c("gene.id", "marker.id", "lod")], 
                      by.x = "gene.id", 
                      by.y = "gene.id")
        
        temp$gene.pos <- (temp$start + temp$end) / 2.0
        colnames(temp)[3] <- "gene.chr"

        return (merge(x    = markers[, c("marker.id", "chr", "pos")],
                      y    = temp[, c("gene.id", "symbol", "gene.chr", "gene.pos", "lod", "marker.id")], 
                      by.x = "marker.id", 
                      by.y = "marker.id"))
    } else if (tolower(ds$datatype) == "protein") {
        temp <- merge(x    = ds$annot.protein[,c("protein.id", "gene.id", "symbol", "chr", "start", "end")], 
                      y    = peaks[, c("protein.id", "marker.id", "lod")], 
                      by.x = "protein.id", 
                      by.y = "protein.id")

        temp$gene.pos <- (temp$start + temp$end) / 2.0
        colnames(temp)[4] <- "gene.chr"

        return (merge(x    = markers[, c("marker.id", "chr", "pos")],
                      y    = temp[, c("protein.id", "gene.id", "symbol","gene.chr", "gene.pos", "lod", "marker.id")],
                      by.x = "marker.id", 
                      by.y = "marker.id"))                      
    } else if(isPheno(ds)) {
        temp <- merge(x    = ds$annot.pheno[,c("data.name", "short.name", "description")], 
                      y    = peaks[, c("data.name", "marker.id", "lod")], 
                      by.x = "data.name", 
                      by.y = "data.name")
        
        return (merge(x    = markers[, c("marker.id", "chr", "pos")],
                      y    = temp, 
                      by.x = "marker.id", 
                      by.y = "marker.id"))
    } else {
        stop("invalid datatype")
    }
}    


#' Get the LOD peaks
#' 
#' @param dataset the dataset identifier
#' @param interactive.covariate the interactive covarariate, if NULL than 
#' the additive is used
#' 
#' @return a data.frame with the following columns: marker, chr, pos
#' and depending upons dataset$dataType the following columns:
#' mRNA = gene_id, symbol, gene_chrom, middle, lod
#' protein = protein_id, gene_id, symbol, gene_chrom, middle, lod
#' phenotype = data_name, short_name, description, lod
#' 
GetLODPeaksAll <- function(dataset) {
    # get the dataset
    ds = GetDataSet(dataset)
    
    if (gtools::invalid(ds)) {
        stop(paste0("dataset not found '", dataset, "'"))
    }
    
    # get the LOD peaks for each covarint
    additive <- GetLODPeaks(dataset)
    
    peaks <- list(additive = additive)

    for (i in 1:nrow(ds$covar.info)) {
        inf <- ds$covar.info[i, ]
        
        if (inf$interactive) {
            tmp.peaks <- GetLODPeaks(dataset, inf$sample.column)
            peaks[[inf$sample.column]] <- tmp.peaks
        }
    } 
    
    peaks
}    


#' Get the correlation.
#' 
#' @param dataset the dataset identifier
#' @param id the identifier
#' @param dataset.correlate the dataset to correlate to
#' @param max.items maximum number of items
#' 
#' @return a data.frame with the following columns: 
#' 
PerformCorrelation <- function(dataset, id, 
                               dataset.correlate = NULL, max.items = NULL) {
    # get the dataset
    ds <- GetDataSet(dataset)

    if (gtools::invalid(ds)) {
        stop(paste0("dataset not found '", dataset, "'"))
    }

    # get the data
    data <- GetData(dataset)
    
    if (gtools::invalid(data)) {
        stop(paste0("data not found for '", dataset, "'"))
    }
    
    # get the index
    idx <- which(colnames(data) == id)
    
    if (gtools::invalid(idx)) {
        stop(paste0("id not found '", id, "'"))
    }
    
    # get the dataset to correlate to
    dataset.correlate <- nvl(dataset.correlate, dataset)
    ds.correlate = GetDataSet(dataset.correlate)

    if (gtools::invalid(ds.correlate)) {
        stop(paste0("dataset.correlate not found '", dataset.correlate, "'"))
    }
    
    data.correlate <- GetData(dataset.correlate)
    
    if (gtools::invalid(data.correlate)) {
        stop(paste0("data not found for '", dataset.correlate, '"'))
    }
    
    samples <- intersect(rownames(data), rownames(data.correlate))
    data <- data[samples, ]
    data.correlate <- data.correlate[samples, ]    

    pcor <- cor(data[, idx], data.correlate, use = "pair")
    pcor <- pcor[1, order(abs(pcor), decreasing = TRUE)]
    pcor <- pcor[1:nvlInteger(max.items, length(pcor))]

    if (tolower(ds.correlate$datatype) == "mrna") {
        return (as_tibble(data.table(cor    = pcor, 
                           id     = names(pcor), 
                           symbol = ds.correlate$annot.mrna$symbol[match(names(pcor), ds.correlate$annot.mrna$gene.id)],
                           chr    = ds.correlate$annot.mrna$chr[match(names(pcor), ds.correlate$annot.mrna$gene.id)],
                           start  = ds.correlate$annot.mrna$start[match(names(pcor), ds.correlate$annot.mrna$gene.id)],
                           end    = ds.correlate$annot.mrna$end[match(names(pcor), ds.correlate$annot.mrna$gene.id)])))
    } else if (tolower(ds.correlate$datatype) == "protein") {
        return (as_tibble(data.table(cor     = pcor, 
                           id      = names(pcor), 
                           gene_id = ds.correlate$annot.protein$gene.id[match(names(pcor), ds.correlate$annot.protein$protein.id)],
                           symbol  = ds.correlate$annot.protein$symbol[match(names(pcor), ds.correlate$annot.protein$protein.id)],
                           chr     = ds.correlate$annot.protein$chr[match(names(pcor), ds.correlate$annot.protein$protein.id)],
                           start   = ds.correlate$annot.protein$start[match(names(pcor), ds.correlate$annot.protein$protein.id)],
                           end     = ds.correlate$annot.protein$end[match(names(pcor), ds.correlate$annot.protein$protein.id)])))
    } else if (isPheno(ds.correlate)) {
        return (as_tibble(data.table(cor    = pcor, 
                           id     = names(pcor))))
    }
}


#' Get the correlation.
#' 
#' @param dataset the dataset identifier
#' @param id the identifier
#' @param dataset.correlate the dataset to correlate to
#' @param id.correlate the identifier from the correlate dataset
#' 
#' @return a data.frame with the following columns: 
#' 
GetCorrelationPlotData <- function(dataset, id, dataset.correlate, id.correlate) {
    # get the dataset
    ds <- GetDataSet(dataset)
    
    if (gtools::invalid(ds)) {
        stop(paste0("dataset not found '", dataset, "'"))
    }

    # get the data
    data <- GetData(dataset)
    
    if (gtools::invalid(data)) {
        stop(paste0("data not found for '", dataset, "'"))
    }
    
    # get the index
    idx <- which(colnames(data) == id)
    
    if (gtools::invalid(idx)) {
        stop(paste0("id not found '", id, "'"))
    }
    
    ds.correlate <- GetDataSet(dataset.correlate)
    
    if (gtools::invalid(ds.correlate)) {
        stop(paste0("dataset not found '", dataset.correlate, "'"))
    }
    
    # get the data
    data.correlate <- GetData(dataset.correlate)
    
    if (gtools::invalid(data.correlate)) {
        stop(paste0("data not found for '", dataset.correlate, "'"))
    }

    # get the index
    idx.correlate <- which(colnames(data.correlate) == id.correlate)
    
    if (gtools::invalid(idx.correlate)) {
        stop(paste0("id not found: ", id.correlate))
    }

    samples <- intersect(rownames(data), rownames(data.correlate))
    data <- data[samples, ]
    data.correlate <- data.correlate[samples, ]
    
    id.column <- match('mouse.id', tolower(colnames(ds$annot.samples)))
    samples.idx <- which(ds$annot.samples[[id.column]] == samples)
    
    # get the covar factors data
    sample.info <- list()
    dt <- list()
    for (s in ds$covar.info$sample.column) {
        stopifnot(!is.null(ds$annot.samples[[s]]))
        sample.info[[s]] <- ds$annot.samples[samples.idx, ][[s]]

        if (is.factor(ds$annot.samples[[s]])) {
            dt[[s]] <- mixedsort(levels(ds$annot.samples[[s]]))
        } else {
            dt[[s]] <- mixedsort(unique(ds$annot.samples[[s]]))
        }
    }
    
    ret <- list(dataset           = dataset,
                id                = id,                    
                dataset.correlate = dataset.correlate,
                id.correlate      = id.correlate,
                datatypes         = dt,
                data              = as_tibble(data.frame(mouse.id        = rownames(data), 
                                                         x               = data[,idx], 
                                                         y               = data.correlate[, idx.correlate],
                                                         sample.info,
                                                         stringsAsFactors = FALSE)))

    ret
}


