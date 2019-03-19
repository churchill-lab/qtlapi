# #############################################################################
#
# Load the libraries will we need
#
# #############################################################################

library(data.table)
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

if (exists("debugMode")) {
    debugMode <- get("debugMode")
} else {
    debugMode <- FALSE
}

if (debugMode) {
    print('DEBUG MODE: Make sure the data is loaded and db.file is defined')
} else {
    print("Finding the data file...")
    dataFiles <- list.files("/app/qtlapi/data", "*.RData", full.names=TRUE)

    if (length(dataFiles) == 0) {
        stop("There needs to be an .RData file in /app/qtlapi/data")
    } else if (length(dataFiles) > 1) {
        stop("There needs to be only 1 .RData file in /app/qtlapi/data")
    }

    print(paste0("Loading the data file: ", dataFiles))

    load(dataFiles, .GlobalEnv)

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
    if (tolower(dataset$datatype) == "phenotype") {
        return(TRUE)
    }

    FALSE
}


#' Check and subset the dataset to match up with the genoprobs and kinship
#' matrix.
#' 
#' @param dataset a list containg all values for a data set
#' 
#' @return a list containing the subset of genoprobs, data, K, covar
#'
SynchronizeSamples <- function(dataset, genoprobs, K) {
    if (isPheno(dataset)) {
        temp <- dataset$pheno        
    } else {
        temp <- dataset$rankz
    }

    samples <- intersect(rownames(temp), rownames(genoprobs[[1]]))
    samples <- intersect(samples, rownames(K[[1]]))
    samples <- intersect(samples, rownames(dataset$covar))
    samples <- sort(samples)
    
    data <- temp[samples, , drop = FALSE]
    covar <- dataset$covar[samples, , drop = FALSE]

    for(i in 1:length(genoprobs)) {
        genoprobs[[i]] <- genoprobs[[i]][samples, , ]
        K[[i]] <- K[[i]][samples, ]
    }

    list(data     = data, 
        covar     = covar, 
        genoprobs = genoprobs, 
        K         = K)
}


#' Log the request and time.
#' 
#' @param req the request object
#' @param elapsed the elapsed time
#' 
TrackTime <- function(req, elapsed) {
    print(paste0("[", Sys.time(), "] [",
                 req$REMOTE_ADDR, "] [",
                 req$REQUEST_METHOD, "] [",
                 req$PATH_INFO, 
                 req$QUERY_STRING, "] [", 
                 sprintf("%3f", elapsed), "]"))
}


# #############################################################################
#
# Data methods
#
# #############################################################################


#' Get the dataset by id (a string).
#' 
#' @param id the dataset id as a string
#' 
#' @return the dataset object
#' 
GetDataSet <- function(id) {
    if (exists(id)) {
        return(get(id))
    }

    NULL
}

#' Get the interactive covariates for a dataset.
#' 
#' @param id the dataset id as a string
#' 
#' @return the interactive covariates
#' 
GetIntCovar <- function(dataset) {
    # get the dataset
    ds = GetDataSet(dataset)

    if (gtools::invalid(ds)) {
        stop(paste0("dataset not found: ", dataset))
    }
    
    if (gtools::invalid(ds$covar.factors)) {
        return(NULL)
    }

    ds$covar.factors[which(!is.na(ds$covar.factors$covar.name)),]$column.name
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

        if (ds$datatype == "mRNA") {
            annotations <- list(ids = ds$annots$gene_id)
        } else if(ds$datatype == "protein") {
            annotations <- list(ids = data.frame(protein_id = ds$annots$protein_id, 
                                                 gene_id    = ds$annots$gene_id))
        } else if(ds$datatype == "phenotype") {
            annotSubset <- ds$annots[which(ds$annots$omit == FALSE), ]
            annotations <- data.frame(dataName     = annotSubset$data_name,
                                      shortName    = annotSubset$short_name,
                                      desc         = annotSubset$description,
                                      category     = annotSubset$category,
                                      units        = annotSubset$units,
                                      isNumeric    = annotSubset$is_numeric,
                                      isDate       = annotSubset$is_date,
                                      isFactor     = annotSubset$is_factor,
                                      isCovar      = annotSubset$is_covar,
                                      isPheno      = annotSubset$is_pheno,
                                      isDerived    = annotSubset$is_derived,
                                      factorLevels = annotSubset$factor_levels,
                                      useCovar     = annotSubset$use_covar)
        }

        dsInfo <- list(id             = d, 
                       annotations    = annotations, 
                       displayName    = nvl(ds$display.name, d), 
                       dataType       = ds$datatype, 
                       ensemblVersion = nvl(ds$ensembl.version, ''),
                       covarFactors   = ds$covar.factors)

        ret <- c(ret, list(dsInfo))
    }

    list(dataSets       = ret, 
         ensemblVersion = ensembl.version)
}


#' Perform the LOD scan.
#' 
#' @param dataset the dataset identifier
#' @param id the identifier
#' @param intCovar the interactive covariate
#' @param nCores number of cores to use (0=ALL)
#' 
#' @return a data.table with the following columns: id, chr, pos, lod
#' 
GetLODScan <- function(dataset, id, intCovar = NULL, nCores = 0) {
    # get the dataset
    ds = GetDataSet(dataset)

    if (gtools::invalid(ds)) {
        stop(paste0("dataset not found: ", dataset))
    }

    # get the index and data 
    if (ds$datatype == "mRNA") {
        idx <- which(ds$annots$gene_id == id)
        data <- ds$rankz
    } else if (ds$datatype == "protein") {
        idx <- which(ds$annots$protein_id == id)
        data <- ds$rankz
    } else if (ds$datatype == "phenotype") {
        idx <- which(ds$annots$data_name == id)
        data <- ds$pheno
    } else {
        stop("invalid datatype")
    }
    
    if (gtools::invalid(idx)) {
        stop(paste0("id not found: ", id))
    }

    # make sure nCores is appropriate  
    numCores = nvlInteger(nCores, 0)
 
    # set the interactiveCovariates, to be used in scan as intcovar
    interactiveCovariate <- NULL

    if (isPheno(ds)) {
        covarStr <- strsplit(ds$annots$use_covar[idx], ":")[[1]]
        covarStr <- paste0("~", paste0(covarStr, collapse = "+"))
        covar <- model.matrix(as.formula(covarStr), data = ds$samples)[, -1, drop = FALSE]
    } else {
        covar <- NULL
        
        if (!gtools::invalid(ds$covar)) {
            covar <- ds$covar    
        }
    }
    
    if (!gtools::invalid(intCovar)) {
        if (intCovar %not in% ds$covar.factors$column.name) {
            stop(sprintf("covar: %s not found in %s$covar.factors", intCovar, dataset))
        }

        n <- ds$covar.factors[ds$covar.factors$column.name == intCovar, ]
        interactiveCovariate <- ds$covar[, n$covar.name, drop = FALSE]
    }
    
    # perform the scan using QTL2, 
    # - addcovar should always be ALL covars
    # - intcovar should be just the intCovar column
    temp <- scan1(genoprobs = genoprobs,
                  kinship   = K,
                  pheno     = data[, idx, drop = FALSE], 
                  addcovar  = covar, 
                  intcovar  = interactiveCovariate,
                  cores     = numCores,
                  reml      = TRUE)

    # construct a 2 dimensional array of data with id, chr, pos, lod as columns
    tempDT <- data.table(id  = markers$marker,
                         chr = markers$chr,
                         pos = markers$pos,
                         temp)

    # setting colnames to NULL removes the names in the JSON and return an array
    colnames(tempDT)[4] <- "lod"
    
    tempDT
}


#' Perform the LOD scan for each "value" of the interactive covariate.
#' 
#' @param dataset the dataset identifier
#' @param id the identifier
#' @param intCovar the interactive covariate
#' @param chrom the chromosome
#' @param nCores number of cores to use (0=ALL)
#' 
#' @return a data.table with the following columns: id, chr, pos, lod
#' 
GetLODScanBySample <- function(dataset, id, intCovar, 
                               chrom = NULL, nCores = 0) {
    # get the dataset
    ds = GetDataSet(dataset)

    if (gtools::invalid(ds)) {
        stop(paste0("dataset not found: ", dataset))
    }

    # get the index and data 
    if (ds$datatype == "mRNA") {
        idx <- which(ds$annots$gene_id == id)
        data <- ds$rankz
    } else if (ds$datatype == "protein") {
        idx <- which(ds$annots$protein_id == id)
        data <- ds$rankz
    } else if (ds$datatype == "phenotype") {
        idx <- which(ds$annots$data_name == id)
        data <- ds$pheno
    } else {
        stop("invalid datatype")
    }

    if (gtools::invalid(idx)) {
        stop(paste0("id not found: ", id))
    }

    # make sure nCores is appropriate  
    numCores = nvlInteger(nCores, 0)
 
    # create a list with each element is indexed by the covar and the 
    # values are the corresponding unique values in the samples table
    covarCategories <- list()
    for (f in ds$covar.factors$column.name) {
        stopifnot(!is.null(ds$samples[[f]]))
        if (is.factor(ds$samples[[f]])) {
            covarCategories[[f]] <- mixedsort(levels(ds$samples[[f]]))
        } else {
            covarCategories[[f]] <- mixedsort(unique(ds$samples[[f]]))
        }
    }

    # retreive all the samples and the value for intCovar
    sampleValues <- ds$samples[, intCovar, drop = FALSE]

    ret <- list()

    # loop through intCovar's unique values
    for (x in covarCategories[[intCovar]]) {
        # samplesIdx will contain ONLY the samples that match x
        samplesIdx <- rownames(sampleValues)[sampleValues[, 1] == x]

        # set covariates
        if (isPheno(ds)) {
            covarStr <- strsplit(ds$annots$use_covar[idx], ":")[[1]]
            covarStr <- paste0("~", paste0(covarStr, collapse = "+"))
            covar <- model.matrix(as.formula(covarStr), data = ds$samples)[samplesIdx, -1, drop = FALSE]
        } else {
            # subset the covars
            covar <- ds$covar[samplesIdx, , drop = FALSE]
        }

        # get the interactive covariate
        n <- ds$covar.factors[ds$covar.factors$column.name == intCovar, ]

        # exclude covar columns that contain it's name
        covar <- covar[, -which(grepl(n$covar.name, colnames(covar)))]

        if (gtools::invalid(chrom)) {
            temp <- (scan1(genoprobs = genoprobs[samplesIdx, ],
                           kinship   = K[samplesIdx, samplesIdx],
                           pheno     = data[samplesIdx, idx, drop = FALSE],
                           addcovar  = covar, 
                           cores     = numCores,
                           reml      = TRUE))
            
            tempMarkers <- markers
        } else {
            temp <- scan1(genoprobs = genoprobs[samplesIdx, chrom],
                          kinship   = K[[chrom]][samplesIdx, samplesIdx],
                          pheno     = data[samplesIdx, idx, drop = FALSE],
                          addcovar  = covar, 
                          cores     = numCores,
                          reml      = TRUE)

            tempMarkers <- markers[which(markers$chr == chrom), ]
        }

        # construct a 2 dimensional array of data
        tempDT <- data.table(id  = tempMarkers$marker, 
                             chr = tempMarkers$chr, 
                             pos = tempMarkers$pos, 
                             temp)

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
#' @param intCovar the interactive covariate
#' @param blup whether or not to perform BLUP
#' @param center whether or not to center the data
#' @param nCores number of cores to use (0=ALL)
#' 
#' @return a data.table with the following columns: id, chr, pos, and A-H
#' 
GetFoundercoefs <- function(dataset, id, chrom, intCovar = NULL, 
                            blup = FALSE, center = TRUE, nCores = 0) {
    # get the dataset
    ds = GetDataSet(dataset)
    
    if (gtools::invalid(ds)) {
        stop(paste0("dataset not found: ", dataset))
    }

    # get the index 
    if (ds$datatype == "mRNA") {
        idx <- which(ds$annots$gene_id == id)
        data <- ds$rankz
    } else if (ds$datatype == "protein") {
        idx <- which(ds$annots$protein_id == id)
        data <- ds$rankz
    } else if (ds$datatype == "phenotype") {
        idx <- which(ds$annots$data_name == id)
        data <- ds$pheno
    } else {
        stop("invalid datatype")
    }

    if (gtools::invalid(idx)) {
        stop(paste0("id not found: ", id))
    }

    # make sure the chromosome data exists
    if (gtools::invalid(K[[chrom]])) {
        stop(paste0("chrom not found: ", chrom))
    }

    # make sure ncores is appropriate  
    numCores = nvlInteger(nCores, 0)

    # set covariates
    if (isPheno(ds)) {
        covarStr <- strsplit(ds$annots$use_covar[idx], ":")[[1]]
        covarStr <- paste0("~", paste0(covarStr, collapse = "+"))
        covar <- model.matrix(as.formula(covarStr), data = ds$pheno)[, -1, drop = FALSE]
    } else {
        covar <- NULL
        
        if (!gtools::invalid(ds$covar)) {
            covar <- ds$covar    
        }
    }

    ret <- list()

    if (gtools::invalid(intCovar)) {
        if (toBoolean(blup)) {
            temp <- scan1blup(genoprobs = genoprobs[, chrom],
                              pheno     = data[, idx, drop = FALSE],
                              kinship   = K[[chrom]],
                              addcovar  = covar,
                              cores     = numCores)
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

        ret[['additive']] <- data.table(id = names(map[[chrom]]), 
                                        chr = chrom, 
                                        pos = map[[chrom]], 
                                        temp[,LETTERS[1:8]])
    } else {
        if (intCovar %not in% ds$covar.factors$column.name) {
            stop(sprintf("covar: %s not found in %s$covar.factors", intCovar, dataset))
        }

        # create a list with each element is indexed by the covar and the
        # values are the corresponding unique values in the samples table
        covarCategories <- list()
        for (f in ds$covar.factors$column.name) {
            stopifnot(!is.null(ds$samples[[f]]))
            if (is.factor(ds$samples[[f]])) {
                covarCategories[[f]] <- mixedsort(levels(ds$samples[[f]]))
            } else {
                covarCategories[[f]] <- mixedsort(unique(ds$samples[[f]]))
            }
        }

        # retrive all the samples and the value for intCovar
        sampleValues <- ds$samples[, intCovar, drop = FALSE]

        # loop through intCovar's unique values
        for (x in covarCategories[[intCovar]]) {
            # samplesIdx will contain ONLY the samples that match x
            samplesIdx <- rownames(sampleValues)[sampleValues[, 1] == x]

            # subset the covars
            subsetCovar <- covar[samplesIdx, , drop = FALSE]

            # get the interactive covariate
            n <- ds$covar.factors[ds$covar.factors$column.name == intCovar, ]

            # exclude covar columns that contain it's name
            subsetCovar <- subsetCovar[, -which(grepl(n$covar.name, colnames(subsetCovar)))]

            if (toBoolean(blup)) {
                temp <- scan1blup(genoprobs = genoprobs[samplesIdx, chrom],
                                  pheno     = data[samplesIdx, idx, drop = FALSE],
                                  kinship   = K[[chrom]][samplesIdx, samplesIdx],
                                  addcovar  = subsetCovar[samplesIdx, , drop = FALSE],
                                  cores     = numCores)
            } else {
                temp <- scan1coef(genoprobs = genoprobs[samplesIdx, chrom],
                                  pheno     = data[samplesIdx, idx, drop = FALSE],
                                  kinship   = K[[chrom]][samplesIdx, samplesIdx],
                                  addcovar  = subsetCovar[samplesIdx, , drop = FALSE])
            }

            if (toBoolean(center)) {
                temp[, LETTERS[1:8]] <- 
                    temp[, LETTERS[1:8]] - rowMeans(temp[, LETTERS[1:8]], na.rm = TRUE)
            }

            ret[[toString(x)]] <- data.table(id = names(map[[chrom]]), 
                                             chr = chrom, 
                                             pos = map[[chrom]], 
                                             temp[,LETTERS[1:8]])
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
        stop(paste0("dataset not found: ", dataset))
    }

    # get the index 
    idx <- 0
    if (ds$datatype == "mRNA") {
        idx <- which(ds$annots$gene_id == id)
    } else if (ds$datatype == "protein") {
        idx <- which(ds$annots$protein_id == id)
    } else if (ds$datatype == "phenotype") {
        idx <- which(ds$annots$data_name == id)
    } else {
        stop("invalid datatype")
    }

    if (gtools::invalid(idx)) {
        stop(paste0("id not found: ", id))
    }
    
    if (gtools::invalid(ds$covar.factor)) {
        t <- NULL
    } else {
        t <- list()
        for (f in ds$covar.factors$column.name) {
            stopifnot(!is.null(ds$samples[[f]]))
            if (is.factor(ds$samples[[f]])) {
                t[[f]] <- mixedsort(levels(ds$samples[[f]]))
            } else {
                t[[f]] <- mixedsort(unique(ds$samples[[f]]))
            }
        }
    }

    if (isPheno(ds)) {
        output <- cbind(ds$samples, expression=ds$pheno[, idx])
    } else {
        output <- cbind(ds$samples, expression=ds$rankz[, idx])
    }

    # rename 'mouse.id' to be 'mouse_id' for easier JSON with JavaScript
    colnames(output)[tolower(colnames(output))=="mouse.id"] <- "mouse_id"

    # eliminate the _row column down line for JSON
    rownames(output) <- NULL

    list(data = output, 
         dataTypes = t)
}    


#' Get the mediation data.
#' 
#' @param dataset the dataset identifier
#' @param id the identifier
#' @param mid the marker identifier
#' @param datasetMediate the dataset to mediate against
#' 
#' @return a data.frame with the following columns depending on datatype: 
#'         mRNA = gene_id, symbol, chr, pos, LOD
#'         protein = protein_id, gene_id, symbol, chr, pos, LOD
#'         phenotype = 
#'         
GetMediate <- function(dataset, id, mid, datasetMediate = NULL) {
    # get the dataset
    ds = GetDataSet(dataset)

    if (gtools::invalid(ds)) {
        stop(paste0("dataset not found: ", dataset))
    }

    # get the index 
    idx <- 0
    annot <- NULL
    if (ds$datatype == "mRNA") {
        idx <- which(ds$annots$gene_id == id)
        data <- ds$rankz
    } else if (ds$datatype == "protein") {
        idx <- which(ds$annots$protein_id == id)
        data <- ds$rankz
    } else if (ds$datatype == "phenotype") {
        idx <- which(ds$annots$data_name == id)
        data <- ds$pheno
    } else {
        stop("invalid datatype")
    }

    if (gtools::invalid(idx)) {
        stop(paste0("id not found: ", id))
    }

    # get the dataset we are mediating against
    datasetMediate <- nvl(datasetMediate, dataset)
    dsMediate = GetDataSet(datasetMediate)
    if (gtools::invalid(ds)) {
        stop(paste0("dataset not found: ", datasetMediate))
    }

    # get the annotations
    if (dsMediate$datatype == "mRNA") {
        annot <- dsMediate$annots[,c("gene_id", "symbol", "chr")]
    } else if (dsMediate$datatype == "protein") {
        annot <- dsMediate$annots[,c("protein_id", "gene_id", "symbol", "chr")]
    } else if (dsMediate$datatype == "phenotype") {
        stop("invalid datatype to mediate against")
    } else {
        stop("invalid datatype")
    }

    # get the marker index 
    mrkx <- which(markers$marker == mid)

    if (gtools::invalid(mrkx)) {
        stop(paste0("mid not found: ", mid))
    }

    chrTmp = as.character(markers[mrkx, 2])
    annot$middle_point <- dsMediate$annots$middle
    
    if (!gtools::invalid(dsMediate$covar)) {
        # perform the mediation
        toReturn <- mediation.scan(target     = data[,idx, drop = FALSE],
                                   mediator   = dsMediate$rankz,
                                   annotation = annot,
                                   covar      = dsMediate$covar,
                                   qtl.geno   = genoprobs[[chrTmp]][rownames(dsMediate$rankz), , mid],
                                   verbose    = FALSE)
    } else {
        toReturn <- mediation.scan(target     = data[,idx, drop = FALSE],
                                   mediator   = dsMediate$rankz,
                                   annotation = annot,
                                   qtl.geno   = genoprobs[[chrTmp]][rownames(dsMediate$rankz), , mid],
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
#' @param windowSize the size of the window to scan before and after location
#' @param nCores number of cores to use (0=ALL)
#' 
#' @return a data.frame with the following columns: 
#' snp, chr, pos, alleles, sdp, ensembl_gene, csq, index, interval, on_map, lod      
#' 
GetSnpAssocMapping <- function(dataset, id, chrom, location, windowSize=500000,
                               nCores=0) {
    # get the dataset
    ds = GetDataSet(dataset)

    if (gtools::invalid(ds)) {
        stop(paste0("dataset not found: ", dataset))
    }

    # get the index 
    if (ds$datatype == "mRNA") {
        idx <- which(ds$annots$gene_id == id)
    } else if (ds$datatype == "protein") {
        idx <- which(ds$annots$protein_id == id)
    } else if (ds$datatype == "phenotype") {
        idx <- which(ds$annots$data_name == id)
    } else {
        stop("invalid datatype")
    }

    if (gtools::invalid(idx)) {
        stop(paste0("id not found: ", id))
    }

    sel.chr <- chrom

    loc <- nvlInteger(location, -1)
    if (loc == -1) {
        stop(paste0("location is invalid: ", location))
    }

    pos <- loc / 1000000.0

    window.size <- nvlInteger(windowSize, -1)
    if (window.size == -1) {
        stop(paste0("window is invalid: ", windowSize))
    }

    window.length <- window.size/1000000.0
    window.range <- pos + c(-1,1)*window.length

    # make sure ncores is appropriate  
    numCores = nvlInteger(nCores, 0)
    
    haveSNPS <- FALSE
    tries <- 1

    # extract SNPs from the database
    while (!haveSNPS) {
        myDB <- src_sqlite(db.file, create = FALSE)
        window.snps <- tbl(myDB, sql("SELECT * FROM snps")) %>%
            filter(chr==sel.chr, 
                   pos_Mbp>=window.range[1], 
                   pos_Mbp<=window.range[2]) %>%
            arrange(pos_Mbp) %>%
            collect(n=Inf)
        
        if (NROW(window.snps) > 0) {
            haveSNPS <- TRUE
        } else {
            window.size <- window.size + nvlInteger(windowSize, -1)
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
        covarStr <- strsplit(ds$annots$use_covar[idx], ":")[[1]]
        covarStr <- paste0("~", paste0(covarStr, collapse = "+"))
        covar <- model.matrix(as.formula(covarStr), data = ds$pheno)[, -1, drop = FALSE]
    } else {
        covar <- NULL
        
        if (!gtools::invalid(ds$covar)) {
            covar <- ds$covar    
        }
    }

    # convert allele probs to SNP probs
    snppr <- genoprob_to_snpprob(genoprobs, window.snps)

    # finally the scan
    if (isPheno(ds)) {
        outSnps <- scan1(pheno     = ds$pheno[, idx, drop = FALSE], 
                         kinship   = K[[sel.chr]], 
                         genoprobs = snppr, 
                         addcovar  = covar, 
                         cores     = numCores)
    } else {
        outSnps <- scan1(pheno     = ds$rankz[, idx, drop = FALSE], 
                         kinship   = K[[sel.chr]], 
                         genoprobs = snppr, 
                         addcovar  = covar, 
                         cores     = numCores)
    }
    
    mapTmp <- qtl2:::snpinfo_to_map(window.snps)
    tmp <- qtl2:::expand_snp_results(outSnps, mapTmp, window.snps)

    toReturn <- window.snps
    toReturn$lod <- tmp$lod[, 1]
    
    toReturn
}


#' Get the LOD peaks
#' 
#' @param dataset the dataset identifier
#' @param intCovar the interactive covarariate, if NULL than 
#' the additive is used
#' 
#' @return a data.frame with the following columns: marker, chr, pos
#' and depending upons dataset$dataType the following columns:
#' mRNA = gene_id, symbol, gene_chrom, middle, lod
#' protein = protein_id, gene_id, symbol, gene_chrom, middle, lod
#' phenotype = data_name, short_name, description, lod
#' 
GetLODPeaks <- function(dataset, intCovar = NULL) {
    # get the dataset
    ds = GetDataSet(dataset)

    if (gtools::invalid(ds)) {
        stop(paste0("dataset not found: ", dataset))
    }

    peaks <- NULL
    
    if (gtools::invalid(intCovar)) {
        peaks <- ds$lod.peaks$additive
    } else {
        # find the covar and get the name of the lod peaks
        if (intCovar %in% ds$covar.factors$column.name) {
            n <- ds$covar.factors[ds$covar.factors$column.name == intCovar,]
            peaks <- ds$lod.peaks[[n$lod.peaks]]  
        } else {
            stop(sprintf("covar: %s not found in %s$covar.factors", intCovar, dataset))
        }
    }

    if (gtools::invalid(peaks)) {
        stop(sprintf("no peaks found for covar %s in %s", intCovar, dataset))
    }

    if (ds$datatype == "mRNA") {
        temp <- merge(x    = ds$annots[,c("gene_id", "symbol", "chr", "middle")], 
                      y    = peaks[, c("annot.id", "marker.id", "lod")], 
                      by.x = "gene_id", 
                      by.y = "annot.id")

        colnames(temp)[3] <- "gene_chrom"

        return (merge(x    = markers[, c("marker", "chr", "pos")],
                      y    = temp, 
                      by.x = "marker", 
                      by.y = "marker.id"))                      
    } else if (ds$datatype == "protein") {
        temp <- merge(x    = ds$annots[,c("protein_id", "gene_id", "symbol", "chr", "middle")], 
                      y    = peaks[, c("annot.id", "marker.id", "lod")], 
                      by.x = "protein_id", 
                      by.y = "annot.id")

        colnames(temp)[4] <- "gene_chrom"

        return (merge(x    = markers[, c("marker", "chr", "pos")],
                      y    = temp, 
                      by.x = "marker", 
                      by.y = "marker.id"))                      
    } else if (ds$datatype == "phenotype") {
        temp <- merge(x    = ds$annots[,c("data_name", "short_name", "description")], 
                      y    = peaks[, c("annot.id", "marker.id", "lod")], 
                      by.x = "data_name", 
                      by.y = "annot.id")

        return (merge(x    = markers[, c("marker", "chr", "pos")],
                      y    = temp, 
                      by.x = "marker", 
                      by.y = "marker.id"))
    } else {
        stop("invalid datatype")
    }
}    


#' Get the correlation.
#' 
#' @param dataset the dataset identifier
#' @param id the identifier
#' @param datasetCorrelate the dataset to correlate to
#' @param maxItems maximum number of items
#' 
#' @return a data.frame with the following columns: 
#' 
GetCorrelation <- function(dataset, id, 
                           datasetCorrelate = NULL, maxItems = NULL) {
    # get the dataset
    ds <- GetDataSet(dataset)

    if (gtools::invalid(ds)) {
        stop(paste0("dataset not found: ", dataset))
    }

    data <- NULL

    # get the index 
    if (ds$datatype == "mRNA") {
        idx <- which(ds$annots$gene_id == id)
        data <- as.matrix(ds$rankz) 
    } else if (ds$datatype == "protein") {
        idx <- which(ds$annots$protein_id == id)
        data <- as.matrix(ds$rankz) 
    } else if (ds$datatype == "phenotype") {
        idx <- which(colnames(ds$pheno[,ds$annots$is_pheno == TRUE]) == id)
        data <- as.matrix(ds$pheno[,ds$annots$is_pheno == TRUE]) 
    } else {
        stop("invalid datatype")
    }

    if (gtools::invalid(idx)) {
        stop(paste0("id not found: ", id))
    }

    # get the dataset to correlate to
    datasetCorrelate <- nvl(datasetCorrelate, dataset)

    dsCorrelate = GetDataSet(datasetCorrelate)
    if (gtools::invalid(dsCorrelate)) {
        stop(paste0("datasetCorrelate not found: ", dataset))
    }

    dataCorrelate <- NULL

    if (dsCorrelate$datatype == "mRNA") {
        dataCorrelate <- as.matrix(dsCorrelate$rankz)
    } else if (dsCorrelate$datatype == "protein") {
        dataCorrelate <- as.matrix(dsCorrelate$rankz)
    } else if (dsCorrelate$datatype == "phenotype") {
        dataCorrelate <- as.matrix(dsCorrelate$pheno[,dsCorrelate$annots$is_pheno == TRUE]) 
    } else {
        stop("invalid datatype for datasetCorrelate")
    }

    samples <- intersect(rownames(data), rownames(dataCorrelate))
    data <- data[samples,]
    dataCorrelate <- dataCorrelate[samples,]    

    pcor <- cor(data[,idx], dataCorrelate, use = "pair")
    pcor <- pcor[1, order(abs(pcor), decreasing = TRUE)]
    pcor <- pcor[1:nvlInteger(maxItems, length(pcor))]

    if (dsCorrelate$datatype == "mRNA") {
        return (data.table(cor    = pcor, 
                           id     = names(pcor), 
                           symbol = dsCorrelate$annots$symbol[match(names(pcor), dsCorrelate$annots$gene_id)],
                           chr    = dsCorrelate$annots$chr[match(names(pcor), dsCorrelate$annots$gene_id)],
                           start  = dsCorrelate$annots$start[match(names(pcor), dsCorrelate$annots$gene_id)],
                           end    = dsCorrelate$annots$end[match(names(pcor), dsCorrelate$annots$gene_id)]))
    } else if (dsCorrelate$datatype == "protein") {
        return (data.table(cor     = pcor, 
                           id      = names(pcor), 
                           gene_id = dsCorrelate$annots$gene_id[match(names(pcor), dsCorrelate$annots$protein_id)],
                           symbol  = dsCorrelate$annots$symbol[match(names(pcor), dsCorrelate$annots$protein_id)],
                           chr     = dsCorrelate$annots$chr[match(names(pcor), dsCorrelate$annots$protein_id)],
                           start   = dsCorrelate$annots$start[match(names(pcor), dsCorrelate$annots$protein_id)],
                           end     = dsCorrelate$annots$end[match(names(pcor), dsCorrelate$annots$protein_id)]))
    } else if (dsCorrelate$datatype == "phenotype") {
        return (data.table(cor    = pcor, 
                           id     = names(pcor)))
    }
}


#' Get the correlation.
#' 
#' @param dataset the dataset identifier
#' @param id the identifier
#' @param datasetCorrelate the dataset to correlate to
#' @param idCorrelate the identifier from the correlate dataset
#' 
#' @return a data.frame with the following columns: 
#' 
GetCorrelationPlotData <- function(dataset, id, datasetCorrelate, idCorrelate) {
    # get the dataset
    ds <- GetDataSet(dataset)
    
    if (gtools::invalid(ds)) {
        stop(paste0("dataset not found: ", dataset))
    }

    dsCorrelate <- GetDataSet(datasetCorrelate)
    
    if (gtools::invalid(dsCorrelate)) {
        stop(paste0("dataset not found: ", datasetCorrelate))
    }

    # get the index 
    if (ds$datatype == "mRNA") {
        idx <- which(ds$annots$gene_id == id)
        data <- as.matrix(ds$rankz) 
    } else if (ds$datatype == "protein") {
        idx <- which(ds$annots$protein_id == id)
        data <- as.matrix(ds$rankz) 
    } else if (ds$datatype == "phenotype") {
        idx <- which(colnames(ds$pheno[,ds$annots$is_pheno == TRUE]) == id)
        data <- as.matrix(ds$pheno[,ds$annots$is_pheno == TRUE]) 
    } else {
        stop("invalid datatype")
    }

    if (gtools::invalid(idx)) {
        stop(paste0("id not found: ", id))
    }

    # get the index of the correlate
    if (dsCorrelate$datatype == "mRNA") {
        idxCorrelate <- which(dsCorrelate$annots$gene_id == idCorrelate)
        dataCorrelate <- as.matrix(dsCorrelate$rankz) 
    } else if (dsCorrelate$datatype == "protein") {
        idxCorrelate <- which(dsCorrelate$annots$protein_id == idCorrelate)
        dataCorrelate <- as.matrix(dsCorrelate$rankz) 
    } else if (dsCorrelate$datatype == "phenotype") {
        idxCorrelate <- which(colnames(dsCorrelate$pheno[,dsCorrelate$annots$is_pheno == TRUE]) == idCorrelate)
        dataCorrelate <- as.matrix(dsCorrelate$pheno[,dsCorrelate$annots$is_pheno == TRUE]) 
    } else {
        stop("invalid datatype")
    }

    if (gtools::invalid(idxCorrelate)) {
        stop(paste0("id not found: ", idCorrelate))
    }

    samples <- intersect(rownames(data), rownames(dataCorrelate))
    data <- data[samples,]
    dataCorrelate <- dataCorrelate[samples,]

    # get the covar factors data
    sampleInfo <- list()
    dt <- list()
    for (f in ds$covar.factors$column.name) {
        stopifnot(!is.null(ds$samples[[f]]))
        sampleInfo[[f]] <- ds$samples[samples,][[f]]

        if (is.factor(ds$samples[[f]])) {
            dt[[f]] <- mixedsort(levels(ds$samples[[f]]))
        } else {
            dt[[f]] <- mixedsort(unique(ds$samples[[f]]))
        }
    }

    toReturn <- list(dataset          = dataset,
                     id               = id,                    
                     datasetCorrelate = datasetCorrelate,
                     idCorrelate      = idCorrelate,
                     dataTypes        = dt,
                     data             = data.frame(mouse_id        = rownames(data), 
                                                   x               = data[,idx], 
                                                   y               = dataCorrelate[,idxCorrelate],
                                                   sampleInfo,
                                                   stringsAsFactors = FALSE))

    toReturn
}


#' Generate an error response 
#' 
#' @param res the response object
#' @param code the error code, defaults to 400
#' @param message the error message, defaults to "Error"
#'
#' @return a response object that is JSON data containing error information
#'
SetError <- function(res, code = 400, message = "Error") {
    res$status <- code
    res$body <- toJSON(list(error = message), auto_unbox = TRUE)
    res
}


#' Set request header "Access-Control-Allow-Origin" to "*".
#' 
#* @filter access_control
function(req, res){
    res$setHeader("Access-Control-Allow-Origin", "*")
    plumber::forward()
}

