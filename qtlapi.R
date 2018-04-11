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
library(plumber)
library(gtools)


# #############################################################################
#
# Load the main data file
#
# #############################################################################

print("Finding the data file...")

debug <- FALSE

if (debug) {
    print('DEBUG MODE: Make sure datafile is loaded and db.file is defined')
} else {
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
        return (value)
    } else if (is.character(value)) {
        return (toupper(value) %in% c("T", "TRUE", "YES", "Y", "1"))
    } else if (is.numeric(value)) {
        return (value == 1)
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
        return (default)
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
    tryCatch(
        {
            n <- as.numeric(value)
            if ((n %% 1) == 0) {
                return (n)
            } else {
                return (default)
            }
        },
        error = function(cond) {
            return (default)
        },
        warning = function(cond) {
            return (default)
        },
        finally={
            # nothing
        }
    )
}


#' Check dataset to see if the datatype value is "phenotype".
#'
#' @param dataset a list containg all values for a data set
#' 
#' @return TRUE if the datatype is phenotype, FALSE otherwise
#' 
isPheno <- function(dataset) {
    if (dataset$datatype == "phenotype") {
        return (TRUE)
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
        temp <- dataset$expr
    }
    
    samples <- intersect(rownames(temp), rownames(genoprobs[[1]]))
    samples <- intersect(samples, rownames(K[[1]]))
    samples <- intersect(samples, rownames(dataset$covar))
    samples <- sort(samples)
    
    data <- temp[samples,,drop = FALSE]
    covar <- dataset$covar[samples,,drop = FALSE]
    
    for(i in 1:length(genoprobs)) {
        genoprobs[[i]] <- genoprobs[[i]][samples,,]
        K[[i]] <- K[[i]][samples,]
    }
    
    list(data = data, covar = covar, genoprobs = genoprobs, K = K)
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
        return (get(id))
    }
    
    NULL
}


#' Get all "dataset.*" information.
#' 
#' list(dataSets =, ensemblVersion =, numMarkers =)
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
            annotSubset <- ds$annots[which(ds$annots$omit == FALSE),]
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
         ensemblVersion = ensembl.version, 
         numMarkers     = dim(markers)[1])
}


#' Perform the LOD scan.
#' 
#' @param dataset the dataset identifier
#' @param id the identifier
#' @param regressLocal TRUE to regress on local genotype, FALSE to not
#' @param nCores number of cores to use (0=ALL)
#' 
#' @return a data.table with the following columns: id, chr, pos, lod
#' 
GetLODScan <- function(dataset, id, regressLocal=FALSE, nCores=0) {
    # get the dataset
    ds = GetDataSet(dataset)
    if (is.null(ds)) {
        stop(paste0("dataset not found: ", dataset))
    }

    # get the index and data based 
    if (ds$datatype == "mRNA") {
        idx <- which(ds$annots$gene_id == id)
        data <- ds$expr
    } else if (ds$datatype == "protein") {
        idx <- which(ds$annots$protein_id == id)
        data <- ds$expr
    } else if (ds$datatype == "phenotype") {
        idx <- which(ds$annots$data_name == id)
        data <- ds$pheno
    } else {
        stop("invalid datatype")
    }
    
    if (length(idx) == 0) {
        stop(paste0("id not found: ", id))
    }

    # make sure nCores is appropriate  
    numCores = nvlInteger(nCores, 0)
 
    # set covariates
    if (isPheno(ds)) {
        covar_str <- strsplit(ds$annots$use_covar[idx], ":")[[1]]
        covar_str <- paste0("~", paste0(covar_str, collapse = "+"))
        covar <- model.matrix(as.formula(covar_str), data = ds$pheno)[, -1, drop = FALSE]
    } else {
        covar <- ds$covar

        # to regress local genotype, add neareast marker to covariates
        if (toBoolean(regressLocal)) {
            mkr = as.character(markers[ds$annots$nearest.marker.id[idx], 1])
            chr = as.character(markers[ds$annots$nearest.marker.id[idx], 2])
            covar <- cbind(covar, genoprobs[, chr][, -1, mkr])
        }
    }

    # perform the scan using QTL2
    temp <- (scan1(genoprobs = genoprobs,
                   kinship   = K,
                   pheno     = data[, idx, drop = F], 
                   addcovar  = covar, 
                   cores     = numCores,
                   reml      = TRUE))

    # construct a 2 dimensional array of data
    tempDT <- data.table(id  = markers$marker, 
                         chr = markers$chr, 
                         pos = markers$pos, 
                         temp)

    # setting colnames to NULL removes the names in the JSON and return an array
    colnames(tempDT)[4] <- "lod"

    tempDT
}


#' Get the founder coefficients.
#' 
#' @param dataset the dataset identifier
#' @param id the identifier
#' @param chrom the chromosome
#' @param regressLocal TRUE to regress on local genotype, FALSE to not
#' @param blup whether or not to perform BLUP
#' @param center whether or not to center the data
#' @param nCores number of cores to use (0=ALL)
#' 
#' @return a data.table with the following columns: id, chr, pos, and A-H
#' 
GetFoundercoefs <- function(dataset, id, chrom, regressLocal = FALSE, 
                            blup = FALSE, center = TRUE, nCores = 0) {
    # get the dataset
    ds = GetDataSet(dataset)
    if (is.null(ds)) {
        stop(paste0("dataset not found: ", dataset))
    }
    
    # get the index 
    if (ds$datatype == "mRNA") {
        idx <- which(ds$annots$gene_id == id)
        data <- ds$expr
    } else if (ds$datatype == "protein") {
        idx <- which(ds$annots$protein_id == id)
        data <- ds$expr
    } else if (ds$datatype == "phenotype") {
        idx <- which(ds$annots$data_name == id)
        data <- ds$pheno
    } else {
        stop("invalid datatype")
    }
    
    if (length(idx) == 0) {
        stop(paste0("id not found: ", id))
    }
    
    # make sure the chromosome data exists
    if (is.null(K[[chrom]])) {
        stop(paste0("chrom not found: ", chrom))
    }
    
    # make sure ncores is appropriate  
    numCores = nvlInteger(nCores, 0)
    
    # set covariates
    if (isPheno(ds)) {
        covar_str <- strsplit(ds$annots$use_covar[idx], ":")[[1]]
        covar_str <- paste0("~", paste0(covar_str, collapse = "+"))
        covar <- model.matrix(as.formula(covar_str), data = ds$pheno)[, -1, drop = FALSE]
    } else {
        covar <- ds$covar

        # to regress local genotype, add neareast marker to covariates
        if (toBoolean(regressLocal)) {
            mkr = as.character(markers[ds$annots$nearest.marker.id[idx], 1])
            chr = as.character(markers[ds$annots$nearest.marker.id[idx], 2])
            covar <- cbind(covar, genoprobs[, chr][, -1, mkr])
        }
    }
    
    if (toBoolean(blup)) {
        temp <- scan1blup(genoprobs = genoprobs[, chrom],
                          pheno     = data[, idx, drop = F],
                          kinship   = K[[chrom]],
                          addcovar  = covar,
                          cores     = numCores)
    } else {
        temp <- scan1coef(genoprobs = genoprobs[, chrom],
                          pheno     = data[, idx, drop = F],
                          kinship   = K[[chrom]],
                          addcovar  = covar)
    }
    
    if (toBoolean(center)) {
        temp[, LETTERS[1:8]] <- 
            temp[, LETTERS[1:8]] - rowMeans(temp[, LETTERS[1:8]], na.rm = TRUE)
    }
    
    data.table(id = names(map[[chrom]]), 
               chr = chrom, 
               pos = map[[chrom]], 
               temp[,LETTERS[1:8]])
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
    if (is.null(ds)) {
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
    
    if (length(idx) == 0) {
        stop(paste0("id not found: ", id))
    }
    
    # TODO: do we need to create this and pass back everything or maybe a subset?
    # the types of expression data
    t <- list()
    for (f in ds$covar.factors$column.name) {
        stopifnot(!is.null(ds$samples[[f]]))
        if (is.factor(ds$samples[[f]])) {
            t[[f]] <- mixedsort(levels(ds$samples[[f]]))
        } else {
            t[[f]] <- mixedsort(unique(ds$samples[[f]]))
        }
    }
    
    #if (isPheno(ds)) {
    #    # TODO: what is going on here?
    #    output <- list()
    #    phenoFact <- c()
    #    allPheno <- c(grep('^pheno*', apropos('pheno\\.'), value=TRUE))
    #    for (p in allPheno) {
    #        pText <- substr(p, 7, 100)
    #        temp <- cbind(ds$samples, expression=ds$expr[,idx], pheno=pText)
    #        temp[,1] = paste0(temp[,1], "_", pText)
    #        output[[pText]] <- temp
    #        phenoFact = c(phenoFact, pText)
    #    }
    #    
    #    t[["pheno"]] = phenoFact
    #    output <- do.call("rbind", output)
    #    colnames(output)[1] <- ("mouse_id")
    #} else {
    #    # simple to get the expression data and tack on
    #    output <- cbind(ds$samples, expression=ds$expr[,idx])
    #    colnames(output)[1] <- ("mouse_id")
    #}
    
    if (isPheno(ds)) {
        output <- cbind(ds$samples, expression=ds$pheno[, idx])
        colnames(output)[1] <- ("mouse_id")
    } else {
        # simple to get the expression data and tack on
        output <- cbind(ds$samples, expression=ds$expr[, idx])
        colnames(output)[1] <- ("mouse_id")
    }
    
    # elimate the _row column down line for JSON
    rownames(output) <- NULL
    
    list(data = output, dataTypes = t)
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
GetMediate <- function(dataset, id, mid, datasetMediate=NULL) {
    # get the dataset
    ds = GetDataSet(dataset)
    if (is.null(ds)) {
        stop(paste0("dataset not found: ", dataset))
    }

    # get the index 
    idx <- 0
    annot <- NULL
    if (ds$datatype == "mRNA") {
        idx <- which(ds$annots$gene_id == id)
        data <- ds$expr
    } else if (ds$datatype == "protein") {
        idx <- which(ds$annots$protein_id == id)
        data <- ds$expr
    } else if (ds$datatype == "phenotype") {
        idx <- which(ds$annots$data_name == id)
        data <- ds$pheno
    } else {
        stop("invalid datatype")
    }
    
    if (length(idx) == 0) {
        stop(paste0("id not found: ", id))
    }

    # get the dataset we are mediating against
    datasetMediate <- nvl(datasetMediate, dataset)
    dsMediate = GetDataSet(datasetMediate)
    if (is.null(ds)) {
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
    
    if (length(mrkx) == 0) {
        stop(paste0("mid not found: ", mid))
    }
    
    chrTmp = as.character(markers[mrkx, 2])
    annot$middle_point <- dsMediate$annots$middle


    
    # synchronize the samples
    # (list(genoprobs = genoprobs, data = data, K = K, covar = covar))
    #temp <- SynchronizeSamples(dataset   = ds, 
    #                           genoprobs = genoprobs,
    #                           K         = K)

    # perform the mediation
    toReturn <- (mediation.scan(target     = data[,idx, drop=FALSE],
                                mediator   = dsMediate$expr,
                                annotation = annot,
                                covar      = dsMediate$covar,
                                qtl.geno   = genoprobs[[chrTmp]][rownames(dsMediate$expr),,mid],
                                verbose    = FALSE))

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
    if (is.null(ds)) {
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

    if (length(idx) == 0) {
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
        print(paste0("window.range=", window.range))
        myDB <- src_sqlite(db.file, create=FALSE)
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
            window.length <- window.size/1000000.0
            window.range <- pos + c(-1,1)*window.length
            tries <- tries + 1

            if (tries > 10) {
                stop(paste0("Cannot find snps in region: ", window.range))
            }
        }
    }

    print(paste0("window.range=", window.range))
    
    colnames(window.snps)[c(1,3)] = c("snp", "pos")
    window.snps = index_snps(map = map, window.snps)

    # set covariates
    if (isPheno(ds)) {
        covar_str <- strsplit(ds$annots$use_covar[idx], ":")[[1]]
        covar_str <- paste0("~", paste0(covar_str, collapse = "+"))
        covar <- model.matrix(as.formula(covar_str), data = ds$pheno)[, -1, drop = FALSE]
    } else {
        covar <- ds$covar
    }
    
    # convert allele probs to SNP probs
    snppr <- genoprob_to_snpprob(genoprobs, window.snps)

    # finally the scan
    if (isPheno(ds)) {
        outSnps <- scan1(pheno     = ds$pheno[, idx, drop=F], 
                         kinship   = K[[sel.chr]], 
                         genoprobs = snppr, 
                         addcovar  = covar, 
                         cores     = numCores)
    } else {
        outSnps <- scan1(pheno     = ds$expr[, idx, drop=F], 
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
#' 
#' @return a data.frame with the following columns: marker, chr, pos
#' and depending upons dataset$dataType the following columns:
#' mRNA = gene_id, symbol, gene_chrom, middle, lod
#' protein = 
#' phenotype = data_name, short_name, description, lod
#' 
GetLODPeaks <- function(dataset) {
    # get the dataset
    ds = GetDataSet(dataset)
    if (is.null(ds)) {
        stop(paste0("dataset not found: ", dataset))
    }

    if (ds$datatype == "mRNA") {
        temp <- merge(x    = ds$annots[,c("gene_id", "symbol", "chr", "middle")], 
                      y    = ds$lod.peaks[, c("annot.id", "marker.id", "lod")], 
                      by.x = "gene_id", 
                      by.y = "annot.id")

        colnames(temp)[3] <- "gene_chrom"

        return (merge(x    = markers[, c("marker", "chr", "pos")],
                      y    = temp, 
                      by.x = "marker", 
                      by.y = "marker.id"))                      
    } else if (ds$datatype == "protein") {
        temp <- merge(x    = ds$annots[,c("protein_id", "gene_id", "symbol", "chr", "middle")], 
                      y    = ds$lod.peaks[, c("annot.id", "marker.id", "lod")], 
                      by.x = "protein_id", 
                      by.y = "annot.id")

        colnames(temp)[4] <- "gene_chrom"

        return (merge(x    = markers[, c("marker", "chr", "pos")],
                      y    = temp, 
                      by.x = "marker", 
                      by.y = "marker.id"))                      
    } else if (ds$datatype == "phenotype") {
        temp <- merge(x    = ds$annots[,c("data_name", "short_name", "description")], 
                      y    = ds$lod.peaks[, c("annot.id", "marker.id", "lod")], 
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
GetCorrelation <- function(dataset, id, datasetCorrelate=NULL, maxItems=NULL) {
    # get the dataset
    ds <- GetDataSet(dataset)
    if (is.null(ds)) {
        stop(paste0("dataset not found: ", dataset))
    }
    
    data <- NULL
    
    # get the index 
    if (ds$datatype == "mRNA") {
        idx <- which(ds$annots$gene_id == id)
        data <- as.matrix(ds$expr) 
    } else if (ds$datatype == "protein") {
        idx <- which(ds$annots$protein_id == id)
        data <- as.matrix(ds$expr) 
    } else if (ds$datatype == "phenotype") {
        idx <- which(colnames(ds$pheno[,ds$annots$is_pheno == TRUE]) == id)
        data <- as.matrix(ds$pheno[,ds$annots$is_pheno == TRUE]) 
    } else {
        stop("invalid datatype")
    }
    
    if (length(idx) == 0) {
        stop(paste0("id not found: ", id))
    }
    
    # get the dataset to correlate to
    datasetCorrelate <- nvl(datasetCorrelate, dataset)
    
    dsCorrelate = GetDataSet(datasetCorrelate)
    if (is.null(dsCorrelate)) {
        stop(paste0("datasetCorrelate not found: ", dataset))
    }
    
    dataCorrelate <- NULL
    
    if (dsCorrelate$datatype == "mRNA") {
        dataCorrelate <- as.matrix(dsCorrelate$expr)
    } else if (dsCorrelate$datatype == "protein") {
        dataCorrelate <- as.matrix(dsCorrelate$expr)
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
    if (is.null(ds)) {
        stop(paste0("dataset not found: ", dataset))
    }
    
    dsCorrelate <- GetDataSet(datasetCorrelate)
    if (is.null(dsCorrelate)) {
        stop(paste0("dataset not found: ", datasetCorrelate))
    }
    
    # get the index 
    if (ds$datatype == "mRNA") {
        idx <- which(ds$annots$gene_id == id)
        data <- as.matrix(ds$expr) 
    } else if (ds$datatype == "protein") {
        idx <- which(ds$annots$protein_id == id)
        data <- as.matrix(ds$expr) 
    } else if (ds$datatype == "phenotype") {
        idx <- which(colnames(ds$pheno[,ds$annots$is_pheno == TRUE]) == id)
        data <- as.matrix(ds$pheno[,ds$annots$is_pheno == TRUE]) 
    } else {
        stop("invalid datatype")
    }
    
    if (length(idx) == 0) {
        stop(paste0("id not found: ", id))
    }
    
    # get the index of the correlate
    if (dsCorrelate$datatype == "mRNA") {
        idxCorrelate <- which(dsCorrelate$annots$gene_id == idCorrelate)
        dataCorrelate <- as.matrix(dsCorrelate$expr) 
    } else if (dsCorrelate$datatype == "protein") {
        idxCorrelate <- which(dsCorrelate$annots$protein_id == idCorrelate)
        dataCorrelate <- as.matrix(dsCorrelate$expr) 
    } else if (dsCorrelate$datatype == "phenotype") {
        idxCorrelate <- which(colnames(dsCorrelate$pheno[,dsCorrelate$annots$is_pheno == TRUE]) == idCorrelate)
        dataCorrelate <- as.matrix(dsCorrelate$pheno[,dsCorrelate$annots$is_pheno == TRUE]) 
    } else {
        stop("invalid datatype")
    }

    if (length(idxCorrelate) == 0) {
        stop(paste0("id not found: ", idCorrelate))
    }
    
    samples <- intersect(rownames(data), rownames(dataCorrelate))
    data <- data[samples,]
    dataCorrelate <- dataCorrelate[samples,]
 
    toReturn <- list(dataset          = dataset,
                     id               = id,                    
                     datasetCorrelate = datasetCorrelate,
                     idCorrelate      = idCorrelate,
                     data             = data.frame(mouse.id = rownames(data), 
                                                   x        = data[,idx], 
                                                   y        = dataCorrelate[,idxCorrelate]))

    toReturn
}


# #############################################################################
#
# Plumber filters, routes, and utils
#
# #############################################################################

#' @rdname serializers
#' @export
SerializerQTLJSON <- function() {
    function(val, req, res, errorHandler) {
        tryCatch({
            json <- jsonlite::toJSON(val, auto_unbox = TRUE, digits=10)

            res$setHeader("Content-Type", "application/json")
            res$body <- json

            return(res$toResponse())
        }, error = function(e){
            errorHandler(req, res, e)
        })
    }
}

addSerializer("qtlJSON", SerializerQTLJSON)


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


#' The following is used for debugging to show where the requests are coming 
#' from.
#' 
#' @param req the request object
#' 
#* @filter logger
function(req){
    print(paste0("[", Sys.time(), "] [",
                 req$REMOTE_ADDR, "] [",
                 req$REQUEST_METHOD, "] [",
                 req$PATH_INFO, 
                 req$QUERY_STRING, "]"))
    plumber::forward()
}

#' Make sure the server is alive
#'
#' @param req the request object
#' @param res the response object
#'
#' @return 'OK' if alive
#'
#' @serializer unboxedJSON
#* @get /ping
HttpPing <- function(req, res) {
    'OK'
}


#' Get the system information.
#'
#' @param req the request object
#' @param res the response object
#'
#' @return JSON of the system information
#'
#' @serializer unboxedJSON
#* @get /sysinfo
HttpSysInfo <- function(req, res) {
    # start the clock
    ptm <- proc.time()
    
    toReturn <- as.list(Sys.info())
    
    # stop the clock
    elapsed <- proc.time() - ptm
    TrackTime(req, elapsed["elapsed"])
    
    toReturn
}


#' Get dataset information.
#'
#' For now the levels should be either mrna/protein OR pheno.
#'
#' @param req the request object
#' @param res the response object
#'
#' @return JSON object of the options
#'
#' Example Response:
#'    {"dataSets":[{"id":"some_identifier",
#'                  "annotations":[{"dataName":"mouse.id",
#'                                  "shortName":"mouse.id",
#'                                  "desc":"the mouse identifier"}],
#'                  "displayName":"My Data",
#'                  "dataType":"mRNA",
#'                  "ensemblVersion":"90",
#'                  "covarFactors":[{"column.name":"Sex",
#'                                   "display.name":"Sex"},
#'                                  {"column.name":"Generation",
#'                                   "display.name":"Generation"}]
#'                 }]
#'    }
#'
#' @serializer qtlJSON
#' @get /datasets
HttpDatasetInfo <- function(req, res) {
    
    result <- tryCatch(
        {
            ptm <- proc.time()
            
            datasets <- GetDatasetInfo()
            
            elapsed <- proc.time() - ptm
            TrackTime(req, elapsed["elapsed"])
            
            list(result = datasets,
                 time   = elapsed["elapsed"])
        },
        error = function(cond) {
            res$status <- 400
            list(error=jsonlite::unbox(cond$message))
        }
    )

    result
}


#' Perform the LOD scan
#'
#' @param req the request object
#' @param res the response object
#' @param dataset the dataset identifier
#' @param id the identifier
#' @param regressLocal TRUE to regress on local genotype, FALSE to not
#' @param nCores number of cores to use (0=ALL)
#' @param expand TRUE to expand the JSON, FALSE to condense
#'
#' @return JSON data
#'
#' Example of expand=FALSE result:
#'     {"result":[["1_4530778","1",40055,1.688],
#'                ["1_4533435","1",4.5334,1.6709],
#'                ...
#'                ["1_4536092","1",4.5361,1.6539]
#'               ],
#'      "time":5.3
#'     }
#'
#' Example of expand=TRUE result:
#'     {"result":[{"id":"1_4530778","chr":"1","pos":4.5308,"lod":1.688},
#'                {"id":"1_4533435","chr":"1","pos":4.5334,"lod":1.6709},
#'                ...
#'                {"id":"1_4536092","chr":"1","pos":4.5361,"lod":1.6539}
#'               ],
#'      "time":5.3
#'     }
#'     
#' @serializer qtlJSON
#* @get /lodscan
#* @post /lodscan
HttpLODScan <- function(req, res, dataset, id, regressLocal=FALSE, nCores=0, 
                        expand=FALSE) {
    
    result <- tryCatch(
        {
            ptm <- proc.time()
            
            lod <- GetLODScan(dataset      = dataset, 
                              id           = id,
                              regressLocal = regressLocal,
                              nCores       = nCores)
            
            if (!(toBoolean(expand))) {
                # by setting column names to NULL, the result will be a
                # 2 dimensional array
                colnames(lod) <- NULL
            }

            elapsed <- proc.time() - ptm
            TrackTime(req, elapsed["elapsed"])

            list(result = lod,
                 time   = elapsed["elapsed"])
        },
        error = function(cond) {
            res$status <- 400
            list(error=jsonlite::unbox(cond$message))
        }
    )

    result
}


#' Get founder coefficients
#'
#' @param req the request object
#' @param res the response object
#' @param dataset the dataset identifier
#' @param id an identifier
#' @param chrom the chromosome
#' @param regressLocal TRUE to regress local genotype
#' @param blup TRUE to perform Best Linear Unbiased Predictors 
#' @param center TRUE to center around 0
#' @param nCores number of cores to use (0=ALL)
#'
#' @return JSON data
#'
#' Example of expand=FALSE result:
#'     {"result":[["2_4292516","2",4.2925,
#'                 0.012,-0.0267,0.0345,0.1994,
#'                 0.1305,-0.4769,0.1217,0.0055},
#'                {"2_4337139","2",4.3371,
#'                 0.0114,-0.0267,0.0343,0.2,
#'                 0.1307,-0.4769,0.1218,0.0055},
#'                ...
#'                {"2_4369263","2",4.3693,
#'                 0.012,-0.0259,0.0407,0.2006,
#'                 0.1323,-0.4783,0.1134,0.0052}
#'                ],
#'      "time":3.3
#'     }
#'     
#' Example of expand=TRUE result:
#'     {"result":[{"id":"2_4292516","chr":"2","pos":4.2925,
#'                 "A":0.012,"B":-0.0267,"C":0.0345,"D":0.1994,
#'                 "E":0.1305,"F":-0.4769,"G":0.1217,"H":0.0055},
#'                {"id":"2_4337139","chr":"2","pos":4.3371,
#'                 "A":0.0114,"B":-0.0267,"C":0.0343,"D":0.2,
#'                 "E":0.1307,"F":-0.4769,"G":0.1218,"H":0.0055},
#'                ...
#'                {"id":"2_4369263","chr":"2","pos":4.3693,
#'                 "A":0.012,"B":-0.0259,"C":0.0407,"D":0.2006,
#'                 "E":0.1323,"F":-0.4783,"G":0.1134,"H":0.0052}
#'                ],
#'      "time":3.3
#'     }
#'     
#' @serializer qtlJSON
#* @get /foundercoefs
#* @post /foundercoefs
HttpFoundercoefs <- function(req, res, dataset, id, chrom, regressLocal=FALSE, 
                             blup=FALSE, center=TRUE, nCores=0, expand=FALSE) {
    ptm <- proc.time()
    
    result <- tryCatch(
        {
            ptm <- proc.time()
            
            effect <- GetFoundercoefs(dataset      = dataset, 
                                      id           = id,
                                      chrom        = chrom, 
                                      regressLocal = regressLocal, 
                                      blup         = blup, 
                                      center       = center, 
                                      nCores       = nCores)
            
            if (!(toBoolean(expand))) {
                # by setting column names to NULL, the result will be a
                # 2 dimensional array
                colnames(effect) <- NULL
            }
            
            elapsed <- proc.time() - ptm
            TrackTime(req, elapsed["elapsed"])
            
            list(result = effect,
                 time   = elapsed["elapsed"])
        },
        error = function(cond) {
            res$status <- 400
            list(error=jsonlite::unbox(cond$message))
        }
    )

    # return the data
    result
}


#' Get the expression
#'
#' @param req the request object
#' @param res the response object
#' @param dataset the dataset identifier
#' @param id an identifier
#'
#' @return JSON data
#'
#' JSON data is dependent upon covariates.
#'
#' Example:
#'  {"data":[{"mouse_id":"F01","Sex":"F","Generation":"G4","Litter":2,
#'            "Diet":"hf","Coat.Color":"agouti","GenerationLitter":"G4_2",
#'            "expression":0.2243,"_row":"F01"},
#'           ...
#'           {"mouse_id":"F02","Sex":"F","Generation":"G4","Litter":2,
#'            "Diet":"hf","Coat.Color":"black","GenerationLitter":"G4_2",
#'            "expression":0.3498,"_row":"F02"}
#'          ]
#'
#' @serializer qtlJSON
#* @get /expression
#* @post /expression
HttpExpression <- function(req, res, dataset, id) {
    # start the clock
    ptm <- proc.time()
    
    result <- tryCatch(
        {
            ptm <- proc.time()
            
            expression <- GetExpression(dataset = dataset, 
                                        id      = id)
            
            elapsed <- proc.time() - ptm
            TrackTime(req, elapsed["elapsed"])
            
            list(result = expression,
                 time   = elapsed["elapsed"])
        },
        error = function(cond) {
            res$status <- 400
            list(error=jsonlite::unbox(cond$message))
        }
    )
    
    result
}


#' Perform mediation analysis.
#'
#' @param req the request object
#' @param res the response object
#' @param dataset the dataset identifier
#' @param id an identifier
#' @param mid a marker identifier
#' @param datasetMediate the dataset identifier to mediate against
#' @param expand TRUE to expand the JSON, FALSE to condense
#'
#' @return JSON data
#'
#' Example of expand=FALSE result:
#'     {"result":[["ENSMUSG00000000001","Gnai3","3",108.1267,1.3625],
#'                ["ENSMUSG00000000028","Cdc45","16",18.7962,1.605],
#'                 ...
#'                [ENSMUSG00000000037","Scml2","X",161.1877,1.2535]
#'               ],
#'      "time":10.3
#'     }
#'
#' Example of expand=TRUE result:
#'     {"result":[{"gene_id":"ENSMUSG00000000001","symbol":"Gnai3",
#'                 "chr":"3","pos":108.1267,"LOD":1.3625},
#'                {"gene_id":"ENSMUSG00000000028","symbol":"Cdc45",
#'                 "chr":"16","pos":18.7962,"LOD":1.605},
#'                 ...
#'                {"gene_id":"ENSMUSG00000000037","symbol":"Scml2",
#'                 "chr":"X","pos":161.1877,"LOD":1.2535}
#'               ],
#'      "time":10.3
#'     }
#'     
#' @serializer qtlJSON
#* @get /mediate
#* @post /mediate
HttpMediate <- function(req, res, dataset, id, mid, datasetMediate=NULL, expand=FALSE) {
    result <- tryCatch(
        {
            ptm <- proc.time()
            
            mediation <- GetMediate(dataset        = dataset, 
                                    id             = id,
                                    mid            = mid,
                                    datasetMediate = datasetMediate)
            
            if (!(toBoolean(expand))) {
                # by setting column names to NULL, the result will be a
                # 2 dimensional array
                colnames(mediation) <- NULL
            }
            
            elapsed <- proc.time() - ptm
            TrackTime(req, elapsed["elapsed"])
            
            list(result = mediation,
                 time   = elapsed["elapsed"])
        },
        error = function(cond) {
            res$status <- 400
            list(error=jsonlite::unbox(cond$message))
        }
    )
    
    #res$setHeader("Content-Type", "application/json")
    #res$body <- toJSON(result, dataframe = c("values"))
    #res
    
    result
}


#' Perform on SNP Association mapping.
#'
#' @param req the request object
#' @param res the response object
#' @param dataset the dataset identifier
#' @param id unique dentifier
#' @param chrom the chromsome
#' @param location the location in base pairs
#' @param windowSize how many base pairs (before and after) to perform scan
#' @param nCores number of cores to use (0=as many as there is)
#' @param expand TRUE to expand the JSON, FALSE to condense
#'
#' @return JSON data
#'
#' Example of expand=FALSE result:
#'     {"result":[["rs245710663","1",14.5001,
#'                 "G|C",64,:"",
#'                 "intergenic_variant",1,147,
#'                 false,:0.0015],
#'                ["rs263649601","chr":"1","pos":14.5002,
#'                  "C|T",64,"",
#'                  "intergenic_variant",1,147,
#'                  false,0.0015],
#'                 ...
#'                ["rs32279922","chr":"1","pos":14.5004,
#'                  "C|T",:153,:"",
#'                  "intergenic_variant",3,:147,
#'                  false,3.5395e-06]
#'               ],
#'      "time":11.1
#'     }
#'
#' Example of expand=TRUE result:
#'     {"result":[{"snp":"rs245710663","chr":"1","pos":14.5001,
#'                 "alleles":"G|C","sdp":64,"ensembl_gene":"",
#'                 "csq":"intergenic_variant","index":1,"interval":147,
#'                 "on_map":false,"lod":0.0015},
#'                 {"snp":"rs263649601","chr":"1","pos":14.5002,
#'                 "alleles":"C|T","sdp":64,"ensembl_gene":"",
#'                 "csq":"intergenic_variant","index":1,"interval":147,
#'                 "on_map":false,"lod":0.0015},
#'                 ...
#'                 {"snp":"rs32279922","chr":"1","pos":14.5004,
#'                  "alleles":"C|T","sdp":153,"ensembl_gene":"",
#'                  "csq":"intergenic_variant","index":3,"interval":147,
#'                  "on_map":false,"lod":3.5395e-06}
#'               ],
#'      "time":11.1
#'     }
#'     
#' @serializer qtlJSON
#* @get /snpassoc
#* @post /snpassoc
HTTPSnpAssocMapping <- function(req, res, dataset, id, chrom, location, 
                                windowSize=500000, nCores=0, expand=FALSE) {
    # start the clock
    ptm <- proc.time()
    
    result <- tryCatch(
        {
            ptm <- proc.time()
            
            snpAssoc <- GetSnpAssocMapping(dataset    = dataset, 
                                           id         = id,
                                           chrom      = chrom,
                                           location   = location,
                                           windowSize = windowSize,
                                           nCores     = nCores)            
            if (!(toBoolean(expand))) {
                # by setting column names to NULL, the result will be a
                # 2 dimensional array
                colnames(snpAssoc) <- NULL
            }
            
            elapsed <- proc.time() - ptm
            TrackTime(req, elapsed["elapsed"])
            
            list(result = snpAssoc,
                   time = elapsed["elapsed"])
        },
        error = function(cond) {
            res$status <- 400
            list(error=jsonlite::unbox(cond$message))
        }
    )

    result
}


#' Get LOD Peaks for a dataset
#'
#' @param req the request object
#' @param res the response object
#' @param dataset the dataset identifier
#' @param expand TRUE to expand the JSON, FALSE to condense
#'
#' @return JSON data
#' 
#' Example of expand=FALSE result and dataset$dataType=mRNA:
#'     {"result":[["1_100007442","1",100.0074,
#'                 "ENSMUSG00000028028","Alpk1","3",127.7255,6.1147],
#'                ["1_100078094","1",100.0781,
#'                 "ENSMUSG00000024727","Trpm6","19",18.8212,6.048],
#'                ...
#'                ["1_100148722","1",100.1487,
#'                 "ENSMUSG00000045725","Prr15","6",54.3286,6.5079],
#'               ],
#'      "time":4.9
#'     }
#'
#' Example of expand=TRUE result and dataset$dataType=mRNA:
#'     {"result":[{"marker":"1_100007442","chr":"1","pos":100.0074,
#'                 "gene_id":"ENSMUSG00000028028","symbol":"Alpk1",
#'                 "gene_chrom":"3","middle":127.7255,"lod":6.1147},
#'                {"marker":"1_100078094","chr":"1","pos":100.0781,
#'                 "gene_id":"ENSMUSG00000024727","symbol":"Trpm6",
#'                 "gene_chrom":"19","middle":18.8212,"lod":6.048},
#'                ...
#'                {"marker":"1_100148722","chr":"1","pos":100.1487,
#'                 "gene_id":"ENSMUSG00000045725","symbol":"Prr15",
#'                 "gene_chrom":"6","middle":54.3286,"lod":6.5079}
#'               ],
#'      "time":4.9
#'     }
#'
#' @serializer qtlJSON
#* @get /lodpeaks
#* @post /lodpeaks
HttpLODPeaks <- function(req, res, dataset, expand=FALSE) {
    result <- tryCatch(
        {
            ptm <- proc.time()
            
            lodPeaks <- GetLODPeaks(dataset)

            if (!(toBoolean(expand))) {
                # by setting column names to NULL, the result will be a
                # 2 dimensional array
                colnames(lodPeaks) <- NULL
            }
            
            elapsed <- proc.time() - ptm
            TrackTime(req, elapsed["elapsed"])
            
            list(result = lodPeaks,
                 time = elapsed["elapsed"])
        },
        error = function(cond) {
            res$status <- 400
            list(error=jsonlite::unbox(cond$message))
        }
    )

    result
}


#' Perform correlation on a phenotype or gene.
#'
#' @param req the request object
#' @param res the response object
#' @param dataset the dataset identifier
#' @param id unique dentifier
#' @param dataset the dataset identifier to correlate to
#' @param maxItems maximum number of items to return
#'
#' @return JSON data
#'
#' Example of expand=FALSE result:
#'     
#' @serializer qtlJSON
#* @get /correlation
#* @post /correlation
HTTPCorrelation <- function(req, res, dataset, id, datasetCorrelate=NULL, maxItems=NULL) {
    # start the clock
    ptm <- proc.time()
    
    result <- tryCatch(
        {
            ptm <- proc.time()
            
            correlation <- GetCorrelation(dataset          = dataset, 
                                          id               = id,
                                          datasetCorrelate = datasetCorrelate,
                                          maxItems         = maxItems)

            elapsed <- proc.time() - ptm
            TrackTime(req, elapsed["elapsed"])
            
            list(result = correlation,
                 time = elapsed["elapsed"])
        },
        error = function(cond) {
            res$status <- 400
            list(error=jsonlite::unbox(cond$message))
        }
    )

    result
}


#' Perform correlation on a phenotype or gene.
#'
#' @param req the request object
#' @param res the response object
#' @param dataset the dataset identifier
#' @param id unique dentifier
#' @param datasetCorrelate the dataset identifier to correlate to
#' @param idCorrelate id from correlate dataset
#'
#' @return JSON data
#'
#' Example of expand=FALSE result:
#'     
#' @serializer qtlJSON
#* @get /correlationPlotData
#* @post /correlationPlotData
HTTPCorrelationPlotData <- function(req, res, dataset, id, datasetCorrelate, idCorrelate) {
    # start the clock
    ptm <- proc.time()
    
    result <- tryCatch(
        {
            ptm <- proc.time()
            
            correlation <- GetCorrelationPlotData(dataset          = dataset, 
                                                  id               = id,
                                                  datasetCorrelate = datasetCorrelate,
                                                  idCorrelate      = idCorrelate)

            elapsed <- proc.time() - ptm
            TrackTime(req, elapsed["elapsed"])
            
            list(result = correlation,
                 time = elapsed["elapsed"])
        },
        error = function(cond) {
            res$status <- 400
            list(error=jsonlite::unbox(cond$message))
        }
    )

    result
}
