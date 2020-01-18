# #############################################################################
#
# Load the libraries will we need
#
# #############################################################################

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

if (exists('debug_mode')) {
    debug_mode <- get('debug_mode')
} else {
    debug_mode <- FALSE
}

if (debug_mode) {
    print('DEBUG MODE: Make sure the data is loaded and db_file is defined')
} else {
    print('Finding the data file to load...')
    data_files <- list.files('/app/qtlapi/data', 
                             '\\.RData$', 
                             ignore.case = TRUE, 
                             full.names = TRUE)
    
    if (length(data_files) == 0) {
        stop('There needs to be an .RData file in /app/qtlapi/data')
    } else if (length(data_files) > 1) {
        stop('There needs to be only 1 .RData file in /app/qtlapi/data')
    }
    
    print(sprintf('Loading the data file: %s', data_files))
    
    load(data_files, .GlobalEnv)
    
    db_file <- list.files('/app/qtlapi/data', 
                          '\\.sqlite$', 
                          ignore.case =TRUE, 
                          full.names = TRUE)
    
    if (length(db_file) == 0) {
        stop('There needs to be an .sqlite file in /app/qtlapi/data')
    } else if (length(db_file) > 1) {
        stop('There needs to be only 1 .sqlite file in /app/qtlapi/data')
    }
    
    print(sprintf('Using SNP db file: %s', db_file))
}


# #############################################################################
#
# Utility functions
#
# #############################################################################


#' Opposite of %in%  
`%not in%` <- function (x, table) match(x, table, nomatch = 0L) == 0L

#' Convert to boolean/logical
#' 
#' Acceptable TRUE boolean values are TRUE, 1, "T", "TRUE", "YES", "Y", "1"
#'
#' @param value The value to check.
#' 
#' @return TRUE if it is a boolean value and it's value is one of the
#'   acceptable matching inputs.
#' 
#' @noRd
to_boolean <- function(value) {
    if (is.logical(value)) {
        return(value)
    } else if (is.character(value)) {
        return(toupper(value) %in% c('T', 'TRUE', 'YES', 'Y', '1'))
    } else if (is.numeric(value)) {
        return(value == 1)
    }
    
    FALSE
}


#' Check value for validity and return it or a default
#'
#' @param value The value to check.
#' @param default The default value to return.
#' 
#' @return Either value if it is valid or the default value.
#' 
#' @noRd
nvl <- function(value, default) {
    if (gtools::invalid(value)) {
        return(default)
    }
    
    value
}


#' Covert value to numeric if possible
#'
#' @param value The value to check.
#' @param default The default value to return.
#' 
#' @return Either value if it is valid or the default value.
#' 
#' @noRd
nvl_int <- function(value, default) {
    tryCatch({
        n <- as.numeric(value)
        if ((n %% 1) == 0) {
            return(n)
        }
    },
    error = function(cond) {
    },
    warning = function(cond) {
    },
    finally={
    })
    
    default
}


# #############################################################################
#
# Data methods
#
# #############################################################################


#' Get the dataset by id (a string)
#' 
#' @param id A string, either 'dataset.name' or just 'name'.
#' 
#' @return The dataset element.
#' 
get_dataset <- function(id) {
    if (exists(id)) {
        return(get(id))
    } else {
        expanded <- sprintf('dataset.%s', id) 
        if (exists(expanded)) {
            return(get(expanded))
        }
    }
    
    stop(sprintf('dataset "%s" not found', id))
}


#' Check dataset to see if the datatype value is "phenotype"
#'
#' @param dataset Either list containing all values for a dataset or character 
#'   string which is the name of the dataset.
#' 
#' @return TRUE if the datatype is phenotype, FALSE otherwise.
#' 
is_phenotype <- function(dataset) {
    if (typeof(dataset) == 'character') {
        dataset <- get_dataset(dataset)
    }
    
    if (startsWith(tolower(dataset$datatype), 'pheno')) {
        return(TRUE)
    }

    FALSE
}


#' Get the dataset by id (a string)
#' 
#' @param id The dataset id as a string.
#' @param data A string containing which data element from the dataset's data
#'     element.
#' 
#' @return The data element.
#' 
get_data <- function(id, data = NULL) {
    ds <- get_dataset(id)
    
    ret <- NULL
    
    if (is.matrix(ds$data) && (!is.null(data))) {
        stop(sprintf('specified data "%s" not found in "%s"', data, id))
    }
    
    # rz, norm, raw, log, transformed
    if (!gtools::invalid(data)) {
        if (is.matrix(ds$data)) {
            stop(sprintf('only 1 data element in "%s"', id))
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
    
    if (gtools::invalid(ret)) {
        stop(sprintf('unable to find data in dataset "%s"', id))
    }
    
    ret
}

#' Get all "dataset.*" information
#' 
#' This will return a named list of all datasets and the ensmebl version.
#'
#' @return A named list of all the dataset objects, along with the 
#'   ensembl.version.
#' 
get_dataset_info <- function() {
    datasets <- grep('^dataset*', apropos('dataset\\.'), value = TRUE)
    ret <- c()
    
    for (d in datasets) {
        ds <- get(d)
        
        annotations <- list()
        
        if (tolower(ds$datatype) == 'mrna') {
            annotations <- list(ids = ds$annot.mrna$gene.id)
        } else if(tolower(ds$datatype) == 'protein') {
            annotations <- 
                list(ids = tibble(protein.id = ds$annot.protein$protein.id,
                                  gene.id    = ds$annot.protein$gene.id))
        } else if(is_phenotype(ds)) {
            annotations <- 
                as_tibble(ds$annot.phenotype[which(ds$annot.phenotype$omit == FALSE), ])
        }
        
        temp <- list(id              = d, 
                     annotations     = annotations, 
                     display.name    = nvl(ds$display.name, d), 
                     datatype        = ds$datatype, 
                     covar.info      = ds$covar.info)
        
        ret <- c(ret, list(temp))
    }
    
    list(datasets        = ret, 
         ensembl.version = nvl(ensembl.version, NA))
}

#' Get all "dataset.*" statistics.
#' 
#' This will return a named list of all datasets and some statistics.
#'
#' @return A named list of all the dataset objects.
#' 
get_dataset_stats <- function() {
    datasets <- grep('^dataset*', apropos('dataset\\.'), value = TRUE)
    ret <- c()
    
    for (d in datasets) {
        ds <- get(d)
        
        num.annotations <- NA
        
        if (tolower(ds$datatype) == 'mrna') {
            num.annotations <- NROW(ds$annot.mrna)
        } else if(tolower(ds$datatype) == 'protein') {
            num.annotations <- NROW(ds$annot.protein)
        } else if(is_phenotype(ds)) {
            num.annotations <- NROW(ds$annot.phenotype)
        }
        
        temp <- list(id              = d, 
                     display.name    = nvl(ds$display.name, d), 
                     datatype        = ds$datatype,
                     num.annotations = num.annotations,
                     num.samples     = NROW(ds$annot.samples))
        
        ret <- c(ret, list(temp))
    }
    
    ret
}


#' Loop through the datasets to see if the id is used.
#' 
#' @param id The identifier.
#' 
#' @return A named list of all the dataset objects.
#' 
has_annotation <- function(id) {
    datasets <- grep('^dataset*', apropos('dataset\\.'), value = TRUE)
    ret <- c()
    
    for (d in datasets) {
        ds <- get(d)
        data <- get_data(d)
        
        found <- FALSE

        if (id %in% colnames(data)) {
            found <- TRUE
        }

        temp <- list(id      = id,
                     dataset = d,
                     found   = found) 

        ret <- c(ret, list(temp))
    }
    
    ret
}


#' Perform a LOD scan
#' 
#' @param dataset The dataset identifier.
#' @param id The identifier.
#' @param covar The interactive covariate.
#' @param cores number of cores to use (0=ALL).
#' 
#' @return A tibble with the following columns: id, chr, pos, lod.
#' 
get_lod_scan <- function(dataset, id, intcovar = NULL, cores = 0) {
    # get the dataset and the data
    ds <- get_dataset(dataset)
    data <- get_data(dataset)
    
    # check if id exists 
    idx <- which(colnames(data) == id)
    
    if (gtools::invalid(idx)) {
        stop(sprintf('id "%s" not found', id))
    }
    
    # grab the data for the id
    id_data <- data[, id, drop = FALSE]
    
    # make sure num_cores is appropriate  
    num_cores = nvl_int(cores, 0)
    
    if (is_phenotype(ds)) {
        # get the annot.pheno row to get use.covar variable from the annotations
        pheno <- ds$annot.pheno %>% filter(data.name == id)
        
        if (gtools::invalid(pheno)) {
            stop(sprintf('id (data.name) "%s" not found in annot.phenotype', id))
        }
        
        # convert samples to data.frame because QTL2 relies heavily
        # on rownames and colnames, rownames currently are or will
        # soon be deprecated in tibbles
        samples <- as.data.frame(ds$annot.samples)
        
        # set the rownames so scan1 will work, finding a match to some
        # variation of mOuSE.Id
        rownames(samples) <- 
            (samples %>% select(matches('^mouse\\.id$')))[[1]]

        # create a string (model formula) from the use.covar column
        formula_str <- paste0('~', gsub(':', '+', pheno$use.covar))
        
        # [, -1, drop = FALSE] will drop the (Intercept) column
        covar <- model.matrix(as.formula(formula_str), data = samples)
        covar <- covar[, -1, drop = FALSE]
    } else {
        covar <- NULL
        
        if (!gtools::invalid(ds$covar.matrix)) {
            covar <- ds$covar.matrix    
        }
    } 
    
    # set the interactive.covariate, to be used in scan1
    # as scan1(intcovar=interactive_covariate)
    interactive_covariate <- NULL
    
    if (!gtools::invalid(intcovar)) {
        if (intcovar %not in% ds$covar.info$sample.column) {
            stop(sprintf('intcovar "%s" not found in %s$covar.info', 
                         intcovar, 
                         dataset))
        }
        
        # grabbing all the columns from covar (covar.matrix) that
        # match, i.e., "batch" will match "batch2", "BATCH3", etc
        interactive_covariate <-
            covar[, which(grepl(intcovar, colnames(covar), ignore.case = T))]
    }
    
    # perform the scan using QTL2, 
    # - addcovar should always be ALL covars
    # - intcovar should be just the interactive covariate column
    temp <- scan1(genoprobs = genoprobs,
                  kinship   = K,
                  pheno     = id_data, 
                  addcovar  = covar, 
                  intcovar  = interactive_covariate,
                  cores     = num_cores,
                  reml      = TRUE)
    
    # construct a 2 dimensional array of data with id, chr, pos, lod as columns
    # we perform a left join here to make sure that the number of elements match
    left_join(as_tibble(temp, rownames = 'marker.id'), 
              markers, 
              by = 'marker.id') %>% 
        select(id = marker.id, chr, pos, lod = id)
}


#' Perform the LOD scan for each "value" of the interactive covariate.
#' 
#' @param dataset The dataset identifier.
#' @param id The identifier.
#' @param intcovar The interactive covariate, sample.column in covar.info.
#' @param chrom The chromosome.
#' @param cores Number of cores to use (0=ALL).
#' 
#' @return A named list with each name a "sample value" and the element is a 
#'   tibble with the following columns: id, chr, pos, lod.
#' 
get_lod_scan_by_sample <- function(dataset, id, intcovar, chrom, cores = 0) {
    # get the dataset and the data
    ds <- get_dataset(dataset)
    data <- get_data(dataset)
    
    # check if id exists 
    idx <- which(colnames(data) == id)
    
    if (gtools::invalid(idx)) {
        stop(sprintf('id "%s" not found', id))
    }
    
    # make sure the chromosome data exists
    if (gtools::invalid(K[[chrom]])) {
        stop(sprintf('chrom "%s" not found', chrom))
    }
    
    # extract the markers just for the chromosome
    markers_chrom <- markers %>% filter(chr == chrom)
    
    # make sure nCores is appropriate  
    num_cores <- nvl_int(cores, 0)
    
    if (intcovar %not in% ds$covar.info$sample.column) {
        stop(sprintf('covar "%s" not found in %s$covar.info', intcovar, dataset))
    }

    # get all the unique values for the interactive.covar and sort them 
    if (is.factor(ds$annot.samples[[intcovar]])) {
        covar_unique <- gtools::mixedsort(levels(ds$annot.samples[[intcovar]]))
    } else {
        covar_unique <- gtools::mixedsort(unique(ds$annot.samples[[intcovar]]))
    }

    # convert samples to data.frame because QTL2 relies heavily
    # on rownames and colnames, rownames currently are or will
    # soon be deprecated in tibbles
    samples <- as.data.frame(ds$annot.samples)
    
    # set the rownames so scan1 will work, finding a match to some
    # variation of mOuSE.Id
    rownames(samples) <- (samples %>% select(matches('^mouse\\.id$')))[[1]]

    # ret will be a named list of tibbles with LOD scores
    # each name is a unique sample value
    ret <- list()
    
    # loop through the unique values for the intcovar
    for (u in covar_unique) {
        # samples.names will contain ONLY the samples that match x
        # take all samples
        # filter rows by value, i.e. sex = "F"
        # select just need the mouse.id column
        sample_names <- 
            ds$annot.samples %>%                        
            filter(UQ(as.name(intcovar)) == u) %>%  
            select(matches('^mouse\\.id$'))                  
        
        sample_names <- c(sample_names[[1]])

        # set covariates
        if (is_phenotype(ds)) {
            # get the annot.pheno row to get use.covar variable
            pheno <- ds$annot.pheno %>% filter(data.name == id)
            
            if (gtools::invalid(pheno)) {
                stop(sprintf('id "%s" not found in annot.phenotype', id))
            }
            
            # create a string (model formula) from the use.covar column
            formula_str <- paste0('~', gsub(':', '+', pheno$use.covar))
            
            # construct the covar matrix
            covar <- model.matrix.lm(as.formula(formula_str), 
                                     data = samples, 
                                     na.action = na.pass)
            
            # [, -1, drop = FALSE] will drop the (Intercept) column
            covar <- covar[sample_names, -1, drop = FALSE]
        } else {
            # subset the covars, don't need -1 since no (Intercept)
            covar <- ds$covar.matrix[sample_names, , drop = FALSE]
        }
        
        # exclude covar columns that contain it's name 
        # since this is a matrix we cannot use %>% select
        covar <- 
            covar[, -which(grepl(intcovar, colnames(covar), ignore.case = T))]
        
        temp <- scan1(genoprobs = genoprobs[sample_names, chrom],
                      kinship   = K[[chrom]][sample_names, sample_names],
                      pheno     = data[sample_names, id, drop = FALSE],
                      addcovar  = covar, 
                      cores     = num_cores,
                      reml      = TRUE)
        
        # construct a 2 dimensional array of data with id, chr, pos, lod
        ret[[toString(u)]] <- tibble(id  = markers_chrom$marker.id,
                                     chr = markers_chrom$chr,
                                     pos = markers_chrom$pos,
                                     lod = temp[, 1])
    }
    
    ret
}


#' Get the founder coefficients
#' 
#' @param dataset The dataset identifier.
#' @param id The identifier.
#' @param chrom The chromosome.
#' @param intcovar The interactive covariate, sample.column in covar.info.
#' @param blup whether or not to perform BLUP
#' @param center whether or not to center the data
#' @param cores Number of cores to use (0=ALL).
#' 
#' @return A nameds list with each element being a tibble with the following 
#'         columns: id, chr, pos, and A-H
#' 
get_founder_coefficients <- function(dataset, id, intcovar, chrom, 
                                     blup = FALSE, center = TRUE, cores = 0) {
    # get the dataset and data
    ds <- get_dataset(dataset)
    data <- get_data(dataset)
    
    # check if id exists 
    idx <- which(colnames(data) == id)
    
    if (gtools::invalid(idx)) {
        stop(sprintf('id "%s" not found', id))
    }
    
    # make sure the chromosome data exists
    if (gtools::invalid(K[[chrom]])) {
        stop(sprintf('chrom "%s" not found', chrom))
    }
    
    # make sure ncores is appropriate  
    num_cores <- nvl_int(cores, 0)
    
    # set covariates
    if (is_phenotype(ds)) {
        # get the annot.pheno row to get use.covar variable from the annotations
        pheno <- ds$annot.pheno %>% filter(data.name == id)
        
        if (gtools::invalid(pheno)) {
            stop(sprintf('id "%s" not found in annot.phenotype', id))
        }
        
        # convert samples to data.frame because QTL2 relies heavily
        # on rownames and colnames, rownames currently are or will
        # soon be deprecated in tibbles
        samples <- as.data.frame(ds$annot.samples)
        
        # set the rownames so scan1 will work, finding a match to some
        # variation of mOuSE.Id
        rownames(samples) <- (samples %>% select(matches('^mouse\\.id$')))[[1]]

        # create a string (model formula) from the use.covar column
        formula_str <- paste0('~', gsub(':', '+', pheno$use.covar))
        
        # construct the covar matrix
        covar <- model.matrix.lm(as.formula(formula_str), 
                                 data = samples,
                                 na.action = na.pass)
        
        # [, -1, drop = FALSE] will drop the (Intercept) column
        covar <- covar[, -1, drop = FALSE]
    } else {
        covar <- NULL
        
        if (!gtools::invalid(ds$covar.matrix)) {
            covar <- ds$covar.matrix    
        }
    }
    
    ret <- list()
    
    if (gtools::invalid(intcovar)) {
        if (to_boolean(blup)) {
            temp <- scan1blup(genoprobs = genoprobs[, chrom],
                              pheno     = data[, idx, drop = FALSE],
                              kinship   = K[[chrom]],
                              addcovar  = covar,
                              cores     = num_cores)
        } else {
            temp <- scan1coef(genoprobs = genoprobs[, chrom],
                              pheno     = data[, idx, drop = FALSE],
                              kinship   = K[[chrom]],
                              addcovar  = covar)
        }
        
        if (to_boolean(center)) {
            temp[, LETTERS[1:8]] <- 
                temp[, LETTERS[1:8]] - rowMeans(temp[, LETTERS[1:8]], 
                                                na.rm = TRUE)
        }
        
        ret[['additive']] <- as_tibble(data.frame(id = names(map[[chrom]]), 
                                                  chr = chrom, 
                                                  pos = map[[chrom]], 
                                                  temp[,LETTERS[1:8]]))
    } else {
        if (intcovar %not in% ds$covar.info$sample.column) {
            stop(sprintf('covar "%s" not found in %s$covar.info', intcovar, dataset))
        }
        
        # get all the unique values for the intcovar and sort them 
        if (is.factor(ds$annot.samples[[intcovar]])) {
            covar.unique <- 
                gtools::mixedsort(levels(ds$annot.samples[[intcovar]]))
        } else {
            covar.unique <- 
                gtools::mixedsort(unique(ds$annot.samples[[intcovar]]))
        }
        
        # convert samples to data.frame because QTL2 relies heavily
        # on rownames and colnames, rownames currently are or will
        # soon be deprecated in tibbles
        samples <- as.data.frame(ds$annot.samples)
        
        # set the rownames so scan1 will work, finding a match to some
        # variation of mOuSE.Id
        rownames(samples) <- 
            (samples %>% select(matches('^mouse\\.id$')))[[1]]
        
        # loop through the unique values for the interactive.covar
        for (u in covar.unique) {
            # sample_names will contain ONLY the samples that match x
            # take all samples
            # filter rows by value, i.e. sex = "F
            # select just need the mouse.id column
            sample_names <- 
                samples %>%                                   
                filter(UQ(as.name(intcovar)) == u) %>%  
                select(matches('^mouse\\.id$'))                  
            
            sample_names <- c(sample_names[[1]])
            
            # exclude covar columns that contain it's name
            covar_subset <- 
                covar[sample_names, -which(grepl(intcovar, 
                                                 colnames(covar), 
                                                 ignore.case = T))]
            
            if (to_boolean(blup)) {
                temp <- scan1blup(genoprobs = genoprobs[sample_names, chrom],
                                  pheno     = data[sample_names, idx, drop = FALSE],
                                  kinship   = K[[chrom]][sample_names, sample_names],
                                  addcovar  = covar_subset,
                                  cores     = num_cores)
            } else {
                temp <- scan1coef(genoprobs = genoprobs[sample_names, chrom],
                                  pheno     = data[sample_names, idx, drop = FALSE],
                                  kinship   = K[[chrom]][sample_names, sample_names],
                                  addcovar  = covar_subset)
            }
            
            if (to_boolean(center)) {
                temp[, LETTERS[1:8]] <- 
                    temp[, LETTERS[1:8]] - rowMeans(temp[, LETTERS[1:8]], 
                                                    na.rm = TRUE)
            }
            
            ret[[toString(u)]] <- as_tibble(data.frame(id = names(map[[chrom]]), 
                                                       chr = chrom, 
                                                       pos = map[[chrom]], 
                                                       temp[, LETTERS[1:8]]))
        }
    }
    
    ret
}


#' Get the expression data.
#' 
#' @param dataset The dataset identifier.
#' @param id The identifier.
#' 
#' @return A list with elements of tibble and list of the datatypes.
#' 
get_expression <- function(dataset, id) {
    # get the dataset and data
    ds <- get_dataset(dataset)
    data <- get_data(dataset)
    
    # check if id exists 
    idx <- which(colnames(data) == id)
    
    if (gtools::invalid(idx)) {
        stop(sprintf('id "%s" not found', id))
    }
    
    if (gtools::invalid(ds$covar.info)) {
        datatypes <- NULL
    } else {
        datatypes <- list()
        for (i in ds$covar.info$sample.column) {
            stopifnot(!is.null(ds$annot.samples[[i]]))
            if (is.factor(ds$annot.samples[[i]])) {
                datatypes[[i]] <- mixedsort(levels(ds$annot.samples[[i]]))
            } else {
                datatypes[[i]] <- mixedsort(unique(ds$annot.samples[[i]]))
            }
        }
    }
    
    # make sure 'mouse.id' is lowercase when passed back and only pass back the
    # columns we need
    samples <- ds$annot.samples %>% 
               rename(mouse.id = matches('^mouse\\.id$')) %>%
               select(c('mouse.id', names(datatypes)))
    
    # bind the data
    output <- cbind(samples, expression = data[, idx])
    
    list(data = output, 
         datatypes = datatypes)
}    


#' Perform mediation
#' 
#' @param dataset The dataset identifier.
#' @param id The identifier.
#' @param marker_id The id of the marker.
#' @param intcovar The interactive covariate, sample.column in covar.info.
#' @param dataset_mediate The dataset to mediate against.
#' 
#' @return A data.frame with the following columns depending on datatype: 
#'         mRNA = gene_id, symbol, chr, pos, LOD
#'         protein = protein_id, gene_id, symbol, chr, pos, LOD
#'         phenotype = NONE
#'         
get_mediation <- function(dataset, id, marker_id, intcovar, 
                          dataset_mediate = NULL) {
    # get the dataset and data
    ds <- get_dataset(dataset)
    data <- get_data(dataset)
    
    # get the dataset we are mediating against
    ds_mediate <- get_dataset(nvl(dataset_mediate, dataset))
    data_mediate <- get_data(nvl(dataset_mediate, dataset))

    # check if id exists 
    idx <- which(colnames(data) == id)
    
    if (gtools::invalid(idx)) {
        stop(sprintf('id "%s" not found', id))
    }
    
    # get the marker index 
    mrkx <- which(markers$marker.id == marker_id)
    
    if (gtools::invalid(mrkx)) {
        stop(sprintf('marker.id "%s" not found', marker_id))
    }

    # get the annotations
    if (tolower(ds_mediate$datatype) == 'mrna') {
        # grab the annotations, create middle_point, and select what is needed
        annot <-
            inner_join(enframe(colnames(data_mediate), name=NULL),
                       ds_mediate$annot.mrna,
                       by = c('value' = 'gene.id')) %>%
            mutate(middle_point = (start + end) / 2) %>%
            select(gene.id = value, symbol, chr, middle_point)
        
    } else if (tolower(ds_mediate$datatype) == 'protein') {
        # grab the annotations, create middle_point, and select what is needed
        annot <-
            inner_join(enframe(colnames(data_mediate), name=NULL),
                       ds_mediate$annot.protein,
                       by = c('value' = 'protein.id')) %>%
            mutate(middle_point = (start + end) / 2) %>%
            select(protein.id = value, gene.id, symbol, chr, middle_point)

    } else if (is_phenotype(ds_mediate)) {
        stop('invalid datatype to mediate against')
    } else {
        stop('invalid datatype')
    }
    
    chrom = as.character(markers[mrkx, 'chr'])
    
    covar <- NULL
    
    # perform the mediation
    if (!gtools::invalid(ds_mediate$covar.matrix)) {
        covar <- ds_mediate$covar.matrix
    }
        
    mediation.scan(target     = data[, idx, drop = FALSE],
                   mediator   = data_mediate,
                   annotation = annot,
                   covar      = covar,
                   qtl.geno   = genoprobs[[chrom]][rownames(data_mediate), , marker_id],
                   verbose    = FALSE)
}


#' Get the SNP association mapping
#' 
#' @param dataset The dataset identifier.
#' @param id The identifier.
#' @param chrom The chromosome.
#' @param intcovar The interactive covariate, sample.column in covar.info.
#' @param chrom The chromosome.
#' @param location Location on chromosome in base pairs.
#' @param window_size The size of the window to scan before and after location.
#' @param cores Number of cores to use (0=ALL).
#' 
#' @return A data.frame with the following columns: snp, chr, pos, alleles, 
#'   sdp, ensembl_gene, csq, index, interval, on_map, lod.
#' 
get_snp_assoc_mapping <- function(dataset, id, intcovar, chrom, location,
                                  window_size = 500000, cores = 0) {
    # get the dataset and data
    ds <- get_dataset(dataset)
    data <- get_data(dataset)
    
    # check if id exists 
    idx <- which(colnames(data) == id)
    
    if (gtools::invalid(idx)) {
        stop(sprintf('id "%s" not found', id))
    }

    # make sure location is valid    
    if (nvl_int(location, -1) == -1) {
        stop(sprintf('location "%s" is invalid', location))
    }
    
    position <- nvl_int(location, -1) / 1000000.0
    
    if (nvl_int(window_size, -1) == -1) {
        stop(sprintf('window_size is invalid: %d', window_size))
    }
    
    window_size <- nvl_int(window_size, -1)
    window_length <- window_size / 1000000.0
    window_range_start <- max(position - window_length, 0)
    window_range_end <- position + window_length
    
    # make sure nCores is appropriate  
    num_cores = nvl_int(cores, 0)
    
    have_snps <- FALSE
    tries <- 1
    
    # extract SNPs from the database, we allow up to 10 tries to
    # find SNPs in a window
    while (!have_snps) {
        db_snps <- src_sqlite(db_file, create = FALSE)
        window_snps <- 
            tbl(db_snps, 'snps') %>%
            filter(chr     == chrom, 
                   pos_Mbp >= window_range_start,
                   pos_Mbp <= window_range_end) %>%
            arrange(pos_Mbp) %>%
            collect(n = Inf)
        
        if (NROW(window_snps) > 0) {
            have_snps <- TRUE
        } else {
            window_size <- window_size + nvl_int(window_size, -1)
            window_length <- window_size / 1000000.0
            window_range_start <- max(position - window_length, 0)
            window_range_end <- position + window_length
            tries <- tries + 1
            
            if (tries > 10) {
                stop(sprintf('Cannot find snps in region: %s:%f-%f', 
                             chrom, window.range_start, window.range_end))
            }
        }
    }
    
    colnames(window_snps)[c(1,3)] = c('snp', 'pos')
    window_snps = index_snps(map = map, window_snps)
    
    # set covariates
    if (is_phenotype(ds)) {
        # get the annot.pheno row to get use.covar variable from the annotations
        pheno <- ds$annot.pheno %>% filter(data.name == id)
        
        if (gtools::invalid(pheno)) {
            stop(sprintf('id "%s" not found in annot.phenotype', id))
        }
        
        # convert samples to data.frame because QTL2 relies heavily
        # on rownames and colnames, rownames currently are or will
        # soon be deprecated in tibbles
        samples <- as.data.frame(ds$annot.samples)
        
        # set the rownames so scan1 will work, finding a match to some
        # variation of mOuSE.Id
        rownames(samples) <- (samples %>% select(matches('^mouse\\.id$')))[[1]]
        
        # create a string (model formula) from the use.covar column
        formula_str <- paste0('~', gsub(':', '+', pheno$use.covar))
        
        # construct the covar matrix
        covar <- model.matrix.lm(as.formula(formula_str), 
                                 data = samples, 
                                 na.action = na.pass)
        
        # [, -1, drop = FALSE] will drop the (Intercept) column
        covar <- covar[, -1, drop = FALSE]
    } else {
        covar <- NULL
        
        if (!gtools::invalid(ds$covar.matrix)) {
            covar <- ds$covar.matrix    
        }
    }
    
    # convert allele probs to SNP probs
    snp_prob <- genoprob_to_snpprob(genoprobs, window_snps)
    
    # perform the scan using QTL2, 
    # - addcovar should always be ALL covars
    # - intcovar should be just the interactive covariate column
    out_snps <- scan1(pheno     = data[, idx, drop = FALSE], 
                      kinship   = K[[chrom]], 
                      genoprobs = snp_prob, 
                      addcovar  = covar,
                      cores     = num_cores)
    
    map_tmp <- qtl2:::snpinfo_to_map(window_snps)
    tmp <- qtl2:::expand_snp_results(out_snps, map_tmp, window_snps)
    
    ret <- window_snps
    ret$lod <- tmp$lod[, 1]
    
    # set the interactive_covariates, to be used in scan1
    # as scan1(intcovar=interactive.covariate)
    interactive_covariate <- NULL

    if (!gtools::invalid(intcovar)) {
        if (intcovar %not in% ds$covar.info$sample.column) {
            stop(sprintf('covar "%s" not found in %s$covar.info', intcovar, dataset))
        }
        
        # grabbing all the columns from covar (covar.matrix) that
        # match, i.e., "batch" will match "batch2", "BATCH3", etc
        interactive_covariate <- 
            covar[, which(grepl(intcovar, colnames(covar), ignore.case = T))]

        # perform the scan using QTL2, 
        # - addcovar should always be ALL covars
        # - intcovar should be just the interactive covariate column
        out_snps <- scan1(pheno     = data[, idx, drop = FALSE], 
                          kinship   = K[[chrom]], 
                          genoprobs = snp_prob, 
                          addcovar  = covar,
                          intcovar  = interactive_covariate,
                          cores     = num_cores)
    
        map_tmp <- qtl2:::snpinfo_to_map(window_snps)
        tmp <- qtl2:::expand_snp_results(out_snps, map_tmp, window_snps)
        
        ret$lod_intcovar <- tmp$lod[, 1]
    }
    
    ret
}


#' Get the LOD peaks
#' 
#' @param dataset The dataset identifier.
#' @param intcovar The interactive covariate, sample.column in covar.info.
#' 
#' @return a data.frame with the following columns: marker, chr, pos
#' and depending upons dataset$dataType the following columns:
#' mRNA = gene_id, symbol, gene_chrom, middle, lod
#' protein = protein_id, gene_id, symbol, gene_chrom, middle, lod
#' phenotype = data_name, short_name, description, lod
#' 
get_lod_peaks <- function(dataset, intcovar = NULL) {
    # get the dataset
    ds <- get_dataset(dataset)
    
    peaks <- NULL
    
    if (gtools::invalid(intcovar)) {
        peaks <- ds$lod.peaks$additive
    } else {
        # find the covar and get the name of the lod peaks
        if (intcovar %in% ds$covar.info$sample.column) {
            n <- ds$covar.info[ds$covar.info$sample.column == intcovar, ]
            peaks <- ds$lod.peaks[[n$lod.peaks]]  
        } else {
            stop(sprintf('covar "%s" not found in %s$covar.info', intcovar, dataset))
        }
    }
    
    if (gtools::invalid(peaks)) {
        stop(sprintf('no peaks found for covar "%s" in %s', intcovar, dataset))
    }
    
    if (tolower(ds$datatype) == 'mrna') {
        temp <- merge(x    = ds$annot.mrna[, c('gene.id', 'symbol', 'chr', 'start', 'end')], 
                      y    = peaks[, c('gene.id', 'marker.id', 'lod')], 
                      by.x = 'gene.id', 
                      by.y = 'gene.id')
        
        temp$gene.pos <- (temp$start + temp$end) / 2.0
        colnames(temp)[3] <- 'gene.chr'
        
        return (merge(x    = markers[, c('marker.id', 'chr', 'pos')],
                      y    = temp[, c('gene.id', 'symbol', 'gene.chr', 'gene.pos', 'lod', 'marker.id')], 
                      by.x = 'marker.id', 
                      by.y = 'marker.id'))
    } else if (tolower(ds$datatype) == 'protein') {
        temp <- merge(x    = ds$annot.protein[,c('protein.id', 'gene.id', 'symbol', 'chr', 'start', 'end')], 
                      y    = peaks[, c('protein.id', 'marker.id', 'lod')], 
                      by.x = 'protein.id', 
                      by.y = 'protein.id')
        
        temp$gene.pos <- (temp$start + temp$end) / 2.0
        colnames(temp)[4] <- 'gene.chr'
        
        return (merge(x    = markers[, c('marker.id', 'chr', 'pos')],
                      y    = temp[, c('protein.id', 'gene.id', 'symbol','gene.chr', 'gene.pos', 'lod', 'marker.id')],
                      by.x = 'marker.id', 
                      by.y = 'marker.id'))                      
    } else if(is_phenotype(ds)) {
        temp <- merge(x    = ds$annot.phenotype[,c('data.name', 'short.name', 'description')], 
                      y    = peaks[, c('data.name', 'marker.id', 'lod')], 
                      by.x = 'data.name', 
                      by.y = 'data.name')
        
        return (merge(x    = markers[, c('marker.id', 'chr', 'pos')],
                      y    = temp, 
                      by.x = 'marker.id', 
                      by.y = 'marker.id'))
    } else {
        stop('invalid datatype')
    }
}    


#' Get the LOD peaks
#' 
#' @param dataset The dataset identifier.
#' 
#' @return a data.frame with the following columns: marker, chr, pos
#' and depending upons dataset$dataType the following columns:
#' mRNA = gene_id, symbol, gene_chrom, middle, lod
#' protein = protein_id, gene_id, symbol, gene_chrom, middle, lod
#' phenotype = data_name, short_name, description, lod
#' 
get_lod_peaks_all <- function(dataset) {
    # get the dataset
    ds = get_dataset(dataset)
    
    # get the additive LOD peaks
    peaks <- list(additive = get_lod_peaks(dataset))
    
    # get the rest
    for (i in seq(nrow(ds$covar.info))) {
        inf <- ds$covar.info[i, ]
        
        if (inf$interactive) {
            peaks[[inf$sample.column]] <- 
                get_lod_peaks(dataset, inf$sample.column)
        }
    } 
    
    peaks
}    


#' Calculate the resudual matrix
#' 
#' @param variable_matrix The data  matrix for first set.
#' @param adjust_matrix The data matrix for the second set.
#' @param variables_interest List of variables of interest.
#' @param variables_compare List of variables to compare.
#' @param use_qr TRUE to use QR decomposition (FASTER).
#' @param impute TRUE to impute NA data.
#' 
#' @return residual matrix
#' 
calculate_residual_matrix <- function(variable_matrix, 
                                      adjust_matrix, 
                                      variables_interest, 
                                      variables_compare,
                                      use_qr = TRUE,
                                      impute = TRUE) {
    # make sure we have the same samples
    samples <- intersect(rownames(variable_matrix), rownames(adjust_matrix))

    # combine the data    
    data <- 
        as.data.frame(cbind(variable_matrix[samples,], adjust_matrix[samples,]))
    
    if(use_qr) {
        # Fast, but can't handle NAs in y??
        formula_str <- 
            paste('~ + 1 +', paste(variables_compare, collapse = ' + '))
        
        X_0 <- model.matrix.lm(as.formula(formula_str), 
                               data = data,
                               na.action = na.pass)
        
        X_0 <- X_0[samples, ]
        
        y_data <- data[samples, variables_interest, drop = FALSE]
        colnames(y_data) <- variables_interest
        
        if (impute) {
            if (any(is.na(X_0))) {
                X_0 <- missMDA::imputeFAMD(X = X_0)$completeObs
            }
            
            if (any(is.na(y_data))) {
                y_data <- missMDA::imputeFAMD(X = y_data)$completeObs
            }
        } 
        
        qr_0 <- qr(X_0)
        
        residual_matrix <- sapply(seq(variables_interest), function(i) {
            d <- y_data[, variables_interest[i]]
            return(qr.resid(qr_0, d))
        }, simplify = TRUE)
        
        rownames(residual_matrix) <- samples
    }
    else{
        ## Way too slow, use QR trick
        residual_matrix <- sapply(seq(variables_interest), function(i) {
            formula_str <- 
                paste(variables_interest[i], '~', paste(variables_compare, collapse = ' + '))
            fit <- lm(formula(formula_str), data = data)
            return(fit$residuals[samples])
        }, simplify = TRUE)
    }
    
    colnames(residual_matrix) <- variables_interest
    residual_matrix
}



#' Get the correlation.
#' 
#' @param dataset The dataset identifier.
#' @param id The identifier.
#' @param dataset_correlate The dataset to correlate against.
#' @param intcovar The interactive covariate, sample.column in covar.info.
#' @param use_qr TRUE to use QR decomposition (FASTER).
#' 
#' 
#' @return A tibble with the correlation and annotations.
#' 
get_correlation <- function(dataset, id, 
                            dataset_correlate = NULL, intcovar = NULL,
                            use_qr = TRUE) {
    # get the dataset and data
    ds <- get_dataset(dataset)
    data <- get_data(dataset)
    
    # get the dataset to correlate to
    ds_correlate = get_dataset(nvl(dataset_correlate, dataset))
    data_correlate <- get_data(nvl(dataset_correlate, dataset))
    
    # check if id exists 
    idx <- which(colnames(data) == id)
    
    if (gtools::invalid(idx)) {
        stop(sprintf('id "%s" not found', id))
    }

    # make sure we have the same samples
    samples <- intersect(rownames(data), rownames(data_correlate))
    data <- data[samples, ]
    data_correlate <- data_correlate[samples, ]    

    if (gtools::invalid(intcovar)) {
        pcor <- cor(data[, idx], data_correlate, use = 'pair')
    } else {
        interactive_covariate <- 
            colnames(ds$covar.matrix)[grepl(intcovar, 
                                            colnames(ds$covar.matrix), 
                                            ignore.case = T)]
        
        id_residual_matrix <- 
            calculate_residual_matrix(variable_matrix    = data,
                                      adjust_matrix      = ds$covar.matrix, 
                                      variables_interest = c(id),
                                      variables_compare  = interactive_covariate,
                                      use_qr             = use_qr)
            
        residual_matrix <- 
            calculate_residual_matrix(variable_matrix    = data_correlate,
                                      adjust_matrix      = ds$covar.matrix, 
                                      variables_interest = colnames(data_correlate),
                                      variables_compare  = interactive_covariate,
                                      use_qr             = use_qr)
        
        resid_samples <- intersect(rownames(id_residual_matrix), 
                                   rownames(residual_matrix))
        
        pcor <- cor(id_residual_matrix[resid_samples, ], 
                    residual_matrix[resid_samples, ], 
                    use = 'pair')
    }

    # reorder
    pcor <- pcor[1, order(abs(pcor), decreasing = TRUE)]
    
    if (tolower(ds_correlate$datatype) == 'mrna') {
        return (tibble(cor    = pcor, 
                       id     = names(pcor), 
                       symbol = ds_correlate$annot.mrna$symbol[match(names(pcor), ds_correlate$annot.mrna$gene.id)],
                       chr    = ds_correlate$annot.mrna$chr[match(names(pcor), ds_correlate$annot.mrna$gene.id)],
                       start  = ds_correlate$annot.mrna$start[match(names(pcor), ds_correlate$annot.mrna$gene.id)],
                       end    = ds_correlate$annot.mrna$end[match(names(pcor), ds_correlate$annot.mrna$gene.id)]))
    } else if (tolower(ds_correlate$datatype) == 'protein') {
        return (tibble(cor     = pcor,
                       id      = names(pcor),
                       gene_id = ds_correlate$annot.protein$gene.id[match(names(pcor), ds_correlate$annot.protein$protein.id)],
                       symbol  = ds_correlate$annot.protein$symbol[match(names(pcor), ds_correlate$annot.protein$protein.id)],
                       chr     = ds_correlate$annot.protein$chr[match(names(pcor), ds_correlate$annot.protein$protein.id)],
                       start   = ds_correlate$annot.protein$start[match(names(pcor), ds_correlate$annot.protein$protein.id)],
                       end     = ds_correlate$annot.protein$end[match(names(pcor), ds_correlate$annot.protein$protein.id)]))
    } else if (is_phenotype(ds_correlate)) {
        return (tibble(cor    = pcor, 
                       id     = names(pcor)))
    }
}


#' Get the correlation.
#' 
#' @param dataset the dataset identifier
#' @param id the identifier
#' @param dataset.correlate the dataset to correlate to
#' @param id.correlate the identifier from the correlate dataset
#' @param intcovar The interactive covariate, sample.column in covar.info.
#' 
#' @return A named list with the data to plot
#' 
get_correlation_plot_data <- function(dataset, id, 
                                      dataset_correlate, id_correlate,
                                      intcovar = NULL) {
    # get the dataset and data
    ds <- get_dataset(dataset)
    data <- get_data(dataset)
    
    # get the dataset to correlate to
    dataset_correlate <- nvl(dataset_correlate, dataset)
    ds_correlate = get_dataset(dataset_correlate)
    data_correlate <- get_data(dataset_correlate)

    # make sure we have the same samples
    samples <- intersect(rownames(data), rownames(data_correlate))
    data <- data[samples, ]
    data_correlate <- data_correlate[samples, ]    
    
    # check if id exists 
    idx <- which(colnames(data) == id)
    
    if (gtools::invalid(idx)) {
        stop(sprintf('id "%s" not found', id))
    }
    
    # get the index
    idx_correlate <- which(colnames(data_correlate) == id_correlate)
    
    if (gtools::invalid(idx_correlate)) {
        stop(sprintf('id "%s" not found: ', id_correlate))
    }
    
    if (!gtools::invalid(intcovar)) {
        interactive_covariate <- 
            colnames(ds$covar.matrix)[grepl(intcovar, 
                                            colnames(ds$covar.matrix), 
                                            ignore.case = T)]
            
        data <- 
            calculate_residual_matrix(variable_matrix    = data,
                                      adjust_matrix      = ds$covar.matrix, 
                                      variables_interest = c(id),
                                      variables_compare  = interactive_covariate,
                                      use_qr             = TRUE)
        
        data_correlate <- 
            calculate_residual_matrix(variable_matrix    = data_correlate,
                                      adjust_matrix      = ds$covar.matrix, 
                                      variables_interest = colnames(data_correlate),
                                      variables_compare  = interactive_covariate,
                                      use_qr             = TRUE)
        
        x <- data[, 1]
        y <- data_correlate[, idx_correlate]
    } else {
        x <- data[, idx]
        y <- data_correlate[, idx_correlate]
    }
        
    id_column <- match('mouse.id', tolower(colnames(ds$annot.samples)))
    samples <- intersect(rownames(data), rownames(data_correlate))
    samples_idx <- which(ds$annot.samples[[id_column]] %in% samples)
    
    # get the covar factors data
    sample_info <- list()
    dt <- list()
    for (s in ds$covar.info$sample.column) {
        stopifnot(!is.null(ds$annot.samples[[s]]))
        sample_info[[toString(s)]] <- ds$annot.samples[samples_idx, ][[s]]
        
        if (is.factor(ds$annot.samples[[s]])) {
            dt[[toString(s)]] <- 
                gtools::mixedsort(levels(ds$annot.samples[[s]]))
        } else {
            dt[[toString(s)]] <- 
                gtools::mixedsort(unique(ds$annot.samples[[s]]))
        }
    }
    
    ret_data <- 
        as_tibble(data.frame(mouse.id         = rownames(data), 
                             x                = x, 
                             y                = y,
                             sample_info,
                             stringsAsFactors = FALSE))

    list(dataset           = dataset,
         id                = id,
         dataset.correlate = dataset_correlate,
         id.correlate      = id_correlate,
         datatypes         = dt,
         data              = ret_data)
}


