library(dplyr)
library(tibble)

#' Opposite of %in%  
`%not in%` <- function (x, table) match(x, table, nomatch = 0L) == 0L


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
    
    if (startsWith(tolower(dataset$datatype), 'pheno')) {
        return(TRUE)
    }
    
    FALSE
}




#' Use the existing covar.matrix or create it.
#'
#' @param dataset The dataset id as a string.
#' @param id The phenotype identifier.
#'
#' @return The covar element.
#'
get_covar_matrix <- function(dataset, id = NULL) {
    ds <- get_dataset(dataset)
    
    if (exists('covar.matrix', ds)) {
        covar <- ds$covar.matrix
    } else {
        # we can generate covar.matrix if it doesn't exist
        if (is_phenotype(ds)) {
            # get the annot.pheno row to get use.covar variable from the 
            # annotations
            pheno <- ds$annot.pheno %>% filter(data.name == id)

            if (gtools::invalid(pheno)) {
                stop(sprintf("Cannot find phenotype '%s' in '%s'", id, dataset))
            }
            
            # create a string (model formula) from the use.covar column
            formula_str <- paste0("~", gsub(":", "+", pheno$use.covar))
        } else {
            formula_str <- paste0(ds$covar.info$sample.column, collapse="+")
            formula_str <- paste0("~", formula_str)
        }

        # convert samples to data.frame because QTL2 relies heavily
        # on rownames and colnames, rownames currently are or will
        # soon be deprecated in tibbles
        samples <- as.data.frame(ds$annot.samples)

        # set the rownames so scan1 will work, finding a match to some
        # variation of mOuSE.Id
        rownames(samples) <-
            (samples %>% select(matches("^mouse\\.id$")))[[1]]

        # [, -1, drop = FALSE] will drop the (Intercept) column
        covar <- model.matrix.lm(
            as.formula(formula_str), 
            data = samples,
            na.action = na.pass
        )
        
        covar <- covar[, -1, drop = FALSE]        
    }
    
    covar
}


check_dataset <- function(d) {
    # get the dataset
    ds <- NA
    
    if (exists(d)) {
        ds <- get(d)
    }
    
    if (gtools::invalid(ds)) {
        message(paste0("dataset not found '", d, '"'))
    } 
    
    cat("\nChecking dataset '", d, "'\n", sep = '')
    
    if (!is.list(ds)) {
        message(paste0("dataset should be a list, but found '", class(d), "'"))
    }
    
    if ('datatype' %not in% names(ds)) {
        message("dataset must contain 'datatype'")
    }
    
    datatype = ds[['datatype']]
    is.mrna <- FALSE
    is.protein <- FALSE
    is.phenotype <- FALSE
    
    if (tolower(datatype) == 'mrna') {
        is.mrna <- TRUE
    } else if (tolower(datatype) == 'protein') {
        is.protein <- TRUE
    } else if (is_phenotype(ds)) {
        is.phenotype <- TRUE
    } else {
        message(paste0("datatype is '", datatype, "', but should be mRNA, protein, or phenotype"))
    }
    
    ###########################################################################
    #
    # annotations
    #
    ###########################################################################
    cat("Checking annotations\n")
    
    if (is.mrna) {
        if ('annot.mrna' %not in% names(ds)) {
            message(paste0("annot.mrna not found in '", d, "'"))
        }
        
        if (!is_tibble(ds$annot.mrna)) {
            message(paste0("annot.mrna should be a tibble, but found '", class(ds$annot.mrna), "'"))
        }
        
        if (any(duplicated(ds$annot.mrna$gene.id))) {
            message("There are duplicated gene.id annotations in annot.mrna")
        }
    } else if (is.protein) {
        if ('annot.protein' %not in% names(ds)) {
            message(paste0("annot.protein not found in '", d, "'"))
        }
        
        if (!is_tibble(ds$annot.protein)) {
            message(paste0("annot.protein should be a tibble, but found ", class(ds$annot.protein)))
        }
        
        if (any(duplicated(ds$annot.protein$protein.id))) {
            message("There are duplicated protein.id annotations in annot.protein")
        }
    } else if (is.phenotype) {
        if ('annot.phenotype' %not in% names(ds)) {
            message(paste0("annot.phenotype not found in '", d, "'"))
        }
        
        if (!is_tibble(ds$annot.pheno)) {
            message(paste0("annot.phenotype should be a tibble, but found ", class(ds$annot.phenotype)))
        }
        
        if (any(duplicated(ds$annot.phenotype$data.name))) {
            message("There are duplicated data.name annotations in annot.phenotype")
        }
    }
    
    if (is.phenotype) {
        if ('data.name' %not in% names(ds$annot.phenotype)) {
            message('data.name not found in annot.phenotype')
        } 
        if ('short.name' %not in% names(ds$annot.phenotype)) {
            message('short.name not found in annot.phenotype')
        } 
        if ('description' %not in% names(ds$annot.phenotype)) {
            message('description not found in annot.phenotype')
        } 
        if ('units' %not in% names(ds$annot.phenotype)) {
            message('units not found in annot.phenotype')
        } 
        if ('is.id' %not in% names(ds$annot.phenotype)) {
            message('is.id not found in annot.phenotype')
        } 
        if ('category' %not in% names(ds$annot.phenotype)) {
            message('category not found in annot.phenotype')
        }
        if ('R.category' %not in% names(ds$annot.phenotype)) {
            message('R.category not found in annot.phenotype')
        }
        if ('is.numeric' %not in% names(ds$annot.phenotype)) {
            message('is.numeric not found in annot.phenotype')
        }
        if ('is.date' %not in% names(ds$annot.phenotype)) {
            message('is.date not found in annot.phenotype')
        }
        if ('is.factor' %not in% names(ds$annot.phenotype)) {
            message('is.factor not found in annot.phenotype')
        }
        if ('factor.levels' %not in% names(ds$annot.phenotype)) {
            message('factor.levels not found in annot.phenotype')
        }
        if ('is.covar' %not in% names(ds$annot.phenotype)) {
            message('is.covar not found in annot.phenotype')
        }
        if ('is.pheno' %not in% names(ds$annot.phenotype)) {
            message('is.pheno not found in annot.phenotype')
        } 
        if ('is.derived' %not in% names(ds$annot.phenotype)) {
            message('is.derived not found in annot.phenotype')
        }
        if ('omit' %not in% names(ds$annot.phenotype)) {
            message('omit not found in annot.phenotype')
        } 
        if ('use.covar' %not in% names(ds$annot.phenotype)) {
            message('use.covar not found in annot.phenotype')
        }
        
        # one and only 1 ID
        a.id <- ds$annot.phenotype[which(ds$annot.phenotype$is.id == TRUE),]$data.name
        
        if(length(a.id) != 1) {
            message('annot.phenotype$is.id should have 1 and only 1 row set to TRUE')
        }
    } else {
        annot <- NULL
        annot.name <- NULL
        
        if (is.protein) {
            if ('protein.id' %not in% names(ds$annot.protein)) {
                message('protein.id not found in annot.protein')
            }
            annot <- ds$annot.protein
            annot.name <- 'annot.protein'
        } else {
            annot <- ds$annot.mrna
            annot.name <- 'annot.mrna'
        }
        
        if ('gene.id' %not in% names(annot)) {
            message(paste('gene.id not found in ', annot.name))
        } else if ('symbol' %not in% names(annot)) {
            message(paste('symbol not found in ', annot.name))
        } else if ('chr' %not in% names(annot)) {
            message(paste('chr not found in ', annot.name))
        } else if ('start' %not in% names(annot)) {
            message(paste('start not found in ', annot.name))
        } else if ('end' %not in% names(annot)) {
            message(paste('end not found in ', annot.name))
        } else if ('strand' %not in% names(annot)) {
            message(paste('strand not found in ', annot.name))
        } else if ('middle' %not in% names(annot)) {
            message(paste('strand not found in ', annot.name))
        } else if ('nearest.marker.id' %not in% names(annot)) {
            message(paste('nearest.marker.id not found in ', annot.name))
        }
        
        if (any(annot$start > 10000.0)) {
            message(paste0(annot.name, '$start should be in Mbp not bp', sep=''))
        } else if (any(annot$end > 10000.0)) {
            message(paste0(annot.name, '$end should be in Mbp not bp', sep=''))
        } else if (any(annot$middle > 10000.0)) {
            message(paste0(annot.name, '$middle should be in Mbp not bp', sep=''))
            message()
        } 
    }
    
    ###########################################################################
    #
    # samples
    #
    ###########################################################################
    cat("Checking annot.samples\n")
    
    if ('annot.samples' %not in% names(ds)) {
        message('annot.samples not found')
    }
    
    if (!is_tibble(ds$annot.sample)) {
        message(paste0("annot.samples should be a tibble, but found ", class(ds$annot.sample)))
    }
    
    id.column <- match('mouse.id', tolower(colnames(ds$annot.samples)))
    if (gtools::invalid(id.column)) {
        message('mouse.id not found in annot.samples')
    }
    
    if (any(duplicated(ds$annot.sample$mouse.id))) {
        message("There are duplicated mouse.id annotations in annot.samples")
    }
    
    ###########################################################################
    #
    # covar.info
    #
    ###########################################################################
    cat("Checking covar.info\n")
    
    if ('covar.info' %not in% names(ds)) {
        message('covar.info not found')
    }
    
    if (!is_tibble(ds$covar.info)) {
        message(paste0("covar.info should be a tibble, but found '", class(ds$covar.info), "'"))
    }
    
    if ('sample.column' %not in% names(ds$covar.info)) {
        message('sample.column not found in covar.info')
    } else if ('display.name' %not in% names(ds$covar.info)) {
        message('display.name not found in covar.info')
    } else if ('interactive' %not in% names(ds$covar.info)) {
        message('interactive not found in covar.info')
    } else if ('primary' %not in% names(ds$covar.info)) {
        message('primary not found in covar.info')
    } else if ('lod.peaks' %not in% names(ds$covar.info)) {
        message('lod.peaks not found in covar.info')
    }
    
    for(i in 1:nrow(ds$covar.info)) {
        row <- ds$covar.info[i, ]
        
        if (row$sample.column %not in% colnames(ds$annot.samples)) {
            message(paste0("covar.info$sample.column ('", row$sample.column, "') is not a column name in annot.samples"))
        }
        
        if (gtools::invalid(row$display.name)) {
            message("covar.info$display.name needs to have a value")
        }
        
        if (row$interactive) {
            if (gtools::invalid(row$lod.peaks)) {
                message("covar.info$interactive is TRUE, but covar.info$lod.peaks is NA")
            } else {
                # check for existence of lod.peaks
                if (gtools::invalid(ds$lod.peaks[[row$lod.peaks]])) {
                    message(paste0("covar.info$interactive is TRUE, but covar.info$lod.peaks ('", row$lod.peaks, "') is not in lod.peaks"))
                }
            }
        } else {
            if (!gtools::invalid(row$lod.peaks)) {
                message(paste0("covar.info$interactive is FALSE, but covar.info$lod.peaks ('", row$lod.peaks, "') is set"))
            }
        }
    }
    
    if (!any(ds$covar.info$primary)) {
        message("covar.info$primary needs at least 1 value set to TRUE")
    }
    
    
    ###########################################################################
    #
    # display.name
    #
    ###########################################################################
    cat("Checking display.name\n")
    
    if ('display.name' %not in% names(ds)) {
        message('display.name not found, will use dataset name')
    }
    
    ###########################################################################
    #
    # lod.peaks
    #
    ###########################################################################
    cat("Checking to see if we can generate covar.matrix\n")
    
    if (is_phenotype(ds)) {
        phenos <- ds$annot.phenotype %>% filter(is.pheno == TRUE)
        phenos <- phenos$data.name
        
        for (x in phenos) {
            tryCatch(
                {
                    temp <- temp <- get_covar_matrix(d, x)
                },
                error = function(cond) {
                    message('Unable to generate covar.matrix, check covar.info')
                    message(cond$message)
                    message('Check phenotype ', x)
                },
                warning = function(cond) {
                },
                finally = {
                }
            )
        }
        
    } else {
        tryCatch(
            {
                temp <- get_covar_matrix(d)
            },
            error = function(cond) {
                message('Unable to generate covar.matrix, check covar.info')
                message(cond$message)
                message('Check phenotype ', x)
            },
            warning = function(cond) {
            },
            finally = {
            }
        )
    }
            
    
    ###########################################################################
    #
    # lod.peaks
    #
    ###########################################################################
    cat("Checking lod.peaks\n")
    
    if ('lod.peaks' %not in% names(ds)) {
        message('lod.peaks not found')
    }
    
    if (!is.list(ds$lod.peaks)) {
        message(paste0("lod.peaks should be a list, but found ", class(ds$lod.peaks)))
    }
    
    if ('additive' %not in% names(ds$lod.peaks)) {
        message("additive should be an element in lod.peaks")
    }
    
    #
    # lod.peaks
    #
    for (lp in names(ds$lod.peaks)) {
        if (is.phenotype) {
            if (length(setdiff(ds$lod.peaks[[lp]]$data.name, ds$annot.phenotype$data.name))) {
                message(paste0('not all lod.peaks$data.name are in annot.phenotype$data.name'))
            }
        } else if (is.protein) {
            if (length(setdiff(ds$lod.peaks[[lp]]$protein.id, ds$annot.protein$protein.id))) {
                message(paste0('not all lod.peaks$protein.id are in annot.protein$protein.id'))
            }
        } else if (is.mrna) {
            if (length(setdiff(ds$lod.peaks[[lp]]$gene.id, ds$annot.mrna$gene.id))) {
                message(paste0('not all lod.peaks$gene.id are in annot.mrna$gene.id'))
            }
        } 
        
        if (length(setdiff(ds$lod.peaks[[lp]]$marker.id, markers$marker.id))) {
            message(paste0('not all ', d, '$lod.peaks$marker.id are in markers'))
        } 
    }
    
    
    
    ###########################################################################
    #
    # data
    #
    ###########################################################################
    cat("Checking data\n")
    
    num_annot_samples <- nrow(ds$annot.samples)
    
    if ('data' %not in% names(ds)) {
        message('data not found')
    }
    
    if (is.matrix(ds$data)) {
        # check if the data is numeric
        if (!is.numeric(ds$data)) {
            message(paste0(d, '$data matrix is not numeric'))
        }
        
        if (num_annot_samples != nrow(ds$data)) {
            message(paste0('number of samples (', num_annot_samples, ') != ',
                           'number of data rows (', nrow(ds$data), ')'))
        }
        
        
        
    } else if (is.list(ds$data)) {
        data.found <- FALSE
        if (any(c('rz','norm','log','transformed','raw') %in% tolower(names(ds$data)))) {
            data.found <- TRUE
        }
        
        if (!data.found) {
            message("'rz','norm','log','transformed', OR 'raw' NOT found in data, must be ONE of them")
        }
        
        data_list <- get('data', ds)
        data_names <- ls(data_list)
        for (i in 1:length(data_names)) {
            data_to_check <- paste0(d, '$data$', data_names[i])
            temp_data <- get(data_names[i], data_list)
            if (!is.numeric(temp_data)) {
                message(paste0(data_to_check, ' is not numeric'))
            }
            
            if (num_annot_samples != nrow(temp_data)) {
                message(paste0('number of samples (', num_annot_samples, ') != ',
                               'number of data rows (', nrow(temp_data), ')'))
            }
            rm(temp_data)
        }
    }
    
    rm(num_annot_samples)
    
}

perform_checks <- function() {
    # grab the datasets in the environment
    datasets <- grep('^dataset*', apropos('dataset\\.'), value=TRUE)

    if (gtools::invalid(datasets)) {
        message("No datasets found!")
    }
    
    if (!exists('ensembl.version')) {
        message("ensembl.version does not exist")
    }
    
    # Check genoprobs and K.
    if(length(genoprobs) != length(K)) {
        message(paste0("genoprobs (", length(genoprobs), ") and K (", length(K), ") do not have the same length."))
    } else {
        if(any(names(genoprobs) != names(K))) {
            message("names of genoprobs and K do not match.")
        }
        
        rownames.eq <- mapply(function(x, y) { all(rownames(x) == rownames(y)) }, genoprobs, K)
        if(any(rownames.eq == FALSE)) {
            message("sample IDs do not match between genoprobs and K.")
        }
    }
    
    # Check Marker IDs for genoprobs and map
    if(length(genoprobs) != length(map)) {
        message(paste0("genoprobs (", length(genoprobs), ") and map (", length(map), ") do not have the same length."))
    } else {
        rownames.eq <- mapply(function(x, y) { all(dimnames(x)[[3]] == names(y)) }, genoprobs, map)
        if(any(rownames.eq == FALSE)) {
            message("marker names do not match between genoprobs and map.")
        }
    }
    
    # Check dimensions of markers and map.
    map.length = sum(sapply(map, length))
    if(map.length != nrow(markers)) {
        message(paste("Number of rows in markers (", nrow(markers), ") does not equal the number of markers in map (", map.length, ")"))
    }        
    
    for (ds in datasets) {
        check_dataset(ds)
    }
}





