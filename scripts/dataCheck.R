library(gtools)

#' Opposite of %in%  
`%not in%` <- function (x, table) match(x, table, nomatch = 0L) == 0L


GetDatasets <- function() {
    # grab the datasets in the environment
    datasets <- grep('^dataset*', apropos('dataset\\.'), value=TRUE)
    return(datasets)
}


CheckMRNAProtein <- function(ds) {
    print("Checking for mRNA data")

    # annots
    if ('annots' %not in% names(ds)) {
        stop('annots not found')
    }
    
    isProtein <- FALSE
    if (ds$datatype == 'protein') {
        isProtein <- TRUE
    }
    
    if (!is.data.frame(ds$annots)) {
        stop(paste0("annots should be a data.frame, but found ", class(ds$annots)))
    }
    

    # covar
    if ('covar' %not in% names(ds)) {
        stop('covar not found')
    }
    
    if (!is.matrix(ds$covar)) {
        stop(paste0("covar should be a matrix, but found ", class(ds$covar)))
    }
    
    # covar.factors
    if ('covar.factors' %not in% names(ds)) {
        stop('covar.factors not found')
    }
    
    if (!is.data.frame(ds$covar.factors)) {
        stop(paste0("covar.factors should be a data.frame, but found ", class(ds$covar.factors)))
    }
    
    if ('column.name' %not in% names(ds$covar.factors)) {
        stop('column.name not found in covar.factors')
    } else if ('display.name' %not in% names(ds$covar.factors)) {
        stop('display.name not found in annots')
    } else if ('int.covar' %not in% names(ds$covar.factors)) {
        stop('int.covar not found in annots')
    } else if ('lod.peaks' %not in% names(ds$covar.factors)) {
        stop('lod.peaks not found in annots')
    } else if ('covar.name' %not in% names(ds$covar.factors)) {
        stop('covar.name not found in annots')
    }
    
    # display.name
    if ('display.name' %not in% names(ds)) {
        warning('display.name not found, will use dataset name')
    }
    
    # lod.peaks
    if ('lod.peaks' %not in% names(ds)) {
        stop('lod.peaks not found')
    }
    
    if (!is.list(ds$lod.peaks)) {
        stop(paste0("pheno should be a list, but found ", class(ds$lod.peaks)))
    }
    
    if ('additive' %not in% names(ds$lod.peaks)) {
        stop("additive should be an element in lod.peaks")
    }
    
    # pheno
    if ('pheno' %not in% names(ds)) {
        stop('pheno not found')
    }
    
    if (!is.data.frame(ds$pheno)) {
        stop(paste0("pheno should be a data.frame, but found ", class(ds$pheno)))
    }
    
    # samples
    if ('samples' %not in% names(ds)) {
        stop('samples not found')
    }
    
    if (!is.data.frame(ds$samples)) {
        stop(paste0("samples should be a data.frame, but found ", class(ds$samples)))
    }
    
    if ('mouse.id' %not in% names(ds$samples)) {
        stop('mouse.id not found in samples')
    }
    
    #
    # all the names exist, make sure they match up
    #
    
    # covar should be samples (rows) x covariates (columns)
    if (nrow(ds$covar) != nrow(ds$samples)) {
        stop(paste0('The number of rows in covar (', nrow(ds$covar), ') != rows in samples (', nrow(ds$samples), ')'))
    }
    
    # pheno should be samples (rows) x annots (columns)
    if (nrow(ds$pheno) != nrow(ds$samples)) {
        stop(paste0('The number of rows in pheno (', nrow(ds$pheno), ') != rows in samples (', nrow(ds$samples), ')'))
    }
    
    if (ncol(ds$pheno) != nrow(ds$annots)) {
        stop(paste0('The number of columns in pheno (', ncol(ds$pheno), ') != rows in annots (', nrow(ds$annots), ')'))
    }
    
    # samples rownames should be the same as mouse.id in the data
    if (!(all(rownames(ds$samples) == ds$samples[,'mouse.id']))) {
        stop('rownames of samples should be the same as samples$mouse.id')
    }
    
    #
    # some more complex rules
    #
    
    for(i in 1:nrow(ds$covar.factors)) {
        row <- ds$covar.factors[i,]
        
        if (row$column.name %not in% colnames(ds$samples)) {
            stop(paste0('covar.factors$column.name (', row$column.name, ') is not a column name in samples'))
        }
        
        if (!gtools::invalid(row$int.covar)) {
            if (row$int.covar %not in% c('factor', 'numeric')) {
                stop('covar.factors$int.covar needs to be NA, factor, or numeric')
            }
        }
        
        if (!gtools::invalid(row$lod.peaks)) {
            if (row$lod.peaks %not in% names(ds$lod.peaks)) {
                stop(paste0('covar.factors$lod.peaks (', row$lod.peaks, ') is not in lod.peaks'))
            }
        }
        
        if (!gtools::invalid(row$covar.name)) {
            if (row$covar.name %not in% colnames(ds$covar)) {
                stop(paste0('covar.factors$covar.name (', row$covar.name, ') is not a column in covar'))
            }
        }
    }
    
    for (lp in names(ds$lod.peaks)) {
        if (length(setdiff(ds$lod.peaks[[lp]]$annot.id, ds$annots$R_name))) {
            stop(paste0('not all lod.peaks$annot.id are in annots$R_name'))
        } 
        
        if (length(setdiff(ds$lod.peaks[[lp]]$marker_id, markers$marker))) {
            stop(paste0('not all lod.peaks$marker_id are in markers'))
        } 
    }
    
    
    
    
    
}
CheckProtein <- function(ds) {
    print("Checking for protein data")
    
}
CheckPhenotype <- function(ds) {

    
    
    
    
    
    
}

CheckDataset <- function(d) {
    # get the dataset
    ds <- NA
    
    if (exists(d)) {
        ds <- get(d)
    }

    if (gtools::invalid(ds)) {
        stop(paste0("dataset not found: ", d))
    } 
    
    print(paste0("Checking dataset: ", d))
    
    if (!is.list(ds)) {
        stop(paste0("dataset should be a list, but found ", class(d)))
    }
    
    if (!('datatype' %in% names(ds))) {
        stop("dataset must contain datatype")
    }
    
    dataType = ds[['datatype']]
    isMRNA <- FALSE
    isProtein <- FALSE
    isPheno <- FALSE
    
    if (dataType == 'mRNA') {
        isMRNA <- TRUE
    } else if (dataType == 'protein') {
        isProtein <- TRUE
    } else if (dataType == 'phenotype') {
        isPheno <- TRUE
    } else {
        stop(paste0("datatype is ", dataType, ", but should be mRNA, protein, or phenotype"))
    }

    print(paste0("Checking for ", dataType, " data"))
    
    # annots
    if ('annots' %not in% names(ds)) {
        stop('annots not found')
    }
    
    if (!is.data.frame(ds$annots)) {
        stop(paste0("annots should be a data.frame, but found ", class(ds$annots)))
    }
    
    if (isPheno) {
        if ('data_name' %not in% names(ds$annots)) {
            stop('data_name not found in annots')
        } else if ('short_name' %not in% names(ds$annots)) {
            warning('short_name not found in annots')
        } else if ('R_name' %not in% names(ds$annots)) {
            stop('R_name not found in annots')
        } else if ('description' %not in% names(ds$annots)) {
            stop('description not found in annots')
        } else if ('units' %not in% names(ds$annots)) {
            warning('units not found in annots')
        } else if ('is_id' %not in% names(ds$annots)) {
            stop('is_id not found in annots')
        } else if ('category' %not in% names(ds$annots)) {
            warning('category not found in annots')
        } else if ('R_category' %not in% names(ds$annots)) {
            warning('R_category not found in annots')
        } else if ('is_numeric' %not in% names(ds$annots)) {
            stop('is_numeric not found in annots')
        } else if ('is_date' %not in% names(ds$annots)) {
            warning('is_date not found in annots')
        } else if ('is_factor' %not in% names(ds$annots)) {
            stop('is_factor not found in annots')
        } else if ('factor_levels' %not in% names(ds$annots)) {
            warning('factor_levels not found in annots')
        } else if ('is_covar' %not in% names(ds$annots)) {
            stop('is_covar not found in annots')
        } else if ('is_pheno' %not in% names(ds$annots)) {
            stop('is_pheno not found in annots')
        } else if ('is_derived' %not in% names(ds$annots)) {
            warning('is_derived not found in annots')
        } else if ('omit' %not in% names(ds$annots)) {
            warning('omit not found in annots')
        } else if ('use_covar' %not in% names(ds$annots)) {
            stop('use_covar not found in annots')
        }
        
        # one and only 1 ID
        a_id <- ds$annots[which(ds$annots$is_id == TRUE),]$R_name    
        
        if(length(a_id) != 1) {
            stop('annots$is_id should have 1 and only 1 row set to TRUE')
        }
    } else {
        if (isProtein) {
            if ('protein_id' %not in% names(ds$annots)) {
                stop('protein_id not found in annots')
            }
        }
        
        if ('gene_id' %not in% names(ds$annots)) {
            stop('gene_id not found in annots')
        } else if ('symbol' %not in% names(ds$annots)) {
            stop('symbol not found in annots')
        } else if ('chr' %not in% names(ds$annots)) {
            stop('chr not found in annots')
        } else if ('start' %not in% names(ds$annots)) {
            stop('start not found in annots')
        } else if ('end' %not in% names(ds$annots)) {
            warning('end not found in annots')
        } else if ('strand' %not in% names(ds$annots)) {
            stop('strand not found in annots')
        } else if ('middle' %not in% names(ds$annots)) {
            warning('middle not found in annots')
        } else if ('nearest.marker.id' %not in% names(ds$annots)) {
            warning('nearest.marker.id not found in annots')
        }
        
        if (isProtein) {
            if (!(all(rownames(ds$annots) == ds$annots[,'protein_id']))) {
                stop('rownames of annots should be the same as annots$protein_id')
            }
        } else {
            if (!(all(rownames(ds$annots) == ds$samples[,'gene_id']))) {
                stop('rownames of annots should be the same as annots$gene_id')
            }
        }
    }
    
    # covar
    if ('covar' %not in% names(ds)) {
        stop('covar not found')
    }
    
    if (!is.matrix(ds$covar)) {
        stop(paste0("covar should be a matrix, but found ", class(ds$covar)))
    }
    
    # covar.factors
    if ('covar.factors' %not in% names(ds)) {
        stop('covar.factors not found')
    }
    
    if (!is.data.frame(ds$covar.factors)) {
        stop(paste0("covar.factors should be a data.frame, but found ", class(ds$covar.factors)))
    }
    
    if ('column.name' %not in% names(ds$covar.factors)) {
        stop('column.name not found in covar.factors')
    } else if ('display.name' %not in% names(ds$covar.factors)) {
        stop('display.name not found in annots')
    } else if ('int.covar' %not in% names(ds$covar.factors)) {
        stop('int.covar not found in annots')
    } else if ('lod.peaks' %not in% names(ds$covar.factors)) {
        stop('lod.peaks not found in annots')
    } else if ('covar.name' %not in% names(ds$covar.factors)) {
        stop('covar.name not found in annots')
    }
    
    # display.name
    if ('display.name' %not in% names(ds)) {
        warning('display.name not found, will use dataset name')
    }
    
    # lod.peaks
    if ('lod.peaks' %not in% names(ds)) {
        stop('lod.peaks not found')
    }
    
    if (!is.list(ds$lod.peaks)) {
        stop(paste0("pheno should be a list, but found ", class(ds$lod.peaks)))
    }
    
    if ('additive' %not in% names(ds$lod.peaks)) {
        stop("additive should be an element in lod.peaks")
    }
    
    if (isPheno) {
        # pheno
        if ('pheno' %not in% names(ds)) {
            stop('pheno not found')
        }
        
        if (!is.data.frame(ds$pheno)) {
            stop(paste0("pheno should be a data.frame, but found ", class(ds$pheno)))
        }
    } else {
        if ('rankz' %not in% names(ds)) {
            stop('rankz not found')
        }
        
        if (!is.matrix(ds$rankz)) {
            stop(paste0("rankz should be a matrix, but found ", class(ds$rankz)))
        }
        
    }
    
    
    # samples
    if ('samples' %not in% names(ds)) {
        stop('samples not found')
    }
    
    if (!is.data.frame(ds$samples)) {
        stop(paste0("samples should be a data.frame, but found ", class(ds$samples)))
    }
    
    if ('mouse.id' %not in% names(ds$samples)) {
        stop('mouse.id not found in samples')
    }
    
    #
    # all the names exist, make sure they match up
    #
    
    # covar should be samples (rows) x covariates (columns)
    if (nrow(ds$covar) != nrow(ds$samples)) {
        stop(paste0('The number of rows in covar (', nrow(ds$covar), ') != rows in samples (', nrow(ds$samples), ')'))
    }
    
    if (isPheno) {
        # pheno should be samples (rows) x annots (columns)
        if (nrow(ds$pheno) != nrow(ds$samples)) {
            stop(paste0('The number of rows in pheno (', nrow(ds$pheno), ') != rows in samples (', nrow(ds$samples), ')'))
        }
        if (ncol(ds$pheno) != nrow(ds$annots)) {
            stop(paste0('The number of columns in pheno (', ncol(ds$pheno), ') != rows in annots (', nrow(ds$annots), ')'))
        }
    } else {
        # rankz should be samples (rows) by annots (columns)
        if (nrow(ds$rankz) != nrow(ds$samples)) {
            stop(paste0('The number of rows in rankz (', nrow(ds$rankz), ') != rows in samples (', nrow(ds$samples), ')'))
        }
        if (ncol(ds$rankz) != nrow(ds$annots)) {
            stop(paste0('The number of columns in rankz (', ncol(ds$rankz), ') != rows in annots (', nrow(ds$annots), ')'))
        }
    }

    # samples rownames should be the same as mouse.id in the data
    if (!(all(rownames(ds$samples) == ds$samples[,'mouse.id']))) {
        stop('rownames of samples should be the same as samples$mouse.id')
    }
    
    #
    # some more complex rules
    #
    
    for(i in 1:nrow(ds$covar.factors)) {
        row <- ds$covar.factors[i,]
        
        if (row$column.name %not in% colnames(ds$samples)) {
            stop(paste0('covar.factors$column.name (', row$column.name, ') is not a column name in samples'))
        }
        
        if (!gtools::invalid(row$int.covar)) {
            if (row$int.covar %not in% c('factor', 'numeric')) {
                stop('covar.factors$int.covar needs to be NA, factor, or numeric')
            }
        }
        
        if (!gtools::invalid(row$lod.peaks)) {
            if (row$lod.peaks %not in% names(ds$lod.peaks)) {
                stop(paste0('covar.factors$lod.peaks (', row$lod.peaks, ') is not in lod.peaks'))
            }
        }
        
        if (!gtools::invalid(row$covar.name)) {
            if (row$covar.name %not in% colnames(ds$covar)) {
                stop(paste0('covar.factors$covar.name (', row$covar.name, ') is not a column in covar'))
            }
        }
    }
    
    for (lp in names(ds$lod.peaks)) {
        if (isPheno) {
            if (length(setdiff(ds$lod.peaks[[lp]]$annot.id, ds$annots$R_name))) {
                stop(paste0('not all lod.peaks$annot.id are in annots$R_name'))
            }
        } else if (isProtein) {
            if (length(setdiff(ds$lod.peaks[[lp]]$annot.id, ds$annots$protein_id))) {
                stop(paste0('not all lod.peaks$annot.id are in annots$R_name'))
            }
        } else if (isMRNA) {
            if (length(setdiff(ds$lod.peaks[[lp]]$annot.id, ds$annots$gene_id))) {
                stop(paste0('not all lod.peaks$annot.id are in annots$R_name'))
            }
        } 
        
        if (length(setdiff(ds$lod.peaks[[lp]]$marker_id, markers$marker))) {
            stop(paste0('not all lod.peaks$marker_id are in markers'))
        } 
    }


    # Check Sample IDs for genoprobs and K.
    if(length(genoprobs) != length(K)) {
        stop(paste0("genoprobs (", length(genoprobs), ") and K (", length(K), ") do not have the same length."))
    } else {
        if(any(names(genoprobs) != names(K))) {
            stop("names of genoprobs and K do not match.")
        }
        
        rownames.eq <- mapply(function(x, y) { all(rownames(x) == rownames(y)) }, genoprobs, K)
        if(any(rownames.eq == FALSE)) {
            stop("sample IDs do not match between genoprobs and K.")
        }
    }
    
    # Check Marker IDs for genoprobs and map.
    if(length(genoprobs) != length(map)) {
        stop(paste0("genoprobs (", length(genoprobs), ") and map (", length(map), ") do not have the same length."))
    } else {
        rownames.eq <- mapply(function(x, y) { all(dimnames(x)[[3]] == names(y)) }, genoprobs, map)
        if(any(rownames.eq == FALSE)) {
            stop("marker names do not match between genoprobs and map.")
        }
    }
    
    # Check dimensions of markers and map.
    map.length = sum(sapply(map, length))
    if(map.length != nrow(markers)) {
        stop(paste("Number of rows in markers (", nrow(markers), ") does not equal the number of markers in map (", map.length, ")",))
    }    
        
    
    
    

}

PerformChecks <- function() {
    datasets <- GetDatasets()
    
    if (gtools::invalid(datasets)) {
        stop("No datasets found!")
    }
    
    for (ds in datasets) {
        CheckDataset(ds)
    }
    
    
}

PerformChecks()
