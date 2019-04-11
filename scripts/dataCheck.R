
GetDatasets <- function() {
    # grab the datasets in the environment
    datasets <- grep('^dataset*', apropos('dataset\\.'), value=TRUE)
    return(datasets)
}

CheckDataset <- function(d) {
    # get the dataset
    ds <- NA
    
    if (exists(d)) {
        ds <- get(d)
    }

    if (gtools::invalid(ds)) {
        stop(paste0("dataset not found '", d, '"'))
    } 
    
    print(paste0("Checking dataset '", d, "'"))
    
    if (!is.list(ds)) {
        stop(paste0("dataset should be a list, but found '", class(d), "'"))
    }
    
    if ('datatype' %not in% names(ds)) {
        stop("dataset must contain 'datatype'")
    }
    
    datatype = ds[['datatype']]
    is.mrna <- FALSE
    is.protein <- FALSE
    is.pheno <- FALSE
    
    if (tolower(datatype) == 'mrna') {
        is.mrna <- TRUE
    } else if (tolower(datatype) == 'protein') {
        is.protein <- TRUE
    } else if (isPheno(ds)) {
        is.pheno <- TRUE
    } else {
        stop(paste0("datatype is '", datatype, "', but should be mRNA, protein, or phenotype"))
    }

    ###########################################################################
    #
    # annotations
    #
    ###########################################################################
    print("Checking annotations")
    
    if (is.mrna) {
        if ('annot.mrna' %not in% names(ds)) {
            stop(paste0("annot.mrna not found in '", d, "'"))
        }

        if (!is_tibble(ds$annot.mrna)) {
            stop(paste0("annot.mrna should be a tibble, but found '", class(ds$annot.mrna), "'"))
        }
        
        if (any(duplicated(ds$annot.mrna$gene.id))) {
            stop("There are duplicated gene.id annotations in annot.mrna")
        }
    } else if (is.protein) {
        if ('annot.protein' %not in% names(ds)) {
            stop(paste0("annot.protein not found in '", d, "'"))
        }
        
        if (!is_tibble(ds$annot.protein)) {
            stop(paste0("annot.protein should be a tibble, but found ", class(ds$annot.protein)))
        }
        
        if (any(duplicated(ds$annot.protein$protein.id))) {
            stop("There are duplicated protein.id annotations in annot.protein")
        }
    } else if (is.pheno) {
        if ('annot.pheno' %not in% names(ds)) {
            stop(paste0("annot.pheno not found in '", d, "'"))
        }

        if (!is_tibble(ds$annot.pheno)) {
            stop(paste0("annot.pheno should be a tibble, but found ", class(ds$annot.pheno)))
        }
        
        if (any(duplicated(ds$annot.pheno$data.name))) {
            stop("There are duplicated data.name annotations in annot.pheno")
        }
    }
    
    if (is.pheno) {
        if ('data.name' %not in% names(ds$annot.pheno)) {
            stop('data.name not found in annot.pheno')
        } else if ('short.name' %not in% names(ds$annot.pheno)) {
            warning('short.name not found in annot.pheno')
        } else if ('R.name' %not in% names(ds$annot.pheno)) {
            stop('R.name not found in annot.pheno')
        } else if ('description' %not in% names(ds$annot.pheno)) {
            stop('description not found in annot.pheno')
        } else if ('units' %not in% names(ds$annot.pheno)) {
            warning('units not found in annot.pheno')
        } else if ('is.id' %not in% names(ds$annot.pheno)) {
            stop('is.id not found in annot.pheno')
        } else if ('category' %not in% names(ds$annot.pheno)) {
            warning('category not found in annot.pheno')
        } else if ('R.category' %not in% names(ds$annot.pheno)) {
            warning('R.category not found in annot.pheno')
        } else if ('is.numeric' %not in% names(ds$annot.pheno)) {
            stop('is.numeric not found in annot.pheno')
        } else if ('is.date' %not in% names(ds$annot.pheno)) {
            warning('is.date not found in annot.pheno')
        } else if ('is.factor' %not in% names(ds$annot.pheno)) {
            stop('is.factor not found in annot.pheno')
        } else if ('factor.levels' %not in% names(ds$annot.pheno)) {
            warning('factor.levels not found in annot.pheno')
        } else if ('is.covar' %not in% names(ds$annot.pheno)) {
            stop('is.covar not found in annot.pheno')
        } else if ('is.pheno' %not in% names(ds$annot.pheno)) {
            stop('is.pheno not found in annot.pheno')
        } else if ('is.derived' %not in% names(ds$annot.pheno)) {
            warning('is.derived not found in annot.pheno')
        } else if ('omit' %not in% names(ds$annot.pheno)) {
            warning('omit not found in annot.pheno')
        } else if ('use.covar' %not in% names(ds$annot.pheno)) {
            stop('use.covar not found in annot.pheno')
        }
        
        # one and only 1 ID
        a.id <- ds$annot.pheno[which(ds$annot.pheno$is.id == TRUE),]$R.name    
        
        if(length(a.id) != 1) {
            stop('annot.pheno$is.id should have 1 and only 1 row set to TRUE')
        }
    } else {
        annot <- NULL
        annot.name <- NULL
        
        if (is.protein) {
            if ('protein.id' %not in% names(ds$annot.protein)) {
                stop('protein.id not found in annot.protein')
            }
            annot <- ds$annot.protein
            annot.name <- 'annot.protein'
        } else {
            annot <- ds$annot.mrna
            annot.name <- 'annot.mrna'
        }
        
        if ('gene.id' %not in% names(annot)) {
            stop(paste('gene.id not found in ', annot.name))
        } else if ('symbol' %not in% names(annot)) {
            warning(paste('symbol not found in ', annot.name))
        } else if ('chr' %not in% names(annot)) {
            warning(paste('chr not found in ', annot.name))
        } else if ('start' %not in% names(annot)) {
            stop(paste('start not found in ', annot.name))
        } else if ('end' %not in% names(annot)) {
            stop(paste('end not found in ', annot.name))
        } else if ('strand' %not in% names(annot)) {
            warning(paste('strand not found in ', annot.name))
        } else if ('middle' %not in% names(annot)) {
            warning(paste('strand not found in ', annot.name))
        } else if ('nearest.marker.id' %not in% names(annot)) {
            stop(paste('nearest.marker.id not found in ', annot.name))
        }
    }
    
    ###########################################################################
    #
    # samples
    #
    ###########################################################################
    print("Checking annot.samples")
    
    if ('annot.samples' %not in% names(ds)) {
        stop('annot.samples not found')
    }
    
    if (!is_tibble(ds$annot.sample)) {
        stop(paste0("annot.samples should be a tibble, but found ", class(ds$annot.sample)))
    }
    
    id.column <- match('mouse.id', tolower(colnames(ds$annot.samples)))
    if (gtools::invalid(id.column)) {
        stop('mouse.id not found in annot.samples')
    }
    
    
    ###########################################################################
    #
    # covar.info
    #
    ###########################################################################
    print("Checking covar.info")
    
    if ('covar.info' %not in% names(ds)) {
        stop('covar.info not found')
    }
    
    if (!is_tibble(ds$covar.info)) {
        stop(paste0("covar.info should be a tibble, but found '", class(ds$covar.info), "'"))
    }
    
    if ('sample.column' %not in% names(ds$covar.info)) {
        stop('sample.column not found in covar.info')
    } else if ('covar.column' %not in% names(ds$covar.info)) {
        stop('covar.column not found in covar.info')
    } else if ('display.name' %not in% names(ds$covar.info)) {
        stop('display.name not found in covar.info')
    } else if ('interactive' %not in% names(ds$covar.info)) {
        stop('interactive not found in covar.info')
    } else if ('primary' %not in% names(ds$covar.info)) {
        stop('primary not found in covar.info')
    } else if ('lod.peaks' %not in% names(ds$covar.info)) {
        stop('lod.peaks not found in covar.info')
    }
    
    for(i in 1:nrow(ds$covar.info)) {
        row <- ds$covar.info[i, ]
        
        if (row$sample.column %not in% colnames(ds$annot.samples)) {
            stop(paste0("covar.factors$sample.column ('", row$column.name, "') is not a column name in annot.samples"))
        }
        
        if (!any(grepl(row$covar.column, colnames(ds$covar.matrix)))) {
            stop(paste0("covar.info$covar.column ('", row$covar.column, "') is not a match to a column in covar.matrix"))
        }
        
        if (gtools::invalid(row$display.name)) {
            stop("covar.info$display.name needs to have a value")
        }
        
        if (row$interactive) {
            if (gtools::invalid(row$lod.peaks)) {
                stop("covar.info$interactive is TRUE, but covar.info$lod.peaks is NA")
            } else {
                # check for existance of lod.peaks
                if (gtools::invalid(ds$lod.peaks[[row$lod.peaks]])) {
                    stop(paste0("covar.info$interactive is TRUE, but covar.info$lod.peaks ('", row$lod.peaks, "') is not in lod.peaks"))
                }
            }
        } else {
            if (!gtools::invalid(row$lod.peaks)) {
                stop(paste0("covar.info$interactive is FALSE, but covar.info$lod.peaks ('", row$lod.peaks, "') is set"))
            }
        }
    }
    
    if (!any(ds$covar.info$primary)) {
        stop("covar.info$primary needs at least 1 value set to TRUE")
    }
    
    
    ###########################################################################
    #
    # covar.matrix
    #
    ###########################################################################
    print("Checking covar.matrix")
    
    if ('covar.matrix' %not in% names(ds)) {
        stop('covar.matrix not found')
    }
    
    if (!is.matrix(ds$covar.matrix)) {
        stop(paste0("covar.matrix should be a matrix, but found '", class(ds$covar.matrix), "'"))
    }
    
    ###########################################################################
    #
    # display.name
    #
    ###########################################################################
    print("Checking display.name")

    if ('display.name' %not in% names(ds)) {
        warning('display.name not found, will use dataset name')
    }
    
    ###########################################################################
    #
    # lod.peaks
    #
    ###########################################################################
    print("Checking lod.peaks")
    
    if ('lod.peaks' %not in% names(ds)) {
        stop('lod.peaks not found')
    }
    
    if (!is.list(ds$lod.peaks)) {
        stop(paste0("lod.peaks should be a list, but found ", class(ds$lod.peaks)))
    }
    
    if ('additive' %not in% names(ds$lod.peaks)) {
        stop("additive should be an element in lod.peaks")
    }
    
    #
    # lod.peaks
    #
    for (lp in names(ds$lod.peaks)) {
        if (is.pheno) {
            if (length(setdiff(ds$lod.peaks[[lp]]$data.name, ds$annot.pheno$data.name))) {
                stop(paste0('not all lod.peaks$data.name are in annot.pheno$data.name'))
            }
        } else if (is.protein) {
            if (length(setdiff(ds$lod.peaks[[lp]]$protein.id, ds$annot.protein$protein_id))) {
                stop(paste0('not all lod.peaks$protein.id are in annot.protein$protein.id'))
            }
        } else if (is.mrna) {
            if (length(setdiff(ds$lod.peaks[[lp]]$gene.id, ds$annot.mrna$gene.id))) {
                stop(paste0('not all lod.peaks$gene.id are in annot.mrna$gene.id'))
            }
        } 
        
        if (length(setdiff(ds$lod.peaks[[lp]]$marker.id, markers$marker.id))) {
            stop(paste0('not all lod.peaks$marker.id are in markers'))
        } 
    }
    
    
    
    ###########################################################################
    #
    # data
    #
    ###########################################################################
    print("Checking data")
    
    if ('data' %not in% names(ds)) {
        stop('data not found')
    }
    
    if (is.matrix(ds$data)) {
        
    } else if (is.list(ds$data)) {
        data.found <- FALSE
        if (any(c('rz','norm','log','transformed','raw') %in% tolower(names(ds$data)))) {
            data.found <- TRUE
        }
        
        if (!data.found) {
            stop("'rz','norm','log','transformed', OR 'raw' NOT found in data, must be ONE of them")
        }
    }
    
    
    
    

}

PerformChecks <- function() {
    datasets <- GetDatasets()
    
    if (gtools::invalid(datasets)) {
        stop("No datasets found!")
    }
    
    # Check genoprobs and K.
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
    
    # Check Marker IDs for genoprobs and map
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
        stop(paste("Number of rows in markers (", nrow(markers), ") does not equal the number of markers in map (", map.length, ")"))
    }        
    
    for (ds in datasets) {
        CheckDataset(ds)
    }
}

PerformChecks()
