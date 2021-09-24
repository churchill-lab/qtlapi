findPeaks <- function(fName, dataset, start, step, threshold = 6.0, drop = 2, nCores = 0) {
    print(paste0('fName: ', fName))
    print(paste0('dataset: ', dataset))
    print(paste0('start: ', start))
    print(paste0('step: ', step))
    print(paste0('threshold: ', threshold))
    print(paste0('drop: ', drop))
    
    output <- FALSE
    first <- TRUE
    
    # get the dataset
    ds = get_dataset(dataset)

    if (is.null(ds)) {
        stop(paste0("dataset not found: ", dataset))
    }

    # get the index and data
    if (tolower(ds$datatype) == "mrna") {
        total <- dim(ds$annot.mrna)[1]
        data <- ds$annot.mrna[start:min(start+step-1, total), ]$gene.id
    } else if (tolower(ds$datatype) == "protein") {
        total <- dim(ds$annot.protein)[1]
        data <- ds$annot.protein[start:min(start+step-1, total), ]$protein.id
    } else if (is_phenotype(ds)) {
        data <- ds$annot.pheno[ds$annot.pheno$is.pheno==TRUE, ]
        data <- data[data$omit==FALSE, ]
        total <- dim(data)[1]
        data <- data[start:min(start+step-1, total), ]$data.name
    } else {
        stop("invalid datatype")
    }

    fileName <- paste0(start, '_', fName)

    for (id in data) {
        print(paste0('GetLODScan: ', dataset, ', ', id))

        # get the lod scan
        lodScores <- get_lod_scan(dataset, id, cores = nCores, 
                                  filter_threshold = threshold, filter_peak_drop = drop)
        lodScores$scan <- 'additive'
        lodScoresAdditive <- lodScores
        
        covar_columns <- c('additive')
        
        
        # loop through interactive covariates
        for (i in 1:nrow(ds$covar.info)) {
            inf <- ds$covar.info[i, ]
            if (inf$interactive) {
                
                print(paste0('GetLODScan: ', dataset, ', ', id, ', ', inf$sample.column))
                lodScoresCovar <- get_lod_scan(dataset, id, intcovar = inf$sample.column, cores = nCores)
                
                # DO NOT STORE LOD, SCORE LOD DIFF
                temp <- lodScoresAdditive %>% 
                    inner_join(lodScoresCovar, by=c('id' = 'id')) %>% 
                    mutate(lod = lod.y - lod.x) %>% 
                    select(id, chr=chr.x, pos=pos.x, lod)
                
                temp$scan <- inf$sample.column
                covar_columns <- c(covar_columns, inf$sample.column)
                
                lodScores <- lodScores %>% bind_rows(temp)
            }
        }
    
    
        maxLodScores <- lodScores %>% 
            group_by(scan, chr) %>%     # group by
            top_n(1, lod) %>%           # greatest lod score
            filter(lod > 6.0) %>%       # filter
            arrange(chr, pos) %>%       # sort
            data.frame() %>%            # remove the grouping variables stored by tibble
            as_tibble()                 # back to a tibble
        
        newLodScores <- maxLodScores %>% 
            select(id) %>% 
            inner_join(lodScores %>% filter(scan == 'additive'), by=c('id' = 'id')) %>% 
            select(id, additive = lod)
        
        for (i in 1:nrow(ds$covar.info)) {
            inf <- ds$covar.info[i, ]
            if (inf$interactive) {
                
                temp <- newLodScores %>% 
                    select(id) %>% 
                    inner_join(lodScores %>% filter(scan == inf$sample.column), by=c('id' = 'id')) %>%
                    select(id, UQ(as.name(inf$sample.column)) := lod)            
                
                newLodScores <- newLodScores %>% bind_cols(temp[,inf$sample.column])
            }
        }
        
    
        if (nrow(newLodScores) > 0) {
            allLods <- newLodScores %>% 
                add_column(annot.id = id, .before=1) %>%
                rename(marker.id = id)
            
            allLodsFull <- allLods %>%
                inner_join(markers, by = c('marker.id' = 'marker.id')) %>%
                select(annot.id, marker.id, chr, pos, additive)
            
            # get the allele effects only for additive
            ds_data <- get_data(dataset)
            
            for (i in 1:nrow(allLodsFull)) {
                t_score <- allLodsFull[i, ]
                t_chr <- t_score$chr
                t_marker <- t_score$marker.id
                t_probs <- list()
                
                t_probs[[t_chr]] <- genoprobs[[t_chr]][,,t_marker, drop = FALSE]
                
                af <- scan1blup(genoprobs = t_probs,
                                pheno     = ds_data[, id, drop = FALSE],
                                kinship   = K[[t_chr]])
                
                if (first) {
                    first <- FALSE
                    output <- allLods[i, ] %>% bind_cols(as_tibble(t(af[, LETTERS[1:8]])))
                } else {
                    temp_out <- allLods[i, ] %>% bind_cols(as_tibble(t(af[, LETTERS[1:8]])))
                    output <- output %>% bind_rows(temp_out)
                }
            }
        }    
    }
    
    readr::write_csv(output, fileName, col_names = TRUE, append = FALSE)
}


# args:
# 1 == RData file
# 2 == QTLAPI source
# 3 == dataset
# 4 == start,
# 5 == step size
# 6 == number of cores

args <- commandArgs(trailingOnly = TRUE)
print(paste0('Loading: ', args[1]))
load(args[1])
debug_mode <- TRUE
print(paste0('Sourcing: ', args[2]))
source(args[2])
find_peaks(args[3], args[4], strtoi(args[5]), strtoi(args[6]))
print('DONE')


#peaks <- find_peaks('deleteme3', 'dataset.DOheart.mrna',  1, 3, 5)
