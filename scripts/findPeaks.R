findPeaks <- function(fName, dataset, start, step, threshold = 6.0, drop = 2, thresholdX = 6.0, dropX = 2, nCores = 0) {
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
        data  <- ds$annot.mrna[start:min(start+step-1, total), ]$gene.id
    } else if (tolower(ds$datatype) == "protein") {
        total <- dim(ds$annot.protein)[1]
        data  <- ds$annot.protein[start:min(start+step-1, total), ]$protein.id
    } else if (is_phenotype(ds)) {
        data  <- ds$annot.pheno[ds$annot.pheno$is.pheno==TRUE, ]
        data  <- data[data$omit==FALSE, ]
        total <- dim(data)[1]
        data  <- data[start:min(start+step-1, total), ]$data.name
    } else {
        stop("invalid datatype")
    }
    
    fileName <- paste0(start, '_', fName)
    
    first <- TRUE
    
    for (id in data) {
        print(paste0('GetLODScan: ', dataset, ', ', id))
        
        # get the lod scan
        lodScoresAdditive <- get_lod_scan_new(dataset, id, cores = nCores, 
                                           filter_threshold = threshold, filter_peak_drop = drop,
                                           filter_thresholdX = thresholdX, filter_peakdropX = dropX,
                                           scan1_output = TRUE)
        lodScoresAdditive$lod_scores$scan <- 'additive'
        lodScoresAll <- lodScoresAdditive$lod_scores
        
        # find the peaks
        peaksAdditive <- qtl2::find_peaks(lodScoresAdditive$scan1, map, threshold = threshold, drop = drop)
        
        if (gtools::invalid(peaksAdditive)) {
            peaksAdditive <- tibble(
                lodindex = numeric(),
                lodcolumn = character(),
                chr = character(),
                pos = numeric(),
                lod = numeric(),
                scan = character()
            )            
        } else {
            peaksAdditive$scan = 'additive'
        }
        
        peaksAll <- peaksAdditive
        
        # loop through interactive covariates
        for (i in 1:nrow(ds$covar.info)) {
            inf <- ds$covar.info[i, ]
            if (inf$interactive) {
                
                print(paste0('GetLODScan: ', dataset, ', ', id, ', ', inf$sample.column))
                lodScoresCovar <- get_lod_scan_new(dataset, id, intcovar = inf$sample.column, cores = nCores,
                                                filter_threshold = threshold, filter_peak_drop = drop,
                                                filter_thresholdX = thresholdX, filter_peakdropX = dropX,
                                                scan1_output = TRUE)
                lodScoresCovar$lod_scores$scan <- inf$sample.column
                
                # DO NOT STORE LOD, SCORE LOD DIFF
                temp <- lodScoresAll %>% 
                    inner_join(lodScoresCovar$lod_scores, by=c('id' = 'id')) %>% 
                    mutate(lod = lod.y - lod.x) %>% 
                    select(id, chr=chr.x, pos=pos.x, lod)
                
                
                lodScoresAll <- lodScoresAll %>% bind_rows(temp)
                
                # peaks
                temp <- lodScoresCovar$scan1 - lodScoresAdditive$scan1
                peaksCovar <- qtl2::find_peaks(temp, map, threshold = threshold, drop = drop)
                
                if(!gtools::invalid(peaksCovar)) {
                    peaksCovar$scan = inf$sample.column
                    peaksAll <- peaksAll %>% bind_rows(peaksCovar)
                }
            }
        }
        
        if (nrow(peaksAll) > 0) {
            lod_peaks <- inner_join(peaksAll,
                                    markers, 
                                    by = c("chr", "pos")) %>%
                select(annot.id = lodcolumn, marker.id, chr, pos, lod, scan) %>%
                mutate_at(c("lod"), as.numeric) %>% 
                as_tibble
            
            # get the allele effects only for additive
            ds_data <- get_data(dataset)
            
            for (i in 1:nrow(lod_peaks)) {
                t_score <- lod_peaks[i, ]
                if(t_score$scan == 'additive') {
                    t_chr <- t_score$chr
                    t_marker <- t_score$marker.id
                    t_probs <- list()
                    
                    t_probs[[t_chr]] <- genoprobs[[t_chr]][,,t_marker, drop = FALSE]
                    
                    af <- scan1blup(genoprobs = t_probs,
                                    pheno     = ds_data[, id, drop = FALSE],
                                    kinship   = K[[t_chr]])
                    
                    output <- lod_peaks[i, ] %>% bind_cols(as_tibble(t(af[, LETTERS[1:8]])))
                } else {
                    output <- lod_peaks[i, ] %>% bind_cols(tibble(A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0))
                }
                
                readr::write_csv(output, fileName, col_names = first, append = !first)
                first <- FALSE
            }
        }    
    }
}


# args:
# 1 == RData file
# 2 == QTLAPI source
# 3 == dataset
# 4 == start,
# 5 == step size
# 6 == number of cores

#args <- commandArgs(trailingOnly = TRUE)
#print(paste0('Loading: ', args[1]))
#load(args[1])
#debug_mode <- TRUE
#print(paste0('Sourcing: ', args[2]))
#source(args[2])
#findPeaks(args[3], args[4], strtoi(args[5]), strtoi(args[6]))
#print('DONE')


#peaks <- find_peaks('deleteme3', 'dataset.DOheart.mrna',  1, 3, 5)
