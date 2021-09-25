findPeaks <- function(fName, dataset, start, step, threshold = 6.0, peakdrop = 2, thresholdX = 6.0, peakdropX = 2, nCores = 0) {
    print(paste0('fName: ', fName))
    print(paste0('dataset: ', dataset))
    print(paste0('start: ', start))
    print(paste0('step: ', step))
    print(paste0('threshold: ', threshold))
    print(paste0('peakdrop: ', peakdrop))
    
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
                                           filter_threshold = threshold, filter_peak_drop = peakdrop,
                                           filter_thresholdX = thresholdX, filter_peakdropX = peakdropX,
                                           scan1_output = TRUE)
        lodScoresAdditive$lod_scores$scan <- 'additive'
        lodScoresAll <- lodScoresAdditive$lod_scores
        
        # find the peaks
        peaksAdditive <- qtl2::find_peaks(lodScoresAdditive$scan1, map, 
                                          threshold = threshold, peakdrop = peakdrop)
        
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
                                                filter_threshold = threshold, filter_peak_drop = peakdrop,
                                                filter_thresholdX = thresholdX, filter_peakdropX = peakdropX,
                                                scan1_output = TRUE)
                lodScoresCovar$lod_scores$scan <- inf$sample.column
                
                # DO NOT STORE LOD, SCORE LOD DIFF
                temp <- lodScoresAll %>% 
                    inner_join(lodScoresCovar$lod_scores, by=c('id' = 'id')) %>% 
                    mutate(lod = lod.y - lod.x) %>% 
                    select(id, chr=chr.x, pos=pos.x, lod, scan=scan.y)
                
                
                lodScoresAll <- lodScoresAll %>% bind_rows(temp)
                
                # peaks
                temp <- lodScoresCovar$scan1 - lodScoresAdditive$scan1
                peaksCovar <- qtl2::find_peaks(temp, map, threshold = threshold, peakdrop = peakdrop)
                
                if(!gtools::invalid(peaksCovar)) {
                    peaksCovar$scan = inf$sample.column
                    peaksAll <- peaksAll %>% bind_rows(peaksCovar)
                }
            }
        }
        
        if (nrow(peaksAll) > 0) {
            lod_peaks <- left_join(peaksAll,
                                   markers, 
                                   by = c('chr', 'pos')) %>%
                select(annot.id = lodcolumn, marker.id, chr, pos, lod, scan) %>%
                mutate_at(c("lod"), as.numeric) %>% 
                as_tibble
            
            ds_data <- get_data(dataset)
            
            for (i in 1:nrow(lod_peaks)) {
                peak <- lod_peaks[i, ]

                if(gtools::invalid(peak$marker.id)) {
                    mrk_id <- qtl2::find_marker(map, peak$chr, peak$pos)
                    mrk <- markers %>% filter(marker.id == mrk_id)
                    peak$marker.id <- mrk_id
                    peak$pos <- mrk$pos
                    lod_peaks[i, ] <- peak
                }

                # get the allele effects only for additive
                if(peak$scan == 'additive') {
                    temp_probs <- list()
                    temp_probs[[peak$chr]] <- genoprobs[[peak$chr]][,,peak$marker.id, drop = FALSE]
                    
                    
                    af <- scan1blup(genoprobs = temp_probs,
                                    pheno     = ds_data[, id, drop = FALSE],
                                    kinship   = K[[peak$chr]])
                    
                    output <- lod_peaks[i, ] %>% bind_cols(as_tibble(t(af[, LETTERS[1:8]])))
                } else {
                    output <- lod_peaks[i, ] %>% bind_cols(tibble(A=NA,B=NA,C=NA,D=NA,E=NA,F=NA,G=NA,H=NA))
                }
                
                readr::write_csv(output, fileName, col_names = first, append = !first, na = '')
                first <- FALSE
            }
        }    
    }
}


# args:
# 1 == RData file
# 2 == QTLAPI source
# 3 == filename
# 4 == dataset
# 5 == start
# 6 == step size
# 7 == number of cores

args <- commandArgs(trailingOnly = TRUE)
print(paste0('Loading: ', args[1]))
load(args[1])
debug_mode <- TRUE
print(paste0('Sourcing: ', args[2]))
source(args[2])
findPeaks(args[3], args[4], strtoi(args[5]), strtoi(args[6]), nCores=strtoi(args[7]))
print('DONE')

#peaks <- find_peaks('deleteme3', 'dataset.DOheart.mrna',  1, 3, 5)
