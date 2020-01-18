find_peaks <- function(fName, dataset, start, step, nCores = 0) {
    print(paste0('fName: ', fName))
    print(paste0('dataset: ', dataset))
    print(paste0('start: ', start))
    print(paste0('step: ', step))

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
        lodScores <- get_lod_scan(dataset, id, cores = nCores)

        maxLodScores <- lodScores %>% 
                        group_by(chr) %>% 
                        top_n(1, lod) %>%
                        filter(lod > 6.0)

        if (nrow(maxLodScores) > 0) {
            # get the allele effects only for additive
            print('getting data...')
            ds_data <- get_data(dataset)

            for (i in 1:nrow(maxLodScores)) {
                t_score <- maxLodScores[i, ]
                t_chr <- t_score$chr
                t_marker <- t_score$id
                t_probs <- list()
                t_probs[[t_chr]] <- genoprobs[[t_chr]][,,t_marker, drop = FALSE]

                af <- scan1blup(genoprobs = t_probs,
                                pheno     = ds_data[, id, drop = FALSE],
                                kinship   = K[[t_chr]])
    

                output <- data.frame(intCovar  = 'additive', 
                                     annot.id  = id, 
                                     marker.id = maxLodScores$id, 
                                     lod       = maxLodScores$lod,
                                     t(af[, LETTERS[1:8]]))

                write.table(output,
                            fileName,
                            sep = ",", 
                            row.names = FALSE, 
                            col.names = FALSE, 
                            append = TRUE)
            }
        }

        # loop through interactive covariates
        for (i in 1:nrow(ds$covar.info)) {
            inf <- ds$covar.info[i, ]
            if (inf$interactive) {

                print(paste0('GetLODScan: ', dataset, ', ', id, ', ', inf$sample.column))
                lodScoresCovar <- get_lod_scan(dataset, id, intcovar = inf$sample.column, cores = nCores)

                # DO NOT STORE LOD, SCORE LOD DIFF
                temp <- merge(x=lodScores, y=lodScoresCovar[,c("id", "lod")], by="id", all.x=TRUE)
                temp$lod <- temp$lod.y - temp$lod.x
                temp$lod.x <- NULL
                temp$lod.y <- NULL

                maxLodScoresCovar <- temp %>% 
                                     group_by(chr) %>% 
                                     top_n(1, lod) %>%
                                     filter(lod > 6.0)

                if (nrow(maxLodScoresCovar) > 0) {
                    maxLodScoresCovar <- data.frame(intCovar=inf$sample.column, annot.id=id, marker.id=maxLodScoresCovar$id, lod=maxLodScoresCovar$lod)
                    write.table(maxLodScoresCovar,
                                fileName,
                                sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
                }
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

args <- commandArgs(trailingOnly = TRUE)
print(paste0('Loading: ', args[1]))
load(args[1])
debug_mode <- TRUE
print(paste0('Sourcing: ', args[2]))
source(args[2])
find_peaks(args[3], args[4], strtoi(args[5]), strtoi(args[6]))
print('DONE')