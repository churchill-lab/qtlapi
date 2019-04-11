library(RestRserve)


prettyJSON <- function(timestamp, level, logger_name, pid, message) {
    x = to_json(
        list(
            timestamp = format(timestamp, "%Y-%m-%d %H:%M:%OS6"),
            level = as.character(level),
            name = as.character(logger_name),
            pid = as.integer(pid),
            message = message
        )
    )
    cat(jsonlite::prettify(x), file = "", append = TRUE, sep = "\n")
}

# setup app and logging
#logger <- RestRserve::Logger$new(level=TRACE, name="app", FUN=prettyJSON)
logger <- RestRserve::Logger$new(level=TRACE, name="app")
RestRserveApp <- RestRserve::RestRserveApplication$new()
RestRserveApp$logger <- logger


HttpGetDataset <- function(request, response) {
    result <- tryCatch({
        # start the clock
        ptm <- proc.time()
        
        datasets <- GetDatasetInfo()
        
        elapsed <- proc.time() - ptm
        #TrackTime(req, elapsed["elapsed"])
        #logger$trace(paste0('time: ', elapsed["elapsed"]))
        #logger$debug(paste0('time: ', elapsed["elapsed"]))
        #logger$info(paste0('time: ', elapsed["elapsed"]))
        #logger$warning(paste0('time: ', elapsed["elapsed"]))
        #logger$error(paste0('time: ', elapsed["elapsed"]))
        #logger$fatal(paste0('time: ', elapsed["elapsed"]))
        
        response$body <- jsonlite::toJSON(list(result = datasets,
                                               time   = elapsed["elapsed"]),
                                          auto_unbox = TRUE)
    },
    error = function(cond) {
        logger$error("ERROR")
        logger$error(jsonlite::toJSON(list(error=jsonlite::unbox(cond$message))))
        response$status_code <- 400
        response$body <- jsonlite::toJSON(list(error=jsonlite::unbox(cond$message)))
    })
    
    RestRserve::forward()
}



HttpGetLODScan <- function(request, response) {
  #' ---
  #' get: 
  #'   description: Get scan information.
  #'   responses:
  #'     200:
  #'       description: API response
  #'       content:
  #'         application/json:
  #'           schema:
  #'             type: string
  #'             example: 5
  #'     400:
  #'       description: API Error
  #'       content:
  #'         application/json:
  #'           schema:
  #'             type: string
  #'             example: Server error
  #' ---    
    result <- tryCatch({
        # start the clock
        ptm <- proc.time()

        dataset <- request$query[["dataset"]]
        id <- request$query[["id"]]
        interactive.covar <- request$query[["interactive.covar"]]
        cores <- nvlInteger(request$query[["cores"]], 5)
        expand <- request$query[["expand"]]
        
        if (tolower(nvl(interactive.covar, '')) == 'additive') {
            interactive.covar <- NULL
        }
        
        lod <- PerformLODScan(dataset           = dataset, 
                              id                = id,
                              interactive.covar = interactive.covar,
                              cores             = cores)
        
        if (!(toBoolean(expand))) {
            # by setting column names to NULL, the result will be a
            # 2 dimensional array
            colnames(lod) <- NULL
        }
        
        elapsed <- proc.time() - ptm
        #TrackTime(req, elapsed["elapsed"])
        
        response$body <- jsonlite::toJSON(list(result = lod,
                                               time   = elapsed["elapsed"]),
                                          auto_unbox = TRUE)
    },
    error = function(cond) {
        logger$error("ERROR")
        logger$error(jsonlite::toJSON(list(error=jsonlite::unbox(cond$message))))
        response$status_code <- 400
        response$body <- jsonlite::toJSON(list(error=jsonlite::unbox(cond$message)))
    })
    
    RestRserve::forward()
}


HttpGetLODScanSamples <- function(request, response) {
  #' ---
  #' get: 
  #'   description: Get scan by sample information.
  #'   responses:
  #'     200:
  #'       description: API response
  #'       content:
  #'         application/json:
  #'           schema:
  #'             type: string
  #'             example: 5
  #'     400:
  #'       description: API Error
  #'       content:
  #'         application/json:
  #'           schema:
  #'             type: string
  #'             example: Server error
  #' ---    
    result <- tryCatch({
        # start the clock
        ptm <- proc.time()

        dataset <- request$query[["dataset"]]
        id <- request$query[["id"]]
        interactive.covar <- request$query[["interactive.covar"]]
        chrom <- request$query[["chrom"]]
        regress.local <- request$query[["regress.local"]]
        cores <- nvlInteger(request$query[["cores"]], 5)
        expand <- request$query[["expand"]]

        if (tolower(nvl(interactive.covar, 'additive')) == 'additive') {
            stop("intCovar should not be additive or null")
        }
        
        lod <- PerformLODScanBySample(dataset           = dataset, 
                                      id                = id,
                                      chrom             = chrom,
                                      interactive.covar = interactive.covar,
                                      cores             = cores)
        
        if (!(toBoolean(expand))) {
            # by setting column names to NULL, the result will be a
            # 2 dimensional array
            for (element in names(lod)) {
                colnames(lod[[element]]) <- NULL
            }
        }
        
        elapsed <- proc.time() - ptm
        #TrackTime(req, elapsed["elapsed"])
        
        response$body <- jsonlite::toJSON(list(result = lod,
                                               time   = elapsed["elapsed"]),
                                          auto_unbox = TRUE)
    },
    error = function(cond) {
        logger$error("ERROR")
        logger$error(jsonlite::toJSON(list(error=jsonlite::unbox(cond$message))))
        response$status_code <- 400
        response$body <- jsonlite::toJSON(list(error=jsonlite::unbox(cond$message)))
    })
    
    RestRserve::forward()
}

HttpGetFounderCoefsScan <- function(request, response) {
    result <- tryCatch({
        # start the clock
        ptm <- proc.time()
        
        dataset <- request$query[["dataset"]]
        id <- request$query[["id"]]
        chrom <- request$query[["chrom"]]
        interactive.covar <- request$query[["interactive.covar"]]
        blup <- request$query[["blup"]]
        center <- request$query[["center"]]
        cores <- nvlInteger(request$query[["cores"]], 5)
        expand <- request$query[["expand"]]

        if (tolower(nvl(interactive.covar, '')) == 'additive') {
            intCovar <- NULL
        }
        
        effect <- PerformFounderCoefsScan(dataset           = dataset, 
                                          id                = id,
                                          chrom             = chrom, 
                                          interactive.covar = interactive.covar,
                                          blup              = blup, 
                                          center            = center, 
                                          cores             = cores)
        
        if (!(toBoolean(expand))) {
            # by setting column names to NULL, the result will be a
            # 2 dimensional array
            
            for (element in names(effect)) {
                colnames(effect[[element]]) <- NULL
            }
        }
        
        elapsed <- proc.time() - ptm
        #TrackTime(req, elapsed["elapsed"])
        
        response$body <- jsonlite::toJSON(list(result = effect,
                                               time   = elapsed["elapsed"]),
                                          auto_unbox = TRUE)
    },
    error = function(cond) {
        logger$error("ERROR")
        logger$error(jsonlite::toJSON(list(error=jsonlite::unbox(cond$message))))
        response$status_code <- 400
        response$body <- jsonlite::toJSON(list(error=jsonlite::unbox(cond$message)))
    })
    
    RestRserve::forward()
}

HttpGetExpression <- function(request, response) {
    result <- tryCatch({
        # start the clock
        ptm <- proc.time()
        
        dataset <- request$query[["dataset"]]
        id <- request$query[["id"]]

        expression <- GetExpression(dataset = dataset, 
                                    id      = id)
        
        elapsed <- proc.time() - ptm
        #TrackTime(req, elapsed["elapsed"])
        
        response$body <- jsonlite::toJSON(list(result = expression,
                                               time   = elapsed["elapsed"]),
                                          auto_unbox = TRUE)
    },
    error = function(cond) {
        logger$error("ERROR")
        logger$error(jsonlite::toJSON(list(error=jsonlite::unbox(cond$message))))
        response$status_code <- 400
        response$body <- jsonlite::toJSON(list(error=jsonlite::unbox(cond$message)))
    })
    
    RestRserve::forward()
}

HttpGetMediation <- function(request, response) {
    result <- tryCatch({
        # start the clock
        ptm <- proc.time()
        
        dataset <- request$query[["dataset"]]
        id <- request$query[["id"]]
        marker.id <- request$query[["marker.id"]]
        dataset.mediate <- request$query[["dataset.mediate"]]
        expand <- request$query[["expand"]]

        mediation <- PerformMediation(dataset         = dataset, 
                                      id              = id,
                                      marker.id       = marker.id,
                                      dataset.mediate = dataset.mediate)
        
        if (!(toBoolean(expand))) {
            # by setting column names to NULL, the result will be a
            # 2 dimensional array
            colnames(mediation) <- NULL
        }
        
        elapsed <- proc.time() - ptm
        #TrackTime(req, elapsed["elapsed"])
        
        response$body <- jsonlite::toJSON(list(result = mediation,
                                               time   = elapsed["elapsed"]),
                                          auto_unbox = TRUE)
    },
    error = function(cond) {
        logger$error("ERROR")
        logger$error(jsonlite::toJSON(list(error=jsonlite::unbox(cond$message))))
        response$status_code <- 400
        response$body <- jsonlite::toJSON(list(error=jsonlite::unbox(cond$message)))
    })
    
    RestRserve::forward()
}

HttpGetSNPAssocMapping <- function(request, response) {
    result <- tryCatch({
        # start the clock
        ptm <- proc.time()
        
        dataset <- request$query[["dataset"]]
        id <- request$query[["id"]]
        chrom <- request$query[["chrom"]]
        location <- request$query[["location"]]
        window.size <- nvlInteger(request$query[["window.size"]], 500000)
        cores <- nvlInteger(request$query[["cores"]], 5)
        expand <- request$query[["expand"]]

        snp.assoc <- PerformSNPAssocMapping(dataset     = dataset, 
                                            id          = id,
                                            chrom       = chrom,
                                            location    = location,
                                            window.size = window.size,
                                            cores       = cores)
        
        if (!(toBoolean(expand))) {
            # by setting column names to NULL, the result will be a
            # 2 dimensional array
            colnames(snp.assoc) <- NULL
        }
        
        elapsed <- proc.time() - ptm
        #TrackTime(req, elapsed["elapsed"])
        
        response$body <- jsonlite::toJSON(list(result = snp.assoc,
                                               time   = elapsed["elapsed"]),
                                          auto_unbox = TRUE)
    },
    error = function(cond) {
        logger$error("ERROR")
        logger$error(jsonlite::toJSON(list(error=jsonlite::unbox(cond$message))))
        response$status_code <- 400
        response$body <- jsonlite::toJSON(list(error=jsonlite::unbox(cond$message)))
    })
    
    RestRserve::forward()
}

HttpGetLODPeaks <- function(request, response) {
    result <- tryCatch({
        # start the clock
        ptm <- proc.time()
        
        dataset <- request$query[["dataset"]]
        interactive.covar <- request$query[["interactive.covar"]]
        expand <- request$query[["expand"]]

        lod.peaks <- GetLODPeaks(dataset, interactive.covar)
        
        if (!(toBoolean(expand))) {
            # by setting column names to NULL, the result will be a
            # 2 dimensional array
            colnames(lod.peaks) <- NULL
        }
        
        elapsed <- proc.time() - ptm
        #TrackTime(req, elapsed["elapsed"])

        response$body <- jsonlite::toJSON(list(result = lod.peaks,
                                               time   = elapsed["elapsed"]),
                                          auto_unbox = TRUE)

    },
    error = function(cond) {
        logger$error("ERROR")
        logger$error(jsonlite::toJSON(list(error=jsonlite::unbox(cond$message))))
        response$status_code <- 400
        response$body <- jsonlite::toJSON(list(error=jsonlite::unbox(cond$message)))
    })
    
    RestRserve::forward()
}


HttpGetLODPeaksAll <- function(request, response) {
    result <- tryCatch({
        # start the clock
        ptm <- proc.time()
        
        dataset <- request$query[["dataset"]]
        expand <- request$query[["expand"]]

        # get the LOD peaks for each covarint
        peaks <- GetLODPeaksAll(dataset)
        
        if (!(toBoolean(expand))) {
            for (n in names(peaks)) {
                colnames(peaks[[n]]) <- NULL
            }
        }
        
        elapsed <- proc.time() - ptm
        #TrackTime(req, elapsed["elapsed"])
        
        response$body <- jsonlite::toJSON(list(id     = dataset,
                                               result = peaks,
                                               time   = elapsed["elapsed"]),
                                          auto_unbox = TRUE)
    },
    error = function(cond) {
        logger$error("ERROR")
        logger$error(jsonlite::toJSON(list(error=jsonlite::unbox(cond$message))))
        response$status_code <- 400
        response$body <- jsonlite::toJSON(list(error=jsonlite::unbox(cond$message)))
    })
    
    RestRserve::forward()
}

HttpGetCorrelation <- function(request, response) {
    result <- tryCatch({
        # start the clock
        ptm <- proc.time()

        dataset <- request$query[["dataset"]]
        id <- request$query[["id"]]
        dataset.correlate <- request$query[["dataset.correlate"]]
        max.items <- request$query[["max.items"]]
        
        correlation <- PerformCorrelation(dataset           = dataset, 
                                          id                = id,
                                          dataset.correlate = dataset.correlate,
                                          max.items         = max.items)
        
        elapsed <- proc.time() - ptm
        #TrackTime(req, elapsed["elapsed"])
        
        response$body <- jsonlite::toJSON(list(result = correlation,
                                               time   = elapsed["elapsed"]),
                                          auto_unbox = TRUE)
    },
    error = function(cond) {
        logger$error("ERROR")
        logger$error(jsonlite::toJSON(list(error=jsonlite::unbox(cond$message))))
        response$status_code <- 400
        response$body <- jsonlite::toJSON(list(error=jsonlite::unbox(cond$message)))
    })
    
    RestRserve::forward()
}

HttpGetCorrelationPlot <- function(request, response) {
    result <- tryCatch({
        # start the clock
        ptm <- proc.time()
        
        dataset <- request$query[["dataset"]]
        id <- request$query[["id"]]
        dataset.correlate <- request$query[["dataset.correlate"]]
        id.correlate <- request$query[["id.correlate"]]

        correlation <- GetCorrelationPlotData(dataset           = dataset, 
                                              id                = id,
                                              dataset.correlate = dataset.correlate,
                                              id.correlate      = id.correlate)
        
        elapsed <- proc.time() - ptm
        #TrackTime(req, elapsed["elapsed"])
        
        response$body <- jsonlite::toJSON(list(result = correlation,
                                               time   = elapsed["elapsed"]),
                                          auto_unbox = TRUE)
    },
    error = function(cond) {
        logger$error("ERROR")
        logger$error(jsonlite::toJSON(list(error=jsonlite::unbox(cond$message))))
        response$status_code <- 400
        response$body <- jsonlite::toJSON(list(error=jsonlite::unbox(cond$message)))
    })
    
    RestRserve::forward()
}


RestRserveApp$add_get(path = "/datasets", FUN = HttpGetDataset, add_head = FALSE)
RestRserveApp$add_get(path = "/lodscan", FUN = HttpGetLODScan, add_head = FALSE)
RestRserveApp$add_get(path = "/lodscansamples", FUN = HttpGetLODScanSamples, add_head = FALSE)
RestRserveApp$add_get(path = "/foundercoefs", FUN = HttpGetFounderCoefsScan, add_head = FALSE)
RestRserveApp$add_get(path = "/expression", FUN = HttpGetExpression, add_head = FALSE)
RestRserveApp$add_get(path = "/mediate", FUN = HttpGetMediation, add_head = FALSE)
RestRserveApp$add_get(path = "/snpassoc", FUN = HttpGetSNPAssocMapping, add_head = FALSE)
RestRserveApp$add_get(path = "/lodpeaks", FUN = HttpGetLODPeaks, add_head = FALSE)
RestRserveApp$add_get(path = "/lodpeaksall", FUN = HttpGetLODPeaksAll, add_head = FALSE)
RestRserveApp$add_get(path = "/correlation", FUN = HttpGetCorrelation, add_head = FALSE)
RestRserveApp$add_get(path = "/correlationplot", FUN = HttpGetCorrelationPlot, add_head = FALSE)
RestRserveApp$add_static(path = "/qtl.RData", file_path = data.files[1])
#RestRserveApp$add_openapi(path = "/openapi.yaml", file_path = "openapi.yaml")
#RestRserveApp$add_swagger_ui(path = "/swagger", 
#                   path_openapi = "/openapi.yaml", 
#                   path_swagger_assets = "/__swagger__")


RestRserveApp$add_openapi()
RestRserveApp$add_swagger_ui("/")

logger$info("READY")

