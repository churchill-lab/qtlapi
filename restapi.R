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


HttpDatasetGet <- function(request, response) {
  #' ---
  #' description: Get dataset information.
  #' responses:
  #'   200:
  #'     description: API response
  #'     content:
  #'       application/json:
  #'         schema:
  #'           type: string
  #'           example: 5
  #'   400:
  #'     description: API Error
  #'     content:
  #'       application/json:
  #'         schema:
  #'           type: string
  #'           example: Server error
  #' ---    

    result <- tryCatch({
        # start the clock
        ptm <- proc.time()
        
        datasets <- GetDatasetInfo()
        
        elapsed <- proc.time() - ptm
        #TrackTime(req, elapsed["elapsed"])
        
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



HttpLODScanGet <- function(request, response) {
    result <- tryCatch({
        # start the clock
        ptm <- proc.time()

        dataset <- request$query[["dataset"]]
        id <- request$query[["id"]]
        intCovar <- request$query[["intCovar"]]
        nCores <- nvlInteger(request$query[["nCores"]], 5)
        expand <- request$query[["expand"]]
        
        print(paste0('intCovar', intCovar))

        if (tolower(nvl(intCovar, '')) == 'additive') {
            intCovar <- NULL
        }
        
        print(paste0('intCovar', intCovar))

        lod <- GetLODScan(dataset      = dataset, 
                          id           = id,
                          intCovar     = intCovar,
                          nCores       = nCores)
        
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


HttpLODScanSamplesGet <- function(request, response) {
    result <- tryCatch({
        # start the clock
        ptm <- proc.time()

        dataset <- request$query[["dataset"]]
        id <- request$query[["id"]]
        intCovar <- request$query[["intCovar"]]
        chrom <- request$query[["chrom"]]
        reegressLocal <- request$query[["regressLocal"]]
        nCores <- nvlInteger(request$query[["nCores"]], 5)
        expand <- request$query[["expand"]]

        if (tolower(nvl(intCovar, 'additive')) == 'additive') {
            stop("intCovar should not be additive or null")
        }
        
        lod <- GetLODScanBySample(dataset      = dataset, 
                                  id           = id,
                                  intCovar     = intCovar,
                                  chrom        = chrom,
                                  nCores       = nCores)
        
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

HttpFoundercoefsGet <- function(request, response) {
    result <- tryCatch({
        # start the clock
        ptm <- proc.time()
        
        dataset <- request$query[["dataset"]]
        id <- request$query[["id"]]
        chrom <- request$query[["chrom"]]
        intCovar <- request$query[["intCovar"]]
        blup <- request$query[["blup"]]
        center <- request$query[["center"]]
        nCores <- nvlInteger(request$query[["nCores"]], 5)
        expand <- request$query[["expand"]]

        if (tolower(nvl(intCovar, '')) == 'additive') {
            intCovar <- NULL
        }
        
        effect <- GetFoundercoefs(dataset      = dataset, 
                                  id           = id,
                                  chrom        = chrom, 
                                  intCovar     = intCovar,
                                  blup         = blup, 
                                  center       = center, 
                                  nCores       = nCores)
        
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

HttpExpressionGet <- function(request, response) {
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

HttpMediateGet <- function(request, response) {
    result <- tryCatch({
        # start the clock
        ptm <- proc.time()
        
        dataset <- request$query[["dataset"]]
        id <- request$query[["id"]]
        mid <- request$query[["mid"]]
        datasetMediate <- request$query[["datasetMediate"]]
        expand <- request$query[["expand"]]


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

HttpSnpAssocMappingGet <- function(request, response) {
    result <- tryCatch({
        # start the clock
        ptm <- proc.time()
        
        dataset <- request$query[["dataset"]]
        id <- request$query[["id"]]
        chrom <- request$query[["chrom"]]
        location <- request$query[["location"]]
        windowSize <- nvlInteger(request$query[["windowSize"]], 500000)
        nCores <- nvlInteger(request$query[["nCores"]], 5)
        expand <- request$query[["expand"]]

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
        #TrackTime(req, elapsed["elapsed"])
        
        response$body <- jsonlite::toJSON(list(result = snpAssoc,
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

HttpLODPeaksGet <- function(request, response) {
    result <- tryCatch({
        # start the clock
        ptm <- proc.time()
        
        dataset <- request$query[["dataset"]]
        intCovar <- request$query[["intCovar"]]
        expand <- request$query[["expand"]]

        lodPeaks <- GetLODPeaks(dataset, intCovar)
        
        if (!(toBoolean(expand))) {
            # by setting column names to NULL, the result will be a
            # 2 dimensional array
            colnames(lodPeaks) <- NULL
        }
        
        elapsed <- proc.time() - ptm
        #TrackTime(req, elapsed["elapsed"])

        response$body <- jsonlite::toJSON(list(result = lodPeaks,
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


HttpLODPeaksAllGet <- function(request, response) {
    result <- tryCatch({
        # start the clock
        ptm <- proc.time()
        
        dataset <- request$query[["dataset"]]
        expand <- request$query[["expand"]]

        # get the LOD peaks for each covarint
        additivePeaks <- GetLODPeaks(dataset)
        
        if (!(toBoolean(expand))) {
            colnames(additivePeaks) <- NULL
        }
        
        peaks <- list(additive = additivePeaks)
        intCovars <- GetIntCovar(dataset)
        
        for (ic in intCovars) {
            icPeaks <- GetLODPeaks(dataset, ic)
            
            if (!(toBoolean(expand))) {
                colnames(icPeaks) <- NULL
            }
            
            peaks <- c(list(icPeaks), peaks)
            names(peaks)[1] = ic
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

HttpCorrelationGet <- function(request, response) {
    result <- tryCatch({
        # start the clock
        ptm <- proc.time()

        dataset <- request$query[["dataset"]]
        id <- request$query[["id"]]
        datasetCorrelate <- request$query[["datasetCorrelate"]]
        maxItems <- request$query[["maxItems"]]
        
        correlation <- GetCorrelation(dataset          = dataset, 
                                      id               = id,
                                      datasetCorrelate = datasetCorrelate,
                                      maxItems         = maxItems)
        
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

HttpCorrelationPlotDataGet <- function(request, response) {
    result <- tryCatch({
        # start the clock
        ptm <- proc.time()
        
        dataset <- request$query[["dataset"]]
        id <- request$query[["id"]]
        datasetCorrelate <- request$query[["datasetCorrelate"]]
        idCorrelate <- request$query[["idCorrelate"]]

        correlation <- GetCorrelationPlotData(dataset          = dataset, 
                                              id               = id,
                                              datasetCorrelate = datasetCorrelate,
                                              idCorrelate      = idCorrelate)
        
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


RestRserveApp$add_get(path = "/datasets", FUN = HttpDatasetGet)
RestRserveApp$add_get(path = "/lodscan", FUN = HttpLODScanGet)
RestRserveApp$add_get(path = "/lodscansamples", FUN = HttpLODScanSamplesGet)
RestRserveApp$add_get(path = "/foundercoefs", FUN = HttpFoundercoefsGet)
RestRserveApp$add_get(path = "/expression", FUN = HttpExpressionGet)
RestRserveApp$add_get(path = "/mediate", FUN = HttpMediateGet)
RestRserveApp$add_get(path = "/snpassoc", FUN = HttpSnpAssocMappingGet)
RestRserveApp$add_get(path = "/lodpeaks", FUN = HttpLODPeaksGet)
RestRserveApp$add_get(path = "/lodpeaksall", FUN = HttpLODPeaksAllGet)
RestRserveApp$add_get(path = "/correlation", FUN = HttpCorrelationGet)
RestRserveApp$add_get(path = "/correlationPlotData", FUN = HttpCorrelationPlotDataGet)





#RestRserveApp$add_openapi(path = "/openapi.yaml", file_path = "openapi.yaml")
#RestRserveApp$add_swagger_ui(path = "/swagger", 
#                   path_openapi = "/openapi.yaml", 
#                   path_swagger_assets = "/__swagger__")


RestRserveApp$add_openapi()
RestRserveApp$add_swagger_ui("/")

logger$info("READY")

