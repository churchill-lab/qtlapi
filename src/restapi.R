library(RestRserve)
library(jsonlite)

pretty_JSON <- function(timestamp, level, logger_name, pid, message) {
    x = to_json(
        list(
            timestamp = format(timestamp, "%Y-%m-%d %H:%M:%OS6"),
            level = as.character(level),
            name = as.character(logger_name),
            pid = as.integer(pid),
            message = message
        )
    )
    cat(prettify(x), file = "", append = TRUE, sep = "\n")
}

simple_csv <- function(timestamp, level, logger_name, pid, message) {
    x = sprintf("[%s][%s][%s][%d][\"%s\"]",
                format(timestamp, "%Y-%m-%d %H:%M:%OS6"),
                as.character(level),
                as.character(logger_name),
                as.integer(pid),
                message)
    cat(x, file = "", append = TRUE, sep = "\n")
}

# setup app and logging
logger <- RestRserve::Logger$new(level=TRACE, name="app", FUN=simple_csv)
#logger <- RestRserve::Logger$new(level=TRACE, name="app")
RestRserveApp <- RestRserve::RestRserveApplication$new()
RestRserveApp$logger <- logger


http_get_dataset <- function(request, response) {
    result <- tryCatch({
        # start the clock
        ptm <- proc.time()
    
        datasets <- get_dataset_info()
    
        elapsed <- proc.time() - ptm
        
        logger$info(paste0('http_get_dataset - time: ', elapsed["elapsed"]))
        
        #logger$trace(paste0('time: ', elapsed["elapsed"]))
        #logger$debug(paste0('time: ', elapsed["elapsed"]))
        #logger$info(paste0('time: ', elapsed["elapsed"]))
        #logger$warning(paste0('time: ', elapsed["elapsed"]))
        #logger$error(paste0('time: ', elapsed["elapsed"]))
        #logger$fatal(paste0('time: ', elapsed["elapsed"]))
    
        response$body <- toJSON(list(result = datasets,
                                     time   = elapsed["elapsed"]),
                                auto_unbox = TRUE)
    },
    error = function(cond) {
        logger$error(sprintf("ERROR: http_get_dataset - %s", cond$message))
        response$status_code <- 400
        response$body <- toJSON(list(method = 'http_get_dataset',
                                     error  = cond$message),
                                auto_unbox = TRUE)
    })
  
    RestRserve::forward()
}


http_get_lodscan <- function(request, response) {
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
        intcovar <- request$query[["intcovar"]]
        cores <- nvl_int(request$query[["cores"]], 5)
        expand <- nvl(request$query[["expand"]], 'FALSE')
        
        if (tolower(nvl(intcovar, '')) %in% c('', 'additive')) {
            intcovar <- NULL
        }
        
        lod <- get_lod_scan(dataset  = dataset, 
                            id       = id,
                            intcovar = intcovar,
                            cores    = cores)
        
        if (!(to_boolean(expand))) {
            # by setting column names to NULL, the result will be a
            # 2 dimensional array
            colnames(lod) <- NULL
        }
        
        elapsed <- proc.time() - ptm
  
        logger$info(paste0('http_get_lodscan - time: ', elapsed["elapsed"]))
        
        response$body <- toJSON(list(result = lod,
                                     time   = elapsed["elapsed"]),
                                auto_unbox = TRUE)
    },
    error = function(cond) {
        logger$error(sprintf("ERROR: http_get_lodscan - %s", cond$message))
        response$status_code <- 400
        response$body <- toJSON(list(method = 'http_get_lodscan',
                                     error  = cond$message),
                                auto_unbox = TRUE)
    })
      
    RestRserve::forward()
}


http_get_lodscan_samples <- function(request, response) {
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
        intcovar <- request$query[["intcovar"]]
        chrom <- request$query[["chrom"]]
        #regress.local <- request$query[["regress.local"]]
        cores <- nvl_int(request$query[["cores"]], 5)
        expand <- nvl(request$query[["expand"]], 'FALSE')
        
        if (tolower(nvl(intcovar, '')) %in% c('', 'additive')) {
            stop("intcovar should not be additive or null")
        }
        
        lod <- get_lod_scan_by_sample(dataset  = dataset, 
                                      id       = id,
                                      chrom    = chrom,
                                      intcovar = intcovar,
                                      cores    = cores)
        
        if (!(to_boolean(expand))) {
            # by setting column names to NULL, the result will be a
            # 2 dimensional array
            for (element in names(lod)) {
                colnames(lod[[element]]) <- NULL
            }
        }
        
        elapsed <- proc.time() - ptm

        logger$info(paste0('http_get_lodscan_samples - time: ', elapsed["elapsed"]))
        
        response$body <- toJSON(list(result = lod,
                                     time   = elapsed["elapsed"]),
                                auto_unbox = TRUE)
    },
    error = function(cond) {
        logger$error(sprintf("ERROR: http_get_lodscan_samples - %s", cond$message))
        response$status_code <- 400
        response$body <- toJSON(list(method = 'http_get_lodscan_samples',
                                     error  = cond$message),
                                auto_unbox = TRUE)
    })
    
    RestRserve::forward()
}

http_get_foundercoefficients <- function(request, response) {
    result <- tryCatch({
        # start the clock
        ptm <- proc.time()
        
        dataset <- request$query[["dataset"]]
        id <- request$query[["id"]]
        chrom <- request$query[["chrom"]]
        intcovar <- request$query[["intcovar"]]
        blup <- nvl(request$query[["blup"]], 'FALSE')
        center <- nvl(request$query[["center"]], "TRUE")
        cores <- nvl_int(request$query[["cores"]], 5)
        expand <- nvl(request$query[["expand"]], 'FALSE')
        
        if (tolower(nvl(intcovar, '')) %in% c('', 'additive')) {
            intcovar <- NULL
        }
        
        effect <- get_founder_coefficients(dataset  = dataset, 
                                           id       = id,
                                           chrom    = chrom, 
                                           intcovar = intcovar,
                                           blup     = blup, 
                                           center   = center, 
                                           cores    = cores)
        
        if (!(to_boolean(expand))) {
            # by setting column names to NULL, the result will be a
            # 2 dimensional array
            
            for (element in names(effect)) {
                colnames(effect[[element]]) <- NULL
            }
        }
        
        elapsed <- proc.time() - ptm
  
        logger$info(paste0('http_get_foundercoefficients - time: ', elapsed["elapsed"]))
        
        response$body <- toJSON(list(result = effect,
                                     time   = elapsed["elapsed"]),
                                auto_unbox = TRUE)
    },
    error = function(cond) {
        logger$error(sprintf("ERROR: http_get_foundercoefficients - %s", cond$message))
        response$status_code <- 400
        response$body <- toJSON(list(method = 'http_get_foundercoefficients',
                                     error  = cond$message),
                                auto_unbox = TRUE)
    })
    
    RestRserve::forward()
}

http_get_expression <- function(request, response) {
    result <- tryCatch({
      # start the clock
      ptm <- proc.time()
      
      dataset <- request$query[["dataset"]]
      id <- request$query[["id"]]
      
      expression <- get_expression(dataset = dataset, 
                                   id      = id)
      
      # eliminate the _row column down line for JSON
      rownames(expression$data) <- NULL

      elapsed <- proc.time() - ptm
      
      logger$info(paste0('http_get_expression - time: ', elapsed["elapsed"]))
      
      response$body <- toJSON(list(result = expression,
                                   time   = elapsed["elapsed"]),
                              auto_unbox = TRUE)
    },
    error = function(cond) {
        logger$error(sprintf("ERROR: http_get_expression - %s", cond$message))
        response$status_code <- 400
        response$body <- toJSON(list(method = 'http_get_expression',
                                     error  = cond$message),
                                auto_unbox = TRUE)
    })
    
    RestRserve::forward()
}

http_get_mediation <- function(request, response) {
    result <- tryCatch({
        # start the clock
        ptm <- proc.time()
        
        dataset <- request$query[["dataset"]]
        id <- request$query[["id"]]
        marker_id <- request$query[["marker_id"]]
        dataset_mediate <- request$query[["dataset_mediate"]]
        intcovar <- request$query[["intcovar"]]
        expand <- nvl(request$query[["expand"]], 'FALSE')

        if (tolower(nvl(intcovar, '')) %in% c('', 'none')) {
            intcovar <- NULL
        }

        mediation <- get_mediation(dataset         = dataset, 
                                   id              = id,
                                   marker_id       = marker_id,
                                   dataset_mediate = dataset_mediate,
                                   intcovar        = intcovar)
        
        if (!(to_boolean(expand))) {
            # by setting column names to NULL, the result will be a
            # 2 dimensional array
            colnames(mediation) <- NULL
        }

        elapsed <- proc.time() - ptm

        logger$info(paste0('http_get_mediation - time: ', elapsed["elapsed"]))

        response$body <- toJSON(list(result = mediation,
                                     time   = elapsed["elapsed"]),
                                auto_unbox = TRUE)
    },
    error = function(cond) {
        logger$error(sprintf("ERROR: http_get_mediation - %s", cond$message))
        response$status_code <- 400
        response$body <- toJSON(list(method = 'http_get_mediation',
                                     error  = cond$message),
                                auto_unbox = TRUE)
    })
    
    RestRserve::forward()
}

http_get_snp_assoc_mapping <- function(request, response) {
    result <- tryCatch({
        # start the clock
        ptm <- proc.time()
        
        dataset <- request$query[["dataset"]]
        id <- request$query[["id"]]
        chrom <- request$query[["chrom"]]
        location <- request$query[["location"]]
        window_size <- nvl_int(request$query[["window_size"]], 500000)
        intcovar <- request$query[["intcovar"]]
        cores <- nvl_int(request$query[["cores"]], 5)
        expand <- nvl(request$query[["expand"]], 'FALSE')
        
        if (tolower(nvl(intcovar, '')) %in% c('', 'additive')) {
          intcovar <- NULL
        }
        
        snp_assoc <- get_snp_assoc_mapping(dataset      = dataset, 
                                            id          = id,
                                            chrom       = chrom,
                                            location    = location,
                                            window_size = window_size,
                                            intcovar    = intcovar,
                                            cores       = cores)
        
        if (!(to_boolean(expand))) {
            # by setting column names to NULL, the result will be a
            # 2 dimensional array
            colnames(snp_assoc) <- NULL
        }
        
        elapsed <- proc.time() - ptm

        logger$info(paste0('http_get_snp_assoc_mapping - time: ', elapsed["elapsed"]))
        
        response$body <- toJSON(list(result = snp_assoc,
                                     time   = elapsed["elapsed"]),
                                auto_unbox = TRUE)
    },
    error = function(cond) {
        logger$error(sprintf("ERROR: http_get_snp_assoc_mapping - %s", cond$message))
        response$status_code <- 400
        response$body <- toJSON(list(method = 'http_get_snp_assoc_mapping',
                                     error  = cond$message),
                                auto_unbox = TRUE)
    })
    
    RestRserve::forward()
}

http_get_lod_peaks <- function(request, response) {
    result <- tryCatch({
        # start the clock
        ptm <- proc.time()
        
        dataset <- request$query[["dataset"]]
        intcovar <- request$query[["intcovar"]]
        expand <- nvl(request$query[["expand"]], 'FALSE')
        
        if (tolower(nvl(intcovar, '')) %in% c('', 'additive')) {
            intcovar <- NULL
        }
        
        lod_peaks <- get_lod_peaks(dataset, intcovar)
        
        if (!(to_boolean(expand))) {
            # by setting column names to NULL, the result will be a
            # 2 dimensional array
            colnames(lod_peaks) <- NULL
        }
        
        elapsed <- proc.time() - ptm
        
        logger$info(paste0('http_get_lod_peaks - time: ', elapsed["elapsed"]))
        
        response$body <- toJSON(list(result = lod_peaks,
                                     time   = elapsed["elapsed"]),
                                auto_unbox = TRUE)
    },
    error = function(cond) {
        logger$error(sprintf("ERROR: http_get_lod_peaks - %s", cond$message))
        response$status_code <- 400
        response$body <- toJSON(list(method = 'http_get_lod_peaks',
                                     error  = cond$message),
                                auto_unbox = TRUE)
    })
    
    RestRserve::forward()
}


http_get_lod_peaks_all <- function(request, response) {
    result <- tryCatch({
      # start the clock
      ptm <- proc.time()
      
      dataset <- request$query[["dataset"]]
      expand <- nvl(request$query[["expand"]], 'FALSE')
      
      # get the LOD peaks for each covarint
      peaks <- get_lod_peaks_all(dataset)
      
      if (!(to_boolean(expand))) {
          for (n in names(peaks)) {
              colnames(peaks[[n]]) <- NULL
          }
      }
      
      elapsed <- proc.time() - ptm
      
      logger$info(paste0('http_get_lod_peaks_all - time: ', elapsed["elapsed"]))
      
      response$body <- toJSON(list(id     = dataset,
                                   result = peaks,
                                   time   = elapsed["elapsed"]),
                              auto_unbox = TRUE)
    },
    error = function(cond) {
        logger$error(sprintf("ERROR: http_get_lod_peaks_all - %s", cond$message))
        response$status_code <- 400
        response$body <- toJSON(list(method = 'http_get_lod_peaks_all',
                                     error  = cond$message),
                                auto_unbox = TRUE)
    })
    
    RestRserve::forward()
}

http_get_correlation <- function(request, response) {
    result <- tryCatch({
        # start the clock
        ptm <- proc.time()
        
        dataset <- request$query[["dataset"]]
        id <- request$query[["id"]]
        dataset_correlate <- request$query[["dataset_correlate"]]
        intcovar <- request$query[["intcovar"]]
        max_items <- nvl_int(request$query[["max_items"]], 10000)
        
        if (tolower(nvl(intcovar, '')) %in% c('', 'none')) {
            intcovar <- NULL
        }
        
        correlation <- get_correlation(dataset           = dataset, 
                                       id                = id,
                                       dataset_correlate = dataset_correlate,
                                       intcovar          = intcovar)
        
        correlation <- correlation[1:nvl_int(max.items, length(correlation))]        

        elapsed <- proc.time() - ptm
        
        logger$info(paste0('http_get_correlation - time: ', elapsed["elapsed"]))
        
        response$body <- toJSON(list(result = correlation,
                                     time   = elapsed["elapsed"]),
                                auto_unbox = TRUE)
    },
    error = function(cond) {
        logger$error(sprintf("ERROR: http_get_correlation - %s", cond$message))
        response$status_code <- 400
        response$body <- toJSON(list(method = 'http_get_correlation',
                                     error  = cond$message),
                                auto_unbox = TRUE)
    })
    
    RestRserve::forward()
}

http_get_correlation_plot_data <- function(request, response) {
    result <- tryCatch({
        # start the clock
        ptm <- proc.time()
        
        dataset <- request$query[["dataset"]]
        id <- request$query[["id"]]
        dataset_correlate <- request$query[["dataset_correlate"]]
        id_correlate <- request$query[["id_correlate"]]
        intcovar <- request$query[["intcovar"]]
        
        if (tolower(nvl(intcovar, '')) %in% c('', 'none')) {
            intcovar <- NULL
        }
        
        correlation <- 
          get_correlation_plot_data(dataset           = dataset,
                                    id                = id,
                                    dataset_correlate = dataset_correlate,
                                    intcovar          = intcovar,
                                    id_correlate      = id_correlate)
        
        elapsed <- proc.time() - ptm
        
        logger$info(paste0('http_get_correlation_plot_data - time: ', elapsed["elapsed"]))
        
        response$body <- toJSON(list(result = correlation,
                                     time   = elapsed["elapsed"]),
                                auto_unbox = TRUE)
    },
    error = function(cond) {
        logger$error(sprintf("ERROR: http_get_correlation_plot_data - %s", cond$message))
        response$status_code <- 400
        response$body <- toJSON(list(method = 'http_get_correlation_plot_data',
                                     error  = cond$message),
                                auto_unbox = TRUE)
    })
    
    RestRserve::forward()
}

RestRserveApp$add_get(path = "/datasets", FUN = http_get_dataset, add_head = FALSE)
RestRserveApp$add_get(path = "/lodscan", FUN = http_get_lodscan, add_head = FALSE)
RestRserveApp$add_get(path = "/lodscansamples", FUN = http_get_lodscan_samples, add_head = FALSE)
RestRserveApp$add_get(path = "/foundercoefs", FUN = http_get_foundercoefficients, add_head = FALSE)
RestRserveApp$add_get(path = "/expression", FUN = http_get_expression, add_head = FALSE)
RestRserveApp$add_get(path = "/mediate", FUN = http_get_mediation, add_head = FALSE)
RestRserveApp$add_get(path = "/snpassoc", FUN = http_get_snp_assoc_mapping, add_head = FALSE)
RestRserveApp$add_get(path = "/lodpeaks", FUN = http_get_lod_peaks, add_head = FALSE)
RestRserveApp$add_get(path = "/lodpeaksall", FUN = http_get_lod_peaks_all, add_head = FALSE)
RestRserveApp$add_get(path = "/correlation", FUN = http_get_correlation, add_head = FALSE)
RestRserveApp$add_get(path = "/correlationplot", FUN = http_get_correlation_plot_data, add_head = FALSE)

#RestRserveApp$add_openapi(path = "/openapi.yaml", file_path = "openapi.yaml")
#RestRserveApp$add_swagger_ui(path = "/swagger", 
#                   path_openapi = "/openapi.yaml", 
#                   path_swagger_assets = "/__swagger__")


RestRserveApp$add_openapi()
RestRserveApp$add_swagger_ui("/")

logger$info("READY")

