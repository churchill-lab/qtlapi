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

qtlapiExists <- TRUE

if (!exists("RestRserveApp")) {
    # setup app and logging
    #logger <- RestRserve::Logger$new(level=TRACE, name="app", FUN=prettyJSON)
    logger <- RestRserve::Logger$new(level=TRACE, name="app")
    RestRserveApp <- RestRserve::RestRserveApplication$new()
    RestRserveApp$logger <- logger
    qtlapiExists <- FALSE
}


HttpRawAnnots <- function(request, response) {
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

        ds <- GetDataSet(dataset)

        if (tolower(ds$datatype) == "mrna") {
            columns <- colnames(ds$annot.mrna)
            data <- ds$annot.mrna
        } else if(tolower(ds$datatype) == "protein") {
            columns <- colnames(ds$annot.protein)
            data <- ds$annot.protein
        } else if(isPhenotype(ds)) {
            columns <- colnames(ds$annot.phenotype)
            data <- ds$annot.phenotype
        }

        elapsed <- proc.time() - ptm
        #TrackTime(req, elapsed["elapsed"])
        
        response$body <- jsonlite::toJSON(list(columns = columns,
                                               data    = data,
                                               time    = elapsed["elapsed"]),
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
HttpRawCovarInfo <- function(request, response) {
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

        ds <- GetDataSet(dataset)

        columns <- colnames(ds$covar.info)
        data <- ds$covar.info


        elapsed <- proc.time() - ptm
        #TrackTime(req, elapsed["elapsed"])
        
        response$body <- jsonlite::toJSON(list(columns = columns,
                                               data    = data,
                                               time    = elapsed["elapsed"]),
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

HttpRawAnnotSamples <- function(request, response) {
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

        ds <- GetDataSet(dataset)

        columns <- colnames(ds$annot.samples)
        data <- ds$annot.samples


        elapsed <- proc.time() - ptm
        #TrackTime(req, elapsed["elapsed"])
        
        response$body <- jsonlite::toJSON(list(columns = columns,
                                               data    = data,
                                               time    = elapsed["elapsed"]),
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
RestRserveApp$add_get(path = "/raw/annots", FUN = HttpRawAnnots, add_head = FALSE)
RestRserveApp$add_get(path = "/raw/covarinfo", FUN = HttpRawCovarInfo, add_head = FALSE)
RestRserveApp$add_get(path = "/raw/samples", FUN = HttpRawAnnotSamples, add_head = FALSE)



if (!qtlapiExists) {

    #RestRserveApp$add_openapi(path = "/openapi.yaml", file_path = "openapi.yaml")
    #RestRserveApp$add_swagger_ui(path = "/swagger", 
    #                   path_openapi = "/openapi.yaml", 
    #                   path_swagger_assets = "/__swagger__")


    RestRserveApp$add_openapi()
    RestRserveApp$add_swagger_ui("/")
}

#logger$info("READY")

