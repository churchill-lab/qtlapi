get_ds <- function(id) {
    if (exists(id)) {
        return(get(id))
    } else {
        expanded <-paste0('dataset.', id) 
        if (exists(expanded)) {
            return(get(expanded))
        }
    }
    
    stop(sprintf("dataset '%s' not found", id))
}

get_d <- function(id, data = NULL) {
    ds <- get_ds(id)
    
    ret <- NULL
    
    if (is.matrix(ds$data) && (!is.null(data))) {
        stop(sprintf("specified data '%s' not found in '%s'", data, id))
    }
    
    if (is.matrix(ds$data)) {
        ret <- ds$data
    } else {
        if (!is.null(ds$data$rz)) {
            ret <- ds$data$rz
        } else if (!is.null(ds$data$norm)) {
            ret <- ds$data$norm
        } else if (!is.null(ds$data$raw)) {
            ret <- ds$data$raw
        } else if (!is.null(ds$data$log)) {
            ret <- ds$data$log
        } else if (!is.null(ds$data$transformed)) {
            ret <- ds$data$transformed
        }
    }
    ret
}

generate_shell_scripts <- function(file_name_data, outdir) {
    print(paste0('Loading: ', file_name_data))
    load(file_name_data, .GlobalEnv)
    fname <- gsub('[_.]', '', basename(file_name_data))

    for (ds_name in grep("^dataset.*", names(.GlobalEnv), value=TRUE)) {
        ds <- get_ds(ds_name)
        print(paste0('Generating script for: ', ds_name))

        cols <- ncol(get_d(ds_name))
        batch_size <- min(100, cols)

        file_name_base <- paste0(fname, '.', gsub('dataset.', '', ds_name))
        file_name_sh <- paste0(file_name_base, '.sh')
        file_name_out <- paste0(file_name_base, '.out.log')
        file_name_error <- paste0(file_name_base, '.error.log')
        
        cat('#!/bin/bash', file=file_name_sh, sep="\n")
        
        for (x in seq(1, cols, batch_size)) {
            # build the string
            command <- 'sbatch '
            command <- paste0(command, '-o ', file_name_out, ' ')
            command <- paste0(command, '-e ', file_name_error, ' ')
            command <- paste0(command, '--job-name=', file_name_base, '_', x, ' ')
            command <- paste0(command, '-v FINDPEAKSR=', G_FINDPEAKSR, ' ')
            command <- paste0(command, '-v DATAFILE=', file_name_data, ' ')
            command <- paste0(command, '-v QTLAPISOURCE=', G_QTLAPISOURCE, ' ')
            command <- paste0(command, '-v OUTPUTFILE=', file_name_base, '.PEAKS.out ')
            command <- paste0(command, '-v DATASET=', ds_name, ' ')
            command <- paste0(command, '-v START=', x, ' ')
            command <- paste0(command, '-v STEP=', batch_size, ' ')
            command <- paste0(command, G_FINDPEAKSSH)
            
            cat(command, file=file_name_sh, sep="\n", append = T)
        }
    }
}


#Rscript generatePeaksShellScript.R datafile.Rdata envFile
print('Generating shell scripts...')
args <- commandArgs(trailingOnly = TRUE)
outdir <- getwd()

if (length(args) != 2) {
    print('Rscript generatePeaksShellScript.R datafile.Rdata envFile')
}

print(paste0('Sourcing ', normalizePath(args[2]), ' ...'))
source(normalizePath(args[2]))

generate_shell_scripts(normalizePath(args[1]), outdir)



