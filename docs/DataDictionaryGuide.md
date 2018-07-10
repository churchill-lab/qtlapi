---
title: "Data Dictionary Definitions"
author: Gary Churchill
output: 
  github_document:
    md_extensions: +raw_html+markdown_in_html_blocks
---

# Data Dictionary Definitions
(version 2: June 15, 2017 – GAC)

This is a summary of the minimal sufficient data descriptions required to convert a data file into a generic R/qtl2 data environment for QTL-browser. Each column in the Dictionary is a data *descriptor* with definitions detailed below. Each row in the Dictionary corresponds to a column in the data file. The data file should have a single row of column headers that are referred to as *data_name*.  All columns in the data file should be defined in the Dictionary.

`Note. Blank fields are interpreted as FALSE for logical and NA for all other types.`

***data_name***: this is the name that the investigator of the data column. It must match the column header in the data file and is not case sensitive. 
> Required; must be unique.

`Note. Each column in the data file should have a header corresponding to *data_name*.  However, if a column header in the data file does not match any data_name in the dictionary they will be ignored. Matching must be exact but is not case sensitive. Conversely, any *data_name* in the dictionary that does not correspond to a column header in the data file will not be used in the R environment. These conventions allow to us subset data either by deleting columns from the data file with changing the dictionary or by deleting row entries in the data dictionary without changing the data file. Although matching dictionary data_names to columns in the data file is case insensitive, the case of *data_name* value will be preserved.`

***short_name***: It should be an abbreviation or short string that uniquely identifies the data column and is used to label axes on graphs or displayed in browser tools. If no *short_name* is provided, *data_name* will be substituted.
> Optional.

***R_name***: this is the name given to the data column in the R environment.  It should be short, begin with a letter, and no special characters except “.” or “_” which are useful when variable names contain data, e.g.  BW_6  BW_12  and BW_18 are bodyweights of the same mouse at different ages.  The R_name can be parsed in R, e.g. with tidyr::separate(). If “NA” then omit must be set to TRUE and the data column will be skipped. 
* Required; must be unique.

***description***: this is a description field for the data column.  It should be a clear and complete but concise description of the data.
> Required. 

***units***: units of measurement could be part of description, but when units are needed, e.g. for data conversions, they should be explicit here.  Use standard abbreviations, e.g. “cm” or “ug”.
> Optional. 

***is_id***: a logical variable that defines which data column is used to identify the animal or experimental unit. It must have exactly one element TRUE.
> Required.

***category***: this is a grouping variable for data columns. In large data sets, variables can be organized into meaningful groups, e.g., “Hematology” and “Body Composition”.  This is useful for organizing pull-down menus as well understanding data structure.
> Optional.

`Note.  Use the reserved category name “ID” for variables that identify the animal/unit. These will not be added to pull-down menu and will be retained unless omit is TRUE.`

***R_category***: a short R-compatible variable name for each category.
> Optional.

***is_numeric***: an indicator, TRUE if the data are numeric values.
> Required.

***is_date***: an indicator, TRUE if the data are dates.
> Optional.

***is_factor***: an indicator, TRUE if the data are discrete factor levels.
> Required.  

`Note. If is_numeric, is_date, and is_factor are all FALSE, the data column will be treated as a text field. Is_factor and is_numeric can both be TRUE but the factor levels should be given as integers, e.g., for mice samples at 3 ages use 6:12:18. Data that are listed as both numeric and factor, will be treated as numeric if used as phenotpyes and will be treated as a factor if used as a covariate.`

***factor_levels***: a string of factor level names with “:” as separator The factor level names must match the entries in the data column, e.g., Female:Male or F:M. The order of factor levels will match the order given here. If missing, factor levels are computed and ordered alphabetically using strings in data column.
> Optional.

***is_covar***: an indicator, TRUE if data column is to be used as a covariate.
> Required.

***is_pheno***: an indicator, TRUE if data column is to be used as a phenotype.
> Required.

`Note. Data columns that are covariates or phenotypes are placed in different data structures in the R/qtl2 environment.  If both is_covar and is_pheno are true the data column will be represented in both data structures. See above note about factors and numeric variables.`

***is_derived***: an indicator that is TRUE if the data column is calculated from other data columns. Maybe useful to indicate if a phenotype is computed as a ratio or linear combination of other phenotypes.
> Optional.

***omit***: an indicator that is TRUE if we do not want to include this data column in the R environment. Useful for any fields that we want to retain in the data but are not otherwise used.
> Optional.

***use_covar***: is a colon-separated “:”  list of covariates that will used as default in genome scans. It is only used if is_pheno=TRUE; it must be NA if is_pheno=FALSE. Names in the colon-separated list must match the R_data_name to be used as a covariate.
> Optional.

