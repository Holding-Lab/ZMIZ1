#' Metadata for RNA-seq time course
#'
#'.
#' @docType data
#' @usage annotationTable
#'
#' @keywords datasets
#'
#' @examples
#' annotationTable
"annotationTable"

#' RNA-seq time course data
#'
#'.
#' @docType data
#' @usage rawcounts
#'
#' @keywords datasets
#'
#' @examples
#' rawcounts
"rawcounts"

#' First half of the TCGA Breast Cancer Expression Matrix
#'
#'.
#' @docType data
#' @usage expmat_brca_a
#'
#' @keywords datasets
#'
#' @examples
#' expmat_brca_a
"expmat_brca_a"

#' Second half of the TCGA Breast Cancer Expression Matrix
#'
#'.
#' @docType data
#' @usage expmat_brca_b
#'
#' @keywords datasets
#'
#' @examples
#' expmat_brca_b
"expmat_brca_b"

#' TCGA samples subtype metadata
#'
#'.
#' @docType data
#' @usage subtypes_brca
#'
#' @keywords datasets
#'
#' @examples
#' subtypes_brca
"subtypes_brca"

#' VIPER Activity Matrix for TCGA BRCA patient samples.
#'
#'.
#' @docType data
#' @usage vipermat_brca
#'
#' @keywords datasets
#'
#' @examples
#' vipermat_brca
"vipermat_brca"


#' METABRIC Breast Cancer Expression Matrix
#'
#'.
#' @docType data
#' @usage expmat_metabric
#'
#' @keywords datasets
#'
#' @examples
#' expmat_metabric
"expmat_metabric"


#' METABRIC samples subtype metadata
#'
#'.
#' @docType data
#' @usage subtypes_metabric
#'
#' @keywords datasets
#'
#' @examples
#' subtypes_metabric
"subtypes_metabric"

#' VIPER Activity Matrix for METABRIC patient samples.
#'
#'.
#' @docType data
#' @usage vipermat_metabric
#'
#' @keywords datasets
#'
#' @examples
#' vipermat_metabric
"vipermat_metabric"

#' MSigDB gene sets.
#'
#'.
#' @docType data
#' @usage msigdb
#'
#' @keywords datasets
#'
#' @examples
#' msigdb
"msigdb"

#' MSigDB gene sets.
#'
#'.
#' @docType data
#' @usage proteinData
#'
#' @keywords datasets
#'
#' @examples
#' proteinData
"proteinData"

## Make global
utils::globalVariables(c("annotationTable", "expmat_brca_a", "expmat_brca_b", "msigdb", "rawcounts"))



