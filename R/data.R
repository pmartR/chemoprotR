# label_setup --------------------------------------------------------

#' Labeled HTP Chemoproteomics Setup File
#'
#' These data are the original formatting of the high-throughput (HTP) chemoproteomics
#' file after loading in the data from excel. 
#'
#' @format A List object of Length 7
#' \describe{
#'   \item{T_alias}{a data frame with dimensions \eqn{22 \times 3} pertaining to alias data}
#'   \item{T_data}{a data frame with dimensions \eqn{2890 \times 22} pertaining to peptide-level data}
#'   \item{T_Data_Package_Analysis_Jobs}{a data frame with dimensions \eqn{4 \times 25} pertaining to T data package analysis jobs data
#'   where each row corresponds to a different job}
#'   \item{T_Reporter_Ions_Typed}{a data frame with dimensions \eqn{28188 \times 19} data pertaining to the reported intensity for each of the 
#'   ionizations}
#'   \item{Mage}{a data frame with dimensions \eqn{39} pertaining to Mage information}
#'   \item{f_data}{a data frame with dimensions \eqn{22 \times 6} pertaining to feature data.
#'   Each row corresponds to a sample with a column giving the unique sample identifiers found in e_data
#'   column names, the injection order, the group information, and the batch information.}
#'   \item{Protein_collection_data}{a data frame with dimensions \eqn{5564 \times 8} pertaining to protein collection data}
#' }
#' @rdname label_setup
#' @name label_setup
NULL

# label_pmart --------------------------------------------------------

#' Cleaned Labeled HTP pmartR-friendly Object
#'
#' These data originates from the label_setup dataset and pertain to values from Job 2160486. However, the data
#' has been transformed and filtered such that the data is compatible with pmartR
#'
#' @format An  isobaricpepData object (see
#'   \code{\link[pmartR]{as.isobaricpepData}} for details)
#' \describe{
#'   \item{e_data}{a \eqn{12665 \times 12} data frame of expression data.
#'   Each row corresponds to expression data for each isobarically labeled peptide}
#'   \item{f_data}{a data frame with \eqn{1 \times 6} data frame of feature data.
#'   Each row corresponds to a sample with a column giving the unique sample identifiers found in e_data
#'   column names, the grouping, job number, replicate, plex, and ionization}
#'   \item{e_meta}{a data frame with \eqn{12665 \times 14} data frame of biomolecule information.
#'   Each row corresponds to data pertaining to biomolecule information for each isobarically labeled peptide}
#' }
#' @rdname label_pmart
#' @name label_pmart
NULL

# labelfree_setup --------------------------------------------------------

#' Label-Free HTP Chemoproteomics Setup File
#'
#' These data are the original formatting of the high-throughput (HTP) chemoproteomics
#' file after loading in the data from excel. 
#'
#' @format A List object of Length 4
#' \describe{
#'   \item{Metadata}{a data frame with dimensions \eqn{6 \times 27} pertaining to metadata where each row
#'   corresponds to a different job}
#'   \item{f_data}{a data frame with \eqn{6 \times 4} data frame of feature data. Each row corresponds to a sample with a column giving the unique sample identifiers found in e_data
#'   column names, the injection order, the group information, and the batch information.}
#'   \item{Mage}{a data frame with dimensions \eqn{132732 \times 39} pertaining to mage information, and stores
#'   the Total Ion Intensity values}
#'   \item{Protein collection list}{a data frame with dimensions \eqn{243824 \times 8} pertaining to protein collection data}
#' }
#' @rdname labelfree_setup
#' @name labelfree_setup
NULL

# labelfree_pmart --------------------------------------------------------

#' Cleaned Labeled HTP pmartR-friendly Object
#'
#' These data originates from the labelfree_setup dataset. However, the data
#' has been transformed and filtered such that the data is compatible with pmartR
#'
#' @format A pepData object (see
#'   \code{\link[pmartR]{as.pepData}} for details)
#' \describe{
#'   \item{e_data}{a \eqn{294655 \times 10} data frame of expression data.
#'   Each row corresponds to expression data for each label-free peptide}
#'   \item{f_data}{a data frame with \eqn{9 \times 5} data frame of feature data.
#'   Each row corresponds to a sample with a column giving the unique sample identifiers found in e_data
#'   column names, the grouping, job, replicate, and job number}
#'   \item{e_meta}{a \eqn{294655 \times 16} data frame of expression data.
#'   Each row corresponds to data pertaining to biomolecule information for each label-free peptide}
#' }
#' @rdname labelfree_pmart
#' @name labelfree_pmart
NULL