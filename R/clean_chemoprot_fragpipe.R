#' Clean Labeled Chemoproteomics Data for MSGF/MASIC Processing
#' 
#' This function takes labeled chemoproteomics data from an excel format 
#' and returns an S3 object compatible with another R package, pmartR, that 
#' can then be ran in a markdown file for further filtering, normalization,
#' and statistical analyses. This function is used for data processed using
#' MSGF/MASIC.
#' 
#' @param m List containing information from each of the different TMT excel sheets
#' @param study_df 1 x n data.frame specifying the names of the n different excel sheets
#' in the TMT excel file and what function to which they correspond. At the time the function
#' was written, n = 9
#' @param mage_cols 1 x n data.frame specifying the n mage column names that are used
#' for analyses. At the time the function was written, n = 6
#' @param analysis_cols 1 x n data.frame specifying the n package analysis column names
#' that are used for analyses. At the time the function was written, n = 2
#' @param reporter_cols 1 x n data.frame specifying the n reporter ions typed column names
#' that are used for analyses. At the time the function was written, n = 2
#' @param fdata_cols 1 x n data.frame specifying the n f_data column names that are used
#' for analyses. At the time the function was written, n = 6 
#' @param protein_collection_cols 1 x n data.frame specifying the n protein collection
#' column names that are used for analyses. At the time the function was written, n = 1
#'  
#' @return List of isobaricpepObjects that have been cleaned and ready for further labeled chemoproteomics analysis
#' 
#' @examples
#' library(chemoprotR)
#' data("fragpipe_m")
#' data("fragpipe_study_df")
#' 
#' @author Damon Leach, Vlad Petyuk
#' 
#' @export
#' 
clean_chemoprot_fragpipe <- function(m, study_df){
  # check that ReferenceIntensity in m
  if(!"ReferenceIntensity" %in% names(fData(m))){
    stop(paste0("The MSnSet object m must have the following column it's fData object: ReferenceIntensity"))
  }
  
  # check that column names are in study_df
  if(any(!(c("SampleID","ProbeGroup","ExperimentName") %in% names(study_df)))){
    stop (paste0("The data.frame study_df requires the following column names in the
                   data: SampleID, ProbeGroup, and ExperimentName"))
  }
  
  # check that m and study_df have the same values
  # i.e. number of cols in exprs(m) match study_df rows
  if(nrow(pData(m)) != nrow(study_df)){
    stop(paste0("The number of samples in the object m and study_df do not match"))
  }
  
  if(any(!(study_df$SampleID %in% pData(m)$sample_name))){
    stop (paste0("There are sample names in study_df and m that do not match"))
  }
  
  # create edata
  e_data <- sweep(exprs(m), 1, fData(m)$ReferenceIntensity, FUN = "+")
  e_data <- e_data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Protein")
  
  # create fdata
  f_data <- study_df
  
  # create emeta
  e_meta <- fData(m) %>%
    dplyr::rename(Protein = featureName)
  
  # split fdata by ExperimentName
  fragpipe_nest <- f_data %>%
    dplyr::group_by(ExperimentName) %>%
    tidyr::nest() %>%
    dplyr::rename(f_data = data) %>%
    dplyr::mutate(f_data = purrr::map2(ExperimentName, f_data, function(met, fdata){
      fdata <- fdata %>%
        dplyr::mutate(ExperimentName = met)
      return(fdata)
    }))
  
  # split edata by ExperimentName
  fragpipe_nest <- fragpipe_nest %>%
    dplyr::mutate(e_data = purrr::map2(ExperimentName, f_data, function(met, fdata){
      edata <- e_data %>%
        dplyr::select(Protein, dplyr::contains(met))
      return(edata)
    })) %>%
    dplyr::mutate(e_meta = purrr::map(ExperimentName,function(met){
      e_meta
    }))
  
  fragpipe_nest <- fragpipe_nest %>%
    dplyr::mutate(pmartObj = purrr::pmap(list(e_data,f_data,e_meta),function(edat,fdat,emet){
      prodata <- as.proData(e_data = edat,
                            f_data = fdat,
                            e_meta = emet,
                            edata_cname = "Protein",
                            fdata_cname = "SampleID",
                            emeta_cname = "Protein",
                            data_scale = "log2",
                            is_normalized = FALSE)
      return(prodata)
    }))
  
  pmart_list <- fragpipe_nest$pmartObj
  names(pmart_list) <- fragpipe_nest$ExperimentName
  
  return(pmart_list)
}