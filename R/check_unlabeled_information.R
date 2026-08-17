#' Check Unlabeled Chemoproteomics Data
#' 
#' This function
#' 
#' @param unlabeled_info unlabeled information that will be passed onto rmarkdown vignette template
#'  
#' @return printed message if no errors occurred
#' 
#' @author Damon Leach, Kelly Stratton
#' 
#' @export
#' 
check_unlabeled_information <- function(unlabeled_info){
  
  if(!is.list(unlabeled_info)){
    stop( paste0("unlabeled_info must be a named list"))
  }
  
  if(any(!c("pmartObj","msgf","fdata_info","mage_info","norm_info","protein_info","data_name",
            "rollup_info","outlier_info","comparisons_info") %in% names(unlabeled_info))){
    stop (paste0("unlabeled_info must be a named list with the following named elements:
                   pmartObj,msgf,fdata_info,mage_info,norm_info,protein_info,data_name,
             rollup_info,outlier_info,comparisons_info"))
  }
  
  # pmartObj
  if(class(unlabeled_info$pmartObj) != "pepData"){
    stop( paste0("The element pmartObj in unlabeled_info must be an omics object of class pepData"))
  }
  
  # msgf
  if(class(unlabeled_info$msgf) != "numeric"){
    stop( paste0("The element msgf in unlabeled_info must be numeric"))
  }
  if(unlabeled_info$msgf < 0){
    stop( paste0("The element msgf in unlabeled_info must be greater than 0"))
  }
  
  # fdata_info
  if(class(unlabeled_info$fdata_info) != "data.frame"){
    stop( paste0("The element fdata_info in unlabeled_info must be a data.frame"))
  }
  if(nrow(unlabeled_info$fdata_info) != 1){
    stop( paste0("The element fdata_info in unlabeled_info must be a data.frame with 1 row"))
  }
  if(any(!c("sampleID_name","group_name","rep_name","job_name") %in% colnames(unlabeled_info$fdata_info))){
    stop (paste0("The element fdata_info in unlabeled_info must be a data.frame with the column names:
                   sampleID_name,group_name,rep_name,job_name"))
  }
  if(any(!unlabeled_info$fdata_info[1,] %in% names(unlabeled_info$pmartObj$f_data))){
    stop (paste0("The values in the data.frame fdata_info must be column names in the f_data
                   element of pmartObj"))
  }
  
  # mage_info
  if(class(unlabeled_info$mage_info) != "data.frame"){
    stop( paste0("The element mage_info in unlabeled_info must be a data.frame"))
  }
  if(nrow(unlabeled_info$mage_info) != 1){
    stop( paste0("The element mage_info in unlabeled_info must be a data.frame with 1 row"))
  }
  if(any(!c("job_name","scan_name","peptide_name","protein_name","qvalue_name","pepqvalue_name",
            "total_ion_intensity_name","peak_area_name","msgf_specprob_name") %in% colnames(unlabeled_info$mage_info))){
    stop (paste0("The element mage_info in unlabeled_info must be a data.frame with the column names:
                   job_name,scan_name,peptide_name,protein_name,qvalue_name,pepqvalue_name,
             total_ion_intensity_name,peak_area_name,msgf_specprob_name"))
  }
  mage_filt1 <- unlabeled_info$mage_info %>% dplyr::select(-c(peak_area_name,job_name))
  mage_filt2 <- unlabeled_info$mage_info %>% dplyr::select(job_name)
  if(any(!mage_filt1[1,] %in% names(unlabeled_info$pmartObj$e_meta))){
    stop (paste0("The values in the data.frame mage_info must be column names in the e_meta
                element of pmartObj (aside from peak_area_name and job_name)"))
  }
  if(any(!mage_filt2[1,] %in% names(unlabeled_info$pmartObj$f_data))){
    stop (paste0("The value job_name in the data.frame mage_info must be a column name in the f_data
                element of pmartObj"))
  }
  
  # norm_info
  if(class(unlabeled_info$norm_info) != "data.frame"){
    stop( paste0("The element norm_info in unlabeled_info must be a data.frame"))
  }
  if(nrow(unlabeled_info$norm_info) != 1){
    stop( paste0("The element norm_info in unlabeled_info must be a data.frame with 1 row"))
  }
  if(any(!c("norm_fn","backtransform") %in% colnames(unlabeled_info$norm_info))){
    stop (paste0("The element norm_info in unlabeled_info must be a data.frame with the column names:
                   norm_fn, backtransform"))
  }
  if(!unlabeled_info$norm_info$norm_fn %in% c("mean","median")){
    stop( paste0("The argument norm_fn in norm_info can only take the following values mean or median"))
  }
  if(!unlabeled_info$norm_info$backtransform %in% c(TRUE, FALSE)){
    stop( paste0("The argument backtransform in norm_info can only take the following values TRUE or FALSE"))
  }
  
  # protein_info
  if(class(unlabeled_info$protein_info) != "data.frame"){
    stop( paste0("The element protein_info in unlabeled_info must be a data.frame"))
  }
  if(nrow(unlabeled_info$protein_info) != 1){
    stop( paste0("The element protein_info in unlabeled_info must be a data.frame with 1 row"))
  }
  if(any(!c("proteinName_name","proteinCollectionID_name","proteinCollection_name","description_name",
            "referenceID_name","residueCount_name","monotopicMass_name","proteinID_name") %in% colnames(unlabeled_info$protein_info))){
    stop (paste0("The element protein_info in unlabeled_info must be a data.frame with the column names:
                   proteinName_name,proteinCollectionID_name,proteinCollection_name,description_name,
                   referenceID_name,residueCount_name,monotopicMass_name,proteinID_name"))
  }
  prot_filt <- unlabeled_info$protein_info[1,] %>% dplyr::select(-proteinName_name)
  if(any(!prot_filt[1,] %in% names(unlabeled_info$pmartObj$e_meta))){
    stop (paste0("The values in the data.frame protein_info must be column names in the e_meta
                   element of pmartObj (except for proteinName_name)"))
  }
  
  # data_name
  if(class(unlabeled_info$data_name) != "character"){
    stop( paste0("The element data_name in unlabeled_info must be a character string"))
  }
  
  # rollup_info
  if(class(unlabeled_info$rollup_info) != "data.frame"){
    stop( paste0("The element rollup_info in unlabeled_info must be a data.frame"))
  }
  if(nrow(unlabeled_info$rollup_info) != 1){
    stop( paste0("The element rollup_info in unlabeled_info must be a data.frame with 1 row"))
  }
  if(any(!c("rollup_method","centering_fn") %in% colnames(unlabeled_info$rollup_info))){
    stop (paste0("The element rollup_info in unlabeled_info must be a data.frame with the column names:
                   rollup_method, centering_fn"))
  }
  if(!unlabeled_info$rollup_info$rollup_method %in% c("rollup","rrollup","summation")){
    stop( paste0("The argument rollup_method in rollup_info can only take the following values rollup, rrollup, or summation"))
  }
  if(!unlabeled_info$rollup_info$centering_fn %in% c("none","median","mean")){
    stop( paste0("The argument backtransform in rollup_info can only take the following values none, median, or mean"))
  }
  if(unlabeled_info$rollup_info$rollup_method != "summation" & unlabeled_info$rollup_info$centering_fn == "none"){
    stop (paste0("The option none for centering_fn is only applicable when rollup_method is summation"))
  }
  
  # outlier_info
  if(length(unlabeled_info$outlier_info) > 0){
    if(class(unlabeled_info$outlier_info) != "character"){
      stop( paste0("The element outlier_info in unlabeled_info must be a character string if not NULL"))
    }
    if(!unlabeled_info$outlier_info %in% unlabeled_info$pmartObj$f_data$job_number){
      stop (paste0("Outliers must be listed in the format of 'job_' and the corresonding job number"))
    }
  }
  
  # comparisons_info
  if(length(unlabeled_info$comparisons_info) > 0){
    if(class(unlabeled_info$comparisons_info) != "data.frame"){
      stop( paste0("The element comparisons_info in unlabeled_info must be a data.frame if not NULL"))
    }
    if(any(!colnames(unlabeled_info$comparisons_info) %in% c("Control","Test"))){
      stop (paste0("If not NULL, the data.frame comparisons_info must take the column names of Control and Test exactly"))
    }
    unique_groups = unique(unlabeled_info$pmartObj$f_data[[unlabeled_info$fdata_info$group_name]])
    if(any(!unlabeled_info$comparisons_info$Control %in% unique_groups)){
      stop (paste0("Some control comparisons in comparisons_df are not a grouping value in the f_data"))
    }
    if(any(!unlabeled_info$comparisons_info$Test %in% unique_groups)){
      stop (paste0("Some test comparisons in comparisons_df are not a grouping value in the f_data"))
    }
    if(any(unlabeled_info$comparisons_info$Test == unlabeled_info$comparisons_info$Control)){
      stop (paste0("Some test and control comparisons are matching. Please remove those comparisons"))
    }
  }
  
  print("No errors detected")
}