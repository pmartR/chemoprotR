#' Clean Unlabeled Chemoproteomics Data
#' 
#' This function takes unlabeled chemoproteomics data from an excel format 
#' and returns an S3 object compatible with another R package, pmartR, that 
#' can then be ran in a markdown file for further filtering, normalization,
#' and statistical analyses
#' 
#' @param dat_list List containing information from each of the different TMT excel sheets
#' @param tab_names 1 x n data.frame specifying the names of the n different excel sheets
#' in the TMT excel file required for analysis. At the time the function was written,
#' n = 4
#' @param mage_cols 1 x n data.frame specifying the n mage column names that are used
#' for analyses. At the time that the function was written, n = 9
#' @param fdata_cols 1 x n data.frame specifying the n f_data column names that are used
#' for analyses. At the time the function was written, n = 9
#' @param protein_collection_cols 1 x n data.frame specifying the protein collection
#' column names that are used for analyses. At the time the function was written, n = 1
#'  
#' @return List of pepObjects that have been cleaned and ready for further unlabeled chemoproteomics analysis
#' 
#' @examples
#' library(chemoprotR)
#' data("labelfree_setup")
#' tb_name <- data.frame(mage = "Mage",
#'                         metadata = "Metadata",
#'                         protein_collection = "Protein collection list",
#'                         fdata = "f_data")
#' mg_col <- data.frame(job_name = "Job",
#'                         scan_name = "Scan",
#'                         peptide_name = "Peptide",
#'                         protein_name = "Protein",
#'                         qvalue_name = "QValue",
#'                         pepqvalue_name = "PepQValue",
#'                         total_ion_intensity_name = "TotalIonIntensity",
#'                         peak_area_name = "PeakArea",
#'                         msgf_specprob_name = "MSGF_SpecProb")
#' f_col <- data.frame(sampleID_name = "SampleID",
#'                          group_name = "Grouping",
#'                          rep_name = "biological_replicate",
#'                          job_name = "Job")
#' pro_col <- data.frame(proteinName_name = "protein_name",
#'                       proteinCollectionID_name = "protein_collection_id",
#'                       proteinCollection_name = "protein_collection",
#'                       description_name = "description",
#'                       referenceID_name = "reference_id",
#'                       residueCount_name = "residue_count",
#'                       monotopicMass_name = "monoisotopic_mass",
#'                       proteinID_name = "protein_id")
#' normalization_info <- data.frame(norm_fn = "mean",
#'                                  backtransform = TRUE)
#' htp_pmart_cleaned <- clean_chemoprot_unlabeled(labelfree_setup,tb_name,mg_col,f_col,pro_col)
#' 
#' @author Damon Leach, Kelly Stratton
#' 
#' @export
#' 
clean_chemoprot_unlabeled <- function(dat_list,tab_names,
                            mage_cols,fdata_cols,
                            protein_collection_cols){
  
  ####################### run through potential warnings #######################
  
  # check that dat_list is a list
  if(!is.list(dat_list)){
    stop (paste("dat_list must be a list and not a data.frame"))
  }
  
  # check that names in tab_names are names of dat_list
  if(sum(!names(dat_list) %in% tab_names) > 0){
    badegg = which(!tab_names %in% names(dat_list))
    stop (paste0("The following argument in tab_names is not a name found in dat_list: ", tab_names[badegg], ". "))
  }
  
  # check mage columns
  if(sum(!(mage_cols) %in% names(dat_list[[tab_names$mage]]))){
    bad_mage = which(!mage_cols %in% names(dat_list[[tab_names$mage]]))
    stop (paste0("The following argument in mage_cols is not a column name found in ", tab_names$mage, ": ", 
                 mage_cols[bad_mage]))
  }
  
  # check fdata columns
  if(sum(!(fdata_cols) %in% names(dat_list[[tab_names$fdata]]))){
    bad_fdata = which(!fdata_cols %in% names(dat_list[[tab_names$fdata]]))
    stop (paste0("The following argument in fdata_cols is not a column name found in ", tab_names$fdata, ": ", 
                 fdata_cols[bad_fdata]))
  }
  
  # # check Ionization samples are the same name as in reporter ions
  #   if(sum(!dat_list[[tab_names$fdata]][[fdata_cols$ionization_name]] %in% colnames(dat_list[[tab_names$reporter_ions]])) > 0){
  #   bad_ionization = which(!dat_list[[tab_names$fdata]][[fdata_cols$ionization_name]] %in% colnames(dat_list[[tab_names$reporter_ions]]))
  #   bad_ionization_name = fdata[[fdata_cols$ionization_name]][bad_ionization]
  #   stop (paste0("The names of the ionization modes must match the columns pertaining to those ionization modes in
  #                T_Reporter_Ions_Typed exactly. If T_Reporter_Ions_Typed lists them as 'Ion_",bad_ionization_name[1],"', they cannot
  #                be listed as '",bad_ionization_name[1],"'in the f_data."))
  # }
  
  # check protein collection columns
  if(sum(!(protein_collection_cols) %in% names(dat_list[[tab_names$protein_collection]]))){
    bad_protein_collection = which(!protein_collection_cols %in% names(dat_list[[tab_names$protein_collection]]))
    stop (paste0("The following argument in protein_collection_cols is not a column name found in ", tab_names$protein_collection, ": ",
                 protein_collection_cols[bad_protein_collection]))
  }
  
  ######################## run the actual processing now #######################
  
  # find which number in the list is the reporter cols sheet
  # find which number in the list is the Mage sheet
  mage_num = which(names(dat_list) == tab_names$mage)
  # find which number in the list is the T_alias sheet
  fdata_num = which(names(dat_list) == tab_names$fdata)
  # find which number in the list is the fdata sheet
  protCollection_num = which(names(dat_list) == tab_names$protein_collection)
  
  # remove unnecessary columns from mage
  # create a string job name variable
  mage_package <- dat_list[[mage_num]]
  mage_package$job_number <- paste0("job_",mage_package[[mage_cols$job_name]])
  mage_package <- mage_package %>%
    # only retain the columns we need
    dplyr::select(job_number,mage_cols$scan_name,mage_cols$peptide,mage_cols$protein,mage_cols$qvalue,
                  mage_cols$pepqvalue,mage_cols$total_ion_intensity,mage_cols$peak_area,mage_cols$msgf_specprob_name)
  
  # remove * from peptides and replace with ""
  mage_package[[mage_cols$peptide_name]] <- gsub("*","",as.character(mage_package[[mage_cols$peptide_name]]), fixed = TRUE)
  
  # make sure that we have unique identifier for each row
  mage_package$ID <- seq.int(nrow(mage_package))
  # make the data wider
  mage_spread <- tidyr::spread(mage_package,job_number,mage_cols$peak_area_name,fill = 0)
  
  # create unique peptide and protein names
  # set up names for unique peptide and unique protein column names
  uniquepep = paste0("Unique_",mage_cols$peptide_name)
  uniqueprot = paste0("Unique_",mage_cols$protein_name)
  mage_spread <- mage_spread %>%
    dplyr::mutate(!!uniquepep := make.unique(mage_spread[[mage_cols$peptide_name]],sep = "__"),
                  !!uniqueprot := make.unique(mage_spread[[mage_cols$protein_name]],sep = "__"))
  
  # save all of the job names
  job_names <- unique(dat_list[[tab_names$mage]][[mage_cols$job_name]])
  
  # now we create the pmart object datasets
  # edata
  edata <- mage_spread %>%
    dplyr::select(!!uniquepep,dplyr::starts_with("job"))
  # emeta
  # find which column in protein collection information corresponds to protein in mage spread
  whichColProtName = which(names(dat_list[[protCollection_num]]) == protein_collection_cols$proteinName_name)
  protCollect = dat_list[[protCollection_num]]
  
  protCollect <- protCollect %>%
    dplyr::mutate(dplyr::across(which(names(protCollect) == protein_collection_cols$proteinName_name),as.character))
  mage_spread <- mage_spread %>%
    dplyr::mutate(dplyr::across(which(names(mage_spread) == mage_cols$protein_name),as.character))

  names(protCollect)[whichColProtName] <- mage_cols$protein_name
  # now left join them together
  emeta <- mage_spread %>%
    dplyr::select(!dplyr::starts_with("job")) %>%
    dplyr::left_join(protCollect, by = mage_cols$protein_name)
  
  # fdata
  fdata <- dat_list[[tab_names$fdata]]
  fdata$job_number <- paste0("job_",fdata[[fdata_cols$job_name]])
  
  # create the pmart object for each job
  # unlike labeled there is only one dataset
  pmart_object <- pmartR::as.pepData(e_data = edata, edata_cname = uniquepep,
                                     f_data = fdata, fdata_cname = "job_number",
                                     e_meta = emeta, emeta_cname = uniqueprot)
  return(pmart_object)
}
