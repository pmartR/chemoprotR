#' Clean Labeled Chemoproteomics Data
#' 
#' This function takes labeled chemoproteomics data from an excel format 
#' and returns an S3 object compatible with another R package, pmartR, that 
#' can then be ran in a markdown file for further filtering, normalization,
#' and statistical analyses
#' 
#' @param dat_list List containing information from each of the different TMT excel sheets
#' @param tab_names 1 x n data.frame specifying the names of the n different excel sheets
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
#' data("label_pmart")
#' tb_name <- data.frame(analysis_jobs = "T_Data_Package_Analysis_Jobs",
#'                         reporter_ions = "T_Reporter_Ions_Typed",
#'                         mage = "Mage",
#'                         protein_collection = "Protein_collection_data",
#'                         fdata = "f_data")
#' mg_col <- data.frame(job_name = "Job",
#'                         scan_name = "Scan",
#'                         peptide_name = "Peptide",
#'                         protein_name = "Protein",
#'                         qvalue_name = "QValue",
#'                         msgf_specprob_name = "MSGF_SpecProb")
#' anlys_col <- data.frame(job_name = "Job",
#'                             dataset_id_name = "Dataset_ID")
#' rep_col <- data.frame(dataset_id_name = "Dataset",
#'                             scan_name = "ScanNumber")
#' f_col <- data.frame(sampleID_name = "SampleID",
#'                          group_name = "Grouping",
#'                          job_name = "Job",
#'                          replicate_name = "Replicate",
#'                          plex_name = "Plex",
#'                          ionization_name = "Ionization")
#' pro_col <- data.frame(proteinName_name = "protein_name")
#' htp_clean <- clean_chemoprot_labeled(label_setup,tb_name,mg_col,anlys_col,rep_col,f_col,pro_col)
#' 
#' @author Damon Leach, Kelly Stratton
#' 
#' @export
#' 
clean_chemoprot_labeled <- function(dat_list,tab_names,
                                    mage_cols,analysis_cols,reporter_cols, fdata_cols,
                                    protein_collection_cols){
  
  ####################### run through potential warnings #######################
  
  # check that dat_list is a list
  if(!is.list(dat_list) | is.data.frame(dat_list)){
    stop (paste("dat_list must be a list and not a data.frame"))
  }
  
  # check that names in tab_names are names of dat_list
  if(sum(which(!tab_names %in% names(dat_list))) > 0){
    badegg = which(!tab_names %in% names(dat_list))
    stop (paste0("The following argument in tab_names is not a name found in dat_list: ", tab_names[badegg], ". "))
  }
  
  # check mage columns
  if(sum(!(mage_cols) %in% names(dat_list[[tab_names$mage]])) > 0){
    bad_mage = which(!mage_cols %in% names(dat_list[[tab_names$mage]]))
    stop (paste0("The following argument in mage_cols is not a column name found in ", tab_names$mage, ": ", 
                 mage_cols[bad_mage]))
  }
  
  # check analysis columns
  if(sum(!(analysis_cols) %in% names(dat_list[[tab_names$analysis_jobs]])) > 0){
    bad_analysis = which(!analysis_cols %in% names(dat_list[[tab_names$analysis_jobs]]))
    stop (paste0("The following argument in analysis_cols is not a column name found in ", tab_names$analysis_jobs, ": ", 
                 analysis_cols[bad_analysis]))
  }
  
  # check report columns
  if(sum(!(reporter_cols) %in% names(dat_list[[tab_names$reporter_ions]])) > 0){
    bad_reporter = which(!reporter_cols %in% names(dat_list[[tab_names$reporter_ions]]))
    stop (paste0("The following argument in reporter_cols is not a column name found in ", tab_names$reporter_ions, ": ", 
                 reporter_cols[bad_reporter]))
  }
  
  # check fdata columns
  if(sum(!(fdata_cols) %in% names(dat_list[[tab_names$fdata]])) > 0){
    bad_fdata = which(!fdata_cols %in% names(dat_list[[tab_names$fdata]]))
    stop (paste0("The following argument in fdata_cols is not a column name found in ", tab_names$fdata, ": ", 
                 fdata_cols[bad_fdata]))
  }
  
  # check Ionization samples are the same name as in reporter ions
  if(sum(!dat_list[[tab_names$fdata]][[fdata_cols$ionization_name]] %in% colnames(dat_list[[tab_names$reporter_ions]])) > 0){
    bad_ionization = which(!dat_list[[tab_names$fdata]][[fdata_cols$ionization_name]] %in% colnames(dat_list[[tab_names$reporter_ions]]))
    bad_ionization_name = dat_list[[tab_names$fdata]][[fdata_cols$ionization_name]][bad_ionization]
    stop (paste0("The names of the ionization modes must match the columns pertaining to those ionization modes in
                 T_Reporter_Ions_Typed exactly. If T_Reporter_Ions_Typed lists them 
                 as 'Ion_",bad_ionization_name[1],"', they cannot
                 be listed as '",bad_ionization_name[1],"'in the f_data."))
  }
  
  # check protein collection columns
  if(sum(!(protein_collection_cols) %in% names(dat_list[[tab_names$protein_collection]])) > 0){
    bad_analysis = which(!protein_collection_cols %in% names(dat_list[[tab_names$protein_collection]]))
    stop (paste0("The following argument in protein_collection_cols is not a column name found in ", tab_names$protein_collection, ": ", 
                 protein_collection_cols[bad_analysis]))
  }
  
  ######################## run the actual processing now #######################
  
  # find which number in the list is the T_Data_Package_Analysis_Jobs sheet
  analysis_jobs_num = which(names(dat_list) == tab_names$analysis_jobs)
  # find which number in the list is the Mage sheet
  mage_num = which(names(dat_list) == tab_names$mage)
  # find which number in the list is the T_Reporter_Ions_Typed sheet
  reporter_ions_num = which(names(dat_list) == tab_names$reporter_ions)
  # find which number in the list is the T_alias sheet
  fdata_num = which(names(dat_list) == tab_names$fdata)
  # create the fdata object using t_alias information
  fdata <- dat_list[[fdata_num]]
  # find which number in the list is the fdata sheet
  protCollection_num = which(names(dat_list) == tab_names$protein_collection)
  
  # find the unique ion names
  ion_names <- unique(fdata[[fdata_cols$ionization_name]])
  
  # merge mage and analysis jobs data frames together (both have "Job" columns)
  mage_package <- dat_list[[analysis_jobs_num]] %>%
    dplyr::left_join(dat_list[[mage_num]], by = stats::setNames(mage_cols$job_name, analysis_cols$job_name),multiple = "all") %>%
    # only retain the columns we need
    dplyr::select(analysis_cols$job_name,analysis_cols$dataset_id_name,mage_cols$qvalue_name,mage_cols$msgf_specprob_name,
                  mage_cols$peptide_name,mage_cols$protein_name,mage_cols$scan_name)
  # add a new column that matches datasetID and scan (so we can merge with reporter ions)
  mage_package$DatasetScan = paste0(mage_package[[analysis_cols$dataset_id_name]],"_",mage_package[[mage_cols$scan_name]])
  
  # edit reporter_ions using "T_Reporter_Ions_Typed" sheet
  reporter_ions_scan <- dat_list[[reporter_ions_num]]
  # add a dataset and scan number column to match mage_package
  reporter_ions_scan$DatasetScan = paste0(reporter_ions_scan[[reporter_cols$dataset_id_name]],"_",reporter_ions_scan[[reporter_cols$scan_name]])
  
  # merge that information with mage package dataset
  # inner join is used as we want the intersection between both datasets
  merged_df <- dplyr::inner_join(mage_package,reporter_ions_scan, by = "DatasetScan")
  
  # clean up the dataframe by removing "*" from peptides and replace with ""
  merged_df[[mage_cols$peptide_name]] <- gsub("*","",as.character(merged_df[[mage_cols$peptide_name]]), fixed = TRUE)
  # only retain the columns that we need ('Job', 'QValue', 'Peptide', 'Protein', 'Ion Values')
  # we also add columns to create unique peptides and proteins for pmart-friendly setup
  # set up names for unique peptide and unique protein column names
  uniquepep = paste0("Unique_",mage_cols$peptide_name)
  uniqueprot = paste0("Unique_",mage_cols$protein_name)
  Combined_df <- merged_df %>%
    dplyr::select(analysis_cols$job_name,mage_cols$qvalue_name,mage_cols$peptide_name,
                  mage_cols$protein_name,mage_cols$msgf_specprob_name,dplyr::all_of(ion_names)) %>%
    # fix peptide, protein names
    # sample identifier column name
    dplyr::mutate(!!uniquepep := make.unique(merged_df[[mage_cols$peptide_name]],sep = "__"),
                  !!uniqueprot := make.unique(merged_df[[mage_cols$protein_name]],sep = "__"))
  
  # now we create the pmart datasets
  # create the edata object from Combined_df (unique peptide and ion values)
  edata <- Combined_df %>%
    dplyr::select(dplyr::all_of(uniquepep),analysis_cols$job_name,dplyr::all_of(ion_names))
  
  # create the emeta object from the non Ion values in Combined_df
  # find which column in protein collection information corresponds to protein in mage spread
  whichColProtName = which(names(dat_list[[protCollection_num]]) == protein_collection_cols$proteinName_name)
  protCollect = dat_list[[protCollection_num]]
  
  # ensure protein names are characters (not numbers)
  protCollect <- protCollect %>%
    dplyr::mutate(dplyr::across(which(names(protCollect) == protein_collection_cols$proteinName_name),as.character))
  merged_df <- merged_df %>%
    dplyr::mutate(dplyr::across(which(names(merged_df) == mage_cols$protein_name),as.character))
  
  names(protCollect)[whichColProtName] <- mage_cols$protein_name
  emeta <- Combined_df %>%
    dplyr::select(!dplyr::all_of(ion_names)) %>%
    dplyr::left_join(protCollect, by = mage_cols$protein_name)
  
  # create the pmart object for each job
  pmart_objects <- list()
  for(i in 1:length(unique(unlist(emeta[analysis_cols$job_name])))){
    # find which job we are working with
    job_name <- unique(unlist(emeta[analysis_cols$job_name]))[i]
    
    # subset the jobs to only that job name
    # fancy way to call 
    fdata_mini <- fdata %>% dplyr::filter(!!as.symbol(fdata_cols$job_name) == job_name)
    edata_mini <- edata %>% dplyr::filter(!!as.symbol(analysis_cols$job_name) == job_name) %>% dplyr::select(-!!as.symbol(analysis_cols$job_name))
    uniquepep_col <- which(colnames(edata_mini) == uniquepep)
    sample_ordering <- match(colnames(edata_mini)[-uniquepep_col],fdata_mini[[fdata_cols$ionization_name]])
    colnames(edata_mini)[-uniquepep_col] <- fdata_mini[[fdata_cols$sampleID_name]][sample_ordering]
    emeta_mini <- emeta %>% dplyr::filter(!!as.symbol(analysis_cols$job_name) == job_name)
    pmart_objects[[i]] <- pmartR::as.isobaricpepData(e_data = edata_mini, edata_cname = uniquepep,
                                                     f_data = fdata_mini, fdata_cname = fdata_cols$sampleID_name,
                                                     e_meta = emeta_mini, emeta_cname = uniqueprot)
  }
  names(pmart_objects) <- paste0("Job_",unique(unlist(emeta[analysis_cols$job_name])))
  # return the list of different pmart objects
  return(pmart_objects)
}
