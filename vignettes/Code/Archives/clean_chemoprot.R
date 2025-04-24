# this function will later be in its own file but lives here currently
# takes in the dat_list object and all the sheet names
clean_chemoprot <- function(dat_list,tab_names,
                            mage_cols,analysis_cols,reporter_cols){
  # find which number in the list is the T_Data_Package_Analysis_Jobs sheet
  analysis_jobs_num = which(names(dat_list) == tab_names$analysis_jobs)
  # find which number in the list is the Mage sheet
  mage_num = which(names(dat_list) == tab_names$mage)
  # find which number in the list is the T_Reporter_Ions_Typed sheet
  reporter_ions_num = which(names(dat_list) == tab_names$reporter_ions)
  # find which number in the list is the T_alias sheet
  fdata_num = which(names(dat_list) == tab_names$fdata)
  
  # merge mage and analysis jobs data frames together (both have "Job" columns)
  mage_package <- dat_list[[analysis_jobs_num]] %>%
    dplyr::left_join(dat_list[[mage_num]], by = stats::setNames(mage_cols$job_name, analysis_cols$job_name)) %>%
    # only retain the columns we need
    dplyr::select(analysis_cols$job_name,analysis_cols$dataset_id_name,analysis_cols$qvalue_name,analysis_cols$msgf_specprob_name,
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
  Combined_df <- merged_df %>%
    dplyr::select(analysis_cols$job_name,analysis_cols$qvalue_name,mage_cols$peptide_name,
                  mage_cols$protein_name,analysis_cols$msgf_specprob_name,dplyr::starts_with("Ion")) %>%
    dplyr::mutate(UniquePeptide = make.unique(Peptide,sep = "__"),
                  UniqueProtein = make.unique(Protein,sep = "__"))
  
  # now we create the pmart datasets
  # create the edata object from Combined_df (unique peptide and ion values)
  edata <- Combined_df %>%
    dplyr::select(UniquePeptide,analysis_cols$job_name,dplyr::starts_with("Ion"))
  
  # create the emeta object from the non Ion values in Combined_df
  emeta <- Combined_df %>%
    dplyr::select(!dplyr::starts_with("Ion"))
  
  # create the fdata object using t_alias information
  fdata <- dat_list[[fdata_num]]

  # create the pmart object for each job
  pmart_objects <- list()
  for(i in 1:length(unique(unlist(emeta[analysis_cols$job_name])))){
    # find which job we are working with
    job_name <- unique(unlist(emeta[analysis_cols$job_name]))[i]
    
    # subset the jobs to only that job name
    # fancy way to call 
    fdata_mini <- fdata %>% dplyr::filter(!!as.symbol(analysis_cols$job_name) == job_name)
    edata_mini <- edata %>% dplyr::filter(!!as.symbol(analysis_cols$job_name) == job_name) %>% dplyr::select(-!!as.symbol(analysis_cols$job_name))
    emeta_mini <- emeta %>% dplyr::filter(!!as.symbol(analysis_cols$job_name) == job_name)
    pmart_objects[[i]] <- pmartR::as.isobaricpepData(e_data = edata_mini, edata_cname = "UniquePeptide",
                                                     f_data = fdata_mini, fdata_cname = "SampleID",
                                                     e_meta = emeta_mini, emeta_cname = "UniqueProtein")
  }
  names(pmart_objects) <- paste0("Job_",unique(unlist(emeta[analysis_cols$job_name])))
  # return the list of different pmart objects
  return(pmart_objects)
}
