# this function will later be in its own file but lives here currently
# takes in the dat_list object and all the sheet names
clean_chemoprot_labeled <- function(dat_list,tab_names,
                            mage_cols,analysis_cols,reporter_cols, fdata_cols,
                            normalization_info){
  
  ####################### run through potential warnings #######################
  
  # check that dat_list is a list
  if(!is.list(dat_list)){
    stop (paste("dat_list must be a list"))
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
  
  # check analysis columns
  if(sum(!(analysis_cols) %in% names(dat_list[[tab_names$analysis_jobs]]))){
    bad_analysis = which(!analysis_cols %in% names(dat_list[[tab_names$analysis_jobs]]))
    stop (paste0("The following argument in analysis_cols is not a column name found in ", tab_names$analysis_jobs, ": ", 
                 analysis_cols[bad_analysis]))
  }
  
  # check report columns
  if(sum(!(reporter_cols) %in% names(dat_list[[tab_names$reporter_ions]]))){
    bad_reporter = which(!reporter_cols %in% names(dat_list[[tab_names$reporter_ions]]))
    stop (paste0("The following argument in reporter_cols is not a column name found in ", tab_names$reporter_ions, ": ", 
                 reporter_cols[bad_reporter]))
  }
  
  # check fdata columns
  if(sum(!(fdata_cols) %in% names(dat_list[[tab_names$fdata]]))){
    bad_fdata = which(!fdata_cols %in% names(dat_list[[tab_names$fdata]]))
    stop (paste0("The following argument in fdata_cols is not a column name found in ", tab_names$fdata, ": ", 
                 fdata_cols[bad_fdata]))
  }
  
  # check normalization info
  # check that reference name is a value within fdata of
  if(!normalization_info$reference_name %in% dat_list[[tab_names$fdata]][[fdata_cols$group_name]]){
    stop (paste0("reference_name must be a value within the column ", fdata_cols$group_name, " in ", tab_names$fdata))
  }
  
  # check that norm_fn is mean or median
  if(!normalization_info$norm_fn %in% c("mean","median")) {
    stop (paste("norm_fn can only take the value of 'mean' or 'median'."))
  }
  
  # check that backtransform is logical
  if(!is.logical(normalization_info$backtransform)){
    stop (paste("The value backtransform must be a logical value (either TRUE or FALSE)."))
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
  emeta <- Combined_df %>%
    dplyr::select(!dplyr::all_of(ion_names))
  
  # create the pmart object for each job
  pmart_objects <- list()
  for(i in 1:length(unique(unlist(emeta[analysis_cols$job_name])))){
    # find which job we are working with
    job_name <- unique(unlist(emeta[analysis_cols$job_name]))[i]
    
    # subset the jobs to only that job name
    # fancy way to call 
    fdata_mini <- fdata %>% dplyr::filter(!!as.symbol(analysis_cols$job_name) == job_name)
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
