# this function will later be in its own file but lives here currently
# takes in the dat_list object and all the sheet names
clean_chemoprot_unlabeled <- function(dat_list,tab_names,
                            mage_cols,fdata_cols,protein_collection_cols){
  
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
  
  # check fdata columns
  if(sum(!(fdata_cols) %in% names(dat_list[[tab_names$fdata]]))){
    bad_fdata = which(!fdata_cols %in% names(dat_list[[tab_names$fdata]]))
    stop (paste0("The following argument in fdata_cols is not a column name found in ", tab_names$fdata, ": ", 
                 fdata_cols[bad_fdata]))
  }
  
  # check protein collection columns
  if(sum(!(protein_collection_cols) %in% names(dat_list[[tab_names$protein_collection]]))){
    bad_protein_collection = which(!protein_collection_cols %in% names(dat_list[[tab_names$protein_collection]]))
    stop (paste0("The following argument in protein_collection_cols is not a column name found in ", tab_names$protein_collection, ": ",
                 protein_collection_cols[bad_protein_collection]))
  }
  
  # check normalization info
  # check that norm_fn is mean or median
  if(!normalization_info$norm_fn %in% c("mean","median")) {
    stop (paste("norm_fn can only take the value of 'mean' or 'median'."))
  }
  
  # check that backtransform is logical
  if(!is.logical(normalization_info$backtransform)){
    stop (paste("The value backtransform must be a logical value (either TRUE or FALSE)."))
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
    dplyr::select(job_number,mage_cols$msgf_specprob,mage_cols$peptide,mage_cols$protein,mage_cols$qvalue,
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
