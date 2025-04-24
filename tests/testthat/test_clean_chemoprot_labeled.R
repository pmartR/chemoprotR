context('Clean Labeled Chemoproteomics Data')

test_that('cleans labeled chemoproteomics to be able to use with pmartR functions',{
  
  # Load the reduced peptide data frames ---------------------------------------
  labeled_data <- readRDS(system.file('testdata',
                                       'example_htp_labeled.RDS',
                                       package = 'chemoprotR'))
  
  # Run through potential warnings ---------------------------------------------
    
  # tab_info
  tab_info <- data.frame(talias = "T_alias",
                          tdata = "T_Data",
                          analysis_jobs = "T_Data_Package_Analysis_Jobs",
                          reporter_ions = "T_Reporter_Ions_Typed",
                          mage = "Mage",
                          protein_collection = "Sheet1",
                          fdata = "f_data")
  # mage_info
  mage_info <- data.frame(job_name = "Job",
                          scan_name = "Scan",
                          peptide_name = "Peptide",
                          protein_name = "Protein",
                          qvalue_name = "QValue",
                          msgf_specprob_name = "MSGF_SpecProb")
  # analysis_info
  analysis_info <- data.frame(job_name = "Job",
                              dataset_id_name = "Dataset_ID")
  # reporter_info
  reporter_info <- data.frame(dataset_id_name = "Dataset",
                              scan_name = "ScanNumber")
  # fdata_info
  fdata_info <- data.frame(sampleID_name = "SampleID",
                           group_name = "Grouping",
                           job_name = "Job",
                           replicate_name = "Replicate",
                           plex_name = "Plex",
                           ionization_name = "Ionization")
  # prot_info
  prot_info <- data.frame(proteinName_name = "Protein_Name",
                                        proteinCollectionID_name = "Protein_Collection_ID",
                                        proteinCollection_name = "Protein_Collection",
                                        description_name = "Description",
                                        referenceID_name = "Reference_ID",
                                        residueCount_name = "Residue_Count",
                                        monotopicMass_name = "Monoisotopic_Mass",
                                        proteinID_name = "Protein_ID")
  
  # must be a list and not data.frame
  expect_error(clean_chemoprot_labeled(dat_list = labeled_data$T_Data,
                                       tab_names = tab_info,
                                       mage_cols = mage_info,
                                       analysis_cols = analysis_info,
                                       reporter_cols = reporter_info,
                                       fdata_cols = fdata_info,
                                       protein_collection_cols = prot_info),
               "dat_list must be a list and not a data.frame")
  
  # tab_info must be the names of dat_list
  expect_error(clean_chemoprot_labeled(dat_list = labeled_data,
                                       tab_names = mage_info,
                                       mage_cols = mage_info,
                                       analysis_cols = analysis_info,
                                       reporter_cols = reporter_info,
                                       fdata_cols = fdata_info,
                                       protein_collection_cols = prot_info),
               "The following argument in tab_names is not a name found in")
  
  # mage_info must be column names in Mage
  expect_error(clean_chemoprot_labeled(dat_list = labeled_data,
                                       tab_names = tab_info,
                                       mage_cols = tab_info,
                                       analysis_cols = analysis_info,
                                       reporter_cols = reporter_info,
                                       fdata_cols = fdata_info,
                                       protein_collection_cols = prot_info),
               "The following argument in mage_cols is not a column name found in")
  
  # analysis_info must be column names in analysis
  expect_error(clean_chemoprot_labeled(dat_list = labeled_data,
                                       tab_names = tab_info,
                                       mage_cols = mage_info,
                                       analysis_cols = tab_info,
                                       reporter_cols = reporter_info,
                                       fdata_cols = fdata_info,
                                       protein_collection_cols = prot_info),
               "The following argument in analysis_cols is not a column name found in")
  
  # reporter_info must be column names in reporter
  expect_error(clean_chemoprot_labeled(dat_list = labeled_data,
                                       tab_names = tab_info,
                                       mage_cols = mage_info,
                                       analysis_cols = analysis_info,
                                       reporter_cols = analysis_info,
                                       fdata_cols = fdata_info,
                                       protein_collection_cols = prot_info),
               "The following argument in reporter_cols is not a column name found in")

  # fdata_info must be column names in fdata
  expect_error(clean_chemoprot_labeled(dat_list = labeled_data,
                                       tab_names = tab_info,
                                       mage_cols = mage_info,
                                       analysis_cols = analysis_info,
                                       reporter_cols = reporter_info,
                                       fdata_cols = reporter_info,
                                       protein_collection_cols = prot_info),
               "The following argument in fdata_cols is not a column name found in")
  
  # ionization samples must match between reporter ions and fdata
  labeled_data2 <- labeled_data
  labeled_data2$f_data$Ionization <- stringr::str_remove(labeled_data2$f_data$Ionization,"Ion_")
  expect_error(clean_chemoprot_labeled(dat_list = labeled_data2,
                                       tab_names = tab_info,
                                       mage_cols = mage_info,
                                       analysis_cols = analysis_info,
                                       reporter_cols = reporter_info,
                                       fdata_cols = fdata_info,
                                       protein_collection_cols = prot_info),
               "The names of the ionization modes must match the columns pertaining to those ionization modes")
  
  # protein collection information
  expect_error(clean_chemoprot_labeled(dat_list = labeled_data,
                                       tab_names = tab_info,
                                       mage_cols = mage_info,
                                       analysis_cols = analysis_info,
                                       reporter_cols = reporter_info,
                                       fdata_cols = fdata_info,
                                       protein_collection_cols = fdata_info),
               "The following argument in")
  
  # run the actual thing now  --------------------------------------------------
  cleaning_labeled <- clean_chemoprot_labeled(dat_list = labeled_data,
                                               tab_names = tab_info,
                                               mage_cols = mage_info,
                                               analysis_cols = analysis_info,
                                               reporter_cols = reporter_info,
                                               fdata_cols = fdata_info,
                                               protein_collection_cols = prot_info)
  
  # Run manually ---------------------------------------------------------------
  
  # find the unique ion names
  ion_names <- unique(labeled_data$f_data$Ionization)
  
  # merge mage and analysis jobs data frames together (both have "Job" columns)
  mage_package <- labeled_data$T_Data_Package_Analysis_Jobs %>%
    dplyr::left_join(labeled_data$Mage, by = "Job",multiple = "all") %>%
    # only retain the columns we need
    dplyr::select(Job,Dataset_ID,QValue,MSGF_SpecProb,Peptide,Protein,Scan)
  # add a new column that matches datasetID and scan (so we can merge with reporter ions)
  mage_package$DatasetScan = paste0(mage_package$Dataset_ID,"_",mage_package$Scan)
  
  # edit reporter_ions using "T_Reporter_Ions_Typed" sheet
  reporter_ions_scan <- labeled_data$T_Reporter_Ions_Typed
  # add a dataset and scan number column to match mage_package
  reporter_ions_scan$DatasetScan = paste0(reporter_ions_scan$Dataset,"_",reporter_ions_scan$ScanNumber)
  
  # merge that information with mage package dataset
  # inner join is used as we want the intersection between both datasets
  merged_df <- dplyr::inner_join(mage_package,reporter_ions_scan, by = "DatasetScan")
  
  # clean up the dataframe by removing "*" from peptides and replace with ""
  merged_df$Peptide <- gsub("*","",as.character(merged_df$Peptide), fixed = TRUE)
  # only retain the columns that we need ('Job', 'QValue', 'Peptide', 'Protein', 'Ion Values')
  # we also add columns to create unique peptides and proteins for pmart-friendly setup
  # set up names for unique peptide and unique protein column names
  Combined_df <- merged_df %>%
    dplyr::select(Job,QValue,Peptide,Protein,MSGF_SpecProb,dplyr::all_of(ion_names)) %>%
    dplyr::mutate(Unique_Peptide = make.unique(merged_df$Peptide,sep = "__"),
                  Unique_Protein = make.unique(merged_df$Protein,sep = "__"))
  
  # now we create the pmart datasets
  # create the edata object from Combined_df (unique peptide and ion values)
  edata <- Combined_df %>%
    dplyr::select(Unique_Peptide,Job,dplyr::starts_with("Ion"))
  
  # create the emeta object from the non Ion values in Combined_df
  # find which column in protein collection information corresponds to protein in mage spread
  protCollect = labeled_data$Sheet1
  names(protCollect)[3] <- "Protein"
  emeta <- Combined_df %>%
    dplyr::select(-dplyr::starts_with("Ion")) %>%
    dplyr::left_join(protCollect, by = "Protein")
  
  # create the pmart object for each job
  pmart_objects <- list()
  for(i in 1:length(unique(unlist(emeta$Job)))){
    # find which job we are working with
    job_name <- unique(unlist(emeta$Job))[i]
    
    # subset the jobs to only that job name
    # fancy way to call 
    fdata_mini <- labeled_data$f_data %>% dplyr::filter(Job == job_name)
    edata_mini <- edata %>% dplyr::filter(Job == job_name) %>% dplyr::select(-Job)
    uniquepep_col <- which(colnames(edata_mini) == "Unique_Peptide")
    sample_ordering <- match(colnames(edata_mini)[-uniquepep_col],fdata_mini$Ionization)
    colnames(edata_mini)[-uniquepep_col] <- fdata_mini$SampleID[sample_ordering]
    emeta_mini <- emeta %>% dplyr::filter(Job == job_name)
    pmart_objects[[i]] <- pmartR::as.isobaricpepData(e_data = edata_mini, edata_cname = "Unique_Peptide",
                                                     f_data = fdata_mini, fdata_cname = "SampleID",
                                                     e_meta = emeta_mini, emeta_cname = "Unique_Protein")
  }
  names(pmart_objects) <- paste0("Job_",unique(unlist(emeta$Job)))
  
  # compare results ------------------------------------------------------------
  # manual and nonmanual should be the same edata
  expect_identical(cleaning_labeled$Job_2044574,pmart_objects$Job_2044574)
  expect_identical(cleaning_labeled$Job_2044575,pmart_objects$Job_2044575)
  expect_identical(cleaning_labeled$Job_2044576,pmart_objects$Job_2044576)
  
})

