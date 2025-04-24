context('Clean Unlabeled Chemoproteomics Data')

test_that('cleans unlabeled chemoproteomics to be able to use with pmartR functions',{
  
  # Load the reduced peptide data frames ---------------------------------------
  unlabeled_data <- readRDS(system.file('testdata',
                                      'example_htp_unlabeled.RDS',
                                      package = 'chemoprotR'))
  
  # Run through potential warnings ---------------------------------------------
  
  # add in sheet names
  tab_info <- data.frame(mage = "Mage",
                          metadata = "Metadata",
                          protein_collection = "Protein_Collection",
                          fdata = "f_data")
  
  mage_info <- data.frame(job_name = "Job",
                          scan_name = "Scan",
                          peptide_name = "Peptide",
                          protein_name = "Protein",
                          qvalue_name = "QValue",
                          pepqvalue_name = "PepQValue",
                          total_ion_intensity_name = "TotalIonIntensity",
                          peak_area_name = "PeakArea",
                          msgf_specprob_name = "MSGF_SpecProb")
  
  fdata_info <- data.frame(sampleID_name = "SampleID",
                           group_name = "Grouping",
                           rep_name = "Replicate",
                           job_name = "Job")
  
  prot_info <- data.frame(proteinName_name = "Protein_Name",
                                        proteinCollectionID_name = "Protein_Collection_ID",
                                        proteinCollection_name = "Protein_Collection",
                                        description_name = "Description",
                                        referenceID_name = "Reference_ID",
                                        residueCount_name = "Residue_Count",
                                        monotopicMass_name = "Monoisotopic_Mass",
                                        proteinID_name = "Protein_ID")
  
  
  # must be a list and not data.frame
  expect_error(clean_chemoprot_unlabeled(dat_list = unlabeled_data$Metatdata,
                                         tab_names = tab_info,
                                         mage_cols = mage_info,
                                         fdata_cols = fdata_info,
                                         protein_collection_cols = prot_info),
               "dat_list must be a list and not a data.frame")

  # tab_info must be the names of dat_list
  expect_error(clean_chemoprot_unlabeled(dat_list = unlabeled_data,
                                         tab_names = mage_info,
                                         mage_cols = mage_info,
                                         fdata_cols = fdata_info,
                                         protein_collection_cols = prot_info),
               "The following argument in tab_names is not a name found in")
  
  # mage_info must be column names in Mage
  expect_error(clean_chemoprot_unlabeled(dat_list = unlabeled_data,
                                         tab_names = tab_info,
                                         mage_cols = tab_info,
                                         fdata_cols = fdata_info,
                                         protein_collection_cols = prot_info),
               "The following argument in mage_cols is not a column name found in")
  
  # fdata_info must be column names in fdata
  expect_error(clean_chemoprot_unlabeled(dat_list = unlabeled_data,
                                         tab_names = tab_info,
                                         mage_cols = mage_info,
                                         fdata_cols = mage_info,
                                         protein_collection_cols = prot_info),
               "The following argument in fdata_cols is not a column name found in")
  
  # protein collection information
  expect_error(clean_chemoprot_unlabeled(dat_list = unlabeled_data,
                                         tab_names = tab_info,
                                         mage_cols = mage_info,
                                         fdata_cols = fdata_info,
                                         protein_collection_cols = fdata_info),
               "The following argument in")
  
  # run the actual thing now  --------------------------------------------------

  cleaning_unlabeled <- clean_chemoprot_unlabeled(dat_list = unlabeled_data,
                                                  tab_names = tab_info,
                                                  mage_cols = mage_info,
                                                  fdata_cols = fdata_info,
                                                  protein_collection_cols = prot_info)
  
  # Run manually ---------------------------------------------------------------
  
  # remove unnecessary columns from mage
  # create a string job name variable
  mage_package <- unlabeled_data$Mage
  mage_package$job_number <- paste0("job_",mage_package$Job)
  mage_package <- mage_package %>%
    # only retain the columns we need
    dplyr::select(job_number,Scan,Peptide,Protein,QValue,
                  PepQValue,TotalIonIntensity,PeakArea,MSGF_SpecProb)
  
  # remove * from peptides and replace with ""
  mage_package$Peptide <- gsub("*","",as.character(mage_package$Peptide), fixed = TRUE)
  
  # make sure that we have unique identifier for each row
  mage_package$ID <- seq.int(nrow(mage_package))
  # make the data wider
  mage_spread <- tidyr::spread(mage_package,job_number,PeakArea,fill = 0)
  
  # create unique peptide and protein names
  # set up names for unique peptide and unique protein column names
  mage_spread <- mage_spread %>%
    dplyr::mutate(Unique_Peptide = make.unique(mage_spread$Peptide,sep = "__"),
                  Unique_Protein = make.unique(mage_spread$Protein,sep = "__"))
  
  # save all of the job names
  job_names <- unique(unlabeled_data$Mage$Job)
  
  # now we create the pmart object datasets
  # edata
  edata <- mage_spread %>%
    dplyr::select(Unique_Peptide,dplyr::starts_with("job"))
  # emeta
  # find protein collection information
  protCollect = unlabeled_data$Protein_Collection
  # now left join them together
  emeta <- mage_spread %>%
    dplyr::select(!dplyr::starts_with("job")) %>%
    dplyr::left_join(protCollect, by = c("Protein" = "Protein_Name"))
  
  # fdata
  fdata <- unlabeled_data$f_data
  fdata$job_number <- paste0("job_",fdata$Job)
  
  # create the pmart object for each job
  # unlike labeled there is only one dataset
  pmart_object <- pmartR::as.pepData(e_data = edata, edata_cname = "Unique_Peptide",
                                     f_data = fdata, fdata_cname = "job_number",
                                     e_meta = emeta, emeta_cname = "Unique_Protein")
  
  # compare results ------------------------------------------------------------
  # manual and nonmanual should be the same edata
  expect_identical(cleaning_unlabeled,pmart_object)
})
