context('Sum Peptide Redundancy')

test_that('redundancy_sum sums up the peptides to a the peptide level as intended',{
  
  # Load the reduced peptide data frames ---------------------------------------
  labeled_pmart <- readRDS(system.file('testdata',
                                       'htp_labeled_pmart.RDS',
                                       package = 'chemoprotR'))
  # Run through potential warnings ---------------------------------------------
  # need information for the other variable names
  prot_collect <- data.frame(proteinName_name = "Protein_Name", proteinCollectionID_name = "Protein_Collection_ID",
                             proteinCollection_name = "Protein_Collection", description_name = "Description",
                             referenceID_name = "Reference_ID", residueCount_name = "Residue_Count",
                             monotopicMass_name = "Monoisotopic_Mass",proteinID_name = "Protein_ID")
  
  mage <- data.frame(job_name = "Job",scan_name = "Scan",peptide_name = "Peptide",
                     protein_name = "Protein",qvalue_name = "QValue",
                     msgf_specprob_name = "MSGF_SpecProb")
  
  # what if we enter in parameters wrong though?
  # pmartObj is not actually a pmartObj
  expect_error(redundancy_sum(labeled_pmart$e_data,protein_collect_names = prot_collect,
                              mage_names = mage,labeled = TRUE),
               "pmartObj must be of class 'pepData', 'proData', 'metabData', 'lipidData', or 'nmrData'")
  
  
  # protein collect names
  # needs to be a data.frame
  expect_error(redundancy_sum(labeled_pmart,protein_collect_names = c("Protein_Name","Prote_Collection_ID"),
                              mage_names = mage,labeled = TRUE),
               "protein_collect_names must be a 1xn data.frame")
  
  # data.frame must only be one row
  prot_collect_manyRows = data.frame(Name = data.frame(Names = c("proteinName_name","proteinCollectionID_name",
                                                                 "proteinCollection_name","description_name",
                                                                 "referenceID_name", "residueCount_name","monotopicMass_name",
                                                                 "proteinID_name"),
                                                       Values = c("Protein_Name","Protein_Collection_ID",
                                                                  "Protein_Collection","Description",
                                                                  "Reference_ID", "Residue_Count",
                                                                  "Monoisotopic_Mass","Protein_ID")))
  expect_error(redundancy_sum(labeled_pmart,protein_collect_names = prot_collect_manyRows,
                              mage_names = mage,labeled = TRUE),
               "protein_collect_names must be a 1xn data.frame")
  
  # we need specific wording for the names in the data.frame
  prot_collect2 = prot_collect
  colnames(prot_collect2)[1] = "protName"
  expect_error(redundancy_sum(labeled_pmart,protein_collect_names = prot_collect2,
                              mage_names = mage,labeled = TRUE),
               "At least one column name in protein_collect_names is missing or mislabeled")
  
  # mage names
  # needs to be a data.frame
  expect_error(redundancy_sum(labeled_pmart,protein_collect_names = prot_collect,
                              mage_names = c("Job","Scan","Peptide"),labeled = TRUE),
               "mage_names must be a 1xn data.frame")
  
  # data.frame must only be one row
  mage_manyRows = data.frame(Name = data.frame(Names = c("job_name","scan_name",
                                                         "peptide_name","protein_name",
                                                         "qvalue_name", "msgf_specprob_name"),
                                               Values = c("Job","Scan",
                                                          "Peptide","Protein",
                                                          "QValue", "MSGF_SpecProb")))
  expect_error(redundancy_sum(labeled_pmart,protein_collect_names = prot_collect,
                              mage_names = mage_manyRows,labeled = TRUE),
               "mage_names must be a 1xn data.frame")
  
  # we need specific wording for the names in the data.frame
  mage2 = mage
  colnames(mage2)[3] = "peptidezName"
  expect_error(redundancy_sum(labeled_pmart,protein_collect_names = prot_collect,
                              mage_names = mage2,labeled = TRUE),
               "At least one column name in mage_names is missing or mislabeled")
  
  # labeled not logical 
  expect_error(redundancy_sum(labeled_pmart,protein_collect_names = prot_collect,
                              mage_names = mage,labeled = "yes"),
               "labeled must be a logical argument")
  
  
  # run the actual thing now  --------------------------------------------------
  pmart_summed = redundancy_sum(pmartObj = labeled_pmart,protein_collect_names = prot_collect,
                                mage_names = mage,labeled = TRUE)
  
  # do this manually too
  # need to convert to abundance
  labeled_pmart <- pmartR::edata_transform(labeled_pmart,"abundance")
  # select the columns we care about - edata cname and mage peptide name
  pep_uniquepep <- labeled_pmart$e_meta %>%
    dplyr::select(Unique_Peptide,Peptide)
  
  # sum up redundancy for edata
  edatUnique <- labeled_pmart$e_data %>%
    dplyr::left_join(pep_uniquepep) %>%
    dplyr::group_by(Peptide) %>%
    dplyr::summarise(dplyr::across(dplyr::where(is.numeric), sum,na.rm=T))
  
  # now update emeta object
  emetUnique <- labeled_pmart$e_meta %>% 
    dplyr::select(-c(Unique_Peptide,MSGF_SpecProb,QValue,Unique_Protein)) %>%
    dplyr::distinct()
  
  #recombine back into pmart object
  labeled_summed <- pmartR::as.pepData(e_data = edatUnique, edata_cname = "Peptide",
                                       f_data = labeled_pmart$f_data, fdata_cname = "SampleID",
                                       e_meta = emetUnique, emeta_cname = "Protein")
  
  # convert back to log2
  labeled_summed <- pmartR::edata_transform(labeled_summed,"log2")
  
  
  # compare the two results
  # edata should be the length of unique peptides and same number of columns
  expect_identical(nrow(labeled_summed$e_data),length(unique(labeled_pmart$e_meta$Peptide)))
  expect_identical(ncol(labeled_summed$e_data),ncol(labeled_pmart$e_data))
  # fdata should stay the same
  expect_identical(dim(labeled_summed$f_data),dim(labeled_pmart$f_data))
  # emeta should be length of unique peptides and 4 less than original
  # as we remove unique_peptide, msgf_specprob, qvalue, and unique_protein
  expect_identical(nrow(labeled_summed$e_meta),length(unique(labeled_pmart$e_meta$Peptide)))
  expect_equal(ncol(labeled_summed$e_meta),ncol(labeled_pmart$e_meta)-4)
  
  # Run manually and compare results -------------------------------------------
  # manual and nonmanual should be the same edata
  expect_equal(labeled_summed$e_data,pmart_summed$e_data)
  
  # Check that attributes update accordingly -----------------------------------
  
  # Inspecticate the attributes of the bc_combat data frame.
  # all of these should be the same between old and new data objects
  
  # isobaric information
  expect_equal(attr(pmart_summed,"isobaric_info"),attr(labeled_pmart,"isobaric_info"))
  
  # filters
  expect_equal(attr(pmart_summed,"filters"),attr(labeled_pmart,"filters"))
  
  # class
  expect_equal(attr(pmart_summed,"class"),attr(labeled_pmart,"class"))
  
  # group_DF
  expect_equal(attr(pmart_summed,"group_DF"),attr(labeled_pmart,"group_DF"))
 
  # these aspects of data_info should be the same
  expect_equal(attr(pmart_summed,"data_info")$data_scale,attr(labeled_summed,"data_info")$data_scale)
  expect_equal(attr(pmart_summed,"data_info")$norm_info,attr(labeled_summed,"data_info")$norm_info)
  expect_equal(attr(pmart_summed,"data_info")$data_types,attr(labeled_summed,"data_info")$data_types)
  expect_equal(attr(pmart_summed,"data_info")$batch_info,attr(labeled_summed,"data_info")$batch_info)
  })

