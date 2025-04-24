context('Protein Rollup Summation')

test_that('protein rollup for summation method runs as expected',{
  
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
  
  # run protein redundancy
  pmart_no_redundancy = redundancy_sum(pmartObj = labeled_pmart,protein_collect_names = prot_collect,
                                mage_names = mage,labeled = TRUE)
  # group designation
  pmart_no_redundancy <- pmartR::group_designation(pmart_no_redundancy,main_effects = "Grouping")
  # normalize data
  pmart_norm <- htp_normalize(pmartObj = pmart_no_redundancy,norm_function = "mean", backtransform_val = TRUE)
  
  # what if we enter parameters wrong
  # must be pepData
  expect_error(sum_pep_quant(pmart_norm$e_data, top3 = F),
               "pepData must be an object of class 'pepData")
  # top3 must be logical argument
  expect_error(sum_pep_quant(pmart_norm, top3 = "False"),
               "top3 must be a logical argument")

  # run the actual thing now (sum_pep_quant) -----------------------------------
  pmart_summation <- sum_pep_quant(pmart_norm, top3 = F)
  
  # do this manually too
  # convert to abundance first
  pmart_abund <- pmartR::edata_transform(pmart_norm,"abundance")
  pmart_summation_manual <- pmart_abund$e_data %>%
    dplyr::left_join(pmart_norm$e_meta %>% dplyr::select(Peptide,Protein), by = "Peptide") %>%
    dplyr::group_by(Protein) %>%
    dplyr::summarise(dplyr::across(dplyr::where(is.numeric),sum,na.rm=TRUE)) %>%
    # convert back to log
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), log2))
  pmart_summation_manual[pmart_summation_manual == "-Inf"] <- NA
  
  # compare by first arranging them by protein name
  pmart_summation_df <- pmart_summation$e_data %>% dplyr::arrange(Protein) %>% data.frame()
  pmart_summation_manual_df <- pmart_summation_manual %>% dplyr::arrange(Protein) %>% data.frame()
  
  # these two data.frames should now be equal to each other
  expect_identical(pmart_summation_df,pmart_summation_manual_df)
  
  # what if top3 = TRUE
  pmart_summation_top3 <- sum_pep_quant(pmart_norm, top3 = T)
  
  # do this manually too
  pmart_summation_manual_top3 <- pmart_abund$e_data %>%
    dplyr::left_join(pmart_norm$e_meta %>% dplyr::select(Peptide,Protein), by = "Peptide") %>%
    dplyr::select(-Peptide) %>%
    dplyr::group_by(Protein) %>%
    dplyr::mutate(dplyr::across(
      .cols = -dplyr::any_of("Protein"),
      .fns = function(x) if(all(is.na(x))) NA else 
        sum(rev(sort(x))[1:(min(3, length(x)))], na.rm = T)
    )) %>%
    dplyr::distinct() %>%
    data.frame() %>%
    dplyr::relocate(Protein, .before = dplyr::everything()) %>%
    # convert back to log
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), log2))
  pmart_summation_manual_top3[pmart_summation_manual_top3 == "-Inf"] <- NA
  
  # compare by first arranging them by protein name
  pmart_summation_top3_df <- pmart_summation_top3$e_data %>% dplyr::arrange(Protein) %>% data.frame()
  pmart_summation_manual_top3_df <- pmart_summation_manual_top3 %>% dplyr::arrange(Protein) %>% data.frame()
  
  # these two data.frames should now be equal to each other
  expect_identical(pmart_summation_top3_df,pmart_summation_manual_top3_df)
  
  # run the actual thing now (sum_pep_function) --------------------------------
  # this assumes everything is on the same scale (so keep it at abundance across the board)
  pmart_summation_internal <- sum_pep_function(pmart_abund, top3 = F)
  
  # do this manually too
  # convert to abundance first
  #pmart_abund <- pmartR::edata_transform(pmart_norm,"abundance")
  pmart_summation_manual_internal <- pmart_abund$e_data %>%
    dplyr::left_join(pmart_norm$e_meta %>% dplyr::select(Peptide,Protein), by = "Peptide") %>%
    dplyr::group_by(Protein) %>%
    dplyr::summarise(dplyr::across(dplyr::where(is.numeric),sum,na.rm=TRUE))
  pmart_summation_manual_internal[pmart_summation_manual_internal == 0] <- NA
  
  # compare by first arranging them by protein name
  pmart_summation_internal_df <- pmart_summation_internal$e_data %>% dplyr::arrange(Protein) %>% data.frame()
  pmart_summation_manual_internal_df <- pmart_summation_manual_internal %>% dplyr::arrange(Protein) %>% data.frame()
  
  # these two data.frames should now be equal to each other
  expect_identical(pmart_summation_internal_df,pmart_summation_manual_internal_df)
  
  # what if top3 = TRUE
  pmart_summation_top3_internal <- sum_pep_function(pmart_abund, top3 = T)
  
  # do this manually too
  pmart_summation_manual_top3_internal <- pmart_abund$e_data %>%
    dplyr::left_join(pmart_norm$e_meta %>% dplyr::select(Peptide,Protein), by = "Peptide") %>%
    dplyr::select(-Peptide) %>%
    dplyr::group_by(Protein) %>%
    dplyr::mutate(dplyr::across(
      .cols = -dplyr::any_of("Protein"),
      .fns = function(x) if(all(is.na(x))) NA else 
        sum(rev(sort(x))[1:(min(3, length(x)))], na.rm = T)
    )) %>%
    dplyr::distinct() %>%
    data.frame() %>%
    dplyr::relocate(Protein, .before = dplyr::everything())
  pmart_summation_manual_top3_internal[pmart_summation_manual_top3_internal == 0] <- NA
  
  # compare by first arranging them by protein name
  pmart_summation_top3_internal_df <- pmart_summation_top3_internal$e_data %>% dplyr::arrange(Protein) %>% data.frame()
  pmart_summation_manual_top3_internal_df <- pmart_summation_manual_top3_internal %>% dplyr::arrange(Protein) %>% data.frame()
  
  # these two data.frames should now be equal to each other
  expect_identical(pmart_summation_top3_internal_df,pmart_summation_manual_top3_internal_df)
  
  # Check that attributes update accordingly -----------------------------------
  # names
  expect_equal(attributes(pmart_summation)$names,attributes(pmart_norm)$names)
  # cnames
  expect_equal(attributes(pmart_summation)$cnames,list(edata_cname = "Protein",emeta_cname = "Protein",fdata_cname = "SampleID",techrep_cname = NULL))
  
  # data_info
  expect_equal(attributes(pmart_summation)$data_info$data_scale_orig,"abundance")
  expect_equal(attributes(pmart_summation)$data_info$data_scale,"log2")
  expect_equal(attributes(pmart_summation)$data_info$norm_info$is_normalized,TRUE)
  expect_equal(attributes(pmart_summation)$data_info$num_edata,nrow(pmart_summation$e_data))
  expect_equal(attributes(pmart_summation)$data_info$num_miss_obs,sum(is.na(pmart_summation$e_data)))
  expect_equal(attributes(pmart_summation)$data_info$prop_missing,sum(is.na(pmart_summation$e_data[,-1]))/((sum(is.na(pmart_summation$e_data[,-1])) + sum(!is.na(pmart_summation$e_data[,-1])))))
  expect_equal(attributes(pmart_summation)$data_info$num_samps,nrow(pmart_summation$f_data))
  expect_equal(attributes(pmart_summation)$data_info$data_types,NULL)
  expect_equal(attributes(pmart_summation)$data_info$batch_info$is_bc,FALSE)
  
  # meta_info
  expect_equal(attributes(pmart_summation)$meta_info,attributes(pmart_norm)$meta_info)
  
  # filters
  # these get removed for proData after being converted from pepData
  expect_equal(attributes(pmart_summation)$filters,list())
  
  # pro_quant_info
  expect_equal(attributes(pmart_summation)$pro_quant_info$method,"sum")
  
  # class (move from pepData to proData)
  expect_equal(attributes(pmart_norm)$class,c("isobaricpepData","pepData"))
  expect_equal(attributes(pmart_summation)$class,"proData")
  
  # group_DF
  expect_equal(attributes(pmart_summation)$group_DF,attributes(pmart_norm)$group_DF)
  
  # now for pep sum function
  expect_equal(length(pmart_summation_internal),2)
  expect_equal(names(pmart_summation_internal),c("e_data","e_meta"))
  })

