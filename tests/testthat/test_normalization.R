context('HTP Normalization')

test_that('htp_normalize normalizes the data as expected',{
  
  # Load the reduced peptide data frames ---------------------------------------
  labeled_pmart <- readRDS(system.file('testdata',
                   'htp_labeled_pmart.RDS',
                   package = 'chemoprotR'))
  
  # Run through potential warnings ---------------------------------------------
  
  # we do not have group designation yet so this should fail
  expect_error(htp_normalize(labeled_pmart,norm_function = "mean",
                             backtransform_val = TRUE),
               "pmartObj object must undergo group_designation")
  
  # add in group designation
  labeled_pmart <- pmartR::group_designation(labeled_pmart, main_effects = "Grouping")
  
  # what if we enter in parameters wrong though?
  # pmartObj is not actually a pmartObj
  expect_error(htp_normalize(labeled_pmart$e_data,norm_function = "mean",
                             backtransform_val = TRUE),
               "pmartObj must be of class 'pepData', 'proData', 'metabData', 'lipidData', or 'nmrData'")
  
  # norm function not mean or median
  expect_error(htp_normalize(labeled_pmart,norm_function = "mode",
                             backtransform_val = TRUE),
               "norm_function can only take the values 'mean' or 'median'")
  
  expect_error(htp_normalize(labeled_pmart,norm_function = 1,
                             backtransform_val = TRUE),
               "norm_function must be a character string indicating 'mean' or 'median'")
  
  # backtransform is not logical
  expect_error(htp_normalize(labeled_pmart,norm_function = "median",
                             backtransform_val = "yes"),
               "backtransform_val must be a logical argument")

  # Check dimensions of results  -----------------------------------------------
  
  # SCENARIO 1: BACKTRANSFORM = TRUE, NORM_FN = MEAN
  # run the function
  labeled_normalized <- htp_normalize(pmartObj = labeled_pmart,norm_function = "mean",
                                      backtransform_val = TRUE)
  
  # check edata, fdata, and emeta all stay the same dimensions as before
  expect_identical(dim(labeled_normalized$e_data),dim(labeled_pmart$e_data))
  expect_identical(dim(labeled_normalized$f_data),dim(labeled_pmart$f_data))
  expect_identical(dim(labeled_normalized$e_meta),dim(labeled_pmart$e_meta))
  
  # Run manually and compare results -------------------------------------------
  
  norm_glob_data <- pmartR::normalize_global(omicsData = labeled_pmart, subset_fn = "all",
                                             norm_fn = "mean", apply_norm = TRUE, backtransform = FALSE)
  # obtain group designation information
  group_info <- attr(labeled_pmart,"group_DF")
  
  # find fdata_cname
  fdat_cname <- pmartR::get_fdata_cname(labeled_pmart)
  # find edata_cname
  edat_cname <- pmartR::get_edata_cname(labeled_pmart)
  # find which column has the edata_cname
  edat_cname_col <- which(colnames(labeled_pmart$e_data) == edat_cname)
  # find the ordering
  ordering <- match(colnames(labeled_pmart$e_data)[-edat_cname_col],group_info[[fdat_cname]])
  
  # create a data frame with that is in proper ordering of samples
  pg_ordered <- labeled_pmart$e_data[-edat_cname_col][,ordering]
  # create a blank data frame to store updated values
  normalized_dat <- data.frame(matrix(data = NA, ncol = ncol(pg_ordered),nrow = nrow(pg_ordered)))
  # add in colname and rowname information for later use
  rownames(normalized_dat) <- labeled_pmart$e_data$Peptide
  colnames(normalized_dat) <- colnames(pg_ordered)
  
  # find center value using mean
  center_val = apply(labeled_pmart$e_data[,-edat_cname_col],2,mean,na.rm=T)
  
  # create a data frame that obtains the center value
  center_val_df <- data.frame(center_val) %>%
    tibble::rownames_to_column(var = fdat_cname)
  # find the max value for each comparison (this is the one we will backtransform by)
  max_group_val <- group_info %>%
    dplyr::left_join(center_val_df) %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(Backtransform = max(center_val))
  # add this max value info into the group information
  group_info <- group_info %>%
    dplyr::left_join(max_group_val)
  # find the proper order between edata and group information
  proper_order <- order(colnames(norm_glob_data$e_data[,-edat_cname_col]),group_info[[fdat_cname]])
  
  # create a data frame that is just the proper max value for each specific sample
  centerM <- data.frame(matrix(group_info$Backtransform,nrow=dim(norm_glob_data$e_data[,-edat_cname_col])[1],ncol=dim(norm_glob_data$e_data[,-edat_cname_col])[2], byrow=TRUE))
  colnames(centerM) <- group_info[[fdat_cname]][proper_order]
  
  # obtain just the unique max values that we care about for each of the groups to be put in attributes
  norm_info_attr <- group_info %>% dplyr::select(SampleID,Backtransform)
  
  # now adjust the normalized data to be transformed
  pmart_htp_norm <- norm_glob_data
  pmart_htp_norm$e_data[,-edat_cname_col] <- centerM + pmart_htp_norm$e_data[,-edat_cname_col]
  
  # these values should be the same
  expect_equal(pmart_htp_norm$e_data,labeled_normalized$e_data)
  
  # Check that attributes update accordingly -----------------------------------
  
  # Inspecticate the attributes of the bc_combat data frame.
  # all of these should be the same between old and new data objects
  # check names
  expect_equal(attr(labeled_normalized,"names"),attr(labeled_pmart,"names"))
  
  # check cnames
  expect_equal(attr(labeled_normalized,"cnames"),attr(labeled_pmart,"cnames"))
  
  # isobaric information
  expect_equal(attr(labeled_normalized,"isobaric_info"),attr(labeled_pmart,"isobaric_info"))
  
  # check.names
  expect_equal(attr(labeled_normalized,"check.names"),attr(labeled_pmart,"check.names"))
  
  # meta_info
  expect_equal(attr(labeled_normalized,"meta_info"),attr(labeled_pmart,"meta_info"))
  
  # filters
  expect_equal(attr(labeled_normalized,"filters"),attr(labeled_pmart,"filters"))
  
  # class
  expect_equal(attr(labeled_normalized,"class"),attr(labeled_pmart,"class"))
  
  # group_DF
  expect_equal(attr(labeled_normalized,"group_DF"),attr(labeled_pmart,"group_DF"))
  
  # data info should be different though with regards to norm_type
  expect_equal(attributes(labeled_normalized)$data_info$data_scale_orig,
               attributes(labeled_pmart)$data_info$data_scale_orig)
  expect_equal(attributes(labeled_normalized)$data_info$data_scale,
               attributes(labeled_pmart)$data_info$data_scale)
  expect_equal(attributes(labeled_normalized)$data_info$num_edata,
               attributes(labeled_pmart)$data_info$num_edata)
  expect_equal(attributes(labeled_normalized)$data_info$num_miss_obs,
               attributes(labeled_pmart)$data_info$num_miss_obs)
  expect_equal(attributes(labeled_normalized)$data_info$prop_missing,
               attributes(labeled_pmart)$data_info$prop_missing)
  expect_equal(attributes(labeled_normalized)$data_info$num_samps,
               attributes(labeled_pmart)$data_info$num_samps)
  expect_equal(attributes(labeled_normalized)$data_info$data_types,
               attributes(labeled_pmart)$data_info$data_types)
  expect_equal(attributes(labeled_normalized)$data_info$batch_info,
               attributes(labeled_pmart)$data_info$batch_info)
  # norm info will differ though
  expect_error(attributes(labeled_normalized)$data_info$norm_info == 
               attributes(labeled_pmart)$data_info$norm_info)
  
  # labeled_normalized should have updated information regarding normalization
  expect_equal(attributes(labeled_normalized)$data_info$norm_info$is_normalized,
               TRUE)
  expect_equal(attributes(labeled_normalized)$data_info$norm_info$norm_type,
               "global")
  expect_equal(attributes(labeled_normalized)$data_info$norm_info$subset_fn,
               "all")
  expect_equal(attributes(labeled_normalized)$data_info$norm_info$norm_fn,
               "mean")
  expect_equal(attributes(labeled_normalized)$data_info$norm_info$n_features_calc,
               nrow(labeled_normalized$e_data))
  expect_equal(attributes(labeled_normalized)$data_info$norm_info$prop_features_calc,
               1)
  expect_equal(attributes(labeled_normalized)$data_info$norm_info$params$bt_location,
               norm_info_attr)

  # SCENARIO 2: BACKTRANSFORM = FALSE, NORM_FN = MEDIAN
  # run the function
  labeled_normalized <- htp_normalize(pmartObj = labeled_pmart,norm_function = "median",
                                      backtransform_val = FALSE)
  
  # check edata, fdata, and emeta all stay the same dimensions as before
  expect_identical(dim(labeled_normalized$e_data),dim(labeled_pmart$e_data))
  expect_identical(dim(labeled_normalized$f_data),dim(labeled_pmart$f_data))
  expect_identical(dim(labeled_normalized$e_meta),dim(labeled_pmart$e_meta))
  
  # Run manually and compare results -------------------------------------------
  
  norm_glob_data <- pmartR::normalize_global(omicsData = labeled_pmart, subset_fn = "all",
                                             norm_fn = "median", apply_norm = TRUE, backtransform = FALSE)
  
  # these values should be the same
  expect_equal(norm_glob_data$e_data,labeled_normalized$e_data)
  
  # Check that attributes update accordingly -----------------------------------
  
  # Inspecticate the attributes of the bc_combat data frame.
  # all of these should be the same between old and new data objects
  # check names
  expect_equal(attr(labeled_normalized,"names"),attr(labeled_pmart,"names"))
  
  # check cnames
  expect_equal(attr(labeled_normalized,"cnames"),attr(labeled_pmart,"cnames"))
  
  # isobaric information
  expect_equal(attr(labeled_normalized,"isobaric_info"),attr(labeled_pmart,"isobaric_info"))
  
  # check.names
  expect_equal(attr(labeled_normalized,"check.names"),attr(labeled_pmart,"check.names"))
  
  # meta_info
  expect_equal(attr(labeled_normalized,"meta_info"),attr(labeled_pmart,"meta_info"))
  
  # filters
  expect_equal(attr(labeled_normalized,"filters"),attr(labeled_pmart,"filters"))
  
  # class
  expect_equal(attr(labeled_normalized,"class"),attr(labeled_pmart,"class"))
  
  # group_DF
  expect_equal(attr(labeled_normalized,"group_DF"),attr(labeled_pmart,"group_DF"))
  
  # data info should be different though with regards to norm_type
  expect_equal(attributes(labeled_normalized)$data_info$data_scale_orig,
               attributes(labeled_pmart)$data_info$data_scale_orig)
  expect_equal(attributes(labeled_normalized)$data_info$data_scale,
               attributes(labeled_pmart)$data_info$data_scale)
  expect_equal(attributes(labeled_normalized)$data_info$num_edata,
               attributes(labeled_pmart)$data_info$num_edata)
  expect_equal(attributes(labeled_normalized)$data_info$num_miss_obs,
               attributes(labeled_pmart)$data_info$num_miss_obs)
  expect_equal(attributes(labeled_normalized)$data_info$prop_missing,
               attributes(labeled_pmart)$data_info$prop_missing)
  expect_equal(attributes(labeled_normalized)$data_info$num_samps,
               attributes(labeled_pmart)$data_info$num_samps)
  expect_equal(attributes(labeled_normalized)$data_info$data_types,
               attributes(labeled_pmart)$data_info$data_types)
  expect_equal(attributes(labeled_normalized)$data_info$batch_info,
               attributes(labeled_pmart)$data_info$batch_info)
  # norm info will differ though
  expect_error(attributes(labeled_normalized)$data_info$norm_info == 
                 attributes(labeled_pmart)$data_info$norm_info)
  
  # labeled_normalized should have updated information regarding normalization
  expect_equal(attributes(labeled_normalized)$data_info$norm_info$is_normalized,
               TRUE)
  expect_equal(attributes(labeled_normalized)$data_info$norm_info$norm_type,
               "global")
  expect_equal(attributes(labeled_normalized)$data_info$norm_info$subset_fn,
               "all")
  expect_equal(attributes(labeled_normalized)$data_info$norm_info$norm_fn,
               "median")
  expect_equal(attributes(labeled_normalized)$data_info$norm_info$n_features_calc,
               nrow(labeled_normalized$e_data))
  expect_equal(attributes(labeled_normalized)$data_info$norm_info$prop_features_calc,
               1)
  expect_equal(attributes(labeled_normalized)$data_info$norm_info$params$bt_location,
               NULL)
  
})

