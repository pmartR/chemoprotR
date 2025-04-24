htp_normalize <- function(pmart_obj,norm_function = "mean", backtransform_val = TRUE){
  
  # check that original_pmart is of appropriate class #
  if (!inherits(pmart_obj, c("pepData", "proData", "metabData", "lipidData",
                             "nmrData"))) {
    
    stop (paste("pmart_obj must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', or 'nmrData'",
                sep = ' '))
  }
  
  # check that original pmart has group designation
  if (is.null(attr(pmart_obj,"group_DF"))){
    
    stop (paste("pmart_obj object must undergo group_designation"))
    
  }
  
  # check that norm_function is a character string that is either mean or median
  if (!is.character(norm_function)) {
    
    stop (paste("norm_function must be a character string indicating 'mean' or 'median'"))
    
  }
  
  if (! norm_function %in% c("mean","median")) {
    
    stop (paste("norm_function can only take the values 'mean' or 'median'"))
    
  }
  
  # check that backtransform_val is a logical argument
  
  if(!is.logical(backtransform_val)) {
    
    stop (paste("backtransform_val must be a logical argument"))
  }
  
  # now run normalize global on the data
  # we use the norm function that the user prefers
  # we do not backtransform - we add this back in later if desired
  norm_glob_data <- pmartR::normalize_global(omicsData = pmart_obj, subset_fn = "all",
                                             norm_fn = norm_function, apply_norm = TRUE, backtransform = FALSE)
  
  # if backtransform_val is FALSE we can just return the data as it is now
  # otherwise we have to calculate the backtransform value back
  if(backtransform_val == FALSE){
    return(norm_glob_data)
  } else {
    # obtain group designation information
    group_info <- attr(pmart_obj,"group_DF")
    
    # find fdata_cname
    fdat_cname <- pmartR::get_fdata_cname(pmart_obj)
    # find edata_cname
    edat_cname <- pmartR::get_edata_cname(pmart_obj)
    # find which column has the edata_cname
    edat_cname_col <- which(colnames(pmart_obj$e_data) == edat_cname)
    # find the ordering
    ordering <- match(colnames(pmart_obj$e_data)[-edat_cname_col],group_info[[fdat_cname]])
    
    # create a data frame with that is in proper ordering of samples
    pg_ordered <- pmart_obj$e_data[-edat_cname_col][,ordering]
    # create a blank data frame to store updated values
    normalized_dat <- data.frame(matrix(data = NA, ncol = ncol(pg_ordered),nrow = nrow(pg_ordered)))
    # add in colname and rowname information for later use
    rownames(normalized_dat) <- pmart_obj$e_data$Peptide
    colnames(normalized_dat) <- colnames(pg_ordered)
    
    # determine whether we are finding the mean or median as the centering value to backtransform by
    if(norm_function == "mean"){
      center_val = apply(pmart_obj$e_data[,-edat_cname_col],2,mean,na.rm=T)
    } else {
      center_val = apply(pmart_obj$e_data[,-edat_cname_col],2,median,na.rm=T)
    }
    
    # create a data frame that obtains the center value
    center_val_df <- data.frame(center_val) %>%
      tibble::rownames_to_column(var = fdat_cname)
    # find the max value for each comparison (this is the one we will backtransform by)
    max_group_val <- group_info %>%
      dplyr::left_join(center_val_df) %>%
      dplyr::group_by(Group) %>%
      dplyr::summarise(max_val = max(center_val))
    # add this max value info into the group information
    group_info <- group_info %>%
      dplyr::left_join(max_group_val)
    # find the proper order between edata and group information
    proper_order <- order(colnames(norm_glob_data$e_data[,-edat_cname_col]),group_info[[fdat_cname]])
    
    # create a data frame that is just the proper max value for each specific sample
    centerM <- data.frame(matrix(group_info$max_val,nrow=dim(norm_glob_data$e_data[,-edat_cname_col])[1],ncol=dim(norm_glob_data$e_data[,-edat_cname_col])[2], byrow=TRUE))
    colnames(centerM) <- group_info[[fdat_cname]][proper_order]
    
    # obtain just the unique max values that we care about for each of the groups to be put in attributes
    distinctVals <- group_info %>% dplyr::select(Group,max_val) %>% dplyr::distinct()
    
    # now adjust the normalized data to be transformed
    pmart_htp_norm <- norm_glob_data
    pmart_htp_norm$e_data[,-edat_cname_col] <- centerM + pmart_htp_norm$e_data[,-edat_cname_col]
    
    # update attributes
    attributes(pmart_htp_norm) <- attributes(pmart_obj)
    attr(pmart_htp_norm,"data_info")$norm_info$is_normalized <- TRUE
    attr(pmart_htp_norm,"data_info")$norm_info$norm_type <- "htp_chemoproteomics"
    attr(pmart_htp_norm,"data_info")$norm_info$subset_fn <- "all"
    attr(pmart_htp_norm,"data_info")$norm_info$subset_params <- NULL
    attr(pmart_htp_norm,"data_info")$norm_info$norm_fn <- norm_function
    attr(pmart_htp_norm,"data_info")$norm_info$n_features_calc <- nrow(pmart_obj$e_data)
    attr(pmart_htp_norm,"data_info")$norm_info$prop_features_calc <- 1
    attr(pmart_htp_norm,"data_info")$norm_info$params$norm_scale <- NULL
    attr(pmart_htp_norm,"data_info")$norm_info$params$norm_location <- NULL
    attr(pmart_htp_norm,"data_info")$norm_info$params$bt_scale <- NULL
    attr(pmart_htp_norm,"data_info")$norm_info$params$bt_location <- distinctVals
    
    # return back transformed normalized data
    return(pmart_htp_norm)
  }
}
