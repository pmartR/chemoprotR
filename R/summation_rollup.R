#' Sum Peptide Helper Function
#' 
#' This function aids in the process of summing up peptide abundance values
#' for rollup via summation
#' 
#' @param pepData an omicsData object of the class 'pepData'
#' @param top3 logical indicating whether to use top 3 most abundant values or all values
#'  
#' @return list with information regarding edata and emeta
#'
#' @author Rachel Richardson, Damon Leach, Kelly Stratton
#' 
#' @export
#' 
sum_pep_function <- function(pepData, top3) {
  # check that pepData is of appropriate class #
  if (!inherits(pepData, "pepData")) {
    stop("pepData must be an object of class 'pepData'")
  }
  
  # check that a protein mapping is provided #
  if (is.null(pepData$e_meta)) {
    stop(paste("A mapping to proteins must be provided in order to use the",
               "protein_filter function.",
               sep = " "
    ))
  }
  
  # Fish out the e_data, f_data, and e_meta column names corresponding to the
  # peptide, sample, and protein IDs.
  pep_id <- attr(pepData, "cnames")$edata_cname
  samp_id = attr(pepData, "cnames")$fdata_cname
  pro_id <- attr(pepData, "cnames")$emeta_cname
  
  # Quantitate the heck out of the peptides!
  if(top3){
    
    res <- merge(
      x = pepData$e_meta[, c(pep_id, pro_id)],
      y = pepData$e_data,
      by = pep_id,
      all.x = FALSE,
      all.y = TRUE
    ) %>%
      dplyr::select(-dplyr::sym(pep_id)) %>%
      dplyr::group_by(!!dplyr::sym(pro_id)) %>%
      dplyr::mutate(dplyr::across(
        .cols = -dplyr::any_of(pro_id),
        .fns = function(x) if(all(is.na(x))) NA else 
          sum(rev(sort(x))[1:(min(3, length(x)))], na.rm = T)
      )) %>%
      dplyr::distinct() %>%
      data.frame(check.names = FALSE)
    
  } else {
    
    res <- merge(
      x = pepData$e_meta[, c(pep_id, pro_id)],
      y = pepData$e_data,
      by = pep_id,
      all.x = FALSE,
      all.y = TRUE
    ) %>%
      dplyr::select(-dplyr::sym(pep_id)) %>%
      dplyr::group_by(!!dplyr::sym(pro_id)) %>%
      dplyr::mutate(dplyr::across(
        .cols = -dplyr::any_of(pro_id),
        .fns = function(x) if(all(is.na(x))) NA else sum(x, na.rm = T)
      )) %>%
      dplyr::distinct() %>%
      data.frame(check.names = FALSE)
    
  }
  
  return(
    list(
      e_data = res,
      e_meta = NULL
    )
  )
}

#' Sum Peptide Quantitation
#' 
#' This function sums up peptide abundance value for rollup
#' 
#' @param pepData an omicsData object of the class 'pepData'
#' @param top3 logical indicating whether to use top 3 most abundant values or all values (default is FALSE)
#'  
#' @return an omicsData object of the class 'proData'
#'
#' @author Rachel Richardson, Damon Leach, Kelly Stratton
#' 
#' @export
#' 
sum_pep_quant <- function(pepData, top3 = F) {
  
  # checks
  if (!inherits(pepData, "pepData")) {
    stop("pepData must be an object of class 'pepData'")
  }
  
  # check that a protein mapping is provided #
  if (is.null(pepData$e_meta)) {
    stop(paste("A mapping to proteins must be provided in order to use the",
               "protein_filter function.",
               sep = " "
    ))
  }
  
  # check that top3 is a logical argument
  if(!is.logical(top3)) {
    stop(paste("top3 must be a logical argument (either TRUE or FALSE)"))
  }
  
  # Extract attribute info to be used throughout the function ------------------
  
  # Pull out column names from e_data, f_data, and e_meta.
  edata_cname <- attr(pepData, "cnames")$edata_cname
  fdata_cname <- attr(pepData, "cnames")$fdata_cname
  emeta_cname <- attr(pepData, "cnames")$emeta_cname
  
  # Extricate e_data column name index.
  edata_cname_id <- which(names(pepData$e_data) == edata_cname)
  
  # Grab more attributes that will be used at some point somewhere.
  data_scale <- pmartR::get_data_scale(pepData)
  is_normalized <- attr(pepData, "data_info")$norm_info$is_normalized
  
  pepData <- pmartR::edata_transform(pepData, "abundance")
  # Quantitate the heck out of the peptides ------------------------------------
  results <- sum_pep_function(pepData, top3)
  
  # Update e_meta after quantitation -------------------------------------------
  results$e_meta <- pepData$e_meta
  
  # Check if isoformRes is NULL. results$e_meta will be updated differently
  # depending on whether isoformRes is present.
  # Update e_meta with peptide counts.
  results$e_meta <- results$e_meta %>%
    dplyr::group_by(!!dplyr::sym(emeta_cname)) %>%
    dplyr::mutate(peps_per_pro = dplyr::n()) %>%
    # Only qrollup will cause n_peps_used != peps_per_protein when isoformRes
    # is NULL.
    dplyr::mutate(
      n_peps_used = if ("n_peps_used" %in% colnames(results$e_meta)) {
        n_peps_used
      } else {
        peps_per_pro
      }
    ) %>%
    # Move n_pep_used to the end. This line will only make a change to the
    # order of the columns when qrollup is selected.
    dplyr::relocate(n_peps_used, .after = dplyr::last_col()) %>%
    # Keep the mapping variable and new columns created plus any columns
    # specified by the user.
    dplyr::select(dplyr::any_of(c(
      emeta_cname, "peps_per_pro",
      "n_peps_used"
    ))) %>%
    # Only keep distinct combinations of the columns that are kept.
    dplyr::distinct() %>%
    data.frame(check.names = FALSE)
  
  # The following runs when isoformRes is present. In this case n_peps_used
  # will be calculated based on protein isoform instead of protein (which
  # includes all isoforms).
  
  # Create a proData object ----------------------------------------------------
  
  prodata <- pmartR::as.proData(
    e_data = results$e_data,
    f_data = pepData$f_data,
    e_meta = results$e_meta,
    edata_cname = emeta_cname,
    fdata_cname = fdata_cname,
    emeta_cname = emeta_cname,
    data_scale = "abundance",
    is_normalized = is_normalized
  )
  
  # Update proData attributes --------------------------------------------------
  
  # Update the original data scale for the proData object. This needs to be
  # manually updated because the original data scale for the proData object will
  # be set to the current data scale of the pepData object when the as.proData
  # function was called. If these two scales are different this is the only way
  # to set the original data scale for the proData object to the original data
  # scale of the pepData object.
  attr(prodata, "data_info")$data_scale_orig <- pmartR::get_data_scale_orig(pepData)
  
  # Update the group_DF attribute (if it exists). This attribute will be "reset"
  # when the as.proData function is called in the rollup functions. It will need
  # to be manually updated to reflect anything done to the peptide data before
  # protein_quant.
  attr(prodata, "group_DF") <- attr(pepData, "group_DF")
  
  # Update the pro_quant_info attribute to reflect which rollup method was used.
  attr(prodata, "pro_quant_info")$method <- "sum"
  
  prodata <- pmartR::edata_transform(prodata, "log2")
  
  return(prodata)
}