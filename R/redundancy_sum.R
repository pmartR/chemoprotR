#' Peptide Redundancy Summation
#' 
#' This function takes an S3 object of class pepData that contains either
#' labeled or unlabeled chemoproteomics data and sums up peptides to remove
#' redundancies such that each row in the e_data corresponds to a single unique 
#' peptide
#' 
#' @param pmartObj S3 omicsData object of the class 'pepData'
#' @param mage_names 1 x n data.frame specifying the n mage column names that are used
#' for analyses
#' @param protein_collect_names 1 x n data.frame specifying the n protein collection
#' column names that are used for analyses
#' @param labeled logical argument specifying whether the peptide level data
#' is labeled or unlabeled (TRUE indicating labeled and FALSE indicating unlabeled)
#'  
#' @return S3 object that has been summed to remove redundancies
#' 
#' @examples
#' library(chemoprotR)
#' data(label_pmart)
#' prot_collect <- data.frame(proteinName_name = "protein_name")
#' mage_info <- data.frame(job_name = "Job",
#'                         scan_name = "Scan",
#'                         peptide_name = "Peptide",
#'                         protein_name = "Protein",
#'                         qvalue_name = "QValue",
#'                         pepqvalue_name = "PepQValue",
#'                         total_ion_intensity_name = "TotalIonIntensity",
#'                         peak_area_name = "PeakArea",
#'                         msgf_specprob_name = "MSGF_SpecProb")
#' htp_no_redundancy <- redundancy_sum(pmartObj = label_pmart,mage_names = mage_info,
#'                                     protein_collect_names = prot_collect,labeled = TRUE)
#'                                     
#' @author Damon Leach, Kelly Stratton
#' 
#' @export
#' 
redundancy_sum <- function(pmartObj,mage_names,protein_collect_names,labeled = TRUE){
  # check that original pmart is of appropriate class #
  if (!inherits(pmartObj, c("pepData", "proData", "metabData", "lipidData",
                            "nmrData"))) {
    
    stop (paste("pmartObj must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', or 'nmrData'",
                sep = ' '))
  }
  
  # check that protein_collect is data.frame with important info we care about
  if(!is.data.frame(protein_collect_names)){
    stop (paste("protein_collect_names must be a 1xn data.frame"))
  }
  
  # make sure there is only one row
  if(nrow(protein_collect_names) != 1){
    stop (paste("protein_collect_names must be a 1xn data.frame"))
  }
  
  # make sure we have the right names we care about
  prot_names <- c("proteinName_name")
  if(sum(!prot_names %in% colnames(protein_collect_names)) > 0){
    stop( paste("At least one column name in protein_collect_names is missing or mislabeled."))
  }
  
  # check that mage_names is data.frame with important info we care about
  if(!is.data.frame(mage_names)){
    stop (paste("mage_names must be a 1xn data.frame"))
  }
  
  # make sure there is only one row
  if(nrow(mage_names) != 1){
    stop (paste("mage_names must be a 1xn data.frame"))
  }
  
  # make sure we have the right names we care about
  mage_check_names <- c("job_name","scan_name","peptide_name",
                        "protein_name","qvalue_name","msgf_specprob_name")
  if(sum(!mage_check_names %in% colnames(mage_names)) > 0){
    stop( paste("At least one column name in mage_names is missing or mislabeled."))
  }
  
  # check that labeled is logical
  if(!is.logical(labeled)) {
    stop (paste("labeled must be a logical argument"))
  }
  
  ########## now proceed to running code ###########
  # check to see if we are on log scale
  pmartObjOld <- pmartObj
  og_scale = attributes(pmartObj)$data_info$data_scale
  if(og_scale != "abundance"){
    # convert to abundance
    pmartObj <- pmartR::edata_transform(pmartObj,"abundance")
  }
  
  # save values
  edata_cname = pmartR::get_edata_cname(pmartObj)
  fdata_cname = pmartR::get_fdata_cname(pmartObj)
  emeta_cname = pmartR::get_emeta_cname(pmartObj)
  
  # check that labeled is logical TRUE or FALSE
  if(labeled == TRUE){
    # select the columns we care about - edata cname and mage peptide name
    pep_uniquepep <- pmartObj$e_meta %>%
      dplyr::select(!!as.symbol(edata_cname),!!as.symbol(mage_names$peptide_name))
    
    # sum up redundancy for edata
    edatUnique <- suppressWarnings(pmartObj$e_data %>%
      dplyr::left_join(pep_uniquepep) %>%
      dplyr::group_by(!!as.symbol(mage_names$peptide_name)) %>%
      dplyr::summarise(dplyr::across(dplyr::where(is.numeric), sum,na.rm=T)))
    
    # now update emeta object
    emetUnique <- pmartObj$e_meta %>% 
      dplyr::select(-c(!!as.symbol(edata_cname),!!as.symbol(mage_names$msgf_specprob_name),
                       !!as.symbol(mage_names$qvalue_name),!!as.symbol(emeta_cname))) %>%
      dplyr::distinct()
    
    #recombine back into pmart object
    
    # is this still isobaric normalization?
    htp_no_redundancy <- pmartR::as.pepData(e_data = edatUnique, edata_cname = mage_names$peptide_name,
                                            f_data = pmartObj$f_data, fdata_cname = fdata_cname,
                                            e_meta = emetUnique, emeta_cname = mage_names$protein_name)
    
  } else {
    
    # select the columns we care about
    pep_uniquepep <- pmartObj$e_meta %>%
      dplyr::select(!!as.symbol(edata_cname),!!as.symbol(mage_names$peptide_name))
    
    # sum up redundancy for edata
    edatUnique <- suppressWarnings(pmartObj$e_data %>%
      dplyr::left_join(pep_uniquepep) %>%
      dplyr::group_by(!!as.symbol(mage_names$peptide_name)) %>%
      dplyr::summarise(dplyr::across(dplyr::where(is.numeric), sum,na.rm=T)))
    
    # add on protein collection info back in
    add_on_emeta <- pmartObj$e_meta[,c(which(colnames(pmartObj$e_meta) == mage_names$protein_name),
                                       which(colnames(pmartObj$e_meta) %in% protein_collect_names))] %>%
      dplyr::distinct()
    
    # now update emeta object
    emetUnique <- pmartObj$e_meta %>%
      dplyr::group_by(!!as.symbol(mage_names$peptide_name),!!as.symbol(mage_names$protein_name)) %>%
      dplyr::count() %>%
      dplyr::rename(pep_count = n) %>%
      dplyr::left_join(add_on_emeta, by = mage_names$protein_name)
    
    # create pmart object with no redundancy
    htp_no_redundancy <- pmartR::as.pepData(e_data = edatUnique, edata_cname = mage_names$peptide_name,
                                            f_data = pmartObj$f_data, fdata_cname = fdata_cname,
                                            e_meta = emetUnique, emeta_cname = mage_names$protein_name)
  }
  
  # convert back to original scale if not log2
  if(og_scale != "abundance"){
    htp_no_redundancy <- pmartR::edata_transform(htp_no_redundancy,og_scale)
  }
  
  # update attributes
  attr(htp_no_redundancy,"class") <- attr(pmartObjOld,"class")
  attr(htp_no_redundancy,"filters") <- attr(pmartObjOld,"filters")
  attr(htp_no_redundancy,"group_DF") <- attr(pmartObjOld,"group_DF")
  attr(htp_no_redundancy,"data_info")$data_scale = attr(pmartObjOld,"data_info")$data_scale
  attr(htp_no_redundancy,"data_info")$norm_info = attr(pmartObjOld,"data_info")$norm_info
  attr(htp_no_redundancy,"data_info")$data_types = attr(pmartObjOld,"data_info")$data_types
  attr(htp_no_redundancy,"data_info")$batch_info = attr(pmartObjOld,"data_info")$batch_info
  if(!is.null(attr(pmartObjOld,"isobaric_info"))){
    attr(htp_no_redundancy,"isobaric_info") = attr(pmartObjOld,"isobaric_info")
  }
  
  # return pmart obj
  return(htp_no_redundancy)
}