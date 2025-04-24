# UNLABELED ####################################################################
# step 1: load in libraries and data ###########################################
# load in necessary libraries
library(pmartR)
library(readxl)
#library(chemoprotR)

# find the data directory
datadir = paste0(here::here(),"/Data")

# step 2: enter in the arguments ###############################################
# data is the name of xlsx file that we are loading in

# unlabeled
# data is the name of xlsx file that we are loading in
data = "LabelFree_Practice_Dataset.xlsx"

# enter in the msgf value
msgf = 2.7113E-8

# add in sheet names
tab_names <- data.frame(mage = "Mage",
                        metadata = "Metadata",
                        protein_collection = "Protein_Collection",
                        fdata = "f_data")

mage_cols <- data.frame(job_name = "Job",
                        scan_name = "MSGF_SpecProb",
                        peptide_name = "Peptide",
                        protein_name = "Protein",
                        qvalue_name = "QValue",
                        pepqvalue_name = "PepQValue",
                        total_ion_intensity_name = "TotalIonIntensity",
                        peak_area_name = "PeakArea",
                        msgf_specprob_name = "MSGF_SpecProb")

fdata_cols <- data.frame(sampleID_name = "SampleID",
                         group_name = "Grouping",
                         rep_name = "Replicate",
                         job_name = "Job")

protein_collection_cols <- data.frame(proteinName_name = "Protein_Name",
                                      proteinCollectionID_name = "Protein_Collection_ID",
                                      proteinCollection_name = "Protein_Collection",
                                      description_name = "Description",
                                      referenceID_name = "Reference_ID",
                                      residueCount_name = "Residue_Count",
                                      monotopicMass_name = "Monoisotopic_Mass",
                                      proteinID_name = "Protein_ID")

normalization_info <- data.frame(norm_fn = "mean",
                                 backtransform = TRUE)

# step 3: create pmart object ##################################################
# put xlsx object into more R friendly formatting
sheet_names <- excel_sheets(paste0(datadir,"/",data))
# this creates a list with each element being a sheet by its corresponding name
dat_list <- lapply(sheet_names, function(x) {          # Read all sheets to list
  as.data.frame(read_excel(paste0(datadir,"/",data), sheet = x)) } )
names(dat_list) <- sheet_names

# we also want the protein collection subset information

# step 4: run cleaning chemoprot function ######################################
# create pmartObj
htp_pmart_cleaned <- clean_chemoprot_unlabeled(dat_list,tab_names,mage_cols,fdata_cols,protein_collection_cols)
# save the results and msgf/job information
unlabeled_information <- list(pmartObj = htp_pmart_cleaned,msgf = msgf,
                              fdata_info = fdata_cols,mage_info = mage_cols, norm_info = normalization_info,
                              protein_info = protein_collection_cols)
saveRDS(unlabeled_information,paste0(datadir,"/unlabeled_information.RDS"))

# step 5: render markdown and results ##########################################
# run this for loop to generate reports for each job
rmarkdown::render(input = paste0(here::here(),"/Code/report_unlabeled.Rmd"),
                  output_file = paste0("report_unlabeled",".html"))
