# front end script
# damon leach, kelly stratton, gerard lomas
# april 4, 2023

# UNLABELED ####################################################################
# step 1: load in libraries and data ###########################################
# load in necessary libraries
library(pmartR)
library(readxl)
library(chemoprotR)
library(here)

here::i_am("Code/front_end_script_unlabeled.R")

# step 2: enter in the arguments ###############################################
# mydata is the name of xlsx file that we are loading in

# unlabeled
# mydata is the name of xlsx file that we are loading in
mydata = "60473_Practice_Unlabeled.xlsx"

# enter in the msgf value
msgf = 6.76378E-9

# add in sheet names
tab_names <- data.frame(mage = "Mage",
                        metadata = "Metadata",
                        protein_collection = "Protein collection list",
                        fdata = "f_data")

mage_cols <- data.frame(job_name = "Job",
                        scan_name = "Scan",
                        peptide_name = "Peptide",
                        protein_name = "Protein",
                        qvalue_name = "QValue",
                        pepqvalue_name = "PepQValue",
                        total_ion_intensity_name = "TotalIonIntensity",
                        peak_area_name = "PeakArea",
                        msgf_specprob_name = "MSGF_SpecProb")

fdata_cols <- data.frame(sampleID_name = "SampleID",
                         group_name = "Grouping",
                         rep_name = "biological_replicate",
                         job_name = "Job")

protein_collection_cols <- data.frame(proteinName_name = "protein_name",
                                      proteinCollectionID_name = "protein_collection_id",
                                      proteinCollection_name = "protein_collection",
                                      description_name = "description",
                                      referenceID_name = "reference_id",
                                      residueCount_name = "residue_count",
                                      monotopicMass_name = "monoisotopic_mass",
                                      proteinID_name = "protein_id")

normalization_info <- data.frame(norm_fn = "mean",
                                 backtransform = TRUE)

# rollup_method can take the values: rollup, rrollup, or summation
# centering_fn can take the values: mean, median, (or none if using summation)
rollup_info <- data.frame(rollup_method = "summation",
                          centering_fn = "none")

# list outliers to be removed (recommend leaving this blank for first pass at data
# and then running the report a second time with outliers to remove if there are any)
outlier_samples <- c()

# step 3: create pmart object ##################################################
# put xlsx object into more R friendly formatting
sheet_names <- excel_sheets(here("Data",mydata))
# this creates a list with each element being a sheet by its corresponding name
dat_list <- lapply(sheet_names, function(x) {          # Read all sheets to list
  as.data.frame(read_excel(here("Data",mydata), sheet = x)) } )
names(dat_list) <- sheet_names
# remove xlsx from data_name
data_name <- stringr::str_remove(mydata,".xlsx")

# step 4: run cleaning chemoprot function ######################################
# create pmartObj
htp_pmart_cleaned <- clean_chemoprot_unlabeled(dat_list,tab_names,mage_cols,fdata_cols,
                                               protein_collection_cols)
# save the results and msgf/job information
unlabeled_information <- list(pmartObj = htp_pmart_cleaned,msgf = msgf,
                              fdata_info = fdata_cols,mage_info = mage_cols, norm_info = normalization_info,
                              protein_info = protein_collection_cols,data_name = mydata, rollup_info = rollup_info,
                              outlier_info = outlier_samples)
saveRDS(unlabeled_information,here("Data","unlabeled_information.RDS"))

# remove xlsx from data_name
# step 5: render markdown and results ##########################################
# run this for loop to generate reports for each job
rmarkdown::render(input = here("Code","report_unlabeled.Rmd"),
                  output_file = paste0("report_unlabeled_",data_name,".html"))
