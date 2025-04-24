## -----------------------------------------------------------------------------
library(readxl)
library(chemoprotR)
library(here)

## -----------------------------------------------------------------------------
here::i_am("Code/front_end_script_labeled.R")

## ----eval = FALSE-------------------------------------------------------------
# # data is the name of xlsx file that we are loading in
# data = "TMT_Practice_Dataset.xlsx"

## ----echo = FALSE-------------------------------------------------------------
data("label_setup")
data = label_setup

## -----------------------------------------------------------------------------
# enter the msgf number that is important for this analysis
msgf = 2.7113E-8

# add in the sheet names
tab_names <- data.frame(talias = "T_alias",
                          analysis_jobs = "T_Data_Package_Analysis_Jobs",
                          reporter_ions = "T_Reporter_Ions_Typed",
                          mage = "Mage",
                          protein_collection = "Protein_collection_data",
                          fdata = "f_data")

# add in the column name information
mage_cols <- data.frame(job_name = "Job",
                        scan_name = "Scan",
                        peptide_name = "Peptide",
                        protein_name = "Protein",
                        qvalue_name = "QValue",
                        msgf_specprob_name = "MSGF_SpecProb")

analysis_cols <- data.frame(job_name = "Job",
                            dataset_id_name = "Dataset_ID")

reporter_cols <- data.frame(dataset_id_name = "Dataset",
                            scan_name = "ScanNumber")

fdata_cols <- data.frame(sampleID_name = "SampleID",
                    group_name = "Grouping",
                    job_name = "Job",
                    replicate_name = "Replicate",
                    plex_name = "Plex",
                    ionization_name = "Ionization")

protein_collection_cols <- data.frame(proteinName_name = "protein_name")

normalization_info <- data.frame(reference_name = "Reference",
                                 norm_fn = "mean",
                                 backtransform = TRUE)

## ----eval = FALSE-------------------------------------------------------------
# # put xlsx object into more R friendly formatting
# mydata <- read_excel(here("Data",data))
# sheet_names <- excel_sheets(here("Data",data))
# # this creates a list with each element being a sheet by its corresponding name
# dat_list <- lapply(sheet_names, function(x) {          # Read all sheets to list
#   as.data.frame(read_excel(here("Data",data),sheet = x))})
# names(dat_list) <- sheet_names

## ----echo = FALSE-------------------------------------------------------------
dat_list <- data

## -----------------------------------------------------------------------------
# step 4: run cleaning chemoprot function ######################################
# create pmartObj
htp_pmart_cleaned <- clean_chemoprot_labeled(dat_list,tab_names,mage_cols,analysis_cols,
                                        reporter_cols,fdata_cols,protein_collection_cols)
# save the results and msgf/job information
labeled_information <- list(pmartObj = htp_pmart_cleaned,msgf = msgf,fdata_info = fdata_cols,
                            analysis_info = analysis_cols,mage_info = mage_cols,norm_info = normalization_info,
                            protein_info = protein_collection_cols)

## ----eval = FALSE-------------------------------------------------------------
# saveRDS(labeled_information,paste0(here("Data",data),"labeled_information.RDS"))

## ----eval = FALSE-------------------------------------------------------------
# # run this for loop to generate reports for each job
# for(i in 1:length(htp_pmart_cleaned)){
#   rmarkdown::render(input = here("Code","report_labeled.Rmd"),
#                     output_file = paste0("report_labeled_job_",unique(htp_pmart_cleaned[[i]]$e_meta$Job),".html"))
# }

