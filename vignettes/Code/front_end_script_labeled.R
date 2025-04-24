# front end script
# damon leach, kelly stratton, gerard lomas
# january 3, 2023

# LABELED ######################################################################
# step 1: load in libraries and data ###########################################
# load in necessary libraries
library(pmartR)
library(readxl)
#library(chemoprotR)

# find the data directory
datadir = paste0(here::here(),"/Data")

# step 2: enter in the arguments ###############################################
# labeled
# data is the name of xlsx file that we are loading in
data = "TMT_Practice_Dataset.xlsx"

# enter the msgf number that is important for this analysis
msgf = 2.7113E-8

# add in the sheet names
tab_names <- data.frame(talias = "T_alias",
                          tdata = "T_Data",
                          analysis_jobs = "T_Data_Package_Analysis_Jobs",
                          reporter_ions = "T_Reporter_Ions_Typed",
                          mage = "Mage",
                          sheet1 = "Sheet1",
                          fdata = "f_data",
                        sheet7 = "Sheet7")

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

normalization_info <- data.frame(reference_name = "Reference",
                                 norm_fn = "mean",
                                 backtransform = TRUE)

# step 3: create pmart object ##################################################
# put xlsx object into more R friendly formatting
sheet_names <- excel_sheets(paste0(datadir,"/",data))
# this creates a list with each element being a sheet by its corresponding name
dat_list <- lapply(sheet_names, function(x) {          # Read all sheets to list
  as.data.frame(read_excel(paste0(datadir,"/",data), sheet = x)) } )
names(dat_list) <- sheet_names

# step 4: run cleaning chemoprot function ######################################
# create pmartObj
htp_pmart_cleaned <- clean_chemoprot_labeled(dat_list,tab_names,mage_cols,analysis_cols,
                                     reporter_cols,fdata_cols,normalization_info)
# save the results and msgf/job information
labeled_information <- list(pmartObj = htp_pmart_cleaned,msgf = msgf,fdata_info = fdata_cols,
                            analysis_info = analysis_cols,mage_info = mage_cols,norm_info = normalization_info)
saveRDS(labeled_information,paste0(datadir,"/labeled_information.RDS"))

# step 5: render markdown and results ##########################################
# run this for loop to generate reports for each job
for(i in 1:length(htp_pmart_cleaned)){
  rmarkdown::render(input = paste0(here::here(),"/Code/report_labeled.Rmd"),
                    output_file = paste0("report_labeled_job_",unique(htp_pmart_cleaned[[i]]$e_meta$Job),".html"))
}