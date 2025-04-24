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
                          fdata = "f_data")

# add in the column name information
mage_cols <- data.frame(job_name = "Job",
                        scan_name = "Scan",
                        peptide_name = "Peptide",
                        protein_name = "Protein")

analysis_cols <- data.frame(job_name = "Job",
                            dataset_id_name = "Dataset_ID",
                            qvalue_name = "QValue",
                            msgf_specprob_name = "MSGF_SpecProb")
reporter_cols <- data.frame(dataset_id_name = "Dataset",
                            scan_name = "ScanNumber")

fdata <- data.frame()

#unlabeled and labeled fdata to send to them

# step 3: create pmart object ##################################################
# put xlsx object into more R friendly formatting
sheet_names <- excel_sheets(paste0(datadir,"/",data))
# this creates a list with each element being a sheet by its corresponding name
dat_list <- lapply(sheet_names, function(x) {          # Read all sheets to list
  as.data.frame(read_excel(paste0(datadir,"/",data), sheet = x)) } )
names(dat_list) <- sheet_names

# step 4: run cleaning chemoprot function ######################################
# create pmartObj
htp_pmart <- clean_chemoprot_labeled(dat_list,tab_names,mage_cols,analysis_cols,
                                     reporter_cols,fdata_cols)
# save the results and msgf/job information
pmart_and_msgf <- list(pmartObj = htp_pmart,msgf = msgf)
saveRDS(pmart_and_msgf,paste0(datadir,"/pmart_and_msgf.RDS"))

# step 5: render markdown and results ##########################################
# run this for loop to generate reports for each job
for(i in 1:length(htp_pmart)){
  rmarkdown::render(input = paste0(here::here(),"/Code/report_labeled.Rmd"),
                    output_file = paste0("report_labeled_job_",unique(htp_pmart$e_meta$Job),".html"))
}

# UNLABELED ####################################################################
# step 1: load in libraries and data ###########################################
# load in necessary libraries
library(pmartR)
library(readxl)
#library(chemoprotR)

# step 2: enter in the arguments ###############################################
# data is the name of xlsx file that we are loading in

# unlabeled
# data is the name of xlsx file that we are loading in
data = "LabelFree_Practice_Dataset.xlsx"

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
                         group_name = "Group",
                         rep_name = "Replicate",
                         job_name = "Job")

# step 3: create pmart object ##################################################
# put xlsx object into more R friendly formatting
sheet_names <- excel_sheets(paste0(datadir,"/",data))
# this creates a list with each element being a sheet by its corresponding name
dat_list <- lapply(sheet_names, function(x) {          # Read all sheets to list
  as.data.frame(read_excel(paste0(datadir,"/",data), sheet = x)) } )
names(dat_list) <- sheet_names

# we also want the protein collection subset information
protein_collection_dat <- dat_list[[tab_names$protein_collection]]

# step 4: run cleaning chemoprot function ######################################
# create pmartObj
htp_pmart <- clean_chemoprot_unlabeled(dat_list,tab_names,mage_cols,
                                       metadata_cols,fdata_cols)
# save the results and msgf/job information
pmart_msgf_protein <- list(pmartObj = htp_pmart,msgf = msgf,prot_collection = protein_collection_dat)
saveRDS(pmart_msgf_protein,paste0(datadir,"/pmart_msgf_protein.RDS"))

# step 5: render markdown and results ##########################################
# run this for loop to generate reports for each job
rmarkdown::render(input = paste0(here::here(),"/Code/report_unlabeled.Rmd"),
                    output_file = paste0("report_unlabeled",".html"))









