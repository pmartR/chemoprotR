# front end script
# damon leach, kelly stratton, gerard lomas
# april 4 2023

# LABELED ######################################################################
# step 1: load in libraries and data ###########################################
# load in necessary libraries
library(pmartR)
library(readxl)
library(chemoprotR)
library(here)

here::i_am("Code/front_end_script_labeled.R")

# step 2: enter in the arguments ###############################################
# labeled
mydata = "60473_Practice_Labeled.xlsx"

# enter the msgf number that is important for this analysis
msgf = 6.76378E-09

# add in the sheet names
tab_names <- data.frame(analysis_jobs = "T_Data_Package_Analysis_Jobs",
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
                                 norm_fn = "median",
                                 backtransform = TRUE)

# rollup_method can take the values: rollup, rrollup, or summation
# centering_fn can take the values: mean, median, (or none if using summation)
rollup_info <- data.frame(rollup_method = "summation",
                          centering_fn = "none")

# list outliers to be removed (recommend leaving this blank for first pass at data
# and then running the report a second time with outliers to remove if there are any)
outlier_samples <- c()

# isobaric information
# can be the value "Yes" or "No"
run_isobaric <- "Yes"

# step 3: create pmart object ##################################################
# put xlsx object into more R friendly formatting
sheet_names <- excel_sheets(here("Data",mydata))
# this creates a list with each element being a sheet by its corresponding name
dat_list <- lapply(sheet_names, function(x) {          # Read all sheets to list
  as.data.frame(read_excel(here("Data",mydata),sheet = x))})
names(dat_list) <- sheet_names
# remove xlsx from data_name
data_name <- stringr::str_remove(mydata,".xlsx")

# step 4: run cleaning chemoprot function ######################################
# create pmartObj
htp_pmart_cleaned <- clean_chemoprot_labeled(dat_list,tab_names,mage_cols,analysis_cols,
                                     reporter_cols,fdata_cols,protein_collection_cols)
# save the results and msgf/job information
labeled_information <- list(pmartObj = htp_pmart_cleaned,msgf = msgf,fdata_info = fdata_cols,
                            analysis_info = analysis_cols,mage_info = mage_cols,norm_info = normalization_info,
                            protein_info = protein_collection_cols,data_name = mydata,
                            rollup_info = rollup_info, outlier_info = outlier_samples,
                            run_isobaric = run_isobaric)
saveRDS(labeled_information,here("Data","labeled_information.RDS"))

# step 5: render markdown and results ##########################################
num_jobs = length(htp_pmart_cleaned)
job_names = names(htp_pmart_cleaned)
# run this for loop to generate reports for each job
for(i in 1:num_jobs){
  rmarkdown::render(input = here("Code","report_labeled.Rmd"),
                    output_file = paste0("report_labeled_",data_name,"_",job_names[i],".html"))

}
