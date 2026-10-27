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
library(MSnSet.utils)

here::i_am("Code/front_end_script_fragpipe.R")

# step 2: enter in the arguments ###############################################
study_design = "study_design_adj.csv"
# study_design = "study_design.csv"
mydata = "ratio_gene_None.tsv"

# # normalization preferences
# # norm_fn can take the following values: "mean" or "median"
# # backtransform can take the following values: TRUE or FALSE
normalization_info <- data.frame(norm_fn = "mean",
                                 backtransform = TRUE)

# # rollup preferences
# # rollup_method can take the values: rollup, rrollup, or summation
# # centering_fn can take the values: mean, median, (or none if using summation)
# rollup_info <- data.frame(rollup_method = "summation",
#                           centering_fn = "none")

# list outliers to be removed (recommend leaving this blank for first pass at data
# and then running the report a second time with outliers to remove if there are any)
outlier_samples <- c()

# step 3: create pmart object ##################################################
# put xlsx object into more R friendly formatting
fragpipe_obj <- read_FragPipe_TMT(here("Data",mydata))
# this creates a list with each element being a sheet by its corresponding name
study_df <- read.csv(here("Data",study_design))

# step 4: run cleaning chemoprot function ######################################
# create pmartObj
htp_pmart_cleaned <- clean_chemoprot_fragpipe(m = fragpipe_obj,
                                              study_df = study_df)
# save the results and msgf/job information
fragpipe_information <- list(pmartObj = htp_pmart_cleaned,
                             norm_info = normalization_info,
                             data_name = mydata,
                             rollup_info = rollup_info,
                             outlier_info = outlier_samples)
saveRDS(fragpipe_information,here("Data","fragpipe_information.RDS"))

# step 5: render markdown and results ##########################################
num_jobs = length(htp_pmart_cleaned)
job_names = names(htp_pmart_cleaned)
# run this for loop to generate reports for each jobx
for(i in 1:num_jobs){
  rmarkdown::render(input = here("Code","report_fragpipe.Rmd"),
                    output_file = paste0("report_fragpipe_",job_names[i],".html"))
}
