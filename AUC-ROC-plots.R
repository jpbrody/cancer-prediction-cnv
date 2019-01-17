#
# This program generates the ROCurves and computes AUC for cross validation predictions for different cancer types
# using TCGA data.
#  TCGA data is read in through a call to Google's bigquery. 
#  ISB (seattle) hosts a copy of the relevant TCGA data at https://bigquery.cloud.google.com/dataset/isb-cgc:TCGA_hg38_data_v0

# J. Brody jpbrody@uci.edu April 2018, https://www.biorxiv.org/content/early/2018/04/17/303339
# Released under GPL-3.0 license
#####################################################################################



#
# This will check for needed packages.  If not installed, it will download and install them.
#
list.of.packages <- c("dplyr","h2o","bigrquery","ggplot2","tidyr","httpuv","readr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


# Load packages
library("dplyr")
library("h2o")
library("bigrquery")
library("ggplot2")
library("tidyr")
library("httpuv")
library("readr")




# This variable 'project' is the Google cloud project ID.  

project<- "test-project-cgc"
# You need to change it to *your* Google cloud project ID
# if you don't have one, go to http://bigquery.cloud.google.com/ and make one
# or see this page for help https://support.google.com/cloud/answer/6158840?hl=en
# The project ID is not secret.  You will need to authenticate yourself to google with your secret password.
# Google provides 1 TB of free queries per month.  This program only has one query that will consume less than 1 GB,
# so, it will use less than 1/1000th of your monthly quota.

# the first time you execute a Bigquery, a browser window will open and ask you to authenticate yourself with Google. 
# the authentication token will be cached, so you don't need to do it every time.



topcnvs<-30  # number of CNVs to use in the model

#start the h2o instance
h2o.init(nthreads=4)

# h2o runs a local web server that you can check out.  
# goto http://localhost:54321 you can access and manipulate datafiles there
# However, you don't have to access it to properly run this program.  
#  Just FYI.


#list of  all  TCGA cancers
tissues<-c("UCEC","BLCA","PRAD","BRCA","OV","SARC","GBM","SKCM","HNSC","PAAD","LUSC","KIRP","LGG","LUAD","STAD","THCA","LIHC","COAD","CESC","PCPG","MESO","ESCA","READ","TGCT","KIRC","THYM","UVM","ACC","UCS","DLBC","CHOL","KICH")
#tissues<-c("BRCA","OV","GBM","COAD")
numbers<-seq(1,1) # repeat x times

#create a dataframe for results
tissue_auc <- data.frame(matrix(ncol=length(tissues), nrow = length(numbers)))


#label the dataframe
colnames(tissue_auc)<-tissues
rownames(tissue_auc)<-numbers


# construct the query
# This query grabs the appropriate data from the TCGA data set hosted on bigquery.
#
#  The query first generates a temporary table "topCNVs" of the top copy number variations ranked by how many are recorded..
#  A copy number variation is identified by a chromosome, start_pos, and end_pos.  It has a value segment_mean.
#  This query was initially "SELECT chromosome, start_pos, end_pos, count(*) . . ."
#  But I found that there were several large copy number variations that had the same endpoint, and that
# grouping these by the endpoint gave substantially better results.  So the query now looks like:
#  "SELECT chromosome, end_pos, count(*) . . .


BigQuerySQL <- paste("#standardSQL
                     with topCNVs as",
# First build the topCNVS temp table                     
                     "(SELECT",
                     "chromosome AS chromosome,",
                     "end_pos,",
                     "COUNT(*) AS total",
                     "FROM",
                     "`isb-cgc.TCGA_hg38_data_v0.Copy_Number_Segment_Masked` cn",
                     "JOIN",
                     "`isb-cgc.TCGA_bioclin_v0.Biospecimen` bs",
                     "ON",
                     "cn.sample_barcode = bs.sample_barcode",
                     "WHERE",
                     "bs.sample_type ='10'",
                     "GROUP BY",
                     "chromosome,",
                     "end_pos",
                     "ORDER BY",
                     "total DESC",
                     "limit ", topcnvs,
                     ")",
# now select the data we want
                     "SELECT",
                     "cn.project_short_name as project_name,",
                     "bs.sample_barcode,",
                     "clin.gender as gender,",
                     "clin.age_at_diagnosis as age,",
                     "clin.days_to_death as days_to_death,",
                     "cn.chromosome as chromosome,",
                     "cn.end_pos as end_pos,",
                     "avg(segment_mean) as segment_mean",
                     "FROM",
                     "`isb-cgc.TCGA_hg38_data_v0.Copy_Number_Segment_Masked` cn",
                     "JOIN",
                     "`isb-cgc.TCGA_bioclin_v0.Biospecimen` bs",
                     "ON",
                     "cn.sample_barcode = bs.sample_barcode",
                     "JOIN",
                     "topCNVs top ON (top.chromosome=cn.chromosome and top.end_pos=cn.end_pos )",
                     "JOIN",
                     "`isb-cgc.TCGA_bioclin_v0.Clinical` clin",
                     "ON",
                     "bs.case_barcode = clin.case_barcode",
                     "WHERE",
                     "bs.sample_type ='10'",
                     "group by project_name,bs.sample_barcode,gender,age,days_to_death,chromosome,end_pos"
)


#Now that we built the query, Get the data
copy_number<-query_exec(BigQuerySQL,project,max_pages=Inf, use_legacy_sql = FALSE)

#this line combines the CNV identifier into something like 12_8748575.
# we will  use this identifier as a column header when we transform the data from 
# long to wide format.
copy_number1 <- unite(copy_number, "chromosome_end", c("chromosome","end_pos"))



#loop through all the different types of cancers, called tissues
for (i in 1:length(tissues)){
    tissue<-tissues[i]
    
# this line converts from long to wide format    
    copy_number2 <- spread(copy_number1,chromosome_end,segment_mean)

# all the cases of not $tissue cancer are relabelled as Normal.    
    copy_number2$tissue<-ifelse((copy_number2$project_name == paste('TCGA-',tissue,sep="")),tissue,'Normal')
    

# For gender specific cancers, include only that gender in the normal/control group.    
    if ((tissue=="BRCA") | (tissue=="OV") | (tissue=="CESC") | (tissue=="UCEC") | (tissue=="UCS")){
      tissue_copynumber <- copy_number2 %>%  filter(gender=='FEMALE') %>% select(project_name,tissue,everything())
    } else
    { 
      if ((tissue=="PRAD") | (tissue=="TGCT") ) {
        tissue_copynumber <- copy_number2 %>%  filter(gender=='MALE') %>% select(project_name,tissue,everything()) 
      }
      else {
        tissue_copynumber <- copy_number2 %>%   select(project_name,tissue,everything())
      }
    }
    
# h20 requires the target of the training to be a factor.    
    tissue_copynumber$tissue<-as.factor(tissue_copynumber$tissue)
    tissue_copynumber$gender<-as.factor(tissue_copynumber$gender)
    tissue_copynumber$project_name<-as.factor(tissue_copynumber$project_name)
    #tissue_copynumber <- select(tissue_copynumber,-project_name,-sample_barcode,-days_to_death,-age)
    tissue_copynumber <- select(tissue_copynumber,-project_name,-sample_barcode,-days_to_death)

    # train with all data, we'll cross validate
    # you could make the split here to  train/test/validate datasets
    train <- tissue_copynumber 
    
    #load data into h2o
    train.hex<- as.h2o(train, destination_frame = "train.hex")  
    
    
    model <- h2o.gbm(x=2:ncol(train),
                     y=1,
                     training_frame = train.hex,
                     nfolds=10,
                     ntrees=100)

########################################################################
#  You can use h2o to generate predictions with this bit of code.
#
    # model <- h2o.gbm(x=2:ncol(train),
    #                  y=1,
    #                  training_frame = train.hex,
    #                  nfolds=10,
    #                  ntrees=100,
    #                  keep_cross_validation_predictions = TRUE)   
    # 
    # #get the predictions
    # cvpreds<-as.data.frame(h2o.getFrame(model@model[["cross_validation_holdout_predictions_frame_id"]][["name"]]))
    # 
    # #combine with the original data
    # predictions<-cbind(cvpreds,train)
    # 
#
#    The dataframe predictions will now have columns that contains the prediction, the patient barcode, age, etc
##########################################################################
    
#record the AUC in the dataset
    tissue_auc$tissue=h2o.auc(model, train=FALSE, xval=TRUE)
   
# plot out the ROC.  We type out the tissue and AUC at the top of the ROC.
    plot(h2o.performance(model,train=FALSE, xval=TRUE),type='roc',main=paste(tissue, tissue_auc$tissue))
   
# Also print the tissue and AUC out to the console.    
    message(paste(tissue, tissue_auc$tissue))
}
  




