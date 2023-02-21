

if (!require('readr')) install.packages('readr'); library('readr')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('hexbin')) install.packages('hexbin'); library('hexbin')
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('ukbtools')) install.packages('ukbtools'); library('ukbtools')
if (!require('h2o')) install.packages('h2o'); library('h2o')




#####
# First read in all the data.
#####




alldata<-tibble(xindex=1:10000)
#loop through all chromosomes
mycolumn=1
for( chr in 1:22){
  for(subfile in 0:3){
  mycolumn<-mycolumn+1
  getfile<-paste("D:/jpbrody/UKB/l2r/chr",chr,"_10k0",subfile, sep="")
  chrtestdata<-read_table2(getfile,col_names = FALSE)
  chrtestdata[ chrtestdata == -1] <- NA                          #according to the spec, value of -1 is NA
  chrcolsums<-colSums(chrtestdata,na.rm=TRUE)/nrow(chrtestdata)  #this is the heart of the "analysis". You could change it
  chrcolsums<-as.data.frame(chrcolsums)
  rm(chrtestdata)                                                #if you don't remove these, R uses lots of memory quickly
  
  alldata<- add_column(alldata,chrcolsums$chrcolsums)
  newname<-paste("chr",chr,"_",subfile,sep="")
  names(alldata)[mycolumn]<-newname
    }
}

#Do the same thing for cancer data
allbcdata<-tibble(xindex=1:1535)                              # I hardcoded this for the number of breast cancer cases, bad form.
mycolumn=1
for( chr in 1:22){
  for(subfile in 0:3){
    mycolumn<-mycolumn+1
    getfile<-paste("D:/jpbrody/UKB/l2r/chr",chr,"_breastcancer0",subfile, sep="")
    chrtestdata<-read_table2(getfile,col_names = FALSE)
    chrtestdata[ chrtestdata == -1] <- NA   #according to the spec, value of -1 is NA
    chrcolsums<-colSums(chrtestdata,na.rm=TRUE)/nrow(chrtestdata) #this is the heart of the "analysis" you could change it
    chrcolsums<-as.data.frame(chrcolsums)
    rm(chrtestdata)
    
    allbcdata<- add_column(allbcdata,chrcolsums$chrcolsums)
    newname<-paste("chr",chr,"_",subfile,sep="")
    names(allbcdata)[mycolumn]<-newname
  }
}

#get the breast cancer mapping information.  Maps rownumber to patientID.
bcmap<-read.csv("D:\\jpbrody\\UKB\\l2r\\breastcancermapping.csv")  # this file is created by makecutscripts.R
allbcdata$PatientIDs<-bcmap$patientID


# Read in the patient ID numbers and line them up with the l2r data
#
# The patientIDSall.fam file is obtained from this command "./ukbgene l2r -c1 -m"
# You can use any chromosome "-cN" they all return the same file.  I renamed the file it returns to "patientIDSall.fam"
#
#
patientIDs <- read_table2("D:/jpbrody/UKB/famfiles/patientIDsall.fam", 
                          col_names = FALSE)


patientIDs$rowno<- 1:nrow(patientIDs)


alldata$PatientIDs<-head(patientIDs$X1,nrow(alldata))
alldata<- select(alldata,PatientIDs,everything())
allbcdata<-select(allbcdata,PatientIDs,everything())


#
#  These just combine the four subparts of chromosomes 1,2 and 3 for some of the correlation plots below.
#
alldata$chr1<-alldata$chr1_0+alldata$chr1_1+alldata$chr1_2+alldata$chr1_3
alldata$chr2<-alldata$chr2_0+alldata$chr2_1+alldata$chr2_2+alldata$chr2_3
alldata$chr3<-alldata$chr3_0+alldata$chr3_1+alldata$chr3_2+alldata$chr3_3

alldata$chr10<-alldata$chr10_0+alldata$chr10_1+alldata$chr10_2+alldata$chr10_3
alldata$chr20<-alldata$chr20_0+alldata$chr20_1+alldata$chr20_2+alldata$chr20_3
ggplot(alldata, aes(x=chr1))+geom_histogram(binwidth=0.005, color="black",fill="white")+xlim(-0.3,0.2)+ylim(0,400)
ggplot(alldata, aes(x=chr2))+geom_histogram(binwidth=0.005, color="black",fill="white")+xlim(-0.3,0.2)+ylim(0,400)

# Some plots to look at correlations in the data
 ggplot(alldata, aes(x=chr15_0, y=chr16_0) ) +
   geom_hex() +    theme_bw()
 
 ggplot(alldata, aes(x=chr1, y=chr20) ) +
   geom_hex() +    theme_bw()
 
 
 ggplot(alldata, aes(x=chr15_0, y=chr16_0) ) +
   geom_point() + geom_density2d()+    theme_bw()
 
 ################
 #  These four are used for a figure in the paper:
 #
 ggplot(alldata, aes(x=chr1, y=chr2) ) +
   geom_point() + geom_hex()+ theme_bw()
 
 ggplot(alldata, aes(x=chr2, y=chr3) ) +
   geom_point() + geom_hex()+ theme_bw()
 
 ggplot(alldata, aes(x=chr1_0, y=chr2_0) ) +
   geom_point() + geom_hex()+  theme_bw()  + ylim(-0.1,0.1) + xlim(-0.1,0.1)

  ggplot(alldata, aes(x=chr1_0, y=chr1_1) ) +
   geom_hex() +    theme_bw() 
#
#  
#########################
  
  
 ggplot(alldata, aes(x=chr2_3, y=chr1_1) ) +
   geom_point() + geom_hex()+ theme_bw()
 

 
 ggplot(alldata, aes(x=chr1_2, y=chr1_3) ) +
   geom_hex() +    theme_bw()

 ggplot(alldata, aes(x=chr15_1, y=chr21_1) ) + geom_point() +
   geom_density2d() +   theme_bw() 
# 
# ggplot(alldata, aes(x=Chr17, y=Chr22) ) +
#   geom_hex() + xlim(-1500,1500) +ylim(-1500,1500)+
#   theme_bw()

#compute correlations

corrdata<-alldata %>% select(contains("Chr"))


# This is a nice table of all the cross correlations between different regions of the genome
GLV_corrs<-round(cor(corrdata),2)






#####
#Read in the ukb data.  This has patientID and male/female cancer status etc.
my_ukb_data <- ukb_df("ukb29183",path="D:/jpbrody/UKB")

my_ukb_data<- select(my_ukb_data,eid,sex_f31_0_0,
                     year_of_birth_f34_0_0,
                     cancer_code_selfreported_f20001_0_0,
                 type_of_cancer_icd9_f40013_0_0,
                 type_of_cancer_icd9_f40013_1_0,
                 behaviour_of_cancer_tumour_f40012_0_0)

####
# Merge the patient UKB data with the l2r data
my_data<-merge(alldata,my_ukb_data, by.x="PatientIDs", by.y="eid")

my_bcdata<-merge(allbcdata,my_ukb_data, by.x="PatientIDs", by.y="eid")



#  This builds up the control data, women who don't have cancer.
#
#
 controldata <- my_data %>%  
   filter(sex_f31_0_0=='Female') %>% 
   filter(is.na(type_of_cancer_icd9_f40013_0_0) &
            is.na(type_of_cancer_icd9_f40013_1_0) &
            is.na(cancer_code_selfreported_f20001_0_0) &
            is.na(behaviour_of_cancer_tumour_f40012_0_0)
          ) %>%  #include only those cases that don't say anything that looks like cancer 
   select(PatientIDs,
          sex=sex_f31_0_0,
          birthyear=year_of_birth_f34_0_0,
          cancer_selfreported=cancer_code_selfreported_f20001_0_0, 
          behavior=behaviour_of_cancer_tumour_f40012_0_0,
          icd9_0=type_of_cancer_icd9_f40013_0_0,
          icd9_1= type_of_cancer_icd9_f40013_1_0,
          contains("Chr")
          )

 controldata<-select(controldata,-chr1,-chr2,-chr3,-chr10,-chr20) # this removes these three columns that I don't need, but I added for a figure/table in the paper.
 
 #
 #  This is the dataset of women with breast cancer.
 #
 #
 cancerdata <- my_bcdata %>%  
   filter(sex_f31_0_0=='Female') %>% 
      select(PatientIDs,
          sex=sex_f31_0_0,
          birthyear=year_of_birth_f34_0_0,
          cancer_selfreported=cancer_code_selfreported_f20001_0_0, 
          behavior=behaviour_of_cancer_tumour_f40012_0_0,
          icd9_0=type_of_cancer_icd9_f40013_0_0,
          icd9_1= type_of_cancer_icd9_f40013_1_0,
          contains("Chr")
   )
 

 # this will return all records that are in cancerdata, but not controldata
 # it should be exactly equal to cancerdata.  If it isn't there are some repeated records
 # in both datasets.  Can't do that.  
leak_check<-anti_join(cancerdata,controldata,by="PatientIDs")   
stopifnot(identical(leak_check,cancerdata)) # some of the cancer data leaked into the control data
 
 
 
 
  
# all the cases of not breast cancer are relabelled as Normal.    
controldata$label<-"Normal"

#all the cases of breast cancer are labelled that way
cancerdata$label<-"BreastCancer"


#Now we combine the two datasets into one for machine learning.
trainingdata<-rbind(controldata,cancerdata)

#Just putting the label in the first column for convenience.
trainingdata<- trainingdata %>% select(label,everything())





# h20 requires unordered factors.    
trainingdata$sex<-factor(trainingdata$sex, ordered=FALSE)
trainingdata$cancer_selfreported<-factor(trainingdata$cancer_selfreported, ordered=FALSE)
trainingdata$icd9_0<-factor(trainingdata$icd9_0,ordered=FALSE)
trainingdata$icd9_1<-factor(trainingdata$icd9_1,ordered=FALSE)
trainingdata$label<-factor(trainingdata$label,ordered = FALSE)
trainingdata$behavior<-factor(trainingdata$behavior,ordered=FALSE)
# train with all data, we'll cross validate
# you could make the split here to  train/test/validate datasets

#I'm going to use 85% of training data to train and withold the rest for testing afterwards
train<-trainingdata %>% sample_frac(0.85)
test <- trainingdata %>% anti_join(train, by="PatientIDs")

h2o.init(nthreads=4)

#load data into h2o
train.hex<- as.h2o(train, destination_frame = "train.hex")  
test.hex<- as.h2o(test, destination_frame = "test.hex")
#
#  I usually stop here and goto http://localhost:54321/flow/index.html 
#  h2o runs a local webserver on port 54321, it offers a nice little interface.
#  you can run the AutoML from the web browser there and it has some nice features so that 
# you can monitor the progress of the training.
#

#
#
#
model <- h2o.automl(x=9:ncol(train),
                 y=1,
                 training_frame = train.hex,
                 nfolds=5,
                 max_runtime_secs = 40000,
                 keep_cross_validation_predictions = TRUE)

winner<- model@leader
#record the AUC in the dataset
auc=h2o.auc(winner, train=FALSE, xval=TRUE)



# plot out the ROC.  We type out the tissue and AUC at the top of the ROC.
plot(h2o.performance(winner,train=FALSE, xval=TRUE),type='roc',main='Breast Cancer')

#ROC with ggplot
perf_automl <- h2o.performance(winner, newdata = test.hex)
automl_curve_dat <- data.frame(perf_automl@metrics$thresholds_and_metric_scores) %>%
  select(c(tpr, fpr))
ggplot() +
  geom_line(data=automl_curve_dat,aes(x=fpr, y=tpr), color="blue", size=3, linetype="solid") +
  geom_segment(
    aes(x = 0, y = 0, xend = 1, yend = 1),
    linetype = "dotted",
    color = "grey50"
  ) +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  ggtitle("ROC Curve AUC= ",perf_automl@metrics$AUC) +
  theme_bw()

h2o.explain(model,newdata = test.hex)

#print the entire leaderboard
lb <- model@leaderboard
print(lb, n = nrow(lb))

# get the predictions on the test set
# 
test_predictions <- as.data.frame(h2o.predict(model, test.hex))

#combine with the original data
predictions<-cbind(test_predictions,test)

#rank by tissue predicted score ascending
predictions$rank<- rank(-predictions[,"BreastCancer"])
predictions<- select(predictions, rank, label, BreastCancer, Normal, PatientIDs, birthyear)


# This builds a dataframe for percentile of genetic risk score vs odds ratio

 num_cancers<-nrow(filter(predictions,label=="BreastCancer"))
 total_samples<-nrow(predictions)



ordataframe <- predictions %>% 
  mutate(decile = ntile(BreastCancer, 10)) %>% 
  group_by(decile,label) %>% 
  tally() %>%
  pivot_wider(names_from = label, values_from = n)

ordataframe$OR <- (ordataframe$BreastCancer/ordataframe$Normal)/(num_cancers/(total_samples-num_cancers))

# These compute the top/bottom of the 95% CI
ordataframe$tci<-exp(log(ordataframe$OR)+
                       1.96*sqrt(
                         (1/(ordataframe$BreastCancer+1))+
                           (1/(ordataframe$Normal+1))+
                           (1/sum(ordataframe$BreastCancer))+
                           (1/sum(ordataframe$Normal))
                       ))
ordataframe$bci<-exp(log(ordataframe$OR)-
                       1.96*sqrt(
                         (1/(ordataframe$BreastCancer+1))+
                           (1/(ordataframe$Normal+1))+
                           (1/sum(ordataframe$BreastCancer))+
                           (1/sum(ordataframe$Normal))
                       ))

# 8/31/2021                     top model was stacked ensemble  all 0.765
# after automl 3600 seconds     top non-stacked model was GBM 0.69

# Reran automl 36000 seconds    top model was stacked ensemble all auc=0.837
#                               top model non-stacked was Deeplearning 0.778

#automl  44000 seconds top model stacked ensemble all auc=0.828
#                      top model non-stacked Deeplearning0.767

#automl  88000 seconds top model stacked ensemble all 0.841
#                      top model non-stacked Deeplearning 0.799

#automl  88000 seconds top model stacked ensemble all 0.846
#                      top model non-stacked Deeplearning 0.811

# automl 40000 seconds top model stacked ensemble all 0.830
#automl               top model non-stacked deeplearning 0.788