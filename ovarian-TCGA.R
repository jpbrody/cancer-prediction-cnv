#
#
# This is the Google cloud project ID.  Change it to *your* project ID, or it won't work
# When this program makes the call to Bigquery to retrieve the TCGA data, it will prompt you to 
# login to your google account.

project<- "test-project-cgc"

# Google charges for Bigquery access, but gives you 1TB of queries free per month.
# The query here will consume about 250 MB, or 1/4000 th of your monthly quota.



library("dplyr")
library("h2o")
library("bigrquery")
library("ggplot2")
library("tidyr")
library("openssl")


topcnvs<-25  # number of features to use in the model

#start the h2o instance
h2o.init(nthreads=2)




#list of  all  TCGA cancers
#tissues<-c("UCEC","BLCA","PRAD","BRCA","OV","SARC","GBM","SKCM","HNSC","PAAD","LUSC","KIRP","LGG","LUAD","STAD","THCA","LIHC","COAD","CESC","PCPG","MESO","ESCA","READ","TGCT","KIRC","THYM","UVM","ACC","UCS","DLBC","CHOL","KICH")
tissues<-c("OV")
numbers<-seq(1,5) # repeat 5 times
#create a dataframe for results
heritability_results <- data.frame(matrix(ncol=length(tissues), nrow = length(numbers)))
tissue_auc <- data.frame(matrix(ncol=length(tissues), nrow = length(numbers)))
#label the dataframe
colnames(tissue_auc)<-tissues
rownames(tissue_auc)<-numbers


BigQuerySQL_start_end <- paste("#standardSQL
                     with topCNVs as",
                         "(SELECT",
                         "chromosome AS chromosome,",
                         "start_pos,",
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
                         "start_pos,",
                         "end_pos",
                         "ORDER BY",
                         "total DESC",
                         "limit ", topcnvs,
                         ")",
                         "SELECT",
                         "cn.project_short_name as project_name,",
                         "bs.sample_barcode,",
                         "clin.gender as gender,",
                         "clin.age_at_diagnosis as age,",
                         "clin.race as race,",
                         "clin.ethnicity as ethnicity,",
                         "clin.days_to_death as days_to_death,",
                         "cn.chromosome as chromosome,",
                         "cn.start_pos as start_pos,",
                         "cn.end_pos as end_pos,",
                         "avg(segment_mean) as segment_mean",
                         "FROM",
                         "`isb-cgc.TCGA_hg38_data_v0.Copy_Number_Segment_Masked` cn",
                         "JOIN",
                         "`isb-cgc.TCGA_bioclin_v0.Biospecimen` bs",
                         "ON",
                         "cn.sample_barcode = bs.sample_barcode",
                         "JOIN",
                         #"topCNVs top ON (top.chromosome=cn.chromosome and top.end_pos=cn.end_pos )",
                         "topCNVs top ON (top.chromosome=cn.chromosome and top.start_pos=cn.start_pos and top.end_pos=cn.end_pos )",
                         "JOIN",
                         "`isb-cgc.TCGA_bioclin_v0.Clinical` clin",
                         "ON",
                         "bs.case_barcode = clin.case_barcode",
                         "WHERE",
                         "bs.sample_type ='10'",
                         "group by project_name,bs.sample_barcode,gender,race,ethnicity,age,days_to_death,chromosome,start_pos, end_pos"
                         #"group by project_name,bs.sample_barcode,gender,age,days_to_death,chromosome,end_pos"
)


#Get the data
copy_number_start_end<-query_exec(BigQuerySQL_start_end,project,max_pages=Inf, use_legacy_sql = FALSE)

#manipulate it into the proper form
copy_number1_start_end <- unite(copy_number_start_end, "chromosome_start_end", c("chromosome","start_pos", "end_pos"))


for (j in 1:length(numbers)) 
{
   
  copy_number2_start_end <- spread(copy_number1_start_end,chromosome_start_end,segment_mean)

  # Starting simple, just start_end
  copy_number2<-copy_number2_start_end

  i<-1  
  tissue<-tissues[i]

  copy_number2$tissue<-ifelse((copy_number2$project_name == paste('TCGA-',tissue,sep="")),tissue,'Normal')
    
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
    
    
    #
    tissue_copynumber$tissue<-as.factor(tissue_copynumber$tissue)
    tissue_copynumber$gender<-as.factor(tissue_copynumber$gender)
    
    #remove a bunch of columns we don't need
     tissue_copynumber <- select(tissue_copynumber,-project_name,-sample_barcode,-days_to_death,-age,-race,-ethnicity)
   
    # train with all data, we'll cross validate
    train <- tissue_copynumber 
    
    #load data into h2o
    train.hex<- as.h2o(train, destination_frame = "train.hex")  
    
    
    model <- h2o.gbm(x=3:ncol(train),
                     y=1,
                     training_frame = train.hex,
                     nfolds=10,
                     ntrees=100,
                     keep_cross_validation_predictions = TRUE)   
    
    #get the predictions
    cvpreds<-as.data.frame(h2o.getFrame(model@model[["cross_validation_holdout_predictions_frame_id"]][["name"]]))
    
    #combine with the original data
    predictions<-cbind(cvpreds,train)
    
    #rank by tissue predicted score ascending
    predictions$rank<-rank(-predictions[,tissue])
    predictions<- select(predictions, rank, everything())
    
   #total number of OV cancers
    num_cancers <- nrow(predictions[predictions$tissue==tissue,])
    total_samples <- nrow(predictions)
    heritability<-nrow(subset(predictions[predictions$tissue==tissue,], rank>0.5*total_samples))
    
    heritability_results[j,tissue]=1-(heritability*2/num_cancers)
    
    message(paste(tissue, heritability_results[j,tissue], j))
    
    }

    
num_cancers<-789
total_samples<-1578

# This builds a dataframe for percentile of genetic risk score vs odds ratio
ordataframe<- predictions %>% mutate(decile = ntile(OV, 5)) %>% group_by(decile) %>% 
  summarize(Normal=summary(tissue)[["Normal"]], 
            OV=summary(tissue)[["OV"]], 
            OR=(summary(tissue)[["OV"]]/summary(tissue)[["Normal"]])/(num_cancers/(total_samples-num_cancers)),
            count = n())
#
# Compute 95% confidence intervals (TCI,BCI) top of ci and bottom of ci
# based on statperls https://www.ncbi.nlm.nih.gov/books/NBK431098/
ordataframe$tci<-exp(log(ordataframe$OR)+
                       1.96*sqrt(
                         (1/(ordataframe$OV+1))+
                           (1/(ordataframe$Normal+1))+
                           (1/sum(ordataframe$OV))+
                           (1/sum(ordataframe$Normal))
                       ))
ordataframe$bci<-exp(log(ordataframe$OR)-
                       1.96*sqrt(
                         (1/(ordataframe$OV+1))+
                           (1/(ordataframe$Normal+1))+
                           (1/sum(ordataframe$OV))+
                           (1/sum(ordataframe$Normal))
                       ))

# These lines for building the graph
# ordataframe$xaxis<-sequence(50)*2
# ggplot(ordataframe,aes(x=xaxis,y=OR))+geom_point()+geom_errorbar(aes(ymin=bci, ymax=tci), width=.2,position=position_dodge(0.05)) 





tissue_auc_excel<-t(tissue_auc)

#plot the results
# library(ggthemes)
# tissue_auc$num<-seq.int(nrow(tissue_auc))
# tissue_auc_long <- gather(tissue_auc,tissue,value,1:4)
# theme_wsj()
# ggplot(tissue_auc_long, aes(x=num,y=value))+geom_point(aes(color=tissue))+geom_smooth(aes(color=tissue))+
#   labs(x="Number of sites incoporated into the model",y="Area Under the ROC curve (AUC)", colour="Type of cancer")+theme_minimal()


#plots for paper
#histogram of chr1:
ggplot(copy_number2_start_end,aes(x=`1_3301765_247650984`))+geom_histogram(binwidth=0.0005)+xlim(-0.04,0.04)
ggplot(copy_number2_start_end,aes(x=`17_1074619_82959812`))+geom_histogram(binwidth=0.0005)+xlim(-0.04,0.04)
ggplot(copy_number2_start_end,aes(x=`6_1011760_170596889`))+geom_histogram(binwidth=0.0005)+xlim(-0.04,0.04)
ggplot(copy_number2_start_end,aes(x=`13_18874255_114226675`))+geom_histogram(binwidth=0.0005)+xlim(-0.04,0.04)


#ROC plot
# Model performance 1 time split data 80/20 train/test
smp_size <- floor(0.80 * nrow(train))
train_ind <- sample(seq_len(nrow(train)), size = smp_size)

trainroc <- train[train_ind, ]
testroc <- train[-train_ind, ]

trainroc.hex<- as.h2o(trainroc, destination_frame = "trainroc.hex")  
testroc.hex<- as.h2o(testroc, destination_frame = "testroc.hex")  

# This bit sets the positive class to be cancer:
trainroc.hex["tissue"] <- h2o.relevel(x = trainroc.hex["tissue"], y = "OV")
testroc.hex["tissue"] <- h2o.relevel(x = testroc.hex["tissue"], y = "OV")

model <- h2o.gbm(x=3:ncol(trainroc),
                 y=1,
                 training_frame = trainroc.hex,
                 nfolds=5,
                 ntrees=100)   
perf <- h2o.performance(model, newdata = testroc.hex)




# Extract info for ROC curve
curve_dat <- data.frame(perf@metrics$thresholds_and_metric_scores) %>%
  select(c(tpr, fpr))

# Plot ROC curve
ggplot(curve_dat, aes(x = fpr, y = tpr)) +
  geom_point() +
  geom_line() +
  geom_segment(
    aes(x = 0, y = 0, xend = 1, yend = 1),
    linetype = "dotted",
    color = "grey50"
  ) +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  ggtitle("ROC Curve AUC= ",perf@metrics$AUC) +
  theme_bw()