library(readxl)
library(dada2)
library(DECIPHER)
library(phangorn)
library(phyloseq)
library(gridExtra)
library(dplyr)
library(plyr)                                     #Check dada2 sequences obtained from otu table
library(ggplot2)
library(tidyverse)
library(csv)
library(decontam)
library(metagenomeSeq)
library(splinectomeR)
library(reshape2)
library(fastDummies)    ###########to create dummy variables#####################
library(tibble)
library(tidyr)
library(vegan)
library(gdata) ##########merging columns into one
library(igraph)
library(HMP)
library(dendextend)
library(ppcor)     #########spearman partial correlation coefficient
library(rbin)
library(stats)
library(ggpubr)
library(nloptr)
library(coda)
library(gss)
library(remotes)
library(devtools)
library(SpiecEasi)
library(randomForest)
library(caret)
library(RColorBrewer) ############color palatte
library(plotrix)
library(robCompositions)
library(Maaslin2)
library(psych)    ##################for correlation test###################
library(rmcorr)
library(clr)
library(taxonomizr)
library(gplots) 
library(compositions) 
###########REMOVE OUTLIERS USING Z-SCORE METHOD##############################

df<-read.csv('C:/Akanksha/LAB WORK/Cytokines/cytokines_UNC_04102021.csv',header = TRUE, stringsAsFactors = FALSE)
#df<- na.omit(df)

df<-df[1:60,]

x<-df[,19:39]
x_new<- sapply(x, function(x) replace(x, is.na(x), mean(x, na.rm = TRUE)))

#Checking normality of cytokines concentration

x_log<-x+1
x_log<-log(x_log)
x_log<- lapply(x_log, function(x) replace(x, is.na(x), mean(x, na.rm = TRUE)))
x_log<-as.data.frame(x_log)

y<-df$GAD_7
visit_all<-df$Visit
BMI<-df$BMI

cy_GAD<-cbind.data.frame(x_log,y,visit_all)

write.csv(cy_GAD,
          file="Cytokinevalues_log.csv")
################################################################################################################

#NEED TO GENERATE HEATMAP BETWEEN INDIVIDUAL CYTOKINES ABUNDANCES#################

#Correlation matrix for cytokines using normal values;if missing replaced by mean

cy_corr<-corr.test(x_new,method = "spearman")

corr_values<-cy_corr$r

#melted_corr_values <- melt(corr_values)
#melted_corr_values<-as.matrix(melted_corr_values)
#library(ggplot2)
#ggplot(data = melted_corr_values, aes(x=Var1, y=Var2, fill=value)) + 
#  geom_tile(color = "white")+scale_fill_gradient2(low = "chocolate", high = "Darkgreen", mid = "white", 
 #                                                 midpoint = 0, limit = c(-1,1), space = "Lab", 
#                                                  name="Spearman\nCorrelation")



heatmap.2(corr_values,col =colorRampPalette(c("Darkgreen","white","chocolate")),
          density.info="none",key =TRUE,trace="none",cexCol = 1,cexRow = 1)

#pheatmap(corr_values)
#library(pheatmap)
###############################################################################################################################################
#For Visit V1

df_new_1<-df[(df$Visit=="V1"),]

x_v1<-df_new_1[,19:39]
x_v1_log<-x_v1+1
x_v1_log<-log(x_v1_log)
#Replacing mean in log values
x_v1_log[] <- lapply(x_v1_log, function(x) replace(x, is.na(x), mean(x, na.rm = TRUE)))
y_v1<-df_new_1$GAD_7

age_v1<-df_new_1$Age
bmi_v1<-df_new_1$BMI
education_v1<-df_new_1$Years.of.Education

df_new_1$WhitesVsNOn_whites<-rep("Non-White", nrow(df_new_1))
df_new_1$WhitesVsNOn_whites[(df_new_1$Race==0)]<-"White"
df_new_1_Race<-df_new_1[!(df_new_1$SampleName=="S_701.521.MPI003" | df_new_1$SampleName=="S_721.506.MPI005"),]
WhitesVsNOn_white_v1<-df_new_1_Race$WhitesVsNOn_whites
WhitesVsNOn_white_v1_dummy<-dummy_cols(WhitesVsNOn_white_v1)
WhitesVsNOn_white_v1_dummy<-subset (WhitesVsNOn_white_v1_dummy, select = -c(.data))
WhitesVsNOn_white_v1_dummy<-matrix(unlist(WhitesVsNOn_white_v1_dummy),ncol = 2)



x_v1_Race<-df_new_1_Race[,19:39]
x_v1_log_Race<-x_v1_Race+1
x_v1_log_Race<-log(x_v1_log_Race)
#Replacing mean in log values
x_v1_log_Race[] <- lapply(x_v1_log_Race, function(x) replace(x, is.na(x), mean(x, na.rm = TRUE)))




n1<-21                   ###########use sapply instead of lapply###########(coursera course)

#Tried statistically adjusting(correcting) but still no cytokine is significant with respect to GAD-7 score.

my_lms_v1 <- lapply(1:n1, function(x) lm(x_v1_log[,x]~y_v1+bmi_v1))
coeff_v1<-sapply(my_lms_v1, coef)
summaries_v1 <- lapply(my_lms_v1, summary)
p_v1<-lapply(summaries_v1, function(x) x$coefficients[2,4])
p_v1_matrix <- matrix(unlist(p_v1), ncol = 1)
p_v1_round<-round(p_v1_matrix,digits=4)
p_v1_matrix <-cbind(colnames(x_v1_log),p_v1_round)
p_v1_matrix<-data.frame(p_v1_matrix)
colnames(p_v1_matrix)<-c("names","pvalue")
p_v1_filt<-p_v1_matrix[which(p_v1_matrix[,"pvalue"]<"0.1"),]

my_lms_age_v1<-lapply(1:n1, function(x) lm(x_v1_log[,x] ~ age_v1))
coeff_v1_age<-sapply(my_lms_age_v1, coef)
summaries_v1_age <- lapply(my_lms_age_v1, summary)
p_v1_age<-lapply(summaries_v1_age, function(x) x$coefficients[2,4])
p_v1_age_matrix <- matrix(unlist(p_v1_age), ncol = 1)
p_v1_age_round<-round(p_v1_age_matrix,digits=4)
p_v1_age_matrix <-cbind(colnames(x_v1_log),p_v1_age_round)
p_v1_age_matrix<-data.frame(p_v1_age_matrix)
colnames(p_v1_age_matrix)<-c("names","pvalue")
p_v1_age_filt<-p_v1_age_matrix[which(p_v1_age_matrix[,"pvalue"]<"0.05"),]

my_lms_education_v1<-lapply(1:n1, function(x) lm(x_v1_log[,x] ~ education_v1))
coeff_v1_education<-sapply(my_lms_education_v1, coef)
summaries_v1_education <- lapply(my_lms_education_v1, summary)
p_v1_education<-lapply(summaries_v1_education, function(x) x$coefficients[2,4])
p_v1_education_matrix <- matrix(unlist(p_v1_education), ncol = 1)
p_v1_education_round<-round(p_v1_education_matrix,digits=4)
p_v1_education_matrix <-cbind(colnames(x_v1_log),p_v1_education_round)
p_v1_education_matrix<-data.frame(p_v1_education_matrix)
colnames(p_v1_education_matrix)<-c("names","pvalue")
p_v1_education_filt<-p_v1_education_matrix[which(p_v1_education_matrix[,"pvalue"]<"0.05"),]

my_lms_race_v1<-lapply(1:n1, function(x) lm(x_v1_log_Race[,x] ~ WhitesVsNOn_white_v1_dummy))
coeff_v1_race<-sapply(my_lms_race_v1, coef)
summaries_v1_race<- lapply(my_lms_race_v1, summary)
p_v1_race<-lapply(summaries_v1_race, function(x) x$coefficients[2,4])
p_v1_race_matrix <- matrix(unlist(p_v1_race), ncol = 1)
p_v1_race_round<-round(p_v1_race_matrix,digits=4)
p_v1_race_matrix <-cbind(colnames(x_v1_log),p_v1_race_round)
p_v1_race_matrix<-data.frame(p_v1_race_matrix)
colnames(p_v1_race_matrix)<-c("names","pvalue")
p_v1_race_filt<-p_v1_race_matrix[which(p_v1_race_matrix[,"pvalue"]<"0.05"),]


#####Race ----- IL10--------------------------------

#################################################################################################################################################################################################################
#For visit V2

df_new_2<-df[(df$Visit=="V2"),]
#df_new_2<-df_new_2 %>%select_if(~ !any(is.na(.)))
#df_v2<-df_new[(df$Visit=="V2"),]
#df_v2<-df_v2 %>%select_if(~ !any(is.na(.)))

x_v2<-df_new_2 [,19:39]
x_v2_log<-x_v2+1
x_v2_log<-log(x_v2_log)
x_v2_log[]<- lapply(x_v2_log, function(x) replace(x, is.na(x), mean(x, na.rm = TRUE)))

y_v2<-df_new_2$GAD_7

y_v2_log<-log(y_v2+1)
age_v2<-df_new_2$Age
bmi_v2<-df_new_2$BMI
education_v2<-df_new_2$Years.of.Education

df_new_2$WhitesVsNOn_whites<-rep("Non-White", nrow(df_new_2))
df_new_2$WhitesVsNOn_whites[(df_new_2$Race==0)]<-"White"
WhitesVsNOn_white_v2<-df_new_2$WhitesVsNOn_whites
WhitesVsNOn_white_v2_dummy<-dummy_cols(WhitesVsNOn_white_v2)
WhitesVsNOn_white_v2_dummy<-subset (WhitesVsNOn_white_v2_dummy, select = -c(.data))
WhitesVsNOn_white_v2_dummy<-matrix(unlist(WhitesVsNOn_white_v2_dummy),ncol = 2)




n2<-21

my_lms_v2<- lapply(1:n2, function(x) lm(x_v2_log[,x] ~ y_v2+bmi_v2))
coeff_v2<-sapply(my_lms_v2, coef)
summaries_v2<- lapply(my_lms_v2, summary)
p_v2<-lapply(summaries_v2, function(x) x$coefficients[2,4])
p_v2_matrix<- matrix(unlist(p_v2), ncol = 1)
p_v2_round<-round(p_v2_matrix,digits=4)
p_v2_matrix <-cbind(colnames(x_v2_log),p_v2_round)
p_v2_matrix<-data.frame(p_v2_matrix)
colnames(p_v2_matrix)<-c("names","pvalue")
p_v2_filt<-p_v2_matrix[which(p_v2_matrix[,"pvalue"]<"0.1"),]

my_lms_age_v2<-lapply(1:n2, function(x) lm(x_v2_log[,x] ~ age_v2))
coeff_v2_age<-sapply(my_lms_age_v2, coef)
summaries_v2_age <- lapply(my_lms_age_v2, summary)
p_v2_age<-lapply(summaries_v2_age, function(x) x$coefficients[2,4])
p_v2_age_matrix <- matrix(unlist(p_v2_age), ncol = 1)
p_v2_age_round<-round(p_v2_age_matrix,digits=4)
p_v2_age_matrix <-cbind(colnames(x_v2_log),p_v2_age_round)
p_v2_age_matrix<-data.frame(p_v2_age_matrix)
colnames(p_v2_age_matrix)<-c("names","pvalue")
p_v2_age_filt<-p_v2_age_matrix[which(p_v2_age_matrix[,"pvalue"]<"0.05"),]

my_lms_education_v2<-lapply(1:n2, function(x) lm(x_v2_log[,x] ~ education_v2))
coeff_v2_education<-sapply(my_lms_education_v2, coef)
summaries_v2_education <- lapply(my_lms_education_v2, summary)
p_v2_education<-lapply(summaries_v2_education, function(x) x$coefficients[2,4])
p_v2_education_matrix <- matrix(unlist(p_v2_education), ncol = 1)
p_v2_education_round<-round(p_v2_education_matrix,digits=4)
p_v2_education_matrix <-cbind(colnames(x_v2_log),p_v2_education_round)
p_v2_education_matrix<-data.frame(p_v2_education_matrix)
colnames(p_v2_education_matrix)<-c("names","pvalue")
p_v2_education_filt<-p_v2_education_matrix[which(p_v2_education_matrix[,"pvalue"]<"0.05"),]

my_lms_race_v2<-lapply(1:n2, function(x) lm(x_v2_log[,x] ~ WhitesVsNOn_white_v2_dummy))
coeff_v2_race<-sapply(my_lms_race_v2, coef)
summaries_v2_race<- lapply(my_lms_race_v2, summary)
p_v2_race<-lapply(summaries_v2_race, function(x) x$coefficients[2,4])
p_v2_race_matrix <- matrix(unlist(p_v2_race), ncol = 1)
p_v2_race_round<-round(p_v2_race_matrix,digits=4)
p_v2_race_matrix <-cbind(colnames(x_v2_log),p_v2_race_round)
p_v2_race_matrix<-data.frame(p_v2_race_matrix)
colnames(p_v2_race_matrix)<-c("names","pvalue")
p_v2_race_filt<-p_v2_race_matrix[which(p_v2_race_matrix[,"pvalue"]<"0.05"),]


#####Anxiety scores-----IL5

###########################################################################################################################################################################
#For visit 4

df_new_4<-df[(df$Visit=="V4"),]

df_new_4_education<-df_new_4[!(df_new_4$Participant=="24"),] ############With respect to education only . Removed subject MPI 24

x_v4<-df_new_4[,19:39]
x_v4_log<-x_v4+1
x_v4_log<-log(x_v4_log)
x_v4_log[]<- lapply(x_v4_log, function(x) replace(x, is.na(x), mean(x, na.rm = TRUE)))


x_v4_education<-df_new_4_education[,19:39]
x_v4_log_education<-x_v4_education+1
x_v4_log_education<-log(x_v4_log_education)
x_v4_log_education[]<- lapply(x_v4_log_education, function(x) replace(x, is.na(x), mean(x, na.rm = TRUE)))


y_v4<-df_new_4$GAD_7
age_v4<-df_new_4$Age
bmi_v4<-df_new_4$BMI
education_v4<-df_new_4_education$Years.of.Education

df_new_4$WhitesVsNOn_whites<-rep("Non-White", nrow(df_new_4))
df_new_4$WhitesVsNOn_whites[(df_new_4$Race==0)]<-"White"
WhitesVsNOn_white_v4<-df_new_4$WhitesVsNOn_whites
WhitesVsNOn_white_v4_dummy<-dummy_cols(WhitesVsNOn_white_v4)
WhitesVsNOn_white_v4_dummy<-subset (WhitesVsNOn_white_v4_dummy, select = -c(.data))
WhitesVsNOn_white_v4_dummy<-matrix(unlist(WhitesVsNOn_white_v4_dummy),ncol = 2)

n4<-21

my_lms_v4<- lapply(1:n4, function(x) lm(x_v4_log[,x] ~ y_v4+bmi_v4))
coeff_v4<-sapply(my_lms_v4, coef)
summaries_v4 <- lapply(my_lms_v4, summary)
p_v4<-lapply(summaries_v4, function(x) x$coefficients[2,4])
p_v4_matrix <- matrix(unlist(p_v4), ncol = 1)
p_v4_round<-round(p_v4_matrix,digits=4)
p_v4_matrix <-cbind(colnames(x_v4_log),p_v4_round)
p_v4_matrix<-data.frame(p_v4_matrix)
colnames(p_v4_matrix)<-c("names","pvalue")
p_v4_filt<-p_v4_matrix[which(p_v4_matrix[,"pvalue"]<"0.1"),]

my_lms_age_v4<-lapply(1:n4, function(x) lm(x_v4_log[,x] ~ age_v4))
coeff_v4_age<-sapply(my_lms_age_v4, coef)
summaries_v4_age <- lapply(my_lms_age_v4, summary)
p_v4_age<-lapply(summaries_v4_age, function(x) x$coefficients[2,4])
p_v4_age_matrix <- matrix(unlist(p_v4_age), ncol = 1)
p_v4_age_round<-round(p_v4_age_matrix,digits=4)
p_v4_age_matrix <-cbind(colnames(x_v4_log),p_v4_age_round)
p_v4_age_matrix<-data.frame(p_v4_age_matrix)
colnames(p_v4_age_matrix)<-c("names","pvalue")
p_v4_age_filt<-p_v4_age_matrix[which(p_v4_age_matrix[,"pvalue"]<"0.05"),]

my_lms_education_v4<-lapply(1:n4, function(x) lm(x_v4_log_education[,x] ~ education_v4))
coeff_v4_education<-sapply(my_lms_education_v4, coef)
summaries_v4_education <- lapply(my_lms_education_v4, summary)
p_v4_education<-lapply(summaries_v4_education, function(x) x$coefficients[2,4])
p_v4_education_matrix <- matrix(unlist(p_v4_education), ncol = 1)
p_v4_education_round<-round(p_v4_education_matrix,digits=4)
p_v4_education_matrix <-cbind(colnames(x_v4_log),p_v4_education_round)
p_v4_education_matrix<-data.frame(p_v4_education_matrix)
colnames(p_v4_education_matrix)<-c("names","pvalue")
p_v4_education_filt<-p_v4_education_matrix[which(p_v4_education_matrix[,"pvalue"]<"0.05"),]

my_lms_race_v4<-lapply(1:n4, function(x) lm(x_v4_log[,x] ~ WhitesVsNOn_white_v4_dummy))
coeff_v4_race<-sapply(my_lms_race_v4, coef)
summaries_v4_race<- lapply(my_lms_race_v4, summary)
p_v4_race<-lapply(summaries_v4_race, function(x) x$coefficients[2,4])
p_v4_race_matrix <- matrix(unlist(p_v4_race), ncol = 1)
p_v4_race_round<-round(p_v4_race_matrix,digits=4)
p_v4_race_matrix <-cbind(colnames(x_v4_log),p_v4_race_round)
p_v4_race_matrix<-data.frame(p_v4_race_matrix)
colnames(p_v4_race_matrix)<-c("names","pvalue")
p_v4_race_filt<-p_v4_race_matrix[which(p_v4_race_matrix[,"pvalue"]<"0.05"),]

#####Anxiety scores----GM_CSF,IL2


#For 1st visit

df_new_1$log_IL10<-x_v1_log$IL10

#With Respect to Race

#fig_1 <- ggplot(df_new_1, aes(x=as.factor(WhitesVsNOn_whites), y=log_IL10,fill=as.factor(WhitesVsNOn_whites))) + 
#  geom_boxplot(width=0.75,lwd=0.5)+
#  labs(x="Race", y = "log IL10(pg/ml)")+ggtitle("Visit 1")
#fig_1<-fig_1 + scale_fill_brewer(palette="Dark2")+guides(fill=guide_legend(title='Race'))+theme_bw()+theme(aspect.ratio = 1)
#fig_1<-fig_1+theme(text = element_text(size = 20))
#plot(fig_1)

#For 2nd visit

df_new_2$log_IL5<-x_v2_log$IL5

#With respect to GAD-7

#ADD CONFIDENCE INTERVAL TO THE PLOT#############

fig_1 <- ggplot(df_new_2, aes(x=GAD_7, y=log_IL5)) + 
  geom_point(colour = "Maroon1",size=5)+geom_smooth(method = "lm",color= "black",size=2)+
  labs(x="GAD-7 score", y = "log IL5(pg/ml)")+ggtitle("Visit 2")+theme_bw()+theme(aspect.ratio = 1)
fig_1<-fig_1+theme(text = element_text(size = 20))
plot(fig_1)

#For Postpartum visit

df_new_4$log_GM_CSF<-x_v4_log$GM_CSF
#df_new_4$IL12_p70<-x_v4_log$IL12_p70
#df_new_4$IL13<-x_v4_log$IL13
df_new_4$log_IL2<-x_v4_log$IL2
#df_new_4$log_IFNg<-x_v4_log$IFNg

#With respect to GAD-7

fig_2 <- ggplot(df_new_4, aes(x=GAD_7, y=log_GM_CSF)) + 
  geom_point(colour = "darkgreen",size=5)+geom_smooth(method = "lm",color= "black",size=2)+
  labs(x="GAD-7 score", y = "log GM-CSF(pg/ml)")+ggtitle("Visit 4")+theme_bw()+theme(aspect.ratio = 1)
fig_2<-fig_2+theme(text = element_text(size = 20))
plot(fig_2)

fig_3<- ggplot(df_new_4, aes(x=GAD_7, y=log_IL2)) + 
 geom_point(colour = "chocolate",size=5)+geom_smooth(method = "lm",color="black",size = 2)+
 labs(x="GAD-7 score", y = "log IL2(pg/ml)")+ggtitle("Visit 4")+theme_bw()+theme(aspect.ratio = 1)
fig_3<-fig_3+theme(text = element_text(size = 20))
plot(fig_3)

#With respect to Age

#fig_5<- ggplot(df_new_4, aes(x=Age, y=log_IFNg)) + 
 # geom_point(colour = "darkgreen",size=4)+geom_smooth(method = "lm",color="black",size=2)+
 # labs(x=" Age", y = "log IFNg(pg/ml)")+ggtitle("Visit 4")+theme_bw()+theme(aspect.ratio = 1)
#fig_5<-fig_5+theme(text = element_text(size = 20))
#plot(fig_5)


#######################Metagenomeseq################################################################

path<-'C:/Akanksha/LAB WORK/Cytokines'
list.files(path)

fnFs<-sort(list.files(path, pattern = "_L001_R1_001.fastq"))
fnRs<-sort(list.files(path, pattern = "_L001_R2_001.fastq"))

##Extracting sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

sample.names<-paste("S_",sample.names,sep="")

# Place filtered files in filtered/ subdirectory
filtfwd_df <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtrvs_df<- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtfwd_df) <- sample.names
names(filtrvs_df) <- sample.names

out<-filterAndTrim(fnFs,filtfwd_df,fnRs,filtrvs_df,verbose=TRUE,trimLeft=20,multithread=FALSE)

errF <- learnErrors(filtfwd_df, nbases=1e8, multithread=TRUE,randomize=TRUE,verbose=1)
errR <- learnErrors(filtrvs_df, nbases=1e8, multithread=TRUE,randomize=TRUE,verbose=1)

dadaFwd <- dada(filtfwd_df, err=errF, multithread = TRUE)

dadaRvs <- dada(filtrvs_df, err=errR, multithread = TRUE)
dadaFwd[[45]]


mergers <- mergePairs(dadaFwd, filtfwd_df, dadaRvs, filtrvs_df, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)
# Inspect distribution of sequence lengths
dim(seqtab)         
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)   
sum(seqtab.nochim)/sum(seqtab)

for (i in 1:nrow(seqtab.nochim)) {
  for (j in 1:ncol(seqtab.nochim)) {
    if(seqtab.nochim[i,j]<=10){
      seqtab.nochim[i,j]<-0
    }
    
  }
  
}

seqtabnochim1<-seqtab.nochim[rowSums(seqtab.nochim)>0,]
seqtabnochim2<-seqtabnochim1[,colSums(seqtabnochim1)>0]

dim(seqtabnochim2)  

taxa <- assignTaxonomy(seqtabnochim2, "C:/Akanksha/LAB WORK/Cytokines/Silva/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
taxa_1<- addSpecies(taxa,"C:/Akanksha/LAB WORK/Cytokines/Silva/silva_species_assignment_v138.fa.gz")

taxonomy<-array(taxa_1)
taxonomy0<-data.frame(taxonomy)

mapfile0<-read.csv('C:/Akanksha/LAB WORK/Cytokines/cytokines_UNC_04102021.csv', header = TRUE, stringsAsFactors = FALSE)

samdf<-mapfile0[1:60, ]
samdf1 <- data.frame(samdf[,-1], row.names = samdf[,1])
samdf2<-sample_data(samdf1)

seqs<-getSequences(seqtabnochim2)
names(seqs)<-seqs
alignment<-AlignSeqs(DNAStringSet(seqs), anchor = NA)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR_1 <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                      rearrangement = "stochastic", control = pml.control(trace = 0))

ps<-phyloseq(tax_table(taxa_1),sample_data(samdf2),otu_table(seqtabnochim2,taxa_are_rows = FALSE), phy_tree(fitGTR_1$tree))
rank_names(ps)
table(tax_table(ps)[, "Phylum"], exclude = NULL) 

ps0<- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))  

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps0),
               MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

prevdf_1 = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(ps0),
                      tax_table(ps0))

plyr::ddply(prevdf_1, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

filterPhyla = c("Cyanobacteria","Fusobacteriota","Synergistota","Euryarchaeota","Patescibacteria")  

##Removing these two phyla because their mean is too low and also its sum


tax_table_ps0<-data.frame(tax_table(ps0))
tax_table_ps0 <-tax_table_ps0[(!(tax_table_ps0$Phylum=="Cyanobacteria") & !(tax_table_ps0$Phylum=="Euryarchaeota")& 
                                 !(tax_table_ps0$Phylum=="Fusobacteriota")&
                                 !(tax_table_ps0$Phylum=="Synergistota")& 
                                 !(tax_table_ps0$Phylum=="Patescibacteria")),]


#rownames(tax_table_ps0)<-paste("Seq",1:nrow(tax_table_ps0), sep="_")
tax_table_ps0<-as.matrix(tax_table_ps0)

otu_table_ps0<-data.frame(otu_table(ps0))
#colnames(otu_table_ps0)<-paste("Seq",1:ncol(otu_table_ps0), sep="_")
otu_table_ps0<-t(otu_table_ps0)
otu_table_ps0<-data.frame(otu_table_ps0)
seqs_ps0<-rownames(tax_table_ps0)
otu_table_ps0<-otu_table_ps0[c(seqs_ps0),]
otu_table_ps0<-t(otu_table_ps0)

samdf3<-data.frame(sample_data(ps0))

ps1<-phyloseq(tax_table(tax_table_ps0),sample_data(samdf2),otu_table(otu_table_ps0, taxa_are_rows = FALSE))


#ps1 = subset_taxa(ps0, !Phylum %in% filterPhyla)

# Subset to the remaining phyla
prevdf_2 = subset(prevdf_1, Phylum %in% get_taxa_unique(ps1, "Phylum"))
#rownames(prevdf_2)<-paste("Seq",1:nrow(prevdf_2), sep="_")

# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.1 * nsamples(ps1)

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf_2)[(prevdf_2$Prevalence >= prevalenceThreshold)]

##########ps2-------->keeping taxa which is more than 5% of total samples####################################
ps2 = prune_taxa(keepTaxa, ps1)

ps3 = tax_glom(ps2, "Genus", NArm = TRUE)
#####################################################################################################################################################################
#Showing relationship between cytokine and microbiome only for cytokines that are significant with respect to anxiety scores
set.seed(100)

x_cyto_meta<-cbind(sample.names,x_log)
x_cyto_meta<- data.frame(x_cyto_meta[,-1], row.names = x_cyto_meta[,1])
x_cyto_meta<-sample_data(x_cyto_meta)
otu_table_cleanup<-data.frame(otu_table(ps2))
colnames(otu_table_cleanup)<-paste("Seq",1:ncol(otu_table_cleanup), sep="_")
otu_table_cleanup<-otu_table_cleanup[,colSums(otu_table_cleanup)>0]
otu_table_cleanup<-otu_table_cleanup[rowSums(otu_table_cleanup)>0,]
tax_split<-as.matrix(tax_table(ps2))
rownames(tax_split)<-paste("Seq",1:nrow(tax_split), sep="_")

final_filtered<-phyloseq(tax_table(tax_split),sample_data(x_cyto_meta),otu_table(otu_table_cleanup,taxa_are_rows = FALSE))

metaResults<-c()

for (icyto in 1:ncol(x_cyto_meta)){
  print(icyto)
  metadata<-sample_data(final_filtered)
  matrixData<-matrix(otu_table(final_filtered),ncol=ncol(otu_table(final_filtered)))
  rownames(matrixData)<-rownames(otu_table(final_filtered))
  colnames(matrixData)<-colnames(otu_table(final_filtered))
  matrixData<-t(matrixData)
  featureData =data.frame(tax_table(final_filtered))
  otus_metagenomeSeq<-newMRexperiment(matrixData, 
                                      phenoData = AnnotatedDataFrame(metadata), 
                                      featureData = AnnotatedDataFrame(featureData))
  otus_metagenomeSeq_filter = filterData(otus_metagenomeSeq,present=1, depth = 1) 
  p = cumNormStatFast(otus_metagenomeSeq_filter)
  cytokinesData<-cumNorm(otus_metagenomeSeq_filter, p =p)
  cytoVariable=pData(cytokinesData)[,icyto]
  
  mod_cyto= model.matrix(~cytoVariable+BMI) 
  settings = zigControl(maxit = 10, verbose = TRUE) 
  fit_cyto = fitZig(obj = cytokinesData,mod = mod_cyto,useCSSoffset = FALSE, control = settings)
  mean_eff_size=median(calculateEffectiveSamples(fit_cyto))
  otus_metagenomeSeq_filter = filterData(otus_metagenomeSeq, present = round(mean_eff_size,0), depth = 1) 
  p = cumNormStatFast(otus_metagenomeSeq_filter)
  cytokinesData<-cumNorm(otus_metagenomeSeq_filter, p =p) 
  
  #Define the normalisation factor
  norm.factor <- normFactors(cytokinesData)
  norm.factor <- log2(norm.factor/median(norm.factor)+1)    #####min value is 0.0061 and max value is 1.462 
  cytoVariable=pData(cytokinesData)[,icyto]
  mod_cyto=model.matrix(~cytoVariable+BMI+norm.factor) 
  settings = zigControl(maxit = 10, verbose = TRUE) 
  fit_cyto = fitZig(obj = cytokinesData, mod = mod_cyto, useCSSoffset =FALSE, control = settings)
  zigFit = fit_cyto@fit 
  finalMod = fit_cyto@fit
  fit_cyto = eBayes(zigFit)
  topTable(fit_cyto)
  
  pvalues_all<-apply(cbind(fit_cyto$p.value[,"cytoVariable"]),2,function(x){p.adjust(as.numeric(x),method="fdr")})
  
  otus_mG<-cbind.data.frame(fit_cyto$coefficients[,"cytoVariable"], pvalues_all)
  colnames(otus_mG)<-c("cytoVariable","cytoVariable.1")
  
  number_sig_seq<-sum(apply(otus_mG[,(ncol(otus_mG)/2):ncol(otus_mG)],1,function(x){ ifelse(any(as.numeric(x)<0.1), TRUE, FALSE)}))
  otus_mG_filtered<-otus_mG[apply(otus_mG[,(ncol(otus_mG)/2):ncol(otus_mG)],1,function(x){ifelse(any(as.numeric(x)<0.1), TRUE, FALSE)}),]

  number_sig_seq_icyto<-sum(apply(otus_mG[,grep("cytoVariable", colnames(otus_mG))],1,function(x){ifelse(as.numeric(x[2])<0.1, TRUE, FALSE)}))
  otus_mG_filtered_icyto<-otus_mG[apply(otus_mG[,grep("cytoVariable", colnames(otus_mG))],1,function(x){ifelse(as.numeric(x[2])<0.1, TRUE, FALSE)}),]
  number_sig_seq<-nrow(otus_mG_filtered)
  
  if (number_sig_seq>0){
    if (!is.null(dim(otus_mG_filtered))){
      tax_table_sig_otu<-merge(tax_table(final_filtered), otus_mG_filtered, by=0)
      
      tax_table_sig_otu_icyto<-merge(tax_table(final_filtered), otus_mG_filtered_icyto, by=0)
      tax_table_sig_otu_icyto$Genus<-as.character(tax_table_sig_otu_icyto$Genus)
      tax_table_sig_otu_icyto<-tax_table_sig_otu_icyto[order(tax_table_sig_otu_icyto[,"cytoVariable"], decreasing = TRUE),]
      tax_table_sig_otu_icyto$Row.names<-factor(tax_table_sig_otu_icyto$Row.names, 
                                                levels=tax_table_sig_otu_icyto$Row.names) 
      metaResults<-rbind(metaResults,cbind(rep(icyto, nrow(tax_table_sig_otu_icyto)),
                                       as.character(tax_table_sig_otu_icyto$Row.names), 
                                       tax_table_sig_otu_icyto$Genus, tax_table_sig_otu_icyto$cytoVariable, 
                                       tax_table_sig_otu_icyto$cytoVariable.1))  
      colnames(metaResults)<-c("Cytokine_number","Seq_no","Genus","Coefficient","Pvalue")
 
    }
  }  
}


sign<-sapply(as.numeric(metaResults[,"Coefficient"]),sign)
metaResults<-cbind.data.frame(metaResults,sign)
metaResults<-na.omit(metaResults)
metaResults$cytokine_names<-rep("ITAC",nrow(metaResults))
metaResults$cytokine_names[(metaResults$Cytokine_number==2)]<-"GM-CSF"
metaResults$cytokine_names[(metaResults$Cytokine_number==3)]<-"Fracktalkine"
#metaResults$cytokine_names[(metaResults$Cytokine_number==4)]<-"IFNg"
#metaResults$cytokine_names[(metaResults$Cytokine_number==5)]<-"IL10"
metaResults$cytokine_names[(metaResults$Cytokine_number==6)]<-"MIP3a"
#metaResults$cytokine_names[(metaResults$Cytokine_number==7)]<-"IL12_p70"
metaResults$cytokine_names[(metaResults$Cytokine_number==8)]<-"IL13"
metaResults$cytokine_names[(metaResults$Cytokine_number==9)]<-"IL17A"
#metaResults$cytokine_names[(metaResults$Cytokine_number==11)]<-"IL2"
#metaResults$cytokine_names[(metaResults$Cytokine_number==12)]<-"IL21"
#metaResults$cytokine_names[(metaResults$Cytokine_number==13)]<-"IL4"
metaResults$cytokine_names[(metaResults$Cytokine_number==14)]<-"IL23"
metaResults$cytokine_names[(metaResults$Cytokine_number==15)]<-"IL5"
#metaResults$cytokine_names[(metaResults$Cytokine_number==16)]<-"IL6"
metaResults$cytokine_names[(metaResults$Cytokine_number==17)]<-"IL7"
metaResults$cytokine_names[(metaResults$Cytokine_number==19)]<-"MIP1a"
#metaResults$cytokine_names[(metaResults$Cytokine_number==20)]<-"MIP1B"
metaResults$cytokine_names[(metaResults$Cytokine_number==21)]<-"TNFa"

write.table(metaResults,"C:/Akanksha/LAB WORK/Cytokines/cytokines_microbiome.txt",
            sep="\t", col.names =TRUE,row.names = FALSE)

####Not significant with respect to IL6,IL1B, IL8

################For Visit 1 #################################################
set.seed(100)
sample.names_v1<-df_new_1$SampleName

x_cyto_meta_v1<-x_cyto_meta[c(sample.names_v1),]
otu_table_cleanup_v1<-otu_table_cleanup[c(sample.names_v1),]

otu_table_cleanup_v1<-otu_table_cleanup_v1[,colSums(otu_table_cleanup_v1)>0]
otu_table_cleanup_v1<-otu_table_cleanup_v1[rowSums(otu_table_cleanup_v1)>0,]

#otu_table_cleanup_v1<-t(otu_table_cleanup_v1)

otu_table_cleanup_v1<-as.matrix(otu_table_cleanup_v1)
tax_split_v1<-as.matrix(tax_split[match(colnames(otu_table_cleanup_v1),rownames(tax_split)),])

final_filtered_v1<-phyloseq(tax_table(tax_split_v1),sample_data(x_cyto_meta_v1),otu_table(otu_table_cleanup_v1,taxa_are_rows = FALSE))

metaResults_v1<-c()


for (icyto_v1 in 1:ncol(x_cyto_meta_v1)){
  print(icyto_v1)
  metadata_v1<-sample_data(final_filtered_v1)
  matrixData_v1<-matrix(otu_table(final_filtered_v1),ncol=ncol(otu_table(final_filtered_v1)))
  rownames(matrixData_v1)<-rownames(otu_table(final_filtered_v1))
  colnames(matrixData_v1)<-colnames(otu_table(final_filtered_v1))
  matrixData_v1<-t(matrixData_v1)
  featureData_v1 =data.frame(tax_table(final_filtered_v1))
  otus_metagenomeSeq_v1<-newMRexperiment(matrixData_v1, 
                                      phenoData = AnnotatedDataFrame(metadata_v1), 
                                      featureData = AnnotatedDataFrame(featureData_v1))
  otus_metagenomeSeq_filter_v1 = filterData(otus_metagenomeSeq_v1,present=1, depth = 1) 
  p_v1 = cumNormStatFast(otus_metagenomeSeq_filter_v1)
  cytokinesData_v1<-cumNorm(otus_metagenomeSeq_filter_v1, p =p_v1)
  cytoVariable_v1=pData(cytokinesData_v1)[,icyto_v1]
  
  mod_cyto_v1= model.matrix(~cytoVariable_v1+bmi_v1) 
  settings = zigControl(maxit = 10, verbose = TRUE) 
  fit_cyto_v1 = fitZig(obj = cytokinesData_v1,mod = mod_cyto_v1,useCSSoffset = FALSE, control = settings)
  mean_eff_size_v1 =median(calculateEffectiveSamples(fit_cyto_v1))
  otus_metagenomeSeq_filter_v1 = filterData(otus_metagenomeSeq_v1, present = round(mean_eff_size_v1,0), depth = 1) 
  p_v1 = cumNormStatFast(otus_metagenomeSeq_filter_v1)
  cytokinesData_v1<-cumNorm(otus_metagenomeSeq_filter_v1, p =p_v1) 
  
  #Define the normalisation factor
  norm.factor_v1 <- normFactors(cytokinesData_v1)
  norm.factor_v1 <- log2(norm.factor_v1/median(norm.factor_v1)+1)    #####min value is 0.0061 and max value is 1.462 
  cytoVariable_v1=pData(cytokinesData_v1)[,icyto_v1]
  mod_cyto_v1=model.matrix(~cytoVariable_v1+bmi_v1+norm.factor_v1) 
  settings = zigControl(maxit = 10, verbose = TRUE) 
  fit_cyto_v1 = fitZig(obj = cytokinesData_v1, mod = mod_cyto_v1, useCSSoffset =FALSE, control = settings)
  zigFit_v1 = fit_cyto_v1@fit 
  finalMod_v1 = fit_cyto_v1@fit
  fit_cyto_v1 = eBayes(zigFit_v1)
  topTable(fit_cyto_v1)
 
  pvalues_v1<-apply(cbind(fit_cyto_v1$p.value[,"cytoVariable_v1"]),2,function(x){p.adjust(as.numeric(x),method="fdr")})
  
  otus_mG_v1<-cbind.data.frame(fit_cyto_v1$coefficients[,"cytoVariable_v1"], pvalues_v1)
 
  colnames(otus_mG_v1)<-c("cytoVariable_v1","cytoVariable_v1.1")
  
  number_sig_seq_v1<-sum(apply(otus_mG_v1[,(ncol(otus_mG_v1)/2):ncol(otus_mG_v1)],1,function(x){ ifelse(any(as.numeric(x)<0.05), TRUE, FALSE)}))
  otus_mG_filtered_v1<-otus_mG_v1[apply(otus_mG_v1[,(ncol(otus_mG_v1)/2):ncol(otus_mG_v1)],1,function(x){ifelse(any(as.numeric(x)<0.05), TRUE, FALSE)}),]
  
  number_sig_seq_icyto_v1<-sum(apply(otus_mG_v1[,grep("cytoVariable_v1", colnames(otus_mG_v1))],1,function(x){ifelse(as.numeric(x[2])<0.05, TRUE, FALSE)}))
  otus_mG_filtered_icyto_v1<-otus_mG_v1[apply(otus_mG_v1[,grep("cytoVariable_v1", colnames(otus_mG_v1))],1,function(x){ifelse(as.numeric(x[2])<0.05, TRUE, FALSE)}),]
  number_sig_seq_v1<-nrow(otus_mG_filtered_v1)
 
  if (number_sig_seq_v1>0){
    if (!is.null(dim(otus_mG_filtered_v1))){
      tax_table_sig_otu_v1<-merge(tax_table(final_filtered_v1), otus_mG_filtered_v1, by=0)
      
      tax_table_sig_otu_icyto_v1<-merge(tax_table(final_filtered_v1), otus_mG_filtered_icyto_v1, by=0)
      tax_table_sig_otu_icyto_v1$Genus<-as.character(tax_table_sig_otu_icyto_v1$Genus)
      tax_table_sig_otu_icyto_v1<-tax_table_sig_otu_icyto_v1[order(tax_table_sig_otu_icyto_v1[,"cytoVariable_v1"], decreasing = TRUE),]
      tax_table_sig_otu_icyto_v1$Row.names<-factor(tax_table_sig_otu_icyto_v1$Row.names, 
                                                levels=tax_table_sig_otu_icyto_v1$Row.names) 
      metaResults_v1<-rbind(metaResults_v1,cbind(rep(icyto_v1, nrow(tax_table_sig_otu_icyto_v1)),
                                           as.character(tax_table_sig_otu_icyto_v1$Row.names), 
                                           tax_table_sig_otu_icyto_v1$Genus, tax_table_sig_otu_icyto_v1$cytoVariable_v1, 
                                           tax_table_sig_otu_icyto_v1$cytoVariable_v1.1))  
      colnames(metaResults_v1)<-c("Cytokine_number","Seq_no","Genus","Coefficient","Pvalue")
      
    }
  }  
}

sign_v1<-sapply(as.numeric(metaResults_v1[,"Coefficient"]),sign)
metaResults_v1<-cbind.data.frame(metaResults_v1,sign_v1)
metaResults_v1<-na.omit(metaResults_v1)


metaResults_v1$cytokine_names<-rep("ITAC",nrow(metaResults_v1))


metaResults_v1$cytokine_names[(metaResults_v1$Cytokine_number==2)]<-"GM-CSF"
metaResults_v1$cytokine_names[(metaResults_v1$Cytokine_number==3)]<-"Fracktalkine"
metaResults_v1$cytokine_names[(metaResults_v1$Cytokine_number==4)]<-"IFNg"

metaResults_v1$cytokine_names[(metaResults_v1$Cytokine_number==5)]<-"IL10"
metaResults_v1$cytokine_names[(metaResults_v1$Cytokine_number==7)]<-"IL12_p70"
metaResults_v1$cytokine_names[(metaResults_v1$Cytokine_number==9)]<-"IL17A"
metaResults_v1$cytokine_names[(metaResults_v1$Cytokine_number==12)]<-"IL21"
metaResults_v1$cytokine_names[(metaResults_v1$Cytokine_number==14)]<-"IL23"
metaResults_v1$cytokine_names[(metaResults_v1$Cytokine_number==15)]<-"IL5"

metaResults_v1$cytokine_names[(metaResults_v1$Cytokine_number==17)]<-"IL7"

metaResults_v1$cytokine_names[(metaResults_v1$Cytokine_number==18)]<-"IL8"
metaResults_v1$cytokine_names[(metaResults_v1$Cytokine_number==19)]<-"MIP1a"
metaResults_v1$cytokine_names[(metaResults_v1$Cytokine_number==20)]<-"MIP1B"

write.table(metaResults_v1,"C:/Akanksha/LAB WORK/Cytokines/cytokines_microbiome_V1.txt",
            sep="\t", col.names =TRUE,row.names = FALSE)

##########################Visit 2 ########################################
set.seed(100)
sample.names_v2<-df_new_2$SampleName

x_cyto_meta_v2<-x_cyto_meta[c(sample.names_v2),]
otu_table_cleanup_v2<-otu_table_cleanup[c(sample.names_v2),]

otu_table_cleanup_v2<-otu_table_cleanup_v2[,colSums(otu_table_cleanup_v2)>0]
otu_table_cleanup_v2<-otu_table_cleanup_v2[rowSums(otu_table_cleanup_v2)>0,]

#otu_table_cleanup_v1<-t(otu_table_cleanup_v1)

otu_table_cleanup_v2<-as.matrix(otu_table_cleanup_v2)
tax_split_v2<-as.matrix(tax_split[match(colnames(otu_table_cleanup_v2),rownames(tax_split)),])

final_filtered_v2<-phyloseq(tax_table(tax_split_v2),sample_data(x_cyto_meta_v2),otu_table(otu_table_cleanup_v2,taxa_are_rows = FALSE))

metaResults_v2<-c()

for (icyto_v2 in 1:ncol(x_cyto_meta_v2)){
  print(icyto_v2)
  metadata_v2<-sample_data(final_filtered_v2)
  matrixData_v2<-matrix(otu_table(final_filtered_v2),ncol=ncol(otu_table(final_filtered_v2)))
  rownames(matrixData_v2)<-rownames(otu_table(final_filtered_v2))
  colnames(matrixData_v2)<-colnames(otu_table(final_filtered_v2))
  matrixData_v2<-t(matrixData_v2)
  featureData_v2 =data.frame(tax_table(final_filtered_v2))
  otus_metagenomeSeq_v2<-newMRexperiment(matrixData_v2, 
                                         phenoData = AnnotatedDataFrame(metadata_v2), 
                                         featureData = AnnotatedDataFrame(featureData_v2))
  otus_metagenomeSeq_filter_v2 = filterData(otus_metagenomeSeq_v2,present=1, depth = 1) 
  p_v2 = cumNormStatFast(otus_metagenomeSeq_filter_v2)
  cytokinesData_v2<-cumNorm(otus_metagenomeSeq_filter_v2, p =p_v2)
  cytoVariable_v2=pData(cytokinesData_v2)[,icyto_v2]
  
  mod_cyto_v2= model.matrix(~cytoVariable_v2+bmi_v2) 
  settings = zigControl(maxit = 10, verbose = TRUE) 
  fit_cyto_v2 = fitZig(obj = cytokinesData_v2,mod = mod_cyto_v2,useCSSoffset = FALSE, control = settings)
  mean_eff_size_v2 =median(calculateEffectiveSamples(fit_cyto_v2))
  otus_metagenomeSeq_filter_v2 = filterData(otus_metagenomeSeq_v2, present = round(mean_eff_size_v2,0), depth = 1) 
  p_v2 = cumNormStatFast(otus_metagenomeSeq_filter_v2)
  cytokinesData_v2<-cumNorm(otus_metagenomeSeq_filter_v2, p =p_v2) 
  
  #Define the normalisation factor
  norm.factor_v2 <- normFactors(cytokinesData_v2)
  norm.factor_v2 <- log2(norm.factor_v2/median(norm.factor_v2)+1)    #####min value is 0.044 and max value is 1.501
  cytoVariable_v2=pData(cytokinesData_v2)[,icyto_v2]
  mod_cyto_v2=model.matrix(~cytoVariable_v2+bmi_v2+norm.factor_v2) 
  settings = zigControl(maxit = 10, verbose = TRUE) 
  fit_cyto_v2 = fitZig(obj = cytokinesData_v2, mod = mod_cyto_v2, useCSSoffset =FALSE, control = settings)
  zigFit_v2 = fit_cyto_v2@fit 
  finalMod_v2 = fit_cyto_v2@fit
  fit_cyto_v2 = eBayes(zigFit_v2)
  topTable(fit_cyto_v2)
  
  pvalues_v2<-apply(cbind(fit_cyto_v2$p.value[,"cytoVariable_v2"]),2,function(x){p.adjust(as.numeric(x),method="fdr")})
  
  otus_mG_v2<-cbind.data.frame(fit_cyto_v2$coefficients[,"cytoVariable_v2"], pvalues_v2)
  
  colnames(otus_mG_v2)<-c("cytoVariable_v2","cytoVariable_v2.1")
  
  number_sig_seq_v2<-sum(apply(otus_mG_v2[,(ncol(otus_mG_v2)/2):ncol(otus_mG_v2)],1,function(x){ ifelse(any(as.numeric(x)<0.01), TRUE, FALSE)}))
  otus_mG_filtered_v2<-otus_mG_v2[apply(otus_mG_v2[,(ncol(otus_mG_v2)/2):ncol(otus_mG_v2)],1,function(x){ifelse(any(as.numeric(x)<0.01), TRUE, FALSE)}),]
  
  number_sig_seq_icyto_v2<-sum(apply(otus_mG_v2[,grep("cytoVariable_v2", colnames(otus_mG_v2))],1,function(x){ifelse(as.numeric(x[2])<0.01, TRUE, FALSE)}))
  otus_mG_filtered_icyto_v2<-otus_mG_v2[apply(otus_mG_v2[,grep("cytoVariable_v2", colnames(otus_mG_v2))],1,function(x){ifelse(as.numeric(x[2])<0.01, TRUE, FALSE)}),]
  number_sig_seq_v2<-nrow(otus_mG_filtered_v2)
  
  if (number_sig_seq_v2>0){
    if (!is.null(dim(otus_mG_filtered_v2))){
      tax_table_sig_otu_v2<-merge(tax_table(final_filtered_v2), otus_mG_filtered_v2, by=0)
      
      tax_table_sig_otu_icyto_v2<-merge(tax_table(final_filtered_v2), otus_mG_filtered_icyto_v2, by=0)
      tax_table_sig_otu_icyto_v2$Genus<-as.character(tax_table_sig_otu_icyto_v2$Genus)
      tax_table_sig_otu_icyto_v2<-tax_table_sig_otu_icyto_v2[order(tax_table_sig_otu_icyto_v2[,"cytoVariable_v2"], decreasing = TRUE),]
      tax_table_sig_otu_icyto_v2$Row.names<-factor(tax_table_sig_otu_icyto_v2$Row.names, 
                                                   levels=tax_table_sig_otu_icyto_v2$Row.names) 
      metaResults_v2<-rbind(metaResults_v2,cbind(rep(icyto_v2, nrow(tax_table_sig_otu_icyto_v2)),
                                                 as.character(tax_table_sig_otu_icyto_v2$Row.names), 
                                                 tax_table_sig_otu_icyto_v2$Genus, tax_table_sig_otu_icyto_v2$cytoVariable_v2, 
                                                 tax_table_sig_otu_icyto_v2$cytoVariable_v2.1))  
      colnames(metaResults_v2)<-c("Cytokine_number","Seq_no","Genus","Coefficient","Pvalue")
      
    }
  }  
}


sign_v2<-sapply(as.numeric(metaResults_v2[,"Coefficient"]),sign)
metaResults_v2<-cbind.data.frame(metaResults_v2,sign_v2)
metaResults_v2<-na.omit(metaResults_v2)

metaResults_v2$cytokine_names<-rep("ITAC",nrow(metaResults_v2))

metaResults_v2$cytokine_names[(metaResults_v2$Cytokine_number==3)]<-"Fracktalkine"

metaResults_v2$cytokine_names[(metaResults_v2$Cytokine_number==5)]<-"IL10"
metaResults_v2$cytokine_names[(metaResults_v2$Cytokine_number==6)]<-"MIP3a"


metaResults_v2$cytokine_names[(metaResults_v2$Cytokine_number==7)]<-"IL12_p70"
metaResults_v2$cytokine_names[(metaResults_v2$Cytokine_number==10)]<-"IL1b"
metaResults_v2$cytokine_names[(metaResults_v2$Cytokine_number==13)]<-"IL4"


metaResults_v2$cytokine_names[(metaResults_v2$Cytokine_number==14)]<-"IL23"

metaResults_v2$cytokine_names[(metaResults_v2$Cytokine_number==17)]<-"IL7"
metaResults_v2$cytokine_names[(metaResults_v2$Cytokine_number==20)]<-"MIP1B"

write.table(metaResults_v2,"C:/Akanksha/LAB WORK/Cytokines/cytokines_microbiome_V2.txt",
            sep="\t", col.names =TRUE,row.names = FALSE)
###########################Visit 4########################################################################
set.seed(100)
df_4<-df[(df$Visit=="V4"),]
sample.names_v4<-df_4$SampleName

BMI_v4<-df_4$BMI

x_cyto_meta_v4<-x_cyto_meta[c(sample.names_v4),]
otu_table_cleanup_v4<-otu_table_cleanup[c(sample.names_v4),]

otu_table_cleanup_v4<-otu_table_cleanup_v4[,colSums(otu_table_cleanup_v4)>0]
otu_table_cleanup_v4<-otu_table_cleanup_v4[rowSums(otu_table_cleanup_v4)>0,]

#otu_table_cleanup_v1<-t(otu_table_cleanup_v1)

otu_table_cleanup_v4<-as.matrix(otu_table_cleanup_v4)
tax_split_v4<-as.matrix(tax_split[match(colnames(otu_table_cleanup_v4),rownames(tax_split)),])

final_filtered_v4<-phyloseq(tax_table(tax_split_v4),sample_data(x_cyto_meta_v4),otu_table(otu_table_cleanup_v4,taxa_are_rows = FALSE))

metaResults_v4<-c()

for (icyto_v4 in 1:ncol(x_cyto_meta_v4)){
  print(icyto_v4)
  metadata_v4<-sample_data(final_filtered_v4)
  matrixData_v4<-matrix(otu_table(final_filtered_v4),ncol=ncol(otu_table(final_filtered_v4)))
  rownames(matrixData_v4)<-rownames(otu_table(final_filtered_v4))
  colnames(matrixData_v4)<-colnames(otu_table(final_filtered_v4))
  matrixData_v4<-t(matrixData_v4)
  featureData_v4 = data.frame(tax_table(final_filtered_v4))
  otus_metagenomeSeq_v4<-newMRexperiment(matrixData_v4, 
                                         phenoData = AnnotatedDataFrame(metadata_v4), 
                                         featureData = AnnotatedDataFrame(featureData_v4))
  otus_metagenomeSeq_filter_v4 = filterData(otus_metagenomeSeq_v4,present=1, depth = 1) 
  p_v4 = cumNormStatFast(otus_metagenomeSeq_filter_v4)
  cytokinesData_v4<-cumNorm(otus_metagenomeSeq_filter_v4, p =p_v4)
  cytoVariable_v4 = pData(cytokinesData_v4)[,icyto_v4]
  
  mod_cyto_v4 = model.matrix(~cytoVariable_v4+BMI_v4) 
  settings = zigControl(maxit = 10, verbose = TRUE) 
  fit_cyto_v4 = fitZig(obj = cytokinesData_v4,mod = mod_cyto_v4,useCSSoffset = FALSE, control = settings)
  mean_eff_size_v4 =median(calculateEffectiveSamples(fit_cyto_v4))
  otus_metagenomeSeq_filter_v4 = filterData(otus_metagenomeSeq_v4, present = round(mean_eff_size_v4,0), depth = 1) 
  p_v4 = cumNormStatFast(otus_metagenomeSeq_filter_v4)
  cytokinesData_v4<-cumNorm(otus_metagenomeSeq_filter_v4, p =p_v4) 
  
  #Define the normalisation factor
  norm.factor_v4 <- normFactors(cytokinesData_v4)
  norm.factor_v4 <- log2(norm.factor_v4/median(norm.factor_v4)+1)    #####min value is 0.0061 and max value is 1.462 
  cytoVariable_v4 = pData(cytokinesData_v4)[,icyto_v4]
  mod_cyto_v4 = model.matrix(~cytoVariable_v4+BMI_v4+norm.factor_v4) 
  settings = zigControl(maxit = 10, verbose = TRUE) 
  fit_cyto_v4 = fitZig(obj = cytokinesData_v4, mod = mod_cyto_v4, useCSSoffset =FALSE, control = settings)
  zigFit_v4 = fit_cyto_v4@fit 
  finalMod_v4 = fit_cyto_v4@fit
  fit_cyto_v4 = eBayes(zigFit_v4)
  topTable(fit_cyto_v4)
  
  pvalues_v4<-apply(cbind(fit_cyto_v4$p.value[,"cytoVariable_v4"]),2,function(x){p.adjust(as.numeric(x),method="fdr")})
  
  otus_mG_v4<-cbind.data.frame(fit_cyto_v4$coefficients[,"cytoVariable_v4"], pvalues_v4)
  
  colnames(otus_mG_v4)<-c("cytoVariable_v4","cytoVariable_v4.1")
 
  number_sig_seq_v4<-sum(apply(otus_mG_v4[,(ncol(otus_mG_v4)/2):ncol(otus_mG_v4)],1,function(x){ ifelse(any(as.numeric(x)<0.05), TRUE, FALSE)}))
  otus_mG_filtered_v4<-otus_mG_v4[apply(otus_mG_v4[,(ncol(otus_mG_v4)/2):ncol(otus_mG_v4)],1,function(x){ifelse(any(as.numeric(x)<0.05), TRUE, FALSE)}),]
  
  number_sig_seq_icyto_v4<-sum(apply(otus_mG_v4[,grep("cytoVariable_v4", colnames(otus_mG_v4))],1,function(x){ifelse(as.numeric(x[2])<0.05, TRUE, FALSE)}))
  otus_mG_filtered_icyto_v4<-otus_mG_v4[apply(otus_mG_v4[,grep("cytoVariable_v4", colnames(otus_mG_v4))],1,function(x){ifelse(as.numeric(x[2])<0.05, TRUE, FALSE)}),]
  number_sig_seq_v4<-nrow(otus_mG_filtered_v4)
  
  if (number_sig_seq_v4>0){
    if (!is.null(dim(otus_mG_filtered_v4))){
      tax_table_sig_otu_v4<-merge(tax_table(final_filtered_v4), otus_mG_filtered_v4, by=0)
      
      tax_table_sig_otu_icyto_v4<-merge(tax_table(final_filtered_v4), otus_mG_filtered_icyto_v4, by=0)
      tax_table_sig_otu_icyto_v4$Genus<-as.character(tax_table_sig_otu_icyto_v4$Genus)
      tax_table_sig_otu_icyto_v4<-tax_table_sig_otu_icyto_v4[order(tax_table_sig_otu_icyto_v4[,"cytoVariable_v4"], decreasing = TRUE),]
      tax_table_sig_otu_icyto_v4$Row.names<-factor(tax_table_sig_otu_icyto_v4$Row.names, 
                                                   levels=tax_table_sig_otu_icyto_v4$Row.names) 
      metaResults_v4<-rbind(metaResults_v4,cbind(rep(icyto_v4, nrow(tax_table_sig_otu_icyto_v4)),
                                                 as.character(tax_table_sig_otu_icyto_v4$Row.names), 
                                                 tax_table_sig_otu_icyto_v4$Genus, tax_table_sig_otu_icyto_v4$cytoVariable_v4, 
                                                 tax_table_sig_otu_icyto_v4$cytoVariable_v4.1))  
      colnames(metaResults_v4)<-c("Cytokine_number","Seq_no","Genus","Coefficient","Pvalue")
      
    }
  }  
}

sign_v4<-sapply(as.numeric(metaResults_v4[,"Coefficient"]),sign)
metaResults_v4<-cbind.data.frame(metaResults_v4,sign_v4)
metaResults_v4<-na.omit(metaResults_v4)

metaResults_v4$cytokine_names<-rep("ITAC",nrow(metaResults_v4))
metaResults_v4$cytokine_names[(metaResults_v4$Cytokine_number==2)]<-"GM-CSF"
metaResults_v4$cytokine_names[(metaResults_v4$Cytokine_number==4)]<-"IFNg"

metaResults_v4$cytokine_names[(metaResults_v4$Cytokine_number==5)]<-"IL10"

metaResults_v4$cytokine_names[(metaResults_v4$Cytokine_number==7)]<-"IL12_p70"
metaResults_v4$cytokine_names[(metaResults_v4$Cytokine_number==9)]<-"IL17A"
metaResults_v4$cytokine_names[(metaResults_v4$Cytokine_number==10)]<-"IL1b"
metaResults_v4$cytokine_names[(metaResults_v4$Cytokine_number==11)]<-"IL2"


metaResults_v4$cytokine_names[(metaResults_v4$Cytokine_number==13)]<-"IL4"
metaResults_v4$cytokine_names[(metaResults_v4$Cytokine_number==14)]<-"IL23"
metaResults_v4$cytokine_names[(metaResults_v4$Cytokine_number==15)]<-"IL5"

metaResults_v4$cytokine_names[(metaResults_v4$Cytokine_number==17)]<-"IL7"
metaResults_v4$cytokine_names[(metaResults_v4$Cytokine_number==19)]<-"MIP1a"
metaResults_v4$cytokine_names[(metaResults_v4$Cytokine_number==20)]<-"MIP1B"



write.table(metaResults_v4,"C:/Akanksha/LAB WORK/Cytokines/cytokines_microbiome_V4.txt",
            sep="\t", col.names =TRUE,row.names = FALSE)
