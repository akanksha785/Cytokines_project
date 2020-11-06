library(readxl)
library(tidyverse)
library(standardize)
library(ggplot2)
library(ggpubr)
library(gginnards)
library(fastDummies)
library(metagenomeSeq)
library(corrplot)
library(complexHeatmap)
library(pheatmap)


df<-read.csv('C:/Akanksha/LAB WORK/Cytokines/cytokines_UNC_10142020.csv',header = TRUE, stringsAsFactors = FALSE)
#df<- na.omit(df)

df<-df[1:63,]
#df<-df%>%select(Ethnicity:Race)
df_new<-dummy_cols(df,select_columns = c('Ethnicity','Race'))

x<-df_new %>% select(ITAC_pg_ml:TNFa_pg_ml)
x_matrix<-data.matrix(x)
y<-df_new$GAD_7
x_std<-scale(x) 
x_std<-data.frame(x_std)
#x_std<- dummy.data.frame(x_std)
#y_std<-scale(y)
#y_std<-data.frame(y_std)
#x_log<-log(x)
GAD_7<-df_new$GAD_7

df

#For Visit V1

df_new_1<-df[(df$Visit=="V1"),]
df_new_1<-df_new_1 %>%select_if(~ !any(is.na(.)))

df_v1<-df_new[(df$Visit=="V1"),]
df_v1<-df_v1 %>%select_if(~ !any(is.na(.)))

x_v1<-df_v1 %>% select(ITAC_pg_ml:TNFa_pg_ml)
#x_v1_log<-log(x_v1)
y_v1<-df_v1$GAD_7
age_v1<-df_v1$Age
education_v1<-df_v1$Years.of.Education
#ethnicity_v1<-df_v1$Ethnicity
#race_v1<-df_v1$Race
#education_v1<-df_v1%>%select(Years.of.Education_13:Years.of.Education_20)
#education_v1_13<-education_v1$Years.of.Education_13
#education_v1_16<-education_v1$Years.of.Education_16
#education_v1_17<-education_v1$Years.of.Education_17
#education_v1_18<-education_v1$Years.of.Education_18
#education_v1_19<-education_v1$Years.of.Education_19
#education_v1_20<-education_v1$Years.of.Education_20
#education_v1<-data.frame(education_v1)
ethnicity_v1<-df_v1$Ethnicity
ethnicity_v1_dummy<-df_v1%>%select(Ethnicity_0:Ethnicity_1)
ethnicity_v1_dummy<-matrix(unlist(ethnicity_v1_dummy),ncol = 2)
#ethnicity_v1_0<-ethnicity_v1$Ethnicity_0
#ethnicity_v1_1<-ethnicity_v1$Ethnicity_1
#ethnicity_v1<-cbind.data.frame(ethnicity_v1_0,ethnicity_v1_1)
race<-df_v1$Race
race_v1_dummy<-df_v1%>%select(Race_0:Race_3)
race_v1_dummy<-matrix(unlist(race_v1_dummy),ncol=4)
#race_v1_0<-race_v1$Race_0
#race_v1_1<-race_v1$Race_1
#race_v1_2<-race_v1$Race_2
#race_v1_3<-race_v1$Race_3
type_v1<-df_v1$Type


n1<-18

my_lms_v1 <- lapply(1:n1, function(x) lm(x_v1[,x] ~ y_v1))
coeff_v1<-sapply(my_lms_v1, coef)
summaries_v1 <- lapply(my_lms_v1, summary)
p_v1<-lapply(summaries_v1, function(x) x$coefficients[2,4])
p_v1_matrix <- matrix(unlist(p_v1), ncol = 1)
#p_v1_round<-round(p_v1_matrix,digits=2)
p_v1_matrix <-cbind(colnames(x_v1),p_v1_matrix)
p_v1_matrix<-data.frame(p_v1_matrix)
colnames(p_v1_matrix)<-c("names","pvalue")
p_v1_filt<-p_v1_matrix[which(p_v1_matrix[,"pvalue"]<="0.05"),]


my_lms_age_v1<-lapply(1:n1, function(x) lm(x_v1[,x] ~ age_v1))
coeff_v1_age<-sapply(my_lms_age_v1, coef)
summaries_v1_age <- lapply(my_lms_age_v1, summary)
p_v1_age<-lapply(summaries_v1_age, function(x) x$coefficients[2,4])
p_v1_age_matrix <- matrix(unlist(p_v1_age), ncol = 1)
p_v1_age_matrix <-cbind(colnames(x_v1),p_v1_age_matrix)
p_v1_age_matrix<-data.frame(p_v1_age_matrix)
colnames(p_v1_age_matrix)<-c("names","pvalue")
p_v1_age_filt<-p_v1_age_matrix[which(p_v1_age_matrix[,"pvalue"]<="0.05"),]


my_lms_education_v1<-lapply(1:n1, function(x) lm(x_v1[,x] ~ education_v1))
coeff_v1_education<-sapply(my_lms_education_v1, coef)
summaries_v1_education <- lapply(my_lms_education_v1, summary)
p_v1_education<-lapply(summaries_v1_education, function(x) x$coefficients[2,4])
p_v1_education_matrix <- matrix(unlist(p_v1_education), ncol = 1)
p_v1_education_matrix <-cbind(colnames(x_v1),p_v1_education_matrix)
p_v1_education_matrix<-data.frame(p_v1_education_matrix)
colnames(p_v1_education_matrix)<-c("names","pvalue")
p_v1_education_filt<-p_v1_education_matrix[which(p_v1_education_matrix[,"pvalue"]<="0.05"),]


my_lms_ethnicity_v1<-lapply(1:n1, function(x) lm(x_v1[,x] ~ ethnicity_v1))
coeff_v1_ethnicity<-sapply(my_lms_ethnicity_v1, coef)
summaries_v1_ethnicity <- lapply(my_lms_ethnicity_v1, summary)
p_v1_ethnicity<-lapply(summaries_v1_ethnicity, function(x) x$coefficients[2,4])
p_v1_ethnicity_matrix <- matrix(unlist(p_v1_ethnicity), ncol = 1)
p_v1_ethnicity_matrix <-cbind(colnames(x_v1),p_v1_ethnicity_matrix )
p_v1_ethnicity_matrix <-data.frame(p_v1_ethnicity_matrix )
colnames(p_v1_ethnicity_matrix)<-c("names","pvalue")
p_v1_ethnicity_filt<-p_v1_ethnicity_matrix[which(p_v1_ethnicity_matrix[,"pvalue"]<="0.05"),]

my_lms_race_v1<-lapply(1:n1, function(x) lm(x_v1[,x] ~ race_v1))
coeff_v1_race<-sapply(my_lms_race_v1, coef)
summaries_v1_race<- lapply(my_lms_race_v1, summary)
p_v1_race<-lapply(summaries_v1_race, function(x) x$coefficients[2,4])
p_v1_race_matrix <- matrix(unlist(p_v1_race), ncol = 1)
p_v1_race_matrix <-cbind(colnames(x_v1),p_v1_race_matrix)
p_v1_race_matrix<-data.frame(p_v1_race_matrix)
colnames(p_v1_race_matrix)<-c("names","pvalue")
p_v1_race_filt<-p_v1_race_matrix[which(p_v1_race_matrix[,"pvalue"]<="0.05"),]

#For visit V2

df_new_2<-df[(df$Visit=="V2"),]
df_new_2<-df_new_2 %>%select_if(~ !any(is.na(.)))

df_v2<-df[(df$Visit=="V2"),]
df_v2<-df_v2 %>%select_if(~ !any(is.na(.)))


x_v2<-df_v2 %>% select(ITAC_pg_ml:TNFa_pg_ml)
#x_v2_log<-log(x_v2)
y_v2<-df_v2$GAD_7
age_v2<-df_v2$Age
education_v2<-df_v2$Years.of.Education
#ethnicity_v2<-df_v2$Ethnicity
ethnicity_v2<-df_v2%>%select(Ethnicity_0:Ethnicity_1)
ethnicity_v2<-matrix(unlist(ethnicity_v2),ncol = 2)

#race_v2<-df_v2$Race
race_v2<-df_v2%>%select(Race_0:Race_3)
race_v2<-matrix(unlist(race_v2),ncol = 4)
type_v2<-df_v2$Type
#education_v2<-df_v2$Years.of.Education
#education_v2<-df_v2%>%select(Years.of.Education_13:Years.of.Education_20)
#education_v2_13<-education_v2$Years.of.Education_13
#education_v2_16<-education_v2$Years.of.Education_16
#education_v2_17<-education_v2$Years.of.Education_17
#education_v2_18<-education_v2$Years.of.Education_18
#education_v2_19<-education_v2$Years.of.Education_19
#education_v2_20<-education_v2$Years.of.Education_20
#
#ethnicity_v2_0<-ethnicity_v2$Ethnicity_0
#ethnicity_v2_1<-ethnicity_v2$Ethnicity_1
#race_v2_0<-race_v2$Race_0
#race_v2_1<-race_v2$Race_1
#race_v2_2<-race_v2$Race_2
#race_v2_3<-race_v2$Race_3


n2<-18

my_lms_v2<- lapply(1:n2, function(x) lm(x_v2[,x] ~ y_v2))
coeff_v2<-sapply(my_lms_v2, coef)
summaries_v2<- lapply(my_lms_v2, summary)
p_v2<-lapply(summaries_v2, function(x) x$coefficients[2,4])
p_v2_matrix <- matrix(unlist(p_v2), ncol = 1)
p_v2_matrix <-cbind(colnames(x_v2),p_v2_matrix)
p_v2_matrix<-data.frame(p_v2_matrix)
colnames(p_v2_matrix)<-c("names","pvalue")
p_v2_filt<-p_v2_matrix[which(p_v2_matrix[,"pvalue"]<="0.05"),]

my_lms_age_v2<-lapply(1:n2, function(x) lm(x_v2[,x] ~ age_v2))
coeff_v2_age<-sapply(my_lms_age_v2, coef)
summaries_v2_age <- lapply(my_lms_age_v2, summary)
p_v2_age<-lapply(summaries_v2_age, function(x) x$coefficients[2,4])
p_v2_age_matrix <- matrix(unlist(p_v2_age), ncol = 1)
p_v2_age_matrix <-cbind(colnames(x_v2),p_v2_age_matrix)
p_v2_age_matrix<-data.frame(p_v2_age_matrix)
colnames(p_v2_age_matrix)<-c("names","pvalue")
p_v2_age_filt<-p_v2_age_matrix[which(p_v2_age_matrix[,"pvalue"]<="0.05"),]

my_lms_education_v2<-lapply(1:n2, function(x) lm(x_v2[,x] ~ education_v2))
coeff_v2_education<-sapply(my_lms_education_v2, coef)
summaries_v2_education <- lapply(my_lms_education_v2, summary)
p_v2_education<-lapply(summaries_v2_education, function(x) x$coefficients[2,4])
p_v2_education_matrix <- matrix(unlist(p_v2_education), ncol = 1)
p_v2_education_matrix <-cbind(colnames(x_v2),p_v2_education_matrix)
p_v2_education_matrix<-data.frame(p_v2_education_matrix)
colnames(p_v2_education_matrix)<-c("names","pvalue")
p_v2_education_filt<-p_v2_education_matrix[which(p_v2_education_matrix[,"pvalue"]<="0.05"),]

my_lms_ethnicity_v2<-lapply(1:n2, function(x) lm(x_v2[,x] ~ ethnicity_v2))
coeff_v2_ethnicity<-sapply(my_lms_ethnicity_v2, coef)
summaries_v2_ethnicity <- lapply(my_lms_ethnicity_v2, summary)
p_v2_ethnicity<-lapply(summaries_v2_ethnicity, function(x) x$coefficients[2,4])
p_v2_ethnicity_matrix <- matrix(unlist(p_v2_ethnicity), ncol = 1)
p_v2_ethnicity_matrix <-cbind(colnames(x_v2),p_v2_ethnicity_matrix )
p_v2_ethnicity_matrix <-data.frame(p_v2_ethnicity_matrix )
colnames(p_v2_ethnicity_matrix)<-c("names","pvalue")
p_v2_ethnicity_filt<-p_v2_ethnicity_matrix[which(p_v2_ethnicity_matrix[,"pvalue"]<="0.05"),]

my_lms_race_v2<-lapply(1:n2, function(x) lm(x_v2[,x] ~ race_v2))
coeff_v2_race<-sapply(my_lms_race_v2, coef)
summaries_v2_race<- lapply(my_lms_race_v2, summary)
p_v2_race<-lapply(summaries_v2_race, function(x) x$coefficients[2,4])
p_v2_race_matrix <- matrix(unlist(p_v2_race), ncol = 1)
p_v2_race_matrix <-cbind(colnames(x_v2),p_v2_race_matrix)
p_v2_race_matrix<-data.frame(p_v2_race_matrix)
colnames(p_v2_race_matrix)<-c("names","pvalue")
p_v2_race_filt<-p_v2_race_matrix[which(p_v2_race_matrix[,"pvalue"]<="0.05"),]

#For visit 4

df_new_4<-df[(df$Visit=="V4"),]
df_new_4<-df_new_4 %>%select_if(~ !any(is.na(.)))


df_v4<-df[(df$Visit=="V4"),]
df_v4<-df_v4 %>%select_if(~ !any(is.na(.)))

x_v4<-df_v4 %>% select(ITAC_pg_ml:TNFa_pg_ml)
#x_v4_log<-log(x_v4)
y_v4<-df_v4$GAD_7
age_v4<-df_v4$Age
education_v4<-df_v4$Years.of.Education
#ethnicity_v4<-df_v4$Ethnicity
ethnicity_v4<-df_v4%>%select(Ethnicity_0:Ethnicity_1)
ethnicity_v4<-matrix(unlist(ethnicity_v4),ncol=2)
#race_v4<-df_v4$Race
race_v4<-df_v4%>%select(Race_0:Race_3)
race_v4<-matrix(unlist(race_v4),ncol=4)
type_v4<-df_v4$Type
#education_v4<-df_v4$Years.of.Education
#education_v4<-df_v4%>%select(Years.of.Education_13:Years.of.Education_20)
#education_v4_13<-education_v4$Years.of.Education_13
#education_v4_16<-education_v4$Years.of.Education_16
#education_v4_17<-education_v4$Years.of.Education_17
#education_v4_18<-education_v4$Years.of.Education_18
#education_v4_19<-education_v4$Years.of.Education_19
#education_v4_20<-education_v4$Years.of.Education_20

#ethnicity_v4_0<-ethnicity_v4$Ethnicity_0
#ethnicity_v4_1<-ethnicity_v4$Ethnicity_1

#race_v4_0<-race_v4$Race_0
#race_v4_1<-race_v4$Race_1
#race_v4_2<-race_v4$Race_2
#race_v4_3<-race_v4$Race_3


n4<-17

my_lms_v4<- lapply(1:n4, function(x) lm(x_v4[,x] ~ y_v4))
coeff_v4<-sapply(my_lms_v4, coef)
summaries_v4 <- lapply(my_lms_v4, summary)
p_v4<-lapply(summaries_v4, function(x) x$coefficients[2,4])
p_v4_matrix <- matrix(unlist(p_v4), ncol = 1)
#p_v4_round<-formatC(p_v4_matrix, digits = 2, format ="fg")
p_v4_matrix <-cbind(colnames(x_v4),p_v4_matrix)
p_v4_matrix<-data.frame(p_v4_matrix)
colnames(p_v4_matrix)<-c("names","pvalue")
p_v4_filt<-p_v4_matrix[which(p_v4_matrix[,"pvalue"]<="0.05"),]

my_lms_age_v4<-lapply(1:n4, function(x) lm(x_v4[,x] ~ age_v4))
coeff_v4_age<-sapply(my_lms_age_v4, coef)
summaries_v4_age <- lapply(my_lms_age_v4, summary)
p_v4_age<-lapply(summaries_v4_age, function(x) x$coefficients[2,4])
p_v4_age_matrix <- matrix(unlist(p_v4_age), ncol = 1)
p_v4_age_matrix <-cbind(colnames(x_v4),p_v4_age_matrix)
p_v4_age_matrix<-data.frame(p_v4_age_matrix)
colnames(p_v4_age_matrix)<-c("names","pvalue")
p_v4_age_filt<-p_v4_age_matrix[which(p_v4_age_matrix[,"pvalue"]<="0.05"),]

my_lms_education_v4<-lapply(1:n4, function(x) lm(x_v4[,x] ~ education_v4))
coeff_v4_education<-sapply(my_lms_education_v4, coef)
summaries_v4_education <- lapply(my_lms_education_v4, summary)
p_v4_education<-lapply(summaries_v4_education, function(x) x$coefficients[2,4])
p_v4_education_matrix <- matrix(unlist(p_v4_education), ncol = 1)
p_v4_education_matrix <-cbind(colnames(x_v4),p_v4_education_matrix)
p_v4_education_matrix<-data.frame(p_v4_education_matrix)
colnames(p_v4_education_matrix)<-c("names","pvalue")
p_v4_education_filt<-p_v4_education_matrix[which(p_v4_education_matrix[,"pvalue"]<="0.05"),]

my_lms_ethnicity_v4<-lapply(1:n4, function(x) lm(x_v4[,x] ~ ethnicity_v4))
coeff_v4_ethnicity<-sapply(my_lms_ethnicity_v4, coef)
summaries_v4_ethnicity <- lapply(my_lms_ethnicity_v4, summary)
p_v4_ethnicity<-lapply(summaries_v4_ethnicity, function(x) x$coefficients[2,4])
p_v4_ethnicity_matrix <- matrix(unlist(p_v4_ethnicity), ncol = 1)
p_v4_ethnicity_matrix <-cbind(colnames(x_v4),p_v4_ethnicity_matrix )
p_v4_ethnicity_matrix <-data.frame(p_v4_ethnicity_matrix )
colnames(p_v4_ethnicity_matrix)<-c("names","pvalue")
p_v4_ethnicity_filt<-p_v4_ethnicity_matrix[which(p_v4_ethnicity_matrix[,"pvalue"]<="0.05"),]

my_lms_race_v4<-lapply(1:n4, function(x) lm(x_v4[,x] ~ race_v4))
coeff_v4_race<-sapply(my_lms_race_v4, coef)
summaries_v4_race<- lapply(my_lms_race_v4, summary)
p_v4_race<-lapply(summaries_v4_race, function(x) x$coefficients[2,4])
p_v4_race_matrix <- matrix(unlist(p_v4_race), ncol = 1)
p_v4_race_matrix <-cbind(colnames(x_v4),p_v4_race_matrix)
p_v4_race_matrix<-data.frame(p_v4_race_matrix)
colnames(p_v4_race_matrix)<-c("names","pvalue")
p_v4_race_filt<-p_v4_race_matrix[which(p_v4_race_matrix[,"pvalue"]<="0.05"),]

#For 1st visit

df_new_1<-as.data.frame(df_new_1)
Race1<-df_new_1$Race
IL_6_pg_ml<-df_new_1$IL6_pg_ml
IL_8_pg_ml<-df_new_1$IL8_pg_ml

t1<-ggplot(df_new_1, aes(x = factor(Race1),y=IL_6_pg_ml,color=factor(Race1))) + geom_boxplot(width=0.5,lwd=1)+geom_jitter(width=0.3)+scale_x_discrete(name = "Race")+scale_y_continuous(name = "IL6(pg/ml)")
#t1<-t1+legend("bottomright",c("White", "Black or African American","American Indian or Alaska Native","Asian"),cex=.8)
plot(t1)

t2<-ggplot(df_new_1, aes(x = factor(Race1),y=IL_8_pg_ml,color=factor(Race1))) + geom_boxplot(width=0.5,lwd=1)+geom_jitter(width=0.3)+scale_x_discrete(name = "Race")+scale_y_continuous(name = "IL8(pg/ml)")
plot(t2)

#For 2nd visit

h3<-df_v2$IL5_pg_ml# with respect to GADscore

t3<-ggplot(df_v2, aes(x= y_v2, y=h3))+geom_point(size=3)+geom_smooth(method = "lm",col="firebrick",size=3)+xlab("GAD7 score")+ylab("IL5(pg/ml)")
plot(t3)

df_new_2<-as.data.frame(df_new_2)
Race2<-df_new_2$Race
IL_8_pg_ml<-df_new_2$IL8_pg_ml#RACE
TNFa_pg_ml<-df_new_2$TNFa_pg_ml#RACE


t4<-ggplot(df_new_2, aes(x=factor(Race2), y=IL_8_pg_ml,color=factor(Race2))) + geom_boxplot(width=0.5,lwd=1)+geom_jitter(width=0.3)+scale_x_discrete(name = "Race")+scale_y_continuous(name = "IL8(pg/ml)")
#t4<-t4+guides(fill=guide_legend(title="Type"))
plot(t4)


t5<-ggplot(df_new_2, aes(x=factor(Race2), y=TNFa_pg_ml,color=factor(Race2))) + geom_boxplot(width=0.5,lwd=1)+scale_x_discrete(name = "Race")+geom_jitter(width=0.3)+scale_y_continuous(name = "TNFa(pg/ml)")
#t5<-t5+guides(fill=guide_legend(title="Type"))
plot(t5)

#For Postpartum visit

df_new_4<-as.data.frame(df_new_4)
education<-df_new_4$Years.of.Education
Race4<-df_new_4$Race
ITAC_pg_ml<-df_new_4$ITAC_pg_ml #education
IL13_pg_ml<-df_new$IL13_pg_ml#race

t6<-ggplot(df_new_4, aes(x=factor(education_v4), y=ITAC_pg_ml)) + geom_point(size=3)+scale_x_discrete(name = "Education")+scale_y_continuous(name = "ITAC(pg/ml)")
#t6<-t6+guides(fill=guide_legend(title="Type"))
plot(t6)


t7<-ggplot(df_new_4, aes(x=factor(Race4), y=IL13_pg_ml,color=factor(Race4))) + geom_boxplot(width=0.5,lwd=1)+geom_jitter(width=0.3)+scale_x_discrete(name = "Race")+scale_y_continuous(name = "IL13(pg/ml)")
plot(t7)

#####################################################################################################################
#Correlation matrix for cytokines

res<-cor(x_matrix)

#Removed some cytokines since they have missing values in the correlation matrix

x_matrix_new<-subset(x_matrix,select=-c(IL1b_pg_ml,IL2_pg_ml,MIP1a_pg_ml,MIP1b_pg_ml)) 

res <- cor(x_matrix_new) # using the default pearson method
res<-round(res, 2)

#Heatmap for cytokines

pheatmap(res)

res1<-corr()