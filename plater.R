#load libraries
library(tidyverse)
library(readxl)
library(ggplot2)
library(broom)
library(intubate)
library(magrittr)
library(ggbeeswarm)
library(plater)
library(stringr)
library(readr)

#clear environment
rm(list=ls())

#plan experiment
# five variables: gene, drug, concentration in nM, TGFb Rx, template present
#create plater.csv file of plate plan
#read in plater.csv that reflects plate experiment design
plate<- read_plate("plater.csv")


#Read in excel file of results
df<- read_excel("TGF-9-2-16.xlsx", sheet=1)
df<- df[,1:6] #select columns with raw data results (not calculated)
names(df)<- c("Wells","label","gene","num","thresh","Ct") #rename
options(tibble.print_max=Inf)
df<-filter(df, !is.na(label)) #remove empty rows
df$Ct<- as.numeric(df$Ct) #make Ct numeric
df$num<-str_replace(df$num,"TGF-9-2-16_","") # clean up numbers

#change well labels from A1 to A01 so all are 3 characters long, facilitating join
df$Wells[as.numeric(str_sub(df$Wells,2,3))<10]<-str_replace(df$Wells,str_sub(df$Wells,2,2), paste0("0",str_sub(df$Wells,2,2)))[as.numeric(str_sub(df$Wells,2,3))<10]

#now do the full join
full<-full_join(df,plate) %>% select(-label,-num)

#now summarize replicates
full %>%
  group_by(gene, TGFbeta, drug, conc, template) %>% 
  summarize (mCt=mean(Ct))

# try a facet graph
full %>% 
  filter(template==1, gene=="col1A1") %>% 
ggplot(aes(x=as.factor(conc), y=Ct)) + geom_quasirandom() + facet_grid(drug ~ TGFbeta)
