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
# three variables: drug, concentration in nM, TGFb Rx
#create plater.csv file of plate plan
#read in plater.csv
plate<- read_plate("plater.csv")


#Read in excel file of results
df<- read_excel("TGF-9-2-16.xlsx", sheet=1)
df<- df[,1:6] #select columns with raw data results (not calculated)
names(df)<- c("Wells","label","gene","num","thresh","Ct") #rename
options(tibble.print_max=Inf)
df<-filter(df, !is.na(label)) #remove empty rows
df$Ct<- as.numeric(df$Ct) #make Ct numeric
df$num<-str_replace(df$num,"TGF-9-2-16_","") # clean up numbers

df$Wells[as.numeric(str_sub(df$Wells,2,3))<10]<-str_replace(df$Wells,str_sub(df$Wells,2,2), paste0("0",str_sub(df$Wells,2,2)))[as.numeric(str_sub(df$Wells,2,3))<10]

full<-full_join(df,plate) %>% select(-label,-num)
