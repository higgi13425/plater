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
# six variables: gene, drug, concentration in nM, TGFb Rx, template present, samplenum
#create plater.csv file of plate plan
#read in plater.csv that reflects plate experiment design
plate<- read_plate("plater.csv")

#Read in excel file of results
df<- read_excel("TGF-9-2-16.xlsx", sheet="simulated raw file", skip=6)
df<- df[,1:6] #select columns with raw data results (not calculated)
names(df)<- c("Wells","label","gene","samplenum","thresh","Ct") #rename
options(tibble.print_max=Inf)
df<-filter(df, !is.na(label)) #remove empty rows
df$Ct<- as.numeric(df$Ct) #make Ct numeric
df$samplenum<-str_replace(df$samplenum,"TGF-9-2-16_","") # clean up numbers

#change well labels from A1 to A01 so all are 3 characters long, facilitating join
df$Wells[as.numeric(str_sub(df$Wells,2,3))<10]<-str_replace(df$Wells,str_sub(df$Wells,2,2), paste0("0",str_sub(df$Wells,2,2)))[as.numeric(str_sub(df$Wells,2,3))<10]

#confirm NTC were Undetermined
if (sum(df$samplenum=="NTC") == sum(is.na(df$Ct))) {
  cat("All NTC wells are Undetermined. \nValid Plate")
} else {
  (cat("NTC wells with signal. \nInvalid Plate"))
}

# Now remove NTC wells from both files
df<-filter(df, df$samplenum != "NTC")
plate <- filter(plate, plate$template==1)

df$samplenum <- as.integer(df$samplenum) #change to integer to match plater file

#now do the full join
full<-full_join(df,plate) %>% select(-label, -template)

# check that the plate makes sense
view_plate(full, well_ids_column = "Wells", columns_to_display = c("gene", "drug", "conc", "TGFbeta"))

#identify drug list
drug_list <- distinct(full, drug) %>% filter(drug!=0)

# now replicate UN wells for each drug, so that have appropriate 
# drug conc=0 controls for each graph
full_un <- full %>% filter(drug==0, conc==0, TGFbeta==0)

for(x in 1:nrow(drug_list)){
  new<- full_un
  new$drug<-drug_list$drug[x]
  new$Wells <- "new un"
  full<- rbind(full, new)
}

# now replicate UN wells for each drug, so that have appropriate 
# drug conc=0 controls for each graph
full_t <- full %>% filter(drug==0, conc==0, TGFbeta==1)

for(x in 1:nrow(drug_list)){
  new<- full_t
  new$drug<-drug_list$drug[x]
  new$Wells <- "new T"
  full<- rbind(full, new)
}


#now summarize replicates by mean Ct
full %>%
  group_by(gene, TGFbeta, drug, conc) %>% 
  summarize (mCt=mean(Ct))

# try a facet graph
#first label levels of TGFbeta stimulus
full$TGFbeta <- as.factor(full$TGFbeta)
levels(full$TGFbeta) <- c("no stimulus", "TGFbeta 1 mM")

full %>% 
  filter(gene=="col1A1", drug!=0) %>% 
ggplot(aes(x=as.factor(conc), y=Ct)) + geom_quasirandom() + facet_grid(drug ~ as.factor(TGFbeta)) + 
    labs(x="Drug Concentration (nM)") + 
  theme(axis.title = element_text(family="Arial", face="bold", size=16))
