#the purpose of this script is to process the data in the variant table output by
#the viral_variant_caller pipeline for downstream analysis/figures

#load libraries
library(tidyverse)
library(openxlsx)
library(readxl)

#raw output of viral_variant_caller pipeline is variant_summary.xlsx
#rename with the date of the pipeline run
#import variant table

df <- read.xlsx("variant_summary_0.001_cutoff_7.13.21.xlsx")

####combining replicates -- average them, unless one is NA, then don't keep it
#get rid of the noncoding regions, gene = "."
#make data long and skinny
#removing Cats 2, 3, 4 from the mix -- will still leave the cat #s the same, as they do not have a contact cat & are previously published

df2 <- df %>%
  filter(gene!=".") %>%
  mutate(VOC = if_else(grepl("VOC", featureid), "VOC", "not_VOC")) %>%
  select(!c("Cat_2", "Cat_2_R", "Cat_3","Cat_3_R", "Cat_4", "Cat_4_R")) %>%
  pivot_longer(!c(reference_sequence, position, gene, codon,indel,variant, 
                  reference_base, variant_base, effect, featureid, VOC), 
               names_to= "dataset_ID", 
               values_to ="frequency")

#make a new column labeling replicates (anything with _R is replicate 2, without it gets a 1)
df2$replicate = str_extract(df2$dataset_ID, "_R")
df2$replicate[is.na(df2$replicate)] <- 1
df2$replicate[df2$replicate == "_R"] <- 2

#important note: Cat 5 (Cat 70) was originally amplified with ARTIC v2 primers 
#there is no replicate Cat_5_R
#all other samples, both replicates were sequenced with ARTIC v3 primers

#now delete _R from sample names because that data is encoded in replicate column
df2$dataset_ID <- gsub("_R", "", df2$dataset_ID)

#right now, there is a 0 if we didn't detect the variant at >3% frequency, and an NA if we didn't have enough read coverage to know if we actually didn't detect it, or it didn't show up due to low coverage 
#for the purpose of combining replicates, we will treat the 0s and NAs the same, as NA
df2[df2 == "0"] <- NA

#what we want is to throw out any data where we didn't detect the variant in both replicates, whatever the reason for this may be
#calculate mean of variants, if there is an NA in either value, the mean function will always output NA 
df.all <- df2 %>%
 group_by(reference_sequence, position, gene,codon, indel, variant, 
          reference_base, variant_base, effect, VOC, featureid, dataset_ID) %>%
  summarise(mean = mean(frequency)) %>% 
  pivot_wider(id_cols = c(reference_sequence, position, gene, indel,variant, 
                          reference_base, variant_base, effect, VOC), 
                              names_from=dataset_ID, 
                              values_from=c(mean),
                              names_sort = T)

#I noticed that snpeff/snpsift have added some random commas and periods into gene names and variant names
df.all$gene <- gsub(",.", "", df.all$gene)
df.all$variant <- gsub(",.", "", df.all$variant)

#get rid of variant rows that just have NAs (bc we didn't detect the variant in both replicates)
delete.na <- function(df.all, n=0) {
  df.all[rowSums(is.na(df.all)) <= n,]
}
df.all <- delete.na(df.all,25) #25 is max NAs allowed for 26 datasets (23 cats + 3 stock passages)

#for our final excel spreadsheet, it will be nice to be able to look at variants in cats and in viral stocks separately
#subset into data frames for cats and viral stocks
df.cats <- df.all %>% 
  select(reference_sequence,position, gene, indel,variant, reference_base, variant_base, effect, starts_with("C"))

df.inoc <- df.all %>%
  select(reference_sequence,position, gene, indel,variant, reference_base, variant_base, effect,starts_with("P"))

#use a function to delete any rows in which there are NAs for all cats
delete.na.cats <- function(df.cats, n=0) {
  df.cats[rowSums(is.na(df.cats)) <= n,]
}
df.cats <- delete.na.cats(df.cats, 22) #22 is max NAs allowed for 23 cats

#and do it for inoculums
delete.na.inoc <- function(df.inoc, n=0) {
  df.inoc[rowSums(is.na(df.inoc)) <= n,]
}
df.inoc <- delete.na.inoc(df.inoc, 2) #2 is max NAs allowed for 3 passages

##now lets filter for variants that were not detected at all in the inoculum viral stocks
#filter to keep all data where Passage 1, 2 and 3 viral stock columns have NAs, and name this df.all
df.filtered <-
  df.all %>%
  filter_at(vars(starts_with("Passage")), all_vars(is.na(.)))

#so there were 610 variants found across all datasets (df.all)
#557 of these were found in cats, and 170 of them were in the viral stocks
#and there were 440 variants found in cats that were not found in the inoculum

df.filtered %>%
  group_by(indel) %>%
  count() #280 SNVs, 160 indels

df.filtered %>%
  filter(indel==F) %>%
  group_by(effect) %>%
  count() #103 synonymous, 177 NS

#save variant table as excel sheet, with sheets for each species
wb <- createWorkbook("variant_summary_0.001_processed.xlsx")
addWorksheet(wb, "variant_table_all_means")
addWorksheet(wb, "cats")
addWorksheet(wb, "viral_stock_inoculum")
addWorksheet(wb, "all_variants_not_in_inoculum")
writeData(wb, "variant_table_all_means", df.all,borders="all")
writeData(wb, "cats", df.cats,borders="all")
writeData(wb, "viral_stock_inoculum", df.inoc, borders="all")
writeData(wb, "all_variants_not_in_inoculum", df.filtered,borders="all")
saveWorkbook(wb, "variant_summary_0.001_processed.xlsx", overwrite = TRUE)

#now make tables with just contact cats 

#Cohort A: P3 adult cats
df.contact01 <- df.all %>% 
  select(reference_sequence,position, gene, indel,variant, reference_base, 
         variant_base, effect, Passage_3, Cat_5, Cat_6)

delete.na.contact1 <- function(df.contact01, n=0) {
  df.contact01[rowSums(is.na(df.contact01)) <= n,]
}
df.contact01 <- delete.na.contact1(df.contact01, 2)

#Cohort B: P2 adult cats
df.contact02 <- df.all %>% 
  select(reference_sequence,position, gene, indel,variant, reference_base, 
         variant_base, effect, Passage_2, 
         Cat_7, Cat_8, Cat_9, Cat_10, Cat_11, Cat_12)

delete.na.contact2 <- function(df.contact02, n=0) {
  df.contact02[rowSums(is.na(df.contact02)) <= n,]
}
df.contact02 <- delete.na.contact2(df.contact02, 6)

#Cohort C: P3 juvenile cats
df.contact03 <- df.all %>% 
  select(reference_sequence,position, gene, indel,variant, reference_base, 
         variant_base, effect, Passage_3, 
         Cat_22, Cat_23, Cat_24, Cat_25, Cat_26)

delete.na.contact3 <- function(df.contact03, n=0) {
  df.contact03[rowSums(is.na(df.contact03)) <= n,]
}
df.contact03 <- delete.na.contact3(df.contact03, 5)

#save excel sheet with contact cats
wb <- createWorkbook("all_contact_and_cohoused_cats.xlsx")
addWorksheet(wb, "cohort_1")
addWorksheet(wb, "cohort_2")
addWorksheet(wb, "cohort_3")
writeData(wb, "cohort_1", df.contact01,borders="all")
writeData(wb, "cohort_2", df.contact02,borders="all")
writeData(wb, "cohort_3", df.contact03,borders="all")
saveWorkbook(wb, "all_contact_and_cohoused_cats.xlsx", overwrite = TRUE)


##finally, let's make the 0.03 cutoff variant table from this data
#start with df.all and filter to variants with AF>0.0299

df0.03 <- df.all %>%
  pivot_longer(!c(reference_sequence, position, gene,indel,variant,
                  reference_base, variant_base, effect, VOC), 
               names_to= "dataset_ID", 
               values_to ="frequency")  %>%
  filter(frequency>0.0299) %>% 
  pivot_wider(id_cols = c(reference_sequence, position, gene, indel,variant,
                          reference_base, variant_base, effect, VOC), 
                                          names_from=dataset_ID, 
                                          values_from=c(frequency),
                                          names_sort = T)

#count the observations of these 129 variants found at >3%
df.all %>%
  pivot_longer(!c(reference_sequence, position, gene,indel,variant,
                  reference_base, variant_base, effect, VOC), 
               names_to= "dataset_ID", 
               values_to ="frequency")%>%
  filter(frequency>0.0299) %>% 
  filter(!dataset_ID %in% c("Passage_1", "Passage_2", "Passage_3"))%>%
  na.omit(df.all)
#306 observations of the 129 variants found in cats
#most common non-cell culture variant is T556P observed in 15/23 cats

#for a final excel spreadsheet, it will be nice to be able to look at variants in cats and in viral stocks separately
#subset into data frames for cats and viral stocks
df.cats <- df0.03 %>% 
  select(reference_sequence,position, gene, indel,variant, 
         reference_base, variant_base, effect, starts_with("C"))

df.inoc <- df0.03 %>%
  select(reference_sequence,position, gene, indel,variant, 
         reference_base, variant_base, effect,starts_with("P"))

#use a function to delete any rows in which there are NAs for all cats
delete.na.cats <- function(df.cats, n=0) {
  df.cats[rowSums(is.na(df.cats)) <= n,]
}
df.cats <- delete.na.cats(df.cats, 22) #22 is max NAs allowed for 23 cats

#and do it for inoculums
delete.na.inoc <- function(df.inoc, n=0) {
  df.inoc[rowSums(is.na(df.inoc)) <= n,]
}
df.inoc <- delete.na.inoc(df.inoc, 2) #2 is max NAs allowed for 3 passages

##now lets filter for variants that were not detected at all in the inoculum viral stocks
#filter to keep all data where Passage 1, 2 and 3 viral stock columns have NAs
df.filtered <- df0.03 %>%
  filter_at(vars(starts_with("Passage")), all_vars(is.na(.)))

#so final information here is: 129 across all datasets, and 129 in cats
#19 in the inoculums
#110 in cats but not found in the inoculums

#how many were SNVs vs. SVs?
df0.03 %>%
  group_by(effect) %>%
  count() %>%
  arrange(desc(n))
#111 SNVs, 18 SVs
#80 missense (nonsyn SNV) #31 synonymous

#now lets make this an excel sheet
wb <- createWorkbook("variant_summary_0.03_processed.xlsx")
addWorksheet(wb, "all_means")
addWorksheet(wb, "cats")
addWorksheet(wb, "viral_stock_inoculum")
addWorksheet(wb, "all_variants_not_in_inoculum")
writeData(wb, "all_means", df0.03,borders="all")
writeData(wb, "cats", df.cats,borders="all")
writeData(wb, "viral_stock_inoculum", df.inoc, borders="all")
writeData(wb, "all_variants_not_in_inoculum", df.filtered,borders="all")
saveWorkbook(wb, "variant_summary_0.03_processed.xlsx", overwrite = TRUE)

