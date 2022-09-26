## counting up variants

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

df.inoc <- df.all %>%
  select(reference_sequence,position, gene, indel,variant, reference_base, variant_base, effect,starts_with("P"))

#and do it for inoculums
delete.na.inoc <- function(df.inoc, n=0) {
  df.inoc[rowSums(is.na(df.inoc)) <= n,]
}
df.inoc <- delete.na.inoc(df.inoc, 2) #2 is max NAs allowed for 3 passages

#there were 610 variants found across all datasets (df.all)
#557 of these were found in cats, and 170 of them were in the viral stocks
#and there were 440 variants found in cats that were not found in the inoculum


#make a 3% cat table to combine with 0.1% inoculum table

df0.03cats <- df.all %>%
  pivot_longer(!c(reference_sequence, position, gene,indel,variant,
                  reference_base, variant_base, effect, VOC), 
               names_to= "dataset_ID", 
               values_to ="frequency")  %>%
  filter(frequency>0.0299) %>% 
  pivot_wider(id_cols = c(reference_sequence, position, gene, indel,variant,
                          reference_base, variant_base, effect, VOC), 
              names_from=dataset_ID, 
              values_from=c(frequency),
              names_sort = T) %>%
 select(reference_sequence,position, gene, indel,variant, reference_base, variant_base, effect, starts_with("C"))

#use a function to delete any rows in which there are NAs for all cats
delete.na.cats <- function(df0.03cats, n=0) {
  df0.03cats[rowSums(is.na(df0.03cats)) <= n,]
}
df0.03cats <- delete.na.cats(df0.03cats, 22) #22 is max NAs allowed for 23 cats


#merge 3% cats with 0.1% inoculum

dfmerged <- merge(df0.03cats, df.inoc) %>%
  select(-c(VOC, reference_sequence, indel)) #for clarity

dfmerged <- inner_join(df0.03cats, df.inoc)

#there are 52 variants at 3% or greater in cats that were present in inoculum at 0.1% or greater
#as compared with 18 variants that were present in inoculum at 3% or greater

#so 66/118 or 56% of variants were found in cats at 3% or greater but were not found in the inoculum at 0.1% or greater
118-52



##now lets filter for variants that were not detected at all in the inoculum viral stocks
#filter to keep all data where Passage 1, 2 and 3 viral stock columns have NAs, and name this df.all
df.filtered <-
  dfmerged %>%
  filter_at(vars(starts_with("Passage")), all_vars(is.na(.)))



#look at Cohort B shared variants
#first find everying in Cats 10, 11 and 12 at >3%, 
#then take those and go down to 0.1% and see what was in Cats, 7,8,9

contact <- df0.03cats %>%
  select(position, variant,
         Cat_10, Cat_11, Cat_12)

delete.na.contact <- function(contact, n=0) {
  contact[rowSums(is.na(contact)) <= n,]
}

contact <- delete.na(contact, 2) #2 is max NAs allowed for 3 cats

variants <- contact$variant

# now go down to 0.1% in the entire cohort B
cohortB <- df.all %>%
  select(position, variant, 
         Passage_2, Cat_7, Cat_8, Cat_9, 
         Cat_10, Cat_11, Cat_12)

delete.na <- function(cohortB, n=0) {
  cohortB[rowSums(is.na(cohortB)) <= n,]
}

cohortB <- delete.na(cohortB, 6) #6 is max NAs allowed for 6 cats and 1 P2


#now compare

cohortB01 <- cohortB %>% filter(variant %in% variants)

#number variants in Cats 10, 11, 12 are 7, 5, 6

# Cat 7: 6/7, 4/5, 3/6
# Cat 8: 4/7, 4/5, 2/6
# Cat 9: 5/7, 3/5, 2/6

