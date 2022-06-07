##the purpose of this script is to identify variants that were 
#shared by a lot of cats,
##or reached >50% frequency,
##but weren't in the inoculum at >0.1%

#load libraries
library(tidyverse)
library(openxlsx)
library(readxl)
library(DT)

#starting with the 3% cutoff variant table processed and iVar filtered
df.all <- read_xlsx("/Users/_lbashor/Dropbox/SARS-CoV-2 cat manuscript/cat ms results/R analysis cats/2_ivar_cats/variant_summary_0.03_iVar_processed.xlsx")

#lets grab out all VOC variants (9 of them)
df.VOC <- df.all %>%
  filter(VOC=="VOC") 

#now let's get shared variants
##first, filter out any variants detected in cell culture viral stocks at >3%
#keep all data where Passage 1, 2 and 3 viral stock columns have NAs
df.filtered <-
  df.all %>%
  filter_at(vars(starts_with("Passage")), all_vars(is.na(.)))

#count how many cats each variant was found in-- 
#filter out cell culture viral stock columns
#I made the cutoff that a variant is shared among at least 3 cats
#(most of the ones shared by 2 cats are contact variants)
df.shared <- df.filtered %>%
  select(!c("Passage_1", "Passage_2", "Passage_3")) %>%
  pivot_longer(!c(reference_sequence, position, gene, indel,variant, 
                  reference_base, variant_base, effect, VOC), 
               names_to= "dataset_ID", 
               values_to ="frequency") %>%
  na.omit(df.filtered) %>%
  group_by(variant, gene) %>%
  count() %>%
  arrange(desc(n)) %>%
  filter(n>2)

head(df.shared, 10)

#how many cats was each VOC found in?
df.VOC.shared <- df.all %>%
  select(!c("Passage_1", "Passage_2", "Passage_3")) %>%
  pivot_longer(!c(reference_sequence, position, gene, indel,variant, 
                  reference_base, variant_base, effect, VOC), 
               names_to= "dataset_ID", 
               values_to ="frequency") %>%
  na.omit(df.filtered) %>%
  group_by(variant, gene) %>%
  count() %>%
  arrange(desc(n)) %>%
  filter(variant %in% df.VOC$variant)

#so now we've got 9 VOC variants in df.VOC, 
#and 9 shared variants (not found in the inoculum at >3%) in df.shared

#our final emergent variant type is variants reaching >50%
#(again not found in the inoc at >3%)

#time to pull out variants greater than 50% (specifically for variants not detected above 0.1%)
df05 <- df.filtered %>%
  pivot_longer(!c(reference_sequence, position, gene, indel, variant, 
                  reference_base, variant_base, effect, VOC), 
               names_to= "dataset_ID", 
               values_to ="frequency") %>%
  filter(frequency >= 0.499) %>%
  pivot_wider(id_cols = c(reference_sequence, position, gene, indel, variant, 
                          reference_base, variant_base, effect, VOC), 
              names_from=dataset_ID, 
              values_from=c(frequency),
              names_sort = T)

#so here we have 13 variants found at >=50% and not found in the inoc at >3%

#finally, we want to combine df.VOC, df.shared and df.0.5 
#to make one emergent variant table
#bc these dfs might not include every observation of a given variant, 
#I'll pull out just the variant names first
vectorVOC <- df.VOC$variant
vectorshared <- df.shared$variant
vector05 <- df05$variant

#December 2021 what about omicron!?!?
#we had some crossover variants here-- a new one is P13S, P13L in omicron
omicron <- c("G204_209del", "P13S", "H69R", "D614G", "H655Y", "E484D")

emergent_vector <- c(vectorVOC, vectorshared, vector05, omicron)

#only keep unique variant names
emergent_vector <- unique(emergent_vector)

emergent_variants <- df.all %>%
  ungroup() %>%
  filter(variant %in% emergent_vector) %>%
  select(-c("reference_sequence", "indel",  "VOC",
            "Passage_1", "Passage_2", "Passage_3")) %>%
  relocate(variant) %>%
  arrange(-desc(position)) %>%
  relocate(c(Cat_5, Cat_6, Cat_7, 
             Cat_8, Cat_9, Cat_10), .after=Cat_1)

#create an interactive data table for exploration
datatable(emergent_variants, options = list(pageLength = 30))

#now lets make this an excel sheet
wb <- createWorkbook("emergent_variants.xlsx")
addWorksheet(wb, "emergent_variants")
writeData(wb, "emergent_variants", emergent_variants, borders="all")
saveWorkbook(wb, "emergent_variants.xlsx", overwrite = TRUE)
