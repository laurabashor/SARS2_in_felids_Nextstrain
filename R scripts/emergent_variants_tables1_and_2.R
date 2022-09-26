### identify emergent variants, including:

#(1) variants not detected at any level in viral stocks (>0.1%) that reached 50% AF
# or higher or were at positions of interest (VOC) >>> table 1 in ms

# (2) everything else interesting: variants detected in viral stocks but not above 3% (so they started very low) 
# and increased to >50% in cats, or were found in 3 or more cats, or were 
# at positions of interest (VOC) >>> table 2 in ms

#load libraries
library(tidyverse)
library(openxlsx)
library(readxl)

# table 1 ############################### 


# #start with a processed 0.001 (0.1%) variant table
df01 <- read_xlsx("/Users/_lbashor/Dropbox/SARS-CoV-2 cat manuscript/cat ms results/R analysis cats/1_initial_processing_0.001_cats/variant_summary_0.001_processed.xlsx")


# filter for variants that were not detected at all in the inoculum viral stocks
# filter to keep all data where Passage 1, 2 and 3 viral stock columns have NAs

df.filtered <-
  df01 %>%
  filter_at(vars(starts_with("Passage")), all_vars(is.na(.))) %>%
  select(-c("Passage_1", "Passage_2", "Passage_3"))

# now df.filtered has all variants not detected in viral stocks at any level (>0.1%)

# pull out VOC that reached >3%
vectorVOC <- df.filtered %>%
  filter(VOC=="VOC") %>%
  pivot_longer(!c(reference_sequence, position, gene, indel, variant, 
                  reference_base, variant_base, effect, VOC), 
               names_to= "dataset_ID", 
               values_to ="frequency") %>%
  filter(frequency >= 0.0299) %>%
  pull(variant)

#also prep vector of relevant omicron variants
omicron <- c("G204_209del", "P13S", "H69R", "D614G", "H655Y", "E484D")

# pull out variants that reached >50%:
vector50 <- df.filtered %>%
  pivot_longer(!c(reference_sequence, position, gene, indel, variant, 
                  reference_base, variant_base, effect, VOC), 
               names_to= "dataset_ID", 
               values_to ="frequency") %>%
  filter(frequency >= 0.499) %>%
  pivot_wider(id_cols = c(reference_sequence, position, gene, indel, variant, 
                          reference_base, variant_base, effect, VOC), 
              names_from=dataset_ID, 
              values_from=c(frequency),
              names_sort = T) %>%
  pull(variant)

#finally go back to original table and pull out all these variants
#(using vectors with variant names because the filtered dfs might not include every observation of a given variant)
combined <- c(vectorVOC, vector50, omicron) %>%
  unique()

# now make table 1:

table1 <- df.filtered %>%
  filter(variant %in% combined) %>%
  select(-c("reference_sequence", "indel",  "VOC")) %>%
  relocate(variant) %>%
  arrange(-desc(position)) %>%
  relocate(c(Cat_5, Cat_6, Cat_7, 
             Cat_8, Cat_9, Cat_10), .after=Cat_1)

#make a list of all these table 1 possible de novo variants:

emergent_denovo_variants <- table1 %>%
  pull(variant)

# table 2 ###############################
# now table 2: pulling out any other variants of interest and seeing how many cats they were found in

#start with the 3% cutoff variant table processed and iVar-filtered
df.all <- read_xlsx("/Users/_lbashor/Dropbox/SARS-CoV-2 cat manuscript/cat ms results/R analysis cats/2_ivar_cats/variant_summary_0.03_iVar_processed.xlsx")

#lets pull out the VOC variants from here
VOC <- df.all %>%
  filter(VOC=="VOC") %>%
  pull(variant)

#next pull out variants that reached >50% in cats but weren't in viral stocks at >3%

#first, filter out any variants detected viral stocks at >3%
#keeping all data where Passage 1, 2 and 3 viral stock columns have NAs
df.filtered2 <-
  df.all %>%
  filter_at(vars(starts_with("Passage")), all_vars(is.na(.)))

#pull out variants greater than 50%
v50 <- df.filtered2 %>%
  pivot_longer(!c(reference_sequence, position, gene, indel, variant, 
                  reference_base, variant_base, effect, VOC), 
               names_to= "dataset_ID", 
               values_to ="frequency") %>%
  filter(frequency >= 0.499) %>%
  pivot_wider(id_cols = c(reference_sequence, position, gene, indel, variant, 
                          reference_base, variant_base, effect, VOC), 
              names_from=dataset_ID, 
              values_from=c(frequency),
              names_sort = T) %>%
  pull(variant)

#now look at variants detected in 3 or more cats that weren't in the inoculum at >3%
#(going with 3 because that means that its more than just contact transmission of a variant between two cats
#(even if it is just contact transmission from the same cat to two others (ie Cat 7 to both Cat 10,11) that still means that the variant was transmitted successfully multiple times)

#first count how many cats each variant was found in, then filter to 3+
shared <- df.filtered2 %>%
  pivot_longer(!c(reference_sequence, position, gene, indel,variant, 
                  reference_base, variant_base, effect, VOC), 
               names_to= "dataset_ID", 
               values_to ="frequency") %>%
  na.omit() %>%
  group_by(variant, gene) %>%
  count() %>%
  arrange(desc(n)) %>%
  filter(n>2) %>%
  pull(variant)

emergent_remaining_variants <- c(VOC, v50, omicron, shared) %>%
  unique()

# remove any variants that are going in table 1 because they might be de novo
emergent_remaining_variants <- emergent_remaining_variants[! emergent_remaining_variants %in% emergent_denovo_variants]

#get all info needed for table 2, starting with df01 so we get all observations down to 0.1% of any of these variants
table2_all <- df01 %>%
  filter(variant %in% emergent_remaining_variants) %>%
  filter(variant != "L37fs", ##this variant removed due to ambiguity
         position != "27379") %>% #can see that one of the D61fs variants was only in 1 cat, remove that one based on its position 
  select(-c("reference_sequence", "indel",  "VOC")) %>%
  relocate(variant) %>%
  arrange(-desc(position)) %>%
  relocate(c(Cat_5, Cat_6, Cat_7, 
             Cat_8, Cat_9, Cat_10), .after=Cat_1)

#get the counts of how many cats it was observed in
cat_counts <- table2_all %>%
  select(-c("Passage_1", "Passage_2", "Passage_3")) %>%
  pivot_longer(!c(position, gene, variant, reference_base, variant_base, effect), 
                           names_to= "dataset_ID", 
                           values_to ="frequency") %>%
  na.omit() %>%
  group_by(variant, gene, position) %>%
  count() %>%
  arrange(position)

table2 <- table2_all %>%
  select(-c("Passage_1", "Passage_2", "Passage_3")) %>%
  pivot_longer(!c(position, gene, variant, reference_base, variant_base, effect), 
               names_to= "dataset_ID", 
               values_to ="frequency") %>%
  na.omit() %>%
  pivot_wider(id_cols = c(position, gene,variant, 
                           reference_base, variant_base, effect),
              names_from=dataset_ID, 
              values_from=c(frequency),
              names_sort = T) %>%
  select(1:3, 6) %>%
  left_join(cat_counts)
  

#now lets export these tables as an excel workbook:

wb <- createWorkbook("emergent_variants20220904.xlsx")
addWorksheet(wb, "table1")
writeData(wb, "table1", table1, borders="all")
addWorksheet(wb, "table2")
writeData(wb, "table2", table2, borders="all")
saveWorkbook(wb, "emergent_variants20220904.xlsx", overwrite = TRUE)