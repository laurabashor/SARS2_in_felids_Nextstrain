#is sequencing coverage related to the number of variants or the viral load?
#seq coverage was calculated as mean/median of both tech. rep.s of a sample

#RESULT:
#no significant relationship between viral load and sequencing coverage
#no significant relationship between seq coverage and #variants

#load packages
library(tidyverse)
library(ggforce)
library(RColorBrewer)
library(ggpubr)
library(readxl)
library(readr)
library(broom)

#start with data frame of cleaned, combined replicates
df <- read_xlsx("/Users/_lbashor/Dropbox/SARS-CoV-2 cat manuscript/cat ms results/R analysis cats/2_ivar_cats/variant_summary_0.03_iVar_processed.xlsx")
df_meta <- read_csv("/Users/_lbashor/Dropbox/SARS-CoV-2 cat manuscript/cat ms results/metadata/cat_metadata.csv")
df_depth <- read_csv("Median_and_mean_depths.csv")

df_depth <- df_depth %>%
  rename(dataset_ID = cat_ID) %>%
  select(-virus)

#make data long and skinny and label animal species and infection method (contact cats)
df2 <- df %>%
  pivot_longer(!c(reference_sequence, VOC, position, gene, indel,variant, reference_base, variant_base, effect), 
               names_to= "dataset_ID", values_to ="frequency") %>%
  na.omit(df2) %>%
  mutate(infection_method = ifelse(dataset_ID %in% c("Cat_6","Cat_10", "Cat_11", "Cat_12", "Cat_26"), 
                                   "contact", 
                                   ifelse(dataset_ID %in% c("Passage_1", "Passage_2", "Passage_3"), 
                                          "inoculum", "direct_inoculation"))) %>%
  mutate(species = case_when(grepl("C", dataset_ID) ~ "Cats",
                             grepl("P", dataset_ID) ~"Vero"))

#how many variants were detected in each cat?
df_all <- df2 %>% 
  group_by(species, dataset_ID, infection_method) %>%
  summarize(number_of_variants=n()) 

df_all <- inner_join(df_all, df_depth)

ggplot(df_all, aes(x=median_depth, y=number_of_variants)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()

median_model <- lm(number_of_variants ~ median_depth, df_all)
glance(median_model) #p = 0.928, R2=0.00035

ggplot(df_all, aes(x=mean_depth, y=number_of_variants)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()

mean_model <- lm(number_of_variants ~ mean_depth, df_all)
glance(mean_model) #p = 0.590, R2=0.012


##what about viral titer? here we just need cats

df_cats <- df_all %>%
  filter(species == "Cats")

#add in the metadata 
df_cats <- merge(df_cats, df_meta)

#look at titer and depth
ggplot(df_cats, aes(x=log10(nasal_pfu + 0.001), y=median_depth)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()

median_pfu_model <- lm(median_depth ~ log10(nasal_pfu + 0.001), df_cats)
glance(median_pfu_model) #p = 0.560, R2 = 0.016

ggplot(df_cats, aes(x=log10(nasal_pfu + 0.001), y=mean_depth)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()

mean_pfu_model <- lm(mean_depth ~ log10(nasal_pfu + 0.001), df_cats)
glance(mean_pfu_model) #p = 0.0764, R2 = 0.142, almost lower coverage with higher pfu!

##so nothing is significant
