#the purpose of this script is to look at variant allele frequencies for the first and second technical sequencing replicates

#load libraries
library(tidyverse)
library(openxlsx)
library(readxl)

#import variant table (3% cutoff)
df <- read_xlsx("/Users/_lbashor/Dropbox/SARS-CoV-2 cat manuscript/cat ms results/Pipeline results 7.13.21/variant_summary_7.13.21.xlsx")

#make data long and skinny
df2 <- df %>%
  pivot_longer(!c(reference_sequence, position, gene, codon,indel,variant, reference_base, variant_base, effect,featureid), 
               names_to= "dataset_ID", 
               values_to ="frequency")

#make a new column labeling replicates (anything with _R is replicate 2, without it gets a 1)
df2$replicate = str_extract(df2$dataset_ID, "_R")
df2$replicate[is.na(df2$replicate)] <- 1
df2$replicate[df2$replicate == "_R"] <- 2

#delete _R from dataset_ID names
df2$dataset_ID <- gsub("_R", "", df2$dataset_ID)

#now make it wide and have a column for rep1 AF and rep2 AF
df_wide <- df2 %>% pivot_wider(id_cols = c(reference_sequence, position, gene, indel,variant, reference_base, variant_base, effect, dataset_ID), 
                               names_from=replicate,
                               names_prefix="rep_",
                               values_from=frequency,
                               names_sort = T)



#now plot the replicates

reps <- ggplot(df_wide, aes(x=rep_1, y=rep_2))+
  geom_point()+
  geom_vline(xintercept = 0.03, color="blue")+
  geom_hline(yintercept = 0.03, color="blue")+
  labs(x="Replicate 1", y="Replicate 2")+
  theme_classic()

print(reps)

pdf("replicate_correlation.pdf")
reps
dev.off()

#you can use geom_text(label=variant) to 
#identify the really disparate replicates, which include S650F and T7I for one cat
#or just label one variant to see what its up to

ggplot(df_wide, aes(x=rep_1, y=rep_2, label=variant))+
  geom_point()+
  geom_text()+
  geom_vline(xintercept = 0.03, color="blue")+
  geom_hline(yintercept = 0.03, color="blue")+
  labs(x="Replicate 1", y="Replicate 2")+
  theme_classic()
