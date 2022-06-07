##the purpose of this script is to plot proportions of the nucleotide substitution types

#load libraries
library(tidyverse)
library(openxlsx)
library(readxl)

#start with cleaned up variant table, iVar processed
df <- read_xlsx("/Users/_lbashor/Dropbox/SARS-CoV-2 cat manuscript/cat ms results/R analysis cats/2_ivar_cats/variant_summary_0.03_iVar_processed.xlsx")

#make data long and skinny, label species, and get rid of insertions/deletions
#just looking at single nucleotide substitutions
df2 <- df %>%
  pivot_longer(!c(reference_sequence, position, gene, indel,variant, 
                  reference_base, variant_base, effect, VOC), 
               names_to= "dataset_ID", values_to ="frequency") %>%
  filter(indel==FALSE) %>%
  mutate(species = case_when(grepl("C", dataset_ID) ~ "Cats",
                             grepl("P", dataset_ID) ~"Vero")) %>%
  na.omit(df2)

#paste together reference and variant bases
df2$substitution <- paste(df2$reference_base, df2$variant_base, sep=">")

#now lets count each type of substitution
as.factor(df2$substitution)

substitution_count <- df2 %>% 
  group_by(substitution, species) %>%
  summarise(count = length(substitution))%>%
  mutate(perc=(count/sum(count)))%>%
  select(substitution, species, count)%>%
  group_by(species)%>%
  mutate(percentage_within_species=count/sum(count)) %>%
  arrange(desc(percentage_within_species))

head(substitution_count)

# C>T 30.4% in Vero and 31.4% in cats

#also want to do this without dealing with species info 
substitution_count_all <- df2 %>% 
  filter(indel==FALSE) %>%
  group_by(substitution) %>%
  summarise(count = length(substitution)) %>%
  mutate(perc=(count/sum(count))) %>%
  arrange(desc(perc))

head(substitution_count_all) #31.3% C>T overall

#plot the substitutions across cat and inoculum samples

ggplot(substitution_count_all, aes(x=reorder(substitution,-(count/sum(count))), y=(count/sum(count)), fill=substitution))+
  geom_bar(stat="identity")+
  scale_y_continuous(labels=scales::percent)+
  labs(x="", y="")+
  scale_fill_brewer(palette="Paired")+
  theme_classic()+
  theme(axis.text=element_text(size=16), 
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18), 
        legend.position="none")

#plot this by species
substitution_count$substitution <- factor(substitution_count$substitution, 
                                          levels= c("C>T", "A>G", "A>C","G>A", 
                                                    "T>A","T>G" ,"T>C","G>C", 
                                                    "G>T","C>A" ,  "A>T"))
levels(substitution_count$substitution)

pdf("substitution_percentages_cats_vero.pdf", width=12, height=8)
ggplot(substitution_count, aes(x=substitution, y=(percentage_within_species), fill=substitution))+
  geom_bar(stat="identity")+
  scale_y_continuous(labels=scales::percent)+
  labs(x= "", y="")+
  scale_fill_brewer(palette="Paired")+
  facet_wrap(~species)+
  theme_classic()+
  theme(axis.text.y=element_text(size=16),
        strip.text=element_text(size=18),
        axis.text.x=element_text(angle=60, vjust=0.5,size=16),
        legend.position="none")
dev.off()

#plot just C>T and compare species

ggplot(data = (substitution_count %>% filter(substitution=="C>T")),
               aes(x=reorder(species,-percentage_within_species), y=percentage_within_species, fill=species))+
  scale_y_continuous(labels=scales::percent)+
  geom_bar(stat="identity")+
  theme_classic()+
  scale_fill_brewer(palette="Paired")+
  labs(y="", x="", title="Percentage of C>T Substitutions")+
  theme(text=element_text(size=18), legend.position="none")
# dev.off()