## what is driving variation in the number of variants among cat samples?
#the purpose of this script is to visualize this and do basic counting-- for stats, see number_variants_stats.R script

#load packages
library(tidyverse)
library(ggforce)
library(RColorBrewer)
library(ggpubr)
library(readxl)
library(GGally)

#start with data frame of cleaned data, combined replicates, and iVar-filtered
df <- read_xlsx("/Users/_lbashor/Dropbox/SARS-CoV-2 cat manuscript/cat ms results/R analysis cats/2_ivar_cats/variant_summary_0.03_iVar_processed.xlsx")

#load in additional metadata
df_meta <- read_csv("/Users/_lbashor/Dropbox/SARS-CoV-2 cat manuscript/cat ms results/metadata/cat_metadata.csv")

### clean up data #################
#make data long and skinny and label animal species and infection method (contact cats)
df2 <- df %>%
  pivot_longer(!c(reference_sequence, position, VOC, gene, indel,variant, reference_base, variant_base, effect), 
               names_to= "dataset_ID", values_to ="frequency") %>%
  na.omit(df2) %>%
  mutate(infection_method = ifelse(dataset_ID %in% c("Cat_6","Cat_10", "Cat_11", "Cat_12", "Cat_26"), 
                                   "contact", 
                                   ifelse(dataset_ID %in% c("Passage_1", "Passage_2", "Passage_3"), 
                                          "inoculum", "direct_inoculation"))) %>%
  mutate(species = case_when(grepl("C", dataset_ID) ~ "Cats",
                             grepl("P", dataset_ID) ~"Vero"))

#count the number of variants
df_all <- df2 %>% 
  group_by(species, dataset_ID, infection_method) %>%
  summarize(number_of_variants=n()) 

#make this df just for cats
df_cats <- df_all %>%
  filter(species == "Cats")

#add in the metadata
df.meta <- merge(df_cats, df_meta)

#for some analysis we'll need to remove contact cats as we know they are different
df.direct <- df.meta %>%
  filter(infection_method=="direct_inoculation")

### ggpairs #################
df.pairs <- df.direct %>%
  select(number_of_variants, cohort, nasal_pfu, dpi, dose_pfu)
ggpairs(df.pairs)

### distribution of # of variants ##########
hist <- ggplot(df.meta, aes(x=number_of_variants)) +
  geom_histogram(binwidth=3, position="dodge", 
                 color="grey", fill="lightgrey", size=0.3) + 
  geom_vline(aes(xintercept=median(number_of_variants)), 
             color="grey", linetype="dashed")+
  labs(x="Number of variants", y="Frequency (number of cats)") +
  theme_classic()+
  theme(axis.text=element_text(size=10),
        axis.title.y=element_text(size=12),
        axis.title.x=element_text(size=12),
        legend.position="none")

hist

#choose some colors
colors <- c("forestgreen","pink","gold1","#E64B35FF", "skyblue2","#3C5488FF")

infection_method <- ggplot(df.meta, aes(x=infection_method, y=number_of_variants))+
  geom_boxplot(outlier.shape=NA, lwd=0.3)+
  geom_jitter(aes(color=factor(cohort)), 
              position=position_jitter(0.2), alpha=0.8, size=2)+
  labs(x="Infection method", y="Number of variants", color = "cohort")+
  ylim(c(0,35))+
  scale_x_discrete(labels=c("contact" = "contact", 
                            "direct_inoculation" = "direct inoculation"))+
  scale_color_manual(labels = c("A", "B", "C"), 
                     values=c("#00A087FF","#E64B35FF", "#8491B4FF"))+
  theme_classic()+
  theme(axis.text=element_text(size=10), 
        axis.title.y=element_text(size=12),
        axis.title.x=element_text(size=12),
        legend.title = element_text(size=12),
        legend.position = "none")

infection_method

### plot dose and #variants ######################

df.direct$dose_level <- factor(df.direct$dose_level, 
                               levels=c("low", "medium", 
                                        "med-high", "high"))

dose_level <- ggplot(df.direct, aes(x=dose_level, y=number_of_variants)) +
                       geom_boxplot(outlier.shape=NA, lwd=0.3)+
                       geom_jitter(aes( color=factor(cohort),
                                        shape=dpi), 
                                   position=position_jitter(0.2), 
                                   alpha=0.8, size=2)+
                       labs(x="Dose level", y="Number of variants", color = "cohort", 
                            shape = "DPI") +
  scale_color_manual(labels = c("A", "B", "C"), 
                     values=c("#00A087FF","#E64B35FF", "#8491B4FF")) +
  ylim(c(0,35))+
  theme_classic()+
  theme(axis.text=element_text(size=10), 
        axis.title.y=element_text(size=12),
        axis.title.x=element_text(size=12),
        legend.position = "bottom",
        legend.background=element_rect(color="black"))
  
print(dose_level)

#lets look at dose as numeric values as well

df.direct$dose_pfu <- as.numeric(df.direct$dose_pfu)

dose <- ggplot(df.direct, aes(x=log10(dose_pfu), 
                              y=number_of_variants, 
                              color=factor(cohort),
                              shape=dpi)) +
  geom_point(alpha=0.8)+
  stat_smooth(method = "lm", col = "black", fill=NA)+
  scale_color_manual(values=c("#00A087FF","#8491B4FF" ,"#E64B35FF"))+
  labs(y="Number of variants", x="Dose log10(pfu)", color = "cohort")+
  theme_classic()+
  theme(text=element_text(size=18))

print(dose)

### plot viral titer and #variants ######################

#exploratory plot with all the nasal lavage data
#not including contact cats bc we know they have fewer variants due to bottleneck/something

number_variants_titer <- ggplot(df.direct, 
            aes(x=log10(nasal_pfu+0.001), y=number_of_variants, 
                shape = dpi, color=factor(cohort)), 
       alpha=0.8, size=2)+
  geom_point()+
  scale_color_manual(labels = c("A", "B", "C"), 
                     values=c("#00A087FF","#E64B35FF", "#8491B4FF")) +
  labs(x="Viral titer (log pfu/mL)", y= "Number of variants", 
       color = "Cohort")+
  theme_classic()+
  theme(axis.text=element_text(size=10), 
        axis.title.y=element_text(size=12),
        axis.title.x=element_text(size=12),
        legend.position = "top",
        legend.background=element_rect(color="black"))

number_variants_titer

## figure S4
pdf("number_of_variants_titer_direct_inoculated.pdf", width=5, height=4)
number_variants_titer
dev.off()

### final plot ############
leg <- get_legend(dose_level)
number_of_variants <- ggarrange(hist, infection_method, dose_level,
                                ncol=3, labels="auto", legend.grob = leg)

number_of_variants

pdf("number_of_variants_5_25_22.pdf", width=10, height=4)
number_of_variants
dev.off()

###count the number of variants in each gene

df %>% 
  select(position, gene, variant, 
         effect, starts_with("C")) %>%
  group_by(gene) %>%
  summarize(number_of_variants=n()) %>%
  arrange(desc(number_of_variants)) %>%
  head()

#31 variants in spike, that's 26.3%
31/118
#spike: 3822 nucleotides out of 29727 (https://www.cdc.gov/sars/lab/sequence.html)
3822/29727
#so the spike is 12.9% of genome

#19 variants in nsp3, that's 16.1%
19/118
#nsp3: 5835
5835/29727
#so nsp3 is 19.6% of genome


