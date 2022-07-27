#the purpose of this script is to visualize the selection analysis output by SNPGenie

#SNPGenie generated a population_summary file and product_results file
#it did this for each sequencing replicate for each cat
#these files were generated individually for each sample so I combined and cleaned up the messy data with the selection_cleanup.R script

#load packages
library(tidyverse)
library(ggpubr)
library(openxlsx)
library(broom)

#load in cleaned up files
setwd("/Users/_lbashor/Dropbox/SARS-CoV-2 cat manuscript/cat ms results/R analysis cats/7_selection_cats")
df_pop <- read_csv("population_summary_cleaned.csv")
df_prod <- read_csv("gene_product_cleaned.csv")

#for now, let's not look at cell culture passages
df_pop <- df_pop %>%
  filter(!grepl("Passage", dataset_ID))

df_prod <- df_prod %>%
  filter(!grepl("Passage", dataset_ID))

#first let's look at these data at the population level#########################

#get the mean of the two technical replicates for each cat
pop <- df_pop %>%
  mutate(piN = as.numeric(piN), 
         piS = as.numeric(piS)) %>%
  select(dataset_ID, replicate, piN, piS) %>%
  group_by(dataset_ID) %>%
  summarize(mean_piN = mean(piN), 
            mean_piS = mean(piS))

#calculate piN/piS (name the column piN_piS) and add logical column
#also add info on inoculation (contact or direct)
pop <- pop %>%
  mutate(piN_piS = (mean_piN/mean_piS)) %>% 
  mutate(piN_greater = ifelse((mean_piN > mean_piS), TRUE, FALSE)) %>%
  mutate(infection_method = ifelse(dataset_ID %in% 
                                     c("Cat_6","Cat_10", "Cat_11", 
                                       "Cat_12", "Cat_26"), 
                                   "contact", 
                                   "direct_inoculation"))
  
##now do t-tests to compare nonsynonymous nucleotide diversity to synonymous

#first, look at the difference between piN and piS at a population level
t.test(x=pop$mean_piN, y=pop$mean_piS, paired=TRUE) #p=0.8601, NS

#how many cats had piN>piS overall? 
pop %>%
  count(piN_greater) #15 cats, the other 8 had piS>piN

pop %>%
  filter(!piN_greater) #Cats 1, 13, 15, 16, 17, 20, 25, 26

#all of these piS > piN were directly inoculated cats except for Cat 26

#visualize it
plot1 <- ggpaired(pop, cond1 = "mean_piN", cond2 = "mean_piS",
                  color = "condition", line.color = "gray", line.size = 0.4,
                  palette="npg", ggtheme = theme_classic(), xlab=FALSE, ylab=FALSE, legend="none", ylim=c(0, 0.00032))+
  stat_compare_means(method="t.test", paired=TRUE, size=4, 
                     vjust=-2, label="p.format")

plot1

#there seem to be a couple cats with really dramatic differences
#Cat 16 had the highest piS and most dramatic difference,
#the other two big differences of piS > piN are Cat 16 and Cat 1

#what about contact vs. directly inoculated cats?
#not seeing a pattern here
pop_contact <- pop %>%
  filter(infection_method == "contact") 

t.test(x=pop_contact$mean_piN, y=pop_contact$mean_piS, paired=TRUE) 
#p=0.3198, NS

pop_direct <- pop %>%
  filter(infection_method != "contact") 

t.test(x=pop_direct$mean_piN, y=pop_direct$mean_piS, paired=TRUE) 
#p=0.9265, NS

#visualize it
plot1b <- ggpaired(pop, cond1 = "mean_piN", cond2 = "mean_piS",
                   color = "condition", line.color = "gray", line.size = 0.4,
                   palette="npg", ggtheme = theme_classic(), xlab=FALSE, ylab=FALSE, legend="none", ylim=c(0, 0.00032))+
  stat_compare_means(method="t.test", paired=TRUE, size=4, 
                     vjust=-2, label="p.format") +
  facet_wrap(~infection_method)

plot1b

#what about different cohorts?

df_meta <- read_csv("/Users/_lbashor/Dropbox/SARS-CoV-2 cat manuscript/cat ms results/metadata/cat_metadata.csv")

pop_meta <- merge(pop, df_meta)

cohortA <- pop_meta %>%
  filter(cohort == "A")
cohortB <- pop_meta %>%
  filter(cohort == "B")
cohortC <- pop_meta %>%
  filter(cohort == "C")

t.test(x=cohortA$mean_piN, y=cohortA$mean_piS, paired=TRUE) #P=0.7195
t.test(x=cohortB$mean_piN, y=cohortB$mean_piS, paired=TRUE) #P=0.01843
t.test(x=cohortC$mean_piN, y=cohortC$mean_piS, paired=TRUE) #P=0.5308

plot1c <- ggpaired(pop_meta, cond1 = "mean_piN", cond2 = "mean_piS",
                   color = "condition", line.color = "gray", line.size = 0.4,
                   # label = "dataset_ID",
                   palette="npg", ggtheme = theme_classic(), xlab=FALSE, ylab=FALSE, legend="none", ylim=c(0, 0.00032))+
  stat_compare_means(method="t.test", paired=TRUE, size=4, 
                     vjust=-2, label="p.format") +
  facet_wrap(~cohort)

plot1c

#interestingly, at the summary level piN > piS significantly for Cohort B
#they were infected with P2 virus, three direct and three contact cats

pop_meta %>% filter(cohort == "B")

#now let's look at these data by gene product ##################################

#get the mean of the two technical replicates for each cat
prod <- df_prod %>%
  mutate(piN = as.numeric(piN), 
         piS = as.numeric(piS)) %>%
  select(dataset_ID, product, replicate, piN, piS) %>%
  group_by(dataset_ID, product) %>%
  summarize(mean_piN = mean(piN), 
            mean_piS = mean(piS))

#calculate piN/piS (name the column piN_piS) and add logical columns as well
prod <- prod %>%
  mutate(piN_piS = (mean_piN/mean_piS)) %>% 
  mutate(piN_greater = ifelse((mean_piN > mean_piS), TRUE, FALSE))

#ttest time

#group the data by gene, and then do a t-test comparing piN to piS
prod_summary <- prod %>%
  group_by(product) %>%
  nest() %>%
  mutate(t_test = map(data, ~t.test(x = .x$mean_piN, 
                                    y = .x$mean_piS, 
                                    paired=TRUE)))

#ungroup the data into a summary table
unnested_summary <- prod_summary %>%
  mutate(summary = map(t_test, ~ tidy(.x, drop=TRUE))) %>%
  select(-data, -t_test) %>%
  unnest(summary) %>%
  ungroup()

options(digits=7)
print(unnested_summary)

unnested_summary %>% 
  filter(p.value < 0.05)
#looks super significant in spike
#also significant in orf1ab

#18 out of 23 cats had piN>piS in Spike
#which cats had piS>piN?
prod %>%
  filter(product == "S") %>%
  filter(!piN_greater)

prod %>%
  filter(product == "orf1ab") %>%
  filter(!piN_greater) #14 cats had piS>piN, 9 had piN>piS

#Cats 1, 11, 12, 17, and 26 had piS > piN in Spike
#that's three out of five contact cats

#interesting as piN>piS overall for Cats 11 & 12 and all cats in Cohort B

prod <- prod %>%
  mutate(infection_method = ifelse(dataset_ID %in% 
                                     c("Cat_6","Cat_10", "Cat_11", 
                                       "Cat_12", "Cat_26"), 
                                   "contact", 
                                   "direct_inoculation"))

prod_contact <- prod %>%
  filter(infection_method == "contact")

prod_direct <- prod %>%
  filter(infection_method != "contact")

summary_direct <- prod_direct %>%
  group_by(product) %>%
  nest() %>%
  mutate(t_test = map(data, ~t.test(x = .x$mean_piN, 
                                    y = .x$mean_piS, 
                                    paired=TRUE)))

unnested_summary_direct <- summary_direct %>%
  mutate(summary = map(t_test, ~ tidy(.x, drop=TRUE))) %>%
  select(-data, -t_test) %>%
  unnest(summary) %>%
  ungroup()

unnested_summary_direct %>% 
  filter(p.value < 0.05)

prod_direct %>%
  filter(product == "orf1ab") %>%
  filter(!piN_greater) #12 cats with piS greater, 6 with piN greater 

prod <- prod %>%
  mutate('mean expression(piN)' = mean_piN)

#now we can plot some of this

plot2 <- ggpaired(prod, cond1 = "mean_piN", cond2 = "mean_piS",
                  color = "condition", 
                  line.color = "gray", 
                  line.size = 0.4, palette="npg",
                  xlab=FALSE, ylab=FALSE, 
                  font.tickslab = c(12, "plain", "black"),
                  font.legend = c(12, "plain", "black"),
                  legend.title = "",
                  ylim=c(0,0.0025), ggtheme = theme_classic())+
    stat_compare_means(method="t.test", paired=TRUE, size=4, 
                     vjust=-2, label="p.format") +
  scale_color_manual(labels = c('mean_piN' = expression(pi*N),
                              'mean_piS' = expression(pi*S)),
                       values = c("#E64B35FF", "skyblue2")) +
  scale_x_discrete(labels = c('mean_piN' = expression(pi*N),
                              'mean_piS' = expression(pi*S)))

plot2 <- facet(plot2, facet.by = "product", 
               ncol = 5, panel.labs.font = list(face = "bold", size = 12))
        
plot2

#if we look at just directly inoculated cats, patterns are the same
#if we do this for just contact cats, nothing is significant

plot2c <- ggpaired(prod_contact, cond1 = "mean_piN", cond2 = "mean_piS",
                  color = "condition", label = "dataset_ID",
                  line.color = "gray", 
                  line.size = 0.4, palette="npg",
                  xlab=FALSE, ylab=FALSE, legend="none", 
                  ylim=c(0,0.0012), ggtheme = theme_classic())+
  facet_wrap(~product, ncol = 5) +
  stat_compare_means(method="t.test", paired=TRUE, size=2, 
                     vjust=-2, label="p.format")

plot2c


#final plots/sheets for these results#################################
#(1) piN and piS for all datasets at a population level,
#and by infection method and by cohort 

#(2) piN vs piS by gene across all datasets

#(1)
plot1
plot1b
plot1c

#(2) 
plot2

pdf("selection_population_level.pdf", onefile=F)
plot1
dev.off()

pdf("selection_infection_method.pdf", onefile=F)
plot1b
dev.off()

pdf("selection_cohort.pdf", onefile=F)
plot1c
dev.off()

pdf("selection_gene_level.pdf", onefile=F)
plot2
dev.off()

#export dataframes as an excel sheet for supplemental tables
wb <- createWorkbook("supplemental_selection_tables.xlsx")
addWorksheet(wb, "population_level")
addWorksheet(wb, "gene_product_level")
writeData(wb, "population_level", df_pop, borders="all")
writeData(wb, "gene_product_level", df_prod,borders="all")
saveWorkbook(wb, "supplemental_selection_tables.xlsx", overwrite = TRUE)


