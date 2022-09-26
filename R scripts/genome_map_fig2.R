###the purpose of this script is to look at the number and type of viral variants, and 
#to visualize them along the SARS2 genome

#load packages
library(tidyverse)
library(ggforce)
library(RColorBrewer)
library(ggpubr)
library(readxl)

#start with data frame of cleaned data, combined replicates, and iVar-filtered
df <- read_xlsx("/Users/_lbashor/Dropbox/SARS-CoV-2 cat manuscript/cat ms results/R analysis cats/2_ivar_cats/variant_summary_0.03_iVar_processed.xlsx")

#start off by counting variants! (now that we have filtered for variants not detected by iVar)

#look at variants in cats and in viral stocks separately
df.cats <- df %>% 
  select(reference_sequence,position, gene, indel,variant, 
         reference_base, variant_base, effect, starts_with("C"))

df.inoc <- df %>%
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
df.filtered <- df %>%
  filter_at(vars(starts_with("Passage")), all_vars(is.na(.)))

#so final information here is: 118 variants across all datasets, and 118 found in cats
#18 in the inoculums 
#100 in cats but not found in the inoculums (so all inoculum variants did transfer to cats)

##what's the overall allele frequency distribution?
hist(df2$frequency)

#what types of variants?
df %>%
  group_by(effect) %>%
  count() %>%
  arrange(desc(n))
# 18 structural variants, 100 SNVs (70 N, 30 S)

# missense: 63 + 6 + 1 = 70 nonsynonymous SNVs
# synonymous: 30 synonymous SNVs

# structural: 9 frameshifts, 6 disruptive inframe deletions, 
# 1 disrupt inframe insert, 1 conserv inframe del & 1 conserv inframe insert

#make data long and skinny 
df2 <- df %>%
  pivot_longer(!c(reference_sequence, position, gene,indel,
                  variant, reference_base, variant_base, VOC, effect), 
               names_to= "dataset_ID", values_to ="frequency")

#want to make a plot showing all variants and their effects
#how many variants were SNVs vs SVs?
df2 %>%
  na.omit(df2) %>%
  count(indel) #246 SNV observations, 65 SV observations

#first need to simplify the effects
df2 %>%
  na.omit(df2) %>%
  filter(indel == TRUE) %>%
  count(effect) #this gives # of observations of each effect

#indels can be 8 things:
  #conservative inframe deletion 2
  #conservative_inframe_insertion,custom 7
  # disruptive_inframe_deletion, 11
  # disruptive_inframe_deletion,custom 1
  #disruptive_inframe_insertion, 3
  # frameshift 14
  # frameshift_variant&stop_gained 16
  # frameshift_variant&stop_lost&splice_region 11

#we can condense that to three things: inframe deletion, inframe insertion, frameshift

df2 %>%
  na.omit(df2) %>%
  filter(indel==FALSE) %>%
  count(effect) #this gives # of observations of each effect

#SNVs can be 4 things:
#missense 184
#missense_variant,custom 39
#missense_variant,custom,custom 3
#synonymous 47

#we can condense that to two things: nonsynonymous_SNV and synonymous_SNV
df3 <- df2 %>%
  na.omit(df2) %>%
  mutate(effect2 = ifelse(effect %in% c("missense","missense_variant,custom", "missense_variant,custom,custom"), 
                  "nonsynonymous_SNV",
                          ifelse(effect %in% c("synonymous"),
                                 "synonymous_SNV",
                                 ifelse(effect %in% c("conservative_inframe_deletion", "conservative_inframe_deletion,conservative_inframe_deletion", "disruptive_inframe_deletion", "disruptive_inframe_deletion,custom", "disruptive_inframe_deletion,disruptive_inframe_deletion"),
                                        "inframe_deletion",
                                        ifelse(effect %in% c("conservative_inframe_insertion,custom", "disruptive_inframe_insertion"), 
                                               "inframe_insertion",
                                               ifelse(effect %in% c("frameshift","frameshift_variant,frameshift","frameshift_variant&stop_gained","frameshift_variant&stop_lost&splice_region"),
                                                      "frameshift", "noncoding_SNV"))))),
         inoculum = ifelse(dataset_ID %in% c("Passage_2", 
                                             "Cat_7"  , "Cat_8" ,"Cat_9"  , 
                                             "Cat_10", "Cat_11" , "Cat_12"), "P2", 
                           ifelse(dataset_ID == "Passage_1", "P1", "P3")))

df3 %>%
  count(effect2)

#here are the counts of variant observations

# 1 frameshift           41
# 2 inframe_deletion     14
# 3 inframe_insertion    10
# 4 nonsynonymous_SNV   200
# 5 synonymous_SNV       46

#now we can plot this

#load in nice axis labels
labels <- scan("labels.txt", character(),quote="")

#choose some colors
colors <- c("forestgreen","pink","gold1","#E64B35FF", "skyblue2","#3C5488FF")
colors2 <- c("forestgreen","pink","#E64B35FF", "skyblue2","#3C5488FF")

##plot variants: split panels into the P2 or P3 inoculated cats #################

P3 <- df3 %>%
  filter(dataset_ID %in% c("Passage_3",
                           "Cat_1", "Cat_5" ,"Cat_6",
                           "Cat_13" , "Cat_14", "Cat_15" ,
                           "Cat_16" , "Cat_17", "Cat_18",
                           "Cat_19" , "Cat_20" , "Cat_21" ,
                           "Cat_22", "Cat_23",  "Cat_24"  ,
                           "Cat_25" , "Cat_26"))

P2 <- df3 %>%
  filter(dataset_ID %in% c("Passage_2",
                           "Cat_7"  , "Cat_8" ,"Cat_9"  ,
                           "Cat_10", "Cat_11" , "Cat_12"))

df4 <- df3 %>% 
  filter(inoculum != "P1")

P2_plot <- ggplot(df4) +
  geom_point(data=(df4 %>%
                     filter(!is.na(frequency)) %>%
                     filter(dataset_ID %in% c("Passage_2", "Cat_7", 
                                              "Cat_8", "Cat_9", "Cat_10",  
                                              "Cat_11" ,  "Cat_12"))), 
             aes(y=fct_relevel(dataset_ID, 
                               c("Passage_2",  
                                 "Cat_7"  , "Cat_8" ,  "Cat_9"  , "Cat_10",  
                                 "Cat_11" ,  "Cat_12" )),
                 x=position, 
                 col=factor(effect2, 
                            levels=c("frameshift", "inframe_deletion", 
                                     "inframe_insertion", "synonymous_SNV", 
                                     "nonsynonymous_SNV")), 
                 shape=indel), 
             size=3, alpha=0.8, stroke=FALSE)  +
  geom_vline(aes(xintercept=21563))+
  geom_vline(aes(xintercept=25384))+
  scale_y_discrete(labels=c("P2",  
                            "Cat 7"  , "Cat 8" ,  "Cat 9"  , "Cat 10*",  
                            "Cat 11*" ,  "Cat 12*")) +
  scale_color_manual(name="Variant effect",
                     labels=c("Frameshift",
                              "Inframe deletion",
                              "Inframe insertion",
                              "Synonymous SNV",
                              "Nonsynonymous SNV"),
                     values = colors) +
  scale_shape(name="Variant type", labels=c("SNV", "SV"))+
  theme_bw()+
  theme_classic()+
  theme(legend.background=element_rect(color="black"),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 14), 
        legend.text  = element_text(size = 12),
        legend.key.size = unit(0.8, "lines"))


P3_plot <- ggplot(P3) +
    geom_point(data=(P3 %>%
                       filter(!is.na(frequency))),
               aes(y=fct_relevel(dataset_ID, 
                                 c("Passage_3",
                                   "Cat_1", "Cat_5" ,"Cat_6",
                                   "Cat_13" , "Cat_14", "Cat_15" ,
                                   "Cat_16" , "Cat_17", "Cat_18",
                                   "Cat_19" , "Cat_20" , "Cat_21" ,
                                   "Cat_22", "Cat_23",  "Cat_24"  ,
                                   "Cat_25" , "Cat_26")),
                   x=position, 
                   col=factor(effect2, 
                              levels=c("frameshift", "inframe_deletion", 
                                       "inframe_insertion", "synonymous_SNV", 
                                       "nonsynonymous_SNV")), 
                   shape=indel), 
               size=3, alpha=0.8, stroke=FALSE)  +
    geom_vline(aes(xintercept=21563))+
    geom_vline(aes(xintercept=25384))+
    scale_y_discrete(labels=c("P3",
                              "Cat 1", "Cat 5" ,"Cat 6*",
                              "Cat 13" , "Cat 14", "Cat 15" ,
                              "Cat 16" , "Cat 17", "Cat 18",
                              "Cat 19" , "Cat 20" , "Cat 21" ,
                              "Cat 22", "Cat 23",  "Cat 24"  ,
                              "Cat 25" , "Cat 26*"))+
    scale_color_manual(name="Variant effect", 
                       labels=c("Frameshift", 
                                "Inframe deletion",
                                "Inframe insertion",
                                "Synonymous SNV",
                                "Nonsynonymous SNV"),
                       values=colors2)+
    scale_shape(name="Variant type", labels=c("SNV", "SV"))+
    theme_bw()+
    labs(x="Position in genome (nt)") + 
    theme_classic()+
    theme(axis.title.x = element_text(size=16),
          legend.position="none",
          axis.title.y=element_blank(),
          axis.text=element_text(size=14))

P2_P3 <- ggarrange(P2_plot, P3_plot, ncol=1, heights = c(1,2.57), 
                   legend = "right", common.legend = TRUE, 
                   labels="auto", hjust=-0.2) +
  theme(plot.margin = margin(1,0.5,0.5,2, "cm")) 

ragg::agg_tiff("fig2.tiff", width = 10, height = 8, 
               units = "in", res = 500, compression = "lzw")
P2_P3
dev.off()
