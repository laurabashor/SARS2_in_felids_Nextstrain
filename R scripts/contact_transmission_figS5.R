### plotting donor vs recipient allele frequencies to visualize contact transmission of variants

#load libraries
library(tidyverse)
library(readxl)
library(ggpubr)

#set working directory
setwd("/Users/_lbashor/Dropbox/SARS-CoV-2 cat manuscript/cat ms results/R analysis cats/8_donor_recipient_AF")

# 3% cleaned, processed iVar filtered variant table
df <- read_xlsx("/Users/_lbashor/Dropbox/SARS-CoV-2 cat manuscript/cat ms results/R analysis cats/2_ivar_cats/variant_summary_0.03_iVar_processed.xlsx")

##make data frames that have all the variants that were in either cat at >3%

#cohort A
df.contact01 <- df %>% 
  select(reference_sequence,position, gene, indel,
         variant, reference_base, variant_base, effect, 
        Cat_5, Cat_6)

delete.na.contact1 <- function(df.contact01, n=0) {
  df.contact01[rowSums(is.na(df.contact01)) <= n,]
}
df.contact01 <- delete.na.contact1(df.contact01, 1)

#cohort B
df.contact02 <- df %>% 
  select(reference_sequence,position, gene, indel,
         variant, reference_base, variant_base, effect, 
         Cat_7, Cat_10)

delete.na.contact2 <- function(df.contact02, n=0) {
  df.contact02[rowSums(is.na(df.contact02)) <= n,]
}
df.contact02 <- delete.na.contact2(df.contact02, 1) 

df.contact03 <- df %>% 
  select(reference_sequence,position, gene, indel,
         variant, reference_base, variant_base, effect, 
         Cat_7, Cat_11)

delete.na.contact3 <- function(df.contact03, n=0) {
  df.contact03[rowSums(is.na(df.contact03)) <= n,]
}
df.contact03 <- delete.na.contact3(df.contact03, 1)

df.contact04 <- df %>% 
  select(reference_sequence,position, gene, indel,
         variant, reference_base, variant_base, effect, 
         Cat_7, Cat_12)

delete.na.contact4 <- function(df.contact04, n=0) {
  df.contact04[rowSums(is.na(df.contact04)) <= n,]
}
df.contact04 <- delete.na.contact4(df.contact04, 1)

#cohort C
df.contact05 <- df %>% 
  select(reference_sequence,position, gene, indel,
         variant, reference_base, variant_base, effect, 
         Cat_22, Cat_26)

delete.na.contact5 <- function(df.contact05, n=0) {
  df.contact05[rowSums(is.na(df.contact05)) <= n,]
}
df.contact05 <- delete.na.contact5(df.contact05, 1)

###make block plots

p.1 <- ggplot(df.contact01 %>% pivot_longer(!c(reference_sequence, position, 
                                       gene, indel,variant, reference_base, 
                                       variant_base, effect), 
                                   names_to= "dataset_ID", 
                                   values_to ="frequency"),
              aes(x=fct_reorder(variant, desc(position)), 
                  y=dataset_ID, fill=frequency))+
  geom_tile()+
  labs(fill="Allele frequency")+
  scale_fill_gradient(low = "#fee0d2",high = "#de2d26",
                      na.value = "white",limits=c(0,1))+
  scale_y_discrete(labels = c("Cat 5", "Cat 6")) +
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.title= element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_text(size=10))+
  coord_flip()

p.1

p.2 <- ggplot(df.contact02 %>% pivot_longer(!c(reference_sequence, position, 
                                       gene, indel,variant, reference_base, 
                                       variant_base, effect), 
                                    names_to= "dataset_ID", 
                                    values_to ="frequency"),
              aes(x=fct_reorder(variant, desc(position)), 
                  y=fct_relevel(dataset_ID, "Cat_7", "Cat_10"), fill=frequency)) +
  geom_tile()+
  labs(fill="Allele frequency")+
  scale_fill_gradient(low = "#fee0d2",high = "#de2d26",
                      na.value = "white",limits=c(0,1))+
  scale_y_discrete(labels = c("Cat 7", "Cat 10")) +
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.title= element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_text(size=10))+
  coord_flip()

print(p.2)

p.3 <- ggplot(df.contact03 %>% pivot_longer(!c(reference_sequence, position, 
                                       gene, indel,variant, reference_base, 
                                       variant_base, effect), 
                                    names_to= "dataset_ID", 
                                    values_to ="frequency"),
              aes(x=fct_reorder(variant, desc(position)), 
                  y=fct_relevel(dataset_ID, "Cat_7", "Cat_11"), fill=frequency)) +
  geom_tile()+
  labs(fill="Allele frequency")+
  scale_fill_gradient(low = "#fee0d2",high = "#de2d26",
                      na.value = "white",limits=c(0,1))+
  scale_y_discrete(labels = c("Cat 7", "Cat 11")) +
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.title= element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_text(size=10))+
  coord_flip()

print(p.3)

p.4 <- ggplot(df.contact04 %>% pivot_longer(!c(reference_sequence, position, 
                                       gene, indel,variant, reference_base, 
                                       variant_base, effect), 
                                    names_to= "dataset_ID", 
                                    values_to ="frequency"),
              aes(x=fct_reorder(variant, desc(position)), 
                  y=fct_relevel(dataset_ID, "Cat_7", "Cat_12"), fill=frequency)) +
  geom_tile()+
  labs(fill="Allele frequency")+
  scale_fill_gradient(low = "#fee0d2",high = "#de2d26",
                      na.value = "white",limits=c(0,1))+
  scale_y_discrete(labels = c("Cat 7", "Cat 12")) +
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.title= element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_text(size=10))+
  coord_flip()

print(p.4)


p.5 <- ggplot(df.contact05 %>% pivot_longer(!c(reference_sequence, position, 
                                      gene, indel,variant, reference_base, 
                                      variant_base, effect), 
                                   names_to= "dataset_ID", 
                                   values_to ="frequency"),
              aes(x=fct_reorder(variant, desc(position)), 
                  y=dataset_ID, fill=frequency)) +
  geom_tile()+
  labs(fill="Allele frequency")+
  scale_fill_gradient(low = "#fee0d2",high = "#de2d26",
                      na.value = "white",limits=c(0,1))+
  scale_y_discrete(labels = c("Cat 22", "Cat 26")) +
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.title= element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_text(size=10))+
  coord_flip()

print(p.5)

pdf("Donor_recipient_tile_plots_0.03.pdf", width=12, height=10, onefile=FALSE)
ggarrange(p.1,p.2,p.3,p.4,p.5, ncol=3, nrow=2, labels="AUTO", 
          common.legend=TRUE, legend="top")
dev.off()

pdf("Donor_recipient_tile_plots_0.03_long.pdf", width=8, height=13, onefile=FALSE)
ggarrange(p.1,p.2,p.3,p.4,p.5, ncol=2, nrow=3, labels="AUTO", 
          common.legend=TRUE, legend="top")
dev.off()


#scatterplots
p1 <- ggplot(df.contact01, aes(Cat_5,Cat_6))+
  geom_point()+
  coord_fixed(xlim=c(0,1), ylim=c(0,1))+
  xlab("AF in Cat 5")+
  ylab("AF in Cat 6")+
  theme_classic()

print(p1)

p2 <- ggplot(df.contact02, aes(Cat_7,Cat_10))+
  geom_point()+
  coord_fixed(xlim=c(0,1), ylim=c(0,1))+
  xlab("AF in Cat 7")+
  ylab("AF in Cat 10")+
  theme_classic()

print(p2)

p3 <- ggplot(df.contact03, aes(Cat_7,Cat_11))+
  geom_point()+
  coord_fixed(xlim=c(0,1), ylim=c(0,1))+
  xlab("AF in Cat 7")+
  ylab("AF in Cat 11")+
  theme_classic()

print(p3)

p4 <- ggplot(df.contact04, aes(Cat_7,Cat_12))+
  geom_point()+
  coord_fixed(xlim=c(0,1), ylim=c(0,1))+
  xlab("AF in Cat 7")+
  ylab("AF in Cat 12")+
  theme_classic()

print(p4)

p5 <- ggplot(df.contact05, aes(Cat_22,Cat_26))+
  geom_point()+
  coord_fixed(xlim=c(0,1), ylim=c(0,1))+
  xlab("AF in Cat 22")+
  ylab("AF in Cat 26")+
  theme_classic()

print(p5)

pdf("Donor_recipient_scatter_plots_0.001.pdf", width=11, height=8)
ggarrange(p1,p2,p3,p4,p5, ncol=3, nrow=2, labels="AUTO")
dev.off()

