#can we show that higher titer relates to transmission?

#were viral loads higher in the contact exposed cats than their proposed donor?
    #A: not exactly: in two cats yes, in three others no
    #A: not significant when you average contact pfus you get 1165, donor pfus average 1900

#or than the mean titer of viral inoculum exposed cats?  
    #A: yes, but it wasn't significant by t-test. 
    # average contact titer: 1165, direct inoculation: 305

#can we compare the titers of three donor cats to infer if 
#there is a minimally infectious dose or just the most infectious donor in the room

#load libraries
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(readxl)

#start with data frame of cleaned data, combined replicates, iVar filtered
df <- read_xlsx("/Users/_lbashor/Dropbox/SARS-CoV-2 cat manuscript/cat ms results/R analysis cats/2_ivar_cats/variant_summary_0.03_iVar_processed.xlsx")
#and metadata table
df_meta <- read_csv("/Users/_lbashor/Dropbox/SARS-CoV-2 cat manuscript/cat ms results/metadata/cat_metadata.csv")


#create a contact column for donor, nondonor, or recipient cats
df_meta <- df_meta %>%
  mutate(contact = ifelse(dataset_ID %in% 
                            c("Cat_5","Cat_7", "Cat_22"), 
                          "donor", ifelse(dataset_ID %in% 
                                            c("Cat_6","Cat_10","Cat_11", 
                                              "Cat_12", "Cat_26"),
                                          "recipient", "nondonor")))

#create an infection method column for contact or direct inoc cats
df_meta <- df_meta %>%
  mutate(inf_method = ifelse(dataset_ID %in% 
                               c("Cat_6","Cat_10","Cat_11","Cat_12", "Cat_26"),
                                          "contact", "direct_inoculation"))

####################################################################

#compare donor vs. contact/recipient -- it seems that there is not a clear difference 
#between donor pfu and recipient pfu
df3 <- df_meta %>%
  filter(contact != "nondonor")

df4 <- df3 %>%
  filter(cohort != "A") #because their dpi was weird

var.test(df3$nasal_pfu~df3$contact) #no sig dif between variances
t.test(df3$nasal_pfu~df3$contact, var.equal = T) #p=0.4818

var.test(df4$nasal_pfu~df4$contact) #sig dif between variances
t.test(df4$nasal_pfu~df4$contact) #p=0.374

####################################################################

#compare contact vs. directly inoculated for all cats
var.test(df_meta$nasal_pfu ~ df_meta$inf_method) #sig dif between variances
t.test(df_meta$nasal_pfu ~ df_meta$inf_method) #p=0.3728

#not sure if this is best to do, might need to filter out cohort A because samples were from 1dpi
#then we have two donors, and a lot of nondonors, and four contact cats
meta_3dpi <- df_meta %>%
  filter(cohort != "A")

#compare contact and directly inoculated 3dpi cats
var.test(meta_3dpi$nasal_pfu ~ meta_3dpi$inf_method) #sig dif between variances
t.test(meta_3dpi$nasal_pfu ~ meta_3dpi$inf_method) #p=0.4113
#no difference in nasal pfu between contact and directly inoculated cats
#mean in contact is 1165, mean in direct inoculation is 305


####################################################################

#compare donors and nondonors within cohorts? 
#(just keep P2 cats (Cat 7, 8, 9) and P3 high dose group (Cat 22-25))
df2 <- meta_3dpi %>%
  filter(contact != "recipient") %>%
  filter(dose_pfu != 10500) %>%
  filter(dose_pfu != 465000)
  
var.test(df2$nasal_pfu~df2$contact) #sig dif between variances
t.test(df2$nasal_pfu~df2$contact) #p=0.2048
#donor mean is 1900, nondonor mean is 116
#no significant difference, maybe because we only can look at these two donors?
#even if we include cohort A it's still the same result


dfB <- df2 %>%
  filter(cohort=="B")

dfC <- df_meta %>%
  filter(cohort=="C")

dfA <- df_meta %>%
  filter(cohort == "A",
         dataset_ID == "Cat_5")

pA <- ggplot(dfA) +
  geom_bar(aes(x=reorder(dataset_ID, nasal_pfu), y=log10(nasal_pfu+0.001)), 
           stat = "identity",
           fill = "forestgreen") +
  labs(y = "Viral titer (log10pfu)", x="") +
  theme_bw()

pA

pB <- ggplot(dfB) +
  geom_bar(aes(x=reorder(dataset_ID, nasal_pfu), y=log10(nasal_pfu+0.001)), 
           stat = "identity",
           fill = "skyblue") +
  labs(y = "Viral titer (log10pfu)", x="") +
  theme_bw()

pC <- ggplot(dfC) +
  geom_bar(aes(x=reorder(dataset_ID, nasal_pfu), y=log10(nasal_pfu+0.001)), 
           stat = "identity",
           fill = "darkorange") +
  labs(y = "", x="") +
  theme_bw()

ggarrange(pB, pC)


dfABC <- df_meta %>%
  filter(dataset_ID %in% c("Cat_5", 
                           "Cat_7", "Cat_8", "Cat_9",
                           "Cat_22", "Cat_23", "Cat_24", "Cat_25"))

ggplot(dfABC) +
  geom_bar(aes(x=reorder(dataset_ID, nasal_pfu), 
               y=log10(nasal_pfu), 
               fill=cohort), 
           stat = "identity") +
  scale_fill_manual(values = c("forestgreen","skyblue", "darkorange")) +
  labs(y = "", x="") +
  theme_bw() + 
  facet_wrap(~cohort)

#in two cats recipient titer was greater, in three cats donor was greater

#Cat 6 > Cat 5 but close
#Cat 10 > Cat 7

#Cats 11 & 12 < Cat 7
#Cat 26 < Cat 22
