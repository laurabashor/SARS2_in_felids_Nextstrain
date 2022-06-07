### just the stats for looking at numbers of variants (no plots)
library(tidyverse)
library(broom)

### load and clean data ###################
#start with data frame of cleaned data, combined replicates, and iVar-filtered
df <- read_xlsx("/Users/_lbashor/Dropbox/SARS-CoV-2 cat manuscript/cat ms results/R analysis cats/2_ivar_cats/variant_summary_0.03_iVar_processed.xlsx")

#load in additional metadata
df_meta <- read_csv("/Users/_lbashor/Dropbox/SARS-CoV-2 cat manuscript/cat ms results/metadata/cat_metadata.csv")

#make data long and skinny and label animal species and infection method (contact cats)
df2 <- df %>%
  pivot_longer(!c(reference_sequence, position, gene, indel,variant, reference_base, variant_base, effect), 
               names_to= "dataset_ID", values_to ="frequency") %>%
  na.omit(df2) %>%
  mutate(infection_method = ifelse(dataset_ID %in% c("Cat_6","Cat_10", "Cat_11", "Cat_12", "Cat_26"), 
                                   "contact", 
                                   ifelse(dataset_ID %in% c("Passage_1", "Passage_2", "Passage_3"), 
                                          "inoculum", "direct_inoculation"))) %>%
  mutate(species = case_when(grepl("C", dataset_ID) ~ "Cats",
                             grepl("P", dataset_ID) ~"Vero"))

### counting variant richness & cats vs inoculum #############################
#what is the mean/median number of variants overall?
df_all <- df2 %>% 
  group_by(species, dataset_ID, infection_method) %>%
  summarize(number_of_variants=n()) 

mean(df_all$number_of_variants) #11.96 overall
median(df_all$number_of_variants) #10 overall

#what about the inoculum samples?
#P1 had 7 variants, P2 had 15, and P3 had 8
df_all %>%
  filter(species == "Vero") %>%
  group_by(species) %>%
  summarize(mean = mean(number_of_variants), # 10
            median = median(number_of_variants)) # 8 

#make this df just for cats
df_cats <- df_all %>%
  filter(species == "Cats")

df_cats %>%
  group_by(species) %>%
  summarize(mean = mean(number_of_variants), # 12.2
            median = median(number_of_variants)) # 10

range(df_cats$number_of_variants) # 4 to 33

se <- function(x) {sqrt(var(x)/length(x))}
se(df_cats$number_of_variants)
sd(df_cats$number_of_variants)

#so mean is 12 +/- 1.4 se or 6.8 sd

#is there a difference between cats and inoculum?
var.test(df_all$number_of_variants ~ df_all$species) #not sig, variance is equal
t.test(df_all$number_of_variants~df_all$species, var.equal=TRUE) 
#no significant difference, p=0.5914

### number of variants and infection method #######################

##now compare the mean number of variants by infection method
#mean in contact cats is 6.4, mean in direct is 13.8
df_cats %>% 
  group_by(infection_method) %>% 
  summarize(mean = mean(number_of_variants),
            median = median(number_of_variants),
            se = se(number_of_variants))

var.test(df_cats$number_of_variants~df_cats$infection_method) 
#there is a sig difference between the variances
#can go ahead and use t.test, which does Welch's t-test assuming unequal variances
t.test(df_cats$number_of_variants~df_cats$infection_method) #significant difference p=0.0002725
#unpaired because they are all from different individual cats


### make df.meta and df.direct ############

#add in the metadata
df.meta <- merge(df_cats, df_meta)

#for some analysis we need to remove contact cats as we know they are different
df.direct <- df.meta %>%
  filter(infection_method=="direct_inoculation")

### number of variants and cohort ############

##look at cohort and number of variants for everyone, contact and direct
df.meta$cohort <- as.factor(df.meta$cohort)
aov.cohort <- aov(number_of_variants ~ cohort, df.meta)
summary(aov.cohort)
#p=0.0807, indicating there isn't a significant difference in the #variants among cohorts 
#R2=0.145
summary(lm(number_of_variants ~ cohort, df.meta))

##look at cohort and number of variants for just directly inoculated cats
df.direct$cohort <- as.factor(df.direct$cohort)
aov.cohort2 <- aov(number_of_variants ~ cohort, df.direct)
summary(aov.cohort2)

summary(lm(number_of_variants ~ cohort, df.direct))
#p=0.154, indicating there isn't a significant difference in the #variants among cohorts 
#R2=0.12

### number of variants and cat age and sex ############
##cat sex and number of variants
var.test(number_of_variants ~ sex, df.direct) # equal variances
t.test(df.direct$number_of_variants ~ df.direct$sex, var.equal = T)
#p=0.4777, not significant

##age and number of variants
var.test(number_of_variants ~ age, df.direct)
#no significant diff between variances, equal variances t-test
t.test(df.direct$number_of_variants ~ df.direct$age, var.equal=TRUE)
#p=0.2253, no difference in number of variants by cat age (adult or juvenile)

### number of variants and dpi ##########################

var.test(number_of_variants ~ dpi, df.direct) #equal variances
t.test(number_of_variants ~ dpi, df.direct, var.equal = T) #p = 0.04998
#mean of 1dpi: 5, mean of 3dpi: 14.9375

### number of variants and dose ############

#now look at dose and number of variants-- df.direct

df.direct$dose_level <- factor(df.direct$dose_level, 
                               levels=c("low", "medium", "medium-high", "high"))

summary(lm(number_of_variants ~ dose_level, df.direct)) %>%
  tidy()
#p = 0.0205 for low-high, close to sig for low-medium 0.0554
#overall p value 0.02168, R2=0.3775

summary(aov(number_of_variants ~ dose_level, df.direct)) #same p = 0.0217

## dose as numeric (log-transformed pfu)

df.direct$dose_pfu <- as.numeric(df.direct$dose_pfu)

summary(lm(number_of_variants ~ log10(dose_pfu), df.direct))
# significant p = 0.03265, R2 = 0.2181

#if we add in dpi? model fits even better
summary(lm(number_of_variants ~ log10(dose_pfu) + dpi, df.direct))
#even more significant, p = 0.00697, R2=0.4155

### number of variants and viral titer ###################

# look at log transformation of nasal titer data
pfu.v.all <- df.meta$nasal_pfu
hist(pfu.v.all)
hist(log10(pfu.v.all))

pfu.v <- df.direct$nasal_pfu
hist(log10(pfu.v))

#linear model with log transformation
summary(lm(number_of_variants ~ log10(nasal_pfu+0.001), df.direct))
#p=0.3025 not significant

#if we add in dpi?
summary(lm(number_of_variants ~ log10(nasal_pfu+0.001) + dpi, df.direct))
#wow! suddenly p=0.02087, and p=0.0484 for titer, p=0.0108 for dpi

#SO we know that number of variants is affected by viral titer and dpi together

#and it's not just that viral titer is determined by dpi, because that is not significant
summary(lm(log10(nasal_pfu+0.001) ~ dpi, df.direct))

var.test(log10(nasal_pfu+0.001) ~ dpi, df.direct)
t.test(log10(nasal_pfu+0.001) ~ dpi, df.direct, var.equal = T)

### correlation between dose and nasal pfu ###########
summary(lm(log10(nasal_pfu+0.001) ~ log10(dose_pfu), df.direct)) 
summary(lm(log10(nasal_pfu+0.001) ~ log10(dose_pfu) + dpi, df.direct)) 
#p=0.0057, R2=0.3501 without dpi
#p=0.007257 overall and R2 = 0.4123 with dpi

#just dose and titer, not significant
summary(lm(number_of_variants ~ log10(nasal_pfu+0.001) + dose_pfu, df.direct))

### put all of these together: #variants with dose, titer, dpi #############
summary(lm(number_of_variants ~ log10(nasal_pfu+0.001) + dpi + dose_pfu, df.direct))

#it seems like dpi is driving some of this?

### infection success ##########
df3 <- df3 %>%
  mutate(infection_success = ifelse(nasal_pfu>0, "yes", "no"))

var.test(df3$number_of_variants ~ df3$infection_success) 
#there is no sig difference between the variances

t.test(df3$number_of_variants ~ df3$infection_success) #p=0.1196