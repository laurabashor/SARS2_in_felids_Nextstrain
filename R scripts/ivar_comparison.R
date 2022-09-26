## the purpose of this script is to compare the variant table output by iVar to the one output by the viral_variant_caller pipeline 

#note: need to add in Cat 5 after processing replicates for the rest

#load libraries
library(tidyverse)
library(openxlsx)
library(readxl)
library(data.table)

# load in and clean up data ####################################
setwd("/Users/_lbashor/Dropbox/SARS-CoV-2 cat manuscript/cat ms results/R analysis cats/2_ivar_cats")
files  <- list.files(path = "./ivar_variant_tables", pattern = "*.tsv", full.names=T) #list all variant tables

tables.list <- lapply(files, read_tsv, skip=1) #read all the variant tables into a list 
#skip the first row as iVar output the sample name into the column headers

names <- c("dataset_ID", "position", "reference_base", "variant_base",
           "gene","reference_CODON","reference_AA","variant_CODON",
           "variant_AA", "reference_depth_rep1", "REF_RV_rep1" ,
          "REF_QUAL_rep1" , "variant_depth_rep1", "ALT_RV_rep1", 
          "ALT_QUAL_rep1", "variant_frequency_rep1", "TOTAL_DP_rep1", 
          "PVAL_rep1", "PASS_rep1", "reference_depth_rep2","REF_RV_rep2" , 
          "REF_QUAL_rep2", "variant_depth_rep2",  "ALT_RV_rep2",
           "ALT_QUAL_rep2", "variant_frequency_rep2", "TOTAL_DP_rep2",
          "PVAL_rep2", "PASS_rep2")

tables.list<- lapply(tables.list, setNames, names) #add column names

ivar_df <- rbindlist(tables.list) #now stick all the tables together

ivar_df$dataset_ID <- gsub(".filtered.tsv MN985325", "", ivar_df$dataset_ID) #clean up datasetID name

#let's explore iVar's PASS/FAIL categories
ivar_df_fail <- ivar_df %>%
  filter(PASS_rep1 != T) 

ivar_df_fail %>%
  count(PASS_rep2) #one TRUE in rep2 where rep1 was FALSE

ivar_df_fail2 <- ivar_df %>%
  filter(PASS_rep2 != T)

ivar_df_fail2 %>%
  count(PASS_rep1) #five TRUEs in rep1 where rep2 was FALSE

#I'm not totally sure why there's some variants where one replicate passed but the other didn't

# initial processing to make iVar variant table ##################################################

#moving forward we just want variants that were present and "passed" in both reps
#and then calculate allele frequency as the mean of the two tech reps
#and also filter out noncoding regions (<266, >29674)
ivar_df_pass <- ivar_df %>%
  filter(PASS_rep1 == T) %>%
  filter(PASS_rep2 == T) %>% 
  mutate(frequency = (variant_frequency_rep1+variant_frequency_rep2)/2) %>%
  filter(position>266, 
         position<29674)

#add in Cat_5
cat5 <- read_tsv("Cat_5.tsv") %>%
  filter(PASS == T) %>%
  rename(dataset_ID = 1, position = POS,
         reference_base = REF, variant_base = ALT, 
         frequency = ALT_FREQ, gene = GFF_FEATURE, reference_CODON = REF_CODON,
         reference_AA = REF_AA, variant_CODON = ALT_CODON,
         variant_AA = ALT_AA)

cat5$dataset_ID <- gsub("_ivar.tsv MN985325", "", cat5$dataset_ID) #clean up datasetID name

ivar_df_cat5 <- full_join(ivar_df_pass, cat5, by = c("dataset_ID", "position", 
                                                     "reference_base", 
                                                     "variant_base",
                                                     "frequency", "gene", 
                                                     "reference_CODON",
                                                     "reference_AA",
                                                     "variant_CODON",
                                                     "variant_AA"))

#make a wide variant table ala our other variant table from viral_variant_caller pipeline
ivar <- ivar_df_cat5 %>%
  pivot_wider(id_cols = c(position, reference_base, variant_base, 
                          gene, reference_AA, variant_AA),
              names_from = dataset_ID, 
              values_from = frequency) %>%
  arrange(position) 

#clean up the environment before moving ahead
rm(files, names, tables.list, ivar_df, ivar_df_pass, ivar_df_fail, ivar_df_fail2, cat5, ivar_df_cat5)

# now compare to viral variant caller ######################################

#load in variant table (viral variant caller = vvc)
vvc <- read_xlsx("/Users/_lbashor/Dropbox/SARS-CoV-2 cat manuscript/cat ms results/R analysis cats/1_initial_processing_0.001_cats/variant_summary_0.03_processed.xlsx")

#filter it to just cats, and make sure there's no variants hanging around that were in only the inoculums
vvc <- vvc %>%
  select(-c("Passage_1", "Passage_2", "Passage_3")) %>%
  pivot_longer(!c(reference_sequence, position, gene, indel,variant, 
                  reference_base, variant_base, effect, VOC),
               names_to= "dataset_ID",
               values_to ="frequency") %>%
  filter(!is.na(frequency)) %>%
  pivot_wider(names_from = "dataset_ID", values_from = "frequency")

#let's compare between the iVar table and the vvc table 

#already we can see that we have 129 variants called by vvc, 128 by iVar 
#that seems really close! it could be just one variant different, or it could be more

#one known difference here is that with iVar I did a 3% cutoff separately for each replicate
#whereas for vvc I took the mean of replicates then did the 3% cutoff
#so for the iVar table it was potentially more stringent, as both replicates needed to be above 3%

#now let's pull out where there are different positions between the two tables
#have to do it separately for indels vs snvs

ivar_snvs <- ivar %>%
  filter(!nchar(variant_base) > 1)
vvc_snvs <- vvc %>%
  filter(!indel)

anti_ivar <- anti_join(ivar_snvs, vvc_snvs, by = c("position", "reference_base")) #snv positions in the ivar table not in the vvc table
anti_vvc <- anti_join(vvc_snvs, ivar_snvs, by = c("position", "reference_base")) #snv positions in the vvc table not in the ivar table

#so it's looking like more than two are different
#but what is up with the NINE new variants in Cat 19 according to iVar (SEVEN in orf1ab)

#importantly, there are 16 variants detected by vvc but not iVar-- should we throw them out?

#for some of these, iVar may have dropped them because one rep was below 3%
#in these cases I will KEEP the variant

# -E120K we detected one rep at 2.4%
# -L705L we detected one rep at 2.7%
# -T913I we detected one rep at 2.9%
# -R682W we detected one rep at 2.9%
# -A134V we detected one rep at 2.9%

# -T556P all right around the edge, many reps were 2-3%, 
# BUT it's weird that we found it in so many cats and not at all with iVar
# because it was not at all in iVar, I think we will throw it out

#ultimately we may want THROW OUT: 
# 15106/T556P, 1029/F75S, 2549/D582Y, 4668/S650F, 7988/A1757S
# 10605/P184H, #13035/A4V, 13394/K124E, 
# 15036/K532N, #15168/L576L, #20054/E145G

#is there anything egregious?
#S650F reaches >50% in Cat 9, it's in vvc but not iVar
#this position had very low sequencing coverage 
# and it is in an artic primer binding region

#K124E is also in a primer binding region, which may explain something

ivar_indels <- ivar %>%
  filter(nchar(variant_base) > 1)
vvc_indels <- vvc %>%
  filter(indel)

anti_ivar_indels <- anti_join(ivar_indels, vvc_indels, by = c("position")) #indel positions in the ivar table not in the vvc table
anti_vvc_indels <- anti_join(vvc_indels, ivar_indels, by = c("position")) #indel positions in the vvc table not in the ivar table

#for the one indel not detected with iVar, G89fs, it was found in one rep at 1%
#so it's another artifact of when we apply the 3% cutoff and we can KEEP it

#let's save all this info as an excel spreadsheet

wb <- createWorkbook("variant_comparison_ivar_vvc.xlsx")
addWorksheet(wb, "vvc_variants")
addWorksheet(wb, "ivar_variants")
addWorksheet(wb, "snvs_in_ivar_not_vvc")
addWorksheet(wb, "snvs_in_vvc_not_ivar")
addWorksheet(wb, "indels_in_ivar_not_vvc")
addWorksheet(wb, "indels_in_vvc_not_ivar")
writeData(wb, "vvc_variants",vvc, borders="all")
writeData(wb, "ivar_variants",ivar, borders="all")
writeData(wb, "snvs_in_ivar_not_vvc", anti_ivar, borders = "all")
writeData(wb, "snvs_in_vvc_not_ivar", anti_vvc, borders = "all")
writeData(wb, "indels_in_ivar_not_vvc", anti_ivar_indels, borders = "all")
writeData(wb, "indels_in_vvc_not_ivar", anti_vvc_indels, borders = "all")
saveWorkbook(wb, "variant_comparison_ivar_vvc.xlsx", overwrite = TRUE)

#final thing to do is to create a new variant table that filters out these variants not found in iVar

# 15106/T556P, 1029/F75S, 2549/D582Y, 4668/S650F, 7988/A1757S
# 10605/P184H, #13035/A4V, 13394/K124E, 
# 15036/K532N, #15168/L576L, #20054/E145G

vvc_filtered <- read_xlsx("/Users/_lbashor/Dropbox/SARS-CoV-2 cat manuscript/cat ms results/R analysis cats/1_initial_processing_0.001_cats/variant_summary_0.03_processed.xlsx") %>%
  filter(!variant %in% c("T556P", "F75S", "D582Y", "S650F", 
                         "A1757S", "P184H", "A4V", "K124E", "K532N",
                         "L576L", "E145G"))

wb <- createWorkbook("variant_summary_0.03_iVar_processed.xlsx")
addWorksheet(wb, "variants")
writeData(wb, "variants", vvc_filtered, borders="all")
saveWorkbook(wb, "variant_summary_0.03_iVar_processed.xlsx", overwrite = TRUE)



