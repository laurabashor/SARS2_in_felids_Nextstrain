#the purpose of this script is to clean up the dataset output by SNPGenie

#load libraries
library(tidyverse)
library(openxlsx)

#load in selection analysis tables from SnpGenie analysis
setwd("/Users/_lbashor/Dropbox/SARS-CoV-2 cat manuscript/cat ms results/R analysis cats/selection_cats")
dir<-"./snpgenie_tables/"
subdir<-list.dirs(dir, recursive=F)

#get dir/file names
ff<-do.call(rbind, lapply(subdir, function(x) {
  ff<-list.files(x, "\\.txt$", include.dirs = FALSE, full.names = TRUE)
  data.frame(dir=basename(x), file=basename(ff), 
             fullpath=ff, stringsAsFactors=F)
}))

#helper function from Mr Flick https://gist.github.com/MrFlick/11407464
read.stack <- function(files, ..., select=NULL, extra=NULL, reader=read.table) {
  dd<-data.frame()
  if(!is.null(extra)) {
    stopifnot(is.list(extra))
    stopifnot(all(sapply(extra, length)==length(files)))
  }
  for(i in 1:length(files)) {
    d<-reader(files[i], ...)
    if(!is.null(select)) {
      stopifnot(all(select %in% names(d)))
      d<-d[, select]
    }
    if(!is.null(extra)) {
      d<-do.call(cbind, c(list(d), lapply(extra, '[', i)))
    }
    if(nrow(dd)>0) {
      dd<-rbind(dd, d)
    } else {
      dd<-d
    }
  }
  dd
}

#read pop summary files in as a dataframe and clean up
ff_pop <- ff %>%
  filter(file == "population_summary.txt")

df_pop <- read.stack(ff_pop$fullpath, extra=list(file=ff_pop$file, dir=ff_pop$dir))

#helper function from https://stackoverflow.com/questions/32054368/use-first-row-data-as-column-names-in-r
header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

df_pop <- header.true(df_pop)

names(df_pop)

df_pop$file <- sapply(strsplit(df_pop$file, "/\\s*"), tail, 1)
df_pop$file <- gsub("\\..*","", df_pop$file)

#helpful deletion of all the headers from each table which were introduced as rows
#function from https://stackoverflow.com/questions/39106128/delete-every-evenuneven-row-from-a-dataset
toDelete <- seq(1, nrow(df_pop), 2)
df_pop <- df_pop[ toDelete ,]

df_pop <- df_pop %>%
  mutate(dataset_ID = file) %>%
  select(dataset_ID, everything()) %>%
  select(-c(population_summary.txt, 
            Cat_1_R.MN985325.bam_SNPGenie_Results,
            file)) %>%
  filter(!dataset_ID %in% c("Cat_2", "Cat_2_R",
                            "Cat_3", "Cat_3_R",
                            "Cat_4", "Cat_4_R"))

#read product summary files in as a dataframe and clean up

ff_prod <- ff %>%
  filter(file == "product_results.txt")

df_prod <- read.stack(ff_prod$fullpath, extra=list(file=ff_prod$file, dir=ff_prod$dir))

df_prod <- header.true(df_prod)

names(df_prod)

df_prod$file <- sapply(strsplit(df_prod$file, "/\\s*"), tail, 1)
df_prod$file <- gsub("\\..*","", df_prod$file)

df_prod <- df_prod %>%
  mutate(dataset_ID = file) %>%
  filter(dataset_ID != "file") %>%
  select(dataset_ID, everything()) %>%
  select(-c(product_results.txt, 
            Cat_1_R.MN985325.bam_SNPGenie_Results,
            file)) %>%
  filter(!dataset_ID %in% c("Cat_2", "Cat_2_R",
                            "Cat_3", "Cat_3_R",
                            "Cat_4", "Cat_4_R"))


# make a column indicating replicates
df_pop$replicate = str_extract(df_pop$dataset_ID, "_R")
df_pop$replicate[is.na(df_pop$replicate)] <- 1
df_pop$replicate[df_pop$replicate == "_R"] <- 2

df_prod$replicate = str_extract(df_prod$dataset_ID, "_R")
df_prod$replicate[is.na(df_prod$replicate)] <- 1
df_prod$replicate[df_prod$replicate == "_R"] <- 2

#now create a column "ID_rep" that has the full ID with replicate
#and remove the replicate indicator from the dataset_ID column

df_pop <- df_pop %>% 
  mutate(ID_rep = dataset_ID)

df_prod <- df_prod %>% 
  mutate(ID_rep = dataset_ID) 

df_pop$dataset_ID <- gsub("_R", "", df_pop$dataset_ID)
df_prod$dataset_ID <- gsub("_R", "", df_prod$dataset_ID)


#save cleaned files
write_csv(df_pop, "population_summary_cleaned.csv")
write_csv(df_prod, "gene_product_cleaned.csv")

