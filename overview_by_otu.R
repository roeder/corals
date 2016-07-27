library(readxl)
library(dplyr)
library(tidyr)

### Read in data ---------
its1_raw <- read_excel('Linett_total_ITS1.xlsx')
its2_raw <- read_excel('Linett_total_ITS2.xlsx')

sample_info <- read_excel('ITS1 miseq data with sites.xlsx')[, 2:3]
colnames(sample_info) <- c('diseased', 'sample')

# ITS1: Fixed names are loaded from separate spreadsheet!
samples_its1 <- sample_info$sample

# ITS2: Fix duplicate PCR neg names by appending new indices in order of
# appearance
colnames(its2_raw)[grep('PCR', colnames(its2_raw), ignore.case = T)] <- 
  paste0('PCR_neg_', 1:10)

samples_its2 <- colnames(its2_raw)[28:197]

sample_info <- rbind(sample_info, c("-", 'PCR_neg_10'))
sample_info$diseased[is.na(sample_info$diseased)] <- '-'

source('normalise_sample_names.R')

colnames(its1_raw)[28:197] <- samples_its1
colnames(its2_raw)[28:197] <- samples_its2

dataset_split <- grep('(', sample_info$sample, fixed = T)[1]
sample_info$dataset <- c(rep('original', dataset_split - 1), 
                         rep('extra', (nrow(sample_info) - dataset_split + 1)))

group_overview <- sample_info %>% 
  group_by(dataset, diseased) %>% 
  summarise(samples = n())

n_sick <- sum(sample_info$diseased == '1', na.rm = T)
n_healthy <- sum(sample_info$diseased == '0', na.rm = T)

### Turn into long tables
its1_long <- its1_raw %>% 
  select(OTU, header, one_of(samples_its1)) %>% 
  gather(sample, hits, one_of(samples_its1)) %>% 
  left_join(sample_info, by = 'sample') %>% 
  filter(diseased != '-') 

its1_agg <- its1_long %>% 
  group_by(OTU, diseased) %>% 
  summarise(total_hits = n(),
            total_abundance = sum(hits))