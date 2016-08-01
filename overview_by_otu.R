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
  summarise(total_samples = n())

### Turn into long tables
its1_long <- its1_raw %>% 
  select(OTU, header, one_of(samples_its1)) %>% 
  gather(sample, reads, one_of(samples_its1)) %>% 
  left_join(sample_info, by = 'sample') %>% 
  filter(diseased != '-',
         reads > 1) 

its1_agg <- its1_long %>% 
  group_by(OTU, header, dataset, diseased) %>% 
  summarise(total_hits = n()) %>% 
  left_join(group_overview, by = c("diseased", "dataset")) %>% 
  mutate(normalised_hits = total_hits / total_samples)

its1_genus <- its1_agg %>% 
  mutate(genus = word(header)) %>% 
  group_by(genus, dataset, diseased) %>% 
  summarise(total_hits = mean(total_hits),
            total_samples = mean(total_samples))

its1_genus_hits <- its1_genus %>% 
  select(-total_samples) %>%
  spread(diseased, total_hits, sep = '_hits_', fill = 0)

its1_genus_samples <- its1_genus %>% 
  select(-total_hits) %>%
  spread(diseased, total_samples, sep = '_samples_', fill = 0)

its1_genus_final <- its1_genus_hits %>% 
  left_join(its1_genus_samples, by = c("genus", "dataset")) %>% 
  mutate(freq_0 = round(diseased_hits_0 / diseased_samples_0, 2),
         freq_1 = round(diseased_hits_1 / diseased_samples_1, 2),
         output_0 = paste0(diseased_hits_0, ' (', freq_0, ')'),
         output_1 = paste0(diseased_hits_1, ' (', freq_1, ')'))
