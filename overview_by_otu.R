library(readxl)
library(dplyr)
library(tidyr)

### Read in data ---------
its1_raw <- read_excel('input/Linett_total_ITS1.xlsx')
its2_raw <- read_excel('input/Linett_total_ITS2.xlsx')

sample_info <- read_excel('input/ITS1 miseq data with sites.xlsx')[, 2:3]
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

group_overview$status <- rep(c('not assigned', 'healthy', 'sick'), 2)

### ITS1 by genus ---------
its1_long <- its1_raw %>% 
  select(OTU, header, one_of(samples_its1)) %>% 
  gather(sample, reads, one_of(samples_its1)) %>% 
  left_join(sample_info, by = 'sample') %>% 
  filter(diseased != '-',
         reads > 1) 

its1_genus <- its1_long %>% 
  mutate(genus = word(header, end = 1)) %>% 
  group_by(sample, genus) %>% 
  mutate(occurence = row_number()) %>% 
  filter(occurence == 1) %>% 
  left_join(group_overview[1:3, c(2, 4)], by = 'diseased') %>% 
  ungroup()

genus1 <- its1_genus %>% 
  filter(dataset == 'original') %>% 
  group_by(genus, status) %>% 
  summarise(freq = n()) %>% 
  spread(status, freq, fill = 0) %>% 
  mutate(healthy_total = group_overview$total_samples[group_overview$dataset == "original" & 
                                                        group_overview$status == 'healthy'],
         sick_total = group_overview$total_samples[group_overview$dataset == "original" & 
                                                     group_overview$status == 'sick'],
         healthy_freq = round(healthy / healthy_total, 3),
         sick_freq = round(sick / sick_total, 3))

genus3 <- its1_genus %>% 
  filter(dataset == 'extra') %>% 
  group_by(genus, status) %>% 
  summarise(freq = n()) %>% 
  spread(status, freq, fill = 0) %>% 
  mutate(healthy_total = group_overview$total_samples[group_overview$dataset == "extra" & 
                                                        group_overview$status == 'healthy'],
         sick_total = group_overview$total_samples[group_overview$dataset == "extra" & 
                                                     group_overview$status == 'sick'],
         healthy_freq = round(healthy / healthy_total, 3),
         sick_freq = round(sick / sick_total, 3))

### ITS1 by species ---------
its1_species <- its1_long %>% 
  mutate(species = word(header, end = 2)) %>% 
  group_by(sample, species) %>% 
  mutate(occurence = row_number()) %>% 
  filter(occurence == 1) %>% 
  left_join(group_overview[1:3, c(2, 4)], by = 'diseased') %>% 
  ungroup()

species1 <- its1_species %>% 
  filter(dataset == 'original') %>% 
  group_by(species, status) %>% 
  summarise(freq = n()) %>% 
  spread(status, freq, fill = 0) %>% 
  mutate(healthy_total = group_overview$total_samples[group_overview$dataset == "original" & 
                                                        group_overview$status == 'healthy'],
         sick_total = group_overview$total_samples[group_overview$dataset == "original" & 
                                                     group_overview$status == 'sick'],
         healthy_freq = round(healthy / healthy_total, 3),
         sick_freq = round(sick / sick_total, 3))

species3 <- its1_species %>% 
  filter(dataset == 'extra') %>% 
  group_by(species, status) %>% 
  summarise(freq = n()) %>% 
  spread(status, freq, fill = 0) %>% 
  mutate(healthy_total = group_overview$total_samples[group_overview$dataset == "extra" & 
                                                        group_overview$status == 'healthy'],
         sick_total = group_overview$total_samples[group_overview$dataset == "extra" & 
                                                     group_overview$status == 'sick'],
         healthy_freq = round(healthy / healthy_total, 3),
         sick_freq = round(sick / sick_total, 3))

### ITS2 by genus ---------
its2_long <- its2_raw %>% 
  select(OTU, header, one_of(samples_its2)) %>% 
  gather(sample, reads, one_of(samples_its2)) %>% 
  left_join(sample_info, by = 'sample') %>% 
  filter(diseased != '-',
         reads > 1) 

its2_genus <- its2_long %>% 
  mutate(genus = word(header, end = 1)) %>% 
  group_by(sample, genus) %>% 
  mutate(occurence = row_number()) %>% 
  filter(occurence == 1) %>% 
  left_join(group_overview[1:3, c(2, 4)], by = 'diseased') %>% 
  ungroup()

genus2 <- its2_genus %>% 
  filter(dataset == 'original') %>% 
  group_by(genus, status) %>% 
  summarise(freq = n()) %>% 
  spread(status, freq, fill = 0) %>% 
  mutate(healthy_total = group_overview$total_samples[group_overview$dataset == "original" & 
                                                        group_overview$status == 'healthy'],
         sick_total = group_overview$total_samples[group_overview$dataset == "original" & 
                                                     group_overview$status == 'sick'],
         healthy_freq = round(healthy / healthy_total, 3),
         sick_freq = round(sick / sick_total, 3))

genus4 <- its2_genus %>% 
  filter(dataset == 'extra') %>% 
  group_by(genus, status) %>% 
  summarise(freq = n()) %>% 
  spread(status, freq, fill = 0) %>% 
  mutate(healthy_total = group_overview$total_samples[group_overview$dataset == "extra" & 
                                                        group_overview$status == 'healthy'],
         sick_total = group_overview$total_samples[group_overview$dataset == "extra" & 
                                                     group_overview$status == 'sick'],
         healthy_freq = round(healthy / healthy_total, 3),
         sick_freq = round(sick / sick_total, 3))

### ITS2 by species ---------
its2_species <- its2_long %>% 
  mutate(species = word(header, end = 2)) %>% 
  group_by(sample, species) %>% 
  mutate(occurence = row_number()) %>% 
  filter(occurence == 1) %>% 
  left_join(group_overview[1:3, c(2, 4)], by = 'diseased') %>% 
  ungroup()

species2 <- its2_species %>% 
  filter(dataset == 'original') %>% 
  group_by(species, status) %>% 
  summarise(freq = n()) %>% 
  spread(status, freq, fill = 0) %>% 
  mutate(healthy_total = group_overview$total_samples[group_overview$dataset == "original" & 
                                                        group_overview$status == 'healthy'],
         sick_total = group_overview$total_samples[group_overview$dataset == "original" & 
                                                     group_overview$status == 'sick'],
         healthy_freq = round(healthy / healthy_total, 3),
         sick_freq = round(sick / sick_total, 3))

species4 <- its2_species %>% 
  filter(dataset == 'extra') %>% 
  group_by(species, status) %>% 
  summarise(freq = n()) %>% 
  spread(status, freq, fill = 0) %>% 
  mutate(healthy_total = group_overview$total_samples[group_overview$dataset == "extra" & 
                                                        group_overview$status == 'healthy'],
         sick_total = group_overview$total_samples[group_overview$dataset == "extra" & 
                                                     group_overview$status == 'sick'],
         healthy_freq = round(healthy / healthy_total, 3),
         sick_freq = round(sick / sick_total, 3))

### Output to xlsx files ----------
library(openxlsx)

l_genus <- list("ITS1_orig" = genus1,
                "ITS2_orig" = genus2,
                "ITS1_extra" = genus3,
                "ITS2_extra" = genus4)

write.xlsx(l_genus, file = "output/new_overview_by_genus.xlsx")

l_species <- list("ITS1_orig" = species1,
                "ITS2_orig" = species2,
                "ITS1_extra" = species3,
                "ITS2_extra" = species4)

write.xlsx(l_species, file = "output/new_overview_by_species.xlsx")
