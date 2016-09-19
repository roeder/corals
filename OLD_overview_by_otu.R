### OLD STUPID VERSION - DO NOT USE

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

### Turn ITS1 into long tables ---------
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

### Aggregate ITS1 by genus ---------------
its1_genus <- its1_agg %>% 
  mutate(genus = word(header)) %>% 
  group_by(genus, dataset, diseased) %>% 
  summarise(total_hits = round(mean(total_hits), 2),
            total_samples = round(mean(total_samples), 2))

its1_genus_hits <- its1_genus %>% 
  select(-total_samples) %>%
  spread(diseased, total_hits, sep = '_hits_', fill = 0)

its1_genus_samples <- its1_genus %>% 
  select(-total_hits) %>%
  spread(diseased, total_samples, sep = '_samples_', fill = 0)

its1_genus_freq <- its1_genus_hits %>% 
  left_join(its1_genus_samples, by = c("genus", "dataset")) %>% 
  mutate(freq_h = round(diseased_hits_0 / diseased_samples_0, 2),
         freq_s = round(diseased_hits_1 / diseased_samples_1, 2)) %>% 
  arrange(desc(freq_s))

its1_genus_freq$freq_h[is.na(its1_genus_freq$freq_h)] <- 0
its1_genus_freq$freq_s[is.na(its1_genus_freq$freq_s)] <- 0

# its1_genus_freq <- its1_genus_freq %>% 
#   mutate(output_0 = paste0(diseased_hits_0, ' (', freq_0, ')'),
#          output_1 = paste0(diseased_hits_1, ' (', freq_1, ')'))

names(its1_genus_freq) <- gsub('diseased_', '', names(its1_genus_freq), fixed = T)
names(its1_genus_freq) <- gsub('_0', '_h', names(its1_genus_freq), fixed = T)
names(its1_genus_freq) <- gsub('_1', '_s', names(its1_genus_freq), fixed = T)

### Aggregate ITS1 by species ---------------
its1_species <- its1_agg %>% 
  mutate(species = word(header, end = 2)) %>% 
  group_by(species, dataset, diseased) %>% 
  summarise(total_hits = round(mean(total_hits), 2),
            total_samples = round(mean(total_samples), 2))

its1_species_hits <- its1_species %>% 
  select(-total_samples) %>%
  spread(diseased, total_hits, sep = '_hits_', fill = 0)

its1_species_samples <- its1_species %>% 
  select(-total_hits) %>%
  spread(diseased, total_samples, sep = '_samples_', fill = 0)

its1_species_freq <- its1_species_hits %>% 
  left_join(its1_species_samples, by = c("species", "dataset")) %>% 
  mutate(freq_h = round(diseased_hits_0 / diseased_samples_0, 2),
         freq_s = round(diseased_hits_1 / diseased_samples_1, 2)) %>% 
  arrange(desc(freq_s))

its1_species_freq$freq_h[is.na(its1_species_freq$freq_h)] <- 0
its1_species_freq$freq_s[is.na(its1_species_freq$freq_s)] <- 0

# its1_species_freq <- its1_species_freq %>% 
#   mutate(output_0 = paste0(diseased_hits_0, ' (', freq_0, ')'),
#          output_1 = paste0(diseased_hits_1, ' (', freq_1, ')'))

names(its1_species_freq) <- gsub('diseased_', '', names(its1_species_freq), fixed = T)
names(its1_species_freq) <- gsub('_0', '_h', names(its1_species_freq), fixed = T)
names(its1_species_freq) <- gsub('_1', '_s', names(its1_species_freq), fixed = T)

### Turn ITS2 into long tables ---------
its2_long <- its2_raw %>% 
  select(OTU, header, one_of(samples_its2)) %>% 
  gather(sample, reads, one_of(samples_its2)) %>% 
  left_join(sample_info, by = 'sample') %>% 
  filter(diseased != '-',
         reads > 1) 

its2_agg <- its2_long %>% 
  group_by(OTU, header, dataset, diseased) %>% 
  summarise(total_hits = n()) %>% 
  left_join(group_overview, by = c("diseased", "dataset")) %>% 
  mutate(normalised_hits = total_hits / total_samples)

### Aggregate ITS2 by genus ---------------
its2_genus <- its2_agg %>% 
  mutate(genus = word(header)) %>% 
  group_by(genus, dataset, diseased) %>% 
  summarise(total_hits = round(mean(total_hits), 2),
            total_samples = round(mean(total_samples), 2))

its2_genus_hits <- its2_genus %>% 
  select(-total_samples) %>%
  spread(diseased, total_hits, sep = '_hits_', fill = 0)

its2_genus_samples <- its2_genus %>% 
  select(-total_hits) %>%
  spread(diseased, total_samples, sep = '_samples_', fill = 0)

its2_genus_freq <- its2_genus_hits %>% 
  left_join(its2_genus_samples, by = c("genus", "dataset")) %>% 
  mutate(freq_h = round(diseased_hits_0 / diseased_samples_0, 2),
         freq_s = round(diseased_hits_1 / diseased_samples_1, 2)) %>% 
  arrange(desc(freq_s))

its2_genus_freq$freq_h[is.na(its2_genus_freq$freq_h)] <- 0
its2_genus_freq$freq_s[is.na(its2_genus_freq$freq_s)] <- 0

# its2_genus_freq <- its2_genus_freq %>% 
#   mutate(output_0 = paste0(diseased_hits_0, ' (', freq_0, ')'),
#          output_1 = paste0(diseased_hits_1, ' (', freq_1, ')'))

names(its2_genus_freq) <- gsub('diseased_', '', names(its2_genus_freq), fixed = T)
names(its2_genus_freq) <- gsub('_0', '_h', names(its2_genus_freq), fixed = T)
names(its2_genus_freq) <- gsub('_1', '_s', names(its2_genus_freq), fixed = T)

### Aggregate ITS2 by species ---------------
its2_species <- its2_agg %>% 
  mutate(species = word(header, end = 2)) %>% 
  group_by(species, dataset, diseased) %>% 
  summarise(total_hits = round(mean(total_hits), 2),
            total_samples = round(mean(total_samples), 2))

its2_species_hits <- its2_species %>% 
  select(-total_samples) %>%
  spread(diseased, total_hits, sep = '_hits_', fill = 0)

its2_species_samples <- its2_species %>% 
  select(-total_hits) %>%
  spread(diseased, total_samples, sep = '_samples_', fill = 0)

its2_species_freq <- its2_species_hits %>% 
  left_join(its2_species_samples, by = c("species", "dataset")) %>% 
  mutate(freq_h = round(diseased_hits_0 / diseased_samples_0, 2),
         freq_s = round(diseased_hits_1 / diseased_samples_1, 2)) %>% 
  arrange(desc(freq_s))

its2_species_freq$freq_h[is.na(its2_species_freq$freq_h)] <- 0
its2_species_freq$freq_s[is.na(its2_species_freq$freq_s)] <- 0

# its2_species_freq <- its2_species_freq %>% 
#   mutate(output_0 = paste0(diseased_hits_0, ' (', freq_0, ')'),
#          output_1 = paste0(diseased_hits_1, ' (', freq_1, ')'))

names(its2_species_freq) <- gsub('diseased_', '', names(its2_species_freq), fixed = T)
names(its2_species_freq) <- gsub('_0', '_h', names(its2_species_freq), fixed = T)
names(its2_species_freq) <- gsub('_1', '_s', names(its2_species_freq), fixed = T)

### Output to xlsx files ----------
library(openxlsx)

l_genus <- list("ITS1_orig" = filter(its1_genus_freq, dataset == 'original'),
                "ITS2_orig" = filter(its2_genus_freq, dataset == 'original'),
                "ITS1_extra" = filter(its1_genus_freq, dataset == 'extra'),
                "ITS2_extra" = filter(its2_genus_freq, dataset == 'extra'))

write.xlsx(l_genus, file = "output/overview_by_genus.xlsx")

l_species <- list("ITS1_orig" = filter(its1_species_freq, dataset == 'original'),
                  "ITS2_orig" = filter(its2_species_freq, dataset == 'original'),
                  "ITS1_extra" = filter(its1_species_freq, dataset == 'extra'),
                  "ITS2_extra" = filter(its2_species_freq, dataset == 'extra'))

write.xlsx(l_species, file = "output/overview_by_species.xlsx")