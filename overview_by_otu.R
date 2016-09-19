### Load packages ---------
library(readxl)  # reading xlsx files
library(dplyr)  # data frame manipulation
library(tidyr)  # wide data <-> long data transformations

### Read in data ---------

# Read in raw data from the two big excel sheets
its1_raw <- read_excel('input/Linett_total_ITS1.xlsx')
its2_raw <- read_excel('input/Linett_total_ITS2.xlsx')

# Read in columns 2 and 3 from the excel sheet containing disease status for
# the samples
sample_info <- read_excel('input/ITS1 miseq data with sites.xlsx')[, 2:3]
colnames(sample_info) <- c('diseased', 'sample')  # shorter column names

### Fix sample names -----------
# This part is a bit fuzzy and wasn't optimised to be short or easy on the eyes.
# It's a bunch of code that does the job in terms of cleaning up the sample 
# names.

# Get ITS1 sample names from the sample_info lookup table (it is filled with 
# info on ITS1 samples).
samples_its1 <- sample_info$sample

# Get ITS2 names from column names in the big raw table. 
samples_its2 <- colnames(its2_raw)[28:197]

# The PCR blanks have messed up, repetitive names. We give them unique names
# instead. They now have names from PCR_neg_1 to PCR_neg_10
samples_its2[grep('PCR', samples_its2, ignore.case = T)] <- 
  paste0('PCR_neg_', 1:10)

# Add PCR_neg_10 to the lookup table. (It already has sample names up to and 
# including pcr_neg_9)
sample_info <- rbind(sample_info, c("-", 'PCR_neg_10'))

# Run a script with a bunch of lines that clean up the formatting of sample 
# names: trim away whitespace, convert to lower letters, remove unnecessary 
# spaces and fix typos. This needs to be done for every analysis, so putting 
# those commands in an extra script cleans up the analysis scripts.
source('normalise_sample_names.R')

# Replace the column names that are sample names with the now cleaned up ones.
colnames(its1_raw)[28:197] <- samples_its1
colnames(its2_raw)[28:197] <- samples_its2

### Prepare dataset information ---------------

# Find the first sample name containing a parenthesis character. This is the
# first sample in the "extra" dataset.
dataset_split <- grep('(', sample_info$sample, fixed = T)[1]

# Add a column to the sample info lookup table containing the name of the
# dataset this sample belongs to.
sample_info$dataset <- c(rep('original', dataset_split - 1), 
                         rep('extra', (nrow(sample_info) - dataset_split + 1)))

# Another tweak to the lookup table: set rows without disease information (they
# are NAs right now) to status "-". Now 1 is diseased, 0 healthy, - not assigned
sample_info$diseased[is.na(sample_info$diseased)] <- '-'

# Create a lookup table containing the total number of samples corresponding to 
# original/extra dataset as well as each of the three disease statuses.
group_overview <- sample_info %>% 
  group_by(dataset, diseased) %>% 
  summarise(total_samples = n())

# Add a more descriptive column for the disease status. This will be used for 
# plotting.
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
