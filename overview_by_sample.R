library(readxl)
library(dplyr)
library(tidyr)

### Read in data ---------
its1_raw <- read_excel('Linett_total_ITS1.xlsx')
its2_raw <- read_excel('Linett_total_ITS2.xlsx')

# ITS1: Fixed names are loaded from separate spreadsheet!
samples_its1 <- read_excel('ITS1 miseq data with sites.xlsx')$Samples
colnames(its1_raw)[28:197] <- samples_its1

# ITS2: Fix names by appending letter indices a, b, c & d to duplicate names
# The duplicate names are PCR neg_1 and PCR neg_2
colnames(its2_raw)[grep('neg_1', colnames(its2_raw), ignore.case = T)] <- 
  paste0(colnames(its2_raw)[grep('neg_1', colnames(its2_raw), ignore.case = T)], '_', letters[1:4])

colnames(its2_raw)[grep('neg_2', colnames(its2_raw), ignore.case = T)] <- 
  paste0(colnames(its2_raw)[grep('neg_2', colnames(its2_raw), ignore.case = T)], '_', letters[1:4])

samples_its2 <- colnames(its2_raw)[28:197]

# ITS 1 process
its1_data <- its1_raw %>% 
  select(OTU, header, spread, one_of(samples_its1)) %>% 
  gather(sample, abundance, one_of(samples_its1))

its1_sample_overview <- its1_data %>% 
  group_by(sample) %>% 
  summarise(hits = sum(abundance))

its1_data <- its1_data %>% 
  filter(abundance > 0) %>% 
  arrange(sample, desc(abundance))

file.create('ITS1_by_sample.md')

for (i in 1:nrow(its1_sample_overview)) {
  cat(paste0('### ', samples_its1[i], '\n\n'), file = 'ITS1_by_sample.md', append = T)
  
  if (its1_sample_overview[its1_sample_overview$sample == samples_its1[i], 2] == 0) {
    cat('no hits!\n\n', file = 'ITS1_by_sample.md', append = T)
    next()
  }
  
  cat('OTU | Info | Abundance\n--- | --- | ---\n', file = 'ITS1_by_sample.md', append = T)
  
  sample_hits <- filter(its1_data, sample == samples_its1[i])
  
  for (j in 1:nrow(sample_hits)) {
    output_string <- paste(sample_hits$OTU[j], 
                           substr(sample_hits$header[j], 1, 40), 
                           sample_hits$abundance[j], 
                           sep = ' | ')
    cat(paste0(output_string, '\n'), file = 'ITS1_by_sample.md', append = T)
  }
  
  
  cat('\n', file = 'ITS1_by_sample.md', append = T)
}

# ITS 2 process
its2_data <- its2_raw %>% 
  select(OTU, header, spread, one_of(samples_its2)) %>% 
  gather(sample, abundance, one_of(samples_its2))

its2_sample_overview <- its2_data %>% 
  group_by(sample) %>% 
  summarise(hits = sum(abundance))

its2_data <- its2_data %>% 
  filter(abundance > 0) %>% 
  arrange(sample, desc(abundance))

file.create('ITS2_by_sample.md')

for (i in 1:nrow(its2_sample_overview)) {
  cat(paste0('### ', samples_its2[i], '\n\n'), file = 'ITS2_by_sample.md', append = T)
  
  if (its2_sample_overview[its2_sample_overview$sample == samples_its2[i], 2] == 0) {
    cat('no hits!\n\n', file = 'ITS2_by_sample.md', append = T)
    next()
  }
  
  cat('OTU | Info | Abundance\n--- | --- | ---\n', file = 'ITS2_by_sample.md', append = T)
  
  sample_hits <- filter(its2_data, sample == samples_its2[i])
  
  for (j in 1:nrow(sample_hits)) {
    output_string <- paste(sample_hits$OTU[j], 
                           substr(sample_hits$header[j], 1, 40), 
                           sample_hits$abundance[j], 
                           sep = ' | ')
    cat(paste0(output_string, '\n'), file = 'ITS2_by_sample.md', append = T)
  }
  
  
  cat('\n', file = 'ITS2_by_sample.md', append = T)
}
