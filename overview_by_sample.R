library(readxl)
library(dplyr)
library(tidyr)

### Read in data ---------
its1_raw <- read_excel('Linett_total_ITS1.xlsx')
# its2_data <- read_excel('Linett_total_ITS2.xlsx')

sample_names_its1 <- read_excel('sample_names_ITS1.xlsx')$Samples
colnames(its1_raw)[28:197] <- sample_names_its1

its1_data <- its1_raw %>% 
  select(OTU, header, spread, one_of(sample_names_its1)) %>% 
  gather(sample, abundance, one_of(sample_names_its1))

its1_sample_overview <- its1_data %>% 
  group_by(sample) %>% 
  summarise(hits = sum(abundance))

its1_data <- its1_data %>% 
  filter(abundance > 0) %>% 
  arrange(sample, desc(abundance))

file.create('ITS1_by_sample.md')

for (i in 1:nrow(its1_sample_overview)) {
  cat(paste0('### ', its1_sample_overview$sample[i], '\n\n'), file = 'ITS1_by_sample.md', append = T)
  
  if (its1_sample_overview$hits[i] == 0) {
    cat('no hits!\n\n', file = 'ITS1_by_sample.md', append = T)
    next()
  }
  
  cat('OTU | Info | Abundance\n--- | --- | ---\n', file = 'ITS1_by_sample.md', append = T)
  
  sample_hits <- filter(its1_data, sample == its1_sample_overview$sample[i])
  
  for (j in 1:nrow(sample_hits)) {
    output_string <- paste(sample_hits$OTU[j], 
                           substr(sample_hits$header[j], 1, 40), 
                           sample_hits$abundance[j], 
                           sep = ' | ')
    cat(paste0(output_string, '\n'), file = 'ITS1_by_sample.md', append = T)
  }
  
  
  cat('\n', file = 'ITS1_by_sample.md', append = T)
}