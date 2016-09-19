# We want to fix the character vectors samples_its1 and samples_its2 as well as 
# the data frame column sample_info$sample.

library(stringr)  # package provides the str_trim() function 

# Trim whitespace in the beginning and end, also convert to lower-case.
samples_its1 <- str_trim(tolower(samples_its1))
samples_its2 <- str_trim(tolower(samples_its2))
sample_info$sample <- str_trim(tolower(sample_info$sample))

# Remove extra spaces and underscores.
samples_its1 <- gsub(' (', '(', samples_its1, fixed = T)
samples_its1 <- gsub('_(', '(', samples_its1, fixed = T)
sample_info$sample <- gsub(' (', '(', sample_info$sample, fixed = T)
sample_info$sample <- gsub('_(', '(', sample_info$sample, fixed = T)
samples_its2 <- gsub(' (', '(', samples_its2, fixed = T)
samples_its2 <- gsub(' (', '(', samples_its2, fixed = T) # twice for 100(07.11)

# Fix random typos.
samples_its2 <- gsub('t v', 'tv', samples_its2, fixed = T)
samples_its2 <- gsub(' con ', '_con_', samples_its2, fixed = T)
