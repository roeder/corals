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

### ITS1 data wrangling --------------
# This is where we rearrange the data frame to get a view that we prefer: we 
# want a long table where each row contains the sample name, OTU # and 
# corresponding species information, along with the number of mapped reads.

# Take the raw data, 
its1_long <- its1_raw %>% 
  # select only the relevant columns by name,
  select(OTU, header, spread, one_of(samples_its1)) %>% 
  # then gather the data into a long table with the new colums "sample" (tidyr 
  # calls this the key column) and "abundance" (this is called value column)
  gather(sample, abundance, one_of(samples_its1))

# Note concerning the above: one_of() is specific syntax used to provide a 
# vector of column names to dplyr or tidyr

# Before filtering out rows with zero reads, we want to take a snapshot of how many
# reads map to each sample.
its1_sample_overview <- its1_long %>%
  group_by(sample) %>%
  summarise(hits = sum(abundance))

# Take the long table
its1_long <- its1_long %>% 
  # then filter, so that only rows with at least one read remain
  filter(abundance > 0) %>% 
  # and finally arrange the data first by sample name, then by descending number
  # of reads (for each sample).
  arrange(sample, desc(abundance))

### ITS1 output -------
# This part is a bit intense. It creates the custom output with the tables for
# each sample.
# For each sample name (outer for loop):
#
# - write a nice big heading (declared with ###) followed by two linebreaks (\n)
# - check if there are zero total reads mapping to this sample
#   |- if yes, write "no hits!\n\n" and go to next sample name
#   '- if no, continue
#
# - start this sample's table by writing the column names, a horizontal line and 
#   a linebreak
# - filter out the rows for this sample from the long data table
# - for each row in this filtered set (inner for loop):
#   * create a string that contains this row's OTU #, the first 40 characters of
#     the species description, and the number of mapped reads.
#   * write the created string followed by a linebreak.
#
# - write a final linebreak

# Declare the name of the output markdown file
# Quick overview in RStudio: Help -> Markdown Quick Reference
its1_out <- 'output/ITS1_by_sample.md'

# Create the markdown file
file.create(its1_out)

for (i in 1:nrow(its1_sample_overview)) {
  cat(paste0('### ', samples_its1[i], '\n\n'), file = its1_out, append = T)
  
  if (its1_sample_overview[its1_sample_overview$sample == samples_its1[i], 2] == 0) {
    cat('no hits!\n\n', file = its1_out, append = T)
    next()
  }
  
  cat('OTU | Info | Abundance\n--- | --- | ---\n', file = its1_out, append = T)
  
  sample_hits <- filter(its1_long, sample == samples_its1[i])
  
  for (j in 1:nrow(sample_hits)) {
    output_string <- paste(sample_hits$OTU[j], 
                           substr(sample_hits$header[j], 1, 40), 
                           sample_hits$abundance[j], 
                           sep = ' | ')
    cat(paste0(output_string, '\n'), file = its1_out, append = T)
  }
  
  
  cat('\n', file = its1_out, append = T)
}

### ITS2 data wrangling ------------
# Exactly like with ITS1

its2_long <- its2_raw %>% 
  select(OTU, header, spread, one_of(samples_its2)) %>% 
  gather(sample, abundance, one_of(samples_its2))

its2_sample_overview <- its2_long %>% 
  group_by(sample) %>% 
  summarise(hits = sum(abundance))

its2_long <- its2_long %>% 
  filter(abundance > 0) %>% 
  arrange(sample, desc(abundance))


### ITS2 output -------------
# Exactly like with ITS1
its2_out <- 'output/ITS2_by_sample.md'

file.create(its2_out)

for (i in 1:nrow(its2_sample_overview)) {
  cat(paste0('### ', samples_its2[i], '\n\n'), file = its2_out, append = T)
  
  if (its2_sample_overview[its2_sample_overview$sample == samples_its2[i], 2] == 0) {
    cat('no hits!\n\n', file = its2_out, append = T)
    next()
  }
  
  cat('OTU | Info | Abundance\n--- | --- | ---\n', file = its2_out, append = T)
  
  sample_hits <- filter(its2_long, sample == samples_its2[i])
  
  for (j in 1:nrow(sample_hits)) {
    output_string <- paste(sample_hits$OTU[j], 
                           substr(sample_hits$header[j], 1, 40), 
                           sample_hits$abundance[j], 
                           sep = ' | ')
    cat(paste0(output_string, '\n'), file = its2_out, append = T)
  }
  
  
  cat('\n', file = its2_out, append = T)
}
