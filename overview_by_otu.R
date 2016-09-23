### Load packages ---------
library(readxl)  # reading xlsx files
library(dplyr)  # data frame manipulation
library(tidyr)  # wide data <-> long data transformations
library(openxlsx)  # writing xlsx files
library(stringr)  # word() function for extracting species/genus

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

### Prepare long tables ----------
# First we transform the wide raw data into a long table, so every row contains
# a sample name, OTU id and corresponding header information as well as the 
# number of mapped reads.

# Take the raw, wide data frame,
its1_long <- its1_raw %>% 
  # select only the columns containing OTU id, species information (header), and
  # the columns containing sample names,
  select(OTU, header, one_of(samples_its1)) %>% 
  # then gather the data into a long table with the new colums "sample" (tidyr 
  # calls this the key column) and "reads" (this is called value column),
  gather(sample, reads, one_of(samples_its1)) %>% 
  # then add information about each sample by joining with the sample info 
  # lookup table
  left_join(sample_info, by = 'sample') %>% 
  # We also add the more descriptive disease status information by left-
  # joining with a subset of the group overview lookup table
  left_join(group_overview[1:3, c(2, 4)], by = 'diseased') %>%
  # and finally only keep rows that have at least two mapped reads and a disease
  # status that is 1 or 0 (i.e. not "-")
  filter(diseased != '-',
         reads > 1) 

# In the next block, we
# aggregate the number of hits for each OTU - grouped by disease status. The 
# data frame resulting from this block will directly go into one the output
# excel sheet.

# Take the prepared long table with genus information,
otu_orig <- its1_long %>% 
  # continue only with samples from the original dataset,
  filter(dataset == 'original') %>% 
  # group by OTU and disease status
  group_by(OTU, status) %>% 
  # and aggregate the data according to this grouping: summarise() from dplyr 
  # provides summaries of data groups. Here, we just use n() to count all 
  # possible cases within this grouping. 
  #
  # This changes the data structure substantially: we just have three columns - 
  # OTU number, disease status and number of samples of this disease status 
  # where the OTU shows up (with at least two mapped reads).
  summarise(freq = n()) %>% 
  # In one of the last steps before exporting the table, we turn the long table
  # into a wide table again - wide tables are nice as overview tables and that's 
  # what we want. The spread() function from tidyr does the transformation.
  spread(status, freq, fill = 0) %>% 
  # We add some more information before exporting the data. We grab the total
  # number of samples for this dataset and each disease status from the group
  # overview lookup table
  mutate(healthy_total = group_overview$total_samples[group_overview$dataset == "original" & 
                                                        group_overview$status == 'healthy'],
         sick_total = group_overview$total_samples[group_overview$dataset == "original" & 
                                                     group_overview$status == 'sick'],
         # and calculate frequencies: number of samples with "hits" divided by
         # total number of samples.
         healthy_freq = round(healthy / healthy_total, 3),
         sick_freq = round(sick / sick_total, 3)) %>% 
  # Finally, we also add the header information as extra column again by joining 
  # with a subset of the raw data which only contains the OTU and header 
  # columns.
  left_join(its1_raw[, colnames(its1_raw) == "header" | colnames(its1_raw) == "OTU"], by = "OTU")

# This is the same as above, only with the extra dataset instead of original.
otu_extra <- its1_long %>% 
  filter(dataset == 'extra') %>% 
  group_by(OTU, status) %>% 
  summarise(freq = n()) %>% 
  spread(status, freq, fill = 0) %>% 
  mutate(healthy_total = group_overview$total_samples[group_overview$dataset == "extra" & 
                                                        group_overview$status == 'healthy'],
         sick_total = group_overview$total_samples[group_overview$dataset == "extra" & 
                                                     group_overview$status == 'sick'],
         healthy_freq = round(healthy / healthy_total, 3),
         sick_freq = round(sick / sick_total, 3)) %>% 
  left_join(its1_raw[, colnames(its1_raw) == "header" | colnames(its1_raw) == "OTU"], by = "OTU")

### Output to xlsx files ----------
# We use write.xlsx() from the openxlsx package to export our overviews to excel
# workbooks. First, we create a list containing the data frames we'd like to 
# export. The names of the list members (in quotes) will be the sheet names in 
# the excel file. 
l_otu <- list("original" = otu_orig,
                "extra" = otu_extra)

# Then we just export to the desired file by providing the list object:
write.xlsx(l_otu, file = "output/overview_by_otu.xlsx")
