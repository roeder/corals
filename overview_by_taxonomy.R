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
# a sample name, OTU id and corresponding species information as well as the 
# number of mapped reads.
#
# In subsequent code blocks, we will construct overviews for specific ITS 
# regions and original vs extra dataset based on the long table created here.

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
  # and finally only keep rows that have at least two mapped reads and a disease
  # status that is 1 or 0 (i.e. not "-")
  filter(diseased != '-',
         reads > 1) 

# ITS2 works exactly the same. There will be a lot of copy-pasted blocks with  
# slight differences in this script.
its2_long <- its2_raw %>% 
  select(OTU, header, one_of(samples_its2)) %>% 
  gather(sample, reads, one_of(samples_its2)) %>% 
  left_join(sample_info, by = 'sample') %>% 
  filter(diseased != '-',
         reads > 1) 

### ITS1 by genus ---------
# Next up, we want to aggregate data by genus information (it could just as well
# be species information, the order is not important). Conceptually, what we 
# want to do is find out if an OTU corresponding to a specific genus shows up in
# a sample at least once (with at least two mapped reads, see above).

# Take the long table,
its1_genus <- its1_long %>% 
  # add a column that has the name of the genus by grabbing the first word from
  # the 'header' column, (this is a dumb computer approach - we will also have
  # the genus 'Uncultured' in here)
  mutate(genus = word(header, end = 1)) %>% 
  # then provide some grouping information. From now on, commands that can use 
  # grouping information, will operate on subgroups of data that have the same
  # sample name and the same genus information.
  group_by(sample, genus) %>% 
  # Add a column called 'occurence'. Here we just count how many different OTUs
  # corresponding to the same genus show up (with at least two reads, still) in 
  # one sample. This uses the grouping we just defined.
  # 
  # For example, sample 2.4.h has 23 different Millepora OTUs. So for the rows 
  # with sample == 2.4.h and genus == Millepora, the 'occurence'column will have
  # integers from 1 to 23.
  #
  # Try running this block only up to and including the mutate() command (stop 
  # before the pipe operator %>%) and look at the dataset. It should be pretty
  # intuitive.
  mutate(occurence = row_number()) %>% 
  # As we only want to know if a genus occurs at least once within a sample, we
  # continue with only the first occurences of each genus in a sample.
  filter(occurence == 1) %>% 
  # Finally, we add the more descriptive disease status information by left-
  # joining with a subset of the group overview lookup table
  left_join(group_overview[1:3, c(2, 4)], by = 'diseased') %>%
  # and remove the grouping information again (pretty sure this isn't really
  # necessary. I'm doing it here to be safe because there will be another 
  # grouping later on).
  ungroup()

# The next block is the last unique procedure in this script, the rest will be
# slightly altered copy-pasted versions of things we've already seen. Here, we
# aggregate the number of hits for each genus - grouped by disease status. The 
# data frame resulting from this block will directly go into one of the output
# excel sheets (it's ITS1 data from the original dataset aggregated by genus).

# Take the prepared long table with genus information,
genus1 <- its1_genus %>% 
  # continue only with samples from the original dataset,
  filter(dataset == 'original') %>% 
  # group by genus and disease status
  group_by(genus, status) %>% 
  # and aggregate the data according to this grouping: summarise() from dplyr 
  # provides summaries of data groups. Here, we just use n() to count all 
  # possible cases within this grouping. 
  #
  # This changes the data structure substantially: we just have three columns - 
  # genus name, disease status and number of samples of this disease status 
  # where the genus shows up (with at least OTU and at least two mapped reads).
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
         sick_freq = round(sick / sick_total, 3))

# This is the same as above, only with the extra dataset instead of original.
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
# Aggregation by species works in a way analogous to aggregation by genus. The 
# only difference is we grab the first two words from the header column instead 
# of only the first one.

its1_species <- its1_long %>% 
  mutate(species = word(header, end = 2)) %>%  # end = 2 gets first two words!
  group_by(sample, species) %>% 
  mutate(occurence = row_number()) %>% 
  filter(occurence == 1) %>% 
  left_join(group_overview[1:3, c(2, 4)], by = 'diseased') %>% 
  ungroup()

# Same procedure as seen before! ITS1/original/species
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

# ITS1/extra/species
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
its2_genus <- its2_long %>% 
  mutate(genus = word(header, end = 1)) %>% 
  group_by(sample, genus) %>% 
  mutate(occurence = row_number()) %>% 
  filter(occurence == 1) %>% 
  left_join(group_overview[1:3, c(2, 4)], by = 'diseased') %>% 
  ungroup()

# ITS2/original/genus
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

# ITS2/extra/genus
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

# ITS2/original/species
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

# ITS2/extra/species
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
# We use write.xlsx() from the openxlsx package to export our overviews to excel
# workbooks. First, we create a list containing the data frames we'd like to 
# export. The names of the list members (in quotes) will be the sheet names in 
# the excel file. 
l_genus <- list("ITS1_orig" = genus1,
                "ITS2_orig" = genus2,
                "ITS1_extra" = genus3,
                "ITS2_extra" = genus4)

# Then we just export to the desired file by providing the list object:
write.xlsx(l_genus, file = "output/new_overview_by_genus.xlsx")

# We do the same for the tables aggregated by species.
l_species <- list("ITS1_orig" = species1,
                  "ITS2_orig" = species2,
                  "ITS1_extra" = species3,
                  "ITS2_extra" = species4)

write.xlsx(l_species, file = "output/new_overview_by_species.xlsx")
