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

original_samples <- sample_info$sample[1:dataset_split - 1]
extra_samples <- sample_info$sample[dataset_split:nrow(sample_info)]

control_samples <- c(grep('pcr_neg_', sample_info$sample, fixed = T, value = T),
                     grep('ex_con_', sample_info$sample, fixed = T, value = T))

### PCA 1: ITS1 & original data ---------
df1 <- its1_raw %>% 
  select(OTU, one_of(setdiff(intersect(samples_its1, original_samples), 
                             control_samples))) %>% 
  arrange(OTU)

matrix1 <- t(data.matrix(df1[, -1]))

colnames(matrix1) <- paste0("OTU_", df1$OTU)

pca1 <- prcomp(matrix1, center = F)

pca_var1 <- pca1$sdev ^ 2

cum_var1 <- cumsum(pca_var1 / sum(pca_var1))

### PCA 2: ITS1 & extra data ---------
df2 <- its1_raw %>% 
  select(OTU, one_of(setdiff(intersect(samples_its1, extra_samples), 
                             control_samples))) %>% 
  arrange(OTU)

matrix2 <- t(data.matrix(df2[, -1]))

colnames(matrix2) <- paste0("OTU_", df2$OTU)

pca2 <- prcomp(matrix2, center = F)

pca_var2 <- pca2$sdev ^ 2

cum_var2 <- cumsum(pca_var2 / sum(pca_var2))

### PCA 3: ITS2 & original data ---------
df3 <- its2_raw %>% 
  select(OTU, one_of(setdiff(intersect(samples_its2, original_samples), 
                             control_samples))) %>% 
  arrange(OTU)

matrix3 <- t(data.matrix(df3[, -1]))

colnames(matrix3) <- paste0("OTU_", df3$OTU)

pca3 <- prcomp(matrix3, center = F)

pca_var3 <- pca3$sdev ^ 2

cum_var3 <- cumsum(pca_var3 / sum(pca_var3))

### PCA 4: ITS2 & extra data ---------
df4 <- its2_raw %>% 
  select(OTU, one_of(setdiff(intersect(samples_its2, extra_samples), 
                             control_samples))) %>% 
  arrange(OTU)

matrix4 <- t(data.matrix(df4[, -1]))

colnames(matrix4) <- paste0("OTU_", df4$OTU)

pca4 <- prcomp(matrix4, center = F)

pca_var4 <- pca4$sdev ^ 2

cum_var4 <- cumsum(pca_var4 / sum(pca_var4))

### Plot prep -------------
plot_data1 <- as.data.frame(pca1$x[, 1:2])
plot_data1$sample <- rownames(plot_data1)

plot_data1 <- left_join(plot_data1, sample_info, by = 'sample')

### Plots -----------
library(ggplot2)

ggplot(plot_data1, aes(x = PC1, y = PC2, colour = diseased)) +
  geom_point() +
  ggtitle(paste0("PCA for ITS1 & original data\n", 
                 round(cum_var1[2] * 100, 2), 
                 " % of total variation captured")) +
  guides(colour = guide_legend(title = NULL)) +
  scale_colour_discrete(labels = c("healthy", "sick")) +
  theme_light() 
