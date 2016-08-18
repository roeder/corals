library(readxl)
library(dplyr)
library(tidyr)

### Read in data ---------
its1_raw <- read_excel('Linett_total_ITS1.xlsx')
its2_raw <- read_excel('Linett_total_ITS2.xlsx')

sample_info <- read_excel('ITS1 miseq data with sites.xlsx')[, 1:3]
colnames(sample_info) <- c('site', 'diseased', 'sample')

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

sample_info$site <- str_to_title(sample_info$site)

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

rel_var3 <- pca_var3 / sum(pca_var3)

cum_var3 <- cumsum(rel_var3)

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
plot_data1$individual <- str_sub(plot_data1$sample, end = -3)

plot_data1 <- left_join(plot_data1, sample_info, by = 'sample')

plot_data2 <- as.data.frame(pca2$x[, 1:2])
plot_data2$sample <- rownames(plot_data2)

plot_data2 <- left_join(plot_data2, sample_info, by = 'sample')

plot_data3 <- as.data.frame(pca3$x[, c(1,3)])
plot_data3$sample <- rownames(plot_data3)
plot_data3$individual <- str_sub(plot_data3$sample, end = -3)

plot_data3 <- left_join(plot_data3, sample_info, by = 'sample')

plot_data4 <- as.data.frame(pca4$x[, 1:2])
plot_data4$sample <- rownames(plot_data4)

plot_data4 <- left_join(plot_data4, sample_info, by = 'sample')

### Plots by disease status -----------
library(ggplot2)

ggplot(plot_data1, aes(x = PC1, y = PC2, colour = diseased)) +
  geom_point() +
  ggtitle(paste0("PCA for ITS1 & original data\n", 
                 round(cum_var1[2] * 100, 2), 
                 " % of total variation captured")) +
  guides(colour = guide_legend(title = NULL)) +
  scale_colour_manual(values = c("#0072B2", "#D55E00"),
                      labels = c("healthy", "sick")) +
  theme_light() 

ggplot(plot_data2, aes(x = PC1, y = PC2, colour = diseased)) +
  geom_point() +
  ggtitle(paste0("PCA for ITS1 & extra data\n", 
                 round(cum_var2[2] * 100, 2), 
                 " % of total variation captured")) +
  guides(colour = guide_legend(title = NULL)) +
  scale_colour_manual(values = c("#999999", "#0072B2", "#D55E00"),
                      labels = c("unlabeled", "healthy", "sick")) +
  theme_light() 

ggplot(plot_data3, aes(x = PC1, y = PC3, colour = diseased)) +
  geom_point() +
  ggtitle(paste0("PCA for ITS2 & original data\n", 
                 round((rel_var3[1] + rel_var3[3]) * 100, 2), 
                 " % of total variation captured")) +
  guides(colour = guide_legend(title = NULL)) +
  scale_colour_manual(values = c("#0072B2", "#D55E00"),
                      labels = c("healthy", "sick")) +
  theme_light() 

ggplot(plot_data4, aes(x = PC1, y = PC2, colour = diseased)) +
  geom_point() +
  ggtitle(paste0("PCA for ITS2 & extra data\n", 
                 round(cum_var4[2] * 100, 2), 
                 " % of total variation captured")) +
  guides(colour = guide_legend(title = NULL)) +
  scale_colour_manual(values = c("#999999", "#0072B2", "#D55E00"),
                      labels = c("unlabeled", "healthy", "sick")) +
  theme_light() 

### Plots by site & disease status -----------

ggplot(plot_data1, aes(x = PC1, y = PC2, colour = site, shape = diseased)) +
  geom_point() +
  ggtitle(paste0("PCA for ITS1 & original data\n", 
                 round(cum_var1[2] * 100, 2), 
                 " % of total variation captured")) +
  theme_light() 

ggplot(plot_data2, aes(x = PC1, y = PC2, colour = site, shape = diseased)) +
  geom_point() +
  ggtitle(paste0("PCA for ITS1 & extra data\n", 
                 round(cum_var2[2] * 100, 2), 
                 " % of total variation captured")) +
  theme_light() 

ggplot(plot_data3, aes(x = PC1, y = PC3, colour = site, shape = diseased)) +
  geom_point() +
  ggtitle(paste0("PCA for ITS2 & original data\n", 
                 round((rel_var3[1] + rel_var3[3]) * 100, 2), 
                 " % of total variation captured")) +
  theme_light() 

ggplot(plot_data4, aes(x = PC1, y = PC2, colour = site, shape = diseased)) +
  geom_point() +
  ggtitle(paste0("PCA for ITS2 & extra data\n", 
                 round(cum_var4[2] * 100, 2), 
                 " % of total variation captured")) +
  theme_light() 

### Plots by pairing -------------

ggplot(plot_data1, aes(x = PC1, y = PC2, colour = individual, shape = diseased)) +
  geom_point() +
  ggtitle(paste0("PCA for ITS1 & original data\n", 
                 round(cum_var1[2] * 100, 2), 
                 " % of total variation captured")) +
  theme_light() 

ggplot(plot_data3, aes(x = PC1, y = PC3, colour = individual, shape = diseased)) +
  geom_point() +
  ggtitle(paste0("PCA for ITS2 & original data\n", 
                 round((rel_var3[1] + rel_var3[3]) * 100, 2), 
                 " % of total variation captured")) +
  theme_light() 