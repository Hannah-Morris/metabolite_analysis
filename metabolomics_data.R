#clear working directory:
rm(list=ls())


library(tidyverse)

# Set working directory
setwd("~/Library/CloudStorage/OneDrive-UniversityofKent")

# Load data
df <- read.csv("chenomx_data_all.csv", check.names = FALSE)
meta <- read.csv("gregarine_count.csv", check.names = FALSE)

# Replace blank cells with NA
df[df == ""] <- NA

# Keep metabolite metadata separately
met_info <- df %>%
  dplyr::select(Metabolite, Formula, `HMDB Accession Number`)

# Keep only sample columns
df_data <- df %>%
  dplyr::select(-Metabolite, -Formula, -`HMDB Accession Number`)

# Set metabolite names as row names
rownames(df_data) <- df$Metabolite

# Convert all sample columns to numeric
df_data <- df_data %>%
  mutate(across(everything(), as.numeric))

# Calculate prevalence:
# number of samples where each metabolite is detected (>0 and not NA)
prevalence <- rowSums(!is.na(df_data) & df_data > 0)

# Make prevalence table
prev_df <- data.frame(
  Metabolite = rownames(df_data),
  Prevalence = prevalence,
  row.names = NULL
) %>%
  left_join(met_info, by = "Metabolite") %>%
  arrange(desc(Prevalence))

# View top metabolites
print(head(prev_df, 50))

# Plot top 50 most prevalent metabolites
top_prev <- prev_df %>%
  slice_head(n = 50)

ggplot(top_prev, aes(x = reorder(Metabolite, Prevalence), y = Prevalence)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Top 50 Most Prevalent Metabolites",
    x = "Metabolite",
    y = "Number of Samples Detected"
  ) +
  theme_minimal()

# OPTIONAL: calculate mean abundance across samples
abundance <- rowMeans(df_data, na.rm = TRUE)

abund_df <- data.frame(
  Metabolite = rownames(df_data),
  MeanAbundance = abundance,
  row.names = NULL
) %>%
  left_join(met_info, by = "Metabolite") %>%
  arrange(desc(MeanAbundance))

# View top abundant metabolites
print(head(abund_df, 20))

# Plot top 20 most abundant metabolites
top_abund <- abund_df %>%
  slice_head(n = 20)

ggplot(top_abund, aes(x = reorder(Metabolite, MeanAbundance), y = MeanAbundance)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Top 20 Most Abundant Metabolites",
    x = "Metabolite",
    y = "Mean Concentration"
  ) +
  theme_minimal()




#################### PCA ANALYSIS #############################################

# Transpose data (samples = rows, metabolites = columns)
df_t <- as.data.frame(t(df_data))
df_t$SampleID <- rownames(df_t)

# Replace NA with 0 (required for PCA)
df_t[is.na(df_t)] <- 0


#################### CLEAN METADATA #################################

meta <- meta %>%
  rename(
    SampleID = Sample_ID,
    Host = Species,
    Gregarine_status = Gregarines
  )

meta <- meta %>%
  filter(Host != "Water")

# Join metadata
pca_df <- df_t %>%
  left_join(meta, by = "SampleID")

#################### PCA ANALYSIS #############################################

# Transpose data (samples = rows)
df_t <- as.data.frame(t(df_data))
df_t$SampleID <- rownames(df_t)

# Replace NA with 0
df_t[is.na(df_t)] <- 0

# Join metadata
pca_df <- df_t %>%
  left_join(meta, by = "SampleID")

pca_df <- df_t %>%
  left_join(meta, by = "SampleID") %>%
  filter(
    !is.na(Gregarine_status),
    Gregarine_status %in% c("Positive", "Negative"),
    !is.na(Host),
    !Host %in% c("Unknown", "Water", "NA")
  )

# Extract only numeric metabolite data
pca_numeric <- pca_df %>%
  dplyr::select(where(is.numeric))

# REMOVE columns with zero variance 
pca_numeric <- pca_numeric[, apply(pca_numeric, 2, var) != 0]

# Run PCA
pca <- prcomp(pca_numeric, scale. = TRUE)

# Extract scores
pca_scores <- as.data.frame(pca$x)
pca_scores$SampleID <- pca_df$SampleID
pca_scores$Host <- pca_df$Host
pca_scores$Gregarine_status <- pca_df$Gregarine_status

# Variance explained
var_explained <- summary(pca)$importance[2,]

# Plot
ggplot(pca_scores,
       aes(x = PC1, y = PC2,
           colour = Host,
           shape = Gregarine_status)) +
  geom_point(size = 4) +
  stat_ellipse(aes(group = Host), linetype = 2, linewidth = 0.8) +
  labs(
    title = "PCA of Metabolite Profiles",
    x = paste0("PC1 (", round(var_explained[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2] * 100, 1), "%)")
  ) +
  scale_colour_brewer(palette = "Set1") +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    legend.title = element_text(face = "bold")
  )


###############################################################################
############################# PLS- DA #########################################
###############################################################################
library(mixOmics)
library(tidyverse)
library(dplyr)

# Metabolite matrix
X <- pca_df %>%
  dplyr::select(where(is.numeric)) %>%
  as.matrix()

# Remove zero-variance metabolites
X <- X[, apply(X, 2, var) > 0]

# Response variable
Y <- as.factor(pca_df$Host)

# Check dimensions
dim(X)
table(Y)

# 2- components first:
plsda_model <- plsda(
  X,
  Y,
  ncomp = 2
)

#Score plot:
plotIndiv(
  plsda_model,
  comp = c(1,2),
  group = Y,
  ellipse = TRUE,
  legend = TRUE,
  ind.names = FALSE,
  title = "PLS-DA: Shrimp Species"
)

#Cross-validation
set.seed(123)

perf_plsda <- perf(
  plsda_model,
  validation = "Mfold",
  folds = 5,
  nrepeat = 100,
  progressBar = TRUE
)

#Classification error:
perf_plsda$error.rate

#plot:
plot(perf_plsda)




#### gregarine PLS-DA
Y_inf <- factor(pca_df$Gregarine_status)

plsda_inf <- plsda(
  X,
  Y_inf,
  ncomp = 2
)


plotIndiv(
  plsda_inf,
  comp = c(1,2),
  group = Y_inf,
  ellipse = TRUE,
  legend = TRUE,
  ind.names = FALSE,
  title = "PLS-DA: Gregarine Status"
)





#### VIP scores
vip_scores <- vip(plsda_model)
head(vip_scores)

dim(vip_scores)

str(vip_scores)


vip_df <- data.frame(
  Metabolite = rownames(vip_scores),
  VIP = vip_scores[,1]
)

head(vip_df)


vip_df <- vip_df[order(-vip_df$VIP), ]

head(vip_df, 20)




vip_df %>%
  slice_max(VIP, n = 20) %>%
  ggplot(aes(
    x = reorder(Metabolite, VIP),
    y = VIP
  )) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(
    title = "Top 20 VIP Metabolites",
    x = "",
    y = "VIP Score"
  )

############################Heatmap of VIP metabolites ########################

library(dplyr)
library(pheatmap)



#Select top (10)
top_vips <- vip_df %>%
  arrange(desc(VIP)) %>%
  slice_head(n = 10)

top_metabolites <- top_vips$Metabolite


#Create etabolite matrix
#Extract only top 25
heat_data <- pca_df %>%
  dplyr::select(all_of(top_metabolites))

#Convert to matrix
heat_matrix <- as.matrix(heat_data)

#Sample ID as row name
rownames(heat_matrix) <- pca_df$SampleID



#Z-score scaling (red=higher, blue=lower)
heat_matrix_scaled <- scale(heat_matrix)

#species
annotation_col <- data.frame(
  Species = pca_df$Host
)

rownames(annotation_col) <- pca_df$SampleID

#Species +gregarine status:
annotation_col <- data.frame(
  Species = pca_df$Host,
  Gregarine = pca_df$Gregarine_status
)

rownames(annotation_col) <- pca_df$SampleID


#Group by species
sample_order <- order(pca_df$Host)

heat_matrix_scaled <- heat_matrix_scaled[sample_order, ]

annotation_col <- annotation_col[sample_order, ]

#plot
pheatmap(
  t(heat_matrix_scaled),
  
  annotation_col = annotation_col,
  
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  
  show_colnames = FALSE,
  
  main = "Top 10 VIP Metabolites"
)








###################
### Box plot for top 10 metabolites

library(tidyverse)
library(ggplot2)
ggplot(
  plot_long,
  aes(
    x = Host,
    y = Concentration,
    fill = Host
  )
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15,
    size = 1.2,
    alpha = 0.7
  ) +
  facet_wrap(
    ~ Metabolite,
    scales = "free_y",
    ncol = 2
  ) +
  theme_classic() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )








####################### WITH METADATA #########################################
rm(list = ls())

library(tidyverse)

# Set working directory
setwd("~/Library/CloudStorage/OneDrive-UniversityofKent")

####################### LOAD DATA ##############################################
df <- read.csv("chenomx_data_all.csv", check.names = FALSE)
meta <- read.csv("gregarine_count.csv", check.names = FALSE)

library(tidyverse)

#################### CLEAN METADATA ####################
meta <- meta[, colnames(meta) != ""]
colnames(meta)[1] <- "SampleID"

meta <- meta %>%
  mutate(
    Gregarine_status = case_when(
      Gregarines == "Positive" ~ "Positive",
      Gregarines == "Negative" ~ "Negative",
      TRUE ~ NA_character_
    )
  )

#################### LONG METABOLITE TABLE ####################
df_long <- df %>%
  pivot_longer(
    cols = -c(Metabolite, Formula, `HMDB Accession Number`),
    names_to = "SampleID",
    values_to = "Abundance"
  ) %>%
  mutate(Abundance = as.numeric(Abundance))

#################### JOIN ####################
data <- df_long %>%
  left_join(meta, by = "SampleID")

#################### TOP 20: GREGARINE POSITIVE ####################
top20_greg_pos <- data %>%
  filter(Gregarine_status == "Positive") %>%
  group_by(Metabolite) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 20)

p1 <- ggplot(top20_greg_pos,
             aes(x = reorder(Metabolite, mean_abundance), y = mean_abundance)) +
  geom_col(fill = "#E64B35FF") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Top 20 metabolites in Gregarine-positive samples",
    x = "Metabolite",
    y = "Mean abundance"
  )

#################### TOP 20: GREGARINE NEGATIVE ####################
top20_greg_neg <- data %>%
  filter(Gregarine_status == "Negative") %>%
  group_by(Metabolite) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 20)

p2 <- ggplot(top20_greg_neg,
             aes(x = reorder(Metabolite, mean_abundance), y = mean_abundance)) +
  geom_col(fill = "#4DBBD5FF") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Top 20 metabolites in Gregarine-negative samples",
    x = "Metabolite",
    y = "Mean abundance"
  )


library(patchwork)
p1 + p2



ggplot(
  data %>% filter(Gregarine_status %in% c("Positive", "Negative")),
  aes(x = Gregarine_status, y = Abundance, fill = Gregarine_status)
) +
  geom_boxplot() +
  facet_wrap(~ Metabolite, scales = "free_y") +
  theme_minimal() +
  labs(
    y = "Concentration (µM)",  # adjust unit
    x = ""
  )


meta <- meta %>%
  mutate(
    Blastocystis_status = case_when(
      Blastocystis == "Positive" ~ "Positive",
      Blastocystis == "Negative" ~ "Negative",
      TRUE ~ NA_character_
    )
  )

data <- df_long %>%
  left_join(meta, by = "SampleID")




library(tidyverse)

library(tidyverse)

#################### TOP 20 FROM EACH GROUP ####################

# Top 20 positive
top20_pos <- data %>%
  filter(Gregarine_status == "Positive") %>%
  group_by(Metabolite) %>%
  summarise(mean_pos = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_pos)) %>%
  slice_head(n = 20)

# Top 20 negative
top20_neg <- data %>%
  filter(Gregarine_status == "Negative") %>%
  group_by(Metabolite) %>%
  summarise(mean_neg = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_neg)) %>%
  slice_head(n = 20)

#################### MERGE METABOLITE LIST ####################

# Get combined metabolite list (union of both)
top_metabs <- union(top20_pos$Metabolite, top20_neg$Metabolite)

#################### BUILD FINAL TABLE ####################

greg_table_combined <- data %>%
  filter(Gregarine_status %in% c("Positive", "Negative")) %>%
  filter(Metabolite %in% top_metabs) %>%
  group_by(Metabolite, Gregarine_status) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = Gregarine_status,
    values_from = mean_abundance
  )

#################### ADD LOG2 FOLD CHANGE ####################

greg_table_combined <- greg_table_combined %>%
  mutate(
    log2FC = log2((Positive + 1e-9) / (Negative + 1e-9))
  )

#################### SORT BY EFFECT SIZE ####################

greg_table_combined <- greg_table_combined %>%
  arrange(desc(abs(log2FC)))

greg_table_combined


View(greg_table_combined)


#################### TOP 20: BLASTOCYSTIS POSITIVE ####################
top20_blast_pos <- data %>%
  filter(Blastocystis_status == "Positive") %>%
  group_by(Metabolite) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 20)

p3 <- ggplot(top20_blast_pos,
             aes(x = reorder(Metabolite, mean_abundance), y = mean_abundance)) +
  geom_col(fill = "#00A087FF") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Top 20 metabolites in Blastocystis-positive samples",
    x = "Metabolite",
    y = "Mean abundance"
  )

#################### TOP 20: BLASTOCYSTIS NEGATIVE ####################
top20_blast_neg <- data %>%
  filter(Blastocystis_status == "Negative") %>%
  group_by(Metabolite) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 20)

p4 <- ggplot(top20_blast_neg,
             aes(x = reorder(Metabolite, mean_abundance), y = mean_abundance)) +
  geom_col(fill = "#3C5488FF") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Top 20 metabolites in Blastocystis-negative samples",
    x = "Metabolite",
    y = "Mean abundance"
  )
p3 + p4

