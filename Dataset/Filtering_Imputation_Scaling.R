library(tidyverse)
library(data.table)
library(patchwork)
library(mvdalab)

##########################
## HANDLING MISSINGNESS ##
##########################

# First, we will load in the datasets and plot their properties
microbiome <- fread("Dataset/Original/16S_Edata.csv")
metaproteomics <- fread("Dataset/Original/Metaproteomics.csv")
metab_neg <- fread("Dataset/Original/Metabolomics_Negative.csv")
metab_pos <- fread("Dataset/Original/Metabolomics_Positive.csv")

# Build a function to calculate missingness
calc_missingness <- function(edata, edata_col, view_name) {
  data.table(
    Feature = edata[[edata_col]],
    View = view_name,
    `00wk Missingness` = apply(edata, 1, function(x) {sum(is.na(x)[2:4])}),
    `post-00wk Missingness` = apply(edata, 1, function(x) {sum(is.na(x)[5:13])})
  )
}

# Get all missingness calculations
missingness <- bind_rows(
  calc_missingness(metaproteomics, "Feature", "metaproteomics"),
  calc_missingness(metab_pos, "Feature", "metabolomics positive"),
  calc_missingness(metab_neg, "Feature", "metabolomics negative")
)

# Plot missingness 
missingnessplot <- missingness %>%
  pivot_longer(3:4) %>%
  rename(Type = name, Count = value) %>%
  mutate(
    Count = as.character(Count),
    Kept = map2_chr(Type, Count, function(x, y) {
      if (x == "00wk Missingness" & y >= 2) {
        return("Removed")
      } else if (x == "post-00wk Missingness" & y >= 5) {
        return("Removed")
      } else {return("Kept")}
    })
  ) %>%
  group_by(View, Type, Kept, Count) %>%
  summarize(Frequency = n()) %>%
  mutate(Frequency = log10(Frequency)) %>%
  ggplot(aes(x = Count, y = Frequency, fill = Kept)) +
    geom_bar(stat = "identity") +
    xlab("Number Missing") +
    ylab("Log10 Frequency") +
    scale_fill_manual(values = c("black", "red")) +
    facet_grid(rows = vars(Type), cols = vars(View), scales = "free") +
    theme_bw() +
    theme(legend.title = element_blank())
missingnessplot

# Determine impact of missingness filtering
toFilter <- missingness %>%
  mutate(Filter = `00wk Missingness` >= 2 | `post-00wk Missingness` >= 5)
toFilter %>% 
  group_by(View) %>%
  summarize(
    Count = n() - sum(Filter),
    Filtered = sum(Filter) / n()
  )
bio_filter <- toFilter %>% filter(Filter) %>% select(Feature) %>% unlist()

# Filter
microbiome <- microbiome %>% 
  filter(Feature %in% bio_filter == FALSE) %>%
  data.frame()
metaproteomics <- metaproteomics %>% 
  filter(Feature %in% bio_filter == FALSE) %>%
  data.frame()
metab_pos <- metab_pos %>% 
  filter(Feature %in% bio_filter == FALSE) %>%
  data.frame()
metab_neg <- metab_neg %>% 
  filter(Feature %in% bio_filter == FALSE) %>%
  data.frame()

# Data has already been log transformed and normalized

####################################################
## VISUALIZE DISTRIBUTIONS & CORRELATION PATTERNS ##
####################################################

correlationplot <- rbind(
  microbiome %>% mutate(View = "16S"), 
  metaproteomics %>% mutate(View = "metaproteomics"), 
  metab_pos %>% mutate(View = "metabolomics positive"), 
  metab_neg %>% mutate(View = "metabolomics negative")
) %>%
  relocate(View) %>%
  group_by(View) %>%
  nest() %>%
  mutate(data = map(data, function(x) {
    data <- x[,2:ncol(x)]
    data[is.na(data)] <- 0
    data %>%
      cor() %>%
      data.frame() %>%
      mutate(Sample1 = row.names(.)) %>%
      relocate(Sample1) %>%
      pivot_longer(2:ncol(.)) %>%
      rename(Sample2 = name, Correlation = value) %>%
      mutate(Sample2 = gsub("X", "", Sample2)) 
  })) %>%
  unnest(cols = c(data)) %>%
  mutate(
    Sample1 = gsub("109_|X", "", Sample1),
    Sample2 = gsub("109_", "", Sample2)
  ) %>%
  ggplot(aes(x = Sample1, y = Sample2, fill = Correlation)) + 
    geom_tile() +
    xlab("") + 
    ylab("") + 
    theme(axis.text.x = element_text(angle = 90),
          plot.title = element_text(hjust = 0.5)) +
    facet_wrap(.~View)
correlationplot

##############################################
## IMPUTATION WITH EXPECTATION MAXIMIZATION ##
##############################################

em_impute <- function(omics) {
  imp_res <- imputeEM(data = omics[,-1], impute.ncomps = 3)
  omics[,-1] <- imp_res$Imputed.DataFrames[[1]]
  return(omics)
}
metaproteomics <- em_impute(metaproteomics)
metab_pos <- em_impute(metab_pos)
metab_neg <- em_impute(metab_neg)

######################
## FIX COLUMN NAMES ##
######################

names <- colnames(microbiome) %>% gsub("X", "", .)

microbiome <- data.table(microbiome)
colnames(microbiome) <- names

metaproteomics <- data.table(metaproteomics)
colnames(metaproteomics) <- names

metab_pos <- data.table(metab_pos)
colnames(metab_pos) <- names

metab_neg <- data.table(metab_neg)
colnames(metab_neg) <- names

#####################
## POST IMPUTATION ##
#####################

# For reference, look at correlation patterns after imputation
postimpcorrelationplot <- rbind(
  microbiome %>% mutate(View = "16S"), 
  metaproteomics %>% mutate(View = "metaproteomics"), 
  metab_pos %>% mutate(View = "metabolomics positive"), 
  metab_neg %>% mutate(View = "metabolomics negative")
) %>%
  relocate(View) %>%
  group_by(View) %>%
  nest() %>%
  mutate(data = map(data, function(x) {
    data <- x[,2:ncol(x)]
    data[is.na(data)] <- 0
    data %>%
      cor() %>%
      data.frame() %>%
      mutate(Sample1 = row.names(.)) %>%
      relocate(Sample1) %>%
      pivot_longer(2:ncol(.)) %>%
      rename(Sample2 = name, Correlation = value) %>%
      mutate(Sample2 = gsub("X", "", Sample2)) 
  })) %>%
  unnest(cols = c(data)) %>%
  mutate(
    Sample1 = gsub("109_|X", "", Sample1),
    Sample2 = gsub("109_", "", Sample2)
  ) %>%
  ggplot(aes(x = Sample1, y = Sample2, fill = Correlation)) + 
  geom_tile() +
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(angle = 90),
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(.~View)
postimpcorrelationplot

#################
## MAKE FIGURE ##
#################

Figure1 <- missingnessplot / (correlationplot + postimpcorrelationplot) + plot_annotation(tag_levels = "A")
Figure1

###########
## SCALE ##
###########

# Scale the ranges of each dataset

## Microbiome
microbiome %>%
  pivot_longer(2:ncol(.)) %>%
  mutate(value = scale(value)[,1]) %>%
  pivot_wider(id_cols = Feature) %>%
  fwrite("Dataset/Scaled/16S_Edata.csv", quote = F, row.names = F)

## Metaproteomics
metaproteomics %>%
  pivot_longer(2:ncol(.)) %>%
  mutate(value = scale(value)[,1]) %>%
  pivot_wider(id_cols = Feature) %>%
  fwrite("Dataset/Scaled/Metaproteomics.csv", quote = F, row.names = F)

## Metabolomics Positive ##
metab_pos %>%
  pivot_longer(2:ncol(.)) %>%
  mutate(value = scale(value)[,1]) %>%
  pivot_wider(id_cols = Feature) %>%
  fwrite("Dataset/Scaled/Metabolomics_Positive.csv", quote = F, row.names = F)

## Metabolomics Negative ##
metab_neg %>%
  pivot_longer(2:ncol(.)) %>%
  mutate(value = scale(value)[,1]) %>%
  pivot_wider(id_cols = Feature) %>%
  fwrite("Dataset/Scaled/Metabolomics_Negative.csv", quote = F, row.names = F)
