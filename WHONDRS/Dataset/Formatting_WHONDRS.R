library(tidyverse)
library(data.table)
library(patchwork)
library(pmartR)
library(mvdalab)

setwd("~/Git_Repos/multiomics_integration_tools/WHONDRS/Dataset/Original/")

##########################
## HANDLING MISSINGNESS ##
##########################

clean_names = function(x) {
  x %>% 
    gsub(pattern = "_", replacement = ".", fixed = TRUE) %>% 
    gsub(pattern = "-", replacement = ".", fixed = TRUE)
}

# First, we will load in the datasets and plot their properties
metab = fread("GCMS_Metabolomics_edata.txt")
colnames(metab) = clean_names(colnames(metab))
lip_pos = fread("LCMS_Lipidomics_pos_edata.txt")
colnames(lip_pos) = clean_names(colnames(lip_pos))
lip_neg = fread("LCMS_Lipidomics_neg_edata.txt")
colnames(lip_neg) = clean_names(colnames(lip_neg))

# Load metadata
fdata = fread("Everything_fdata.csv") %>%
  mutate(SampleID = clean_names(SampleID))

# Build a function to calculate missingness - 16 none and 18 surrounding
calc_missingness <- function(edata, viewname) {
  edata %>%
    pivot_longer(-1) %>%
    rename(SampleID = name) %>%
    left_join(fdata) %>%
    group_by(Feature, Vegetation) %>%
    summarize(
      Count = sum(is.na(value))
    ) %>%
    mutate(View = viewname)
}

# Get all missingness calculations
missingness = bind_rows(
  calc_missingness(metab, "metabolomics"),
  calc_missingness(lip_pos, "lipidomics positive"),
  calc_missingness(lip_neg, "lipidomics negative")
)

# Plot missingness 
missingnessplot = missingness %>%
  ungroup() %>%
  group_by(View, Vegetation, Count) %>%
  summarize(Frequency = n()) %>%
  mutate(
    Kept = map2_chr(Vegetation, Count, function(x, y) {
      if (x == "None" & y > 8) {
        return("Removed")
      } else if (x == "Surrounding" & y > 9) {
        return("Removed")
      } else {return("Kept")}
    })
  ) %>%
  ggplot(aes(x = Count, y = Frequency, fill = Kept)) +
    geom_bar(stat = "identity") +
    xlab("Number Missing") +
    ylab("Frequency") +
    scale_fill_manual(values = c("black", "red")) +
    facet_grid(rows = vars(Vegetation), cols = vars(View), scales = "free") +
    theme_bw() +
    theme(legend.title = element_blank())
missingnessplot

# Determine impact of missingness filtering
toFilter <- missingness %>%
  pivot_wider(id_cols = c(Feature, View), names_from = Vegetation, values_from = Count) %>%
  mutate(Filter = None > 8 | Surrounding > 9)
toFilter %>% 
  group_by(View) %>%
  summarize(
    Count = n() - sum(Filter),
    Filtered = sum(Filter) / n()
  )
bio_filter = toFilter %>% filter(Filter) %>% select(Feature) %>% unlist()

# Filter
metab = metab %>% 
  filter(Feature %in% bio_filter == FALSE) %>%
  data.frame()
lip_pos = lip_pos %>% 
  filter(Feature %in% bio_filter == FALSE) %>%
  data.frame()
lip_neg = lip_neg %>% 
  filter(Feature %in% bio_filter == FALSE) %>%
  data.frame()

# Log transform and normalize using pmartR

# P-Value is 0.6788
lipDat1 = as.lipidData(
  e_data = lip_pos, f_data = fdata, edata_cname = "Feature", fdata_cname = "SampleID"
) %>%
  edata_transform(data_scale = "log2") %>%
  group_designation(main_effects = "Vegetation") %>%
  normalize_global("ppp_rip", "mad", list(ppp_rip = list(ppp = 0.5, rip = 0.2)), apply_norm = TRUE) 
plot(dim_reduction(lipDat1))

# P-Value is 0.8630
lipDat2 = as.lipidData(
  e_data = lip_neg, f_data = fdata, edata_cname = "Feature", fdata_cname = "SampleID"
) %>%
  edata_transform(data_scale = "log2") %>%
  group_designation(main_effects = "Vegetation") %>%
  normalize_global("ppp_rip", "mad", list(ppp_rip = list(ppp = 0.5, rip = 0.2)), apply_norm = TRUE) 
plot(dim_reduction(lipDat2))

# P-Value is 0.9724
metabDat = as.metabData(
  e_data = metab, f_data = fdata, edata_cname = "Feature", fdata_cname = "SampleID"
) %>%
  edata_transform(data_scale = "log2") %>%
  group_designation(main_effects = "Vegetation") %>%
  normalize_global("ppp_rip", "mad", list(ppp_rip = list(ppp = 0.5, rip = 0.2)), apply_norm = TRUE) 
plot(dim_reduction(metabDat))

##############################################
## IMPUTATION WITH EXPECTATION MAXIMIZATION ##
##############################################

em_impute <- function(omics) {
  imp_res <- imputeEM(data = omics[,-1], impute.ncomps = 3)
  omics[,-1] <- imp_res$Imputed.DataFrames[[1]]
  return(omics)
}
metabolomics = em_impute(metabDat$e_data)
lipidomics_pos = em_impute(lipDat1$e_data)
lipidomics_neg = em_impute(lipDat2$e_data)

##############################
## SCALE & FIX COLUMN NAMES ##
##############################

fdata = fdata %>%
  mutate(Group = map_chr(Tag, function(x) {
    split = x %>% strsplit("_") %>% unlist()
    paste0(substr(split[1], 1, 1), split[4])
  })) 
fwrite(fdata, "../Scaled/FData.csv", quote = F, row.names = F)

metabolomics %>%
  pivot_longer(-1) %>%
  mutate(value = scale(value)[,1]) %>%
  rename(SampleID = name) %>%
  left_join(fdata %>% select(SampleID, Group)) %>%
  select(-SampleID) %>% 
  pivot_wider(id_cols = Feature, names_from = Group) %>%
  fwrite("../Scaled/Metabolomics.txt", sep = "\t", quote = F, row.names = F)

lipidomics_pos %>%
  pivot_longer(-1) %>%
  mutate(value = scale(value)[,1]) %>%
  rename(SampleID = name) %>%
  left_join(fdata %>% select(SampleID, Group)) %>%
  select(-SampleID) %>% 
  pivot_wider(id_cols = Feature, names_from = Group) %>%
  fwrite("../Scaled/Lipidomics_Positive.txt", sep = "\t", quote = F, row.names = F)

lipidomics_neg %>%
  pivot_longer(-1) %>%
  mutate(value = scale(value)[,1]) %>%
  rename(SampleID = name) %>%
  left_join(fdata %>% select(SampleID, Group)) %>%
  select(-SampleID) %>% 
  pivot_wider(id_cols = Feature, names_from = Group) %>%
  fwrite("../Scaled/Lipidomics_Negative.txt", sep = "\t", quote = F, row.names = F)


