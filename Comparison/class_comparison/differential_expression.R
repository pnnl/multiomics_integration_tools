library(pmartR)
library(tidyverse)
library(data.table)
library(ggrepel)
library(patchwork)

theme_set(theme_bw())

fdata <- data.frame(
  Sample = c("109_00wk_A", "109_00wk_B", "109_00wk_C", 
             "109_04wk_A", "109_04wk_B", "109_04wk_C", 
             "109_08wk_A", "109_08wk_B", "109_08wk_C", 
             "109_12wk_A", "109_12wk_B", "109_12wk_C"),
  Time = c(rep("00wk", 3), rep("post-00wk", 9)),
  Time_edger = c(rep("Before", 3), rep("After", 9)) # EdgeR is picky
)

# 16S---------------------------------------------------------------------------

# Make seqData - no imputation was conducted
seqData <- as.seqData(
  e_data = fread("Dataset/Original/16S_Edata.csv"),
  edata_cname = "Feature",
  f_data = fdata, 
  fdata_cname = "Sample"
)
attr(seqData, "data_info")$data_scale <- "lcpm"

# Group designation
seqData <- group_designation(seqData, main_effects = "Time_edger")

# Run statistics
stats_seq <- diffexp_seq(omicsData = seqData, method = "edgeR")
fwrite(stats_seq, "Comparison/class_comparison/Microbiome_DiffExp.csv", quote = F)

stats_seq_plot <- stats_seq %>%
  mutate(Direction = ifelse(Flag_Before_vs_After == 0, "Not Significant", 
                     ifelse(Flag_Before_vs_After == -1, "Negative", "Positive")),
         Direction = factor(Direction, levels = c("Negative", "Not Significant", "Positive"))) %>%
ggplot(aes(x = Fold_change_Before_vs_After, y = -1 * log10(P_value_Before_vs_After))) +
  geom_label_repel(aes(label = Feature, color = Direction)) + 
  xlim(c(-0.65, 0.65)) +
  scale_color_manual(values = c("Positive" = "blue", "Not Significant" = "gray", "Negative" = "red")) +
  xlab("Fold change") +
  ylab("-Log10 P-Value") +
  ggtitle("16S") +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) 

# Metaproteomics----------------------------------------------------------------

# Load metaproteomics - use non-imputed data
proData <- as.proData(
  e_data = fread("Dataset/Original/Metaproteomics.csv") %>%
    filter(Feature %in% fread("Dataset/Scaled/Metaproteomics.csv")$Feature) %>%
    `colnames<-`(colnames(fread("Dataset/Scaled/Metaproteomics.csv"))),
  edata_cname = "Feature",
  f_data = fdata, 
  fdata_cname = "Sample",
  e_meta = fread("Dataset/EMeta/metaproteomics_emeta.txt") %>% 
    filter(Protein %in% fread("Dataset/Scaled/Metaproteomics.csv")$Feature) %>%
    rename(Feature = Protein) %>%
    mutate(EMeta_CName = 1:nrow(.)),
  emeta_cname = "EMeta_CName",
  data_scale = "log2",
  is_normalized = TRUE
)
proData <- group_designation(proData, main_effects = "Time")

# IMD-ANOVA filter
proData <- applyFilt(imdanova_filter(proData), proData, min_nonmiss_anova = 2)

# Conduct statistics
stats_pro <- imd_anova(proData, test_method = "anova", pval_adjust_a_fdr = "BH")
names_to_change <- grepl("_00wk_vs_post-00wk", colnames(stats_pro))
colnames(stats_pro)[names_to_change] <- gsub("_00wk_vs_post-00wk", "", colnames(stats_pro)[names_to_change])
fwrite(stats_pro, "Comparison/class_comparison/Metaproteomics_DiffExp.csv", quote = F)
stats_pro_plot <- stats_pro %>%
  mutate(Direction = ifelse(Flag_A == 0, "Not Significant", 
                     ifelse(Flag_A == -1, "Negative", "Positive")),
         Direction = factor(Direction, levels = c("Negative", "Not Significant", "Positive"))) %>%
  ggplot(aes(x = Fold_change, y = -1 * log10(P_value_A), color = Direction)) +
    geom_point() +  
    scale_color_manual(values = c("Positive" = "blue", "Not Significant" = "grey", "Negative" = "red")) +
    xlab("Fold change") +
    ylab("-Log10 P-Value") +
    ggtitle("Metaproteomics") +
    theme(
      legend.title = element_blank(), 
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    ) 

# Metabolomics Positive---------------------------------------------------------

# Load metabolomics positive - use non-imputed data
metabPos <- as.metabData(
  e_data = fread("Dataset/Original/Metabolomics_Positive.csv") %>%
    filter(Feature %in% fread("Dataset/Scaled/Metabolomics_Positive.csv")$Feature) %>%
    `colnames<-`(colnames(fread("Dataset/Scaled/Metabolomics_Positive.csv"))),
  edata_cname = "Feature",
  f_data = fdata, 
  fdata_cname = "Sample",
  e_meta = fread("Dataset/EMeta/metabolomics_positive.txt") %>% 
    filter(RT.MZ.ID %in% fread("Dataset/Scaled/Metabolomics_Positive.csv")$Feature) %>%
    rename(Feature = RT.MZ.ID) %>%
    mutate(EMeta_CName = 1:nrow(.)),
  emeta_cname = "EMeta_CName",
  data_scale = "log2",
  is_normalized = TRUE
)
metabPos <- group_designation(metabPos, main_effects = "Time")

# IMD-ANOVA filter
metabPos <- applyFilt(imdanova_filter(metabPos), metabPos, min_nonmiss_anova = 2)

# Conduct statistics
stats_metabp <- imd_anova(metabPos, test_method = "anova", pval_adjust_a_fdr = "BH")
names_to_change <- grepl("_00wk_vs_post-00wk", colnames(stats_metabp))
colnames(stats_metabp)[names_to_change] <- gsub("_00wk_vs_post-00wk", "", colnames(stats_metabp)[names_to_change])
fwrite(stats_metabp, "Comparison/class_comparison/Metabolomics_Positive_DiffExp.csv", quote = F)
stats_metabp_plot <- stats_metabp %>%
  mutate(Direction = ifelse(Flag_A == 0, "Not Significant", 
                     ifelse(Flag_A == -1, "Negative", "Positive")),
         Direction = factor(Direction, levels = c("Negative", "Not Significant", "Positive"))) %>%
  ggplot(aes(x = Fold_change, y = -1 * log10(P_value_A), color = Direction)) +
  geom_point() +  
  scale_color_manual(values = c("Positive" = "blue", "Not Significant" = "grey", "Negative" = "red")) +
  xlab("Fold change") +
  ylab("-Log10 P-Value") +
  ggtitle("Metabolomics Positive") +
  theme(
    legend.title = element_blank(), 
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)) 

# Metabolomics Negative---------------------------------------------------------

# Load metabolomics positive - use non-imputed data
metabNeg <- as.metabData(
  e_data = fread("Dataset/Original/Metabolomics_Negative.csv") %>%
    filter(Feature %in% fread("Dataset/Scaled/Metabolomics_Negative.csv")$Feature) %>%
    `colnames<-`(colnames(fread("Dataset/Scaled/Metabolomics_Negative.csv"))),
  edata_cname = "Feature",
  f_data = fdata, 
  fdata_cname = "Sample",
  e_meta = fread("Dataset/EMeta/metabolomics_negative.txt") %>% 
    filter(RT.MZ.ID %in% fread("Dataset/Scaled/Metabolomics_Negative.csv")$Feature) %>%
    rename(Feature = RT.MZ.ID) %>%
    mutate(EMeta_CName = 1:nrow(.)),
  emeta_cname = "EMeta_CName",
  data_scale = "log2",
  is_normalized = TRUE
)
metabNeg <- group_designation(metabNeg, main_effects = "Time")

# IMD-ANOVA filter
metabNeg <- applyFilt(imdanova_filter(metabNeg), metabNeg, min_nonmiss_anova = 2)

# Conduct statistics
stats_metabn <- imd_anova(metabNeg, test_method = "anova", pval_adjust_a_fdr = "BH")
names_to_change <- grepl("_00wk_vs_post-00wk", colnames(stats_metabn))
colnames(stats_metabn)[names_to_change] <- gsub("_00wk_vs_post-00wk", "", colnames(stats_metabn)[names_to_change])
fwrite(stats_metabn, "Comparison/class_comparison/Metabolomics_Negative_DiffExp.csv", quote = F)
stats_metabn_plot <- stats_metabn %>%
  mutate(Direction = ifelse(Flag_A == 0, "Not Significant", 
                     ifelse(Flag_A == -1, "Negative","Positive")),
         Direction = factor(Direction, levels = c("Negative", "Not Significant", "Positive"))) %>%
  ggplot(aes(x = Fold_change, y = -1 * log10(P_value_A), color = Direction)) +
  geom_point() +  
  scale_color_manual(values = c("Positive" = "blue", "Not Significant" = "grey", "Negative" = "red")) +
  xlab("Fold change") +
  ylab("-Log10 P-Value") +
  ggtitle("Metabolomics Negative") +
  theme(
    legend.title = element_blank(), 
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)) 

# Put information together------------------------------------------------------

# Make combined plot
diffplot <- stats_seq_plot + stats_pro_plot + stats_metabp_plot + stats_metabn_plot +
  plot_annotation(tag_levels = "A")
diffplot
ggsave("Plots/Supplemental_DiffPlot.png", diffplot, width = 8.3, height = 5.72, units = "in")

# Save topK information: View,Feature,AbsoluteWeight,Rank

rbind(
  stats_seq %>% 
    filter(P_value_Before_vs_After <= 0.05) %>% 
    mutate(View = "16S") %>%
    select(View,Feature),
  stats_pro %>% 
    filter(P_value_A <= 0.05) %>% 
    mutate(View = "metaproteomics") %>%
    select(View,Feature),
  stats_metabp %>% 
    filter(P_value_A <= 0.05) %>% 
    mutate(View = "metabolomics positive") %>%
    select(View,Feature),
  stats_metabn %>% 
    filter(P_value_A <= 0.05) %>% 
    mutate(View = "metabolomics negative") %>%
    select(View,Feature)
) %>%
  mutate(AbsoluteWeight=NA,Rank=NA) %>%
  fwrite("Comparison/TopK/DifferentialStats_TopK.csv", quote = F, row.names = F)

# Filter to emeta groups of interest--------------------------------------------

# Save features detected by at least one algorithm
detected <- c(
  fread("Comparison/TopK/DeepIMV_TopK.csv")$Feature,
  fread("Comparison/TopK/DIABLO_TopK.csv")$Feature,
  fread("Comparison/TopK/JACA_TopK.csv")$Feature,
  fread("Comparison/TopK/MOFA_TopK.csv")$Feature,
  fread("Comparison/TopK/SLIDE_TopK.csv")$Feature
) %>%
  unique()

# Make a general plotting schema for differential abundance
make_diff_plot <- function(x) {
  x %>%
    ggplot(aes(x = `Fold Change`, y = -1 * log10(`P Value`), color = Direction)) +
      geom_point() +
      scale_color_manual(values = c("Positive" = "blue", "Not Significant" = "grey", "Negative" = "red")) +
      xlab("Fold change") +
      ylab("-Log10 P-Value") +
      theme(legend.title = element_blank())
}



## Metabolite Classes ##

metab_classes <- c("Organic acids", "Peptides & Derivatives",
                   "Lipids & Derivatives", "Organoheterocyclics", 
                   "Benzenoids", "Carbohydrates & Sugars")

# Load all metabolite information and subset to relevant classes
metab_classes <- rbind(
  stats_metabn,
  stats_metabp
) %>%
  rename(`P Value` = P_value_A, `Fold Change` = Fold_change) %>%
  select(Feature, `P Value`, `Fold Change`) %>%
  filter(Feature %in% detected) %>%
  left_join(
    fread("Comparison/class_comparison/term_labeling_final.txt") %>%
      rename(Feature = RT.MZ.ID, Label = Final.Label) %>%
      select(Feature, Label)
  ) %>% 
  filter(Label %in% metab_classes) %>%
  mutate(
    Label = factor(Label, levels = metab_classes),
    Direction = ifelse(`P Value` > 0.05, "Not Significant",
                ifelse(`Fold Change` > 0, "Positive", "Negative"))
  ) %>%
  make_diff_plot() +
    theme(legend.position = "bottom") +
    facet_wrap(.~Label) 
metab_classes
ggsave("Plots/Supplemental_MetaboliteClassesDiffExp.png", metab_classes, dpi = 300)


## Protein Microbes ##

metaproteomics_emeta <- fread("Dataset/EMeta/metaproteomics_emeta.txt") %>%
  select(Protein, Strain) %>%
  rename(Feature = Protein)

proteomics_microbe_diffexp <- stats_pro %>%
  rename(`P Value` = P_value_A, `Fold Change` = Fold_change) %>%
  select(Feature, `P Value`, `Fold Change`) %>%
  filter(Feature %in% detected) %>%
  left_join(metaproteomics_emeta) %>%
  mutate(Strain = ifelse(Strain == "Sphingopyxisß", "Sphingopyxis", Strain),
         Direction = ifelse(`P Value` > 0.05, "Not Significant",
                            ifelse(`Fold Change` > 0, "Positive", "Negative"))
  ) %>%
  rename(Label = Strain) %>%
  make_diff_plot() +
    theme(legend.position = "bottom") +
    facet_wrap(.~Label) 
proteomics_microbe_diffexp
ggsave("Plots/Supplemental_MetaproteomicsSpeciesDiffPlot.png", proteomics_microbe_diffexp, dpi = 300)



