## This script compares the top k biomolecules from models. 
library(data.table)
library(tidyverse)
library(patchwork)
library(viridis)
library(ggdendro)

theme_set(theme_bw())

## Prepare Multiomics Counts----------------------------------------------------

# Load multiomics data
multiomics <- rbind(
  fread("Dataset/Scaled/16S_Edata.csv") %>% mutate(view = "16S"),
  fread("Dataset/Scaled/Metabolomics_Negative.csv") %>% mutate(view = "metabolomics negative"),
  fread("Dataset/Scaled/Metabolomics_Positive.csv") %>% mutate(view = "metabolomics positive"),
  fread("Dataset/Scaled/Metaproteomics.csv") %>% mutate(view = "metaproteomics")
) %>%
  relocate(view) %>%
  pivot_longer(3:ncol(.)) %>%
  rename(sample = name, feature = Feature) %>%
  mutate(group = map_chr(sample, function(x) {strsplit(x, "_") %>% unlist() %>% head(2) %>% tail(1)}),
         group = ifelse(group == "00wk", group, "post-00wk")) %>%
  select(sample, group, feature, view, value)

# Get counts: 7,286 features - 8 16S, 1752 metab neg, 2800 metab pos, 2726 metaP
multiomics %>%
  select(feature, view) %>%
  unique() %>%
  select(view) %>%
  unlist() %>%
  table()

## Load and make plot-----------------------------------------------------------

diablo <- fread("Comparison/TopK/DIABLO_TopK.csv") %>% mutate(Model = "DIABLO")
jaca <- fread("Comparison/TopK/JACA_TopK.csv") %>% mutate(Model = "JACA")
mofa <- fread("Comparison/TopK/MOFA_TopK.csv") %>% mutate(Model = "MOFA")
multimlp <- fread("Comparison/TopK/MultiMLP_TopK.csv") %>% mutate(Model = "MultiMLP")
slide <- fread("Comparison/TopK/SLIDE_TopK.csv") %>% mutate(Model = "SLIDE")
diffstats <- fread("Comparison/TopK/DifferentialStats_TopK.csv") %>% mutate(Model = "Diff Stats \nw/ PVal Adj")

# Make summary plot
topfeatureplot <- data.table(
  "Model" = c("Diff Stats \nw/ PVal Adj", "DIABLO", "JACA", "MOFA", "MultiMLP", "SLIDE"),
  "Top Feature Proportion" = c(nrow(diffstats) / length(unique(multiomics$feature)),
                               nrow(diablo) / length(unique(multiomics$feature)), 
                               nrow(jaca) / length(unique(multiomics$feature)),
                               nrow(mofa) / length(unique(multiomics$feature)),
                               nrow(multimlp) / length(unique(multiomics$feature)),
                               nrow(slide) / length(unique(multiomics$feature)))
) %>%
  mutate(
    Model = factor(Model, 
                   levels = c("Diff Stats \nw/ PVal Adj",
                              "DIABLO", "JACA", "MOFA", "MultiMLP", "SLIDE")),
    Proportion = round(`Top Feature Proportion`, 3)
  ) %>%
  ggplot(aes(x = Model, y = `Top Feature Proportion`)) +
    geom_bar(stat = "identity") + 
    ylim(c(0,1)) + 
    geom_text(aes(label = Proportion), vjust = -0.3) +
    theme_bw() +
    theme(legend.position = "none") +
    xlab("")
topfeatureplot

# Make proportion of each omics covered
omicscoverage_data <- rbind(
  diffstats %>% select(View, Model) %>% mutate(Model = "Diff Stats w/ PVal Adj"),
  diablo %>% select(View, Model),
  jaca %>% select(View, Model),
  mofa %>% select(View, Model),
  multimlp %>% select(View, Model),
  slide %>% select(View, Model)
) %>%
  group_by(View, Model) %>%
  summarize(Count = n()) %>%
  mutate(
    Proportion = ifelse(View == "16S", Count / 8,
                        ifelse(View == "metabolomics negative", Count / 1752,
                               ifelse(View == "metabolomics positive", Count / 2800, Count / 2726))),
    Proportion = round(Proportion, 3)
  ) %>%
  mutate(Model = factor(Model, levels = c("Diff Stats w/ PVal Adj", "DeepIMV", "DIABLO",
                                          "JACA", "MOFA", "MultiMLP", "SLIDE")))

omicscoverageplot <-  omicscoverage_data %>%
  ggplot(aes(x = View, y = Proportion, fill = View)) +
    geom_bar(stat = "identity") + 
    ylim(c(0,1.1)) + 
    geom_text(aes(label = Proportion), vjust = -0.3) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
    xlab("") +
    facet_wrap(.~Model)
omicscoverageplot

# Visualize overlap
pre_overlap <- rbind(
  diffstats %>% select(Feature, View, Model) %>% mutate(Model = "DiffStats"),
  diablo %>% select(Feature, View, Model),
  jaca %>% select(Feature, View, Model),
  mofa %>% select(Feature, View, Model),
  multimlp %>% select(Feature, View, Model),
  slide %>% select(Feature, View, Model)
) %>%
  group_by(Feature) %>%
  summarize(
    DiffStats = ifelse("DiffStats" %in% Model, 1, 0),
    DIABLO = ifelse("DIABLO" %in% Model, 1, 0),
    JACA = ifelse("JACA" %in% Model, 1, 0),
    MOFA = ifelse("MOFA" %in% Model, 1, 0),
    MultiMLP = ifelse("MultiMLP" %in% Model, 1, 0),
    SLIDE = ifelse("SLIDE" %in% Model, 1, 0),
  ) 

algorithms <- c("DiffStats", "DIABLO", "JACA", "MOFA", "MultiMLP", "SLIDE")

overlap <- expand.grid(algorithms, algorithms) %>%
  data.frame() %>%
  rename(Model1 = Var1, Model2 = Var2) %>%
  mutate(
    `First Count` = map2_dbl(Model1, Model2, function(x, y) {
      first <- pre_overlap$Feature[pre_overlap[[x]] == 1]
      second <- pre_overlap$Feature[pre_overlap[[y]] == 1]
      round(sum(first %in% second))
     }),
    `Second Count` = map_dbl(Model2, function(y) {
      length(pre_overlap$Feature[pre_overlap[[y]] == 1])
    }),
    Label = paste0(`First Count`, "/", `Second Count`),
    `Percent Overlap` = `First Count` / `Second Count` * 100,
    Model1 = ifelse(Model1 == "DiffStats", "Diff Stats w/ PVal Adj", as.character(Model1)),
    Model2 = ifelse(Model2 == "DiffStats", "Diff Stats w/ PVal Adj", as.character(Model2)),
    Model1 = factor(Model1, levels = c("Diff Stats w/ PVal Adj",  
                                       "DIABLO", "JACA", "MOFA", "MultiMLP", "SLIDE") %>% rev()),
    Model2 = factor(Model2, levels = c("Diff Stats w/ PVal Adj", 
                                       "DIABLO", "JACA", "MOFA", "MultiMLP", "SLIDE")),
  )

# Explore overlaps of integration models, focusing on just the largest amount per pair
best_overlap <- overlap %>% 
  filter(Model1 != Model2) %>% 
  filter(Model1 != "Diff Stats w/ PVal Adj" & Model2 != "Diff Stats w/ PVal Adj") %>% 
  mutate(ID = map2_chr(Model1, Model2, function(x,y) {paste(sort(c(x,y)), collapse = " ")})) %>%
  group_by(ID) %>%
  summarize(Best = max(`Percent Overlap`)) %>%
  arrange(-Best) %>%
  ungroup()
mean(best_overlap$Best)

# Build a hierarchical cluster
overlap %>% 
  filter(Model1 != "Diff Stats w/ PVal Adj" & Model2 != "Diff Stats w/ PVal Adj") %>%
  dplyr::select(Model1, Model2, `Percent Overlap`) %>%
  pivot_wider(id_cols = Model1, names_from = Model2, values_from = `Percent Overlap`) %>%
  dplyr::select(-Model1) %>%
  dist() %>%
  hclust() %>%
  ggdendrogram(rotate = TRUE) 

# Calculate a percent overlap with each other
overlapplot <- ggplot(overlap, aes(x = Model2, y = Model1, fill = `Percent Overlap`)) +
  geom_tile() +
  geom_text(aes(label = Label), color = "white") +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("") +
  ylab("")
overlapplot

# Make a topk dataframe
topK <- rbind(diablo, jaca, mofa, multimlp, slide)

# Get counts of the features 
counts <- topK$Feature %>% 
  table(dnn = "Feature") %>%
  data.frame()

# Explore the 2 shared metabolites
fread("Dataset/EMeta/metabolomics_negative.txt") %>%
  filter(RT.MZ.ID %in% unlist(counts[counts$Freq == 5, "Feature"]))


# Plot the count of counts
countplot <- counts$Freq %>%
  table(dnn = "Number") %>%
  data.frame() %>%
  mutate(Number = as.character(Number)) %>%
  mutate(Freq = as.numeric(Freq)) %>%
  ggplot(aes(x = Number, y = Freq)) +
  geom_bar(stat = "identity", color = "black", fill = "steelblue") +
  geom_text(aes(label = Freq), vjust = -0.2) + 
  ylim(c(0, 220)) +
  xlab("Common feature occurrence across models") +
  ylab("Frequency")
countplot

topk_figure <- (topfeatureplot + omicscoverageplot + plot_layout(widths = c(1,2))) / 
  (countplot + overlapplot + plot_layout(widths = c(1,2))) + 
  plot_layout(height = c(1, 1.5)) +
  plot_annotation(tag_levels = "A") 
topk_figure
  
# Two non-significant selections by MultiMLP

features <- pre_overlap %>%
  filter(MultiMLP == 1 & DiffStats == 0)

rbind(
  fread("Comparison/class_comparison/Metabolomics_Negative_DiffExp.csv"),
  fread("Comparison/class_comparison/Metabolomics_Positive_DiffExp.csv")
) %>%
  filter(Feature %in% features$Feature)

rbind(
  fread("Dataset/EMeta/metabolomics_negative.txt"),
  fread("Dataset/EMeta/metabolomics_positive.txt")
) %>%
  filter(RT.MZ.ID %in% features$Feature)

# Extract weights
diablo_v <- fread("Comparison/v_matrix/DIABLO_vmat.csv") %>% dplyr::select(Feature, Weight) %>% 
  mutate(Weight = abs(Weight)) %>% rename(DIABLO = Weight) 
jaca_v <- fread("Comparison/v_matrix/JACA_vmat.csv") %>% dplyr::select(Feature, Weight) %>% 
  mutate(Weight = abs(Weight)) %>% rename(JACA = Weight)
mofa_v <- fread("Comparison/v_matrix/MOFA_vmat.csv") %>% dplyr::filter(Factor == "Factor1") %>% 
  dplyr::select(Feature, Weight) %>% mutate(Weight = abs(Weight)) %>% rename(MOFA = Weight)
multimlp_v <- fread("Comparison/v_matrix/MultiMLP_vmat.csv") %>% dplyr::select(Feature, Weight) %>% 
  mutate(Weight = abs(Weight)) %>% rename(MultiMLP = Weight)
slide_v <- fread("Comparison/v_matrix/SLIDE_vmat.csv") %>% dplyr::select(Feature, Weight) %>% 
  mutate(Weight = abs(Weight)) %>% rename(SLIDE = Weight)

left_join(diablo_v, jaca_v) %>%
  left_join(mofa_v) %>%
  left_join(multimlp_v) %>%
  left_join(slide_v) %>%
  dplyr::select(-Feature) %>%
  cor(method = "spearman")  %>%
  dist() %>%
  hclust() %>%
  ggdendrogram(rotate = TRUE) +
  labs(xlab = "Height")
  






