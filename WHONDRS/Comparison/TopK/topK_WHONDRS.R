## This script compares the top k biomolecules from models. 
library(data.table)
library(tidyverse)
library(patchwork)

theme_set(theme_bw())

## Prepare Multiomics Counts----------------------------------------------------

# Load multiomics data
multiomics <- rbind(
  fread("Dataset/Scaled/Metabolomics.txt") %>% mutate(view = "metabolomics"),
  fread("Dataset/Scaled/Lipidomics_Positive.txt") %>% mutate(view = "lipidomics positive"),
  fread("Dataset/Scaled/Lipidomics_Negative.txt") %>% mutate(view = "lipidomics negative")
) %>%
  relocate(view) %>%
  pivot_longer(3:ncol(.)) %>%
  rename(sample = name, feature = Feature) %>%
  mutate(group = ifelse(grepl("N", sample), "None", "Surrounding")) %>%
  select(sample, group, feature, view, value)

# Get counts: 169, 34 - lipidomics negative, 81 - lipidomics positive, 54 - metabolomics
multiomics %>%
  select(feature, view) %>%
  unique() %>%
  select(view) %>%
  unlist() %>%
  table()

## Load and make plot-----------------------------------------------------------

diablo = fread("Comparison/TopK/DIABLO_TopK.csv") %>% 
  mutate(Model = "DIABLO") %>%
  mutate(Feature = c("DG 34:2", "Cer 39:1;4O", "Cer 34:1;3O", "Cer 34:0;3O", 
                     "PG 32:0", "Cer 36:0;4O", "Cer 42:1;3O", "PE 33:2",
                     "DGTS 38:2", "6-deoxy-D-glucose", "DGCC 36:2", "D-fructose",
                     "2-methyl-3-hydroxybutyric acid", "oxalacetic acid"))
  
  
jaca = fread("Comparison/TopK/JACA_TopK.csv") %>% mutate(Model = "JACA")
mofa = fread("Comparison/TopK/MOFA_TopK.csv") %>% mutate(Model = "MOFA") %>%
  mutate(Feature = gsub("_lipidomics positive", "", Feature, fixed = T))
multimlp = fread("Comparison/TopK/MultiMLP_TopK.csv") %>% mutate(Model = "MultiMLP")
slide = fread("Comparison/TopK/SLIDE_TopK.csv") %>% mutate(Model = "SLIDE") %>%
  mutate(Feature = gsub("_metab", "", Feature, fixed = T))


# Calculate an average rank
top_known = rbind(diablo, jaca, mofa, multimlp, slide) %>%
  group_by(Model) %>%
  mutate(`Scaled Rank` = round(1 - ((Rank - min(Rank)) / max(Rank)), 2)) %>% 
  ungroup() %>%
  select(Model, Feature, `Scaled Rank`) %>%
  pivot_wider(id_cols = Feature, values_from = `Scaled Rank`, names_from = Model) %>%
  pivot_longer(-1) %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>%
  rename(Model = name, `Scaled Rank` = value) %>%
  group_by(Feature) %>%
  mutate(`Mean Scaled Rank` = mean(`Scaled Rank`)) %>%
  arrange(`Mean Scaled Rank`) %>%
  filter(`Scaled Rank` != 0) %>%
  mutate(Count = n()) %>%
  filter(Count != 1)
top_known$Feature = factor(top_known$Feature, levels = unique(top_known$Feature))


feature_rank_plot1 <- ggplot(top_known, aes(x = Model, y = Feature, fill = `Mean Scaled Rank`)) +
  geom_tile() + 
  geom_text(aes(label = `Scaled Rank`)) + 
  scale_fill_gradient(low = "forestgreen", high = "green", na.value = "white") +
  xlab("") +
  ylab("") +
  theme_classic() +
  theme(legend.position = "none")

feature_rank_plot2 <- top_known %>%
  select(Feature, `Mean Scaled Rank`) %>%
  rename(`Mean Rank` = `Mean Scaled Rank`) %>%
  unique() %>%
  mutate(Group = "Mean Rank") %>%
  ggplot(aes(x = Group, y = Feature, fill = `Mean Rank`)) +
  geom_tile() +
  geom_text(aes(label = `Mean Rank`)) + 
  scale_fill_gradient(low = "forestgreen", high = "#1BD11A", na.value = "white") +
  xlab("") + 
  ylab("") +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())

feature_plot <- feature_rank_plot1 + feature_rank_plot2 + plot_layout(widths = c(5, 1))
feature_plot





