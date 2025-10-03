library(tidyverse)
library(data.table)
library(patchwork)

theme_set(theme_bw())

###############
## LOAD DATA ##
###############

# Load top K biomolecules with ranking information. Weights will not be used
# as their values have different meanings per algorithm.
topK <- rbind(
  fread("Comparison/TopK/DIABLO_TopK.csv") %>% mutate(Model = "DIABLO"),
  fread("Comparison/TopK/JACA_TopK.csv") %>% mutate(Model = "JACA"),
  fread("Comparison/TopK/MOFA_TopK.csv") %>% mutate(Model = "MOFA"),
  fread("Comparison/TopK/MultiMLP_TopK.csv") %>% mutate(Model = "MultiMLP"),
  fread("Comparison/TopK/SLIDE_TopK.csv") %>% mutate(Model = "SLIDE")
) %>%
  group_by(Model) %>%
  mutate(`Scaled Rank` = round(1 - ((Rank - min(Rank)) / max(Rank)), 2)) %>% 
  ungroup()

# Add common names
cn_subset <- topK %>%
  right_join(
    rbind(
      fread("Dataset/EMeta_TopK/TopK_Metabolites.txt") %>% 
        rename(Feature = RT.MZ.ID) %>% 
        select(Feature, Identifier),
      fread("Dataset/EMeta_TopK/TopK_Proteins.csv") %>%
        rename(Feature = Protein) %>%
        select(Feature, Identifier),
      data.frame(Feature = unique(topK$Feature[topK$View == "16S"]),
                 Identifier = unique(topK$Feature[topK$View == "16S"]))
    ), by = "Feature", relationship = "many-to-many"
  )

#########################
## Start building plot ##
#########################

# Calculate mean scaled ranks, keep instances of 2 known features across models
top_known <- cn_subset %>%
  group_by(View, Identifier, Model) %>%
  summarize(`Mean Scaled Rank` = mean(`Scaled Rank`)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(View, Identifier), values_from = `Mean Scaled Rank`, names_from = Model) %>%
  pivot_longer(3:7) %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>%
  rename(Model = name, `Scaled Rank` = value) %>%
  group_by(Identifier) %>%
  mutate(`Mean Scaled Rank` = mean(`Scaled Rank`)) %>%
  arrange(`Mean Scaled Rank`) %>%
  filter(`Scaled Rank` != 0) %>%
  mutate(Count = n()) %>%
  filter(Count != 1)
top_known$Identifier <- factor(top_known$Identifier, levels = unique(top_known$Identifier))


feature_rank_plot1 <- ggplot(top_known, aes(x = Model, y = Identifier, fill = `Mean Scaled Rank`)) +
  geom_tile() + 
  geom_text(aes(label = `Scaled Rank`)) + 
  scale_fill_gradient(low = "forestgreen", high = "green", na.value = "white") +
  xlab("") +
  ylab("") +
  theme_classic() +
  theme(legend.position = "none")

feature_rank_plot2 <- top_known %>%
  select(Identifier, `Mean Scaled Rank`) %>%
  rename(`Mean Rank` = `Mean Scaled Rank`) %>%
  unique() %>%
  mutate(Group = "Mean Rank") %>%
  ggplot(aes(x = Group, y = Identifier, fill = `Mean Rank`)) +
    geom_tile() +
    geom_text(aes(label = `Mean Rank`)) + 
    scale_fill_gradient(low = "forestgreen", high = "#1BD11A", na.value = "white") +
    xlab("") + 
    ylab("") +
    theme_classic() + 
    theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())

feature_plot <- feature_rank_plot1 + feature_rank_plot2 + plot_layout(widths = c(5, 1))


