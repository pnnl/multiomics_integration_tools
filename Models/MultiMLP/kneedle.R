library(tidyverse)
library(data.table)
library(kneedle)

############################
## HYPERAPARAMETER TUNING ##
############################

# Visualize hyperparameter tuning
deep_imv_tune <- fread("Models/MultiMLP/MultiMLP_tuning.csv") %>%
  mutate(Size = factor(Size, levels = c("0.25", "0.5", "1", "2", "4"))) %>%
  ggplot(aes(x = Size, y = Loss)) +
    geom_boxplot() +
    xlab("Grid size as an approximate proportion of input data size") +
    geom_jitter(height = 0, width = 0.5) +
    theme_bw()
ggsave("Plots/Supplemental_MultiMLPTuning.png", deep_imv_tune, dpi = 300)
deep_imv_tune

##############
## V MATRIX ##
##############

# Build a v matrix
fread("Models/MultiMLP/MultiMLP_ShapValues.csv") %>% 
  mutate(Factor = "Factor1") %>%
  rename(Weight = `Shapley Value`) %>%
  select(Factor, Weight, Feature, View) %>%
  fwrite("Comparison/v_matrix/MultiMLP_vmat.csv", quote = F, row.names = F)
  
##############################
## CALCULATE SHAPLEY VALUES ##
##############################

# Load shapley values, add rank
shaps <- fread("Models/MultiMLP/MultiMLP_ShapValues.csv") %>%
  mutate(`Absolute Weight` = abs(`Shapley Value`)) %>%
  arrange(-`Absolute Weight`) %>%
  mutate(Rank = 1:nrow(.)) %>%
  select(View, Feature, `Absolute Weight`, Rank)

fwrite(shaps, "Comparison/kneedle_df/MultiMLP_kneedle.csv", quote = F, row.names = F)

#############
## KNEEDLE ##
#############

# Calculate kneedle value
knee <- kneedle(shaps$Rank, shaps$`Absolute Weight`)

# Let's plot
ggplot(shaps, aes(x = Rank, y = `Absolute Weight`, color = View)) +
  geom_point() +
  theme_bw() +
  ggtitle("MultiMLP") + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  geom_hline(yintercept = knee[2], linetype = "dashed")

# Write output
shaps %>%
  filter(Rank <= knee[1]) %>%
  select(View, Feature, `Absolute Weight`, Rank) %>%
  fwrite("Comparison/TopK/MultiMLP_TopK.csv", quote = F, row.names = F)
