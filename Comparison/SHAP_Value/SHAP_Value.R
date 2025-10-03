library(data.table)
library(tidyverse)
library(patchwork)

# Load data
shaps <- rbind(
  fread("Comparison/SHAP_Value/MOFA_SHAP_Value.csv") %>% mutate(Model = "MOFA")
) %>%
  mutate(`Mean Absolute SHAP Value` = round(`Mean Absolute SHAP Value`, 4))
  
# Make function to visualize data by model
viz_shap <- function(theModel) {
  shaps %>%
    filter(Model == theModel) %>%
    arrange(`Mean Absolute SHAP Value`) %>%
    mutate(Factor = factor(Factor, levels = Factor)) %>%
    ggplot(aes(x = Factor, y = `Mean Absolute SHAP Value`)) +
      geom_bar(stat = "identity", color = "black", fill = "steelblue") +
      geom_text(aes(label = `Mean Absolute SHAP Value`), hjust = -0.1) +
      coord_flip() +
      theme_bw() +
      ggtitle(theModel) +
      ylim(c(0,0.25)) + 
      theme(plot.title = element_text(hjust = 0.5))
}

MOFA_shap <- viz_shap("MOFA") 
ggsave("Plots/Supplemental_MOFAShap.png", MOFA_shap, dpi = 300)
MOFA_shap

