library(tidyverse)
library(data.table)

# Function to format a dataframe for the knee plot
format_knee_plot_df <- function(kneedle_path, topk_path) {
  
  # Read kneedle file path
  kneedle <- fread(kneedle_path)
  
  # Read topk file path and get the knee
  topk <- fread(topk_path)$`Absolute Weight` %>% min()
  
  # Get the model name
  model <- kneedle_path %>% strsplit("/") %>% unlist() %>% tail(1) %>% 
    strsplit("_") %>% unlist() %>% head(1)
  
  # Add information to kneedle dataframe
  kneedle %>%
    mutate(
      Knee = topk,
      Model = model
    )

}

# Make plotting data.frame
kneedle <- rbind(
  format_knee_plot_df("Comparison/kneedle_df/DIABLO_kneedle.csv", "Comparison/TopK/DIABLO_TopK.csv"),
  format_knee_plot_df("Comparison/kneedle_df/JACA_kneedle.csv", "Comparison/TopK/JACA_TopK.csv"),
  format_knee_plot_df("Comparison/kneedle_df/MOFA_kneedle.csv", "Comparison/TopK/MOFA_TopK.csv"),
  format_knee_plot_df("Comparison/kneedle_df/MultiMLP_kneedle.csv", "Comparison/TopK/MultiMLP_TopK.csv"),
  format_knee_plot_df("Comparison/kneedle_df/SLIDE_kneedle.csv", "Comparison/TopK/SLIDE_TopK.csv")
)

kneedle_plot <- ggplot(kneedle, aes(x = Rank, y = `Absolute Weight`, color = View)) +
  geom_point() +
  geom_hline(aes(yintercept = Knee), linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom",
        legend.text = element_text(size = 10)) +
  facet_wrap(.~Model, scales = "free")



