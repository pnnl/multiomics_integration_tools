library(data.table)
library(tidyverse)
library(patchwork)

theme_set(theme_bw())

# Let's visualize by large groups and by smaller groups
DIABLO <- fread("Comparison/u_matrix/DIABLO_umat_average.csv") %>%
  mutate(Time = c(rep("00wk", 3), rep("04wk", 3), rep("08wk", 3), rep("12wk", 3))) 
DIABLO$Factor2 <- NA
MOFA <- fread("Comparison/u_matrix/MOFA_umat.csv") %>% 
  rename(Group = group) %>%
  mutate(Time = c(rep("00wk", 3), rep("04wk", 3), rep("08wk", 3), rep("12wk", 3)))
SLIDE <- fread("Comparison/u_matrix/SLIDE_umat.csv") %>% 
  rename(Group = group) %>%
  mutate(Time = c(rep("00wk", 3), rep("04wk", 3), rep("08wk", 3), rep("12wk", 3)))
SLIDE$Factor2 <- NA

plotA <- DIABLO %>%
  ggplot(aes(x = Factor1, y = Factor2, color = Group)) +
  geom_point(size = 5) +
  xlab("Component 1 Average") + 
  ylab("Component 2 Average") +
  ggtitle("DIABLO") +
  scale_color_manual(values = c("orange", "orchid4")) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

plotB <- MOFA %>%
  ggplot(aes(x = Factor1, y = Factor2, color = Group)) +
  geom_point(size = 5) +
  xlab("Component 1") + 
  ylab("Component 2") +
  ggtitle("MOFA") +
  scale_color_manual(values = c("orange", "orchid4")) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

plotC <- SLIDE %>%
  ggplot(aes(x = Factor1, y = Factor2, color = Group)) +
  geom_point(size = 5) +
  xlab("Component 1") + 
  ylab("Component 2") +
  ggtitle("SLIDE") +
  scale_color_manual(values = c("orange", "orchid4")) +
  theme(plot.title = element_text(hjust = 0.5))

plotD <- DIABLO %>%
  ggplot(aes(x = Factor1, y = Factor2, color = Time)) +
  geom_point(size = 5) +
  xlab("Component 1 Average") + 
  ylab("Component 2 Average") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

plotE <- MOFA %>%
  ggplot(aes(x = Factor1, y = Factor2, color = Time)) +
  geom_point(size = 5) +
  xlab("Component 1") + 
  ylab("Component 2") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

plotF <- SLIDE %>%
  ggplot(aes(x = Factor1, y = Factor2, color = Time)) +
  geom_point(size = 5) +
  xlab("Component 1") + 
  ylab("Component 2") +
  theme(plot.title = element_text(hjust = 0.5))

factor_compare <- (plotA + plotB + plotC) / (plotD + plotE + plotF)

