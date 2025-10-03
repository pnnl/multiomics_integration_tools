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
) %>% select(-`Absolute Weight`)

# Add emeta information
microbiome <- fread("Dataset/Scaled/16S_Edata.csv")$Feature %>% unlist()
metab_neg_emeta <- fread("Dataset/EMeta/metabolomics_negative.txt")
metab_pos_emeta <- fread("Dataset/EMeta/metabolomics_positive.txt")
metaproteomics_emeta <- fread("Dataset/EMeta/metaproteomics_emeta.txt")

########################
## 16S REPRESENTATION ##
########################

# Which microbes are represented? Use pivoting to create NAs for ease of visualization
microbe_rank <- topK %>%
  group_by(Model) %>%
  mutate(`Scaled Rank` = round(1 - ((Rank - min(Rank)) / max(Rank)), 2)) %>%
  filter(View == "16S") %>%
  select(Feature, `Scaled Rank`, Model) %>%
  pivot_wider(values_from = `Scaled Rank`, names_from = Model, id_cols = Feature) %>%
  pivot_longer(2:ncol(.)) %>%
  rename(Model = name, `Scaled Rank` = value) 

# Add an average rank weight for all models
microbe_rank_order <- microbe_rank %>%
  mutate(`Scaled Rank` = ifelse(is.na(`Scaled Rank`), 0, `Scaled Rank`)) %>%
  group_by(Feature) %>%
  summarize(`Mean Scaled Rank` = mean(`Scaled Rank`)) %>%
  arrange(`Mean Scaled Rank`)

# Make the microbe plot 
microbe_rank_plot1 <- microbe_rank %>%
  mutate(Feature = factor(Feature, levels = microbe_rank_order$Feature)) %>%
  ggplot(aes(x = Model, y = Feature, fill = `Scaled Rank`)) +
    geom_tile() + 
    geom_text(aes(label = `Scaled Rank`)) + 
    scale_fill_gradient(low = "forestgreen", high = "green", na.value = "white") +
    xlab("") +
    ylab("") +
    theme_classic() + 
    theme(legend.position = "none")

# Make a plot for the average of averages
microbe_rank_plot2 <- microbe_rank_order %>%
  mutate(`Mean Rank` = "Mean Rank",
         `Mean Scaled Rank` = round(`Mean Scaled Rank`, 2),
         Feature = factor(Feature, levels = microbe_rank_order$Feature)) %>%
  ggplot(aes(x = `Mean Rank`, y = Feature, fill = `Mean Scaled Rank`)) +
    geom_tile() +
    geom_text(aes(label = `Mean Scaled Rank`)) + 
    scale_fill_gradient(low = "forestgreen", high = "#1BD11A", na.value = "white") +
    xlab("") + 
    ylab("") +
   theme_classic() + 
    theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())

microbeplot <- microbe_rank_plot1 + microbe_rank_plot2 + plot_layout(widths = c(5, 1))
microbeplot

#############################################
## GENERAL METABOLITE CLASS REPRESENTATION ##
#############################################

# Bind together as there are no duplicate RT.MZ.IDs, filter to only topK features,
# and pull out all general metabolite classes
#rbind(
#  metab_pos_emeta,
#  metab_neg_emeta
#) %>%
#  select(RT.MZ.ID, Compound.ID, Library.ID.Match, MSCC.Compound.Match, 
#         Compound.Class.CANOPUS.Main, 
#         Compound.Class.CANOPUS.UpperClass.1,
#         Compound.Class.CANOPUS.UpperClass.2, Feature.ID.Final.DF) %>%
#  filter(RT.MZ.ID %in% unique(topK$Feature)) %>%
#  fwrite("Comparison/class_comparison/term_labeling.txt", sep = "\t", quote = F, row.names = F)
## CANOPUS was prioritized, followed by MSCC

# Load metabolite types - These were manually collapsed for consistency
general_metabs <- fread("Comparison/class_comparison/term_labeling_final.txt") %>%
  rename(Feature = RT.MZ.ID, Group = Final.Label) %>%
  select(Feature, Group)

# Get mean ranks per group
general_metabs_means <- topK %>%
  filter(Feature %in% general_metabs$Feature) %>%
  mutate(`Scaled Rank` = round(1 - ((Rank - min(Rank)) / max(Rank)), 2),
         `Scaled Rank` = ifelse(`Scaled Rank` == 0, 0.01, `Scaled Rank`)) %>%
  left_join(general_metabs) %>%
  select(Model, `Scaled Rank`, Group) %>%
  group_by(Model, Group) %>%
  summarize(`Mean Scaled Rank` = mean(`Scaled Rank`)) %>%
  pivot_wider(id_cols = Group, names_from = Model, values_from = `Mean Scaled Rank`) %>%
  pivot_longer(2:ncol(.)) %>%
  rename(Model = name, `Mean Scaled Rank` = value) %>%
  mutate(`Mean Scaled Rank` = ifelse(is.na(`Mean Scaled Rank`), 0, `Mean Scaled Rank`))
  
# Get counts 
general_metabs_counts <- general_metabs %>%
  count(Group) %>%
  rename(Count = n) %>%
  left_join(general_metabs_means) %>%
  group_by(Group) %>%
  summarize(
    `Mean Rank` = mean(`Mean Scaled Rank`),
    Count = head(Count, 1)
  ) %>%
  arrange(`Mean Rank`)

# Add group information and take the mean rank per group per model.
# Use pivoting to fill NA values again
metabolite_rank1 <- general_metabs_means %>%
  mutate(
    Group = factor(Group, levels = general_metabs_counts$Group),
    `Mean Scaled Rank` = ifelse(`Mean Scaled Rank` == 0, NA, `Mean Scaled Rank`),
    `Mean Scaled Rank` = round(`Mean Scaled Rank`, 2)
  ) %>%
  ggplot(aes(x = Model, y = Group, fill = `Mean Scaled Rank`)) +
    geom_tile() + 
    geom_text(aes(label = `Mean Scaled Rank`)) + 
    scale_fill_gradient(low = "forestgreen", high = "green", na.value = "white") +
    xlab("") +
    ylab("") +
    theme_classic() +
    theme(legend.position = "none")
metabolite_rank1

# Add mean rank information
metabolite_rank2 <- general_metabs_counts %>%
  drop_na() %>%
  select(Group, `Mean Rank`) %>%
  rename(Rank = `Mean Rank`) %>%
  arrange(Rank) %>%
  mutate(
    Group = factor(Group, levels = general_metabs_counts$Group),
    `Mean Rank` = "Mean Rank"
  ) %>% 
  select(Group, Rank, `Mean Rank`) %>%
  mutate(Rank = round(Rank, 2)) %>%
  ggplot(aes(x = `Mean Rank`, y = Group, fill = Rank)) +
    geom_tile() +
    geom_text(aes(label = Rank)) + 
    scale_fill_gradient(low = "forestgreen", high = "#1BD11A", na.value = "white") +
    xlab("") + 
    ylab("") +
    theme_classic() + 
    theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())
metabolite_rank2

# Add count information
metabolite_rank3 <- general_metabs_counts %>%
  drop_na() %>%
  rename(Freq = Count) %>%
  mutate(
    Group = factor(Group, levels = general_metabs_counts$Group),
    Count = "Count"
  ) %>%
  select(Group, Count, Freq) %>%
  ggplot(aes(x = Count, y = Group, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq)) + 
    scale_fill_gradient(low = "white", high = "white", na.value = "white") +
    xlab("") + 
    ylab("") +
    theme_classic() + 
    theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())

metaboliteplot <- metabolite_rank1 + metabolite_rank2 + metabolite_rank3 + plot_layout(widths = c(4, 1, 1))
metaboliteplot

#####################################
## PROTEIN ACTIVITY REPRESENTATION ##
#####################################

# Make a list to hold the main properties of COG classes
cog_list <- list(
  "C" = "Energy Production/Conversion",
  "D" = "Cell Cycle",
  "E" = "Amino Acid Metabolism/Transport",
  "F" = "Nucleotide Metabolism/Transport",
  "G" = "Carbohydrate Metabolism/Transport",
  "H" = "Coenzyme Metabolism/Transport",
  "I" = "Lipid Metabolism/Transport",
  "J" = "Translation/Ribosome Structure",
  "K" = "Transcription",
  "L" = "Replication/Recombination/Repair",
  "M" = "Cell Wall/Membranes", 
  "N" = "Cell Motility",
  "O" = "PTMs/Chaperones",
  "P" = "Inorganic Ion Metabolism/Transport",
  "Q" = "Metabolite Metabolism/Transport",
  "T" = "Signal Transduction",
  "U" = "Vesicular Transport/Secretion",
  "W" = "Extracellular Structures"
)

# Calculate mean scaled rank per property. Add missing models and use pivoting to add NAs needed.
protein_info <- topK %>%
  group_by(Model) %>%
  mutate(`Scaled Rank` = round(1 - ((Rank - min(Rank)) / max(Rank)), 2),
         `Scaled Rank` = ifelse(`Scaled Rank` == 0, 0.01, `Scaled Rank`)) %>%
  ungroup() %>%
  filter(View == "metaproteomics") %>%
  left_join(metaproteomics_emeta %>% rename(Feature = Protein) %>% select(Feature, COG_category) %>% unique()) %>%
  filter(COG_category %in% c("-", "S", "") == FALSE) %>% # Remove unknown functions
  mutate(
    COG_category = map(COG_category, function(x) {strsplit(x, "") %>% unlist()})
  ) %>%
  unnest(cols = c(COG_category)) %>%
  mutate(
    Property = map_chr(COG_category, function(x) {unlist(cog_list[x])})
  )

protein_count <- protein_info %>%
  count(COG_category) %>%
  mutate(Property = lapply(COG_category, function(x) {cog_list[[x]]}) %>% unlist()) %>%
  rename(Count = n) %>%
  select(Property, Count)

# Save the scaled ranks per model and property
protein_props <- protein_info %>%
  select(Model, `Scaled Rank`, Property) %>%
  rbind( # Add missing models
    c("SLIDE", 0, "Cell Motility")
  ) %>%
  mutate(`Scaled Rank` = as.numeric(`Scaled Rank`)) %>%
  group_by(Model, Property) %>%
  summarize(
    `Mean Scaled Rank` = mean(`Scaled Rank`),
  ) %>%
  ungroup() %>%
  pivot_wider(names_from = Model, values_from = `Mean Scaled Rank`, id_cols = Property) %>%
  pivot_longer(2:ncol(.)) %>%
  rename(Model = name, `Mean Scaled Rank` = value) %>%
  mutate(`Mean Scaled Rank` = ifelse(is.na(`Mean Scaled Rank`), 0, `Mean Scaled Rank`))

# Build order
protein_order <- protein_props %>%
  mutate(`Mean Scaled Rank` = ifelse(is.na(`Mean Scaled Rank`), 0, `Mean Scaled Rank`)) %>%
  group_by(Property) %>%
  summarize(
    `Mean Scaled Rank` = mean(`Mean Scaled Rank`)
  ) %>%
  arrange(`Mean Scaled Rank`) %>%
  left_join(protein_count, by = "Property")
  

protein_rank_plot1 <- protein_props %>%
  mutate(
    Property = factor(Property, levels = protein_order$Property),
    `Mean Scaled Rank` = round(`Mean Scaled Rank`, 2),
    `Mean Scaled Rank` = ifelse(`Mean Scaled Rank` == 0, NA, `Mean Scaled Rank`),
  ) %>%
  ggplot(aes(x = Model, y = Property, fill = `Mean Scaled Rank`)) +
    geom_tile() +
    geom_text(aes(label = `Mean Scaled Rank`)) + 
    scale_fill_gradient(low = "forestgreen", high = "green", na.value = "white") +
    xlab("") + 
    ylab("") +
    theme_classic() + 
    theme(legend.position = "none")
protein_rank_plot1
  
protein_rank_plot2 <- protein_order %>%
  rename(Rank = `Mean Scaled Rank`) %>%
  mutate(
    Property = factor(Property, levels = protein_order$Property),
    `Mean Rank` = "Mean Rank",
    Rank = round(Rank, 2),
    Rank = ifelse(Rank == 0, 0.01, Rank)
  ) %>%
  ggplot(aes(x = `Mean Rank`, y = Property, fill = Rank)) +
  geom_tile() +
  geom_text(aes(label = Rank)) + 
  scale_fill_gradient(low = "forestgreen", high = "#21BA1D", na.value = "white") +
  xlab("") + 
  ylab("") +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())
protein_rank_plot2

# Add count information
protein_rank_plot3 <- protein_count %>%
  rename(Freq = Count) %>%
  mutate(
    Property = factor(Property, levels = protein_order$Property),
    Count = "Count"
  ) %>%
  group_by(Property, Count) %>%
  mutate(Freq = max(Freq, na.rm = T)) %>%
  select(Property, Count, Freq) %>%
  unique() %>%
  ggplot(aes(x = Count, y = Property, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq)) + 
    scale_fill_gradient(low = "white", high = "white", na.value = "white") +
    xlab("") + 
    ylab("") +
    theme_classic() + 
    theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())

proteinplot <- protein_rank_plot1 + protein_rank_plot2 + protein_rank_plot3 + plot_layout(widths = c(4, 1, 1))
proteinplot

#####################
## MAKE FINAL PLOT ##
#####################

Figure4 <- (metaboliteplot) / proteinplot + 
  plot_annotation(tag_levels = list(c("A", "", "", "B", "", "")))
Figure4

##############################
## CLASSES PER TOP PROTEINS ##
##############################

# Get species information
prot_species <- topK %>%
  group_by(Model) %>%
  mutate(`Scaled Rank` = round(1 - ((Rank - min(Rank)) / max(Rank)), 2),
         `Scaled Rank` = ifelse(`Scaled Rank` == 0, 0.01, `Scaled Rank`)) %>%
  ungroup() %>%
  filter(View == "metaproteomics") %>%
  left_join(
    metaproteomics_emeta %>% select(Protein, Strain) %>% rename(Feature = Protein), 
    by = "Feature",
    relationship = "many-to-many"
  ) %>%
  filter(Strain != "") %>%
  mutate(Strain = ifelse(Strain == "Sphingopyxisß", "Shingopyxis", Strain)) %>%
  group_by(Model, Strain) %>%
  summarize(
    `Mean Scaled Rank` = mean(`Scaled Rank`),
  ) %>%
  ungroup() %>%
  rbind(
    c("SLIDE", "Dyadobacter", 0)
  ) %>%
  pivot_wider(names_from = Model, values_from = `Mean Scaled Rank`, id_cols = Strain) %>%
  pivot_longer(2:ncol(.)) %>%
  rename(Model = name, `Mean Scaled Rank` = value) %>%
  mutate(`Mean Scaled Rank` = as.numeric(`Mean Scaled Rank`))

prot_species_count <- metaproteomics_emeta %>%
  filter(Protein %in% topK$Feature) %>%
  count(Strain, name = "Frequency") %>%
  mutate(Strain = ifelse(Strain == "Sphingopyxisß", "Shingopyxis", Strain)) 

prot_species_order <- prot_species %>%
  group_by(Strain) %>%
  mutate(`Mean Scaled Rank` = ifelse(is.na(`Mean Scaled Rank`), 0, `Mean Scaled Rank`)) %>%
  summarize(`Average Rank` = mean(`Mean Scaled Rank`)) %>%
  left_join(prot_species_count) %>%
  arrange(`Average Rank`)
  
# Plot 1: All ranks
prot_species_plot1 <- prot_species %>%
  mutate(
    Strain = factor(Strain, levels = prot_species_order$Strain),
    `Mean Scaled Rank` = round(`Mean Scaled Rank`, 2),
    `Mean Scaled Rank` = ifelse(`Mean Scaled Rank` == 0, NA, `Mean Scaled Rank`),
  ) %>%
  ggplot(aes(x = Model, y = Strain, fill = `Mean Scaled Rank`)) +
    geom_tile() +
    geom_text(aes(label = `Mean Scaled Rank`)) + 
    scale_fill_gradient(low = "forestgreen", high = "green", na.value = "white") +
    xlab("") + 
    ylab("") +
    theme_classic() + 
    theme(legend.position = "none")
prot_species_plot1

# Plot 2: Average rank
prot_species_plot2 <- prot_species_order %>%
  rename(Rank = `Average Rank`) %>%
  mutate(
    Strain = factor(Strain, levels = prot_species_order$Strain),
    `Mean Rank` = "Mean Rank",
    Rank = round(Rank, 2),
    Rank = ifelse(Rank == 0, 0.01, Rank)
  ) %>%
  ggplot(aes(x = `Mean Rank`, y = Strain, fill = Rank)) +
    geom_tile() +
    geom_text(aes(label = Rank)) + 
    scale_fill_gradient(low = "forestgreen", high = "#21BA1D", na.value = "white") +
    xlab("") + 
    ylab("") +
    theme_classic() + 
    theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())
prot_species_plot2

# Plot3: Protein Count
prot_species_plot3 <- prot_species_order %>%
  rename(Freq = Frequency) %>%
  mutate(
    Strain = factor(Strain, levels = prot_species_order$Strain),
    Count = "Count"
  ) %>%
  group_by(Strain, Count) %>%
  mutate(Freq = sum(Freq, na.rm = T)) %>%
  select(Strain, Count, Freq) %>%
  unique() %>%
  ggplot(aes(x = Count, y = Strain, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq)) + 
    scale_fill_gradient(low = "white", high = "white", na.value = "white") +
    xlab("") + 
    ylab("") +
    theme_classic() + 
    theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())
prot_species_plot3

proteinspeciesplot <- prot_species_plot1 + prot_species_plot2 + prot_species_plot3 + plot_layout(widths = c(4, 1, 1))
proteinspeciesplot




