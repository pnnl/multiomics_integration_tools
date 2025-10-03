library(tidyverse)
library(data.table)

# Define group of top ranked
topK <- c("Sphingopyxis", 
          "Streptomyces", 
          "Variovorax", 
          "Dyadobacter", 
          "Neorhizobium",
          "G1_CDS.5066",                # xylF - Streptomyces
          "G7_CDS.3416",                # lpdA - Dyadobacter
          "G12_CDS.1173",               # exaA - Variovorax
          "G1_CDS.2166", "G1_CDS.3045", # gapA - Streptomyces
          "G1_CDS.4718",                # cobN - Streptomyces
          "G1_CDS.6022",                # sucB - Streptomyces
          "G5_CDS.2369", "G5_CDS.2837", # cspA - Neorhizobium
          "1.03_116.034809",            # maleamate
          "1.47_245.100185",            # 3-isobutyl-1-methyxanthine
          "4.4_287.171867",             # ammothamnine
          "1.05_119.034024",            # glycoaldehyde
          "5.78_319.081943",            # demethylsulochrin
          "1.03_103.039633",            # acetoacetate
          "1.01_146.092771",            # 4-guanidinobutanoate
          "1.9_295.115861",             # N-glycosyl-L-asparagine
          "1.04_118.050454",            # N-acetylglycine
          "1.05_98.024382"              # Maleimide
          )

# Create sample group information
fdata <- data.frame(
  name = c("109_00wk_A", "109_00wk_B", "109_00wk_C", 
           "109_04wk_A", "109_04wk_B", "109_04wk_C", 
           "109_08wk_A", "109_08wk_B", "109_08wk_C", 
           "109_12wk_A", "109_12wk_B", "109_12wk_C"),
  Time = c(rep("00wk", 3), rep("post-00wk", 9))
)

# Calculate differences in means per group
mean_diff <- rbind(
  fread("Dataset/Scaled/16S_Edata.csv") %>% pivot_longer(-1),
  fread("Dataset/Scaled/Metabolomics_Negative.csv") %>% pivot_longer(-1),
  fread("Dataset/Scaled/Metabolomics_Positive.csv") %>% pivot_longer(-1),
  fread("Dataset/Scaled/Metaproteomics.csv") %>% pivot_longer(-1)
) %>%
  left_join(fdata) %>%
  group_by(Feature, Time) %>%
  summarize(Mean = mean(value)) %>%
  ungroup() %>%
  pivot_wider(id_cols = Feature, names_from = Time, values_from = Mean) %>%
  mutate(Difference = `post-00wk` - `00wk`,
         AbsDifference = abs(Difference)) %>%
  arrange(-AbsDifference) %>%
  mutate(TopK = Feature %in% topK)

# Now per week
week <- rbind(
  fread("Dataset/Scaled/16S_Edata.csv") %>% pivot_longer(-1),
  fread("Dataset/Scaled/Metabolomics_Negative.csv") %>% pivot_longer(-1),
  fread("Dataset/Scaled/Metabolomics_Positive.csv") %>% pivot_longer(-1),
  fread("Dataset/Scaled/Metaproteomics.csv") %>% pivot_longer(-1)
) %>%
  filter(Feature %in% topK) %>%
  mutate(Week = map_chr(name, function(x) {strsplit(x, "_") %>% unlist() %>% head(2) %>% tail(1)})) %>%
  group_by(Feature, Week) %>%
  summarize(Mean = mean(value)) %>%
  ungroup() %>%
  pivot_wider(id_cols = Feature, names_from = Week, values_from = Mean)
  
# Load weights and look at key players 
keyplayers <- fread("Dataset/chitin_metabolism_keyplayers.csv")

fread("Comparison/TopK/MOFA_TopK.csv") %>% select(`Absolute Weight`) %>% unlist() %>% min()
fread("Comparison/v_matrix/MOFA_vmat.csv") %>% filter(Factor == "Factor1") %>% 
  mutate(Weight = abs(Weight)) %>% arrange(-Weight) %>% right_join(keyplayers %>% rename(Feature = Term)) 

fread("Comparison/TopK/MultiMLP_TopK.csv")  %>% select(`Absolute Weight`) %>% unlist() %>% min()
fread("Comparison/v_matrix/MultiMLP_vmat.csv") %>% 
  mutate(Weight = abs(Weight)) %>% arrange(-Weight) %>% right_join(keyplayers %>% rename(Feature = Term)) 

fread("Comparison/TopK/SLIDE_TopK.csv") %>% select(`Absolute Weight`) %>% unlist() %>% min()
fread("Comparison/v_matrix/SLIDE_vmat.csv") %>% 
  mutate(Weight = abs(Weight)) %>% arrange(-Weight) %>% right_join(keyplayers %>% rename(Feature = Term)) 




  
