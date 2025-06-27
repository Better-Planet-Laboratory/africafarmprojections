# ==============================================================================
# FARM SIZE DATA PROCESSING SCRIPT
# ==============================================================================
# Description: Takes farm size data, checks against FAO, fixes issues and saves output
# ==============================================================================

# Load Required Libraries =====================================================
library(sf)
library(terra)
library(mapview)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)

# Load Data ===================================================================
# Read farm size data
size <- read_sf("input/handoff/SmallholderMap_20201202/SmallholderMap_20201202.shp")

# Read metadata (note: Lowder 2012 not 2016)
meta <- read_xlsx("input/handoff/SmallholderMap_20201202/SmallholderMapping_Metadata.xlsx")

# Data Quality Check ==========================================================
# Check for missing values across farm size categories
size_l <- size %>% 
  pivot_longer(
    cols = A0_1:PN0_20,
    names_to = "name",
    values_to = "value"
  ) %>%
  as.data.frame() %>%
  group_by(NAME_0, name) %>% 
  summarise(n = ifelse(sum(!is.na(value)) > 0, 1, NA), .groups = "drop")

# Create visualization of missing data
miss <- ggplot(size_l, aes(y = NAME_0, x = name, fill = n)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "dark blue", na.value = "gray") +
  theme_minimal() +
  ggtitle("Missing Data Visualization")

# Data Exploration ============================================================
# Note: Some countries have max farm size so proportions weren't calculated
# e.g., DRC listed on FAO as range 2-<5, not >5
# https://www.fao.org/faostat/en/#data/WCAD

# Create maps to visualize data
mapview(size, layer.name = "N0_1", zcol = "N0_1")  # Most sensible but missing data
mapview(size, layer.name = "PN0_2", zcol = "PN0_2")  # Most sensible but missing data

# Data Selection ==============================================================
# Select columns of interest
size <- size %>%
  select(GID_0:NAME_1, N0_1:N20_)

# Load Additional Data ========================================================
# Read farm structure data
#https://www.nature.com/articles/s41893-023-01110-y
#https://github.com/Better-Planet-Laboratory/farm-decline
fn <- readRDS("input/ssp2t.rds")[[4]]

# Load FAO World Census of Agriculture data
#https://www.fao.org/world-census-agriculture/en
FAO <- read.csv("input/WCA_2010.csv")
FAO$Area[FAO$Area == "United Republic of Tanzania"] <- "Tanzania"

# Data Harmonization ==========================================================
# Define farm size categories
item_levels <- c("N0_1", "N1_2", "N2_5", "N5_10", "N10_20", "N20_")

# Harmonize FAO categories to match our classification
FAO <- FAO %>%
  mutate(Item = recode(Item,
                       `Holdings without land` = "N0_1",
                       `Holdings with land size 0-<1` = "N0_1",
                       `Holdings with land size 1-<2` = "N1_2",
                       `Holdings with land size 2-<3` = "N2_5",
                       `Holdings with land size 3-<4` = "N2_5",
                       `Holdings with land size 4-<5` = "N2_5",
                       `Holdings with land size 5-<10` = "N5_10",
                       `Holdings with land size 10-<20` = "N10_20",
                       `Holdings with land size 20-<50` = "N20_",
                       `Holdings with land size 50-<100` = "N20_",
                       `Holdings with land size 100-<200` = "N20_",
                       `Holdings with land size 200-<500` = "N20_",
                       `Holdings with land size 500-<1000` = "N20_",
                       `Holdings with land size 1000-<2500` = "N20_"
  )) %>%
  mutate(Item = factor(Item, levels = item_levels)) %>%
  group_by(Area, Item) %>%
  summarize(Number = sum(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Item, values_from = Number)

# Country-Specific Corrections ================================================

# Fix South Africa data for non-commercial farms
# Reference: 2.6M farms from Aliber study
# https://repository.uwc.ac.za/bitstream/handle/10566/4731/aliber_support_smallholder_farmers_south_africa_2012.pdf
SAr <- c(174.67, 238.20, 119.08, 119.08, 85.55, 29.39)  # Ricciardi gross area (ha) https://www.sciencedirect.com/science/article/pii/S2211912417301293
sizemid <- c(0.5, 1, 2.5, 7.5, 15, 25)  # Midpoint of farm size bins
SAr_count <- SAr / sizemid  # Calculate farm counts
toadd <- round(2600000 / sum(SAr_count) * SAr_count)  # Scale to 2.6M total

# Update South Africa data
FAO[FAO$Area == "South Africa", "N0_1"] <- 0
FAO[FAO$Area == "South Africa", c("N0_1", "N1_2", "N2_5", "N5_10", "N10_20", "N20_")] <- 
  FAO[FAO$Area == "South Africa", c("N0_1", "N1_2", "N2_5", "N5_10", "N10_20", "N20_")] + toadd

# Data Validation ==============================================================
# Check consistency between microdata and FAO data
FAO <- replace(FAO, is.na(FAO), 0)

# Calculate proportions for FAO data
FAOprop <- cbind(FAO[, 1], t(apply(FAO[, 2:7], 1, proportions)))
rowSums(FAOprop[, 2:7])

# Calculate proportions for LUGE data (same countries as FAO)
size.sub <- size %>%
  as.data.frame() %>%
  filter(NAME_0 %in% unique(FAO$Area)) %>%
  select(NAME_0:N20_, -one_of("geometry")) %>%
  group_by(NAME_0) %>%
  summarise(across(N0_1:N20_, sum, na.rm = TRUE), .groups = "drop") %>%
  mutate(row_sums = rowSums(select(., N0_1:N20_))) %>%
  mutate(across(N0_1:N20_, ~ . / row_sums)) %>%
  select(-one_of("row_sums"))

# Prepare data for comparison
FAOprop <- FAOprop %>% pivot_longer(cols = N0_1:N20_)
size.sub <- size.sub %>% pivot_longer(cols = N0_1:N20_)

# Join datasets for comparison
all <- FAOprop %>% 
  left_join(size.sub, by = c("Area" = "NAME_0", "name" = "name"),
            suffix = c("FAO", "LUGE"))

# Create comparison plot
comparison_plot <- ggplot(all, aes(valueLUGE, valueFAO, color = name)) +
  geom_point() +
  facet_wrap(~Area, scales = "free") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "FAO vs LUGE Data Comparison",
       x = "LUGE Values", y = "FAO Values") +
  theme_minimal()

# Data Integration =============================================================
# Prepare South Africa fix
other_countries <- size[size$NAME_0 != "South Africa", ]
size_safix <- size %>%
  filter(NAME_0 == "South Africa") %>%
  mutate_at(vars(starts_with("NAME_0"):starts_with("N20_")), ~NA) %>%
  mutate(geometry = st_union(.)) %>%
  slice(1) %>%
  rbind(other_countries, .)

# Update specific countries with FAO data
sa <- replace(size_safix[size_safix$GID_0 == "ZAF", ], 4:9, 
              FAO[FAO$Area == "South Africa", ][2:7])
sa$GID_0 <- "ZAF"
sa$NAME_0 <- "South Africa"

egypt <- replace(size_safix[size_safix$GID_0 == "EGY", ], 4:9, 
                 FAO[FAO$Area == "Egypt", ][2:7])

congo <- replace(size_safix[size_safix$GID_0 == "COG", ], 4:9, 
                 FAO[FAO$Area == "Congo", ][2:7])

# Combine all fixes
allfix <- rbind(sa, congo, egypt)

size_allfix <- size_safix %>% 
  filter(!(GID_0 %in% c("ZAF", "EGY", "COG"))) %>%
  rbind(allfix)

# South Sudan Correction ======================================================
# South Sudan cannot have more than 1M farms based on:
# - Average household size: 5.9
# - Population in 2000: 6,032,267 (World Bank)
# - Estimated households: 1,022,418
# - 95% in farming = 971,297 farms
# References: FAO South Sudan data, Relief Web ISNA 2023

# Calculate current SSD total
ssd_current_total <- size_allfix %>%
  filter(GID_0 == "SSD") %>%
  as.data.frame() %>%
  select(N0_1:N20_) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  rowSums()

cat("Current SSD total farms:", ssd_current_total, "\n")

# Set target and calculate scaling factor
target_total <- 1000000
scaling_factor <- target_total / ssd_current_total
cat("Scaling factor:", scaling_factor, "\n")

# Apply scaling to South Sudan data
size_allfix <- size_allfix %>%
  mutate(
    N0_1 = ifelse(GID_0 == "SSD", N0_1 * scaling_factor, N0_1),
    N1_2 = ifelse(GID_0 == "SSD", N1_2 * scaling_factor, N1_2),
    N2_5 = ifelse(GID_0 == "SSD", N2_5 * scaling_factor, N2_5),
    N5_10 = ifelse(GID_0 == "SSD", N5_10 * scaling_factor, N5_10),
    N10_20 = ifelse(GID_0 == "SSD", N10_20 * scaling_factor, N10_20),
    N20_ = ifelse(GID_0 == "SSD", N20_ * scaling_factor, N20_)
  )

# Verify correction
ssd_new_total <- size_allfix %>%
  filter(GID_0 == "SSD") %>%
  as.data.frame() %>%
  select(N0_1:N20_) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  rowSums()

cat("New SSD total farms:", round(ssd_new_total), "\n")
cat("Difference from target:", round(ssd_new_total - target_total), "\n")

# Save Output ==================================================================
write_sf(size_allfix, "input/fsfix.shp")

cat("Script completed successfully. Output saved to input/fsfix.shp\n")