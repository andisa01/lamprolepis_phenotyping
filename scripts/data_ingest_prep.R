
###
# Phenotyping
###

# Set up the environment
library(tidyverse)
# For mapping
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
# Extras
source("https://raw.githubusercontent.com/andisa01/andis_utils/main/00_HelperFunctions.R")

# Get the raw iNat records
inat_records <- read.csv("./data/inat_records/observations-684399.csv")

inat_records %>%
  group_by(id) %>%
  filter(n() > 1) # The 'id' column is a unique identifier.

inat_records %>%
  select(
    id,
    image_url
  ) %>%
  .[1:10, ] %>%
  write.csv(
    "./data/temp/inat_records_001_010.csv",
    row.names = FALSE
  ) # Export the iNat records for phenotyping

# Read in the phenotyped records
phenotypes <- read.csv("./data/input/Lamprolepis_smaragdina_phenotyping - observations-684399.csv")

phenotypes %>% glimpse()

phenotypes_classed <- phenotypes %>%
  select(
    id,
    observed_on,
    num_identification_agreements,
    num_identification_disagreements,
    latitude,
    longitude,
    positional_accuracy,
    coordinates_obscured,
    phenotype
  ) %>%
  filter(
    !is.na(phenotype)
  ) %>%
  # Correct some errors
  filter(
    !id %in% c(307240727, 40065742, 232483017),
    !id %in% c(20200757, 325619790), # Possible group 5 error
    !id %in% c(54247271), # Definite group 3 error
    !id %in% c(263605039), # Possible group 3 error
    !id %in% c(160645024), # Possible group 2 error
  ) %>%
  mutate(
    phenotype_num = phenotype,
    phenotype = case_when(
      phenotype == 1 ~ "A_philippinica",
      phenotype == 2 ~ "C_trade_green",
      phenotype == 3 ~ "E_viridipuncta",
      phenotype == 4 ~ "D_brown_spot",
      phenotype == 5 ~ "B_saddleback",
      phenotype == 6 ~ "F_hypermelan",
    )
  )

## Mapping phenotypes

# Get country polygons
world <- ne_countries(scale = "medium", returnclass = "sf")

# Get geo points
pts <- phenotypes_classed %>%
  mutate(.row = row_number()) %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

# Geo points for special points
ref_df <- tribble(
  ~label,                               ~longitude, ~latitude, ~label_x, ~label_y,
  "Kosrae\nLa Coquille holotype",       162.98,      5.31,     164.2,     6.3,
  "Makassar\nIndonesia export hub",     119.43,     -5.14,     121.2,    -4.3,
  "Honiara\nSolomon Islands export hub",159.95,     -9.43,     161.5,   -10.4
)

ref_df <- st_as_sf(ref_df, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

# Geo points for place labels
place_df <- tribble(
  ~label,                  ~longitude, ~latitude,
  "Philippines",            122.0,      12.5,
  "North Sulawesi",         124.8,       1.0,
  "South Sulawesi",         120.3,      -3.2,
  "Lesser Sunda Islands",   121.5,      -8.6,
  "Palau",                  134.5,       7.4,
  "Mariana Islands",        145.5,      15.5,
  "Marshall Islands",       170.5,       7.5,
  "New Guinea",             141.5,      -4.8,
  "Bismarck Archipelago",   151.5,      -4.3,
  "Solomon Islands",        160.2,      -8.0
)

pts_xy <- pts %>%
  mutate(longitude = st_coordinates(.)[,1],
         latitude  = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  select(longitude, latitude)

# Expand extent to include BOTH observations and annotation points
extent_sf <- bind_rows(
  pts %>% select(geometry),
  ann_sf %>% select(geometry)
)

# Set bounding box
bb <- st_bbox(pts)
pad_x <- (bb["xmax"] - bb["xmin"]) * 0.10
pad_y <- (bb["ymax"] - bb["ymin"]) * 0.10
xlim <- c(bb["xmin"] - pad_x, bb["xmax"] + pad_x)
ylim <- c(bb["ymin"] - pad_y, bb["ymax"] + pad_y)

# Plot facetedpanel maps
ggplot() +
  geom_sf(data = world, fill = "grey90", linewidth = 0.3) +
  geom_sf(data = pts, aes(color = factor(phenotype)), size = 3, alpha = 1, pch = 16) +
  # geom_sf_text(data = pts, aes(label = .row, color = factor(phenotype)), size = 5) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(color = "Phenotype", x = NULL, y = NULL) +
  # theme_minimal() +
  facet_wrap(vars(factor(phenotype))) +
  scale_color_manual(values = c("blue4", "#377EB8", "#66A61E", "#E6AB02", "red3", "#984EA3"))
ggsave("./maps_panel.jpg", w = 12, h = 6, dpi = 150)

# Plot main map
ggplot() +
  geom_sf(data = world, fill = "grey90", linewidth = 0.3) +
  geom_sf(data = pts, aes(color = factor(phenotype)), size = 2, alpha = 0.85, pch = 16) +
  # geom_sf_text(data = pts %>% filter(phenotype == 1), aes(label = .row, color = factor(phenotype)), size = 5) +
  # geom_point(
  #   data = ref_df,
  #   aes(x = longitude, y = latitude),
  #   inherit.aes = FALSE,
  #   pch = 13,
  #   size = 5,
  #   color = "black"
  # ) +
  # geom_text(
  #   data = ref_df,
  #   aes(x = longitude, y = latitude, label = label),
  #   inherit.aes = FALSE,
  #   size = 3.5,
  #   hjust = 0,
  #   vjust = 1.5
  # ) +
  # geom_text(
  #   data = place_df,
  #   aes(x = longitude, y = latitude, label = label),
  #   inherit.aes = FALSE,
  #   color = "grey20",
  #   size = 3.1,
  #   check_overlap = TRUE
  # ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(color = "Phenotype", x = NULL, y = NULL) +
  # theme_minimal() +
  theme(legend.potion = "none") +
  scale_color_manual(values = c("blue4", "#377EB8", "#66A61E", "#E6AB02", "red3", "#984EA3"))
ggsave("./map.jpg", w = 12, h = 6, dpi = 150)

###
# Import/Export records ====
###

# Get the main lemis dataset and filter for ETS
# lemis <- read.csv("./data/marshal_etal_2025_lemis_data/lemisDataCorrected_2023-11-11.csv/lemisDataCorrected_2023-11-11.csv")
# 
# lemis_codes <- read.csv("./data/marshal_etal_2025_lemis_data/lemis_codes_2023-04-19.csv")
# 
# lemis_lasm <- lemis %>%
#   filter(genus == "LAMPROLEPIS") %>%
#   left_join(
#     lemis_codes %>%
#       select(
#         country_origin = code,
#         country_origin_name = description
#       )
#   ) %>%
#   left_join(
#     lemis_codes %>%
#       select(
#         country_imp_exp = code,
#         country_imp_exp_name = description
#       )
#   )
# 
# write.csv(lemis_lasm, "./data/lemis_lasm.csv", row.names = FALSE) # Write this out so that I can just work with the smaller dataset

lemis_lasm <- read.csv("./data/lemis_lasm.csv") # Grab the smaller dataset

# How many total ETS imported?
lemis_lasm %>%
  mutate(quantity = as.numeric(quantity)) %>%
  summarize(sum(quantity, na.rm = TRUE))

# How many live ETS imported?
lemis_lasm %>%
  mutate(
    quantity = as.numeric(quantity)
  ) %>%
  filter(quantity != 0) %>%
  filter(import_export == "I") %>%
  filter(description == "LIV") %>%
  summarize(sum(quantity, na.rm = TRUE))

# By year and country of origin
lemis_lasm %>%
  mutate(
    quantity = as.numeric(quantity)
  ) %>%
  filter(quantity != 0) %>%
  filter(import_export == "I") %>%
  filter(description == "LIV") %>%
  group_by(sYear, import_export, country_origin_name, country_imp_exp_name) %>%
  summarize(quantity = sum(quantity, na.rm = TRUE)) %>%
  mutate(country = case_when(
    country_origin_name == "Solomon Islands" | country_imp_exp_name == "Solomon Islands" ~ "Solomon Islands",
    country_origin_name == "Indonesia" | country_imp_exp_name == "Indonesia" ~ "Indonesia",
    TRUE ~ "Other"
  )) %>%
  group_by(sYear, country) %>%
  summarize(quantity = sum(quantity)) %>%
  print(n = nrow(.))

# Proportion of total by country of origin
lemis_lasm %>%
  mutate(
    quantity = as.numeric(quantity)
  ) %>%
  filter(quantity != 0) %>%
  filter(!is.na(quantity)) %>%
  filter(import_export == "I") %>%
  filter(description == "LIV") %>%
  group_by(sYear, import_export, country_origin_name, country_imp_exp_name) %>%
  summarize(quantity = sum(quantity, na.rm = TRUE)) %>%
  mutate(country = case_when(
    country_origin_name == "Solomon Islands" | country_imp_exp_name == "Solomon Islands" ~ "Solomon Islands",
    country_origin_name == "Indonesia" | country_imp_exp_name == "Indonesia" ~ "Indonesia",
    TRUE ~ "Other"
  )) %>%
  group_by(sYear, country) %>%
  summarize(quantity = sum(quantity)) %>%
  group_by(country) %>%
  summarize(quantity = sum(quantity)) %>%
  ungroup() %>%
  mutate(
    total = sum(quantity),
    perc = quantity/total*100
  )

# Plot
lemis_lasm %>%
  mutate(
    quantity = as.numeric(quantity)
  ) %>%
  filter(quantity != 0) %>%
  filter(!is.na(quantity)) %>%
  filter(import_export == "I") %>%
  filter(description == "LIV") %>%
  group_by(sYear, import_export, country_origin_name, country_imp_exp_name) %>%
  summarize(quantity = sum(quantity, na.rm = TRUE)) %>%
  mutate(country = case_when(
    country_origin_name == "Solomon Islands" | country_imp_exp_name == "Solomon Islands" ~ "Solomon Islands",
    country_origin_name == "Indonesia" | country_imp_exp_name == "Indonesia" ~ "Indonesia",
    TRUE ~ "Other"
  )) %>%
  group_by(sYear, country) %>%
  summarize(quantity = sum(quantity)) %>%
  ggplot(aes(x = sYear, fill = country)) +
  geom_bar(aes(y = quantity), stat = "identity") +
  # theme_minimal() +
  theme(legend.position = c(0.2, 0.8)) +
  labs(
    x = "",
    title = "Quantity of live L. smaragdina imported to the US by import year",
    subtitle = "LEMIS data from Marshall et al. 2025",
    caption = "A.Z. Andis Arietta | holotypica.com",
    fill = "Country of origin",
    y = ""
  ) +
  scale_fill_manual(values = c("red3", "grey30","#66A61E")) +
  scale_x_continuous(breaks = seq(2015, 2022))
ggsave("./imports_ets_bars.jpg", w = 8, h = 5, dpi = 150)
