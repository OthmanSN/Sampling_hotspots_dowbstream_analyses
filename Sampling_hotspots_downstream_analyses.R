rm(list = ls())

library(tidyverse)
library(ggalluvial)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)

setwd("C:/Users/Hp/Desktop/NFU2/mismatchA1vsC1")

#-------------------------
# 1. READ FILES
#-------------------------
alluvial_df <- read.csv("dicroglossid_alluvial_input.csv", stringsAsFactors = FALSE)
full_df     <- read.csv("dicroglossid_modelA1_vs_genbank_full_georef.csv", stringsAsFactors = FALSE)
pair_df     <- read.csv("dicroglossid_pair_summary.csv", stringsAsFactors = FALSE)
region_df   <- read.csv("dicroglossid_region_summary.csv", stringsAsFactors = FALSE)

#-------------------------
# 2. CLEANING FUNCTIONS
#-------------------------
clean_taxon <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\s+", " ", x)
  x <- gsub("multistrata", "multistriata", x, ignore.case = TRUE)
  x <- gsub("multistriaat", "multistriata", x, ignore.case = TRUE)
  x <- gsub("rugolosus", "rugulosus", x, ignore.case = TRUE)
  x
}

strip_clade <- function(x) {
  x <- clean_taxon(x)
  x <- gsub("\\s*\\(Clade [IVX]+\\)", "", x, ignore.case = TRUE)
  trimws(x)
}

#-------------------------
# 3. MAIN DATA CLEANING
#-------------------------

# A1 label
if ("A1_full" %in% names(full_df)) {
  full_df$A1_full <- clean_taxon(full_df$A1_full)
} else if ("Species.assignment.of.Model.A1" %in% names(full_df)) {
  full_df$A1_full <- clean_taxon(full_df$Species.assignment.of.Model.A1)
} else if ("a1_assignment" %in% names(full_df)) {
  full_df$A1_full <- clean_taxon(full_df$a1_assignment)
} else {
  stop("Cannot find A1 assignment column in full_df")
}

# GenBank label
if ("GenBank_full" %in% names(full_df)) {
  full_df$GenBank_full <- clean_taxon(full_df$GenBank_full)
} else if ("Taxon.label.in.GenBank" %in% names(full_df)) {
  full_df$GenBank_full <- clean_taxon(full_df$Taxon.label.in.GenBank)
} else if ("genbank_label_original" %in% names(full_df)) {
  full_df$GenBank_full <- clean_taxon(full_df$genbank_label_original)
} else if ("genbank_species" %in% names(full_df)) {
  full_df$GenBank_full <- clean_taxon(full_df$genbank_species)
} else {
  stop("Cannot find GenBank label column in full_df")
}

# stripped species labels
full_df$A1_species      <- strip_clade(full_df$A1_full)
full_df$GenBank_species <- strip_clade(full_df$GenBank_full)

# match status
if ("match_status" %in% names(full_df)) {
  full_df$match_status <- tolower(trimws(as.character(full_df$match_status)))
} else {
  full_df$match_status <- ifelse(full_df$A1_full == full_df$GenBank_full, "match", "mismatch")
}

# country
if ("country" %in% names(full_df)) {
  full_df$country <- trimws(as.character(full_df$country))
} else if ("Country" %in% names(full_df)) {
  full_df$country <- trimws(as.character(full_df$Country))
} else {
  full_df$country <- NA_character_
}
full_df$country[full_df$country %in% c("Verified", "", "Unknown")] <- NA

# locality
if ("locality" %in% names(full_df)) {
  full_df$locality <- trimws(as.character(full_df$locality))
} else if ("Locality" %in% names(full_df)) {
  full_df$locality <- trimws(as.character(full_df$Locality))
} else {
  full_df$locality <- NA_character_
}

# coordinates
if ("lon" %in% names(full_df)) full_df$lon <- as.numeric(trimws(as.character(full_df$lon)))
if ("lat" %in% names(full_df)) full_df$lat <- as.numeric(trimws(as.character(full_df$lat)))

# region_group:
# keep your manually revised values if already present;
# only autofill blanks/missing values
if ("region_group" %in% names(full_df)) {
  full_df$region_group <- trimws(as.character(full_df$region_group))
} else {
  full_df$region_group <- NA_character_
}

missing_region <- is.na(full_df$region_group) | full_df$region_group == ""

full_df$region_group[missing_region] <- case_when(
  grepl("hainan", full_df$locality[missing_region], ignore.case = TRUE) ~ "Hainan (insular)",
  grepl("taiwan", full_df$locality[missing_region], ignore.case = TRUE) |
    full_df$country[missing_region] == "Taiwan" ~ "Taiwan (insular)",
  full_df$country[missing_region] == "Japan" |
    grepl("japan|ryukyu|iriomote", full_df$locality[missing_region], ignore.case = TRUE) ~ "Japan (archipelago)",
  grepl("jiangsu|anhui|henan|shanghai|chongming|yangtze|huai",
        full_df$locality[missing_region], ignore.case = TRUE) ~ "Yangtze–Huai transition",
  full_df$country[missing_region] %in% c("Thailand", "Vietnam", "Laos", "Myanmar", "Cambodia") ~ "Mainland Southeast Asia",
  full_df$country[missing_region] %in% c("Malaysia", "Indonesia", "Philippines") |
    grepl("borneo|sumatra|java|sarawak|sabah",
          full_df$locality[missing_region], ignore.case = TRUE) ~ "Sundaic (Maritime Southeast Asia)",
  full_df$country[missing_region] == "China" ~ "East Asia mainland",
  TRUE ~ "Other / unclear"
)

full_df$region_group <- factor(
  full_df$region_group,
  levels = c(
    "Hainan (insular)",
    "Taiwan (insular)",
    "Japan (archipelago)",
    "Yangtze–Huai transition",
    "East Asia mainland",
    "Mainland Southeast Asia",
    "Sundaic (Maritime Southeast Asia)",
    "Other / unclear"
  )
)

# quick checks
print(sort(unique(full_df$region_group)))
print(table(full_df$match_status, useNA = "ifany"))
print(summary(full_df$lon))
print(summary(full_df$lat))

#-------------------------
# 4. SPECIES ALLUVIAL
#-------------------------
species_df <- full_df %>%
  count(GenBank_species, A1_species, match_status, name = "Freq")

p_species <- ggplot(
  species_df,
  aes(axis1 = GenBank_species, axis2 = A1_species, y = Freq/5)
) +
  geom_alluvium(aes(fill = match_status), width = 0.1, alpha = 0.8) +
  geom_stratum(width = 0.25, fill = "lightgrey", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 0) +
  scale_x_discrete(limits = c("GenBank / C1 label", "Model A1 assignment")) +
  labs(x = NULL, y = "Number of entries", fill = "Status") +
  theme_classic(base_size = 12)

ggsave("Fig_main_species_alluvial_revise_nolabel.pdf", p_species, width = 6, height = 12)
ggsave("Fig_main_species_alluvial_revise_nolabel.png", p_species, width = 6, height = 8, dpi = 300)

#-------------------------
# 5. CLADE ALLUVIAL
#-------------------------
clade_df <- full_df %>%
  count(GenBank_full, A1_full, match_status, name = "Freq")

p_clade <- ggplot(
  clade_df,
  aes(axis1 = GenBank_full, axis2 = A1_full, y = Freq)
) +
  geom_alluvium(aes(fill = match_status), width = 0.2, alpha = 0.8) +
  geom_stratum(width = 0.25, fill = "lightgrey", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 1.5)+
  scale_x_discrete(limits = c("GenBank / C1 label", "Model A1 assignment")) +
  labs(x = NULL, y = "Number of entries", fill = "Status") +
  theme_classic(base_size = 11)

ggsave("FigS_clade_alluvial_revised_label.pdf", p_clade, width = 9, height = 14)
ggsave("FigS_clade_alluvial_revised_label.png", p_clade, width = 9, height = 14, dpi = 300)

#-------------------------
# 6. COUNTRY BARPLOT
#-------------------------
country_df <- full_df %>%
  filter(!is.na(country), country != "") %>%
  count(country, match_status)

p_country_bar <- ggplot(
  country_df,
  aes(x = reorder(country, n), y = n, fill = match_status)
) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = "Number of entries", fill = "Status") +
  theme_classic(base_size = 12)

ggsave("Fig_country_barplot_updated.pdf", p_country_bar, width = 8, height = 6)
ggsave("Fig_country_barplot_updated.png", p_country_bar, width = 8, height = 6, dpi = 300)

#-------------------------
# 7. REGION BARPLOT
#-------------------------
region_df_plot <- full_df %>%
  filter(!is.na(region_group), region_group != "") %>%
  count(region_group, match_status)

p_region_bar <- ggplot(
  region_df_plot,
  aes(x = reorder(region_group, n), y = n, fill = match_status)
) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = "Number of entries", fill = "Status") +
  theme_classic(base_size = 12)

ggsave("Fig_region_barplot_updated.pdf", p_region_bar, width = 8, height = 6)
ggsave("Fig_region_barplot_updated.png", p_region_bar, width = 8, height = 6, dpi = 300)

#-------------------------
# 8. REGION HEATMAP
#-------------------------
heat_region_df <- full_df %>%
  filter(!is.na(region_group), region_group != "", !is.na(A1_species), A1_species != "") %>%
  count(region_group, A1_species, name = "n")

heat_region_df$region_group <- factor(
  heat_region_df$region_group,
  levels = levels(full_df$region_group)
)

p_region_heat <- ggplot(
  heat_region_df,
  aes(x = region_group, y = A1_species, fill = n)
) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "black") +
  labs(x = NULL, y = "Model A1 assignment", fill = "Count") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Fig_region_heatmap_updated.pdf", p_region_heat, width = 10, height = 7)
ggsave("Fig_region_heatmap_updated.png", p_region_heat, width = 10, height = 7, dpi = 300)

#-------------------------
# 9. PAIR HEATMAP
#-------------------------
pair_df$genbank_species <- clean_taxon(pair_df$genbank_species)
pair_df$a1_assignment   <- clean_taxon(pair_df$a1_assignment)

p_pair_heat <- ggplot(pair_df, aes(x = genbank_species, y = a1_assignment, fill = n)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "black") +
  labs(x = "GenBank / C1 label", y = "Model A1 assignment", fill = "Count") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Fig_pair_heatmap.pdf", p_pair_heat, width = 10, height = 7)
ggsave("Fig_pair_heatmap.png", p_pair_heat, width = 10, height = 7, dpi = 300)

#-------------------------
# 10. MISMATCH-ONLY PAIR HEATMAP
#-------------------------
pair_df2 <- pair_df %>%
  filter(match_status == "mismatch")

genbank_order <- pair_df2 %>%
  group_by(genbank_species) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(genbank_species)

a1_order <- pair_df2 %>%
  group_by(a1_assignment) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(a1_assignment)

pair_df2$genbank_species <- factor(pair_df2$genbank_species, levels = genbank_order)
pair_df2$a1_assignment   <- factor(pair_df2$a1_assignment, levels = a1_order)

p_pair_heat_mismatch <- ggplot(pair_df2, aes(x = genbank_species, y = a1_assignment, fill = n)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "black") +
  labs(x = "GenBank / C1 label", y = "Model A1 assignment", fill = "Mismatch count") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10)
  )

ggsave("Fig_pair_heatmap_mismatch.pdf", p_pair_heat_mismatch, width = 10, height = 7)
ggsave("Fig_pair_heatmap_mismatch.png", p_pair_heat_mismatch, width = 10, height = 7, dpi = 300)

#-------------------------
# 11. CHI-SQUARE TEST
#-------------------------
tab <- xtabs(n ~ genbank_species + a1_assignment, data = pair_df)

chisq_res <- chisq.test(tab, simulate.p.value = TRUE, B = 9999)
print(chisq_res)

tab_df <- as.data.frame(as.table(tab))
colnames(tab_df) <- c("GenBank_label", "A1_assignment", "Count")

expected_df <- as.data.frame(as.table(chisq_res$expected))
colnames(expected_df) <- c("GenBank_label", "A1_assignment", "Expected_count")

stdres_df <- as.data.frame(as.table(chisq_res$stdres))
colnames(stdres_df) <- c("GenBank_label", "A1_assignment", "Std_residual")

obs_exp_df <- tab_df %>%
  left_join(expected_df, by = c("GenBank_label", "A1_assignment")) %>%
  mutate(deviation = Count - Expected_count) %>%
  left_join(stdres_df, by = c("GenBank_label", "A1_assignment"))

write.csv(tab_df, "Table_SX_GenBank_vs_A1_contingency.csv", row.names = FALSE)
write.csv(expected_df, "Table_SX_expected_counts.csv", row.names = FALSE)
write.csv(obs_exp_df, "Table_SX_observed_vs_expected.csv", row.names = FALSE)

#-------------------------
# 12. MAP WITH POINTS
#-------------------------
world <- ne_countries(scale = "medium", returnclass = "sf")

full_df_map <- full_df %>%
  filter(!is.na(lon), !is.na(lat))

pts_sf <- st_as_sf(full_df_map, coords = c("lon", "lat"), crs = 4326)

p_points <- ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_sf(data = pts_sf, aes(color = match_status), size = 1.8, alpha = 0.8) +
  coord_sf(xlim = c(90, 140), ylim = c(0, 40), expand = FALSE) +
  labs(color = "Status") +
  theme_classic()

ggsave("Fig_map_match_status.pdf", p_points, width = 8, height = 6)
ggsave("Fig_map_match_status.png", p_points, width = 8, height = 6, dpi = 300)

#-------------------------
# 13. MISMATCH-ONLY MAP
#-------------------------
mismatch_map_df <- full_df %>%
  filter(match_status == "mismatch", !is.na(lon), !is.na(lat))

mismatch_pts_sf <- st_as_sf(mismatch_map_df, coords = c("lon", "lat"), crs = 4326)

p_mismatch_map <- ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_sf(data = mismatch_pts_sf, color = "black", size = 1.8, alpha = 0.8) +
  coord_sf(xlim = c(90, 140), ylim = c(0, 40), expand = FALSE) +
  theme_classic()

ggsave("Fig_map_mismatch_only.pdf", p_mismatch_map, width = 8, height = 6)
ggsave("Fig_map_mismatch_only.png", p_mismatch_map, width = 8, height = 6, dpi = 300)



#-------------------------
# REGION BUBBLE MAP
# size = total records
# fill = mismatch proportion
#-------------------------

region_summary_map <- full_df %>%
  filter(!is.na(region_group), region_group != "",
         region_group != "Other / unclear",
         !is.na(match_status), match_status != "") %>%
  count(region_group, match_status, name = "n") %>%
  group_by(region_group) %>%
  mutate(total_region = sum(n)) %>%
  ungroup() %>%
  filter(match_status == "mismatch") %>%
  mutate(mismatch_prop = n / total_region)

# approximate centroids
region_coords <- data.frame(
  region_group = c(
    "Hainan (insular)",
    "Taiwan (insular)",
    "Japan (archipelago)",
    "Yangtze–Huai transition",
    "East Asia mainland",
    "Mainland Southeast Asia",
    "Sundaic (Maritime Southeast Asia)"
  ),
  lon = c(110.3, 121.0, 138.0, 119.0, 112.0, 101.0, 113.5),
  lat = c(19.2, 23.7, 36.0, 32.0, 35.0, 15.5, 2.5)
)

region_summary_map <- left_join(region_summary_map, region_coords, by = "region_group")

world <- ne_countries(scale = "medium", returnclass = "sf")

p_region_bubble2 <- ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey70", linewidth = 0.2) +
  
  # background raw datapoints (FIRST)
  geom_point(
    data = full_df,
    aes(x = lon, y = lat),
    color = "black",
    size = 0.6,
    alpha = 0.4
  ) +
  
  geom_point(
    data = region_summary_map,
    aes(x = lon, y = lat, size = total_region, fill = mismatch_prop),
    shape = 21,
    color = "grey95",
    stroke = 0.5,
    alpha = 0.65
  ) +
  coord_sf(xlim = c(90, 145), ylim = c(-10, 45), expand = FALSE) +
  scale_fill_gradient(
    low = "#c7e9c0",
    high =  "#1b9e9a",
    name = "Mismatch proportion"
  ) +
  scale_size_continuous(name = "Total records", range = c(3, 12)) +
  theme_classic(base_size = 12)

print(p_region_bubble2)

ggsave("Fig_region_bubble_total_and_mismatch_wide_withrawdatapoint.pdf", p_region_bubble2, width = 8, height = 6)
ggsave("Fig_region_bubble_total_and_mismatch_wide_withrawdatapoint.png", p_region_bubble2, width = 8, height = 6, dpi = 300)

#Save proportion value
region_table <- full_df %>%
  filter(!is.na(region_group), region_group != "",
         region_group != "Other / unclear",
         !is.na(match_status), match_status != "") %>%
  count(region_group, match_status, name = "n") %>%
  group_by(region_group) %>%
  mutate(
    total_region = sum(n),
    proportion = n / total_region
  ) %>%
  ungroup() %>%
  tidyr::pivot_wider(
    names_from = match_status,
    values_from = c(n, proportion),
    values_fill = 0
  )


#rename the column
region_table <- region_table %>%
  rename(
    match_n = n_match,
    mismatch_n = n_mismatch,
    match_prop = proportion_match,
    mismatch_prop = proportion_mismatch
  )


write.csv(region_table,
          "Table_region_match_mismatch_summary.csv",
          row.names = FALSE)



###BIGGER SCALE VERSION OF MISMATCH PROPORTION FOR PUBLICATION
#-------------------------
# REGION BUBBLE MAP
# size = total records
# fill = mismatch proportion
# with background raw datapoints
#-------------------------

library(tidyverse)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)

# version tag so files do not overwrite
version_tag <- "v5_green_withrawpoints_sundaiconly"

#-------------------------
# REGION BUBBLE MAP
# size = total records
# fill = mismatch proportion
# bubble position = median coordinates of raw datapoints in each region
# with small manual offset only for Sundaic region
#-------------------------

# summarize mismatch by region
region_summary_map <- full_df %>%
  filter(!is.na(region_group), region_group != "",
         region_group != "Other / unclear",
         !is.na(match_status), match_status != "") %>%
  count(region_group, match_status, name = "n") %>%
  group_by(region_group) %>%
  mutate(total_region = sum(n)) %>%
  ungroup() %>%
  filter(match_status == "mismatch") %>%
  mutate(mismatch_prop = n / total_region)

# compute region bubble positions from actual datapoints
region_coords <- full_df %>%
  filter(!is.na(region_group), region_group != "",
         region_group != "Other / unclear",
         !is.na(lon), !is.na(lat)) %>%
  group_by(region_group) %>%
  summarise(
    lon = median(lon, na.rm = TRUE),
    lat = median(lat, na.rm = TRUE),
    .groups = "drop"
  )

# join region summaries to median coordinates
region_summary_map <- left_join(region_summary_map, region_coords, by = "region_group")

# small manual offset only for Sundaic / Borneo readability
region_summary_map <- region_summary_map %>%
  mutate(
    lon = case_when(
      region_group == "Sundaic (Maritime Southeast Asia)" ~ lon + 0.6,
      TRUE ~ lon
    ),
    lat = case_when(
      region_group == "Sundaic (Maritime Southeast Asia)" ~ lat + 0.0,
      TRUE ~ lat
    )
  )

# keep only valid raw datapoints
full_df_map <- full_df %>%
  filter(!is.na(lon), !is.na(lat))

# world basemap
world <- ne_countries(scale = "medium", returnclass = "sf")

# plot
p_region_bubble5 <- ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey70", linewidth = 0.2) +
  
  # background raw datapoints
  geom_point(
    data = full_df_map,
    aes(x = lon, y = lat),
    color = "black",
    size = 0.6,
    alpha = 0.4
  ) +
  
  # regional summary bubbles
  geom_point(
    data = region_summary_map,
    aes(x = lon, y = lat, size = total_region, fill = mismatch_prop),
    shape = 21,
    color = "grey95",
    stroke = 0.5,
    alpha = 0.65
  ) +
  
  coord_sf(xlim = c(90, 145), ylim = c(-10, 45), expand = FALSE) +
  
  scale_fill_gradient(
    low = "#c7e9c0",
    high = "#1b9e9a",
    name = "Mismatch proportion"
  ) +
  
  scale_size_continuous(
    name = "Total records",
    range = c(3, 12)
  ) +
  
  labs(x = "lon", y = "lat") +
  theme_classic(base_size = 12)

print(p_region_bubble5)

# save figure
ggsave(
  paste0("Fig_region_bubble_", version_tag, ".pdf"),
  p_region_bubble5,
  width = 8, height = 6
)

ggsave(
  paste0("Fig_region_bubble_", version_tag, ".png"),
  p_region_bubble5,
  width = 8, height = 6, dpi = 300
)

#-------------------------
# SAVE REGION SUMMARY TABLE
#-------------------------

region_table <- full_df %>%
  filter(!is.na(region_group), region_group != "",
         region_group != "Other / unclear",
         !is.na(match_status), match_status != "") %>%
  count(region_group, match_status, name = "n") %>%
  group_by(region_group) %>%
  mutate(
    total_region = sum(n),
    proportion = n / total_region
  ) %>%
  ungroup() %>%
  tidyr::pivot_wider(
    names_from = match_status,
    values_from = c(n, proportion),
    values_fill = 0
  ) %>%
  rename(
    match_n = n_match,
    mismatch_n = n_mismatch,
    match_prop = proportion_match,
    mismatch_prop = proportion_mismatch
  )

write.csv(
  region_table,
  paste0("Table_region_match_mismatch_summary_", version_tag, ".csv"),
  row.names = FALSE
)

#-------------------------
# SAVE MAP INPUT TABLE
#-------------------------

write.csv(
  region_summary_map,
  paste0("Table_region_map_input_", version_tag, ".csv"),
  row.names = FALSE
)

#-------------------------
# SAVE RANGE OF MISMATCH PROPORTION
#-------------------------

range_vals <- range(region_table$mismatch_prop, na.rm = TRUE)

writeLines(
  paste0(
    "Mismatch proportion range: ",
    round(range_vals[1], 3), " – ", round(range_vals[2], 3)
  ),
  con = paste0("Table_region_range_", version_tag, ".txt")
)