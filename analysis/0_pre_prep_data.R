# Prepping data in right format for code
# to reduce mess in clean work through code for paper

box::use(R / dirs[dirs])
box::use(clean = R / clean_data)
box::use(mods = R / for_models)

box::use(ggplot2[...])
box::use(dplyr[...])
box::use(moments)
box::use(sf)
box::use(tidyr)
box::use(magrittr[`%>%`])
box::use(janitor)

# READ IN - Breeding Data -------------------------------------------------

# get all datasets - 1960-2020, 2021, and 2022 and merge them together - save as 1 big dataset

# merge 2021 data
# read in data 1960-2020
breed <- read.csv(file.path(
  dirs$data_raw,
  "LT_breeding_data_greti_bluti_1960_2020.csv"
), na.strings = c("", "NA")) %>%
  janitor::clean_names()
nrow(breed) # 39530

# read in 2021 data
breed2021 <- read.csv(file.path(
  dirs$data_raw,
  "ebmp_broods_2021.csv"
), na.strings = c("", "NA")) %>%
  # add column with year
  dplyr::mutate(year = as.integer(2021)) %>%
  janitor::clean_names()
nrow(breed2021) # 1209

# read in 2022 data
breed2022 <- read.csv(file.path(
  dirs$data_raw,
  "ebmp_broods_2022.csv"
), na.strings = c("", "NA")) %>%
  # add column with year
  dplyr::mutate(year = as.integer(2022)) %>%
  janitor::clean_names()
nrow(breed2022) # 1248


# bind 2
breed <- rbind(breed, breed2021, breed2022)
nrow(breed) # 41987


# clean dataset
breed <- clean$clean_breeding_data(breed)
nrow(breed) # 40254 #40248

write.csv(breed, file = file.path(
  dirs$data_output,
  "LT_breeding_data_greti_bluti_1960_2022.csv"
), row.names = F)


# READ IN - Ringing Data -------------------------------------------------

# get ringing data sets and add together -
# first join 2013-2020 data with 2021 and 2022
ring2 <- read.csv(file.path(dirs$data_raw, "ebmp_database_ringing_record_export_GT&BT_2013-20.csv"),
  na.strings = c("", "NA")
) %>%
  janitor::clean_names()

ring2021 <- read.csv(file.path(dirs$data_raw, "ebmp_ringing_data_2021.csv"),
  na.strings = c("", "NA")
) %>%
  janitor::clean_names()

ring2022 <- read.csv(file.path(dirs$data_raw, "ebmp_ringing_2022.csv"),
  na.strings = c("", "NA")
) %>%
  janitor::clean_names()


nrow(ring2) # 51566
nrow(ring2021) # 8197
nrow(ring2022) # 7469

# bind
ring2 <- rbind(ring2, ring2021, ring2022)
nrow(ring2) # 67232


# READ IN - Ringing Data up to 2013 --------------------------------------------------

# then also add data 1960-2013
ring <- read.csv(file.path(
  dirs$data_raw,
  "legacy_ringing_records_GT&BT_up_to_2013.csv"
), na.strings = c("", "NA")) %>%
  janitor::clean_names() %>%
  # clean with function - removes errors and cleans up dataframe
  clean$clean_ringing_data() %>%
  # keep only some columns
  select(pnum, age, sex, bto_species_code, bto_ring, yr, nb, retrap)
nrow(ring) # 179010

# Clean Ringing data 2013-2022  -------------------------------------

ring2_clean <- ring2 %>%
  # clean with function - removes errors and cleans up dataframe
  clean$clean_ringing_data_2(.) %>%
  mutate(pnum = paste0(date, "1", site)) %>%
  # rename columns
  rename(bto_ring = ring,
                bto_species_code = spec,
                yr = date,
                nb = site,
                retrap = rtype,
                location = place
  ) %>%
  # keep only selected columns
  dplyr::select(pnum, age, sex, bto_species_code, bto_ring, yr, nb, retrap)

nrow(ring2_clean) # 61343

# bind
ring_all <- rbind(ring, ring2_clean)
nrow(ring_all) # 240353

# save
write.csv(ring_all, file = file.path(
  dirs$data_output,
  "ebmp_database_ringing_record_export_GT&BT_all.csv"
), row.names = F)




# SPATIAL DATA ----------------------------------------------

# make combined dataset of all nestbox environment data

# read in nest box data and use some classifcations of habitat
nestbox <- read.csv(file.path(dirs$data_raw, "Nest_box_habitat_data.csv"),
  na.strings = c("", "NA")
) %>%
  janitor::clean_names() %>%
  dplyr::select("box", "section", "x", "y", "edge_edi", "altitude_m", "northness")
nrow(nestbox) # 1019

# oak data - not for all boxes.... leave some with NA
oak <- read.csv(file.path(dirs$data_raw, "Tit_oak_data.csv"), na.strings = c("", "NA")) %>%
  janitor::clean_names() %>%
  dplyr::select("box", "no_trees_75m")
nrow(oak) # 987

box_tree <- merge(nestbox, oak, by = "box", all.x = T)
head(box_tree)
nrow(box_tree) # 1019


# get territory areas ----------------------------------------------

# open wood outline file - extract first polygon and transform
wood_outline <- sf::st_read(file.path(dirs$data_raw, "/maps/perimeter poly with clearings_region.shp"))[1, ] %>%
  sf::st_transform(27700)

# need breeding data with nest box locations
nrow(breed) # 40249

# sort breeding data
breeding_data <- breed %>%
  dplyr::inner_join(., nestbox, by = c("nest_box" = "box")) %>%
  dplyr::filter(species == "g") %>%
  # convert to spatial object
  sf::st_as_sf(coords = c("x", "y"), remove = F, crs = 27700)
nrow(breeding_data) # 17809

# loop to get territory size for individuals within each year

# NEED TO SAVE AND CLOSE AND REOPEN - WASNT GETTING FUNCTION INSIDE FUNCTION***!SECTION


GTIT_allyrs_area <- NULL
for (yr in unique(breeding_data$year)) {
  territories_G <- mods$get_territory_polygons(breeding_data, wood_outline, yr)
  GTIT_allyrs_area <- rbind(GTIT_allyrs_area, territories_G)
}

summary(GTIT_allyrs_area$area_polygon)
nrow(GTIT_allyrs_area) # 17809

# select just some columns and drop geometyr
territory <- GTIT_allyrs_area %>%
  # get rid of geometry
  sf::st_drop_geometry() %>%
  dplyr::mutate(area_polygon = as.numeric(area_polygon)) %>%
  dplyr::select("pnum", "nest_box", "area_polygon")

# merge with other habitat data
box_tree_terr <- merge(box_tree, territory, by.x = "box", by.y = "nest_box", all.y = T) %>%
  dplyr::mutate(area_polygon_sqrt = sqrt(area_polygon))

nrow(box_tree_terr) # 17809

# save
write.csv(box_tree_terr, file = file.path(dirs$data_output, "Habitat_data_nestboxes.csv"), row.names = F)
