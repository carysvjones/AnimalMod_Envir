#Prepping data in right format for code 
#to reduce mess in clean work through code for paper

box::use(.. / R / dirs[dirs])
box::use(clean = ../ R / clean_data)
box::use(mods = ../ R / for_models)

box::use(ggplot2[...])
box::use(dplyr[...])
box::use(moments)
box::use(sf)
box::use(tidyr)
box::use(magrittr[`%>%`])

# READ IN - Breeding Data -------------------------------------------------

#get all datasets - 1960-2020, 2021, and 2022 and merge them together - save as 1 big dataset

#merge 2021 data
#read in data 1960-2020
breed <- read.csv(file.path(dirs$data_raw, 'LT_breeding_data_greti_bluti_1960_2020.csv'), na.strings=c("", "NA"))
nrow(breed) #39530

#read in 2021 data
breed2021 <- read.csv(file.path(dirs$data_raw, 'ebmp_broods_2021.csv'), na.strings=c("", "NA")) %>%
  #add column with year
  dplyr::mutate(year = as.integer(2021))
nrow(breed2021) #1209

#read in 2022 data 
breed2022 <- read.csv(file.path(dirs$data_raw, 'ebmp_broods_2022.csv'), na.strings=c("", "NA")) %>%
  #add column with year
  dplyr::mutate(year = as.integer(2022))
nrow(breed2022) #1248


#bind 2
breed <- rbind(breed, breed2021, breed2022)
nrow(breed) #41987


# clean dataset 
breed <- clean$clean_breeding_data(breed)
nrow(breed) #38017

write.csv(breed, file = file.path(dirs$data_output,'LT_breeding_data_greti_bluti_1960_2022.csv'), row.names = F)


# READ IN - Ringing Data -------------------------------------------------

#get ringing data sets and add together - 
#first join 2013-2020 data with 2021 and 2022
ring2 <- read.csv(file.path(dirs$data_raw, 'ebmp_database_ringing_record_export_GT&BT_2013-20.csv'),
                  na.strings=c("", "NA"))

ring2021 <- read.csv(file.path(dirs$data_raw, 'ebmp_ringing_data_2021.csv'),
                  na.strings=c("", "NA"))

ring2022 <- read.csv(file.path(dirs$data_raw, 'ebmp_ringing_2022.csv'), 
                     na.strings=c("", "NA")) 


nrow(ring2) #51566
nrow(ring2021) #8197
nrow(ring2022) #7469

#bind
ring2 <- rbind(ring2, ring2021, ring2022)
nrow(ring2) #67232


# READ IN - Ringing Data up to 2013 --------------------------------------------------

#then also add data 1960-2013
ring <- read.csv(file.path(dirs$data_raw, 'legacy_ringing_records_GT&BT_up_to_2013.csv'), na.strings=c("", "NA")) %>%
  #clean with function - removes errors and cleans up dataframe
  clean$clean_ringing_data(.) %>%
  #keep only some columns
  select(Pnum, age, sex, bto_species_code, bto_ring, yr, nb, retrap)
nrow(ring) #179010

# Clean Ringing data 2013-2022  -------------------------------------

ring2_clean <- ring2 %>% 
  #clean with function - removes errors and cleans up dataframe
  clean$clean_ringing_data_2(.) %>%
  mutate(Pnum = paste0(Date, '1', Site)) %>%
  #rename columns
  rename(bto_ring = Ring,
                bto_species_code = Spec,
                age = Age,
                sex = Sex,
                yr = Date,
                nb = Site,
                retrap = Rtype,
                location = Place
  ) %>%
  #keep only selected columns
  dplyr::select(Pnum, age, sex, bto_species_code, bto_ring, yr, nb, retrap)

nrow(ring2_clean) #61343

#bind
ring_all <- rbind(ring, ring2_clean)
nrow(ring_all) #240353

#save
write.csv(ring_all, file = file.path(dirs$data_output, 'ebmp_database_ringing_record_export_GT&BT_all.csv'), row.names = F)




# SPATIAL DATA ----------------------------------------------

#make combined dataset of all nestbox environment data

#read in nest box data and use some classifcations of habitat
nestbox <- read.csv(file.path(dirs$data_raw, 'Nest_box_habitat_data.csv'), na.strings=c("", "NA"))
nrow(nestbox) #1019
nestbox <- nestbox %>%
  dplyr::select('Box', 'Section', 'x', 'y', 'edge.EDI.', 'altitude.m.', 'northness')

#oak data - not for all boxes.... leave some with NA
oak <- read.csv(file.path(dirs$data_raw, 'Tit_oak_data.csv'), na.strings=c("", "NA"))
nrow(oak) #987
oak <- oak %>%
  dplyr::select('Box', 'No_trees_75m')

box_tree <- merge(nestbox, oak, by = 'Box', all.x = T)
head(box_tree)
nrow(box_tree) #1019

  
# get territory areas ----------------------------------------------

#open wood outline file - extract first polygon and transform
wood_outline <- sf::st_read(file.path(dirs$data_raw, '/maps/perimeter poly with clearings_region.shp'))[1,] %>%
  sf::st_transform(27700)

#need breeding data with nest box locations 
nrow(breed) #40255

#sort breeding data 
breeding_data <- breed %>%
  dplyr::inner_join(., nestbox, by = c('nest.box' = 'Box')) %>%
  dplyr::filter(Species == 'g') %>%
  #convert to spatial object
  sf::st_as_sf(coords=c("x","y"), remove=F, crs=27700)
nrow(breeding_data) #17812

#loop to get territory size for individuals within each year
GTIT_allyrs_area <- NULL
for (yr in unique(breeding_data$year)){
  territories_G <- get_territory_polygons(breeding_data, wood_outline, yr)
  GTIT_allyrs_area <- rbind(GTIT_allyrs_area, territories_G)
}

summary(GTIT_allyrs_area$area_polygon)
nrow(GTIT_allyrs_area) #17812

#select just some columns and drop geometyr
territory <- GTIT_allyrs_area %>%
  #get rid of geometry 
  sf::st_drop_geometry() %>%
  dplyr::mutate(area_polygon = as.numeric(area_polygon)) %>%
  dplyr::select('Pnum', 'nest.box', 'area_polygon')

#merge with other habitat data
box_tree_terr <- merge(box_tree, territory, by.x = 'Box', by.y = 'nest.box', all.y = T)
nrow(box_tree_terr) #17812

#save
write.csv(box_tree_terr, file = file.path(dirs$data_output,'Habitat_data_nestboxes.csv'), row.names = F)



