#Prepping data in right format for code 
#to reduce mess in clean work through code for paper

box::use(.. / R / dirs[dirs])
box::use(clean = ../ R / clean_data)
box::use(mods = ../ R / for_models)

box::use(ggplot2)
box::use(dplyr)
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
head(breed2021)

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
  dplyr::select(Pnum, age, sex, bto_species_code, bto_ring, yr, nb, retrap)
nrow(ring) #179010

# Clean Ringing data 2013-2022  -------------------------------------

ring2_clean <- ring2 %>% 
  #clean with function - removes errors and cleans up dataframe
  clean$clean_ringing_data_2(.) %>%
  dplyr::mutate(Pnum = paste0(Date, '1', Site)) %>%
  #rename columns
  dplyr::rename(bto_ring = Ring,
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

#open wood outline file
wood.outline <- sf::st_read(file.path(dirs$data_raw, '/maps/perimeter poly with clearings_region.shp'))

#extracting the first polygon because that's the bit of the woods we care about
#(and the polygons go a bit nuts if you don't)
wood.outline <- wood.outline[1,]
#need to transform wood outline
wood.outline <- sf::st_transform(wood.outline, 27700)

#need breeding data with nest box locations 
nrow(breed) #40255

breeding.data <- breed %>%
  dplyr::inner_join(., nestbox, by = c('nest.box' = 'Box'))
nrow(breeding.data) #37620

#get subset of just great tits
breeding.data.G <- subset(breeding.data, Species == 'g' )

#converting it into a spatial object
breeding.data.G <- sf::st_as_sf(breeding.data.G, coords=c("x","y"), remove=F, crs=27700)

#calculating a bounding box for the function that calculates the territory polygons
#can use this for when include blue tits as well later
bbox_polygon_G <- function(x) {
  bb <- sf::st_bbox(x)
  
  p <- matrix(
    c(bb["xmin"], bb["ymin"], 
      bb["xmin"], bb["ymax"],
      bb["xmax"], bb["ymax"], 
      bb["xmax"], bb["ymin"], 
      bb["xmin"], bb["ymin"]),
    ncol = 2, byrow = T
  )
  
  sf::st_polygon(list(p))
}
box <- sf::st_sfc(bbox_polygon_G(breeding.data.G))


#loop to get territory size for individuals within each year
GTIT_allyrs_area <- NULL
for (i in unique(breeding.data$year)){
  #subset 1 year
  breeding.data.G.SUB <- subset(breeding.data.G, year == i)
  
  #get voronoi polygons
  territories_G <- sf::st_voronoi(sf::st_union(breeding.data.G.SUB), box)
  territories_G <- sf::st_intersection(sf::st_cast(territories_G), sf::st_union(wood.outline))
  
  #joining the territory polygons back up with the individuals that bred in them
  territories_G <- sf::st_sf(geom = territories_G)
  territories_G <- sf::st_join(territories_G, breeding.data.G.SUB)
  
  #then run with st_area, it gives the area of each polygon???
  areas_G <- sf::st_area(territories_G$geometry)
  #can I then add the nest box number again?
  territories_G$area_polygon.G <- areas_G
  
  #get rid of geometry
  sf::st_geometry(territories_G) <- NULL
  
  GTIT_allyrs_area <- rbind(GTIT_allyrs_area, territories_G)
}

nrow(GTIT_allyrs_area) #17812

GTIT_allyrs_area$area_polygon.G <- as.numeric(GTIT_allyrs_area$area_polygon.G)

territory <- GTIT_allyrs_area %>%
  dplyr::select('Pnum', 'nest.box', 'area_polygon.G')

#merge with other habitat data
box_tree_terr <- merge(box_tree, territory, by.x = 'Box', by.y = 'nest.box', all.y = T)
nrow(box_tree_terr) #17812

#save
write.csv(box_tree_terr, file = file.path(dirs$data_output,'Habitat_data_nestboxes.csv'), row.names = F)



