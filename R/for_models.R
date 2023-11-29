# DESCRIPTION ──────────────────────────────────────────────────────────────── #

# Functions to start running models

# DEPENDENCIES ─────────────────────────────────────────────────────────────── #

box::use(magrittr[`%>%`])
box::use(dplyr)
box::use(rptR)


# FUNCTIONS ────────────────────────────────────────────────────────────────── #


#' Returns subset of data with years you choose.
#'
#' @param data Wytham breeding data.
#' @return dataset with certain/all years.


choose_data <- function(data, years){
  ifelse(years == 'ALL', 
         output <- data, 
         output <- data %>%
           dplyr::filter(breeding_year %in% years) %>%
           droplevels())
  return(output)
}



#' Gets output of variance components for each model neatly, ready to plots
#' 
#' @param model which model type do want variance components returned for 
#' @return table, with var comp est and SE, as well as relative effect and SE

get_var_comps <- function(model) {
  mod_name <- deparse(substitute(model))
  
  var_comp_names <- switch(
    mod_name,
    LD_basic = c('Vby', 'Va', 'Vpe', 'Vr'),
    HD_basic = c('Vby', 'Va', 'Vpe', 'Vr'),
    LD_NB = c('Vby', 'Vnb', 'Va', 'Vpe', 'Vr'),
    HD_NB = c('Vby', 'Vnb', 'Va', 'Vpe', 'Vr'),
    LD_spat = c('Vby', 'Vspat', 'Va', 'Vpe', 'Vr'),
    HD_spat = c('Vby', 'Vspat', 'Va', 'Vpe', 'Vr'),
    LD_envir = c('Vby', 'Venv', 'Va', 'Vpe', 'Vr'),
    HD_envir = c('Vby', 'Venv', 'Va', 'Vpe', 'Vr'),
    stop("Invalid model specified.")
  )
  
  data <- summary(model)$varcomp %>%
    tibble::rownames_to_column(var = "name") %>%
    base::replace(., 'name', var_comp_names) %>%
    dplyr::rename('Est' = 'component', 'SE' = 'std.error') %>%
    dplyr::mutate(Rel_Est = Est / sum(Est), 
                  Rel_SE = (SE / Est) * Rel_Est) %>%
    dplyr::select(-c('z.ratio', 'bound', '%ch')) %>%
    dplyr::mutate(model_name = as.character(mod_name))
  
  return(data)
}




 #' get heritability and phenotypic variance estimates
#' 
#' @param data as output from get_var_comps function.
#' @param model which model to get estimates for.
#' @return small tibble with Vp, VP within year, h2, h2 within year


######GETS WRONG SE!! FIX THIS

# get_herit <- function(data, model){
#   
#   herits <- tibble::tibble(name = c('Vp', 'Vp_yr', 'h2', 'h2_yr'),
#                  Est = c(#Vp - total phenotypic var
#                    sum(subset(data, model_name == model)$Est, na.rm = T ),
#                    #Vp_yr, within year phenotypic var
#                    sum(subset(data, model_name == model & name != 'Vby')$Est,
#                        na.rm = T),
#                    #h2 - across years
#                    subset(data, model_name == model & name == 'Va')$Est / 
#                      sum(subset(data, model_name == model)$Est, na.rm = T),
#                    #within year h2
#                    subset(data, model_name == model & name == 'Va')$Est / 
#                      sum(subset(data, model_name == model & name != 'Vby')$Est,
#                          na.rm = T)),
#                  SE = c(#get SE for these - first Vp
#                    NA, NA,
#                    #h2 - across years
#                    subset(data, model_name == model & name == 'Va')$SE / 
#                      sum(subset(data, model_name == model)$SE, na.rm = T),
#                    #within year h2
#                    subset(data, model_name == model & name == 'Va')$SE / 
#                      sum(subset(data, model_name == model)$SE)))
#   
#   return(herits)
#   
# }


#' get heritability and phenotypic variance estimates - CORRECT HERIT USING ASREML
#' 
#' @param data as output from get_var_comps function.
#' @param model model to get estimates as asreml V4 output
#' @return small tibble with Vp, VP within year, h2, h2 within year
#' 


get_herit_asreml <- function(data, model, modelname) {
  
  filter_by_model <- subset(data, model_name == modelname)
  
  # Calculate total phenotypic variance (Vp)
  Vp <- sum(filter_by_model$Est, na.rm = TRUE)
  # Calculate within-year phenotypic variance (Vp_yr)
  Vp_yr <- sum(filter_by_model$name != 'Vby', na.rm = TRUE)
  # Calculate heritability across years (h2)
  h2 <- subset(filter_by_model, name == 'Va')[, 2] / Vp
  # Calculate heritability within year (h2_yr)
  h2_yr <- subset(filter_by_model, name == 'Va')[, 2] / Vp_yr
  
  # Choose the appropriate formula based on the number of rows
  if (nrow(filter_by_model) == 4) {
    SE_h2 <- asreml::vpredict(model, h2 ~ V2 / (V1 + V2 + V3 + V4))[, 2]
    SE_h2_yr <- asreml::vpredict(model, h2_wyr ~ V2 / (V1 + V2 + V3 + V4))[, 2]
  } else {
    SE_h2 <- asreml::vpredict(model, h2 ~ V3 / (V1 + V2 + V3 + V4 + V5))[, 2]
    SE_h2_yr <- asreml::vpredict(model, h2_wyr ~ V3 / (V2 + V3 + V4 + V5))[, 2]
  }
    
  # Create a tibble with results
  herits <- tibble::tibble(
    name = c('Vp', 'Vp_yr', 'h2', 'h2_yr'),
    Est = c(Vp, Vp_yr, h2, h2_yr),
    SE = c(NA, NA, SE_h2, SE_h2_yr)
  )
  
  return(herits)
}



#' Get bounding box for calculating territory polygons 
#' 
#' @param breeding_data breeding data with spatial object
#' @return box

bbox_polygon_G <- function(breeding_data) {
  bb <- sf::st_bbox(breeding_data)
  
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



#' Get territory polygons
#' 
#' @param breeding_data
#' @param wood_outline 
#' @param yr year of data
#' @return dataset with areas for each box occupied within that year 

get_territory_polygons <- function(breeding_data, wood_outline, yr) {
  
  #create box 
  box <- sf::st_sfc(mods$bbox_polygon_G(breeding_data))
  
  #find areas
  breeding_data_areas <- breeding_data %>%
    dplyr::filter(if_any(any_of(c("year", "breeding_year")), ~. == yr)) %>%
      #get voronoi polygons
    sf::st_union() %>%
    sf::st_voronoi(box) %>%
    sf::st_cast() %>%
    sf::st_intersection(sf::st_union(wood_outline)) %>%
    #joining the territory polygons back up with the individuals that bred in them
    sf::st_sf() %>%
    sf::st_join(dplyr::filter(breeding_data, if_any(any_of(c("year", "breeding_year")), ~. == yr))) %>%
    #get area in new column 
    dplyr::mutate(area_polygon = sf::st_area(geometry)) 
  
  return(breeding_data_areas)
  
}


