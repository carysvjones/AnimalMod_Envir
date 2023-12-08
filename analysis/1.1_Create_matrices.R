# Create matrices and compare
# Dependencies

box::use(R / dirs[dirs])
box::use(clean = R / clean_data)
box::use(mod = R / for_models)

box::use(readr[read_csv])
box::use(ggplot2[...])
box::use(dplyr[...])
box::use(moments)
box::use(kinship2)
box::use(sf)
box::use(tidyr)
box::use(magrittr[`%>%`])
box::use(Matrix)


gtit_data_sub_nb <- read.csv(
  file.path(
    dirs$data_output,
    "GTIT_data_anim_mod.csv"
  ),
  na.strings = c("", "NA")
)
nrow(gtit_data_sub_nb) # 11661

gtit_data_sub_nb <- gtit_data_sub_nb %>% # don't make breeding year a factor?
  mutate_at(
    c("breeding_year", "Mother", "Fbreeding_age_group", "nest.box"),
    as.factor
  )

# to run with a subset
# put dates as range e.g. 2011:2021, or 'ALL' to get all years of data
breed_data <- gtit_data_sub_nb
mod$choose_data(gtit_data_sub_nb, years = "ALL")
unique(breed_data$breeding_year)
nrow(breed_data) # 11661

# get distinct line for each mother with whole dataset
breed_data_uniq <- breed_data %>%
  distinct(., Mother, .keep_all = TRUE)
# filter(!is.na(alt_mean) & !is.na(No_trees_75m_mean))

nrow(breed_data_uniq) # 7680
length(unique(breed_data_uniq$Mother)) # 7680



# ENVIRONMENT SIM MATRIX  -------------------------------------------------

# As environmental measures may occur on different scales, we first scale
# and center the environmental data.
# This is done by taking away the mean (so the mean becomes 0)
# and dividing by the standard deviation for
# each environmental measure (so the variance becomes 1).
# First gets the relevant columns from the data

env_var <- breed_data_uniq[, c(
  "alt_mean", "northness_mean", "edge_mean",
  "No_trees_75m_mean", "area_polygon_sqrt"
)]

# Centers the data, and scales by standard deviation
env_var <- scale(env_var, center = TRUE, scale = T)
rownames(env_var) <- breed_data_uniq$Mother

# similarity matrix between individuals
# calculates the euclidean distance between each
# individual's environment parameters
env_var_euc <- as.matrix(dist(env_var,
  method = "euclidean",
  diag = TRUE,
  upper = T
))

# then scales so that the values are between 0 and 1
env_euc_scal <- 1 - env_var_euc / max(env_var_euc, na.rm = TRUE)
colnames(env_euc_scal) <-
  rownames(env_euc_scal) <-
  breed_data_uniq$Mother

# p <- as(solve(env_euc_scal), "dgCMatrix")
# p <- as.matrix(env_euc_scal)

# make positive definite so can run in animal model
env_euc_scal_PD <- Matrix::nearPD(env_euc_scal)
env_matrix <- env_euc_scal_PD$mat

saveRDS(env_matrix, file = file.path(dirs$data_output, "env_matrix.rds"))



# SPATIAL MATRIX ---------------------------------------------------------

# make distance matrix
# make coordinates into geometries
spatial_geom <- sf::st_as_sf(breed_data_uniq,
  coords = c("x_mean", "y_mean"),
  crs = 27700
)
# keep just geometry
geom <- spatial_geom[, c("geometry")]
rownames(geom) <- spatial_geom$Mother
# make distance matrix
spat_mat <- as.matrix(sf::st_distance(geom, geom))
colnames(spat_mat) <- rownames(spat_mat) <- spatial_geom$Mother

# can make it to a matrix and array by doing this
spat_mat_conv <- as.numeric(spat_mat)
dim(spat_mat_conv) <- c(
  nrow(spatial_geom),
  nrow(spatial_geom)
) # give it dimensions?
colnames(spat_mat_conv) <- rownames(spat_mat_conv) <- spatial_geom$Mother

# save raw distance matrix
saveRDS(spat_mat_conv, file = file.path(
  dirs$data_output,
  "spatial_matrix_raw.rds"
))

# scale the matrix - so diagonal is 1
spat_mat_conv_sc <- 1 - spat_mat_conv / max(spat_mat_conv, na.rm = TRUE)
spat_mat_conv_sc_PD <- Matrix::nearPD(spat_mat_conv_sc)
spatial_matrix <- spat_mat_conv_sc_PD$mat

# save spatial matrix
saveRDS(spatial_matrix, file = file.path(
  dirs$data_output,
  "spatial_matrix.rds"
))
