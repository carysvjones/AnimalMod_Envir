box::use(.. / R / dirs[dirs])
box::use(clean = .. / R / clean_data)
box::use(mods = .. / R / for_models)

box::use(readr[read_csv])
box::use(ggplot2)
box::use(dplyr)
box::use(moments)
box::use(kinship2)
box::use(sf)
box::use(tidyr)
box::use(nadiv)
box::use(magrittr[`%>%`])




# BREEDING DISPERSAL ----------------------------------------------------------

# SUPPLEMENTARY

gtit_data_sub_nb <- read.csv(file.path(
  dirs$data_output,
  "GTIT_data_anim_mod.csv"
), na.strings = c("", "NA"))


data_doubles <- gtit_data_sub_nb %>%
  group_by(mother) %>%
  filter(n() > 1) %>%
  ungroup()
nrow(data_doubles) # 6568
data_doubles_summ <- data_doubles %>%
  group_by(mother) %>%
  summarise(num = length(mother))
# table of number of breeding attempts
table(data_doubles_summ$num)

# distance between breeding attempts
dists <- data_doubles %>%
  rename(year = breeding_year) %>%
  group_by(mother) %>%
  mutate(
    x2 = dplyr::lead(x),
    y2 = dplyr::lead(y),
    year2 = dplyr::lead(year),
    box2 = dplyr::lead(nest_box)
  ) %>%
  mutate(dist = sqrt((x2 - x)^2 + (y2 - y)^2)) %>%
  ungroup() %>%
  select(mother, year, year2, nest_box, box2, dist)

# DISTANCES
median(dists$dist, na.rm = T) # 60.75m
mean(dists$dist, na.rm = T) # 93.30m
range(dists$dist, na.rm = T) # 0 - 2647.5m

# plot
breed_disp <-
  ggplot(subset(dists, !is.na(dist)), aes(x = dist)) +
  geom_histogram(fill = "#77AADE", colour = "#2F5597") +
  geom_vline(
    xintercept = median(dists$dist, na.rm = T),
    colour = "#2F5597",
    linewidth = 1,
    linetype = "dashed"
  ) +
  scale_x_continuous(breaks = seq(0, 2700, by = 500)) +
  xlab("Breeding dispersal distance (m)") +
  ylab("Count") +
  output_plot_theme
breed_disp

ggsave(
  filename = "./plots/breeding_dispersal_distance.png",
  plot = breed_disp,
  scale = 1,
  width = 15,
  height = 10,
  units = "cm",
  dpi = 350,
  limitsize = FALSE,
  bg = "transparent"
)





# MATRICES -----------------------------------------------------------

# correlation of environmental factors
habitat_data <- read.csv(file.path(
  dirs$data_output,
  "Habitat_data_nestboxes.csv"
), na.strings = "NA") %>%
  janitor::clean_names()

# get one row for each box
habitat_data_uniq <- habitat_data %>%
  distinct(box, .keep_all = TRUE)

# scale and centre altitude_m, northness, edge_edi, no_trees_75m



env_matrix <- readRDS(file.path(
  dirs$data_output,
  "env_matrix.rds"
))

spatial_matrix <- readRDS(file.path(
  dirs$data_output,
  "spatial_matrix.rds"
))
spatial_matrix_raw <- readRDS(file.path(
  dirs$data_output,
  "spatial_matrix_raw.rds"
))


env_matrix[1:10, 1:10]
spatial_matrix[1:10, 1:10]
spatial_matrix_raw[1:10, 1:10]

# subset to test
env_matrix <- round(env_matrix[1:5000, 1:5000], digits = 3)
spat_sub <- round(spatial_matrix[1:10, 1:10], digits = 3)
spat_matrix_raw <- round(spatial_matrix_raw[1:5000, 1:5000], digits = 3)


# PLOTS COMPARING & MANTEL/OTHER TESTS



# Correlation matrices
mantel_matrices <- vegan::mantel(
  env_matrix,
  spatial_matrix,
  method = "spearman",
  permutations = 1,
  na.rm = TRUE
)



# Make to table -----------------------------------------------------------

# Convert the distance matrix to a data frame with row and column names
spatial_raw_df <- as.data.frame(as.table(spatial_matrix_raw))
colnames(spatial_raw_df) <- c("Mother1", "Mother2", "dist")
# remove those with themselves
spatial_raw_df <- spatial_raw_df[
  spatial_raw_df$Mother1 != spatial_raw_df$Mother2,
]
nrow(spatial_raw_df) # 58974720

# 1st try to convert env to df
# regular_matrix <- as.matrix(env_matrix)
# envir_scal_df <- as.data.frame(as.table(regular_matrix))
# colnames(envir_scal_df) <- c("Mother1", "Mother2", "dist")


# For environ matrix - is saved as a different type to save space
# Get row and column names
names_vector <- rownames(env_matrix)
# Create a data frame with three columns: name1, name2, distance
env_df <- expand.grid(Mother1 = names_vector, Mother2 = names_vector)
env_df$env_sim <- as.vector(env_matrix)
# Ensure that name1 and name2 are not the same
env_df <- env_df[env_df$Mother1 != env_df$Mother2, ]
nrow(env_df) # 58974720


head(env_df)
head(spatial_raw_df)

# merge
spat_env_df <- spatial_raw_df %>%
  select(dist) %>%
  cbind(., env_df)
head(spat_env_df)
nrow(spat_env_df) # 58974720

# save
saveRDS(spat_env_df, file = file.path(dirs$data_output, "spat_env_df.rds"))



# SUPP PLOTS -------------------------------------------------------------------

# environmental similarity
envir_sim_hist <-
  ggplot(spat_env_df, aes(x = env_sim)) +
  geom_histogram(aes(y = ..density..),
    bins = 30,
    fill = "#77AADE",
    colour = "#2F5597",
    alpha = 0.7
  ) +
  geom_density(
    fill = "#77AADE",
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = median(spat_env_df$env_sim),
    colour = "#2F5597",
    linetype = "dashed",
    linewidth = 0.7
  ) +
  xlab("Environmental similarity measure") +
  ylab("Count") +
  output_plot_theme


ggsave(
  filename = file.path(dirs$plots, "envir_similarity_histogram.png"),
  plot = envir_sim_hist,
  scale = 1,
  width = 15,
  height = 10,
  units = "cm",
  dpi = 350,
  limitsize = FALSE,
  bg = "transparent"
)



# Compare with distance matrix --------------------------------------------


summary(lm(sim_both$envir ~ sim_both$dist_raw))
cor(sim_both$dist_raw, sim_both$envir)


## density plot
ggplot(spat_env_df, aes(x = dist, y = env_sim)) +
  geom_bin2d() +
  scale_fill_continuous(type = "viridis") +
  output_plot_theme



# for POSTER
comp_mat <- ggplot(spat_env_df, aes(x = dist, y = env_sim)) +
  geom_hex() +
  scale_fill_gradient(low = "#E7E4EF", high = "#4945A0") +
  theme(
    panel.background = element_rect(fill = "#B4C7E7", color = "#B4C7E7"),
    plot.background = element_rect(fill = "#B4C7E7", color = "#B4C7E7"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    # axis.title.x = element_text(colour = "#a4948c"),
    # axis.title.y = element_text(colour = "#a4948c"),
    text = element_text(color = "#2F5597", face = "bold", size = 14),
    axis.text = element_text(color = "#2F5597", face = "bold", size = 10),
    legend.position = "none"
  )


ggsave(
  filename = file.path(dirs$plots, "compare_matrices.png"),
  plot = comp_mat,
  scale = 1,
  width = 26,
  height = 22,
  units = "cm",
  dpi = 350,
  limitsize = FALSE,
  bg = "transparent"
)

# poster plot - opposite colours
comp_mat_opp <- ggplot(spat_env_df, aes(x = dist, y = env_sim)) +
  geom_hex() +
  output_plot_theme +
  scale_fill_gradient(low = "#4945A0", high = "#E7E4EF")

ggsave(
  filename = file.path(dirs$plots, "compare_matrices.png"),
  plot = comp_mat_opp,
  scale = 1,
  width = 26,
  height = 22,
  units = "cm",
  dpi = 350,
  limitsize = FALSE,
  bg = "transparent"
)





# what about make most similar 0
nest_both$dist2 <- 1 - nest_both$dist
nest_both$envir2 <- 1 - nest_both$envir

comp_mat_zero <- ggplot(nest_both, aes(x = dist2, y = envir2)) +
  geom_hex() +
  scale_fill_gradient(low = "#4945A0", high = "#E7E4EF") +
  xlab("Spatial") +
  ylab("Environment") +
  theme(
    panel.background =
      element_rect(fill = "transparent"), # transparent panel bg
    plot.background =
      element_rect(fill = "transparent", color = NA), # transparent plot bg
    panel.grid.major = element_blank(), # remove major gridlines
    panel.grid.minor = element_blank(), # remove minor gridlines
    legend.background =
      element_rect(fill = "transparent"), # transparent legend bg
    legend.box.background =
      element_rect(fill = "transparent"), # transparent legend panel
    # panel.background = element_rect(fill = "#B4C7E7", color = "#B4C7E7"),
    # plot.background = element_rect(fill = "#B4C7E7", color = "#B4C7E7"),
    # panel.grid.major.x = element_blank(),
    # panel.grid.minor.x = element_blank(),
    # panel.grid.minor.y = element_blank(),
    # panel.grid.major.y = element_blank(),
    # axis.title.x = element_text(colour = "#a4948c"),
    # axis.title.y = element_text(colour = "#a4948c"),
    text = element_text(color = "#2F5597", face = "bold", size = 34),
    axis.text = element_text(color = "#2F5597", face = "bold", size = 24)
  )
# legend.position = 'none')


ggsave(
  filename = "./plots/comp_mat_zero_legend.png",
  plot = comp_mat_zero,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 350,
  limitsize = FALSE,
  bg = "transparent"
)


# plot with raw distance and environment similarity with 1 most similar
raw_plot <- ggplot(nest_both, aes(x = dist_raw, y = envir_raw)) +
  geom_hex() +
  scale_fill_gradient(low = "#4945A0", high = "#E7E4EF") +
  xlab("Distance (m)") +
  ylab("Environmental Similarity") +
  theme(
    panel.background =
      element_rect(fill = "transparent"), # transparent panel bg
    plot.background =
      element_rect(fill = "transparent", color = NA), # transparent plot bg
    panel.grid.major = element_blank(), # remove major gridlines
    panel.grid.minor = element_blank(), # remove minor gridlines
    legend.background =
      element_rect(fill = "transparent"), # transparent legend bg
    legend.box.background =
      element_rect(fill = "transparent"), # transparent legend panel
    # panel.background = element_rect(fill = "#B4C7E7", color = "#B4C7E7"),
    # plot.background = element_rect(fill = "#B4C7E7", color = "#B4C7E7"),
    # panel.grid.major.x = element_blank(),
    # panel.grid.minor.x = element_blank(),
    # panel.grid.minor.y = element_blank(),
    # panel.grid.major.y = element_blank(),
    # axis.title.x = element_text(colour = "#a4948c"),
    # axis.title.y = element_text(colour = "#a4948c"),
    text = element_text(color = "#2F5597", face = "bold", size = 34),
    axis.text = element_text(color = "#2F5597", face = "bold", size = 24),
    legend.position = "none"
  )
raw_plot

ggsave(
  filename = "./plots/comp_mat_raw_plot.png",
  plot = raw_plot,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 350,
  limitsize = FALSE,
  bg = "transparent"
)
