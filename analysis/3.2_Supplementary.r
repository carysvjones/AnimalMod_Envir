box::use(R / dirs[dirs])
box::use(clean = R / clean_data)
box::use(mods = R / for_analysis)

box::use(readr[read_csv])
box::use(ggplot2[...])
box::use(dplyr)
box::use(moments)
box::use(kinship2)
box::use(sf)
box::use(tidyr)
box::use(nadiv)
box::use(magrittr[`%>%`])

output_plot_theme <- theme(
  panel.border = element_blank(),
  # make major gridlines grey
  panel.grid.major = element_line(
    colour = "#ece7e7",
    linetype = "dashed"
  ),
  # make just x and y axis white
  axis.line = element_line(colour = "#ece7e7"),
  panel.grid.minor = element_blank(),
  plot.margin = unit(c(.8, .8, .8, .8), "cm"),
  plot.background = element_rect(fill = "white"),
  panel.background = element_rect(fill = "transparent"),
  axis.text.x = element_text(size = 9), # Adjust the size as needed
  axis.text.y = element_text(
    size = 9,
    # face = "bold"
  ),
  axis.title.x = element_text(
    size = 12,
    face = "bold",
    margin = margin(t = 8)
  ),
  axis.title.y = element_text(
    size = 12,
    face = "bold",
    margin = margin(r = 10)
  ),
  legend.title = element_text(
    size = 8,
    # face = "bold",
    hjust = 0.5
  ),
  legend.text = element_text(
    size = 8,
    hjust = 0
  )
)



# BREEDING DISPERSAL ----------------------------------------------------------

# SUPPLEMENTARY

gtit_data_sub_nb <- read.csv(file.path(
  dirs$data_output,
  "GTIT_data_anim_mod.csv"
), na.strings = c("", "NA"))


# get the number of times boxes occupied
nest_box_counts <- gtit_data_sub_nb %>%
  dplyr::group_by(nest_box) %>%
  dplyr::summarise(num = length(nest_box)) %>%
  dplyr::ungroup()
# get median
median(nest_box_counts$num) # 11
# get range
range(nest_box_counts$num) # 1 - 35


# breeding dispersal - when females breed more than once
data_doubles <- gtit_data_sub_nb %>%
  dplyr::group_by(mother) %>%
  dplyr::filter(dplyr::n() > 1) %>%
  dplyr::ungroup()
nrow(data_doubles) # 6568
length(unique(data_doubles$mother)) # 2589
data_doubles_summ <- data_doubles %>%
  dplyr::group_by(mother) %>%
  dplyr::summarise(num = length(mother))
# table of number of breeding attempts
table(data_doubles_summ$num)

# distance between breeding attempts
dists <- data_doubles %>%
  dplyr::rename(year = breeding_year) %>%
  dplyr::group_by(mother) %>%
  dplyr::mutate(
    x2 = dplyr::lead(x),
    y2 = dplyr::lead(y),
    year2 = dplyr::lead(year),
    box2 = dplyr::lead(nest_box)
  ) %>%
  dplyr::mutate(dist = sqrt((x2 - x)^2 + (y2 - y)^2)) %>%
  dplyr::ungroup() %>%
  dplyr::select(mother, year, year2, nest_box, box2, dist)

# DISTANCES
median(dists$dist, na.rm = T) # 60.75m
mean(dists$dist, na.rm = T) # 93.30m
range(dists$dist, na.rm = T) # 0 - 2647.5m

# sum up the number that went 0m
nrow(subset(dists, dist == 0)) # 958

# plot
breed_disp <-
  ggplot(subset(dists, !is.na(dist)), aes(x = dist)) +
  geom_histogram(fill = "#77AADE", colour = "#2F5597") +
  geom_vline(
    xintercept = median(dists$dist, na.rm = T),
    colour = "#2F5597",
    linewidth = 0.7,
    linetype = "dashed"
  ) +
  scale_x_continuous(breaks = seq(0, 2700, by = 500)) +
  xlab("Breeding dispersal distance (m)") +
  ylab("Density") +
  output_plot_theme
breed_disp

ggsave(
  filename = "./plots/breeding_dispersal_distance.pdf",
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

env_matrix <- readRDS(file.path(
  dirs$data_output,
  "env_matrix.rds"
)) %>%
  # replace and values more than 1 to just 1
  replace(., . > 1, 1)

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
env_matrix <- round(env_matrix[1:500, 1:500], digits = 3)
spatial_matrix <- round(spatial_matrix[1:500, 1:500], digits = 3)
spatial_matrix_raw <- round(spatial_matrix_raw[1:500, 1:500], digits = 3)


# PLOTS COMPARING & MANTEL/OTHER TESTS

# Correlation matrices
mantel_matrices <- vegan::mantel(
  env_matrix,
  spatial_matrix_raw,
  method = "spearman",
  permutations = 1,
  na.rm = TRUE
)

# save
saveRDS(results, file = file.path(dirs$data_output, "mantel_output.rds"))


# Make to table -----------------------------------------------------------

# Convert the distance matrix to a data frame with row and column names
spatial_raw_df <- as.data.frame(as.table(spatial_matrix_raw))
colnames(spatial_raw_df) <- c("Mother1", "Mother2", "dist")
# remove those with themselves
spatial_raw_df <- spatial_raw_df[
  spatial_raw_df$Mother1 != spatial_raw_df$Mother2,
]
nrow(spatial_raw_df) # 58974720

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

spat_env_df <- readRDS(file.path(dirs$data_output, "spat_env_df.rds")) %>%
  # replace any values greater thaqn 1 with 1
  dplyr::mutate(env_sim = ifelse(env_sim > 1, 1, env_sim))

# are there cases where Mother1 and Mother2 are the same?
spat_env_df %>%
  dplyr::filter(Mother1 == Mother2) %>%
  nrow() # 0

max(spat_env_df$dist)

# SUPP PLOTS -------------------------------------------------------------------


# environmental similarity

median(spat_env_df$env_sim, na.rm = T)


envir_sim_hist <-
  ggplot(spat_env_df, aes(x = env_sim)) +
  geom_histogram(
    bins = 30,
    fill = "#77AADE",
    colour = "#2F5597"
    # alpha = 0.7
  ) +
  # geom_density(
  #   fill = "#77AADE",
  #   alpha = 0.8
  # ) +
  geom_vline(
    xintercept = median(spat_env_df$env_sim),
    colour = "#2F5597",
    linetype = "dashed",
    linewidth = 0.7
  ) +
  xlab("Breeding environment similarity") +
  ylab("Count") +
  output_plot_theme +
  scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6))
envir_sim_hist

ggsave(
  filename = file.path(dirs$plots, "envir_similarity_histogram.pdf"),
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

summary(lm(spat_env_df$env_sim ~ spat_env_df$dist))
cor(spat_env_df$dist, spat_env_df$env_sim)

# hex plot
comp_mat_opp <- ggplot(spat_env_df, aes(x = dist, y = env_sim)) +
  geom_hex() +
  output_plot_theme +
  scale_fill_gradient(low = "#2F5597", high = "#E7E4EF") +
  xlab("Distance (metres)") +
  ylab("Breeding environment similarity") +
  # name legend
  labs(fill = "Count") +
  theme(
    axis.title.x = element_text(
      size = 16,
      face = "bold",
      margin = margin(t = 8)
    ),
    axis.title.y = element_text(
      size = 16,
      face = "bold",
      margin = margin(r = 10)
    ),
    legend.title = element_text(
      size = 12,
      face = "bold",
      hjust = 0.5
    ),
    legend.text = element_text(
      size = 10,
      hjust = 0
    )
  )

ggsave(
  filename = file.path(dirs$plots, "compare_matrices.pdf"),
  plot = comp_mat_opp,
  scale = 1,
  width = 26,
  height = 22,
  units = "cm",
  dpi = 350,
  limitsize = FALSE,
  bg = "transparent"
)
