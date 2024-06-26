# Dependencies
box::use(R / dirs[dirs])
box::use(clean = R / clean_data)
box::use(mod = R / for_analysis)

box::use(readr[read_csv])
box::use(ggplot2[...])
box::use(dplyr[...])
box::use(moments)
box::use(kinship2)
box::use(sf)
box::use(tidyr)
box::use(magrittr[`%>%`])
box::use(nadiv)
box::use(tibble)
box::use(ggdist[...])
box::use(patchwork[...])
box::use(distributional[...])


output_plot_theme <- theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(.8, .8, .8, .8), "cm"),
    plot.background = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    axis.text.x = element_text(size = 9), # Adjust the size as needed
    axis.text.y = element_text(
        size = 8,
        face = "bold"
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
        face = "bold",
        hjust = 0.5
    ),
    legend.text = element_text(
        size = 8,
        hjust = 0
    )
)

# asreml requires a licence and install
# library(asreml)
# install.packages("/Users/user/Downloads/asreml_4.2.0.276.tgz",
# repos = NULL, type = "source")



# Read in model output -------------------------------------------------------

# Results LD --------------------------------------------------------------

LD_basic <- readRDS(file.path(dirs$model_output, "LD_basic_output.rds"))
LD_NB <- readRDS(file.path(dirs$model_output, "LD_NB_output.rds"))
LD_spat <- readRDS(file.path(dirs$model_output, "LD_spat_output.rds"))
LD_envir <- readRDS(file.path(dirs$model_output, "LD_envir_output.rds"))

HD_basic <- readRDS(file.path(dirs$model_output, "HD_basic_output.rds"))
HD_NB <- readRDS(file.path(dirs$model_output, "HD_NB_output.rds"))
HD_spat <- readRDS(file.path(dirs$model_output, "HD_dist_output.rds"))
HD_envir <- readRDS(file.path(dirs$model_output, "HD_envir_output.rds"))



# basic model outputs
summary(LD_basic)$varcomp
summary(LD_NB)$varcomp
summary(LD_spat)$varcomp
summary(LD_envir)$varcomp

summary(HD_basic)$varcomp
summary(HD_NB)$varcomp
summary(HD_spat)$varcomp
summary(HD_envir)$varcomp




# LAYING DATE OUTPUT  ------------------------------------------------

# laying date results

# variance componenets
LD_varcomps <- rbind(
    mod$get_var_comps(model = LD_basic),
    mod$get_var_comps(model = LD_NB),
    mod$get_var_comps(model = LD_spat),
    mod$get_var_comps(model = LD_envir)
)

LD_varcomps

# heritability
# make new little table - Vp, Vp_year, h2, h2_year

LD_herits <- rbind(
    mod$get_herit_asreml(LD_varcomps, LD_basic, "LD_basic"),
    mod$get_herit_asreml(LD_varcomps, LD_NB, "LD_NB"),
    mod$get_herit_asreml(LD_varcomps, LD_spat, "LD_spat"),
    mod$get_herit_asreml(LD_varcomps, LD_envir, "LD_envir")
) %>%
    dplyr::mutate(
        Rel_Est = NA,
        Rel_SE = NA,
        model_name = c(
            rep("LD_basic", 2),
            rep("LD_NB", 2),
            rep("LD_spat", 2),
            rep("LD_envir", 2)
        )
    )


head(LD_varcomps)
head(LD_herits)

# bind together and to 3 dp - combine all using dplyr
LD_all <- rbind(LD_varcomps, LD_herits) %>%
    dplyr::mutate_at(
        vars(Est, SE, Rel_Est, Rel_SE),
        ~ if_else(is.na(.), "-",
            format(round(., 3), nsmall = 3)
        )
    )
# order them by model
order <- c("LD_basic", "LD_NB", "LD_spat", "LD_env")
# merge with brackets
LD_all_merge <- LD_all %>%
    dplyr::arrange(factor(model_name, levels = order)) %>%
    dplyr::mutate(
        Est = ifelse(SE != "-", paste(Est, " (", SE, ")", sep = ""), Est),
        Rel_Est = ifelse(Rel_SE != "-", paste(Rel_Est, " (", Rel_SE, ")", sep = ""), Rel_Est)
    ) %>%
    dplyr::select(-c(SE, Rel_SE))
# this is table use in results!

write.csv(LD_all_merge, file = file.path(dirs$model_output, "LD_model_output.csv"), row.names = F)


# HATCHING DATE OUTPUT ------------------------------------------------
# variance componenets
HD_varcomps <- rbind(
    mod$get_var_comps(model = HD_basic),
    mod$get_var_comps(model = HD_NB),
    mod$get_var_comps(model = HD_spat),
    mod$get_var_comps(model = HD_envir)
)

HD_varcomps


# heritability
# make new little table - Vp, Vp_year, h2, h2_year
HD_herits <- rbind(
    mod$get_herit_asreml(HD_varcomps, HD_basic, "HD_basic"),
    mod$get_herit_asreml(HD_varcomps, HD_NB, "HD_NB"),
    mod$get_herit_asreml(HD_varcomps, HD_spat, "HD_spat"),
    mod$get_herit_asreml(HD_varcomps, HD_envir, "HD_envir")
) %>%
    dplyr::mutate(
        Rel_Est = NA,
        Rel_SE = NA,
        model_name = c(
            rep("HD_basic", 2),
            rep("HD_NB", 2),
            rep("HD_spat", 2),
            rep("HD_envir", 2)
        )
    )


head(HD_varcomps)
head(HD_herits)

# bind together and to 3 dp - combine all using dplyr

HD_all <- rbind(HD_varcomps, HD_herits) %>%
    mutate_at(
        vars(Est, SE, Rel_Est, Rel_SE),
        ~ if_else(is.na(.), "-",
            format(round(., 3), nsmall = 3)
        )
    )
# order them by model
order <- c("HD_basic", "HD_NB", "HD_spat", "HD_env")
# merge with brackets
HD_all_merge <- HD_all %>%
    dplyr::arrange(factor(model_name, levels = order)) %>%
    mutate(
        Est = ifelse(SE != "-", paste(Est, " (", SE, ")", sep = ""), Est),
        Rel_Est = ifelse(Rel_SE != "-", paste(Rel_Est, " (", Rel_SE, ")", sep = ""), Rel_Est)
    ) %>%
    select(-c(SE, Rel_SE))
# this is table use in results!

write.csv(HD_all_merge, file = file.path(dirs$model_output, "HD_model_output.csv"), row.names = F)




# FIGURES ---------------------------------------------------------------

# Name var comps in order
LD_varcomps$name_plot <- factor(LD_varcomps$name, levels = c("Vr", "Vby", "Vnb", "Vspat", "Venv", "Vpe", "Va"))

LD_varcomps$model_name_plot <-
    factor(LD_varcomps$model_name, levels = c("LD_basic", "LD_NB", "LD_spat", "LD_envir"))


# Set labels
custom_model_labs <- gsub(" ", "\n", c("Basic", "Nestbox", "Spatial Prox.", "Breeding Env.Sim."))
custom_legend_labels <- c(
    expression(V[R]),
    expression(V[BY]),
    expression(V[NB]),
    expression(V[Spatial]),
    expression(V[BreedEnv]),
    expression(V[PE]),
    expression(V[A])
)

LD_var_comp <-
    ggplot(
        LD_varcomps,
        aes(
            x = model_name_plot,
            y = Rel_Est,
            fill = name_plot
        )
    ) +
    geom_col(position = "stack", colour = "black") +
    ylab("Proportion of variance") +
    xlab("Laying Date") +
    labs(fill = "Variance \n component") +
    scale_fill_manual(
        values = c(
            "#CFCFC4", "#43BB99", "#EE8866", "#FEAABB",
            "#99DDFE", "#EEDD88", "#77AADE"
        ),
        labels = custom_legend_labels
    ) +
    scale_x_discrete(labels = custom_model_labs) +
    guides(fill = guide_legend(override.aes = list(color = NA))) +
    output_plot_theme +
    theme(
        axis.ticks = element_blank(),
        axis.ticks.y = element_line(colour = "grey"),
        legend.key = element_rect(color = "transparent"),
        legend.box.margin = margin(2.5, 2.5, 2.5, 2.5),
        legend.box.background = element_rect(color = "black", linewidth = 0.5)
    )

LD_var_comp

ggsave(
    filename = file.path(dirs$plots, "LD_var_comp.png"),
    plot = LD_var_comp,
    scale = 1,
    width = 14,
    height = 18,
    units = "cm",
    dpi = 350,
    limitsize = FALSE,
    bg = "transparent"
)


LD_herits_sub <- LD_herits %>%
    dplyr::filter(name == "h2_yr")

LD_herits_sub$model_name_plot <- c("Basic", "Nestbox", "Spatial", "Breeding env.")
LD_herits_sub$model_name_plot <- factor(LD_herits_sub$model_name_plot,
    levels = c("Breeding env.", "Spatial", "Nestbox", "Basic")
)

LD_herit_plot <-
    ggplot(
        LD_herits_sub,
        aes(
            x = model_name_plot,
            y = Est,
            ydist = distributional::dist_normal(Est, SE)
        )
    ) +
    geom_hline(
        yintercept = c(0.0875, 0.1125, 0.1375, 0.1625, 0.1875, 0.2125, 0.2375),
        linetype = "dashed",
        color = "light grey",
        linewidth = 0.15
    ) +
    geom_errorbar(
        aes(
            ymin = Est - SE,
            ymax = Est + SE
        ),
        width = 0.5,
        colour = "#2F5597"
    ) +
    geom_point(
        colour = "#2F5597",
        size = 5,
        shape = 23,
        fill = "#77AADE"
    ) +
    # stat_eye(colour = '#2F5597',
    #          fill = '#77AADE',
    #          alpha = 0.7,
    #           ) +
    scale_y_continuous(limits = c(0.08, 0.24)) +
    scale_size_area(max_size = 1) +
    coord_flip() +
    ylab("Laying date within-year heritability") +
    output_plot_theme +
    theme(
        panel.grid.major.x = element_line(
            color = "grey",
            linewidth = 0.5
        ),
        panel.grid.minor.x = element_line(
            color = "grey",
            linewidth = 0.25,
            linetype = "dashed"
        ),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "grey", fill = NA, linewidth = 1),
        panel.background = element_rect(fill = "transparent", color = NA)
    )
LD_herit_plot

ggsave(
    filename = file.path(dirs$plots, "LD_herit.png"),
    plot = LD_herit_plot,
    scale = 1,
    width = 30,
    height = 8,
    units = "cm",
    dpi = 350,
    limitsize = FALSE,
    bg = "transparent"
)


# Name var comps in order
HD_varcomps$name_plot <- factor(HD_varcomps$name, levels = c("Vr", "Vby", "Vnb", "Vspat", "Venv", "Vpe", "Va"))

HD_varcomps$model_name_plot <-
    factor(HD_varcomps$model_name, levels = c("HD_basic", "HD_NB", "HD_spat", "HD_envir"))


HD_var_comp <-
    ggplot(
        HD_varcomps,
        aes(
            x = model_name_plot,
            y = Rel_Est,
            fill = name_plot
        )
    ) +
    geom_col(position = "stack", colour = "black") +
    ylab("Proportion of variance") +
    xlab("Hatching Date") +
    labs(fill = "Variance \n component") +
    scale_fill_manual(
        values = c(
            "#CFCFC4", "#43BB99", "#EE8866", "#FEAABB",
            "#99DDFE", "#EEDD88", "#77AADE"
        ),
        labels = custom_legend_labels
    ) +
    scale_x_discrete(labels = custom_model_labs) +
    guides(fill = guide_legend(override.aes = list(color = NA))) +
    output_plot_theme +
    theme(
        axis.ticks = element_blank(),
        axis.ticks.y = element_line(colour = "grey"),
        legend.key = element_rect(color = "transparent"),
        legend.box.margin = margin(2.5, 2.5, 2.5, 2.5),
        legend.box.background = element_rect(color = "black", linewidth = 0.5)
    )

HD_var_comp

ggsave(
    filename = file.path(dirs$plots, "HD_var_comp.png"),
    plot = HD_var_comp,
    scale = 1,
    width = 14,
    height = 18,
    units = "cm",
    dpi = 350,
    limitsize = FALSE,
    bg = "transparent"
)


HD_herits_sub <- HD_herits %>%
    dplyr::filter(name == "h2_yr")

HD_herits_sub$model_name_plot <- c(
    "Basic", "Nestbox", "Spatial",
    "Breeding env."
)
HD_herits_sub$model_name_plot <- factor(HD_herits_sub$model_name_plot,
    levels = c(
        "Breeding env.", "Spatial",
        "Nestbox", "Basic"
    )
)

HD_herit_plot <-
    ggplot(
        HD_herits_sub,
        aes(
            x = model_name_plot,
            y = Est,
            ydist = distributional::dist_normal(Est, SE)
        )
    ) +
    geom_hline(
        yintercept = c(
            0.0875, 0.1125, 0.1375, 0.1625, 0.1875,
            0.2125, 0.2375
        ),
        linetype = "dashed",
        color = "light grey",
        linewidth = 0.15
    ) +
    geom_errorbar(
        aes(
            ymin = Est - SE,
            ymax = Est + SE
        ),
        width = 0.5,
        colour = "#2F5597"
    ) +
    geom_point(
        colour = "#2F5597",
        size = 5,
        shape = 23,
        fill = "#77AADE"
    ) +
    # stat_eye(colour = '#2F5597',
    #          fill = '#77AADE',
    #          alpha = 0.7,
    #           ) +
    # make x axis begin at 0.08
    scale_y_continuous(limits = c(0.08, 0.24)) +
    scale_size_area(max_size = 1) +
    coord_flip() +
    ylab("Hatching date within-year heritability") +
    output_plot_theme +
    theme(
        panel.grid.major.x = element_line(
            color = "grey",
            linewidth = 0.5
        ),
        panel.grid.minor.x = element_line(
            color = "grey",
            linewidth = 0.25,
            linetype = "dashed"
        ),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "grey", fill = NA, linewidth = 1),
        # make whole plot outline transparent
        panel.background = element_rect(fill = "transparent", color = NA)
    )
HD_herit_plot


ggsave(
    filename = file.path(dirs$plots, "HD_herit.png"),
    plot = HD_herit_plot,
    scale = 1,
    width = 30,
    height = 8,
    units = "cm",
    dpi = 350,
    limitsize = FALSE
)




# Compare models -------------------------------------------------------------

# LAY DATE

# Basic vs nest box
2 * (LD_NB$loglik - LD_basic$loglik)
1 - pchisq(2 * (LD_NB$loglik - LD_basic$loglik), 1)

# nest box vs spatial
2 * (LD_spat$loglik - LD_NB$loglik)
1 - pchisq(2 * (LD_spat$loglik - LD_NB$loglik), 1)

# nest box vs environ
2 * (LD_envir$loglik - LD_NB$loglik)
1 - pchisq(2 * (LD_envir$loglik - LD_NB$loglik), 1)

# spatial vs environ
2 * (LD_spat$loglik - LD_envir$loglik)
1 - pchisq(2 * (LD_spat$loglik - LD_envir$loglik), 1)



# HATCH DATE

# Basic vs nest box
2 * (HD_NB$loglik - HD_basic$loglik)
1 - pchisq(2 * (HD_NB$loglik - HD_basic$loglik), 1)

# nest box vs spatial
2 * (HD_spat$loglik - HD_NB$loglik)
1 - pchisq(2 * (HD_spat$loglik - HD_NB$loglik), 1)

# nest box vs environ
2 * (HD_envir$loglik - HD_NB$loglik)
1 - pchisq(2 * (HD_envir$loglik - HD_NB$loglik), 1)

# spatial vs environ
2 * (HD_spat$loglik - HD_envir$loglik)
1 - pchisq(2 * (HD_spat$loglik - HD_envir$loglik), 1)
