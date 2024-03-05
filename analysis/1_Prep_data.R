# Dependencies
box::use(R / dirs[dirs])
box::use(clean = R / clean_data)
box::use(mods = R / for_analysis)

box::use(readr[read_csv])
box::use(ggplot2)
box::use(dplyr)
box::use(moments)
box::use(kinship2)
box::use(sf)
box::use(tidyr)
box::use(nadiv)
box::use(magrittr[`%>%`])


# Script to prepare data for running animal models.
# First creating pedigree by merging together breeding and ringing data.
# Then taking breeding attempt data and getting local environment for each individual.


# PREPARE PEDIGREE

# MERGE BREEDING AND RINGING DATA -----------------------------------------

# input clean breeding attempt data 1960-2022
breed <- read.csv(file.path(
    dirs$data_output,
    "LT_breeding_data_greti_bluti_1960_2022.csv"
), na.strings = c("", "NA"))
nrow(breed) # 40248

# subset to just keep great tits
breed_gtit <- breed %>%
    dplyr::filter(species == "g")
nrow(breed_gtit) # 17996

# quick summaries
nrow(breed_gtit) # 17996
length(unique(breed_gtit$pnum)) # 17996
length(unique(breed_gtit$mother)) # 8531
summary(as.factor(breed_gtit$mother)) # 3830 mothers unknown
length(unique(breed_gtit$father)) # 7215
summary(as.factor(breed_gtit$father)) # 6818


# READ IN - Ringing Data  --------------------------------------------------

ring_all <- read.csv(
    file.path(
        dirs$data_output,
        "ebmp_database_ringing_record_export_GT&BT_all.csv"
    ),
    na.strings = c("", "NA")
)

# just keep great tits and keep those only age 1
ring_gtit <- ring_all %>%
    dplyr::filter(
        bto_species_code == "GRETI",
        age == 1
    )

# quick summaries
nrow(ring_gtit) # 91645
length(unique(ring_gtit$pnum)) # 12711
length(unique(ring_gtit$nb)) # 1121
length(unique(ring_gtit$bto_ring)) # 91645


# MERGE - ring and breed  ---------------------------------------

# merge both datasets by pnum -   #keep all so have all parents in pedigree
breed_ring_gtit <- merge(breed_gtit, ring_gtit, by = "pnum", all = T) %>%
    # keep only some columns
    dplyr::select(bto_ring, mother, father, sex, nest_box) %>%
    # nest box here is nest box the bto_ring bird was born in
    # rename this as use it for regression analysis (supplementary)
    dplyr::rename(nest_box_born_in = nest_box) %>%
    # remove rows where bto_ring, mother and father are NA
    dplyr::filter(!is.na(bto_ring) | !is.na(mother) | !is.na(father))

nrow(breed_gtit) # 17996
nrow(ring_gtit) # 91645

# summaries
nrow(breed_ring_gtit) # 94005
length(unique(breed_ring_gtit$bto_ring)) # 91646
length(unique(breed_ring_gtit$mother)) # 8531
length(unique(breed_ring_gtit$father)) # 7215

# save this
write.csv(breed_ring_gtit,
    file = file.path(dirs$data_output, "GTIT-pedigree-full.csv"), row.names = F
)

# CREATING PEDIGREE ----------------------------------------------------------

# rename columns
ped_gtit_short <- breed_ring_gtit %>%
    dplyr::rename(
        id = bto_ring,
        dam = mother,
        sire = father
    ) %>%
    dplyr::select(id, dam, sire, sex)
head(ped_gtit_short)
nrow(ped_gtit_short) # 94005

# fix the pedigree - adds any missing individuals to id column and orders it correctly

# find any dam or sire that aren't in id column

miss_dam <- unique(ped_gtit_short$dam[!(ped_gtit_short$dam %in%
    ped_gtit_short$id)])
miss_dam_df <- data.frame(id = miss_dam, dam = NA, sire = NA, sex = "F")

miss_sire <- unique(ped_gtit_short$sire[!(ped_gtit_short$sire %in%
    ped_gtit_short$id)])
miss_sire_df <- data.frame(id = miss_sire, dam = NA, sire = NA, sex = "M")

full_ped <- rbind(miss_dam_df, miss_sire_df, ped_gtit_short) %>%
    dplyr::filter(!is.na(id))
nrow(full_ped) # 99239


# this is the overall huge pedigree!
head(full_ped)
nrow(full_ped) # 99239

# so can save this here...
write.csv(full_ped,
    file =
        file.path(dirs$data_output, "GTIT_pedigree_full.csv"), row.names = F
)

# I myself then keep just those that have a laying date recorded, as that's what is useful for my analysis

# make column for LD with dummy values to say if have LD recorded or not - so I can keep only informative individuals
full_ped_LD <- dplyr::mutate(full_ped,
    LD_num =
        ifelse(full_ped$id %in% breed_gtit$father |
            full_ped$id %in% breed_gtit$mother, 1, NA)
)
summary(as.factor(full_ped_LD$LD_num))

# subset those with a LD
ped_gtit_final <- full_ped_LD %>%
    dplyr::filter(!is.na(LD_num)) %>%
    droplevels() %>%
    dplyr::select(id, dam, sire)
nrow(ped_gtit_final) # 15744

ped_gtit_final <- nadiv::prepPed(ped_gtit_final)

write.csv(ped_gtit_final, file = file.path(dirs$data_output, "GTIT-pedigree-for-asreml.csv"), row.names = F)



# PREPARE DATA


# GREAT TITS --------------------------------------------------------------

# BREEDING DATA ---------------------------------------------------

# how many have no mother
nrow(subset(breed_gtit_age_noNA, is.na(mother))) # 3288
nrow(subset(breed_gtit_age_noNA, is.na(father))) # 5789

# have both parents
nrow(subset(breed_gtit_age_noNA, !is.na(mother) & !is.na(father))) # 9260
# have just mother
nrow(subset(breed_gtit_age_noNA, !is.na(mother))) # 12148
nrow(subset(breed_gtit_age_noNA, !is.na(mother) & is.na(father))) # 2888 have no mother but have father
# have jsut father
nrow(subset(breed_gtit_age_noNA, !is.na(father))) # 9647
nrow(subset(breed_gtit_age_noNA, !is.na(father) & is.na(mother))) # 387 have no father but have mother


# ADD AGE  ----------------------------------------------------------------

ring <- read.csv(file.path(dirs$data_output, "ebmp_database_ringing_record_export_GT&BT_all.csv"), na.strings = c("", "NA"))

breed_gtit_age <- breed_gtit %>%
    dplyr::filter(!is.na(mother)) %>%
    clean$get_age(., ring)
table(breed_gtit_age$fem_breed_age)
#    1    2    3    4    5    6    7    8    9
# 7593 3748 1676  742  290   85   26    5    1


length(unique(breed_gtit_age$mother)) # 7791
length(unique(breed_gtit_age$father)) # 6342
length(unique(breed_gtit_age$pnum)) # 12148
summary(breed_gtit_age$april_lay_date) #-5 to 83

nrow(breed_gtit_age) # 14166

# get rid of those with experiment codes, no LD
breed_gtit_age_noNA <- breed_gtit_age %>%
    dplyr::filter(is.na(breed_gtit_age$experiment_codes) &
        !is.na(breed_gtit_age$april_lay_date)) %>%
    dplyr::select(
        "year", "pnum", "section", "april_lay_date", "april_hatch_date",
        "clutch_size", "num_fledglings", "father", "mother", "nest_box",
        "fem_dob", "fem_breed_age"
    )

nrow(breed_gtit_age_noNA) # 12148


# GET RID mother & father KNOWN 2ND BROODS & LATE BROODS ------------------

# use function - decide which ones to use
# should say how many breeding attempts it gets rid of each time
breed_gtit_age_noNA_sub <- clean$remove_late_broods(breed_gtit_age_noNA,
    mothers = T, fathers = T,
    late = F, late_persection = T
)
nrow(breed_gtit_age_noNA_sub) # 11743

# rename year as breeding year - make new column of yearling and adult
breed_gtit_age_noNA_sub <- breed_gtit_age_noNA_sub %>%
    dplyr::rename("breeding_year" = "year") %>%
    dplyr::mutate(Fbreeding_age_group = ifelse(fem_breed_age == 1,
        "Yearling", "Adult"
    ))
summary(as.factor(breed_gtit_age_noNA_sub$Fbreeding_age_group)) # Adult 5420, Yearling 6323
summary(as.factor(breed_gtit_age_noNA_sub$fem_breed_age))




# ADD SPATIAL DATA


# Just open habitat data straight away

habitat_data <- read.csv(file.path(
    dirs$data_output,
    "Habitat_data_nestboxes.csv"
), na.strings = "NA") %>%
    janitor::clean_names()


# now merge with breeding data!
gtit_data_sub_nb <- merge(breed_gtit_age_noNA_sub, habitat_data, by = "pnum")
nrow(gtit_data_sub_nb) # 11658 - only a few hundred missing from before



# SPATIAL MATRIX ----------------------------------------------------------

# FINDING MIDPOINT IN SPACE FOR INDIVIDUALS -------------------------------

# data set of all that have nest box
nrow(gtit_data_sub_nb) # 11658
length(unique(gtit_data_sub_nb$mother)) # 7679

summary(gtit_data_sub_nb$april_lay_date) #-5 to 60, median = 25
summary(gtit_data_sub_nb$april_hatch_date) # 768 NAs
summary(gtit_data_sub_nb$clutch_size) # 42 NAs



# MEAN COORDINATES OF NEST BOX FOR MULTI NESTERS
# first make line for each breeding attempt
# then after individual will have same coords for all rows/breeding attempts
# so can then jsut keep one row per individual

gtit_data_sub_nb <- gtit_data_sub_nb %>%
    dplyr::group_by(mother) %>%
    dplyr::mutate(
        x_mean = mean(x),
        y_mean = mean(y),
        alt_mean = mean(altitude_m, na.rm = T),
        northness_mean = mean(northness, na.rm = T),
        edge_mean = mean(edge_edi, na.rm = T),
        no_trees_75m_mean = mean(no_trees_75m, na.rm = T),
        area_polygon_sqrt_mean = mean(area_polygon_sqrt, na.rm = T)
    ) %>%
    dplyr::ungroup()


# output data
write.csv(gtit_data_sub_nb, file = file.path(dirs$data_output, "GTIT_data_anim_mod.csv"), row.names = F)
