# DESCRIPTION ──────────────────────────────────────────────────────────────── #

# Functions to prepare data for analysis - removes mixed broods, keeps only those
#with lay dates recorded, known Pnums, nesting within the woods and fixes any 
#obvious errors in data

# DEPENDENCIES ─────────────────────────────────────────────────────────────── #

box::use(magrittr[`%>%`])
box::use(dplyr)
box::use(clean = ../ R / clean_data)


# FUNCTIONS ────────────────────────────────────────────────────────────────── #


#' Returns clean breeding data.
#'
#' @param data Wytham breeding data.
#' @return dataset clean. 

clean_breeding_data <- function(data){
  data <- data %>%
    dplyr::mutate(Lay.date = as.Date(data$Lay.date, '%d / %m / %Y'),
                  Hatch.date = as.Date(data$Hatch.date, '%d / %m / %Y'),
                  #make all ring names upper case
                  Mother = toupper(Mother),
                  Father = toupper(Father),
                  #make a column with the nest box code in
                  nest.box = as.factor(substr(Pnum, start = 6, 15))) %>%
    #factorise columns
    dplyr::mutate(across(c('Pnum', 'Section', 'Mixed.species', 'Species', 'Closed', 'Missing.entire.brood', 
                       'Suspected.predation', 'Father', 'Mother', 'Experiment.codes'), factor)) %>%
    #get rid of any mixed species
    dplyr::filter(Mixed.species != TRUE) %>%
    droplevels()
  
  #save this for later - for pedigree want to keep as many as possible                
  #data <- droplevels(subset(data, !is.na(April.lay.date)))
  
  #adjust the section of some boxes
  #some are listed as unknown, but are within the sections of the woods
  #so rename their section with actual section
  #then afterwards can remove those not in the woods
  #pnums of those with unknown section
  unknown_Pnums <- data$Pnum[data$Section == 'unknown']
  #if start with CP, EX, W, SW, O, C, MP, P, then keep, others get rid
  for(i in unknown_Pnums){
    ifelse(
      #for double letter sections
      #if is in one of the double letter sections
      substr(data$nest.box[data$Pnum == i], start = 1,2) %in% c('CP', 'EX',  'SW',  'MP'),
      #if true then paste section into Section column
      data$Section[data$Pnum == i] <- tolower(substr(data$nest.box[data$Pnum == i], start = 1,2)),
      #if not true check if its in single letter section
      ifelse(substr(data$nest.box[data$Pnum == i], start = 1,1) %in% c('W','O', 'C', 'P'),
             #if yes keep and rename section
             data$Section[data$Pnum == i] <- tolower(substr(data$nest.box[data$Pnum == i], start = 1,1)),
             #if not rename with unknown (how do i say just do nothing..?)
             data$Section[data$Pnum == i] <- 'unknown'))
  }
  
  #take out those with section unknown and he (are located outside of the woods)
  data <- data %>% 
    dplyr::filter(Section != 'unknown' & Section != 'he') %>%
  droplevels()

  #NOW FIXING SOME OTHER ERRORS IN DATA
  
  #some individuals recorded as both species
  #20101O253 - has father as blue tit mother as great tit young as blue tits - get rid 
  #20071SW75 - mother and father need to change - mother = T176160, father = NA 
  #20071O74, mother should be V260632, father should be V260844, stay as blue tit 
  data <- data %>%
    #20101O253 - has father as blue tit mother as great tit young as blue tits - get rid 
    dplyr::filter(Pnum != '20101O253' & 
                    #20082W37 listed as blue tit but is actually great tit - something about nest take over so just get rid
                    Pnum != '20082W37') %>%
    droplevels()
  
  #first have to add them as factor levels
  levels(data$Mother) <- c(levels(data$Mother), 'T176160', 'V260632')
  levels(data$Father) <- c(levels(data$Father), 'V260844')
  
  data$Mother[data$Pnum == '20071SW75'] <- 'T176160'
  data$Mother[data$Pnum == '20071O74'] <- 'V260632'
  
  data$Father[data$Pnum == '20071SW75'] <- NA
  data$Father[data$Pnum == '20071O74'] <- 'V260844'

  
  ###GREAT TIT SPECIFIC CLEANING 
  data <- data %>%
    dplyr::filter(
      #get rid of Pnum 20081EX27, 20051EX66
      #mother recorded twice in 1 year lay date same day but in different boxes...
      Pnum != '20081EX27' &
      #same here - of the 2 attempts listed I keep attempt where we know father too
      Pnum != '20051EX66' &
      #remove E985078 - is in there as male and female - 3 times as male so remove one where it is mother
      Pnum != '20081CP35' &
      #remove 20181EX20 - has same ID recorded for mother and father
      Pnum != '20181EX20' &
      #20182EX14 this pnum says no fledglings recorded, and all same fledglings are named under another pnum (20181EX14)
      Pnum != '20182EX14')
      
  #parents are listed as both mother and father?
  #pnums = 
  #TK26998, 20081SW112 - make father NA,  and 20081O75 - make mother NA 
  data$Father[data$Pnum == '20081SW112'] <- NA
  data$Mother[data$Pnum == '20081O75'] <- NA
  
  #X239739, 20101W97 - says X239739 male, and F833150 female, - make male NA
  #         20121EX39 - X239739 is female - get rid of all, because has no father recorded 
  data <- droplevels(subset(data, Pnum != '20121EX39'))
  data$Father[data$Pnum == '20101W97'] <- NA
  
  #listed as male and female - just make NA
  #Y031429 is duplicated
  #20121W56 = father , 20131W56 = mother- listed as male and female get rid - just make both NA?
  data$Father[data$Pnum == '20121W56'] <- NA
  data$Mother[data$Pnum == '20131W56'] <- NA
  
  #these are same mothers 2 breeding attempts in 1 year
  #with less than 20 days inbetween them, but both have fledglings, or first has fledglings - impossible...
  #need almost 2 weeks  to hatch and 3 weeks to fledge usually... 
  #so just cutting all out 
  too_short_int <- c('19651O86', '19821EX58', '20061EX2', '20181C32',
                     '19651O106', '19721O11', '19891W104', '20181SW109')
  data <- droplevels(subset(data, !(Pnum %in% too_short_int)))
  
  #keep 19721O13, 19881W17, 19821EX61, 19891W103A, 20061EX81A, 19971EX47A, 
  #edit 19881W38 to make mother = NA
  data$Mother[data$Pnum == '19881W38'] <- NA
  #edit 19971EX48A, make parents VR26902 and J644096 instead - in breeding data that pnum has 2 diff parents
  #but so second pnum where those parents are recorded with young
  #this pnum has 2 different parents, so change parents in breeding data to them 
  #J644096 is recorded as female at one point so make mother and other father.
  #first have to add them as factor levels
  levels(data$Mother) <- c(levels(data$Mother), 'J644096')
  levels(data$Father) <- c(levels(data$Father), 'VR26902')
  data$Mother[data$Pnum == '19971EX48A'] <- 'J644096'
  data$Father[data$Pnum == '19971EX48A'] <- 'VR26902'
  
  #20121O64 listed as great tit but is actually blue tit 
  data$Species[data$Pnum == '20121O64'] <- 'b'
  
  return(data)
}


#' Returns clean ringing data from old data (up to 2013).
#'
#' @param data Wytham ringing data.
#' @return dataset clean and standardised 

clean_ringing_data <- function(data){
  
  #make all ring names upper case
  data <- data %>%
    dplyr::mutate(#make all ring names upper case
                  bto_ring = toupper(bto_ring)) %>%
    dplyr::rename(Pnum = pnum) %>%
    #make factors
    dplyr::mutate(across(c('id', 'Pnum', 'site', 'retrap', 'age', 'sex', 
                           'bto_species_code', 'bto_ring', 'nb'), factor)) %>%
    #get rid unringed
    dplyr::filter(bto_ring != 'unringed' &
                    bto_ring != 'UNRINGED' &
                    #keep only not retraps - will be removed later when just keep age 1
                    # retrap == 'N' &
                    #keep only age 1 birds - no later!
                    # age == 1 & 
                    #keep those only in main woods
                   location == 'Wytham_Core' &
                   #remove thos with no Pnum
                   Pnum != 'NA' &
                   #just keep years from 1960 - breeding data doesn't go earlier than that 
                   yr > 1959 &
                   #remove UNRRUNT from bto_ring?
                   bto_ring != 'UNRRUNT') %>%
    droplevels()


  #REMOVE SOME GTITS
  #remove some because obvs some errors - maybe come back and work out if can keep?
  #are all duplicated assoc with different pnum each time?
  doubles <- c('VF64204', 'TK27605', 'TK27604', 'TK27508', 'TL39500')
  data <- droplevels(subset(data, !(bto_ring %in% doubles)))
  
  #remove individuals with bto_ring number 400 - can't be correct (Pnum 19981MP49)
  data <- droplevels(subset(data, bto_ring != '400'))
  
  
  #REMOVE SOME BTITS
  #remove bto_rings - 
  remove_btit_rings <- c('V729349', 'E587236', 'E906508', 'E985206', 'KJ72884', 'JN32266')
  #are parents of broods (because have sex) but listed as age 1 so treated as young 
  #so end up as own parents when merged with breeding data 
  data <- droplevels(subset(data, !(bto_ring %in% remove_btit_rings)))
  
  #2 have 2 entries
  #are all duplicated assoc with different pnum each time?
  doubles_btit <- c('F491535', 'T641298')
  #remove some because obvs some errors - maybe come back and work out if can keep?
  data <- droplevels(subset(data, !(bto_ring %in% doubles_btit)))
  
  #GTITs recorded in both early and late ringing data, but with different parents each time 
  #remove in all ringing data 
  #remove some because obvs some errors - same indiv diff parents etc - 
  #same brood id's recorded in 2013 and 2015, with different parents
  #maybe come back and work out if can keep?
  #TX42804, TX42802, TX42806, TX42803, TX42801, TX42805, TX42807, TX42808
  double_reported_twice <- c('TX42804', 'TX42802', 'TX42806', 'TX42803', 'TX42801', 'TX42805', 
                             'TX42807', 'TX42808')
  data <- droplevels(subset(data, !(bto_ring %in% double_reported_twice)))
  
  #this father is listed as age 1, but is not possible, so change to NA
  data$age[data$Pnum == '20002B159' & data$bto_ring == 'VF64204'] <- NA
  
  #get rid of pnum 19981O87 - mother is listed as age 1 and age 6 in 1 year
  #her other pnum is listed in breeding data so keeping that one,
  #just getting rid of this one 
  data <- droplevels(subset(data, Pnum != '19981O87'))
  
  return(data) 
}


#' Returns clean ringing data from new data (post-2013)
#
#' @param data Wytham breeding data.
#' @return  dataset clean and standardised.

clean_ringing_data_2 <- function(data){
  
  #make all ring names upper case
  data <- data %>%
    dplyr::mutate(
      #make all ring names upper case
      Ring = toupper(Ring),
      Date = substr(data$Date, 7, 10)) %>%
    #make factors
    dplyr::mutate(across(c('Spec', 'Sex', 'Place', 'Site', 'Date'), factor)) %>%
    #get rid of 2013 because that is in other data as well
    dplyr::filter(Date != '2013'&
                    #only not retraps
                    # Rtype == 'N' &
                    #keep age 1 only - do later!
                    # Age == 1 &
                    #from Wytham main 
                    Place == 'WYT' &
                    #get rid of those with no Ring
                    !is.na(Ring) & 
                    Ring != 'UNRINGED')
  
  #BTIT doubles in Ring number 
  #are all duplicated assoc with different pnum each time?
  doubles2_btit <- c("AVL4793", "AVL4794", "AVL4795", "AVL4796", "AVL4797", "AVL4798")
  #are one brood ringed age 1 in 2 boxes (EX73 and O257) in same year, surely must be wrong?
  #just remove all 
  data <- droplevels(subset(data, !(Ring %in% c(doubles2_btit))))
  
  #GTITs recorded in both early and late ringing data, but with different parents each time 
  #remove some because obvs some errors - same indiv diff parents etc - 
  #same brood id's recorded in 2013 and 2015, with different parents
  #maybe come back and work out if can keep?
  #TX42804, TX42802, TX42806, TX42803, TX42801, TX42805, TX42807, TX42808
  double_reported_twice <- c('TX42804', 'TX42802', 'TX42806', 'TX42803', 'TX42801', 'TX42805', 
                             'TX42807', 'TX42808')
  data <- droplevels(subset(data, !(Ring %in% double_reported_twice)))
  
  #has ring TS45000 twice , Pnum 20111SW70 (in this brood ring number fits with others) and 
  #20161SW44 (here does not fit) - remove Pnum == '20161SW44'
  data <- droplevels(subset(data, !(Date == '2016' & Ring == 'TS45000')))
  
  return(data)
  
}



#' Add age to all mothers
#' 
#' @param age age dataset.
#' @param data data to take.
#' @return dataset with age of individuals best estimated.

get_age <- function(data, ringing_data){
  
  #Get those we know when born
  #just keep greti
  ring_all_gtit <- ringing_data %>%
    dplyr::filter(bto_species_code == 'GRETI') %>%
    dplyr::mutate(yr = as.integer(yr))
  
  #just get age for those that we know for sure in ringing dataset
  #use BTO age codes
  ages <- ring_all_gtit %>%
    dplyr::mutate(Fem_DOB = 
                    dplyr::case_when(
                      age == '1' ~ yr,
                      age == '1J' ~ yr,
                      age == '3' ~ yr,
                      age == '3J' ~ yr,
                      age == '5' ~ as.integer((yr-1)),
                      age == '6' ~ as.integer((yr-2)))) %>%
    dplyr::group_by(bto_ring) %>%
    dplyr::slice_min(n = 1, order_by = Fem_DOB) %>%  # select the row with the smallest Fem_Age
    #keep only 1 row per individual 
    dplyr::slice_head() %>%
    dplyr::ungroup() %>%
    dplyr::select(bto_ring, Fem_DOB, age)
  
  #merge with breeding data 
  merged <- dplyr::left_join(data, ages, by = c('Mother' = 'bto_ring'))
  
  #get age at breeding attempt for those that know DOB
  merged <- merged %>%
    dplyr::mutate(Fem_breed_age = merged$year - merged$Fem_DOB)
  #table(merged$Fem_breed_age)
  #check and sort for any that are definitely wrong - less than 1 or more than 10, Inf (stand in for unknown) or NA
  false_ageFem <- subset(merged, Fem_breed_age < 1 | 
                           Fem_breed_age > 10 | 
                           is.na(Fem_breed_age))$Mother
  
  #various reasons why missing - some mistakes in data 
  #find these ones earliest breeding year and make their guess of DOB be 1 year before that 
  #assume that when first recorded breeding it is their 1st year actually breeding
  #because usually don't move far to breed would likely have been recorded if they did nest
  #in a box before
  if (length(false_ageFem) > 0){
    for(individual in false_ageFem){
      sub <- merged[!is.na(merged$Mother) & #subset Mother that has false age, so exclude NA mothers
                      merged$Mother == individual, ] 
      #find earliest breeding_year recorded and replace Fem_DOB with that
      merged$Fem_DOB[!is.na(merged$Mother) &
                       merged$Mother == individual] <- rep(min(sub$year) - 1, nrow(sub)) 
    }
  }
  #recalc age 
  merged <- merged %>%
    dplyr::mutate(Fem_breed_age = year - Fem_DOB)
  
  #return this dataset
  return(merged)
  
}



#' Remove second broods, for mothers, fathers, late broods - choose which to remove
#' 
#' @param data data
#' @param mothers T or F to remove known second mothers.
#' @param fathers T or F to remove known second fathers.
#' @param late T or F to remove late broods.
#' @return dataset with removed breeding attempts

remove_late_broods <- function(data, mothers, fathers, late, late_persection){
  
  if(mothers == T){
    #how many mothers have >1 in a given year 
    dupl_M <- data %>%
      dplyr::group_by(year) %>%
      dplyr::filter(clean$isdup(Mother))
    #mothers keep first broods
    dupl_M_first <- dupl_M %>%
      dplyr::group_by(., Mother, year) %>%
      dplyr::slice_min(., order_by = April.lay.date)
    #remove others
    dupl_M_toremove <- dupl_M$Pnum[!(dupl_M$Pnum %in% dupl_M_first$Pnum)]
    
    message(paste(length(unique(dupl_M_toremove))), ' breeding attempts removed as known female 2nd broods')
    
    #remove those Pnums from data
    data <- data[!(data$Pnum %in% dupl_M_toremove) , ]
  }
  
  if(fathers == T){
    #how many mothers have >1 in a given year 
    dupl_F <- data %>%
      dplyr::filter(!is.na(Father)) %>%
      dplyr::group_by(year) %>%
      dplyr::filter(clean$isdup(Father))
    #mothers keep first broods
    dupl_F_first <- dupl_F %>%
      dplyr::group_by(., Father, year) %>%
      dplyr::slice_min(., order_by = April.lay.date)
    #remove others
    dupl_F_toremove <- dupl_F$Pnum[!(dupl_F$Pnum %in% dupl_F_first$Pnum)]
    
    message(paste(length(unique(dupl_F_toremove))), ' breeding attempts removed as known male 2nd broods')
    
    #remove those Pnums from data
    data <- data[!(data$Pnum %in% dupl_F_toremove) , ]
    
  }
  
  if(late == T){
    output <- NULL
    #Find first 5% in each year then add 30
    rem <- data %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(tot = sum(April.lay.date > (sort(April.lay.date)[dplyr::n()*0.05] + 30)))

    data <- data %>%
      dplyr::group_by(year) %>%
      dplyr::filter(April.lay.date <= (sort(April.lay.date)[dplyr::n()*0.05] + 30)) %>%
      dplyr::ungroup()
    
    message(paste(sum(rem$tot), ' breeding attempts removed as late & likely 2nd broods'))
  }
  
  if(late_persection == T){
    output <- NULL
    #Find first 5% in each year then add 30
    rem <- data %>%
      dplyr::group_by(year, Section) %>%
      dplyr::mutate(threshold = sort(April.lay.date)[ceiling(dplyr::n() * 0.05)] + 30) %>%
      dplyr::filter(April.lay.date > threshold) %>%
      dplyr::summarise(tot = dplyr::n(), .groups = "drop")
    
    data <- data %>%
      dplyr::group_by(year, Section) %>%
      dplyr::mutate(threshold = sort(April.lay.date)[ceiling(dplyr::n() * 0.05)] + 30) %>%
      dplyr::filter(April.lay.date <= threshold) %>%
      dplyr::ungroup() %>%
      dplyr::select(-threshold)
    
    message(paste(sum(rem$tot), ' breeding attempts removed as late & likely 2nd broods'))
  }
    
  return(data)
  
}




