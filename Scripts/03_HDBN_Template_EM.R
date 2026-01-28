#############################################################################################
#                                                                                           #
#                                                                                           #
#############           EM learning for data driven CPT Generation                  #########
#                                                                                           #
#                                                                                           #
#############################################################################################


# Read the template HDBN generated from GeNIE
# Try EM learning with biological and environmental data state tables
# Template nodes can be given "learning tables" with rows representing one case of observed evidence 

# INIT -----

  # Clean slate
  rm(list=ls(all=TRUE))
  
  # Libraries initialization
  init1 <- tail(unlist(strsplit(rstudioapi::getActiveDocumentContext()$path, "/")), n = 1)
  source("Scripts/00_Initialisation.R")
  
  # /!\ Exception from the Initialisation.R script
    # Install the package rSMILE:
    # you'll need to download the package manually from BayesFusion:
    # https://download.bayesfusion.com/files.html?category=Academia
    # Copy/paste the .tar.gz file at the root of the Rproj
    # Then install it manually
    # DONE ONCE, THIS PART SHOULD STAY COMMENTED IF THE PACKAGE IS ALREADY INSTALLED
    # install.packages("rSMILE_2.4.0_R_x86_64-pc-linux-gnu.tar.gz", repos = NULL, type = "source")
    library(rSMILE)
    source("License.R")
  
  # eo init ----


  
# CORE MODEL CPTs ----
  
  # Load the model
  net <- Network()
  net$readFile(file.path(pathGra,"Hierarchical_Dynamic_Bayesian_Network_CoralReefs_RORC_Template.xdsl"))    
  
  # Listing nodes
  # getCPT(net, "Fish_Carnivores_Biomass_Station_OBS_XX")
  # getCPT(net, "Invertebrate_Abundance_Station_True_XX")
  
  nodes <- net$getAllNodeIds()
  length(nodes)
  
  # getNodeSpec(net, "Branching_Coral_Cover_Station_TREND")
  # getNodeSpec(net, "Structural_Complexity", temporal_order = 1:3)
  # getNodeSpec(net, "Fish_Diversity", temporal_order = 1:3)
  # getNodeSpec(net, "Sand_and_Silt_Cover_Station_True_XX", temporal_order = 1:3)
  
  lapply(nodes,  function(x) getNodeSpec(net, x))
  
  # getNodeSpec(net, "Live_Coral_Cover")
  
  # eo core model cpts ----
  
  
  
# READ DATA ----
  
  nodes <- net$getAllNodeIds()
  
  # Function to select only relevant data
  keep_evidence_cols <- function(df, nodeIds = nodes) {
    df %>%
      # select(Year, Site, Station, where(~ !is.list(.))) %>%
      # Keep any crossname if present
      # select(any_of(unique(c(names(crossNames), names(crossNames_trend)))))
      # Select by present node in network
      select(any_of(c("Year","Site","Station", nodeIds)))
  }
  
  # Read all data and remap all colnames to fit the hdbn node ids
  # Cross names index for OBS
  crossNames <- c(
    Year = "year",
    Year = "SEASON",
    Year = "Campagne",
    Station = "station",
    Site = "site",
    
    Branching_Coral_Cover_Station_OBS = "RorcHCB",
    Massive_Coral_Cover_Station_OBS = "RorcHCM",
    Other_Forms_Coral_Cover_Station_OBS = "RorcHCO",
    Tabular_Coral_Cover_Station_OBS = "RorcHCT",
    Soft_Coral_Cover_Station_OBS = "RorcSC",
    Algal_Cover_Station_OBS = "RorcFS",
    Sand_and_Silt_Cover_Station_OBS = "RorcSD_SI",
    Debris_Cover_Station_OBS = "RorcRB",
    Rock_and_Slabs_Cover_Station_OBS = "RorcRC",
    Coral_Diversity_Station_OBS = "RorcCoralRichness",
    
    Fish_Diversity_Station_OBS = "RorcFishRichness",
    Fish_Carnivores_Biomass_Station_OBS = "RorcCarnivoreAbund",
    Fish_Herbivores_Biomass_Station_OBS = "RorcHerbivoreAbund",
    Fish_Corallivores_Biomass_Station_OBS = "RorcCorallivoreAbund",
    Fish_Scraper_Biomass_Station_OBS = "RorcScraperAbund",
    
    Invertebrate_Diversity_Station_OBS = "RorcInvRichness",
    Invertebrate_Abundance_Station_OBS = "RorcInvAbund",
    Sea_Urchins_Abundance_Station_OBS = "RorcUrchinAbund"
  )
  
  # For trend
  crossNames_trend <- crossNames
  names(crossNames_trend) <- sub("_OBS$", "_TREND", names(crossNames_trend))
  
  # General and Trend States
  # Coral
  coralG <- read.csv(file.path(pathProSpe,"RORC_Coral_Station_General_States_hdbn.csv"))
  coralG <- coralG %>%
    rename(any_of(crossNames))
  coralG <- keep_evidence_cols(coralG)
  head(coralG)
  
  coralT <- read.csv(file.path(pathProSpe,"RORC_Coral_Station_Trends_States_hdbn.csv"))
  coralT <- coralT %>%
    rename(any_of(crossNames_trend))
  coralT <- keep_evidence_cols(coralT)
  head(coralT)
  
  # Fish
  fishG <- read.csv(file.path(pathProSpe,"RORC_Fish_Station_General_States_hdbn.csv"))
  fishG <- fishG %>%
    rename(any_of(crossNames))
  fishG <- keep_evidence_cols(fishG)
  head(fishG)
  
  fishT <- read.csv(file.path(pathProSpe,"RORC_Fish_Station_Trends_States_hdbn.csv"))
  fishT <- fishT %>%
    rename(any_of(crossNames_trend))
  fishT <- keep_evidence_cols(fishT)
  head(fishT)
  
  # Invertebrates
  invG <- read.csv(file.path(pathProSpe,"RORC_Inv_Station_General_States_hdbn.csv"))
  invG <- invG %>%
    rename(any_of(crossNames))
  invG <- keep_evidence_cols(invG)
  head(invG)
  
  invT <- read.csv(file.path(pathProSpe,"RORC_Inv_Station_Trends_States_hdbn.csv"))
  invT <- invT %>%
    rename(any_of(crossNames_trend))
  invT <- keep_evidence_cols(invT)
  head(invT)
  
  # Environment
    # Oceanic Nino Index
    ninoa <- read.csv(file.path(pathProEnv,"OceanicNinoIndex","Env_Nino_Phase_General_ONI.csv"))
    colnames(ninoa)[grep("State|state",colnames(ninoa))] <- "Env_Nino_Phase_General"
    ninoa <- keep_evidence_cols(ninoa)
    head(ninoa)
    
    # Geomorphology
    geo <- read.csv(file.path(pathProEnv,"Geomorphology","Env_Geomorphology_Station_States_New_Caledonia.csv"))
    colnames(geo)[grep("State|state|Geo",colnames(geo))] <- "Env_Geomorphology"
    geo <- keep_evidence_cols(geo)
    head(geo)
    
    # MPA
    mpa <- read.csv(file.path(pathProEnv,"MPA","Env_MPA_Station_States_New_Caledonia.csv"))
    colnames(mpa)[grep("State|state",colnames(mpa))] <- "Env_MPA"
    mpa <- keep_evidence_cols(mpa)
    head(mpa)
    
    # Gravity
    grav <- read.csv(file.path(pathProEnv,"Gravity","Env_Gravity_Station_States_New_Caledonia.csv"))
    colnames(grav)[grep("State|state|States|states",colnames(grav))] <- "Env_Gravity"
    grav <- keep_evidence_cols(grav)
    head(grav)
    
    # Cyclone season
    cyc <- read.csv(file.path(pathProEnv,"Cyclones","Env_Cyclone_Frequency_General_States_New_Caledonia.csv"))
    cyc <- cyc %>%
      rename(any_of(crossNames_trend))
    colnames(cyc)[grep("State|state|States|states",colnames(cyc))] <- "Env_Cyclone_Frequency_General"
    cyc <- keep_evidence_cols(cyc)
    head(cyc)
    
    # Temperature Site
    temp <- read.csv(file.path(pathProEnv,"Temperature","RORC_Env_Temperature_Regime_Site_States_hdbn.csv"))
    temp <- temp %>%
      rename(any_of(crossNames_trend))
    colnames(temp)[grep("State|state|States|states",colnames(temp))] <- "Env_Site_Temperature_Regime"
    temp <- keep_evidence_cols(temp)
    head(temp) 
    
    # Chla Site
    chla <- read.csv(file.path(pathProEnv,"Chlorophyll_a","RORC_Env_Chlorophyll-a_Regime_Site_States_hdbn.csv"))
    chla <- chla %>%
      rename(any_of(crossNames_trend))
    colnames(chla)[grep("State|state|States|states",colnames(chla))] <- "Env_Site_Chlorophyll_a_Regime"
    chla <- keep_evidence_cols(chla)
    head(chla)
    
    # Heatwaves BAA
    baa <- read.csv(file.path(pathProEnv,"Heatwaves_BAA","Env_Bleaching_Alert_Area_Station_States_New_Caledonia.csv"))
    baa <- baa %>%
      rename(any_of(crossNames_trend))
    colnames(baa)[grep("State|state|States|states",colnames(baa))] <- "Env_Bleaching_Alert_Area"
    baa <- keep_evidence_cols(baa)
    head(baa)
    
    # Cyclones R34
    cyc34 <- read.csv(file.path(pathProEnv,"Cyclones","Env_Cyclone_R34_Station_States_2013_2025_New_Caledonia.csv"))
    cyc34 <- cyc34 %>%
      rename(any_of(crossNames_trend))
    colnames(cyc34)[grep("State|state|States|states",colnames(cyc34))] <- "Env_Cyclone_R34"
    cyc34 <- keep_evidence_cols(cyc34)
    head(cyc34)
    
    # COTS outbreaks
    cots <- read.csv(file.path(pathProEnv,"COTS","Env_COTS_Outbreaks_Station_States_New_Caledonia.csv"))
    cots <- cots %>%
      rename(any_of(crossNames_trend))
    colnames(cots)[grep("State|state|States|states",colnames(cots))] <- "Env_COTS_Outbreak"
    cots <- keep_evidence_cols(cots)
    head(cots)
    
    
  # JOIN All datasets
    
  bio <- coralG %>%
    full_join(coralT, by=c("Year","Site","Station")) %>%
    full_join(fishG,  by=c("Year","Site","Station")) %>%
    full_join(fishT,  by=c("Year","Site","Station")) %>%
    full_join(invG,   by=c("Year","Site","Station")) %>%
    full_join(invT,   by=c("Year","Site","Station"))
  
  
  
  
  
  skeleton <- bio %>%
    select(Year, Site, Station) %>%
    distinct()
  
  
  # Add environment: join by appropriate level
  master <- skeleton %>%
    left_join(geo, by=c("Site","Station")) %>%
    left_join(mpa, by=c("Site","Station")) %>%
    left_join(grav, by=c("Site","Station")) %>%
    left_join(cots, by=c("Year","Site","Station")) %>%
    left_join(baa,  by=c("Year","Site","Station")) %>%
    left_join(cyc34, by=c("Year","Site","Station")) %>%
    left_join(temp, by=c("Year","Site")) %>%
    left_join(chla, by=c("Year","Site")) %>%
    left_join(ninoa, by="Year") %>%
    left_join(cyc,   by="Year") %>%
    left_join(bio,    by=c("Year","Site","Station"))
  
  
  # Check for quality
  # Check duplicates
  stopifnot(nrow(master) == nrow(distinct(master, Year, Site, Station)))
  
  # Check missingness of key env covariates
  master %>% summarise(
    n = n(),
    miss_geomorph = mean(is.na(Env_Geomorphology)),
    miss_mpa      = mean(is.na(Env_MPA)),
    miss_gravity  = mean(is.na(Env_Gravity)),
    miss_nino     = mean(is.na(Env_Nino_Phase_General)),
    miss_cyc      = mean(is.na(Env_Cyclone_Frequency_General)),
    miss_baa      = mean(is.na(Env_Bleaching_Alert_Area)),
    miss_cyc34    = mean(is.na(Env_Cyclone_R34)),
    miss_temp     = mean(is.na(Env_Site_Temperature_Regime)),
    miss_chla     = mean(is.na(Env_Site_Chlorophyll_a_Regime)),
    miss_cots     = mean(is.na(Env_COTS_Outbreak))
    
  )
  
  lapply(master, unique)
  
  check_levels <- function(df, node) {
    lev <- sort(unique(na.omit(df[[node]])))
    net_lev <- net$getOutcomeIds(node)
    setdiff(lev, net_lev)
  }
  
  check_levels(master, "Env_Cyclone_R34")
  check_levels(master, "Env_Bleaching_Alert_Area")
  check_levels(master, "Env_Site_Temperature_Regime")
  check_levels(master, "Env_Site_Chlorophyll_a_Regime")
  
  lapply(colnames(master)[-c(1:3)], check_levels, df = master)
  
  
  
# TRAIN WITH EM ----
  
  # Keep only columns that are node IDs (and optionally keep Year/Site/Station if you want, but they won't match)
  evid <- master %>%
    select(any_of(nodes))  # strict: only network nodes
  
  # Ensure character columns (SMILE expects strings for discrete)
  evid <- evid %>%
    mutate(across(everything(), ~ if (is.factor(.x)) as.character(.x) else .x))
  
  # Optional: convert "" to NA
  evid[evid == ""] <- NA
  
  # Save so the SMILE api can read the file (NA are converted to "" for SMILE)
  write.csv(evid, file.path(pathPro,"Hdbn_Evidence_Full.csv"), na = "", row.names = FALSE)
  
  ds <- DataSet()
  ds$readFile(file.path(pathPro,"Hdbn_Evidence_Full.csv"))
  
  matching <- ds$matchNetwork(net)
  
  # Quick sanity: count matches
  length(matching)
  
  
  fixedIds <- c(
    # exogenous/environment
    "Env_Cyclone_Frequency_General",
    "Env_Nino_Phase_General",
    "Env_Site_Temperature_Regime",
    "Env_Site_Chlorophyll_a_Regime",
    "Env_Site_Conditions",
    "Env_Bleaching_Alert_Area",
    "Env_Cyclone_R34",
    "Env_COTS_Outbreak",
    "Env_Geomorphology",
    "Env_MPA",
    "Env_Gravity",
    
    # engineered / deterministic-ish
    "Cyclone_Impact",
    "Env_Environmental_Shock",
    "Acute_Local_Pressure",
    "Phase_Shift_Window",
    "Coral_Pressure_Index",
    
    # optional but recommended
    "Intrinsic_Resilience"
  )
  
  # Convert to handles (safer than passing strings in some builds)
  fixedHandles <- vapply(fixedIds, net$getNode, integer(1))
  
  
  em <- EM()
  
  # Keep original parameters (no uniformize, no randomize)
  em$setUniformizeParameters(FALSE)
  # em$setRandomizeParameters(FALSE)
  em$setRandomizeParameters(TRUE)
  
  # Stabilize with ESS (start value; adjust if too stiff/too wiggly)
  em$setEqSampleSize(0) #100
  
  # em$set_seed(111)
  # Run learning
  # em$learn(ds, net, matching, fixedNodes = fixedHandles)
  em$learn(ds, net, matching)
  
  
  getCPT(net, "Env_Site_Conditions")
  test <- getCPT(net, "Fish_Biomass")
  
  test <- net$getNodeDefinition("Live_Coral_Cover")
  
  