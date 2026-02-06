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

# --------------------------------------------------------------------------------------------------
# This script is updated to work on a server
# To make it work from the R720:
#  Need to load libraries and paths inside
# Folders are already pre-created on the server
# Data is pre-loaded on the server
# 
# --------------------------------------------------------------------------------------------------


# INIT -----

# Clean slate
rm(list=ls(all=TRUE))

# Paths
# pathPro <- file.path("Data", "01_Processed")
# pathProGra <- file.path("Data", "01_Processed", "Graph")

# Libraries initialization
# init1 <- tail(unlist(strsplit(rstudioapi::getActiveDocumentContext()$path, "/")), n = 1)
# source("Scripts/00_Initialisation.R")

library("dplyr")
library("tibble")
library("purrr")

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
net$readFile(file.path("Rorc_HDBN_Template_Soft_Constrained.xdsl"))    

# Listing nodes
# getCPT(net, "Fish_Carnivores_Biomass_Station_OBS_XX")
# getCPT(net, "Invertebrate_Abundance_Station_True_XX")

nodes <- net$getAllNodeIds()
# length(nodes)

# names(nodes) <- net$getAllNodes()

# getNodeSpec(net, "Branching_Coral_Cover_Station_TREND")
# getNodeSpec(net, "Structural_Complexity", temporal_order = 1:3)
# getNodeSpec(net, "Fish_Diversity", temporal_order = 1:3)
# getNodeSpec(net, "Sand_and_Silt_Cover_Station_True_XX", temporal_order = 1:3)

# lapply(nodes,  function(x) getNodeSpec(net, x))

# getNodeSpec(net, "Live_Coral_Cover")

# eo core model cpts ----



# TRAIN WITH EM ----
  
  # get evidence size
  # This is a value to get to tune ESS
  # ESS = 0 → pure ML
  # ESS = 100 → prior equivalent to 100 records
  # ESS ≫ N → CPTs barely move

  # data <- read.csv(file.path(pathPro,"Hdbn_Evidence_Full.csv"))
  # data <- read.csv(file.path(pathPro,"Hdbn_Evidence_Full_StationYear.csv"), check.names = FALSE)
  data <- read.csv(file.path("Hdbn_Evidence_DBNSlices.csv"), check.names = FALSE)
  # data <- read.csv(file.path(pathPro,"Hdbn_Evidence_DBNSlices_t0t2.csv"), check.names = FALSE)
  series_keys <- c("Site","Station")
  dim(data)  

  count_non_na_by_slice <- function(prefix) {
    sapply(0:2, function(s) {
      cname <- paste0(prefix, "__t", s)
      if (!cname %in% names(data)) return(NA_integer_)
      sum(!is.na(data[[cname]]))
    })
  }
  
  # count_non_na_by_slice("Fish_Scraper_Biomass_Station_TREND")
  # count_non_na_by_slice("Fish_Scraper_Biomass_Station_OBS")
  
  
  # dfNames <- names(data)
  # names(dfNames) <- seq(1:49)
  # nodes
  
  ds <- DataSet()
  # ds$readFile(file.path(pathPro,"Hdbn_Evidence_Full.csv"))
  # ds$readFile(file.path(pathPro,"Hdbn_Evidence_Full_StationYear.csv"))
  ds$readFile(file.path("Hdbn_Evidence_DBNSlices.csv"))
  # ds$readFile(file.path(pathPro,"Hdbn_Evidence_DBNSlices_t0t2.csv"))
  
  
  # ----------------------------
  # C) Build explicit matching with correct slices
  # ----------------------------
  
  # Helper to create a DataMatch object (rSMILE uses reference classes)
  make_match <- function(column_index0, node_handle, slice_index0) {
    m <- DataMatch()
    m$column <- as.integer(column_index0)  # 0-based
    m$node   <- as.integer(node_handle)    # node handle (not node id)
    m$slice  <- as.integer(slice_index0)   # time slice
    m
  }
  
  # Identify dataset columns (exclude keys)
  col_names <- colnames(data)
  data_cols <- setdiff(col_names, series_keys)
  
  # Map column name -> 0-based column index within the DataSet
  # IMPORTANT: DataMatch$column refers to the dataset column index (0-based) in file order.
  col_index0 <- setNames(seq_along(col_names) - 1L, col_names)
  # Compute maximum number of slices across all series
  Tmax <- 12
  # Tmax <- 3
  
  matching <- list()
  
  for (nid in nodes) {
    nh <- net$getNode(nid)
    # Add matches for each slice that exists in the table
    for (s in 0:(Tmax - 1L)) {
      cname <- paste0(nid, "__t", s)
      if (cname %in% col_names) {
        matching[[length(matching) + 1L]] <- make_match(col_index0[[cname]], nh, s)
      }
    }
  }
  
  length(matching)  # should be roughly (#nodes) * (Tmax) minus missing columns
  
  net$setSliceCount(Tmax)
  # net$setNumberOfSlices(Tmax)
  
  
  fixedIds <- c(
    # Observed exogenous/environment
    "Env_Cyclone_Frequency_General",
    "Env_Nino_Phase_General",
    "Env_Site_Temperature_Regime",
    "Env_Site_Chlorophyll_a_Regime",
    "Env_Bleaching_Alert_Area",
    "Env_Cyclone_R34",
    "Env_COTS_Outbreak",
    "Env_Geomorphology",
    "Env_MPA",
    "Env_Gravity",
    
    # Encoder Nodes / engineered / deterministic-ish
    "Env_Cyclone_Impact",
    "Env_Environmental_Shock",
    "Env_Acute_Local_Pressure",
    "Env_Site_Conditions",
    "Phase_Shift_Window",
    "Coral_Pressure_Index",
    "Intrinsic_Resilience"
    
  )
  
  # Fix also OBS and TREND nodes as they should "always be observed" for EM
  # fixedIds <- c(fixedIds, grep("OBS|TREND", nodes, value = TRUE))
  
  # Convert to handles (safer than passing strings in some builds)
  fixedHandles <- vapply(fixedIds, net$getNode, integer(1))
  
  # ----------------------------
  # E) Run EM learning
  # ----------------------------
  # net$BayesianAlgorithmType
  # net$getBayesianAlgorithm()
  # EPIS (recommended stochastic algorithm per BayesFusion)
  # net$setBayesianAlgorithm(net$BayesianAlgorithmType$EPIS_SAMPLING)
  
  # Quality vs runtime: start moderate; increase if estimates are noisy
  # net$setSampleCount(2000)   # try 20k-100k depending on runtime
  
  # Reproducibility
  net$setRandSeed(12345)
  
  em <- EM()
  em$setUniformizeParameters(FALSE)
  em$setRandomizeParameters(FALSE)
  
  # Equivalent sample size (ESS) is your smoothing/regularization knob
  em$setEqSampleSize(0)
  
  # em$learn(ds, net, matching)
  message("Starting learning")
  em$learn(ds, net, matching, fixedNodes = fixedHandles)
  
  
  # Now inspect a transition CPT (t1) for a plate node that truly has slice>0 evidence
  # getCPT(net, "Fish_Scraper_Biomass_Station_TREND")
  
  
  
  
  # matching <- ds$matchNetwork(net)
  # 
  # # Quick sanity: count matches
  # length(matching)
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # em <- EM()
  # 
  # # Keep original parameters (no uniformize, no randomize)
  # em$setUniformizeParameters(FALSE)
  # em$setRandomizeParameters(FALSE)
  # # em$setRandomizeParameters(TRUE)
  # 
  # # Stabilize with ESS (start value; adjust if too stiff/too wiggly)
  # em$setEqSampleSize(0) #100
  # 
  # 
  # 
  # getCPT(net, "Fish_Scraper_Biomass")
  # getCPT(net, "Fish_Scraper_Biomass_Station_OBS")
  # getCPT(net, "Fish_Scraper_Biomass_Station_TREND")
  # 
  # 
  # # em$set_seed(111)
  # # Run learning
  # # em$learn(ds, net, matching, fixedNodes = fixedHandles)
  # em$learn(ds, net, matching)
  # 
  # getCPT(net, "Fish_Scraper_Biomass")
  # getCPT(net, "Fish_Scraper_Biomass_Station_OBS")
  # getCPT(net, "Fish_Scraper_Biomass_Station_TREND")
  # 
  # getCPT(net, "Env_Site_Conditions")
  # getCPT(net, "Live_Coral_Cover")
  # getCPT(net, "Coral_Reef_Ecosystem_Health")
  # getCPT(net, "Coral_Reef_Ecosystem_Health")
  # getCPT(net, "Coral_Reef_Ecosystem_Health")
  # 
  # 
  # test <- getCPT(net, "Coral_Diversity")
  # 
  # test <- net$getNodeDefinition("Live_Coral_Cover")
  
  
  # Save the updated network
  net$writeFile(file.path("Rorc_HDBN_Template_Soft_Constrained_EM.xdsl"))
  message("Learning, model saved.")



