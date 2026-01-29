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







# TRAIN WITH EM ----
  
  
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
  
  




