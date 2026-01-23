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
  net$getTemporalChildren("Branching_Coral_Cover_Station_TREND_XX")
  
  getNodeSpec(net, "Branching_Coral_Cover_Station_TREND_XX")
  getNodeSpec(net, "Structural_Complexity", temporal_order = 1:3)
  getNodeSpec(net, "Fish_Diversity", temporal_order = 1:3)
  getNodeSpec(net, "Sand_and_Silt_Cover_Station_True_XX", temporal_order = 1:3)
  
  lapply(nodes,  function(x) getNodeSpec(net, x))
  
  # Process all Core CPTs and export them to csv
  # This is for later CPT updating using a shiny app 
  
  # eo core model cpts ----
  