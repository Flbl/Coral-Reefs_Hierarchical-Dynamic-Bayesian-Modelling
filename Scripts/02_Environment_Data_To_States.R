#############################################################################################
#                                                                                           #
#                                                                                           #
#############                Convert observed data to states                   #########
#                                                                                           #
#                                                                                           #
#############################################################################################

# This scripts reads, process and convert to states the different environmental forcings used in the model:
# Temperature
# "Heat wave"
# Cyclones
# Chla
# Geomorphologie
# Travel time/Gravity
# COTS outbreak


# INIT -----

  # Clean slate
  rm(list=ls(all=TRUE))
  
  # Libraries initialization
  init1 <- tail(unlist(strsplit(rstudioapi::getActiveDocumentContext()$path, "/")), n = 1)
  source("Scripts/00_Initialisation.R")


  # eo INIT ----


# TEMPERATURE ----
  
  # 5 states : well below, below, average, over, well over
  
  
  
  
  
  
  
# HEATWAVES ----
  
  
  
  
  



