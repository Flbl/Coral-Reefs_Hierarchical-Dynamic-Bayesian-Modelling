#############################################################################################
#                                                                                           #
#                                                                                           #
#############                 Fix needed nodes prior to EM                         ##########
#                                                                                           #
#                                                                                           #
#############################################################################################

# Load the HDBN
# List the nodes that should be fixed manually even before elicitation
# == some nodes are defined to reduce CPT from too many parents.
# == "Encoder node"
# Those nodes are "deterministic-ish" and should not serve EM or the network in general as "garbage solution finder"/find weird EM shortcuts
# It means "Hard coding" assumptions but its okay because they only run on logic
# e.g. If Temp is warm-ish and Chl-a is high-ish → “Above_Usual” likely
# If Temp is cool-ish and Chl-a is low-ish → “Below_Usual” likely
# Otherwise → “Usual” likely

# Here we can list easily the nodes that we want fixed, read their CPT in the shiny app and modify them manually 


# INIT ----

  # Clean slate
  rm(list=ls(all=TRUE))
  
  # Libraries initialization
  init1 <- tail(unlist(strsplit(rstudioapi::getActiveDocumentContext()$path, "/")), n = 1)
  source("Scripts/00_Initialisation.R")
  library(rSMILE)
  source("License.R")


# DATA ----

  # Load the model directly from the genie designed and exported model 
  net <- Network()
  net$readFile(file.path(pathGra,"Hierarchical_Dynamic_Bayesian_Network_CoralReefs_RORC_Template.xdsl"))    
  
  # Get CPT list
  nodes <- net$getAllNodeIds()
  
  
  # Fixing a node means “do not let EM reinterpret this variable."
  # Observed nodes are not protected from learning. However, on directly observed environmental nodes (e.g. Cyclone Season), 
  # its meaningless to try and set CPT because we are not trying to model the stochastic processes that generates environmental forcings.
  # Environmental 'encoder' nodes like Cyclone Impact, Env Site conditions... are to be fixed
  # They are not observed environmental forcing with designed mappings
  # For those, the CPT defines the concept and they must be set intentionally and fixed.
  # Some biological nodes are also encoder nodes (Phase shift window, coral pressure index, Intrinsic Resilience)
  
  toFix <- grep("Env", nodes, value = TRUE)
  toFix <- c(toFix,
             "Phase_Shift_Window",
             "Coral_Pressure_Index",
             "Intrinsic_Resilience")
  
  # Then we have nodes for biological/latent process that could use "soft-constrain" 
  # so EM can't wander into ecologically absurd regimes when data are sparse
  # So here its better to get "ecologically informed" CPT initialization
  # /!\ in this case, not to forget that EM should be done with ESS (equivalent sample size) so the initialization acts like a prior
  # We don't want hard zeros as EM can also struggle with them. Plus there is never a "zero probability" only very low ones
  # Those nodes are 
  # Coral reef ecosystem health
  # Grazing nodes
  # Recovery_capacity
  
  # note: Potentially all latent variables will be soft constrained before EM
  # Or elicited after EM with corrections
  # The goal is to get a less elicitation and more "Expert validation" than corrections
  
  toSoft <- c("Coral_Reef_Ecosystem_Health",
              "Fish_Grazing",
              "Urchin_Grazing",
              "Grazing_pressure",
              "Recovery_Capacity")
  

  
  
  source("Scripts/99_Shiny_CPT_manual_modification_App.R")
  
  run_cpt_app(
    pathProCpt = pathProCpt,
    folder = "01_CPT",
    cpt_patterns = c(toFix, toSoft)
  )

