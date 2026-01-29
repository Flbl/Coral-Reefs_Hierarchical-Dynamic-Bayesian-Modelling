#############################################################################################
#                                                                                           #
#                                                                                           #
#############                      Extract CPT templates                         ############
#                                                                                           #
#                                                                                           #
#############################################################################################


# We load the model constructed in Genie, that has been proofed in 03_HDBN_Template_Evidence
# We extract CPTs in a versioning manner (creating a version/number folder each time this script is launched to save previous attempts)
# CPTs are then to be accessed through their csv or shiny app for manual modification/elicitation (scripts 05,06, 07)
 

# INIT -----

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
  
  
  # Function to extract list wise
  extractCPTs <- function(net, nodes, temporal_order = c(1,2)) {
    lapply(nodes, function(n) {
      message(n)
      list(
        node = n,
        cpt  = getCPT(net, n, temporal_order) #getCPT is a function in 00_Initialisation.R
      )
    })
  }
  
  
  # Create subdir for CPT version folder.
  # That means each time this code is run, if the previous number exists of the folder version, then create a new folder
  # CPT_01
  # CPT_02
  # ...
  
  # Listing existing dirs
  existing <- list.dirs(pathProCpt, full.names = TRUE, recursive = FALSE)
  
  # Transforming dirs into integers after their sequence number (if present, else NA)
  nums <- as.integer(sub("^([0-9]+).*", "\\1", basename(existing)))
  nums <- nums[!is.na(nums)]
  
  # Creating the next integer in line from existing ones
  next_id <- sprintf("%02d", ifelse(length(nums) == 0, 1, max(nums) + 1))
  
  # Creating folder from the next id
  saveFolder <- file.path(pathProCpt, paste0(next_id, "_CPT"))
  dir.create(saveFolder, showWarnings = FALSE)
  
  
  # Getting CPTs
  modelcpts <- extractCPTs(net, nodes)
  
  # Export CPTS to excel
  lapply(modelcpts, function(x){
    # x = modelcpts[[1]]
    write.csv(x$cpt, file.path(file.path(saveFolder), paste0("CPT_",x$node,".csv")), row.names = FALSE)
    
  })
  
  
  
  
  
  
  
  
# Trash -----
  
  #############################################################################################
  #                                                                                           #
  #                                                                                           #
  #############                      Extract CPT templates                         ############
  #                                                                                           #
  #                                                                                           #
  #############################################################################################
  
  
  # We load the model constructed in Genie, that has been proofed in 03_HDBN_Template_Evidence
  # We extract CPTs 
  # We launch the SHINY APP to define fixed nodes manually
  # The typical CPTs to extract include:
  # CPTs from the CORE model
  # CPTs from the chain [variable]_Site_ -> [variable]_Station -> [variable]_Station_OBS / [variable]_Station_Trend
  # CPTs from Environmental archetypes only
  
  
  # # DATA ----
  # 
  # # Load the model
  # net <- Network()
  # net$readFile(file.path(pathProGra,"Rorc_HDBN_Template.xdsl"))    
  # 
  # netOld <- Network()
  # netOld$readFile(file.path(pathGra,"Hierarchical_Dynamic_Bayesian_Network_CoralReefs_RORC_Core_Template.xdsl"))
  # 
  # # Declare CPT list:
  # cpt_registry <- list(
  #   core = grep("Env|Site|Station|Rorc", net$getAllNodeIds(), invert = TRUE, value = TRUE),
  #   
  #   # We take one branch and remove specific names
  #   # hierarchical = c("Site", "Station", "Station_OBS", "Station_Trend"),
  #   hierarchical = grep("RorcHCB.*akai|RorcHCB.*bourail", net$getAllNodeIds(), invert = FALSE, value = TRUE),
  #   
  #   environmental = grep("Env_", netOld$getAllNodeIds(), value = TRUE)
  # )
  # 
  # 
  # # Extract the CPT list:
  # getCPT(net, "Coral_Reef_Ecosystem_Health")
  # getCPT(net, "Fish_Diversity")
  # getCPT(net, "Other_Forms_Coral_Cover")
  # getCPT(netOld, "Env_Temperature_Site_XX")
  # 
  # extractCPTs <- function(net, nodes, temporal_order = c(1,2)) {
  #   lapply(nodes, function(n) {
  #     message(n)
  #     list(
  #       node = n,
  #       cpt  = getCPT(net, n, temporal_order)
  #     )
  #   })
  # }
  # 
  # coreCpts <- extractCPTs(net, cpt_registry$core)
  # 
  # envCpts <- extractCPTs(netOld, cpt_registry$environmental)
  # 
  # hierCpts <- extractCPTs(net, cpt_registry$hierarchical)
  # # Change names to generic
  # hierCpts <- lapply(hierCpts, function(x){
  #   names(x$cpt) <- gsub("akaia|bourail","",names(x$cpt))
  #   names(x$cpt) <- gsub("RorcHCB","RorcVAR",names(x$cpt))
  #   list(node = x$node <- gsub("RorcHCB|akaia|bourail","",x$node),
  #        cpt = x$cpt)
  # })
  # 
  # 
  # # Create subdir for raw exports (later combined with manually modified exports)
  # dir.create(file.path(pathProCpt,"TemplateExports"), showWarnings = FALSE)
  # 
  # # Export CPTS to excel
  # lapply(hierCpts, function(x){
  #   # x = hierCpts[[1]]
  #   write.csv(x$cpt, file.path(file.path(pathProCpt,"TemplateExports"),paste0("CPT_",x$node,".csv")), row.names = FALSE)
  #   
  # })
  # 
  # lapply(coreCpts, function(x){
  #   # x = hierCpts[[1]]
  #   write.csv(x$cpt, file.path(file.path(pathProCpt,"TemplateExports"),paste0("CPT_",x$node,".csv")), row.names = FALSE)
  #   
  # })
  # 
  # lapply(envCpts, function(x){
  #   # x = hierCpts[[1]]
  #   write.csv(x$cpt, file.path(file.path(pathProCpt,"TemplateExports"),paste0("CPT_",x$node,".csv")), row.names = FALSE)
  #   
  # })
  