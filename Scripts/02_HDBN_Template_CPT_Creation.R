#############################################################################################
#                                                                                           #
#                                                                                           #
#############           Construct dynamic bayesian network                 #########
#                                                                                           #
#                                                                                           #
#############################################################################################

# What does this script do:
# Construct core model
# function to create hierarchical nodes given standard structure and adaptation with loaded data


# INIT -----

# Clean slate
rm(list=ls(all=TRUE))

# Libraries initialization
init1 <- tail(unlist(strsplit(rstudioapi::getActiveDocumentContext()$path, "/")), n = 1)
source("Scripts/00_Initialisation.R")


  # Install the package:
  # you'll need to download the package from BayesFusion:
  # https://download.bayesfusion.com/files.html?category=Academia
  # Copy/paste the .tar.gz file at the root of the Rproj
  # Then install it manually
  # DONE ONCE, THIS PART SHOULD STAY COMMENTED IF THE PACKAGE IS ALREADY INSTALLED
  # install.packages("rSMILE_2.4.0_R_x86_64-pc-linux-gnu.tar.gz", repos = NULL, type = "source")
  
library(rSMILE)
source("License.R")
# ?addTemporalArc


# DATA TESTING ----

  # net <- Network()
  # net$readFile(file.path("Data","Genie","Network1.xdsl"))
  # class(net)
  # net$getAllNodes()
  # net$setEvidence(nodeHandle = "Expert_Forecast", outcomeId = "Good")
  # net$updateBeliefs()
  # beliefs <- net$getNodeValue("Success")  
  # for(i in 1:length(beliefs)) {
  #   cat(sprintf("%s = %f\n", net$getOutcomeId("Success", i-1L), beliefs[i]))
  # }  
  

# CORE MODEL CPTs ----
  
  # Load the model
  net <- Network()
  net$readFile(file.path("Data","Genie","Hierarchical_Dynamic_Bayesian_Network_CoralReefs_RORC_Core_Template.xdsl"))    
  
  # Access nodes
  net$getAllNodeIds()
  net$getAllNodes()
  net$getNodeId(c(18))
  unlist(lapply(c(1,2), net$getNodeId))
  net$getTemporalParents("Coral_Diversity",c(1))
    
  # Inspect CPTs
  # node = "Coral_Diversity"
  node = "Coral_Reef_Ecosystem_Health"
  # node = "Branching_Coral_Cover"
  net$getOutcomeIds(node)
  net$getNodeDefinition(node)
  
    # Function to get CPT in dataframe
    getCPT <- function(net, node, temporal_order = c(1,2)) {
      states <- net$getOutcomeIds(node)
      n_states <- length(states)
      
      # parents in the same time slice
      parents_now <- net$getParents(node)
      names(parents_now) <- unlist(lapply(parents_now, net$getNodeId))
      
      # parents from previous slices (temporal)
      parents_temp <- unlist(lapply(temporal_order, function(x) {
        # x = 2
        res <- net$getTemporalParents(node, x)
        if(length(res) == 0) return(NULL)
        resNode <- res[[1]]$handle
        resId <- res[[1]]$id
        resOrder <- res[[1]]$order
        names(resNode) <- paste0(resId, "_t-",resOrder)
        return(resNode)
        
      }))
      
      # merge all parents
      parents <- c(parents_now, parents_temp)
      
      # no parents -> simple prior
      if (length(parents) == 0) {
        probs <- net$getNodeDefinition(node)
        return(data.frame(State = states, Probability = probs))
      }
      
      # list of outcomes for all parents
      parent_states <- lapply(parents, function(p) net$getOutcomeIds(p))
      
      # build all combinations of parent states
      parent_grid <- expand.grid(parent_states, KEEP.OUT.ATTRS = FALSE)
      # names(parent_grid) <- names(parents)
      
      # raw CPT definition from SMILE
      raw <- net$getNodeDefinition(node)
      
      # reshape: each parent combo gets n_states probabilities
      cpt <- cbind(parent_grid, matrix(raw, ncol = n_states, byrow = TRUE))
      names(cpt)[(ncol(cpt)-n_states+1):ncol(cpt)] <- states
      
      return(cpt)
    }
  
  
  
   

    # testing function
    getCPT(net, node)
    
    
  # Process all Core CPTs and export them to csv
  # This is for later CPT updating using a shiny app
    
    

# 
    
          
# Create "automatic" nodes from transects to related master nodes 
  # Read data
  
  # CSMP data
  # Read fish data
  # fishCsmp <- read.csv(file.path("Data","Processed","CSMP_GBR_Fish_Surveys_20182024.csv")) # , col_types = "text"
  # Read coral data
  # coralCsmp <- read.csv(file.path("Data","Processed","CSMP_GBR_Coral_PIT_20182024.csv")) # , col_types = "text"
  
  # RORC
  # Coral
  coralRorc <- read.csv(file.path(pathPro,"RORC_Coral_Station_General_States_hdbn.csv"))
  coralRorc$Station <- gsub(" ","_", coralRorc$Station)
  coralRorc$Site <- gsub(" ","_", coralRorc$Site)
  # Fish
  fishRorc <- read.csv(file.path(pathPro,"RORC_Fish_Station_General_States_hdbn.csv"))
  # Inv
  invRorc <- read.csv(file.path(pathPro,"RORC_Inv_Station_General_States_hdbn.csv"))
  
  
  # Create the nodes
  
  # Function to create CPT node :
  createTemplateCptNode <- function(net, id, name, outcomes, xPos = NULL, yPos = NULL) {
    # Create blank node
    handle <- net$addNode(net$NodeType$CPT, id)
    # Set node name (here same as id but can be more "pretty" with spaces and parenthesis for ex)
    net$setNodeName(handle, name)
    
    # Add position (optional)
    if(!is.null(xPos) &&  !is.null(yPos)) {
      net$setNodePosition(handle, xPos, yPos, 85L, 55L)
    }
    
    # Rename the original two outcomes if standard outcome (which it should be)
    net$setOutcomeId(handle, 0, outcomes[1])
    net$setOutcomeId(handle, 1, outcomes[2])
    
    # add and names outcomes > 2
    if(length(outcomes) > 2) {
      for(outcome in outcomes[3:length(outcomes)]) {
        net$addOutcome(handle, outcome)
      }
      
    }
    
    return(handle)
    
  }
  
  
  # Function that works line by line.
  # Creates the transect node and then if it doesnt exist, create the associated station, site, sector, region node.
  # Then links the transect node hierarchically to the higher level nodes
  # A function to apply to each dataset:
  createTemplateHierarchicalNodes <- function(dataset, 
                                              net, 
                                              structure = c("Region","Station"),
                                              vars = colnames(dataset)[!colnames(dataset) %in% structure]) {
    
    # DEV PART
    # dataset <- coralRorc
    # net <- net
    # structure = c("Region","Sector","Site","Station","Sample")
    # vars = colnames(dataset)[!colnames(dataset) %in% c(structure,"Year","Country")]
    # structure = c("Region","Sector","Site","Station","Sample")
    # vars = c("RorcHCB", "RorcHCM", "RorcHCO", "RorcHCT", "RorcSC", "RorcDC", "RorcSD_SI", "RorcRB", "RorcCoralRichness",
             # "RorcCoralCover")
    # eo DEV PART
    
    # Dataset should be a dataframe in two "parts" : structure columns that will be used to hierarchise nodes based on their order in structure
    # And variable columns that will demultiply the number of nodes in the structure
    # data <- cbind(dataset[,structure], dataset[,vars])
    # data <- data[!duplicated(data),]
    
    # Because we're building a template mode, we don't need all year all samples combination
    # Dataset should only be transformed into unique combinations of the structure, which will be repeated by vars

    data <- unique(dataset[, structure])
    
    # Cache Node ids and append new nodes on each loop instead of reusing net$getAllNodeIds() because its too memory intensive
    nodeList <- net$getAllNodeIds()
    
    # For each variable:
    for (var in vars) {
      # var = vars[1]
      
      cat(var,"\n")
      
      # df <- data[!is.na(data[[var]]),]
      df <- data
      
      # Loop for each dataset line
      for (i in 1:nrow(df)) {
        # i = 1
        cat("row",i,"\n")
        
        # Remember : ids (and not names) need unique identifiers. Samples/transects are not UNIQUE ! found it the hard R-crashing-style way
        region <- paste0(var, "_", "Region_", df$Region[i])
        sector <- paste0(var, "_", "Sector_", df$Sector[i])
        site <- paste0(var, "_", "Site_", df$Site[i])
        station <- paste0(var, "_", "Station_", df$Station[i])
        sample <- paste0(var, "_","Station_", df$Station[i], "_Sample_", df$Sample[i])
        
        # Create nodes if they don't exist
        for (node in c(region, sector, site, station, sample)) {
          # node = region
          if (!(node %in% nodeList)) {
            cat("node : ", node,"\n")
            # Create the node
            createTemplateCptNode(net, gsub("-","_", node), node, c("High","Medium","Low"))
            # Add to node list
            nodeList <- c(nodeList, node)
          }
        
        }
      
      # Add hierarchy arcs (variable-specific) 
        # Check for preexisting links because repeated sampling of transects. 
        # We only have to (actually must to not raise an error) define arcs once for specified nodes
        # However in this loop, region-sector arcs up to stations are asked to be redefined each sample
        # So we check if it already exist we don't need to re-arc it
        if (!(region %in% net$getParentIds(sector))) {
          net$addArc(region, sector)
        }
        
        if (!(sector %in% net$getParentIds(site))) {
          net$addArc(sector, site)
        }
        
        if (!(site %in% net$getParentIds(station))) {
          net$addArc(site, station)
        }
        
        if (!(station %in% net$getParentIds(sample))) {
          net$addArc(station, sample)
        }
        
      # This one normally doesnt exist:
      # net$addArc(station, sample)
        
      }
    
    }
    
  }
  
  
  # Creat hierarchical nodes per variable per dataset
  createTemplateHierarchicalNodes(dataset = coralRorc, 
                          net = net, 
                          structure = c("Region","Sector","Site","Station","Sample"),
                          vars = c("RorcHCB", "RorcHCM", "RorcHCO", "RorcHCT", "RorcSC", "RorcDC", "RorcSD_SI", "RorcRB", "RorcRC", "RorcCoralRichness",
                          "RorcCoralCover"))
  
  length(net$getAllNodeIds())
  
  createTemplateHierarchicalNodes(dataset = fishRorc, 
                                  net = net, 
                                  structure = c("Region","Sector","Site","Station","Sample"),
                                  vars = c("RorcFishRichness","RorcCarnivoreAbund","RorcCorallivoreAbund","RorcHerbivoreAbund"))
  
  length(net$getAllNodeIds())
  
  createTemplateHierarchicalNodes(dataset = invRorc, 
                                  net = net, 
                                  structure = c("Region","Sector","Site","Station","Sample"),
                                  vars = c("RorcInvRichness","RorcInvAbund"))
  
  length(net$getAllNodeIds())
  
  
  
  # Connect the appropriate data type region node to core model related node
  nodes <- net$getAllNodeIds()
  
  # Core nodes
  # Coral:
  # Coral_Diversity <- c("RorcCoralRichness")
  nodes[1:25]
  RegionNodes <- nodes[grep("Region", nodes)] 
  RegionNodes
  
  net$getParentIds("RorcCoralRichness_Region_PS")
  net$getParentIds("Coral_Diversity")
  net$getParentIds("Coral_Reef_Ecosystem_Productivity")
  
  net$addArc
  net$getMaxNodeTemporalOrder("Coral_Diversity")
  net$getMaxNodeTemporalOrder("RorcCoralRichness_Region_PS")
  
  # Coral_Diversity
  for(node in RegionNodes[grep("RorcCoralRichness", RegionNodes)]) {
    net$addArc("Coral_Diversity", node)
  }
  
  # Branching_Coral_Cover
  for(node in RegionNodes[grep("RorcHCB", RegionNodes)]) {
    net$addArc("Branching_Coral_Cover", node)
  }
  
  # Massive_Coral_Cover
  for(node in RegionNodes[grep("RorcHCM", RegionNodes)]) {
    net$addArc("Massive_Coral_Cover", node)
  }
  
  # Tabular_Coral_Cover
  for(node in RegionNodes[grep("RorcHCT", RegionNodes)]) {
    net$addArc("Tabular_Coral_Cover", node)
  }
  
  # Other_Forms_Coral_Cover
  for(node in RegionNodes[grep("RorcHCO", RegionNodes)]) {
    net$addArc("Other_Forms_Coral_Cover", node)
  }
  
  # Soft_Coral_Cover
  for(node in RegionNodes[grep("RorcSC", RegionNodes)]) {
    net$addArc("Soft_Coral_Cover", node)
  }
  
  # Dead_Coral_Cover
  for(node in RegionNodes[grep("RorcDC", RegionNodes)]) {
    net$addArc("Dead_Coral_Cover", node)
  }
  
  # Rock_and_Slabs_Cover
  for(node in RegionNodes[grep("RorcRC", RegionNodes)]) {
    net$addArc("Rock_and_Slabs_Cover", node)
  }
  
  # Debris_Cover
  for(node in RegionNodes[grep("RorcRB", RegionNodes)]) {
    net$addArc("Debris_Cover", node)
  }
  
  # Sand_and_Silt_Cover
  for(node in RegionNodes[grep("RorcSD_SI", RegionNodes)]) {
    net$addArc("Sand_and_Silt_Cover", node)
  }
  
  # Fish_Diversity
  for(node in RegionNodes[grep("RorcFishRichness", RegionNodes)]) {
    net$addArc("Fish_Diversity", node)
  }
  
  # Fish_Herbivores_Biomass
  for(node in RegionNodes[grep("RorcHerbivoreAbund", RegionNodes)]) {
    net$addArc("Fish_Herbivores_Biomass", node)
  }
  
  # Fish_Corallivores_Biomass
  for(node in RegionNodes[grep("RorcCorallivoreAbund", RegionNodes)]) {
    net$addArc("Fish_Corallivores_Biomass", node)
  }
  
  # Fish_Carnivores_Biomass
  for(node in RegionNodes[grep("RorcCarnivoreAbund", RegionNodes)]) {
    net$addArc("Fish_Carnivores_Biomass", node)
  }
  
  # Fish_Scraper_Biomass
  for(node in RegionNodes[grep("RorcCorallivoreAbund", RegionNodes)]) {
    net$addArc("Fish_Scraper_Biomass", node)
  }
  
  # Invertebrate_Diversity
  for(node in RegionNodes[grep("RorcInvRichness", RegionNodes)]) {
    net$addArc("Invertebrate_Diversity", node)
  }
  
  # Invertebrate_Abundance
  for(node in RegionNodes[grep("RorcInvAbund", RegionNodes)]) {
    net$addArc("Invertebrate_Abundance", node)
  }
  
  
  
  
  
  
  
  # EXPORT THE MODEL
  net$writeFile(file.path("Results", "Rorc_HDBN_Template.xdsl"))
  
  # Export the CPTs
  net$getNodeDefinition("RorcInvAbund_Station_menaku_Sample_2")
  net$getParentIds("RorcInvAbund_Station_menaku_Sample_2")
  net$getNodeDefinition("RorcInvAbund_Station_menaku")
  net$getOutcomeIds("RorcInvAbund_Station_menaku")
  
  
  # Extract structure
  # Structure, variables
  structure <- c("")
  
  
  
  
  
  
  
    
# Load CPTs from excel file
  
  
  
  
  
# Trash
  
  # Create "automatic" nodes from transects to related master nodes 
  # Read data
  
  # CSMP data
  # Read fish data
  # fishCsmp <- read.csv(file.path("Data","Processed","CSMP_GBR_Fish_Surveys_20182024.csv")) # , col_types = "text"
  # Read coral data
  # coralCsmp <- read.csv(file.path("Data","Processed","CSMP_GBR_Coral_PIT_20182024.csv")) # , col_types = "text"
  
  # RORC
  # Coral
  coralRorc <- read.csv(file.path("Data", "Species","RORC_Coral_Data_hdbn.csv"))
  coralRorc$Station <- gsub(" ","_", coralRorc$Station)
  coralRorc$Site <- gsub(" ","_", coralRorc$Site)
  # Fish
  fishRorc <- read.csv(file.path("Data", "Species","RORC_Fish_Data_hdbn.csv"))
  # Inv
  invRorc <- read.csv(file.path("Data", "Species","RORC_Inv_Data_hdbn.csv"))
  
  
  # Create the nodes
  
  # Function to create CPT node :
  createTemplateCptNode <- function(net, id, name, outcomes, xPos = NULL, yPos = NULL) {
    # Create blank node
    handle <- net$addNode(net$NodeType$CPT, id)
    # Set node name (here same as id but can be more "pretty" with spaces and parenthesis for ex)
    net$setNodeName(handle, name)
    
    # Add position (optional)
    if(!is.null(xPos) &&  !is.null(yPos)) {
      net$setNodePosition(handle, xPos, yPos, 85L, 55L)
    }
    
    # Rename the original two outcomes if standard outcome (which it should be)
    net$setOutcomeId(handle, 0, outcomes[1])
    net$setOutcomeId(handle, 1, outcomes[2])
    
    # add and names outcomes > 2
    if(length(outcomes) > 2) {
      for(outcome in outcomes[3:length(outcomes)]) {
        net$addOutcome(handle, outcome)
      }
      
    }
    
    return(handle)
    
  }
  
  
  # Function that works line by line.
  # Creates the transect node and then if it doesnt exist, create the associated station, site, sector, region node.
  # Then links the transect node hierarchically to the higher level nodes
  # A function to apply to each dataset:
  createTemplateHierarchicalNodes <- function(dataset, 
                                              net, 
                                              structure = c("Region","Sector","Site","Station","Sample"),
                                              vars = colnames(dataset)[!colnames(dataset) %in% structure]) {
    
    # DEV PART
    # dataset <- coralRorc
    # net <- net
    # structure = c("Region","Sector","Site","Station","Sample")
    # vars = colnames(dataset)[!colnames(dataset) %in% c(structure,"Year","Country")]
    # structure = c("Region","Sector","Site","Station","Sample")
    # vars = c("RorcHCB", "RorcHCM", "RorcHCO", "RorcHCT", "RorcSC", "RorcDC", "RorcSD_SI", "RorcRB", "RorcCoralRichness",
    # "RorcCoralCover")
    # eo DEV PART
    
    # Dataset should be a dataframe in two "parts" : structure columns that will be used to hierarchise nodes based on their order in structure
    # And variable columns that will demultiply the number of nodes in the structure
    # data <- cbind(dataset[,structure], dataset[,vars])
    # data <- data[!duplicated(data),]
    
    # Because we're building a template mode, we don't need all year all samples combination
    # Dataset should only be transformed into unique combinations of the structure, which will be repeated by vars
    
    data <- unique(dataset[, structure])
    
    # Cache Node ids and append new nodes on each loop instead of reusing net$getAllNodeIds() because its too memory intensive
    nodeList <- net$getAllNodeIds()
    
    # For each variable:
    for (var in vars) {
      # var = vars[1]
      
      cat(var,"\n")
      
      # df <- data[!is.na(data[[var]]),]
      df <- data
      
      # Loop for each dataset line
      for (i in 1:nrow(df)) {
        # i = 1
        cat("row",i,"\n")
        
        # Remember : ids (and not names) need unique identifiers. Samples/transects are not UNIQUE ! found it the hard R-crashing-style way
        region <- paste0(var, "_", "Region_", df$Region[i])
        sector <- paste0(var, "_", "Sector_", df$Sector[i])
        site <- paste0(var, "_", "Site_", df$Site[i])
        station <- paste0(var, "_", "Station_", df$Station[i])
        sample <- paste0(var, "_","Station_", df$Station[i], "_Sample_", df$Sample[i])
        
        # Create nodes if they don't exist
        for (node in c(region, sector, site, station, sample)) {
          # node = region
          if (!(node %in% nodeList)) {
            cat("node : ", node,"\n")
            # Create the node
            createTemplateCptNode(net, gsub("-","_", node), node, c("High","Medium","Low"))
            # Add to node list
            nodeList <- c(nodeList, node)
          }
          
        }
        
        # Add hierarchy arcs (variable-specific) 
        # Check for preexisting links because repeated sampling of transects. 
        # We only have to (actually must to not raise an error) define arcs once for specified nodes
        # However in this loop, region-sector arcs up to stations are asked to be redefined each sample
        # So we check if it already exist we don't need to re-arc it
        if (!(region %in% net$getParentIds(sector))) {
          net$addArc(region, sector)
        }
        
        if (!(sector %in% net$getParentIds(site))) {
          net$addArc(sector, site)
        }
        
        if (!(site %in% net$getParentIds(station))) {
          net$addArc(site, station)
        }
        
        if (!(station %in% net$getParentIds(sample))) {
          net$addArc(station, sample)
        }
        
        # This one normally doesnt exist:
        # net$addArc(station, sample)
        
      }
      
    }
    
  }
  
  
  # Creat hierarchical nodes per variable per dataset
  createTemplateHierarchicalNodes(dataset = coralRorc, 
                                  net = net, 
                                  structure = c("Region","Sector","Site","Station","Sample"),
                                  vars = c("RorcHCB", "RorcHCM", "RorcHCO", "RorcHCT", "RorcSC", "RorcDC", "RorcSD_SI", "RorcRB", "RorcRC", "RorcCoralRichness",
                                           "RorcCoralCover"))
  
  length(net$getAllNodeIds())
  
  createTemplateHierarchicalNodes(dataset = fishRorc, 
                                  net = net, 
                                  structure = c("Region","Sector","Site","Station","Sample"),
                                  vars = c("RorcFishRichness","RorcCarnivoreAbund","RorcCorallivoreAbund","RorcHerbivoreAbund"))
  
  length(net$getAllNodeIds())
  
  createTemplateHierarchicalNodes(dataset = invRorc, 
                                  net = net, 
                                  structure = c("Region","Sector","Site","Station","Sample"),
                                  vars = c("RorcInvRichness","RorcInvAbund"))
  
  length(net$getAllNodeIds())
  
  
  
  # Connect the appropriate data type region node to core model related node
  nodes <- net$getAllNodeIds()
  
  # Core nodes
  # Coral:
  # Coral_Diversity <- c("RorcCoralRichness")
  nodes[1:25]
  RegionNodes <- nodes[grep("Region", nodes)] 
  RegionNodes
  
  net$getParentIds("RorcCoralRichness_Region_PS")
  net$getParentIds("Coral_Diversity")
  net$getParentIds("Coral_Reef_Ecosystem_Productivity")
  
  net$addArc
  net$getMaxNodeTemporalOrder("Coral_Diversity")
  net$getMaxNodeTemporalOrder("RorcCoralRichness_Region_PS")
  
  # Coral_Diversity
  for(node in RegionNodes[grep("RorcCoralRichness", RegionNodes)]) {
    net$addArc("Coral_Diversity", node)
  }
  
  # Branching_Coral_Cover
  for(node in RegionNodes[grep("RorcHCB", RegionNodes)]) {
    net$addArc("Branching_Coral_Cover", node)
  }
  
  # Massive_Coral_Cover
  for(node in RegionNodes[grep("RorcHCM", RegionNodes)]) {
    net$addArc("Massive_Coral_Cover", node)
  }
  
  # Tabular_Coral_Cover
  for(node in RegionNodes[grep("RorcHCT", RegionNodes)]) {
    net$addArc("Tabular_Coral_Cover", node)
  }
  
  # Other_Forms_Coral_Cover
  for(node in RegionNodes[grep("RorcHCO", RegionNodes)]) {
    net$addArc("Other_Forms_Coral_Cover", node)
  }
  
  # Soft_Coral_Cover
  for(node in RegionNodes[grep("RorcSC", RegionNodes)]) {
    net$addArc("Soft_Coral_Cover", node)
  }
  
  # Dead_Coral_Cover
  for(node in RegionNodes[grep("RorcDC", RegionNodes)]) {
    net$addArc("Dead_Coral_Cover", node)
  }
  
  # Rock_and_Slabs_Cover
  for(node in RegionNodes[grep("RorcRC", RegionNodes)]) {
    net$addArc("Rock_and_Slabs_Cover", node)
  }
  
  # Debris_Cover
  for(node in RegionNodes[grep("RorcRB", RegionNodes)]) {
    net$addArc("Debris_Cover", node)
  }
  
  # Sand_and_Silt_Cover
  for(node in RegionNodes[grep("RorcSD_SI", RegionNodes)]) {
    net$addArc("Sand_and_Silt_Cover", node)
  }
  
  # Fish_Diversity
  for(node in RegionNodes[grep("RorcFishRichness", RegionNodes)]) {
    net$addArc("Fish_Diversity", node)
  }
  
  # Fish_Herbivores_Biomass
  for(node in RegionNodes[grep("RorcHerbivoreAbund", RegionNodes)]) {
    net$addArc("Fish_Herbivores_Biomass", node)
  }
  
  # Fish_Corallivores_Biomass
  for(node in RegionNodes[grep("RorcCorallivoreAbund", RegionNodes)]) {
    net$addArc("Fish_Corallivores_Biomass", node)
  }
  
  # Fish_Carnivores_Biomass
  for(node in RegionNodes[grep("RorcCarnivoreAbund", RegionNodes)]) {
    net$addArc("Fish_Carnivores_Biomass", node)
  }
  
  # Fish_Scraper_Biomass
  for(node in RegionNodes[grep("RorcCorallivoreAbund", RegionNodes)]) {
    net$addArc("Fish_Scraper_Biomass", node)
  }
  
  # Invertebrate_Diversity
  for(node in RegionNodes[grep("RorcInvRichness", RegionNodes)]) {
    net$addArc("Invertebrate_Diversity", node)
  }
  
  # Invertebrate_Abundance
  for(node in RegionNodes[grep("RorcInvAbund", RegionNodes)]) {
    net$addArc("Invertebrate_Abundance", node)
  }
  
  
  
  
  
  
  
  # EXPORT THE MODEL
  net$writeFile(file.path("Results", "Rorc_HDBN_Template.xdsl"))
  
  
  
  
  # getCPT <- function(net, node) {
  #   # Obtain number of states
  #   states <- net$getOutcomeIds(node)
  #   n_states <- length(states)
  #   
  #   # Obtain parent info
  #   parents <- net$getParents(node)
  #   
  #   if (length(parents) == 0) {
  #     # simple prior distribution
  #     probs <- net$getNodeDefinition(node)
  #     df <- data.frame(State = states, Probability = probs)
  #     return(df)
  #   }
  #   # full CPT for nodes with parents
  #   raw <- net$getNodeDefinition(node)
  #   parent_states <- lapply(parents, function(p) net$getOutcomeIds(p))
  #   
  #   # Cartesian product of parent states → columns
  #   parent_grid <- expand.grid(parent_states, KEEP.OUT.ATTRS = FALSE)
  #   names(parent_grid) <- unlist(lapply(parents, net$getNodeId))
  #   
  #   # reshape: every parent configuration has n_states probabilities
  #   cpt <- cbind(parent_grid, matrix(raw, ncol = n_states, byrow = TRUE))
  #   names(cpt)[(ncol(cpt)-n_states+1):ncol(cpt)] <- states
  #   
  #   return(cpt)
  # }
  
  
  # network_methods <- getRefClass("Network")$methods()
  # print(network_methods)
  # createCptNode = function(net, id, name, outcomes) { #, xPos, yPos
  #   # Dev
  #   # net <- net
  #   # id <-node
  #   # name = node
  #   # outcomes = c("High","Medium","Low")
  #   # eo Dev
  #   
  #   handle <- net$addNode(net$NodeType$CPT, id)
  #   net$setNodeName(handle, name)
  #   # net$setNodePosition(handle, xPos, yPos, 85L, 55L)
  #   outcomesCount <- length(outcomes)
  #   initialOutcomeCount <- net$getOutcomeCount(handle)
  #   i <- 0
  #   for (outcome in outcomes[1:initialOutcomeCount]) {
  #     net$setOutcomeId(handle, i, outcome)
  #     i <- i + 1
  #   }
  #   if ((initialOutcomeCount + 1) <= outcomesCount) {
  #     for (outcome in outcomes[(initialOutcomeCount + 1):(outcomesCount)]) {
  #       net$addOutcome(handle, outcome)
  #     }
  #   }
  #   return(handle)
  # }
  
  
  # createCptNode = function(net, id, name, outcomes, xPos, yPos) {
  #   handle <- net$addNode(net$NodeType$CPT, id)
  #   initialOutcomeCount <- net$getOutcomeCount(handle)
  #   count <- length(outcomes)
  #   net$setNodeName(handle, name)
  #   net$setNodePosition(handle, xPos, yPos, 85L, 55L)
  #   if (!is.null(outcomes)) {
  #     sapply(0:(initialOutcomeCount-1),
  #            function(x) net$setOutcomeId(handle, x, outcomes[x+1]))
  #     if (initialOutcomeCount < count) {
  #       sapply(initialOutcomeCount:(count-1),
  #              function(x) net$addOutcome(handle, outcomes[x+1]))
  #     }
  #   }
  #   return(handle)
  # }



