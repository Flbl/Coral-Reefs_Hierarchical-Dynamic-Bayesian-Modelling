#############################################################################################
#                                                                                           #
#                                                                                           #
#############           Construct dynamic bayesian network                 #########
#                                                                                           #
#                                                                                           #
#############################################################################################

# What does this script do:
# Read the core model made in Genie
# Constructs and exports all CPT from the core model nodes
# Constructs and exports the template for the hierarchy of all biological variables (5 states throughout all variables)
# Constructs and exports the environmental variable CPT
# function to create hierarchical nodes given standard structure and adaptation with loaded data


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
  net$readFile(file.path("Data","Graph","Hierarchical_Dynamic_Bayesian_Network_CoralReefs_RORC_Core_Template.xdsl"))    
  
  # Access nodes
  net$getAllNodeIds()
  net$getAllNodes()
  net$getNodeId(c(18))
  unlist(lapply(c(1,2), net$getNodeId))
  net$getTemporalParents("Coral_Diversity",c(1))
    
  # Inspect CPTs
  # node = "Coral_Diversity"
  # node = "Coral_Reef_Ecosystem_Health"
  # node = "Branching_Coral_Cover"
  node = "Fish_Carnivores_Biomass"
  net$getOutcomeIds(node)
  net$getNodeDefinition(node)
  
    # testing function
    getCPT(net, node)
    
    getCPT(net, "Env_Environmental_Pressure_Station_XX")
    
  # Process all Core CPTs and export them to csv
  # This is for later CPT updating using a shiny app 
    
  # eo core model cpts ----
          

# HIERARCHICAHL NODES GENERATION ----
    
  # Create "automatic" nodes from transects to related master nodes 
  # Read data
  
  # CSMP data
  # Read fish data
  # fishCsmp <- read.csv(file.path("Data","Processed","CSMP_GBR_Fish_Surveys_20182024.csv")) # , col_types = "text"
  # Read coral data
  # coralCsmp <- read.csv(file.path("Data","Processed","CSMP_GBR_Coral_PIT_20182024.csv")) # , col_types = "text"
  
  # RORC
  # Coral
  coralRorc <- read.csv(file.path(pathProSpe,"RORC_Coral_Station_General_States_hdbn.csv"))
  # Fish
  fishRorc <- read.csv(file.path(pathProSpe,"RORC_Fish_Station_General_States_hdbn.csv"))
  # Inv
  invRorc <- read.csv(file.path(pathProSpe,"RORC_Inv_Station_General_States_hdbn.csv"))
  
  
  # Create the nodes
  # Function to create CPT node :
  createTemplateCptNode <- function(net, id, name, outcomes, xPos = NULL, yPos = NULL, temporal = 0) {
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
    
    # if(temporal == TRUE) {
      # Temporal node types are placed in a temporal plate
      # there are 4 states: 
      # 0: CONTEMPORAL => not temporal
      # 1: ANCHOR => Outside the temporal plate with children inside the first time slice of the temporal plate
      # 2: TERMINAL => nodes outside the temporal plate with parents inside the plate, run on the last time slice only
      # 3: PLATE => "The" temporal type: nodes inside the temporal plate with incoming/outgoing temporal arcs
      net$setNodeTemporalType(handle, temporal)
    # }
   
    
    return(handle)
    
  }
  
  
  # Function that works line by line.
  # Creates the transect node and then if it doesnt exist, create the associated station, site, sector, region node.
  # Then links the transect node hierarchically to the higher level nodes
  # A function to apply to each dataset:
  createTemplateHierarchicalNodes <- function(dataset, 
                                              net, 
                                              structure = c("Site","Station"), # Can only be all or part of c("Country","Region","Sector","Site","Station") #removed ,"Sample"
                                              vars = colnames(dataset)[!colnames(dataset) %in% structure],
                                              colExcept = c("Year","Country","Region","Sector","n_samples","RorcCoralCover")) {
    
    # DEV PART
    # dataset <- coralRorc
    # structure = c("Site","Station")
    # colExcept = c("Year","Country","Region","Sector","n_samples","RorcCoralCover", "RorcDC")
    # vars = colnames(dataset)[!colnames(dataset) %in% c(structure,colExcept)]
    # structure = c("Region","Sector","Site","Station","Sample")
    # vars = c("RorcHCB", "RorcHCM", "RorcHCO", "RorcHCT", "RorcSC", "RorcDC", "RorcSD_SI", "RorcRB", "RorcCoralRichness","RorcCoralCover")
    # dataset <- invRorc
    # structure = c("Site","Station")
    # vars = c("RorcInvRichness","RorcInvAbund")
    
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
      # Test Zone
      # var = vars[1]
      # eo tz
      
      cat(var,"\n")
      
      # df <- data[!is.na(data[[var]]),]
      df <- data
      
      # Loop for each dataset line
      for (i in 1:nrow(df)) {
        # i = 1
        cat("row",i,"\n")
        
        # Remember : ids (and not names) need unique identifiers. Samples/transects are not UNIQUE ! found it the hard R-crashing-style way
        if("Country" %in% structure) country <- paste0(var, "_", "Region_", df$Region[i]) else country <- NULL
        if("Region" %in% structure) region <- paste0(var, "_", "Region_", df$Region[i]) else region <- NULL
        if("Sector" %in% structure) sector <- paste0(var, "_", "Sector_", df$Sector[i]) else sector <- NULL
        if("Site" %in% structure) site <- paste0(var, "_", "Site_", df$Site[i]) else site <- NULL
        # Station is subdivided into 3 nodes
        if("Station" %in% structure) station <- c(paste0(var, "_", "Station_", df$Station[i]), paste0(var, "_", "Station_", df$Station[i],"_OBS"),paste0(var, "_", "Station_", df$Station[i],"_Trend")) else station <- NULL
        # if("Sample" %in% structure) sample <- paste0(var, "_","Station_", df$Station[i], "_Sample_", df$Sample[i]) else sample <- NULL
        
        # Create up to site nodes if they don't exist
        for (node in c(country,region, sector, site)) {
          # node = site
          if (!(node %in% nodeList)) {
            cat("node : ", node,"\n")
            # Create the node
            createTemplateCptNode(net, gsub("-","_", node), node, c("Very_Good","Good","Medium","Low","Zero"), temporal = 3)
            # Add to node list
            nodeList <- c(nodeList, node)
          }
        }
        
        # Create station / Station OBS / Station Trend (we lose sample)
        for (node in c(station)) {
          # node = station
          if (!(node %in% nodeList)) {
            cat("node : ", node,"\n")
            # Create the node
            if(length(grep("Trend",node)) == 0){
              createTemplateCptNode(net, gsub("-","_", node), node, c("Very_Good","Good","Medium","Low","Zero"), temporal = 3)
            } else {
              createTemplateCptNode(net, gsub("-","_", node), node, c("Recovering","Stable","Degrading"), temporal = 3)
            }
            # Add to node list
            nodeList <- c(nodeList, node)
          }
        }
        
        
        
      # Add hierarchy arcs (variable-specific) 
        # Check for preexisting links because repeated sampling of transects. 
        # We only have to (actually must to not raise an error) define arcs once for specified nodes
        # However in this loop, region-sector arcs up to stations are asked to be redefined each sample
        # So we check if it already exist we don't need to re-arc it
        # Plus we check if structure level exists
        if(!any(is.null(region),is.null(sector))) {
          if (!(region %in% net$getParentIds(sector))) {
            net$addArc(region, sector)
          }
        }
        
        if(!any(is.null(sector), is.null(site))) {
          if (!(sector %in% net$getParentIds(site))) {
            net$addArc(sector, site)
          }
        }
        
        trueStation <- station[-grep("OBS|Trend",station)]
        
        # Site -> true Station
        if(!any(is.null(site), is.null(station))) {
          if (!(site %in% net$getParentIds(station))) {
            net$addArc(site, trueStation)
          }
        }
        
        # True station to obs
        if(!(trueStation %in% net$getParentIds(station[grep("OBS",station)]))) {
          net$addArc(trueStation, station[grep("OBS",station)])
        }  
        
        # True station to trend
        # True station to obs
        if(!(trueStation %in% net$getParentIds(station[grep("Trend",station)]))) {
          net$addArc(trueStation, station[grep("Trend",station)])
        }  
        
        # Add temporal arc to true station
        net$addTemporalArc(trueStation, trueStation, 1)
        
        # Add temporal arc to true station from trend station
        net$addTemporalArc(station[grep("Trend",station)], trueStation, 1)
        
        
      # This one normally doesnt exist since we removed samples:
      # net$addArc(station, sample)
        # if (!(station %in% net$getParentIds(sample))) {
        #   net$addArc(station, sample)
        # }
        
      }
    
    }
    
  }
  
  
  # Create hierarchical nodes per variable per dataset
  # Coral branches
  createTemplateHierarchicalNodes(dataset = coralRorc, 
                          net = net, 
                          structure = c("Site","Station"),
                          vars = c("RorcHCB", "RorcHCM", "RorcHCO", "RorcHCT", "RorcSC", "RorcSD_SI", "RorcRB", "RorcRC", "RorcCoralRichness"))
  
  length(net$getAllNodeIds())
  
  # Fish branches
  createTemplateHierarchicalNodes(dataset = fishRorc, 
                                  net = net, 
                                  structure = c("Site","Station"),
                                  vars = c("RorcFishRichness","RorcCarnivoreAbund","RorcCorallivoreAbund","RorcHerbivoreAbund","RorcScraperAbund"))
  
  length(net$getAllNodeIds())
  
  # Invertebrate branches
  createTemplateHierarchicalNodes(dataset = invRorc, 
                                  net = net, 
                                  structure = c("Site","Station"),
                                  vars = c("RorcInvRichness","RorcInvAbund"))
  
  length(net$getAllNodeIds())
  
  
  
  # Connect the appropriate data type hierarchical nodes to core model related nodes
  nodes <- net$getAllNodeIds()
  
  
  # Checks
  # Coral:
  # Coral_Diversity <- c("RorcCoralRichness")
  nodes[1:25]
  # RegionNodes <- nodes[grep("Region", nodes)] 
  # RegionNodes
  # net$getParentIds("RorcCoralRichness_Region_PS")
  # net$getMaxNodeTemporalOrder("RorcCoralRichness_Region_PS")
  net$getParentIds("Coral_Diversity")
  net$getChildIds("Coral_Diversity")
  net$getParentIds("Coral_Reef_Ecosystem_Productivity")
  net$getMaxNodeTemporalOrder("Coral_Diversity")
  net$getMaxNodeTemporalOrder("RorcCoralRichness_Site_drueulu")
  
  
  # Get chosen hierarchical rank nodes to connect to main core nodes
  rankNodes <- nodes[grep("Site",nodes)]
  
  # CHOSEN RANK WISE NODE CONNECTION
  # Coral_Diversity
  for(node in rankNodes[grep("RorcCoralRichness", rankNodes)]) {
    net$addArc("Coral_Diversity", node)
  }
  
  # Branching_Coral_Cover
  for(node in rankNodes[grep("RorcHCB", rankNodes)]) {
    net$addArc("Branching_Coral_Cover", node)
  }
  
  # Massive_Coral_Cover
  for(node in rankNodes[grep("RorcHCM", rankNodes)]) {
    net$addArc("Massive_Coral_Cover", node)
  }
  
  # Tabular_Coral_Cover
  for(node in rankNodes[grep("RorcHCT", rankNodes)]) {
    net$addArc("Tabular_Coral_Cover", node)
  }
  
  # Other_Forms_Coral_Cover
  for(node in rankNodes[grep("RorcHCO", rankNodes)]) {
    net$addArc("Other_Forms_Coral_Cover", node)
  }
  
  # Soft_Coral_Cover
  for(node in rankNodes[grep("RorcSC", rankNodes)]) {
    net$addArc("Soft_Coral_Cover", node)
  }
  
  # Dead_Coral_Cover
  # for(node in rankNodes[grep("RorcDC", rankNodes)]) {
  #   net$addArc("Dead_Coral_Cover", node)
  # }
  
  # Rock_and_Slabs_Cover
  for(node in rankNodes[grep("RorcRC", rankNodes)]) {
    net$addArc("Rock_and_Slabs_Cover", node)
  }
  
  # Debris_Cover
  for(node in rankNodes[grep("RorcRB", rankNodes)]) {
    net$addArc("Debris_Cover", node)
  }
  
  # Sand_and_Silt_Cover
  for(node in rankNodes[grep("RorcSD_SI", rankNodes)]) {
    net$addArc("Sand_and_Silt_Cover", node)
  }
  
  # Fish_Diversity
  for(node in rankNodes[grep("RorcFishRichness", rankNodes)]) {
    net$addArc("Fish_Diversity", node)
  }
  
  # Fish_Herbivores_Biomass
  for(node in rankNodes[grep("RorcHerbivoreAbund", rankNodes)]) {
    net$addArc("Fish_Herbivores_Biomass", node)
  }
  
  # Fish_Corallivores_Biomass
  for(node in rankNodes[grep("RorcCorallivoreAbund", rankNodes)]) {
    net$addArc("Fish_Corallivores_Biomass", node)
  }
  
  # Fish_Carnivores_Biomass
  for(node in rankNodes[grep("RorcCarnivoreAbund", rankNodes)]) {
    net$addArc("Fish_Carnivores_Biomass", node)
  }
  
  # Fish_Scraper_Biomass
  for(node in rankNodes[grep("RorcScraperAbund", rankNodes)]) {
    net$addArc("Fish_Scraper_Biomass", node)
  }
  
  # Invertebrate_Diversity
  for(node in rankNodes[grep("RorcInvRichness", rankNodes)]) {
    net$addArc("Invertebrate_Diversity", node)
  }
  
  # Invertebrate_Abundance
  for(node in rankNodes[grep("RorcInvAbund", rankNodes)]) {
    net$addArc("Invertebrate_Abundance", node)
  }
  
  
  # ENVIRONMENTAL VARIABLES
  # Read all env drivers, separate hierarchy levels, then loop over General/Site/station-variable to Create and Connect environmental drivers nodes
  
  # Read all env drivers
  # General wise
  # Cyclone frequency over New Caledonia
  cyc <- read.csv(file.path(pathDat,"01_Processed","Environment","Cyclones","Env_StormCount_General_States_New_Caledonia.csv"), header = TRUE)
  
  # Nino phase:
  ninoa <- read.csv(file.path(pathDat,"01_Processed","Environment","OceanicNinoIndex","Env_General_ONI_New_Caledonia.csv"), header = TRUE)
  
  
  # Site wise
  # Temperature regime
  temp <- read.csv(file.path(pathDat,"01_Processed","Environment","Temperature","RORC_Env_TemperatureRegime_Site_States_hdbn.csv"),header = TRUE)
  
  # Chla regime
  chla <- read.csv(file.path(pathDat,"01_Processed","Environment","Chlorophyll_a","RORC_Env_ChlorophyllaRegime_Site_States_hdbn.csv"),header = TRUE)
  
  
  # Station wise
  # Cyclone R34 influence
  cycR34 <- read.csv(file.path(pathDat,"01_Processed","Environment","Cyclones","Env_R34StormCount_Station_States_New_Caledonia.csv"), header = TRUE)
  
  # Bleaching alert area
  baa <- read.csv(file.path(pathDat,"01_Processed","Environment","Heatwaves_BAA","Env_HeatwavesBAA_Station_States_New_Caledonia.csv"), header = TRUE)
  
  # Cots
  cots <- read.csv(file.path(pathDat,"01_Processed","Environment","COTS","Env_COTS_Outbreaks_Station_States_New_Caledonia.csv"), header = TRUE)
  
  # Gravity
  grav <- read.csv(file.path(pathDat,"01_Processed","Environment","Gravity","Env_Gravity_Station_States_New_Caledonia.csv"), header = TRUE)
  
  
  # Replicating env archetype nodes per spatial unit (General/Station/Site/) per variable
  nodes[grep("Env",nodes)]
  net$getParentIds("Env_Local_Pressure_Station_XX")
  
  # Function to get genie generated env archetype node specifications
  getArchetypeSpec <- function(net, archetype) {
    list(
      outcomes = net$getOutcomeIds(archetype),
      parents  = net$getParentIds(archetype),
      children = net$getChildIds(archetype),
      temporal = net$getNodeTemporalType(archetype)
    )
  }
  
  getArchetypeSpec(net, "Env_Temperature_Site_XX")
  getArchetypeSpec(net, "Env_Bleaching_Alert_Area_Station_XX")
  getArchetypeSpec(net, "Env_Gravity_Station_XX")
  getArchetypeSpec(net, "Env_Local_Pressure_Station_XX")
  getArchetypeSpec(net, "Env_Nino_Phase_General")
  
  
  # Helper functions
  makeSpatialNode <- function(template, station) {
    gsub("XX", station, template)
  }
  
  ensureNode <- function(net, node_id, archetype) {
    if (!(node_id %in% net$getAllNodeIds())) {
      info <- getArchetypeSpec(net, archetype)
      createTemplateCptNode(
        net,
        id = node_id,
        name = node_id,
        outcomes = info$outcomes,
        temporal = info$temporal
      )
    }
  }
  
  ensureArc <- function(net, parent, child) {
    if (!(child %in% net$getChildIds(parent))) {
      net$addArc(parent, child)
    }
  }
  
  # Function to create environmental nodes from model template archetypes and csv data
  # No direct observed csv connect to station/site nodes. So we don't need to wire the arcs to all stations/sites in this function
  createEnvNodesFromCSV <- function(net,
                                    csv,
                                    archetype,
                                    spatial_level = c("Station", "Site"), # , "General" # We remove General as we can already wire the nodes in Genie
                                    # state_col = "State",
                                    station_col = "Station",
                                    site_col = "Site") {
    
    # Test zone
    # net = net
    # csv = baa
    # archetype = "Env_Bleaching_Alert_Area_Station_XX"
    # spatial_level = c("Station")
    # state_col = "State"
    # station_col = "station"
    # site_col = "site"
    
    # net = net
    # csv = temp
    # archetype = "Env_Temperature_Site_XX"
    # spatial_level = c("Site")
    # state_col = "SST_regime_state"
    # station_col = ""
    # site_col = "site"
    
    # net = net
    # csv = cycR34
    # archetype = "Env_Cyclone_R34_Station_XX"
    # spatial_level = c("Station")
    # state_col = "States"
    # station_col = "Station"
    # site_col = ""
    
    # net = net
    # csv = cots
    # archetype = "Env_COTS_Outbreak_Station_XX"
    # spatial_level = c("Station")
    # state_col = "States"
    # station_col = "Station"
    # site_col = ""

    # eo tz
    
    
    # Ensure spatial_level is one of the allowed values
    # This avoids silent errors and standardises downstream logic
    # This runs only when inside the function
    spatial_level <- match.arg(spatial_level)
    
    # # Raise errors for site/station names
    # if (spatial_level == "Site" && !site_col %in% names(csv)) {
    #   stop(sprintf("site_col '%s' not found in CSV", site_col))
    # }
    # 
    # if (spatial_level == "Station" && !station_col %in% names(csv)) {
    #   stop(sprintf("station_col '%s' not found in CSV", station_col))
    # }
    
    
    # Retrieve structural metadata from the archetype environmental node:
    # - its outcomes (not used here directly)
    # - its parents and children (critical for reconnecting structure)
    # - its temporal type (PLATE / ANCHOR / etc.)
    info <- getArchetypeSpec(net, archetype)
    
    # Determine which spatial units the node must be replicated over.
    # This is inferred directly from the CSV content.
    # Station-level drivers → one node per Station
    # Site-level drivers    → one node per Site
    # General drivers       → a single node (no spatial replication)
    units <- switch(
      spatial_level,
      Station = unique(csv[[station_col]]),
      Site    = unique(csv[[site_col]])
      # General = "General"
    )
    
    # Extract the categorical states of the environmental variable.
    # These define the node outcomes (e.g. "No alert","Bleaching alert","Mortality alert" for baa).
    # IMPORTANT:
    # - We only define the outcome space here
    # - CPT values will be learned / set later
    # states <- unique(as.character(csv[[state_col]]))
    
    # Cache all existing node IDs in the network.
    # This avoids repeated expensive calls to net$getAllNodeIds()
    nodeList <- net$getAllNodeIds()
    
    # Iterate over each spatial replication unit
    for (u in units) {
      # u = "akaia"
      # u = "chateaubriand"
      
      
      # Construct the new node ID by replacing the "XX" placeholder
      # used in Genie archetypes with the actual spatial unit.
      new_id <- makeSpatialNode(archetype, u)
      
      # Create the replicated environmental node
      if (!(new_id %in% nodeList)) {
        createTemplateCptNode(
          net,
          id       = new_id,
          name     = new_id,
          outcomes = info$outcomes,
          temporal = info$temporal
        )
        nodeList <- c(nodeList, new_id)
      }
 
      # Reconnect the archetype children (if any)
      for (child in info$children) {
        # child = info$children[1]
        
        # If there is another archetype node in the children we create its name
        child_id <- if (grepl("XX", child)) {
          makeSpatialNode(child, u)
        } else {
          child
        }
        
        # Ensure child exists
        if (!(child_id %in% nodeList)) {
          child_info <- getArchetypeSpec(net, child)
          createTemplateCptNode(
            net,
            id       = child_id,
            name     = child_id,
            outcomes = child_info$outcomes,
            temporal = child_info$temporal
          )
          nodeList <- c(nodeList, child_id)
        }
        
        # Add structural arc
        ensureArc(net, new_id, child_id)
        
        # Add temporal arc for cyclone forcing
        if (grepl("Env_Cyclone_R34", new_id)) {
          net$addTemporalArc(new_id, child_id, 1)
        }
      }
      
      
      # Station level local pressure node logic (specific)
      if (spatial_level == "Station") {
        
        env_pressure_id <- makeSpatialNode(
          "Env_Environmental_Pressure_Station_XX", u
        )
        
        local_pressure_id <- makeSpatialNode(
          "Env_Local_Pressure_Station_XX", u
        )
        
        # Ensure both nodes exist
        ensureNode(net,
                   env_pressure_id,
                   "Env_Environmental_Pressure_Station_XX")
        
        ensureNode(net,
                   local_pressure_id,
                   "Env_Local_Pressure_Station_XX")
        
        # Ensure arc exists
        ensureArc(net, env_pressure_id, local_pressure_id)
      }
    }
  }
  
  
  # Create environmental nodes from csv
  # Initial node number check
  length(net$getAllNodeIds())
  
  # Temperature regime
  createEnvNodesFromCSV(
    net        = net,
    csv        = temp,
    archetype  = "Env_Temperature_Site_XX", # nodes[grep("Temperature",nodes)]
    spatial_level = "Site",
    site_col = "site"
  )
  
  # Checks
  # nodes <- net$getAllNodeIds()
  length(net$getAllNodeIds())
  # nodes[grep("Temperature",nodes)]
  # nodes[grep("Conditio",nodes)]
  
  # Chla regime
  createEnvNodesFromCSV(
    net        = net,
    csv        = chla,
    archetype  = nodes[grep("Chloro",nodes)],
    spatial_level = "Site",
    site_col = "site"
  )
  
  # Checks
  # nodes <- net$getAllNodeIds()
  length(net$getAllNodeIds())
  # nodes[grep("Chloro",nodes)]
  # nodes[grep("Conditio",nodes)]
  # getArchetypeSpec(net, "Env_Site_bourail_Conditions")
  
  # Gravity
  createEnvNodesFromCSV(
    net        = net,
    csv        = grav,
    archetype  = nodes[grep("Grav",nodes)],
    spatial_level = "Station",
    station_col = "Station"
  )
  length(net$getAllNodeIds())
  
  
  # COTS oubreaks
  createEnvNodesFromCSV(
    net        = net,
    csv        = cots,
    archetype  = nodes[grep("COTS",nodes)],
    spatial_level = "Station",
    station_col = "Station"
  )
  length(net$getAllNodeIds())
  
  # BAA
  createEnvNodesFromCSV(
    net        = net,
    csv        = baa,
    archetype  = nodes[grep("Bleaching",nodes)],
    spatial_level = "Station",
    station_col = "station"
  )
  length(net$getAllNodeIds())
  
  # Cyclone R34
  createEnvNodesFromCSV(
    net        = net,
    csv        = cycR34,
    archetype  = "Env_Cyclone_R34_Station_XX",
    spatial_level = "Station",
    station_col = "Station"
  )
  length(net$getAllNodeIds())
  
  
  # Checks
  nodes <- net$getAllNodeIds()
  nodes[grep("Temp",nodes)]
  nodes[grep("Chloro",nodes)]
  nodes[grep("COTS",nodes)]
  nodes[grep("Env",nodes)]
  nodes[grep("Conditio",nodes)]
  nodes[grep("Local",nodes)]
  nodes[grep("Grav",nodes)]
  getArchetypeSpec(net, "Env_Conditions_Site_bourail")
  # getArchetypeSpec(net, "Env_Local_Pressure_Station_port_boise")
  getArchetypeSpec(net, "Env_Local_Pressure_Station_yejele")
  getArchetypeSpec(net, "Env_Gravity_Station_yejele")
  getArchetypeSpec(net, "Env_Local_Pressure_Station_hiengabat")
  getArchetypeSpec(net, "Env_Gravity_Station_hiengabat")
  
  # Delete all archetypal nodes
  todel <- nodes[grep("XX",nodes)]
  lapply(todel, net$deleteNode)
  
  getCPT(net,"Env_Environmental_Pressure_Station_akaia")
  getCPT(net,"Env_Local_Pressure_Station_akaia")
  
  # Get list of environmental nodes connecting to the site/station-variable
  eco_vars <- c(
    "RorcCoralRichness",
    "RorcHCB",
    "RorcHCM",
    "RorcHCT",
    "RorcHCO",
    "RorcSC",
    "RorcRC",
    "RorcRB",
    "RorcSD_SI",
    "RorcFishRichness",
    "RorcHerbivoreAbund",
    "RorcCorallivoreAbund",
    "RorcCarnivoreAbund",
    "RorcScraperAbund",
    "RorcInvRichness",
    "RorcInvAbund"
  )
  
  
  # Loop over the env nodes to create an arc to each site/station wise variables
  # Create a function to link all of this:
  
  linkCompositeDrivers <- function(net, eco_vars) {
    
    nodes <- net$getAllNodeIds()
    
    ## ------------------------------------------------
    ## 1. SITE-LEVEL: Env_Site_<site>_Conditions
    ## ------------------------------------------------
    # site_cond_nodes <- grep("^Env_Site_.*_Conditions$", nodes, value = TRUE)
    site_cond_nodes <- grep("^Env_Conditions_Site_.*", nodes, value = TRUE)
    
    for (cond_node in site_cond_nodes) {
      # cond_node = "Env_Conditions_Site_bourail"
      
      message(cond_node)
      
      # Extract site name
      # site <- sub("^Env_Site_(.*)_Conditions$", "\\1", cond_node)
      site <- sub("^Env_Conditions_Site_", "", cond_node)
      
      for (v in eco_vars) {
        # v = "RorcHCO"
        
        message(v)
        
        target <- paste0(v, "_Site_", site)
        
        # getArchetypeSpec(net, cond_node)
        # ensureArc
        # target %in% net$getChildIds(cond_node)
        ensureArc(net, cond_node, target)
      }
    }
    
    ## ------------------------------------------------
    ## 2. STATION-LEVEL: Env_Local_Pressure_Station_<station>
    ## ------------------------------------------------
    station_press_nodes <- grep("^Env_Local_Pressure_Station_", nodes, value = TRUE)
    
    for (press_node in station_press_nodes) {
      # press_node = "Env_Local_Pressure_Station_akaia"
      
      message(press_node)
      
      # Extract station name
      station <- sub("^Env_Local_Pressure_Station_", "", press_node)
      
      for (v in eco_vars) {
        
        message(v)
        
        target <- paste0(v, "_Station_", station)
        
        ensureArc(net, press_node, target)
      }
    }
  }
  
  # getArchetypeSpec(net, "Env_Local_Pressure_Station_paradis")
  getArchetypeSpec(net, "RorcScraperAbund_Site_bourail")
  # nodes[grep("RorcScraperAbund",nodes)]
  # nodes[grep("Scraper",nodes)]
  
  linkCompositeDrivers(net, eco_vars)
  
  
  # CHECKING TIME
  # Helper Functions to check the model
  # No parents (only specific env nodes should have no parents)
  getNodesWithoutParents <- function(net) {
    nodes <- net$getAllNodeIds()
    nodes[sapply(nodes, function(n) length(net$getParentIds(n)) == 0)]
  }
  roots <- getNodesWithoutParents(net)
  length(roots)
  
  # No children
  getNodesWithoutChildren <- function(net) {
    nodes <- net$getAllNodeIds()
    nodes[sapply(nodes, function(n) length(net$getChildren(n)) == 0)]
  }
  test <- getNodesWithoutChildren(net)
  test[grep("OBS|Trend",test, invert = TRUE)]
  
  # Nodes without parents or children?
  getIsolatedNodes <- function(net) {
    nodes <- net$getAllNodeIds()
    nodes[
      sapply(nodes, function(n)
        length(net$getParentIds(n)) == 0 &&
          length(net$getChildren(n)) == 0
      )
    ]
  }
  getIsolatedNodes(net)
  
  
  # Nodes with only temporal parents or children
  getTemporalOnlyNodes <- function(net) {
    nodes <- net$getAllNodeIds()
    
    nodes[sapply(nodes, function(n) {
      parents <- net$getParentIds(n)
      children <- net$getChildren(n)
      
      has_structural <-
        length(parents) > 0 || length(children) > 0
      
      has_temporal <-
        net$getMaxNodeTemporalOrder(n) != 0
      
      has_temporal && !has_structural
    })]
  }
  getTemporalOnlyNodes(net)
  
  
  getNodesWithEmptyCPT <- function(net) {
    nodes <- net$getAllNodeIds()
    
    nodes[sapply(nodes, function(n) {
      cpt <- try(net$getNodeDefinition(n), silent = TRUE)
      inherits(cpt, "try-error") || length(cpt) == 0
    })]
  }
  getNodesWithEmptyCPT(net)
  
  
  # Outcome consistency
  checkOutcomeConsistency <- function(net, pattern) {
    nodes <- grep(pattern, net$getAllNodeIds(), value = TRUE)
    
    if (length(nodes) < 2) return(invisible(NULL))
    
    ref <- net$getOutcomeIds(nodes[1])
    
    bad <- nodes[sapply(nodes, function(n) {
      !identical(net$getOutcomeIds(n), ref)
    })]
    
    bad
  }
  checkOutcomeConsistency(net, "Env_Temperature_Site_")
  checkOutcomeConsistency(net, "Env_Chlorophylle_a_Site_")
  
  length(net$getAllNodeIds())
  
  # EXPORT THE MODEL
  net$writeFile(file.path(pathDat,"01_Processed","Graph", "Rorc_HDBN_Template.xdsl"))
  
  
  
  
# Trash
  
  # Old function working (remove all child2 stuff)
  # createEnvNodesFromCSV <- function(net,
  #                                   csv,
  #                                   archetype,
  #                                   spatial_level = c("Station", "Site"), # , "General" # We remove General as we can already wire the nodes in Genie
  #                                   # state_col = "State",
  #                                   station_col = "Station",
  #                                   site_col = "Site") {
  #   
  #   # Test zone
  #   # net = net
  #   # csv = baa
  #   # archetype = "Env_Bleaching_Alert_Area_Station_XX"
  #   # spatial_level = c("Station")
  #   # state_col = "State"
  #   # station_col = "station"
  #   # site_col = "site"
  #   
  #   # net = net
  #   # csv = temp
  #   # archetype = "Env_Temperature_Site_XX"
  #   # spatial_level = c("Site")
  #   # state_col = "SST_regime_state"
  #   # station_col = ""
  #   # site_col = "site"
  #   
  #   # net = net
  #   # csv = cycR34
  #   # archetype = "Env_Cyclone_R34_Station_XX"
  #   # spatial_level = c("Station")
  #   # state_col = "States"
  #   # station_col = "Station"
  #   # site_col = ""
  #   
  #   # eo tz
  #   
  #   
  #   # Ensure spatial_level is one of the allowed values
  #   # This avoids silent errors and standardises downstream logic
  #   # This runs only when inside the function
  #   spatial_level <- match.arg(spatial_level)
  #   
  #   # # Raise errors for site/station names
  #   # if (spatial_level == "Site" && !site_col %in% names(csv)) {
  #   #   stop(sprintf("site_col '%s' not found in CSV", site_col))
  #   # }
  #   # 
  #   # if (spatial_level == "Station" && !station_col %in% names(csv)) {
  #   #   stop(sprintf("station_col '%s' not found in CSV", station_col))
  #   # }
  #   
  #   
  #   # Retrieve structural metadata from the archetype environmental node:
  #   # - its outcomes (not used here directly)
  #   # - its parents and children (critical for reconnecting structure)
  #   # - its temporal type (PLATE / ANCHOR / etc.)
  #   info <- getArchetypeSpec(net, archetype)
  #   
  #   # Determine which spatial units the node must be replicated over.
  #   # This is inferred directly from the CSV content.
  #   # Station-level drivers → one node per Station
  #   # Site-level drivers    → one node per Site
  #   # General drivers       → a single node (no spatial replication)
  #   units <- switch(
  #     spatial_level,
  #     Station = unique(csv[[station_col]]),
  #     Site    = unique(csv[[site_col]])
  #     # General = "General"
  #   )
  #   
  #   # Extract the categorical states of the environmental variable.
  #   # These define the node outcomes (e.g. "No alert","Bleaching alert","Mortality alert" for baa).
  #   # IMPORTANT:
  #   # - We only define the outcome space here
  #   # - CPT values will be learned / set later
  #   # states <- unique(as.character(csv[[state_col]]))
  #   
  #   # Cache all existing node IDs in the network.
  #   # This avoids repeated expensive calls to net$getAllNodeIds()
  #   nodeList <- net$getAllNodeIds()
  #   
  #   # Iterate over each spatial replication unit
  #   for (u in units) {
  #     # u = "akaia"
  #     # u = "chateaubriand"
  #     
  #     
  #     # Construct the new node ID by replacing the "XX" placeholder
  #     # used in Genie archetypes with the actual spatial unit.
  #     #
  #     # Examples:
  #     #   Env_Cyclone_R34_Station_XX → Env_Cyclone_R34_Station_jua
  #     #   Env_Chlorophylle_a_Site_XX → Env_Chlorophylle_a_Site_bourail
  #     #   Env_General_Cyclone_Season stays the same
  #     new_id <- switch(
  #       spatial_level,
  #       Station = gsub("XX", u, archetype),
  #       Site    = gsub("XX", u, archetype)
  #       # General = gsub("_XX", "", archetype)
  #     )
  #     
  #     
  #     # If the node already exists (e.g. repeated CSV rows),
  #     # skip creation to avoid SMILE errors
  #     if (new_id %in% nodeList) next
  #     
  #     # Create a new CPT node using the generated ID
  #     # The node is structurally identical to the archetype
  #     handle <- createTemplateCptNode(net, 
  #                                     id = new_id, 
  #                                     name = new_id,
  #                                     outcomes = info$outcomes,
  #                                     temporal = info$temporal
  #     )
  #     
  #     # Update the cached node list to include the newly created node
  #     nodeList <- c(nodeList, new_id)
  #     
  #     # Reconnect the new node to the children of the archetype.
  #     #
  #     # This ensures that:
  #     # - Station/Site-specific environmental drivers influence
  #     #   the correct ecological or composite nodes
  #     # - Existing causal logic in the core model is preserved
  #     for (child in info$children) {
  #       # child = info$children[1]
  #       
  #       # If the child node is also spatially templated (contains "XX"),
  #       # Replicate it to create a new node (check if it exists as well)
  #       # Change its name with the same spatial unit
  #       child_id <- if (grepl("XX", child)) gsub("XX", u, child) else child
  #       # Check if it exists already
  #       # Add the arc but create the child it doesn't exist in the network.
  #       if (child_id %in% nodeList) {
  #         net$addArc(new_id, child_id)
  #         
  #         # Add temporal links
  #         if(grepl("Env_Cyclone_R34", new_id)) net$addTemporalArc(new_id, child_id, 1)
  #         
  #       } else {
  #         info_Child <- getArchetypeSpec(net, child)
  #         
  #         child_handle <- createTemplateCptNode(net, 
  #                                               id = child_id, 
  #                                               name = child_id,
  #                                               outcomes = info_Child$outcomes,
  #                                               temporal = info_Child$temporal
  #         )
  #         # Add the arc
  #         net$addArc(new_id, child_id)
  #         
  #         # Add temporal links
  #         if(grepl("Env_Cyclone_R34", new_id)) net$addTemporalArc(new_id, child_id, 1)
  #         
  #         # Update the cached node list to include the newly created node
  #         nodeList <- c(nodeList, child_id)
  #         
  #         # IF the child is "Env_Local_Pressure_Station_XX", we need to add an arc to it from the children 
  #         if(length(info_Child$children) > 0 && any(grepl("Env_Local_Pressure_Station_XX", info_Child$children))) {
  #           child_id2 <- if (grepl("XX", info_Child$children)) gsub("XX", u, info_Child$children) else info_Child$children
  #           
  #           if (child_id2 %in% nodeList) {
  #             
  #             if(length(net$getChildren(child_id2)) == 0) net$addArc(child_id, child_id2)
  #             
  #           } else {
  #             
  #             # Create the node if it doesn't exist (it should exist depending on the order of nodes generation)
  #             # First fetch info
  #             info_Child2 <- getArchetypeSpec(net, info_Child$children)
  #             
  #             # Then create the node
  #             child2_handle <- createTemplateCptNode(net, 
  #                                                    id = child_id2, 
  #                                                    name = child_id2,
  #                                                    outcomes = info_Child2$outcomes,
  #                                                    temporal = info_Child2$temporal
  #             )
  #             
  #             # And add the arc
  #             net$addArc(child_id, child_id2)
  #             
  #           }
  #         }
  #       }
  #       
  #       # Deal with the link env pressure -> local pressure
  #       
  #     }
  #   }
  # }
  # 
  
  
  # replicateEnvNode <- function(net, archetype, newId) {
  #   spec <- archetypeSpec[[archetype]]
  #   
  #   createTemplateCptNode(
  #     net,
  #     id = newId,
  #     name = newId,
  #     outcomes = spec$outcomes,
  #     temporal = (spec$temporal == 3)
  #   )
  # }
  # # REGION WISE
  # # Coral_Diversity
  # for(node in RegionNodes[grep("RorcCoralRichness", RegionNodes)]) {
  #   net$addArc("Coral_Diversity", node)
  # }
  # 
  # # Branching_Coral_Cover
  # for(node in RegionNodes[grep("RorcHCB", RegionNodes)]) {
  #   net$addArc("Branching_Coral_Cover", node)
  # }
  # 
  # # Massive_Coral_Cover
  # for(node in RegionNodes[grep("RorcHCM", RegionNodes)]) {
  #   net$addArc("Massive_Coral_Cover", node)
  # }
  # 
  # # Tabular_Coral_Cover
  # for(node in RegionNodes[grep("RorcHCT", RegionNodes)]) {
  #   net$addArc("Tabular_Coral_Cover", node)
  # }
  # 
  # # Other_Forms_Coral_Cover
  # for(node in RegionNodes[grep("RorcHCO", RegionNodes)]) {
  #   net$addArc("Other_Forms_Coral_Cover", node)
  # }
  # 
  # # Soft_Coral_Cover
  # for(node in RegionNodes[grep("RorcSC", RegionNodes)]) {
  #   net$addArc("Soft_Coral_Cover", node)
  # }
  # 
  # # Dead_Coral_Cover
  # for(node in RegionNodes[grep("RorcDC", RegionNodes)]) {
  #   net$addArc("Dead_Coral_Cover", node)
  # }
  # 
  # # Rock_and_Slabs_Cover
  # for(node in RegionNodes[grep("RorcRC", RegionNodes)]) {
  #   net$addArc("Rock_and_Slabs_Cover", node)
  # }
  # 
  # # Debris_Cover
  # for(node in RegionNodes[grep("RorcRB", RegionNodes)]) {
  #   net$addArc("Debris_Cover", node)
  # }
  # 
  # # Sand_and_Silt_Cover
  # for(node in RegionNodes[grep("RorcSD_SI", RegionNodes)]) {
  #   net$addArc("Sand_and_Silt_Cover", node)
  # }
  # 
  # # Fish_Diversity
  # for(node in RegionNodes[grep("RorcFishRichness", RegionNodes)]) {
  #   net$addArc("Fish_Diversity", node)
  # }
  # 
  # # Fish_Herbivores_Biomass
  # for(node in RegionNodes[grep("RorcHerbivoreAbund", RegionNodes)]) {
  #   net$addArc("Fish_Herbivores_Biomass", node)
  # }
  # 
  # # Fish_Corallivores_Biomass
  # for(node in RegionNodes[grep("RorcCorallivoreAbund", RegionNodes)]) {
  #   net$addArc("Fish_Corallivores_Biomass", node)
  # }
  # 
  # # Fish_Carnivores_Biomass
  # for(node in RegionNodes[grep("RorcCarnivoreAbund", RegionNodes)]) {
  #   net$addArc("Fish_Carnivores_Biomass", node)
  # }
  # 
  # # Fish_Scraper_Biomass
  # for(node in RegionNodes[grep("RorcCorallivoreAbund", RegionNodes)]) {
  #   net$addArc("Fish_Scraper_Biomass", node)
  # }
  # 
  # # Invertebrate_Diversity
  # for(node in RegionNodes[grep("RorcInvRichness", RegionNodes)]) {
  #   net$addArc("Invertebrate_Diversity", node)
  # }
  # 
  # # Invertebrate_Abundance
  # for(node in RegionNodes[grep("RorcInvAbund", RegionNodes)]) {
  #   net$addArc("Invertebrate_Abundance", node)
  # }
  
  # Create "automatic" nodes from transects to related master nodes 
  # Read data
  
  # CSMP data
  # Read fish data
  # fishCsmp <- read.csv(file.path("Data","Processed","CSMP_GBR_Fish_Surveys_20182024.csv")) # , col_types = "text"
  # Read coral data
  # coralCsmp <- read.csv(file.path("Data","Processed","CSMP_GBR_Coral_PIT_20182024.csv")) # , col_types = "text"
  
  # RORC
  # Coral
  # coralRorc <- read.csv(file.path("Data", "Species","RORC_Coral_Data_hdbn.csv"))
  # coralRorc$Station <- gsub(" ","_", coralRorc$Station)
  # coralRorc$Site <- gsub(" ","_", coralRorc$Site)
  # # Fish
  # fishRorc <- read.csv(file.path("Data", "Species","RORC_Fish_Data_hdbn.csv"))
  # # Inv
  # invRorc <- read.csv(file.path("Data", "Species","RORC_Inv_Data_hdbn.csv"))
  # 
  # 
  # # Create the nodes
  # 
  # # Function to create CPT node :
  # createTemplateCptNode <- function(net, id, name, outcomes, xPos = NULL, yPos = NULL) {
  #   # Create blank node
  #   handle <- net$addNode(net$NodeType$CPT, id)
  #   # Set node name (here same as id but can be more "pretty" with spaces and parenthesis for ex)
  #   net$setNodeName(handle, name)
  #   
  #   # Add position (optional)
  #   if(!is.null(xPos) &&  !is.null(yPos)) {
  #     net$setNodePosition(handle, xPos, yPos, 85L, 55L)
  #   }
  #   
  #   # Rename the original two outcomes if standard outcome (which it should be)
  #   net$setOutcomeId(handle, 0, outcomes[1])
  #   net$setOutcomeId(handle, 1, outcomes[2])
  #   
  #   # add and names outcomes > 2
  #   if(length(outcomes) > 2) {
  #     for(outcome in outcomes[3:length(outcomes)]) {
  #       net$addOutcome(handle, outcome)
  #     }
  #     
  #   }
  #   
  #   return(handle)
  #   
  # }
  # 
  # 
  # # Function that works line by line.
  # # Creates the transect node and then if it doesnt exist, create the associated station, site, sector, region node.
  # # Then links the transect node hierarchically to the higher level nodes
  # # A function to apply to each dataset:
  # createTemplateHierarchicalNodes <- function(dataset, 
  #                                             net, 
  #                                             structure = c("Region","Sector","Site","Station","Sample"),
  #                                             vars = colnames(dataset)[!colnames(dataset) %in% structure]) {
  #   
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
    
  #   data <- unique(dataset[, structure])
  #   
  #   # Cache Node ids and append new nodes on each loop instead of reusing net$getAllNodeIds() because its too memory intensive
  #   nodeList <- net$getAllNodeIds()
  #   
  #   # For each variable:
  #   for (var in vars) {
  #     # var = vars[1]
  #     
  #     cat(var,"\n")
  #     
  #     # df <- data[!is.na(data[[var]]),]
  #     df <- data
  #     
  #     # Loop for each dataset line
  #     for (i in 1:nrow(df)) {
  #       # i = 1
  #       cat("row",i,"\n")
  #       
  #       # Remember : ids (and not names) need unique identifiers. Samples/transects are not UNIQUE ! found it the hard R-crashing-style way
  #       region <- paste0(var, "_", "Region_", df$Region[i])
  #       sector <- paste0(var, "_", "Sector_", df$Sector[i])
  #       site <- paste0(var, "_", "Site_", df$Site[i])
  #       station <- paste0(var, "_", "Station_", df$Station[i])
  #       sample <- paste0(var, "_","Station_", df$Station[i], "_Sample_", df$Sample[i])
  #       
  #       # Create nodes if they don't exist
  #       for (node in c(region, sector, site, station, sample)) {
  #         # node = region
  #         if (!(node %in% nodeList)) {
  #           cat("node : ", node,"\n")
  #           # Create the node
  #           createTemplateCptNode(net, gsub("-","_", node), node, c("High","Medium","Low"))
  #           # Add to node list
  #           nodeList <- c(nodeList, node)
  #         }
  #         
  #       }
  #       
  #       # Add hierarchy arcs (variable-specific) 
  #       # Check for preexisting links because repeated sampling of transects. 
  #       # We only have to (actually must to not raise an error) define arcs once for specified nodes
  #       # However in this loop, region-sector arcs up to stations are asked to be redefined each sample
  #       # So we check if it already exist we don't need to re-arc it
  #       if (!(region %in% net$getParentIds(sector))) {
  #         net$addArc(region, sector)
  #       }
  #       
  #       if (!(sector %in% net$getParentIds(site))) {
  #         net$addArc(sector, site)
  #       }
  #       
  #       if (!(site %in% net$getParentIds(station))) {
  #         net$addArc(site, station)
  #       }
  #       
  #       if (!(station %in% net$getParentIds(sample))) {
  #         net$addArc(station, sample)
  #       }
  #       
  #       # This one normally doesnt exist:
  #       # net$addArc(station, sample)
  #       
  #     }
  #     
  #   }
  #   
  # }
  # 
  # 
  # # Creat hierarchical nodes per variable per dataset
  # createTemplateHierarchicalNodes(dataset = coralRorc, 
  #                                 net = net, 
  #                                 structure = c("Region","Sector","Site","Station","Sample"),
  #                                 vars = c("RorcHCB", "RorcHCM", "RorcHCO", "RorcHCT", "RorcSC", "RorcDC", "RorcSD_SI", "RorcRB", "RorcRC", "RorcCoralRichness",
  #                                          "RorcCoralCover"))
  # 
  # length(net$getAllNodeIds())
  # 
  # createTemplateHierarchicalNodes(dataset = fishRorc, 
  #                                 net = net, 
  #                                 structure = c("Region","Sector","Site","Station","Sample"),
  #                                 vars = c("RorcFishRichness","RorcCarnivoreAbund","RorcCorallivoreAbund","RorcHerbivoreAbund"))
  # 
  # length(net$getAllNodeIds())
  # 
  # createTemplateHierarchicalNodes(dataset = invRorc, 
  #                                 net = net, 
  #                                 structure = c("Region","Sector","Site","Station","Sample"),
  #                                 vars = c("RorcInvRichness","RorcInvAbund"))
  # 
  # length(net$getAllNodeIds())
  # 
  # 
  # 
  # # Connect the appropriate data type region node to core model related node
  # nodes <- net$getAllNodeIds()
  # 
  # # Core nodes
  # # Coral:
  # # Coral_Diversity <- c("RorcCoralRichness")
  # nodes[1:25]
  # RegionNodes <- nodes[grep("Region", nodes)] 
  # RegionNodes
  # 
  # net$getParentIds("RorcCoralRichness_Region_PS")
  # net$getParentIds("Coral_Diversity")
  # net$getParentIds("Coral_Reef_Ecosystem_Productivity")
  # 
  # net$addArc
  # net$getMaxNodeTemporalOrder("Coral_Diversity")
  # net$getMaxNodeTemporalOrder("RorcCoralRichness_Region_PS")
  # 
  # # Coral_Diversity
  # for(node in RegionNodes[grep("RorcCoralRichness", RegionNodes)]) {
  #   net$addArc("Coral_Diversity", node)
  # }
  # 
  # # Branching_Coral_Cover
  # for(node in RegionNodes[grep("RorcHCB", RegionNodes)]) {
  #   net$addArc("Branching_Coral_Cover", node)
  # }
  # 
  # # Massive_Coral_Cover
  # for(node in RegionNodes[grep("RorcHCM", RegionNodes)]) {
  #   net$addArc("Massive_Coral_Cover", node)
  # }
  # 
  # # Tabular_Coral_Cover
  # for(node in RegionNodes[grep("RorcHCT", RegionNodes)]) {
  #   net$addArc("Tabular_Coral_Cover", node)
  # }
  # 
  # # Other_Forms_Coral_Cover
  # for(node in RegionNodes[grep("RorcHCO", RegionNodes)]) {
  #   net$addArc("Other_Forms_Coral_Cover", node)
  # }
  # 
  # # Soft_Coral_Cover
  # for(node in RegionNodes[grep("RorcSC", RegionNodes)]) {
  #   net$addArc("Soft_Coral_Cover", node)
  # }
  # 
  # # Dead_Coral_Cover
  # for(node in RegionNodes[grep("RorcDC", RegionNodes)]) {
  #   net$addArc("Dead_Coral_Cover", node)
  # }
  # 
  # # Rock_and_Slabs_Cover
  # for(node in RegionNodes[grep("RorcRC", RegionNodes)]) {
  #   net$addArc("Rock_and_Slabs_Cover", node)
  # }
  # 
  # # Debris_Cover
  # for(node in RegionNodes[grep("RorcRB", RegionNodes)]) {
  #   net$addArc("Debris_Cover", node)
  # }
  # 
  # # Sand_and_Silt_Cover
  # for(node in RegionNodes[grep("RorcSD_SI", RegionNodes)]) {
  #   net$addArc("Sand_and_Silt_Cover", node)
  # }
  # 
  # # Fish_Diversity
  # for(node in RegionNodes[grep("RorcFishRichness", RegionNodes)]) {
  #   net$addArc("Fish_Diversity", node)
  # }
  # 
  # # Fish_Herbivores_Biomass
  # for(node in RegionNodes[grep("RorcHerbivoreAbund", RegionNodes)]) {
  #   net$addArc("Fish_Herbivores_Biomass", node)
  # }
  # 
  # # Fish_Corallivores_Biomass
  # for(node in RegionNodes[grep("RorcCorallivoreAbund", RegionNodes)]) {
  #   net$addArc("Fish_Corallivores_Biomass", node)
  # }
  # 
  # # Fish_Carnivores_Biomass
  # for(node in RegionNodes[grep("RorcCarnivoreAbund", RegionNodes)]) {
  #   net$addArc("Fish_Carnivores_Biomass", node)
  # }
  # 
  # # Fish_Scraper_Biomass
  # for(node in RegionNodes[grep("RorcCorallivoreAbund", RegionNodes)]) {
  #   net$addArc("Fish_Scraper_Biomass", node)
  # }
  # 
  # # Invertebrate_Diversity
  # for(node in RegionNodes[grep("RorcInvRichness", RegionNodes)]) {
  #   net$addArc("Invertebrate_Diversity", node)
  # }
  # 
  # # Invertebrate_Abundance
  # for(node in RegionNodes[grep("RorcInvAbund", RegionNodes)]) {
  #   net$addArc("Invertebrate_Abundance", node)
  # }
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # # EXPORT THE MODEL
  # net$writeFile(file.path("Results", "Rorc_HDBN_Template.xdsl"))
  # 
  
  
  
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



