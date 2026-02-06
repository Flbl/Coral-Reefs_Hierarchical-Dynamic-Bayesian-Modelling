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
  # net$readFile(file.path(pathGra,"Hierarchical_Dynamic_Bayesian_Network_CoralReefs_RORC_Template.xdsl"))    
  # net$readFile(file.path(pathGra,"Hierarchical_Dynamic_Bayesian_Network_CoralReefs_RORC_Template_Reduced.xdsl"))    
  net$readFile(file.path(pathGra,"Hierarchical_Dynamic_Bayesian_Network_CoralReefs_RORC_Template_noOBSTREND.xdsl"))    
  
  # Get CPT list
  nodes <- net$getAllNodeIds()
  
  
  # Fixing a node means “do not let EM reinterpret this variable."
  # Observed nodes are not protected from learning. However, on directly observed environmental nodes (e.g. Cyclone Season), 
  # its meaningless to try and set CPT because we are not trying to model the stochastic processes that generates environmental forcings.
  # Environmental 'encoder' nodes like Cyclone Impact, Env Site conditions... are to be fixed
  # They are not observed environmental forcing with designed mappings
  # For those, the CPT defines the concept and they must be set intentionally and fixed.
  # Some biological nodes are also encoder nodes (Phase shift window, coral pressure index, Intrinsic Resilience)
  
  # toFix <- grep("Env", nodes, value = TRUE)
  toFix <- c("Env_Site_Conditions","Env_Chronic_Environmental_Stress","Env_Environmental_Shock","Env_Cyclone_Impact","Env_Acute_Local_Pressure")
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
              "Coral_Reef_Ecosystem_Diversity",
              "Fish_Grazing",
              "Urchin_Grazing",
              "Grazing_pressure",
              "Recovery_Capacity")
  

  
  
  source("Scripts/99_Shiny_CPT_manual_modification_App.R")
  
  run_cpt_app(
    pathProCpt = pathProCpt,
    folder = "05_CPT",
    cpt_patterns = c(toFix, toSoft)
  )

  
  # Updating OBS and TREND nodes as encoder nodes to fix them for learning
  obsTrendNodes <- grep("OBS|TREND", nodes, value = TRUE)
  
  run_cpt_app(
    pathProCpt = pathProCpt,
    folder = "05_CPT",
    cpt_patterns = c(obsTrendNodes)
  )

  
# UPDATING MODEL CPTs
  
  # What this does:
  # Reads CPT_*.csv files from a folder
  # For each target node, validates the CSV structure against the network:
  #     * node exists
  #     * probability columns match node outcomes (P_<state>)
  #     * parent columns match node parents (including temporal parents, e.g. _t-1)
  #     * row order matches the SMILE CPT ordering (or gets reordered to match)
  # Converts the CPT table into SMILE's flat probability vector
  # Writes it into the network using net$setNodeDefinition()
  
  # /!\ Assumptions:
  # The CPT CSVs were produced with getCPT() function SO COLUMN NAMES ARE UNCHANGED 
  # Probability column names are P_<state> WHICH EXACTLY MATCH net$getOutcomeIds(node) since its the same structure from getCPT
  # Parent columns are exactly as in getCPT(): current parents by id, temporal as "<id>_t-<lag>"
  
  cptPath <- file.path(pathProCpt,"05_CPT")
  # list.files(cptPath)
  
  test <- read.csv(file.path(cptPath, "CPT_Env_Acute_Local_Pressure_t1.csv"))
  net$getOutcomeIds("Env_Acute_Local_Pressure")
  getCPT(net,"Env_Acute_Local_Pressure")
  getCPT(net,"Env_Environmental_Shock")
  
  # net
  # node = "Acute_Local_Pressure"
  # k = 0
  
  net <- updateNetFromCPTFolder(net, cptPath, verbose = TRUE)
  
  
  
  # Save the new model
  # ""
  # Save the updated network
  # net$writeFile(file.path(pathProGra, "Rorc_HDBN_Template_Soft_Constrained.xdsl"))
  # net$writeFile(file.path(pathProGra, "Rorc_HDBN_Template_Soft_Constrained_Reduced.xdsl"))
  net$writeFile(file.path(pathProGra, "Rorc_HDBN_Template_Soft_Constrained_noOBSTREND.xdsl"))
  
  
# Trash ------
  
  # # Update CPTs in a rSMILE Network from CPT_*.csv files.
  # # - By default updates all nodes for which a CPT CSV exists.
  # # - Optionally restrict to a subset of node IDs.
  # update_cpts_from_csv <- function(net, cptPath, nodeIds = NULL, verbose = TRUE) {
  #   nodeIds = "Coral_Reef_Ecosystem_Health"
  #   
  #   
  #   # 1) list CPT files
  #   files <- list.files(cptPath, pattern = "^CPT_.*\\.csv$", full.names = TRUE)
  #   if (length(files) == 0) stop("No CPT_*.csv files found in: ", cptPath)
  #   
  #   # 2) convert filenames -> nodeIds (CPT_<nodeId>.csv)
  #   file_node <- function(fp) {
  #     x <- basename(fp)
  #     x <- sub("^CPT_", "", x)
  #     x <- sub("\\.csv$", "", x)
  #     x
  #   }
  #   node_from_files <- vapply(files, file_node, character(1))
  #   
  #   # 3) optionally keep only user requested nodes
  #   if (!is.null(nodeIds)) {
  #     keep <- node_from_files %in% nodeIds
  #     files <- files[keep]
  #     node_from_files <- node_from_files[keep]
  #     if (length(files) == 0) stop("None of the requested nodeIds had a CPT file.")
  #   }
  #   
  #   # 4) check nodes exist in the network
  #   net_nodes <- net$getAllNodeIds()
  #   missing <- setdiff(node_from_files, net_nodes)
  #   if (length(missing) > 0) {
  #     stop("These CPT files refer to nodes not in the network: ", paste(missing, collapse = ", "))
  #   }
  #   
  #   # 5) update each node
  #   report <- data.frame(nodeId = node_from_files, file = basename(files), updated = FALSE)
  #   
  #   for (i in seq_along(files)) {
  #     # i = 1
  #     nodeId <- node_from_files[i]
  #     fp <- files[i]
  #     
  #     if (verbose) message("Updating: ", nodeId, "  <-  ", basename(fp))
  #     
  #     # read CSV. check.names=FALSE prevents "t.1" name mangling, but we don't rely on those columns anyway.
  #     df <- read.csv(fp, stringsAsFactors = FALSE, check.names = FALSE)
  #     
  #     # find probability columns (P_<state>)
  #     prob_cols <- grep("^P_", names(df), value = TRUE)
  #     if (length(prob_cols) == 0) stop("No P_* columns in: ", basename(fp))
  #     
  #     # Keep probability columns in the SAME order as node outcomes
  #     states <- net$getOutcomeIds(nodeId)
  #     expected_prob_cols <- paste0("P_", states)
  #     
  #     if (!all(expected_prob_cols %in% prob_cols)) {
  #       stop("Missing some probability columns for node ", nodeId,
  #            ". Expected: ", paste(expected_prob_cols, collapse = ", "))
  #     }
  #     
  #     probs_mat <- as.matrix(df[, expected_prob_cols, drop = FALSE])
  #     
  #     # validate rows sum to 1 (light check)
  #     rs <- rowSums(probs_mat)
  #     if (any(is.na(rs)) || any(abs(rs - 1) > 1e-6)) {
  #       bad <- which(abs(rs - 1) > 1e-6)[1]
  #       stop("Row probabilities do not sum to 1 in ", basename(fp),
  #            " (node ", nodeId, "), bad row: ", bad, ", sum=", rs[bad])
  #     }
  #     
  #     # Flatten into SMILE definition vector: row1 all states, row2 all states, ...
  #     # SMILE expects: P(state1|parents=row1),...,P(stateK|row1), P(state1|row2),..., etc.
  #     definition <- as.numeric(t(probs_mat))
  #     
  #     nodeHandle <- net$getNode(nodeId)
  #     net$setNodeDefinition(nodeHandle, definition)
  #     
  #     report$updated[i] <- TRUE
  #   }
  #   
  #   if (verbose) message("Done. Updated ", sum(report$updated), " CPT(s).")
  #   return(list(net = net, report = report))
  # }
  # 
  # # Run the update
  # out <- update_cpts_from_csv(net, cptPath, nodeIds = c(toFix, toSoft))
  # net <- out$net
  # out$report
  # 
  
  