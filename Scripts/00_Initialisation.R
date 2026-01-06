#############################################################################################
#                                                                                           #
#                                                                                           #
########################      Booting all pipeline prerequisites          ###################
#                                                                                           #
#                                                                                           #
#############################################################################################



# Common paths

pathDat <- file.path("Data")
pathSpe <- file.path("Data", "Species")
pathEnv <- file.path("Data", "Environment")
pathPro <- file.path("Data", "Processed")
pathRes <- file.path("Results")

# Create dirs
dir.create("Data", showWarnings = FALSE)
dir.create(file.path("Data","Processed"), showWarnings = FALSE)
dir.create(file.path("Data","Environment"), showWarnings = FALSE)
dir.create(file.path("Data","Species"), showWarnings = FALSE)


# Load packages ----

# Common Packages

# if (!require("reshape2")) install.packages("reshape2")
# if (!require("plyr")) install.packages("plyr")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tibble")) install.packages("tibble")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("ggnewscale")) install.packages("ggnewscale")
# if (!require("ggrepel")) install.packages("ggrepel")
# if (!require("ggsignif")) install.packages("ggsignif")


# Specific packages for specific scripts

# 01_State_data_conversion 

if(init1 == "01_Data_To_States.R") {
  
  
  # Travel time
  # This package will install pretty much all other needed packages for SIG processing
  if (!require("brms")) install.packages("brms")
  if (!require("mclust")) install.packages("mclust")
  if (!require("purrr")) install.packages("purrr")
  if (!require("lmPerm")) install.packages("lmPerm")
  # if (!require("traveltime")) install.packages("traveltime", repos = c("https://idem-lab.r-universe.dev"))
  # if (!require("mapsapi")) install.packages("mapsapi")
  # if (!require("rSDM")) install.packages("rSDM", repos = c("https://pakillo.r-universe.dev", "https://cloud.r-project.org"))
  # if (!require("concaveman")) install.packages("concaveman")
  # if (!require("nngeo")) install.packages("nngeo")
  
}

if(init1 == "02_HDBN_Template_CPT_Creation.R") {
  
  
  # Travel time
  # This package will install pretty much all other needed packages for SIG processing
  if (!require("brms")) install.packages("brms")
  if (!require("mclust")) install.packages("mclust")
  if (!require("purrr")) install.packages("purrr")
  if (!require("lmPerm")) install.packages("lmPerm")

  
}



# FUNCTIONS

safe_Mclust <- function(x, ...) {
  
  out <- tryCatch(
    {
      mclust::Mclust(x, ...)
    },
    error = function(e) NULL,
    warning = function(w) NULL
  )
  
  return(out)
  
}


# Function to get CPT in dataframe from a graph network from the rSMILE package
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



