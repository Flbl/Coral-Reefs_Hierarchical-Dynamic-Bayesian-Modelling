#############################################################################################
#                                                                                           #
#                                                                                           #
########################      Booting all pipeline prerequisites          ###################
#                                                                                           #
#                                                                                           #
#############################################################################################



# Common paths

pathDat <- file.path("Data")
pathSpe <- file.path("Data","00_Raw","Species")
pathEnv <- file.path("Data","00_Raw", "Environment")
pathGra <- file.path("Data","00_Raw", "Graph")
pathPro <- file.path("Data", "01_Processed")
pathProSpe <- file.path("Data", "01_Processed", "Species")
pathProEnv <- file.path("Data", "01_Processed", "Environment")
pathProGra <- file.path("Data", "01_Processed", "Graph")
pathProCpt <- file.path("Data", "01_Processed", "CPT")
pathRes <- file.path("Results")

# Create dirs
dir.create(pathDat, showWarnings = FALSE)
dir.create(pathSpe, showWarnings = FALSE)
dir.create(pathEnv, showWarnings = FALSE)
dir.create(pathGra, showWarnings = FALSE)
dir.create(pathPro, showWarnings = FALSE)
dir.create(pathProSpe, showWarnings = FALSE)
dir.create(pathProEnv, showWarnings = FALSE)
dir.create(pathProGra, showWarnings = FALSE)
dir.create(pathProCpt, showWarnings = FALSE)
dir.create(pathRes, showWarnings = FALSE)


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


if(init1 == "02_Environment_Data_To_States.R") {
  
  
  if (!require("purrr")) install.packages("purrr")
  if (!require("lubridate")) install.packages("lubridate")
  if (!require("ncdf4")) install.packages("ncdf4")
  if (!require("terra")) install.packages("terra")
  if (!require("stars")) install.packages("stars")
  if (!require("httr2")) install.packages("httr2")
  if (!require("lwgeom")) install.packages("lwgeom")
  if (!require("geosphere")) install.packages("geosphere")
  
  
}


if(init1 == "03_HDBN_Template_Evidence.R") {
  
  if (!require("purrr")) install.packages("purrr")
  
}


if(init1 == "04_HDBN_Extract_CPT.R") {
  
  if (!require("purrr")) install.packages("purrr")
  
}


if(init1 == "05_HDBN_Fixed_Nodes.R") {
  
  if (!require("purrr")) install.packages("purrr")

}


if(init1 == "06_HDBN_EM_learning.R") {
  
  if (!require("purrr")) install.packages("purrr")

}

if(init1 == "07_HDBN_Fixed_EM_Eliciting.R") {
  
  if (!require("purrr")) install.packages("purrr")
  
}

if(init1 == "08_HDBN_Test.R") {
  
  if (!require("purrr")) install.packages("purrr")
  
}


if(init1 == "99_Shiny_CPT_manual_modification_App.R") {
  
  if (!require("purrr")) install.packages("purrr")
  
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
  if (length(parents_now) > 0) {
    names(parents_now) <- vapply(
      parents_now,
      net$getNodeId,
      character(1)
    )
  }
  
  # --- Temporal parents (previous slices) ---
  parents_temp <- list()
  
  for (lag in temporal_order) {
    
    # rSMILE throws an error if the node is not temporal
    res <- tryCatch(
      net$getTemporalParents(node, lag),
      error = function(e) NULL
    )
    
    if (is.null(res) || length(res) == 0) next
    
    # rSMILE can return multiple temporal parents for the same lag
    for (r in res) {
      parent_handle <- r$handle
      parent_name   <- paste0(r$id, "_t-", r$order)
      parents_temp[[parent_name]] <- parent_handle
    }
  }
  
  # merge all parents
  parents <- c(parents_now, parents_temp)
  
  # --- No parents: simple prior ---
  if (length(parents) == 0) {
    probs <- net$getNodeDefinition(node)
    return(
      data.frame(
        State = states,
        Probability = probs,
        row.names = NULL
      )
    )
  }
  
  # list of outcomes for all parents
  parent_states <- lapply(parents, function(p) net$getOutcomeIds(p))
  
  # build all combinations of parent states
  parent_grid <- expand.grid(parent_states, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  names(parent_grid) <- names(parents)
  
  # raw CPT definition from SMILE
  raw <- net$getNodeDefinition(node)
  
  # --- Reshape CPT ---
  cpt_matrix <- matrix(
    raw,
    ncol = n_states,
    byrow = TRUE
  )
  
  colnames(cpt_matrix) <- paste0("P_", states)
  
  cpt <- cbind(parent_grid, cpt_matrix)
  
  return(cpt)
}


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


# Function to get genie generated env node specifications ("archetype" is a kept artifact of previous function version)
# getNodeSpec <- function(net, archetype) {
#   list(
#     nodeName = archetype,
#     outcomes = net$getOutcomeIds(archetype),
#     parents  = net$getParentIds(archetype),
#     children = net$getChildIds(archetype),
#     temporal = net$getNodeTemporalType(archetype)
#   )
# }
archetype = "Structural_Complexity"
archetype = "Branching_Coral_Cover_Station_TREND_XX"
archetype = "Fish_Diversity"
temporal_order = c(1, 2)

getNodeSpec <- function(net, archetype, temporal_order = c(1, 2)) {
  
  safe_temporal <- function(expr) {
    tryCatch(expr, error = function(e) NULL)
  }
  
  # --- Base contemporaneous spec ---
  spec <- list(
    nodeName = archetype,
    nodeHandle = net$getNode(archetype),
    outcomes = net$getOutcomeIds(archetype),
    parents  = net$getParentIds(archetype),
    children = net$getChildIds(archetype),
    temporal = net$getNodeTemporalType(archetype)
  )
  
  # Helper: normalize TemporalInfo output to a list of objects
  normalize_temporal_info <- function(x) {
    if (is.null(x)) return(list())
    # Expected: list of TemporalInfo
    if (is.list(x) && length(x) > 0) return(x)
    # Fallback: a single TemporalInfo (rare wrapper variant)
    list(x)
  }
  
  # Helper: convert list of TemporalInfo to details df
  to_details_df <- function(info_list) {
    if (length(info_list) == 0) {
      return(data.frame(id = character(0), order = integer(0), handle = integer(0)))
    }
    out <- lapply(info_list, function(r) {
      data.frame(
        id     = r$id,
        order  = as.integer(r$order),
        handle = as.integer(r$handle),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, out)
  }
  
  # --- Temporal parents (API supports lag argument) ---
  tp_all <- list()
  for (lag in temporal_order) {
    resP <- safe_temporal(net$getTemporalParents(archetype, lag))
    resP <- normalize_temporal_info(resP)
    if (length(resP) == 0) next
    tp_all <- c(tp_all, resP)
  }
  
  temporalParents_details <- to_details_df(tp_all)
  
  # --- Temporal children (API returns all; filter by $order) ---
  tc_all <- safe_temporal(net$getTemporalChildren(archetype))
  tc_all <- normalize_temporal_info(tc_all)
  
  if (length(tc_all) > 0) {
    tc_all <- Filter(function(r) !is.null(r$order) && as.integer(r$order) %in% temporal_order, tc_all)
  }
  
  temporalChildren_details <- to_details_df(tc_all)
  
  # Optional: stable ordering
  if (nrow(temporalParents_details) > 0) {
    temporalParents_details <- temporalParents_details[order(temporalParents_details$order, temporalParents_details$id), ]
    rownames(temporalParents_details) <- NULL
  }
  if (nrow(temporalChildren_details) > 0) {
    temporalChildren_details <- temporalChildren_details[order(temporalChildren_details$order, temporalChildren_details$id), ]
    rownames(temporalChildren_details) <- NULL
  }
  
  # Human-readable labels (your convention)
  temporalParents_ids <- if (nrow(temporalParents_details) == 0) {
    character(0)
  } else {
    paste0(temporalParents_details$id, "_t-", temporalParents_details$order)
  }
  
  temporalChildren_ids <- if (nrow(temporalChildren_details) == 0) {
    character(0)
  } else {
    paste0(temporalChildren_details$id, "_t+", temporalChildren_details$order)
  }
  
  # Attach
  spec$temporalParents_details  <- temporalParents_details
  spec$temporalChildren_details <- temporalChildren_details
  spec$temporalParents_ids      <- temporalParents_ids
  spec$temporalChildren_ids     <- temporalChildren_ids
  
  return(spec)
}





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


# Coding directly in R a function that calls griddap and its OPeNDAP hyperslab protocol to request data off the NOAA ERDDAP data server
download_erddap_yearly <- function(
    dataset_id,            # ERDDAP dataset identifier, e.g., "noaacrwbaa7dDaily"
    vars,                  # character vector of variable names, e.g., c("bleaching_alert_area","mask")
    start_year,            # first year to download, integer
    end_year,              # last year to download, integer
    start_MMDD = "01-01",  # time MM-DD in case the dataset has different date (e.g. monthly)
    end_MMDD = "12-31",    # time MM-DD in case the dataset has different date (e.g. monthly)
    lat_min,               # minimum latitude bound
    lat_max,               # maximum latitude bound
    lon_min,               # minimum longitude bound
    lon_max,               # maximum longitude bound
    out_dir = ".",         # output directory
    retries = 10,           # number of request retries
    pause_sec = 2          # base pause time between retries
) {
  
  # debugzone
  # dataset_id = "noaacwNPPVIIRSSQchlaDaily"
  # vars = c("chlor_a")
  # start_year = 2013
  # end_year = 2024
  # start_MMDD = "01-01"
  # end_MMDD = "12-31"
  # lat_min = -27
  # lat_max = -14
  # lon_min = 155
  # lon_max = 175
  # out_dir = file.path(pathEnv, "Chlorophyll_a")
  # retries = 5
  # pause_sec = 3
  # eo debugzone
  
  # Base ERDDAP griddap URL (no dataset-specific name yet)
  base_url <- "https://coastwatch.noaa.gov/erddap/griddap/"
  
  # Ensure output directory exists
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Loop through each year from start_year to end_year
  for (year in start_year:end_year) {
    
    # debugzone 
    # year = 2013
    # eo debugzone
    
    # Build ISO time stamps for the start and end of the year
    start_date <- sprintf("%d-%sT12:00:00Z", year, start_MMDD)
    end_date   <- sprintf("%d-%sT12:00:00Z", year, end_MMDD)
    # Build ISO time stamps for the start and end of the year
    # start_date <- sprintf("%d-01-01T12:00:00Z", year)
    # end_date   <- sprintf("%d-12-31T12:00:00Z", year)
    
    # Define output filename for this year's NetCDF file
    outfile <- file.path(out_dir, sprintf("%s_%d.nc", dataset_id, year))
    
    # Create temporary file name for download (so, among others, if there's an error while its not finished its not considered the final file)
    outfileTemp <- file.path(out_dir, sprintf("%s_%d_Temp.nc", dataset_id, year))
    
    message("--------------------------------------------------\nDownloading:")
    message("Dataset: ", dataset_id)
    message("Year: ", year)
    message("Output file: ", outfile)
    
    # Check for existing files to not redownload the same data if some years failed
    if (file.exists(outfile)) {
      message("File already exists for ", year, "; skipping.")
      next
    }
    
    # ---- Build variable query part ----
    # For each variable requested, create a full ERDDAP query slice
    # Structure: var[(time1):1:(time2)][(latmin):1:(latmax)][(lonmin):1:(lonmax)]
    var_query_parts <- vapply(
      vars,
      FUN = function(v) {
        sprintf(
          "%s[(%s):1:(%s)][(%s):1:(%s)][(%s):1:(%s)]",
          v,
          start_date, end_date,
          lat_min, lat_max,
          lon_min, lon_max
        )
      },
      FUN.VALUE = character(1)
    )
    
    # Collapse multiple variable query components into one comma-separated string
    var_query <- paste(var_query_parts, collapse = ",")
    
    # ---- Build full download URL ----
    # Add dataset ID and .nc output extension
    url <- paste0(
      base_url,
      dataset_id,
      ".nc?",
      var_query
    )
    
    # Display the URL if debugging is needed
    # message("Request URL: ", url)
    
    # Adding attempts loop since potentially "req_retry() in httr2 does not wrap streaming to disk (path) properly for retries)."
    for (attempt in 1:retries) {
      # debugzone
      # attempt = 1
      # eo debugzone
      
      startTime <- Sys.time()
      message("Try ", attempt)
      
      # Create the request object
      req <- request(url)
      
      tryCatch(
        {
          message("Requesting data...")
          resp <- req_perform(req, path = outfileTemp, verbosity = 0)
          
          if (resp_status(resp) >= 400) {
            stop("HTTP failure: ", resp_status(resp))
          }
          
          file.rename(outfileTemp,outfile)
          elapsed <- round(as.numeric(Sys.time() - startTime, units = "secs"),1)
          message("Download successful (", elapsed, " s).")
          break
        },
        
        error = function(e) {
          message("Error:", conditionMessage(e))
          
          if(file.exists(outfile)) {
            file.remove(outfile)
          }
          if(attempt < retries) {
            wait <- pause_sec + jitter(attempt)
            message("Retrying in around ",round(wait), " seconds...")
            Sys.sleep(wait)
          } else {
            message("All retries failed for year ",year)
          }
        }
      )
    }
  }
  
  message("All requested years processed.")
}


# Replot extents of a map in the graphical RSTUDIO HUD
replot <- function(x) {
  plot(x)
  pts <- locator(n = 2)
  xlimits <- range(pts$x) #BOTTOM LEFT
  ylimits <- range(pts$y) #TOP RIGHT
  plot(x, xlim = xlimits, ylim = ylimits)
  
}



# Trash
# old getCPfunction
# getCPT <- function(net, node, temporal_order = c(1,2)) {
#   states <- net$getOutcomeIds(node)
#   n_states <- length(states)
#   
#   # parents in the same time slice
#   parents_now <- net$getParents(node)
#   names(parents_now) <- unlist(lapply(parents_now, net$getNodeId))
#   
#   # parents from previous slices (temporal)
#   parents_temp <- unlist(lapply(temporal_order, function(x) {
#     # x = 2
#     res <- net$getTemporalParents(node, x)
#     if(length(res) == 0) return(NULL)
#     resNode <- res[[1]]$handle
#     resId <- res[[1]]$id
#     resOrder <- res[[1]]$order
#     names(resNode) <- paste0(resId, "_t-",resOrder)
#     return(resNode)
#     
#   }))
#   
#   # merge all parents
#   parents <- c(parents_now, parents_temp)
#   
#   # no parents -> simple prior
#   if (length(parents) == 0) {
#     probs <- net$getNodeDefinition(node)
#     return(data.frame(State = states, Probability = probs))
#   }
#   
#   # list of outcomes for all parents
#   parent_states <- lapply(parents, function(p) net$getOutcomeIds(p))
#   
#   # build all combinations of parent states
#   parent_grid <- expand.grid(parent_states, KEEP.OUT.ATTRS = FALSE)
#   # names(parent_grid) <- names(parents)
#   
#   # raw CPT definition from SMILE
#   raw <- net$getNodeDefinition(node)
#   
#   # reshape: each parent combo gets n_states probabilities
#   cpt <- cbind(parent_grid, matrix(raw, ncol = n_states, byrow = TRUE))
#   names(cpt)[(ncol(cpt)-n_states+1):ncol(cpt)] <- states
#   
#   return(cpt)
# }

