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




# Corrected function to get CPT in a datagrame from a graph network with temporal transitions
# Extract DBN CPTs from an rSMILE Network for ALL temporal definitions of a node
# Why this is needed:
# In SMILE/GeNIe DBNs, a plate node with temporal arcs needs multiple CPT "definitions":
#  - t = 0 uses net$getNodeDefinition(node)
#  - t = k (k >= 1) uses net$getNodeTemporalDefinition(node, k)
# The set of parents indexing each definition is different in early slices, so we must use
# net$getUnrolledParents(node, k) to get the correct parent list for definition k.
#
# Plus: The order of CPT dimensions follows the order in which arcs to the node were created.
# So we need to reconstruct the CPT table by iterating over the linear CPT vector in SMILE order,
# Output:
# A named list of data.frames:
#   $t0, $t1, ..., and the last one labeled like "t2plus" (if max order is 2),
# each with parent state columns + P_<outcome> columns.

# For decoding each index into multi-dimensional coordinates, we use a helper function (given in the wrapper pdf)
indexToCoords <- function(index, dimSizes) {
  prod <- 1L
  coords <- integer(length(dimSizes))
  for (i in length(dimSizes):1) {
    coords[i] <- floor(index / prod) %% dimSizes[[i]]
    prod <- prod * dimSizes[[i]]
  }
  coords
}

# Then we can make a function to get a temporal aware, ordered parent list 
getParentsForDefinition <- function(net, node, k) {
  
  # Returns list(handles=..., labels=...) where labels encode temporal order
  if (k == 0) {
    info <- tryCatch(net$getUnrolledParents(node, 0), error=function(e) NULL)
  } else {
    info <- tryCatch(net$getUnrolledParents(node, k), error=function(e) NULL)
  }
  
  if (is.null(info)) {
    if (k == 0) {
      hs <- net$getParents(node)
      labs <- vapply(hs, net$getNodeId, character(1))
      return(list(handles = hs, labels = labs))
    }
    stop("getUnrolledParents failed for node=", node, " k=", k)
  }
  
  handles <- vapply(info, function(x) as.integer(x$handle), integer(1))
  labels  <- vapply(info, function(x) {
    pid <- as.character(x$id)
    ord <- as.integer(x$order)
    if (ord == 0) pid else paste0(pid, "_t-", ord)
  }, character(1))
  
  # IMPORTANT: keep ordering exactly as returned; ensure uniqueness just in case
  labels <- make.unique(labels)
  
  list(handles = handles, labels = labels)
}

# Helper function to get definition
getDefinitionByK <- function(net, node, k) {
  if (k == 0) {
    return(net$getNodeDefinition(node))
  }
  # k >= 1
  return(net$getNodeTemporalDefinition(node, k))
}

# And the function to extract for one temporal definition 
getCPT_one_definition <- function(net, node, k = 0) {
  
  # Parents in the exact order SMILE reports for this temporal definition
  par <- getParentsForDefinition(net, node, k)
  parent_handles <- par$handles
  parent_labels  <- par$labels
  
  child_states <- net$getOutcomeIds(node)
  n_child <- length(child_states)
  
  raw <- getDefinitionByK(net, node, k)
  
  parent_sizes <- if (length(parent_handles) > 0) {
    vapply(parent_handles, net$getOutcomeCount, integer(1))
  } else integer(0)
  
  dimSizes <- c(parent_sizes, n_child)
  expected_len <- prod(dimSizes)
  
  if (length(raw) != expected_len) {
    stop(
      "Definition length mismatch for node=", node, " k=", k,
      " expected=", expected_len, " got=", length(raw),
      " (parent sizes: ", paste(parent_sizes, collapse="x"),
      ", child=", n_child, ")"
    )
  }
  
  # Number of parent configurations (each has n_child probs)
  n_parent_configs <- if (length(parent_sizes) > 0) prod(parent_sizes) else 1L
  
  parent_rows <- vector("list", n_parent_configs)
  prob_mat <- matrix(NA_real_, nrow = n_parent_configs, ncol = n_child)
  
  for (pc in seq_len(n_parent_configs)) {
    # pc = 1 
      
    base_idx0 <- (pc - 1L) * n_child  # 0-based start index for this parent config block
    coords <- indexToCoords(base_idx0, dimSizes)
    
    if (length(parent_handles) > 0) {
      parent_state_ids <- vapply(seq_along(parent_handles), function(j) {
        net$getOutcomeId(parent_handles[j], coords[j])
      }, character(1))
      names(parent_state_ids) <- parent_labels
      parent_rows[[pc]] <- as.list(parent_state_ids)
    } else {
      parent_rows[[pc]] <- list()
    }
    
    prob_mat[pc, ] <- raw[(base_idx0 + 1L):(base_idx0 + n_child)]
  }
  
  parent_df <- if (length(parent_labels) > 0) {
    as.data.frame(do.call(rbind, lapply(parent_rows, as.data.frame)), stringsAsFactors = FALSE)
  } else {
    # IMPORTANT: 1-row empty df so cbind() matches prob_mat rows
    data.frame(row_id = 1, stringsAsFactors = FALSE)[, 0, drop = FALSE]
  }
  
  colnames(prob_mat) <- paste0("P_", child_states)
  
  cbind(parent_df, as.data.frame(prob_mat, check.names = FALSE))
}

# And the final wrapper:
getCPT <- function(net, node, max_k = 2) {
  out <- list()
  
  for (k in 0:max_k) {
    res <- tryCatch(
      getCPT_one_definition(net, node, k = k),
      error = function(e) e
    )
    
    if (inherits(res, "error")) {
      # k==0 should never fail; if it does, propagate
      if (k == 0) stop(res)
      
      # for k>0: stop at first missing temporal definition
      break
    }
    
    out[[paste0("t", k)]] <- res
  }
  
  out
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
# archetype = "Structural_Complexity"
# archetype = "Branching_Coral_Cover_Station_TREND_XX"
# archetype = "Fish_Diversity"
# temporal_order = c(1, 2)




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
# Function to get CPT in dataframe from a graph network from the rSMILE package
# THIS ONE DOESNT CARE FOR TRANSITION CPTS AND RECYCLES FOR TEMPORAL PARENTS
# getCPT <- function(net, node, temporal_order = c(1,2)) {
#   # node = "Coral_Reef_Ecosystem_Health"
#   # temporal_order = c(1,2)
#   
#   states <- net$getOutcomeIds(node)
#   n_states <- length(states)
#   
#   # parents in the same time slice
#   parents_now <- net$getParents(node)
#   if (length(parents_now) > 0) {
#     names(parents_now) <- vapply(
#       parents_now,
#       net$getNodeId,
#       character(1)
#     )
#   }
#   
#   # --- Temporal parents (previous slices) ---
#   parents_temp <- list()
#   
#   for (lag in temporal_order) {
#     
#     # rSMILE throws an error if the node is not temporal
#     res <- tryCatch(
#       net$getTemporalParents(node, lag),
#       error = function(e) NULL
#     )
#     
#     if (is.null(res) || length(res) == 0) next
#     
#     # rSMILE can return multiple temporal parents for the same lag
#     for (r in res) {
#       parent_handle <- r$handle
#       parent_name   <- paste0(r$id, "_t-", r$order)
#       parents_temp[[parent_name]] <- parent_handle
#     }
#   }
#   
#   # merge all parents
#   parents <- c(parents_now, parents_temp)
#   
#   # --- No parents: simple prior ---
#   if (length(parents) == 0) {
#     probs <- net$getNodeDefinition(node)
#     return(
#       data.frame(
#         State = states,
#         Probability = probs,
#         row.names = NULL
#       )
#     )
#   }
#   
#   # list of outcomes for all parents
#   parent_states <- lapply(parents, function(p) net$getOutcomeIds(p))
#   
#   # build all combinations of parent states
#   parent_grid <- expand.grid(parent_states, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
#   names(parent_grid) <- names(parents)
#   
#   # raw CPT definition from SMILE
#   raw <- net$getNodeDefinition(node)
#   
#   # --- Reshape CPT ---
#   cpt_matrix <- matrix(
#     raw,
#     ncol = n_states,
#     byrow = TRUE
#   )
#   
#   colnames(cpt_matrix) <- paste0("P_", states)
#   
#   cpt <- cbind(parent_grid, cpt_matrix)
#   
#   return(cpt)
# }
# 
# getCPT <- function(net, node, max_order = NULL, max_order_search = 10) {
#   
#   # Test zone
#   # node = "Coral_Reef_Ecosystem_Health"
#   # node = "Env_Environmental_Shock"
#   # max_order = NULL
#   # max_order_search = 10
#   # eo tz
#   
#   # --- node outcomes ---
#   states <- net$getOutcomeIds(node)
#   n_states <- length(states)
#   
#   # --- infer maximum temporal order if not provided ---
#   if (is.null(max_order)) {
#     found <- integer(0)
#     for (k in 1:max_order_search) {
#       res <- tryCatch(net$getTemporalParents(node, k), error = function(e) NULL)
#       if (!is.null(res) && length(res) > 0) found <- c(found, k)
#     }
#     max_order <- if (length(found) == 0) 0L else max(found)
#   }
#   
#   # helper to label definitions nicely
#   def_name <- function(k, kmax) {
#     if (k == 0) return("t0")
#     if (k == kmax && kmax >= 2) return(paste0("t", k, "plus")) # applies to t>=k
#     paste0("t", k)
#   }
#   
#   out <- vector("list", length = max_order + 1L)
#   names(out) <- vapply(0:max_order, def_name, character(1), kmax = max_order)
#   
#   for (k in 0:max_order) {
#     # k = 1
#     
#     # --- get the correct parent list for THIS temporal definition ---
#     # rSMILE returns a list of TemporalInfo objects with fields: handle, id, order
#     parent_info <- tryCatch(net$getUnrolledParents(node, k), error = function(e) NULL)
#     
#     # getUnrolledParents doesn't work for initial slice
#     if (is.null(parent_info)) {
#       if (k == 0) {
#         # No temporal information here; just contemporaneous parents
#         parent_handles <- net$getParents(node)
#         parent_labels  <- vapply(parent_handles, net$getNodeId, character(1))
#       } else {
#         stop("getUnrolledParents() failed for node=", node, " order=", k,
#              ". This function requires SMILE wrappers that support temporal definitions.")
#       }
#     } else {
#       # Convert TemporalInfo list to handles + labels (labels must encode temporal order)
#       parent_handles <- vapply(parent_info, function(x) as.integer(x$handle), integer(1))
#       
#       parent_labels <- vapply(parent_info, function(x) {
#         pid <- as.character(x$id)
#         ord <- as.integer(x$order)
#         if (ord == 0) pid else paste0(pid, "_t-", ord)
#       }, character(1))
#       
#       # Ensure uniqueness if something repeats (rare but safe)
#       parent_labels <- make.unique(parent_labels)
#     }
#     
#     # --- get the correct numeric definition vector ---
#     raw <- if (k == 0) {
#       net$getNodeDefinition(node)
#     } else {
#       net$getNodeTemporalDefinition(node, k)
#     }
#     
#     # --- if no parents in this definition: simple prior distribution ---
#     if (length(parent_handles) == 0) {
#       if (length(raw) != n_states) {
#         stop("Node=", node, " order=", k,
#              ": expected ", n_states, " probs (no parents) but got ", length(raw))
#       }
#       out[[def_name(k, max_order)]] <- data.frame(
#         State = states,
#         Probability = raw,
#         row.names = NULL
#       )
#       next
#     }
#     
#     # --- build parent state grid in the correct parent order ---
#     parent_states <- lapply(parent_handles, function(h) net$getOutcomeIds(h))
#     names(parent_states) <- parent_labels
#     
#     parent_grid <- expand.grid(parent_states, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
#     # expand.grid already keeps names(parent_states), so this is optional:
#     names(parent_grid) <- parent_labels
#     
#     n_parent_configs <- nrow(parent_grid)
#     expected_len <- n_parent_configs * n_states
#     
#     if (length(raw) != expected_len) {
#       stop(
#         "Node=", node, " temporalDef=", k, " (", def_name(k, max_order), ")",
#         ": definition length mismatch.\n",
#         "Parents in this definition: ", paste(parent_labels, collapse = ", "), "\n",
#         "Parent configs=", n_parent_configs, ", outcomes=", n_states,
#         " => expected ", expected_len, " numbers, got ", length(raw), ".\n",
#         "This usually means you're reading the wrong temporal definition for the chosen parent set."
#       )
#     }
#     
#     # SMILE stores CPT rows sequentially by parent configuration, each row contains all child outcomes.
#     cpt_matrix <- matrix(raw, ncol = n_states, byrow = TRUE)
#     colnames(cpt_matrix) <- paste0("P_", states)
#     
#     out[[def_name(k, max_order)]] <- cbind(parent_grid, cpt_matrix)
#   }
#   
#   out
# }

