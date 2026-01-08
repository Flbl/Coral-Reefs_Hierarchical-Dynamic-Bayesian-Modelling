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
# COTS outbreak : if even one transect of each station has more than 100 COTS/ha then outbreak 


# INIT -----

  # Clean slate
  rm(list=ls(all=TRUE))
  
  # Libraries initialization
  init1 <- tail(unlist(strsplit(rstudioapi::getActiveDocumentContext()$path, "/")), n = 1)
  source("Scripts/00_Initialisation.R")


  # eo INIT ----


# TEMPERATURE ----
  
  # 5 states : well below, below, average, over, well over
  # We'll use the Coral Temp from the CRW datasets (5km res)
  # There could also be the MUR GHRSST but doesn't cover after 2023 (2km res)
  
  # Coral reef watch SST:
  # original link : https://coastwatch.noaa.gov/erddap/griddap/noaacrwsstDaily.html
  # Use this link to get Dataset ID and grid Variables wanted, then come back here to use the custom R function to request the data through URL
  
  
  download_erddap_yearly(
    dataset_id = "noaacrwsstDaily",
    vars = c("analysed_sst"),  
    start_year = 2013,
    end_year = 2024,
    lat_min = -27,
    lat_max = -14,
    lon_min = 155,
    lon_max = 175,
    out_dir = file.path(pathEnv, "Temperature")
  )
  
  
  

  
  
  
  
  
  
  
# HEATWAVES ----
  # Bleaching alert Area of the Coral Reef Watch
  # Simplifying the stress levels into three states:
  # No alert (first 1-3 levels)
  # Bleaching (Alert 1-2)
  # Mortality (Alert 3-5)
  nc <- nc_open(file.path(pathDat,"Environment","CRW","crw_baa_2013_2018.nc"))
  print(nc)
  names(nc$var)
  ncvar_get(nc, "bleaching_alert_area")
  
  
  
  nc <- nc_open(file.path(pathDat,"Environment","CRW","noaacrwdhwDaily_b396_7599_f8b5.nc"))
  nc <- nc_open(file.path(pathDat,"Environment","CRW","noaacrwdhwDaily_83e9_2084_6dff.nc"))
  
  
  # eo heatwaves ----
  
  

  
  
# COTS OUTBREAKS ----
  # If even one transect (20 meters * 5 meters) has ACA > 100/ha (>1 per transect of 100sqm)
  # See Dumas et al. 2020
  
  
  
  
  


# Trash ----------------------------
  
  # resp <- req_perform(req, path = outfile, verbosity = 3)
  # 
  # # Check HTTP status and report result
  # if (resp_status(resp) >= 400) {
  #   warning("Download failed for year ", year, " (HTTP ", resp_status(resp), ")")
  #   if (file.exists(outfile)) {
  #     file.remove(outfile)
  #   }
  # } else {
  #   message("Downloaded successfully.")
  # }
  
  
  # # Create the request object
  # req <- request(url)
  # 
  # # Add retry policy (exponential backoff)
  # # /!\ maybe doesnt work because "req_retry() in httr2 does not wrap streaming to disk (path) properly for retries)."
  # # So we try adding the for attempts above
  # # req <- req_retry(
  # #   req,
  # #   max_tries = retries
  # #   # backoff = ~ pause_sec * (2 ^ .x)  # 3s, 6s, 12s, etc.
  # #   # backoff = ~ pause_sec + .x
  # #   
  # # )
  # 
  # 
  # # Perform request and stream directly to disk
  # # Add a trycatch for 502 errors serverside
  # 
  # startTime <- Sys.time()
  # 
  # tryCatch(
  #   {
  #     message("Starting request...")
  #     resp <- req_perform(req, path = outfile, verbosity = 3)
  #     
  #     if (resp_status(resp) >= 400) {
  #       stop("HTTP failure: ", resp_status(resp))
  #     }
  #     
  #     elapsed <- round(as.numeric(Sys.time() - startTime, units = "secs"),1)
  #     message("Download successful (", elapsed, " s).")
  #   },
  #   
  #   error = function(e) {
  #     message("Error during download:", conditionMessage(e))
  #     
  #     if(file.exists(outfile)) {
  #       file.remove(outfile)
  #     }
  #   }
  #     
  # )

  
  # download_erddap_yearly <- function(
    #   dataset_id,
  #   var_names,
  #   start_year,
  #   end_year,
  #   lat_min,
  #   lat_max,
  #   lon_min,
  #   lon_max,
  #   out_dir = ".",
  #   retries = 5,
  #   pause_sec = 3
  # ) {
  #   
  #   stopifnot(is.character(dataset_id), length(dataset_id) == 1)
  #   stopifnot(is.character(var_names), length(var_names) >= 1)
  #   
  #   # Build base URL
  #   base_url <- paste0(
  #     "https://coastwatch.noaa.gov/erddap/griddap/",
  #     dataset_id
  #   )
  #   
  #   # ensure output directory
  #   dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  #   
  #   message("Dataset ID: ", dataset_id)
  #   message("Variables: ", paste(var_names, collapse = ", "))
  #   
  #   for (year in start_year:end_year) {
  #     
  #     start_date <- sprintf("%d-01-01T12:00:00Z", year)
  #     end_date   <- sprintf("%d-12-31T12:00:00Z", year)
  #     
  #     outfile <- file.path(out_dir, sprintf("%s_%d.nc", dataset_id, year))
  #     
  #     message("--------------------------------------------------")
  #     message("Year ", year)
  #     message("Output file: ", outfile)
  #     
  #     # build query portion variable-wise
  #     variable_queries <- vapply(
  #       var_names,
  #       function(v) sprintf(
  #         "%s[(%s):1:(%s)][(%s):1:(%s)][(%s):1:(%s)]",
  #         v,
  #         start_date, end_date,
  #         lat_min, lat_max,
  #         lon_min, lon_max
  #       ),
  #       character(1)
  #     )
  #     
  #     query <- paste(variable_queries, collapse = ",")
  #     
  #     url <- paste0(base_url, "?", query)
  #     
  #     # prepare request
  #     req <- request(url)
  #     
  #     # retry with exponential backoff
  #     req <- req_retry(
  #       req,
  #       max_tries = retries,
  #       backoff = ~ pause_sec * (2 ^ .x)
  #     )
  #     
  #     # save to disk
  #     resp <- req_perform(req, path = outfile)
  #     
  #     if (resp_status(resp) >= 400) {
  #       warning("Failed download for year ", year)
  #     } else {
  #       message("Downloaded successfully.")
  #     }
  #   }
  #   
  #   message("Completed all requested years.")
  # }
  # 