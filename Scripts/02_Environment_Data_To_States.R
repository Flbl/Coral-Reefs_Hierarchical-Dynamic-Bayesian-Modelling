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

  
# GET DATA -----
  
  # GET TEMPERATURE ----
  
  # We'll use the Coral Temp from the CRW datasets (5km res)
  # There could also be the MUR GHRSST but doesn't cover after 2023 (2km res)
  # Coral reef watch SST:
  # original link : https://coastwatch.noaa.gov/erddap/griddap/noaacrwsstDaily.html
  # Use this link to get Dataset ID and grid Variables wanted, 
  # then come back here to use the custom R function to request the data through URL (see 00_initialisation.R for the function)
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
  
    # eo get temperature ----
  
  
  # GET HEATWAVES ----

  # Bleaching alert Area of the Coral Reef Watch
  # Same function from ERRDAP, go there:
  # https://coastwatch.noaa.gov/erddap/griddap/noaacrwbaa7dDaily.html

  download_erddap_yearly(
    dataset_id = "noaacrwbaa7dDaily",
    vars = c("bleaching_alert_area"),  
    start_year = 2013,
    end_year = 2024,
    lat_min = -27,
    lat_max = -14,
    lon_min = 155,
    lon_max = 175,
    out_dir = file.path(pathEnv, "Heatwaves_BAA")
  )
    
    # eo get heatwaves ----


  # GET CHL-A ----

  # There is no gap filled/complete chla dataset in NOAA
  # So we use a copernicus dataset
  # Data was downloaded using the Copernicus Marine Toolbox (CLI-Command Line Interface)
  # See the readme

    # eo get chl-a


  # GET GEOMORPHOLOGY ----
  # Dataset from Andrefouet et al.2009 found on Georep.nc
  
    # eo geomorphology ----


  # GET CYCLONES ----

  # This one's a bit tricky
  # There is a global dataset of cyclones: the IBTrACS — Global best-track cyclone archive
  # https://www.ncei.noaa.gov/products/international-best-track-archive
  # There we can download either netcdf or csv or shapefile of the tracks of the cyclone
  # It also has data associated to the track, like type of depression and notably R34, R50 and R64 radii (NE, SE, SW, NW) which represent the sustained wind radius per quadrant for each track point
  # R34 = radius of sustained wind of 34kts
  # This is a common use to determine the intensity of the storm at a given point and its influence
  # "These wind radii represent the spatial footprint of wind hazard, roughly analogous to empirical influence distance metrics."
  # https://link.springer.com/article/10.1007/s13143-022-00274-5?utm_source=chatgpt.com  
  # So we'll build the polygon of influence of each cyclone and cross them with the stations point per year to determine if they lived the cyclone or not

    # eo get cyclones ----


  # GET TRAVEL TIME / GRAVITY ----

  # Done in a separated project. Data available on demande and should be published soon for full availability 
  
    # eo get travel time / gravity -----


  # GET CROWN OF THORNS STARFISH OUTBREAKS ----
  
  # This one comes from 
  
  # Dumas, P., Fiat, S., Durbano, A., Peignon, C., Mou-Tham, G., Ham, J., Gereva, S., Kaku, R., Chateau, O., Wantiez, L., 
  # De Ramon N’Yeurt, A., Adjeroud, M., 2020. Citizen Science, a promising tool for detecting and monitoring outbreaks of 
  # the crown-of-thorns starfish Acanthaster spp. Sci Rep 10, 291. https://doi.org/10.1038/s41598-019-57251-8
  
  # it states that over 100 COTS/ha is an outbreak and states the different outbreaks up to 2020 
  # We'll use the biological dataset that detects COTS and use the threshold to determine if there was an outbreak at the time of the observation or not
  
    # eo get cots oubreaks ----

  




# PROCESS DATA ----

  # TEMPERATURE ----
    
    # 5 states : well below, below, average, over, well over
    # Reading files
    # Creating simple path
    pathTemp <- file.path(pathEnv, "Temperature")
    
    # Getting file name
    sstFileName <- file.path(pathTemp,"noaacrwsstDaily_2013.nc")
    
    # ncdf4 to check
    nc <- nc_open(sstFileName)
    print(nc)
    nc$dim
    
    # Terra brick to access data
    sstR <- rast(sstFileName, subds = "analysed_sst")
    sstR[[1]]
    plot(sstR[[1]])
    plot(sstR[[90]])
    plot(sstR[[180]])
    plot(sstR[[270]])
    
    # Over the ten years data: quantiles 5. Thresholds from there, like the biological data. same method
    # We know different SSTs are highly correlated so anyone will work. let's take the mean of the coldest month why not ?
    
    
  # HEATWAVES ----
    # Bleaching alert Area of the Coral Reef Watch
    # Simplifying the stress levels into three states:
    # No alert (first 1-3 levels)
    # Bleaching (Alert 1-2)
    # Mortality (Alert 3-5)
    
    # Creating simple path
    pathHeat <- file.path(pathEnv, "Heatwaves_BAA")
    
    # Getting file name
    baaFileName <- file.path(pathHeat,"noaacrwbaa7dDaily_2013.nc")
    
    # ncdf4 to check
    nc <- nc_open(baaFileName)
    print(nc)
    nc$dim
    
    # Terra brick to access data
    baaR <- rast(baaFileName, subds = "bleaching_alert_area")
    baaR[[1]]
    plot(baaR[[350]])
    
    
    fix_orientation <- function(r) {
      y <- terra::yFromRow(r)
      if (y[1] > y[length(y)]) {
        r <- terra::flip(r, "vertical")
      }
      r
    }
    
    plot(fix_orientation(baaR[[1]]))
    
    plot(flip(baaR[[350]], direction = "vertical"))
    terra::ext(baaR)
    ext(sstR)
    terra::ylabel(baaR)
    
    nc <- nc_open(file.path(pathDat,"Environment","CRW","crw_baa_2013_2018.nc"))
    print(nc)
    names(nc$var)
    ncvar_get(nc, "bleaching_alert_area")
    
    
    
    nc <- nc_open(file.path(pathDat,"Environment","CRW","noaacrwdhwDaily_b396_7599_f8b5.nc"))
    nc <- nc_open(file.path(pathDat,"Environment","CRW","noaacrwdhwDaily_83e9_2084_6dff.nc"))
    
    
    # eo heatwaves ----
    
    ra <- rast(file.path(pathDat,"Environment","Chlorophyll_a","noaacwNPPVIIRSSQchlaDaily_8938_b7ee_2602.nc"), subds = "chlor_a")
    plot(flip(ra[[2]], direction = "vertical"))
    
    
    
    
    
  # CYCLONES ----
    # Read cyclone tracks
    # cyc <- st_read(file.path(pathEnv, "Cyclones"), layer = "IBTrACS.SP.list.v04r01.points")
    cyc <- st_read(file.path(pathEnv, "Cyclones","IBTrACS.SP.list.v04r01.points.shp"))
    
    # Filtering cyclones within boundaries and between 2013 and 2024
    cyc <- cyc[cyc$year >= 2013,]
    # st_crs(cyc)
    
    # Removing track parts that have NA R34, meaning we consider the depression became "negligible"
    # Potentially add later back points that do not have R34 but USA_WIND values > 34 ?
    quadrantsNames <- c("USA_R34_NE","USA_R34_SE","USA_R34_SW","USA_R34_NW")
    cyc <- cyc[complete.cases(st_drop_geometry(cyc[ , quadrantsNames])), ]
    
    # cyc radii per quadrant is in nautical miles so we convert them to meter for the buffer
    nmToMeter <- 1852
    cyc <- cyc |>
      mutate(
        r34_m_ne = USA_R34_NE*nmToMeter,
        r34_m_se = USA_R34_SE*nmToMeter,
        r34_m_sw = USA_R34_SW*nmToMeter,
        r34_m_nw = USA_R34_NW*nmToMeter
        )
    
    # Function to create quadrants sector polygons
    make_quadrant <- function(center, radius, start_bearing, end_bearing, n = 90) {
      # bearings sequence
      bseq <- seq(start_bearing, end_bearing, length.out = n)
      
      # great-circle points at given radius
      pts <- geosphere::destPoint(center, bseq, radius)
      
      # build polygon closing at center
      pol <- rbind(center, pts, center)
      
      st_polygon(list(as.matrix(pol))) |> st_sfc(crs = 4326)
    }
  
    # Build the four R34 quadrants
    quad_polygons <- cyc |>
      rowwise() |>
      mutate(
        center = list(c(st_coordinates(geometry)[1], st_coordinates(geometry)[2])),
        q_ne = make_quadrant(center, r34_m_ne, 0, 90),
        q_se = make_quadrant(center, r34_m_se, 90, 180),
        q_sw = make_quadrant(center, r34_m_sw, 180, 270),
        q_nw = make_quadrant(center, r34_m_nw, 270, 360)
      )
    
    # Binding quadrants
    quads <- quad_polygons |>
      # st_as_sf() |>
      st_drop_geometry() |>
      select(SID, ISO_TIME, SEASON, NAME, NATURE, q_ne, q_se, q_sw, q_nw) |>
      tidyr::pivot_longer(cols = starts_with("q_"),
                          names_to = "quadrant",
                          values_to = "geometry") |>
      st_as_sf(sf_column_name = "geometry")
    
    
    # Dissolve each point polygon into polygon per cyclone
    storm_polygons <- quads |>
      group_by(SID, SEASON, NAME) |>
      summarise(do_union = TRUE) |>
      st_make_valid()
    
    storm_polygons <- quads |>
      group_by(SID, SEASON, NAME) |>
      summarise(
        geometry = st_union(geometry),  # union of all quadrant polygons per storm
        .groups = "drop"
      ) |>
      st_make_valid()
      
    plot(st_geometry(storm_polygons[c(1:9),]))
    
    # SAVE
    dir.create(file.path(pathDat, "01_Processed","Environment","Cyclones"), showWarnings = FALSE)
    st_write(storm_polygons, dsn = file.path(pathDat, "01_Processed","Environment","Cyclones"), layer = "Cyclones_SP_R34_Influence_2013-2025.shp", driver = "ESRI Shapefile", append = FALSE) 
    
    
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
    
    
    # IBTrACS.SP.v04r01.nc
    # nc <- nc_open(file.path(pathEnv,"Cyclones","IBTrACS.SP.v04r01.nc"))
    # names(nc$var)
    # ra <- rast(file.path(pathEnv,"Cyclones","IBTrACS.SP.v04r01.nc"), subds = "name")
    
    
    