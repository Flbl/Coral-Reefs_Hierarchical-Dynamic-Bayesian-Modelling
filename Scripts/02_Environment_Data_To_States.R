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
    
  # GET ONI INDEX (EL NINO EL NINA PHASES) ----
  
  # We use the noaa historical metrics that we manually enter into a R dataframe to relate yearly main phase from the dataset on their main website
  # Check https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/ONI_v5.php
  # Check https://ggweather.com/enso/oni.htm

    # eo get oni index (el nino el nina phases) ----



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
    start_year = 2002,
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

  # READ BIOLOGICAL DATA FOR STATION ASSIGNATION ----
  # Coral
  coralRorc <- read.csv(file.path("Data", "Species","RORC_Coral_Data_hdbn.csv"))

  # Fish
  fishRorc <- read.csv(file.path("Data", "Species","RORC_Fish_Data_hdbn.csv"))

  # Inv
  invRorc <- read.csv(file.path("Data", "Species","RORC_Inv_Data_hdbn.csv"))

  # Create common station df
  stationsRorc <- st_as_sf(
    merge(
      merge(
        coralRorc[!is.na(coralRorc$Lon),] %>% distinct(Site, Station, Lon, Lat),
        fishRorc[!is.na(fishRorc$Lon),] %>% distinct(Site, Station, Lon, Lat),
        all = TRUE),
      invRorc[!is.na(invRorc$Lon),] %>% distinct(Site, Station, Lon, Lat),
      all = TRUE),
    coords = c("Lon", "Lat"),
    crs = 3163
    )
  
  # Convert to wgs84
  stationsRorc <- st_transform(stationsRorc, 4326)
    
    # eo read biological data for station assignation ----

  
  # ONI INDICATOR FOR EL NINO EL NINA GENERAL ----
  
  # Check https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/ONI_v5.php
  # There are three possible phases: Nino, Neutral, Nina
  # Make a table with rules:
  # If a year only shows one major phase, then its assigned to the phase
  # If there are 2 or more phases in the same year, check values and assign the most important
  # Check https://ggweather.com/enso/oni.htm
  # for appreciation
  
  ninoa <- data.frame(Year = c(2013:2024),
                      State = c(
                        "Neutral",
                        "Neutral",
                        "El Nino",
                        "El Nino",
                        "El Nina",
                        "El Nino",
                        "El Nino",
                        "El Nina",
                        "El Nina",
                        "El Nina",
                        "El Nino",
                        "El Nino"
                      ))
  
  ninoa
  
  dir.create(file.path(pathDat, "01_Processed","Environment","OceaniNinoIndex"), showWarnings = FALSE)
  write.csv(ninoa, file.path(pathDat, "01_Processed","Environment","OceaniNinoIndex","Env_General_ONI_New_Caledonia.csv"), row.names = FALSE)
  

  # TEMPERATURE ----
    
    # Read data, check time series, extract values for each station, get the quantiles as threshold
  
    # Reading files
    # Creating simple path
    pathTemp <- file.path(pathEnv, "Temperature")
  
    # Over the ten years data: quantiles 5. Thresholds from there, like the biological data. same method
    # We know different SSTs (min, mean, max) are highly correlated so anyone will work. annual mean SST will do to provide the background regime of the site/stations
    # 5 states : temperature regime for the year compared to its last 10 years: very cool, slightly cool, usual, slightly warm, very warm
    # "Compared to the last decade at this site, was this year cooler/normal/warmer ?
    # Anomaly magnitude discretized into 5 states:
    # -Very cool regime
    # -Slightly cool regime
    # -Near usual
    # -Slightly warm regime
    # -Very warm regime
    # The anomaly is evaluated through the site quantiles and an ecologically relevant hybrid threshold of 0.3 
    
    # Trying map from purrr package (apply a function to each element of a vector)
    # Reading raw data, compute annual mean on all raster cells and store the annual mean rasters per year inside a list
    yea <- 2002:2024
    annual_mean_sst <- function(year) {
      
      message("Reading file for year ",year)
      r <- rast(file.path(pathTemp, paste0("noaacrwsstDaily_", year, ".nc")), subds = "analysed_sst")
      
      message("calculating annual mean")
      res <- mean(r)
      message("Done")
      
      return(res)
    }
    
    # Applying function
    annualSSTMeanBrick <- map(yea, annual_mean_sst)
    names(annualSSTMeanBrick) <- yea
    # plot(annualSSTMeanBrick[[1]])
    # plot(annualSSTMeanBrick[[4]])
    
    # Now we have a list of single rasters with averaged cell values over their year.
    
    # Extracting station annual SST means
    # Function
    extract_year <- function(year, annual_mean_list, stations_sf) {
      
      r <- annual_mean_list[[as.character(year)]]

      vals <- terra::extract(r, vect(stations_sf), search_radius = 12000) #2 cell large~ 
      
      tibble(
        site = stations_sf$Site,
        station = stations_sf$Station,
        year = year,
        annual_mean = vals[,2]
      )
    }
    
    # Create the Station level annual mean SST time series
    sst_station_year <- map_df(yea, extract_year, 
                                annual_mean_list = annualSSTMeanBrick,
                                stations_sf = stationsRorc)
    
    
    # plot(annualSSTMeanBrick[[1]])
    # replot(annualSSTMeanBrick[[1]])
    # plot(st_geometry(stationsRorc), add = TRUE, pch = 3)
    # text(st_coordinates(stationsRorc),
    #      labels = stationsRorc$Station,
    #      pos = 3, cex = 0.6)
    
    
    # aggregate the station time series to site time series
    sst_site_year <- sst_station_year %>%
      group_by(site, year) %>%
      summarise(
        annual_mean_SST = mean(annual_mean, na.rm = TRUE),
        .groups = "drop"
      )
    
    
    # Now let's compute the 10 yearss prior baseline for each year for each site
    compute_baseline <- function(df) {
      df %>%
        arrange(site, year) %>%
        group_by(site) %>%
        mutate(
          baseline_10yr =
            sapply(seq_along(year), function(i) {
              # years strictly before current year
              past_vals <- annual_mean_SST[max(1, i-10):(i-1)]
              
              if(length(past_vals) < 10) return(NA)  # avoid fragile small windows
              
              mean(past_vals, na.rm = TRUE)
            })
        ) %>%
        ungroup()
    }
    
    # Computing baselines
    sst_site_year <- compute_baseline(sst_site_year)
    
    # Now computing the anomaly of temperature between year observed and 10 year baseline:
    sst_site_year <- sst_site_year %>%
      mutate(
        anomaly = annual_mean_SST - baseline_10yr
      )
    
    # Setting the hybrid delta threshold over which we can trigger a state change (0.3 degrees)
    # Ref of 0.3 degrees on mean annual SST for coral reefs :
    # Coral reef Watch documentation on SST anomalies stating +-0.2°C is accuracy bias 
    # https://coralreefwatch.noaa.gov/product/5km/methodology.php#sst
    # and cite 1-2 paper on 0.3
    # 
    delta <- 0.3
    
    # Getting to quantiles:
    # What's the quantile sequence to get evenly spaced values ?
    # Deciles of appropriate choice for wide central usual region and narrow truly extreme tails
    # That gives 10,30,70,90 instead of 20,40,60,80 for equal quantiles 
    seq(0,1,1/(5))
    
    quantiles_site <- sst_site_year %>%
      group_by(site) %>%
      summarise(
        q10 = quantile(anomaly, 0.1, na.rm = TRUE),
        q30 = quantile(anomaly, 0.3, na.rm = TRUE),
        q70 = quantile(anomaly, 0.7, na.rm = TRUE),
        q90 = quantile(anomaly, 0.9, na.rm = TRUE),
        .groups = "drop"
      )
    
    # And finally getting the states for the temperature data for each site:
    # join quantiles and classify with hybrid rule
    sst_site_year <- sst_site_year %>%
      left_join(quantiles_site, by = "site") %>%
      mutate(
        SST_regime_state = case_when(
          
          # --- VERY COOL REGIME ---
          anomaly < q10 & anomaly <= -delta ~ "Very cool regime",
          
          # --- SLIGHTLY COOL REGIME ---
          anomaly >= q10 & anomaly < q30 ~ "Slightly cool regime",
          anomaly < q10 & abs(anomaly) < delta ~ "Slightly cool regime",
          
          # --- USUAL REGIME ---
          anomaly >= q30 & anomaly <= q70 & abs(anomaly) < delta ~ "Usual regime",
          
          # --- SLIGHTLY WARM REGIME ---
          anomaly > q70 & anomaly <= q90 ~ "Slightly warm regime",
          anomaly > q90 & anomaly < delta ~ "Slightly warm regime",
          
          # --- VERY WARM REGIME ---
          anomaly > q90 & anomaly >= delta ~ "Very warm regime",
          
          TRUE ~ NA_character_
        )
      )
    
    sst_site_year_states <- sst_site_year[sst_site_year$year >2012,]
    
    dir.create(file.path(pathDat, "01_Processed","Environment","Temperature"), showWarnings = FALSE)

    write.csv(sst_site_year_states, file = file.path(pathDat, "01_Processed","Environment","Temperature","RORC_Env_TemperatureRegime_Site_States_hdbn.csv"), row.names = FALSE)
    
    # eo temperature ----
    
    
  # CHLOROPHYL-A ----
    
    # Cholorophyl-a levels here are similar to the sst regime
    # We don't want "high, good or bad" but whether the background regime was unusually low or high compared to the baseline
    # We need a delta to work with the time series to express bias and usual ecological change
    # The Copernicus manual on the dataset states that the ubRMSD (only viable indicator for log transformed chla data) is the right oneand equals 0.34 mg/m³
    # So we'll use the same logic with a delta = 0.34 and then make 5 states with quantiles 10-30-70-90 with rules <q10 & > delta = very unproductive year
    # 0.34 apparently is quite high for coral reef oligotrophy (0.01-0.15~~ to validate) so a change outside the delta could already be considered ("slight" states)
    # Read data, check time series, extract values for each station, get the quantiles as threshold
    
    # Reading files
    # Creating simple path
    pathChla <- file.path(pathEnv, "Chlorophyll_a")
    
    # Read data
    # 2002-2012 period
    r <- rast(file.path(pathChla, paste0("cmems_obs-oc_glo_bgc-plankton_my_l4-gapfree-multi-4km_P1D_multi-vars_155.02E-174.98E_26.98S-14.02S_2002-01-01-2012-12-31.nc")), subds = "CHL")
    # 2013-2024 period
    r2 <- rast(file.path(pathChla, paste0("cmems_obs-oc_glo_bgc-plankton_my_l4-gapfree-multi-4km_P1D_multi-vars_155.02E-174.98E_26.98S-14.02S_2013-01-01-2024-12-31.nc")), subds = "CHL")
    
    # Plot
    # plot(r2[[160]])
    
    # get dates in year-month
    rDates  <- format(time(r), "%Y-%m")
    r2Dates <- format(time(r2), "%Y-%m")
    
    # Compute monthly mean
    rYearsMonthlyMean <- tapp(r, rDates, fun = mean, na.rm = TRUE) 
    # plot(rYearsMonthlyMean[[12]])
    # values(rYearsMonthlyMean)
    r2YearsMonthlyMean <- tapp(r2, r2Dates, fun = mean, na.rm = TRUE) 
    
    # Log-transform monthly CHL because copernicus chl product is log-narmally distributed
    #  + 1e-6 to avoid zero just in case
    rYearsMonthlyMean  <- log10(rYearsMonthlyMean  + 1e-6)
    r2YearsMonthlyMean <- log10(r2YearsMonthlyMean + 1e-6)
    
    # # Get year and month vectors
    # rYear <- year(rDates)
    # rMonth <- month(rDates)
    # r2Year <- year(r2Dates)
    # r2Month <- month(r2Dates)
    # 
    # # split timeseries into yearly list
    # # Baseline (2002-2012)
    # rYearsMonth <- split(r, list(rYear, rMonth))
    # names(rYears) <- as.character(unique(year(rDates)))
    # # Model timeline (2013-2024)
    # r2Years <- split(r2, year(r2Dates))
    # names(r2Years) <- as.character(unique(year(r2Dates)))
    # 
    # # Convert daily data per year into yearly average
    # rYears <- lapply(rYears, mean, na.rm = TRUE) #small data error with only NaN on the 2002/12/31 so we remove it
    # r2Years <- lapply(r2Years, mean)
    
    # Remove the "X" from names to get only Year.Month
    names(rYearsMonthlyMean) <- gsub("X","",names(rYearsMonthlyMean))
    names(r2YearsMonthlyMean) <- gsub("X","",names(r2YearsMonthlyMean))
    
    get_r_year <- function(x) substr(names(x),1,4) #character 1 to 4 = year (e.g. 2002.01)
    # get_r_year(rYearsMonthlyMean[[1]])
    
    # split monthly rasters by year and assign list names
    rYearsMonthlyMeanList  <- split(rYearsMonthlyMean,  get_r_year(rYearsMonthlyMean))
    names(rYearsMonthlyMeanList) <- unique(get_r_year(rYearsMonthlyMean))
    r2YearsMonthlyMeanList <- split(r2YearsMonthlyMean, get_r_year(r2YearsMonthlyMean))
    names(r2YearsMonthlyMeanList) <- unique(get_r_year(r2YearsMonthlyMean))
    
    
    # Function to extract top three year
    extract_top3_year <- function(year, monthly_list, stations_sf, radius = 12000, nbMonth = 3) {
      # Test zone
      # year = "2013"
      # monthly_list = r2YearsMonthlyMeanList
      # stations_sf = stationsRorc
      # radius = 12000
      # nbMonth = 3
      # eo tz
      
      r_year <- monthly_list[[as.character(year)]]
      
      # extract each monthly layer separately (required for search_radius)
      vals_list <- lapply(1:nlyr(r_year), function(i) {
        terra::extract(
          r_year[[i]],
          vect(stations_sf),
          search_radius = radius,
          na.rm = TRUE
        )[,2]   # drop ID column
      })
      
      # bind months column-wise → stations × 12 matrix
      vals_mat <- do.call(cbind, vals_list)
      
      # compute top-3 monthly mean per station
      # apply iterates over a matrix
      # 1 = row-wise (2 = column wise)
      top3_mean <- apply(vals_mat, 1, function(x) {
        # x = one station's 12 monthly values
        # We remove potential missing months (although here we know there aren't any)
        x <- x[!is.na(x)]
        # We enforce a small data requirement (statistics not meaningful if data exist for less than three months over the year) (yeah it's heavy shielded)
        if (length(x) < 3) return(NA_real_) #NA_real_ ensures numeric NA
        # We sort the 12 months by descending magnitude
        # We lost the names/if its the same months over the different stations/years but it doesn't really matter as we state:
        # "We want the mean of the top 3 most productive months", it can be any non consecutive month out of the 12  
        mean(sort(x, decreasing = TRUE)[1:nbMonth])
      })
      
      tibble::tibble(
        site     = stations_sf$Site,
        station = stations_sf$Station,
        year     = as.integer(year),
        chl_top3_log = top3_mean
      )
    }
    
    # Baseline period
    chl_station_year <- purrr::map_df(
      names(rYearsMonthlyMeanList),
      extract_top3_year,
      monthly_list = rYearsMonthlyMeanList,
      stations_sf  = stationsRorc,
      nbMonth = 3
    )
    
    # Model period
    chl_station_year2 <- purrr::map_df(
      names(r2YearsMonthlyMeanList),
      extract_top3_year,
      monthly_list = r2YearsMonthlyMeanList,
      stations_sf  = stationsRorc,
      nbMonth = 3
    )
    
    chl_station_year <- bind_rows(chl_station_year, chl_station_year2)
    
    # aggregate the station time series to site time series
    chl_site_year <- chl_station_year %>%
      group_by(site, year) %>%
      summarise(
        chl_top3_log = mean(chl_top3_log, na.rm = TRUE),
        .groups = "drop"
      )
    
    
    # Now let's compute the 10 yearss prior baseline for each year for each site
    compute_baseline <- function(df) {
      df %>%
        arrange(site, year) %>%
        group_by(site) %>%
        mutate(
          baseline_10yr =
            sapply(seq_along(year), function(i) {
              # years strictly before current year
              past_vals <- chl_top3_log[max(1, i-10):(i-1)]
              
              if(length(past_vals) < 10) return(NA)  # avoid fragile small windows
              
              mean(past_vals, na.rm = TRUE)
            })
        ) %>%
        ungroup()
    }
    
    # Computing baselines
    chl_site_year <- compute_baseline(chl_site_year)
    
    # Now computing the anomaly of CHL between year observed and 10 year baseline:
    chl_site_year <- chl_site_year %>%
      mutate(
        anomaly = chl_top3_log - baseline_10yr
      )
    
    # Setting the hybrid delta threshold over which we can trigger a state change
    # Threshold here is 0.31 as ubRMSD of the GLO chl daily data method (accuracy indicator)
    # /!\ ubRMSD is actually the validation between in situ and satellite. not the accuracy of the model ?
    # Need to check 
    delta <- 0.34
    
    # Getting to quantiles:
    # What's the quantile sequence to get evenly spaced values ?
    # Deciles of appropriate choice for wide central usual region and narrow truly extreme tails
    # That gives 10,30,70,90 instead of 20,40,60,80 for equal quantiles 
    seq(0,1,1/(5))
    
    quantiles_site <- chl_site_year %>%
      group_by(site) %>%
      summarise(
        q10 = quantile(anomaly, 0.1, na.rm = TRUE),
        q30 = quantile(anomaly, 0.3, na.rm = TRUE),
        q70 = quantile(anomaly, 0.7, na.rm = TRUE),
        q90 = quantile(anomaly, 0.9, na.rm = TRUE),
        .groups = "drop"
      )
    
    # And finally getting the states for the temperature data for each site:
    # join quantiles and classify with hybrid rule
    chl_site_year <- chl_site_year %>%
      left_join(quantiles_site, by = "site") %>%
      mutate(
        CHL_regime_state = case_when(
          
          # --- VERY low REGIME ---
          anomaly < q10 & anomaly <= -delta ~ "Very low regime",
          
          # --- SLIGHTLY low REGIME ---
          anomaly >= q10 & anomaly < q30 ~ "Slightly low regime",
          anomaly < q10 & abs(anomaly) < delta ~ "Slightly low regime",
          
          # --- USUAL REGIME ---
          anomaly >= q30 & anomaly <= q70 & abs(anomaly) < delta ~ "Usual regime",
          
          # --- SLIGHTLY high REGIME ---
          anomaly > q70 & anomaly <= q90 ~ "Slightly high regime",
          anomaly > q90 & anomaly < delta ~ "Slightly high regime",
          
          # --- VERY high REGIME ---
          anomaly > q90 & anomaly >= delta ~ "Very high regime",
          
          TRUE ~ NA_character_
        )
      )
    
    chl_site_year_states <- chl_site_year[chl_site_year$year >2012,]
    
    write.csv(chl_site_year_states, file = file.path(pathPro,"RORC_Env_ChlorophyllaRegime_Site_States_hdbn.csv"), row.names = FALSE)
    
    # eo chlorophyll a ----
    
    

    
  # "HEATWAVES" ----
    # Bleaching alert Area of the Coral Reef Watch
    # Simplifying the stress levels into three states:
    # No alert (first 1-3 levels)
    # Bleaching (Alert 1-2)
    # Mortality (Alert 3-5)
    
    # Creating simple path
    pathHeat <- file.path(pathEnv, "Heatwaves_BAA")
    
    # Getting file name
    baaFileName <- file.path(pathHeat,"noaacrwbaa7dDaily_2019.nc")

    # Terra brick to access data
    baaR <- rast(baaFileName, subds = "bleaching_alert_area")
    # baaR[[1]]
    plot(baaR[[1:12]])
    plot(baaR[[181:192]])
    plot(baaR[[90]])
    
    # On linux there are problems of raster orientation but on windows there isn't
    # fix_orientation <- function(r) {
    #   y <- terra::yFromRow(r)
    #   if (y[1] > y[length(y)]) {
    #     r <- terra::flip(r, "vertical")
    #   }
    #   r
    # }
    # plot(baaR[[1]])
    # plot(fix_orientation(baaR[[1]]))

    # for each year, extract the max value
    # Create the year vector
    yea <- 2013:2024

    # Storing yearly data into a list 
    annualDailyBaaBrick <- lapply(yea, function(x) rast(file.path(pathHeat, paste0("noaacrwbaa7dDaily_", x, ".nc")), subds = "bleaching_alert_area"))
    names(annualDailyBaaBrick) <- yea
    
    # Extracting station max baa value for each year
    # Function
    extract_baa_year <- function(year, annualdaily_list, stations_sf, radius = 12000) {
      # Test zone
      # year = "2013"
      # annualdaily_list = annualDailyBaaBrick
      # stations_sf = stationsRorc
      # radius = 12000
      # eo tz
      
      
      r <- annualdaily_list[[as.character(year)]]
      
      # extract each monthly layer separately (required for search_radius)
      vals_list <- lapply(1:nlyr(r), function(i) {
        terra::extract(
          r[[i]],
          vect(stations_sf),
          # fun = max,
          search_radius = radius,
          na.rm = TRUE
        )[,2]   # drop ID column
      })
      
      # bind months column-wise → stations × 12 matrix
      vals_mat <- do.call(cbind, vals_list)
      
      # Extract max
      maxBaa <- apply(vals_mat, 1, function(x) {
        # # x = one station 365-366 daily value
        # # We remove potential missing data
        # x <- x[!is.na(x)]
        # We extract the max
        max(x, na.rm = TRUE)
      })
      

      tibble(
        site = stations_sf$Site,
        station = stations_sf$Station,
        year = year,
        annual_max_BAA = maxBaa
      )
    }
    
    # Create the Station level annual mean SST time series
    baa_station_year <- map_df(yea, extract_baa_year, 
                               annualdaily_list = annualDailyBaaBrick,
                               stations_sf = stationsRorc)
    
    
    
    # Levels for states:
    # No alert (first 1-3 levels) = 0-2
    # Bleaching (Alert 1-2) = 3-4
    # Mortality (Alert 3-5) = 5-7
    
    baa_station_year <- baa_station_year%>%
      mutate(
        State = case_when(
          annual_max_BAA %in% c(0:2) ~ "No alert",
          annual_max_BAA %in% c(3:4) ~ "Bleaching alert",
          annual_max_BAA %in% c(5:7) ~ "Mortality alert",
          TRUE ~ NA_character_
        )
      )
    
    dir.create(file.path(pathDat, "01_Processed","Environment","Heawaves_BAA"), showWarnings = FALSE)
    write.csv(baa_station_year, file.path(pathDat, "01_Processed","Environment","Heawaves_BAA","Env_HeatwavesBAA_Station_States_New_Caledonia.csv"), row.names = FALSE)
    
    
    # eobaaR# eo heatwaves ----
    

    
  # CYCLONES ----
    # Use presence/absence data at the general scale because waves
    # AND station node with R34 because waves+wind = damage
    
    # ---
    # Read raw open data base and extract New Caledonia cyclone tracks (DONE ONCE), then we read directly the new caledonia cyclone file
    # cyc <- st_read(file.path(pathEnv, "Cyclones","IBTrACS.SP.list.v04r01.points.shp"))
    # 
    # # Creating bounding box to filter only New Caledonia
    # ncBbox <- st_bbox(
    #   c(
    #     xmin = 155, ymin = -27, xmax = 175, ymax = -14
    #     ),
    #   crs = st_crs(4326)
    #   )
    # 
    # # Convert to sf
    # ncBboxSf <- st_as_sfc(ncBbox)
    # 
    # # Intersect with cyclones
    # cycNc <- cyc[st_intersects(cyc, ncBboxSf, sparse = FALSE), ]
    # # plot(st_geometry(cycNc))
    # 
    # st_write(cycNc, dsn = file.path(pathDat, "01_Processed","Environment","Cyclones"), layer = "Cyclones_New_Caledonia_IBTrACS.shp", driver = "ESRI Shapefile", append = FALSE) 
    # ---
    
    cyc <- st_read(file.path(pathDat, "01_Processed","Environment","Cyclones","Cyclones_New_Caledonia_IBTrACS.shp"))
    
    # Filtering cyclones within boundaries and between 2013 and 2024
    cyc <- cyc[cyc$year >= 2013,]
    # st_crs(cyc)
    
    
    
    # GET GENERAL COUNT OF STORMS PER YEAR
    
    # Count the number of cyclones  per year and 
    cycCount <- cyc %>%
      distinct(SEASON, NAME) %>%
      mutate(CYCLONE = 1) %>%
      group_by(SEASON) %>%
      summarise(CYCLONE = sum(CYCLONE))
      
    # create states out of it for the general nodes (named "Storms" since it also registers depressions etc)
    cycCount <- cycCount%>%
      mutate(
      State = case_when(
        CYCLONE %in% c(0) ~ "No Storm",
        CYCLONE %in% c(1:4) ~ "Moderate Storm Season",
        CYCLONE %in% c(5:7) ~ "Intense Storm Season",
        TRUE ~ NA_character_
        )
      )
        
    # Save !
    write.csv(cycCount, file.path(pathDat, "01_Processed","Environment","Cyclones","Env_StormCount_General_States_New_Caledonia.csv"), row.names = FALSE)
    
    # Removing track parts that have NA R34, meaning we consider the depression became "negligible"
    # Potentially add later back points that do not have R34 but USA_WIND values > 34 ?
    quadrantsNames <- c("USA_R34_NE","USA_R34_SE","USA_R34_SW","USA_R34_NW")
    cyc <- cyc[complete.cases(st_drop_geometry(cyc[ , quadrantsNames])), ]
    
    
    
    # COMPUTE STATION WISE R34 BUFFER 
    
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
      
    # plot(st_geometry(storm_polygons[c(1:2),]), add = TRUE)
    
    # SAVE
    dir.create(file.path(pathDat, "01_Processed","Environment","Cyclones"), showWarnings = FALSE)
    st_write(storm_polygons, dsn = file.path(pathDat, "01_Processed","Environment","Cyclones"), layer = "Cyclones_SP_R34_Influence_2013-2025.shp", driver = "ESRI Shapefile", append = FALSE) 
    
    # Create states for station wise direct cyclone passage
    # Small crs check just in case
    st_crs(storm_polygons)
    st_crs(stationsRorc)
    
    # Join bot object through intersection
    stations_storms <- st_join(
      stationsRorc,
      storm_polygons,
      join = st_intersects,
      left = TRUE
    )
    
    # Count the number of cyclones per station per year:
    cyclone_R34_counts <- stations_storms %>%
      st_drop_geometry() %>% # drop geometry
      filter(!is.na(SEASON)) %>%        # remove non-intersections
      group_by(Site, Station, SEASON) %>%
      summarise(
        n_cyclones = n(),   # or n(), see note below
        .groups = "drop"
      )
    
    # Since it was an intersection, some years or stations may have been put aside so we build the grid back
    # all_years <- sort(unique(storm_polygons$SEASON))
    
    cyclone_R34_counts <- expand.grid(
      # Site    = unique(stationsRorc$Site),
      Station = unique(stationsRorc$Station),
      SEASON  = all_years
    ) %>%
      left_join(st_drop_geometry(stationsRorc),
                by = c("Station")) %>%
      left_join(cyclone_R34_counts,
                by = c("Site","Station", "SEASON")) %>%
      mutate(n_cyclones = ifelse(is.na(n_cyclones), 0, n_cyclones))
    
    
    cyclone_R34_States <- cyclone_R34_counts %>%
      mutate(
        States = case_when(
          n_cyclones == 0 ~ "No Storm",
          n_cyclones == 1 ~ "1 Heavy Storm",
          n_cyclones >= 2 ~ "Multiple Heavy Storms",
          TRUE ~ NA_character_
        )
      )
    
    # WRITE THE STATES !
    write.csv(cyclone_R34_States, file.path(pathDat, "01_Processed","Environment","Cyclones","Env_R34StormCount_Station_States_New_Caledonia.csv"), row.names = FALSE)
    
    
    
  # GRAVITY ----
    # We'll read the gravity shapefile and distinguish eventually the values using quantiles for 2 to three states
    # (severly, mild, none)
    pathGrav <- file.path(pathEnv, "Gravity")
    grav <- rast(file.path(pathEnv, "Gravity","Human_Gravity_NewCaledonia_2019.tif"))
    plot(grav)
    # Transform to epsg84
    stationsRorc_lambert <- st_transform(stationsRorc, st_crs(grav))
    
    # extract 
    gravStates <- stationsRorc %>%
      st_drop_geometry() %>%
      mutate(
        Gravity = terra::extract(grav, stationsRorc_lambert, search_radius = 500)[,2]
      )
    
    # Create states
    # Taking station's dataset quantiles because the scales are of the station's chart so they'd all be "High" if using the New Caledonia quantiles
    # seq(0,1,1/(3))
    # qt <- quantile(values(grav, na.rm = TRUE), seq(0,1,1/(3)))
    qt <- quantile(gravStates$Gravity, c(1/3,2/3))
    
    gravStates <- gravStates %>%
      mutate(
        States = case_when(
          Gravity < qt[1] ~ "Low",
          Gravity >= qt[1] & Gravity < qt[2] ~ "Moderate",
          Gravity > qt[2] ~ "High"
        )
      )
    
    dir.create(file.path(pathDat, "01_Processed","Environment","Gravity"), showWarnings = FALSE)
    write.csv(gravStates, file.path(pathDat, "01_Processed","Environment","Gravity","Env_Gravity_Station_States_New_Caledonia.csv"), row.names = FALSE)
    
    
    
  
    
    # COTS OUTBREAKS ----
    # ACA > 100/ha = outbreak
    # See Dumas et al. 2020
    # So by conversion to 20*5 = 100sqm, then cots >= 1 ==> outbreak) Not good.
    # So we'll be a bit smoother, we'll count in this case over the 3-4 transects as one because observing only one COTS on 20*5 meters cannot be an indicator of 100COTS/ha
    # Station wise invasion
    # We transform transects into a single transect of x meters (20*nbTransects, usually 4)*5 meters wide
    # Then we compare the surface observed to 100 cots/ha
    
    cots <- read.csv(file.path(pathEnv, "COTS","RORC_Invertebrate.csv"), header = TRUE)
    cots <- cots %>%
      mutate(TransectSurface = 20*5) %>%
      select(Campagne, Site, Station, Transect, TransectSurface, ACA) %>%
      group_by(Campagne, Site, Station) %>%
      summarise(Cots = sum(ACA),
                TransectSurface = sum(TransectSurface),
                .groups = "drop")
    
    # Threshold conversion of 100cots/10 000 sqm
    cots$Threshold <- (cots$TransectSurface*100)/10000
    
    # State
    cots <- cots %>% 
      mutate(
        States = case_when(
          Cots < Threshold ~ "No Outbreak",
          Cots >= Threshold ~ "Outbreak"
        
      ))
    
    # WRITE STATES !
    dir.create(file.path(pathDat, "01_Processed","Environment","COTS"), showWarnings = FALSE)
    write.csv(cots, file.path(pathDat, "01_Processed","Environment","COTS","Env_COTS_Outbreaks_Station_States_New_Caledonia.csv"), row.names = FALSE)
    
    


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
    
    # # Getting file name
    # sstFileName <- file.path(pathTemp,"noaacrwsstDaily_2013.nc")
    # 
    # # ncdf4 to check
    # nc <- nc_open(sstFileName)
    # print(nc)
    # nc$dim
    # 
    # # Terra brick to access data
    # sstR <- rast(sstFileName, subds = "analysed_sst")
    # sstR[[1]]
    # plot(sstR[[1]])
    # plot(sstR[[2]])
    # plot(sstR[[3]])
    # plot(sstR[[4]])
    # plot(sstR[[90]])
    # plot(sstR[[180]])
    # plot(sstR[[270]])
    # 
    