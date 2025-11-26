#############################################################################################
#                                                                                           #
#                                                                                           #
#############                Convert observed data to states                   #########
#                                                                                           #
#                                                                                           #
#############################################################################################


# What does this script do:
  # It reads data that respects a particular format:
  # it has metadata columns with names including one or more of these:
  # c("Year","Country","Region","Sector","Site","Station","Sample")

  # It has biological data with columns starting with the dataset name
  # ex: RorcCoralRichness
  # This helps identifying meta from biological data

  # Later : maybe separate meta from biological obs and keep only one id column to cross with meta
  # (Slightly dangerous due to needing to always update both datasets accordingly)
  # (Probably better like this to put a marker in the column name to identify meta from variables from biological data)
  # (Different datasets can be generated if the data is very heterogeneous (different sites/samplings) to not create many NAs in eventual grid expansions) 

  # This script will read the data and convert each biological variable observation per year into a state (e.g. High, Normal, Low, Zero), depending on the available data 
  # For The General thresholds : Two methods are proposed : GMM to identify 
  # that it takes the whole temporal series of the station and checks for trends and mean value, checks the last value and coordinates a choice based on one or each results depending on data availability
  # 


# INIT -----

# Clean slate
rm(list=ls(all=TRUE))

# Libraries initialization
init1 <- tail(unlist(strsplit(rstudioapi::getActiveDocumentContext()$path, "/")), n = 1)
source("Scripts/00_Initialisation.R")


  # eo INIT ----



# READ DATA ----

# RORC ----

  # Coral
  coralRorc <- read.csv(file.path("Data", "Species","RORC_Coral_Data_hdbn.csv"))
  coralRorc$Station <- gsub(" ","_", coralRorc$Station)
  coralRorc$Site <- gsub(" ","_", coralRorc$Site)
  
  # Fish
  fishRorc <- read.csv(file.path("Data", "Species","RORC_Fish_Data_hdbn.csv"))
  
  # Inv
  invRorc <- read.csv(file.path("Data", "Species","RORC_Inv_Data_hdbn.csv"))


  # eo rorc ----
  
  # eo read data ----
  
# PREPARE DATA ----
  
  # Identify variables
  varCoral <- colnames(coralRorc)[grep("Rorc", colnames(coralRorc))]  
  varFish <- colnames(fishRorc)[grep("Rorc", colnames(fishRorc))]
  varInv <- colnames(invRorc)[grep("Rorc", colnames(invRorc))]
  
  # Aggregate datasets to their means per station
  # Aggregate by station
  # mean Coral
  coralRorcMean <- coralRorc %>%
    # select(-Sample) %>%           # drop the transect column
    group_by(Year, Country, Region, Sector, Site, Station) %>%   # your metadata columns
    summarise(
      n_samples = n_distinct(Sample),
      across(all_of(varCoral), mean),
      .groups = "drop"
      )
  
  # Sd coral (for station wise thresholds later)
  coralRorcSD <- coralRorc %>%
    select(-Sample) %>%           # drop the transect column
    group_by(Year, Country, Region, Sector, Site, Station) %>%   # your metadata columns
    summarise(across(all_of(varCoral), sd), .groups = "drop")
  
  # mean Fish
  fishRorcMean <- fishRorc %>%
    # select(-Sample) %>%           # drop the transect column
    group_by(Year, Country, Region, Sector, Site, Station) %>%   # your metadata columns
    summarise(
      n_samples = n_distinct(Sample),
      across(all_of(varFish), mean),
      .groups = "drop"
    )
  # Sd Fish 
  fishRorcSD <- fishRorc %>%
    select(-Sample) %>%           # drop the transect column
    group_by(Year, Country, Region, Sector, Site, Station) %>%   # your metadata columns
    summarise(across(all_of(varFish), sd), .groups = "drop")
  
  # Mean Inv
  invRorcMean <- invRorc %>%
    # select(-Sample) %>%           # drop the transect column
    group_by(Year, Country, Region, Sector, Site, Station) %>%   # your metadata columns
    summarise(
      n_samples = n_distinct(Sample),
      across(all_of(varInv), mean),
      .groups = "drop"
    )
  # Sd Inv 
  invRorcSD <- invRorc %>%
    select(-Sample) %>%           # drop the transect column
    group_by(Year, Country, Region, Sector, Site, Station) %>%   # your metadata columns
    summarise(across(all_of(varInv), sd), .groups = "drop")
  
  
# GENERAL THRESHOLD ----
  
  # All data is used to calulate data driven thresholds using GMM and quantiles to establish global states of quality
  # This assumes that there is enough data in the temporal series to distinguish the "health" gradient from zero/bad to good/very good
  
  # Create Thresholds from the available data
  # For each variable, make a gaussian Mixture Model and estimate the thresholds. Calculate also the quantiles 
  # Export thresholds to a CSV
  # 1st column is variable name
  # Next columns are threshold1, threshold2, threshold3...
  # Make it threshold number variable
  # (Export two files : "_GMM" for mixture models and "_QTLS" for quantiles)
  # 
  
  # Create Thresholds and save them
  
    # Function
    makeDDThresholds <- function(x, nbState, zeroState = TRUE, plot = TRUE, varName = "Variable", savePlot = FALSE, savePath){
      
      # Test Zone
      # x = coralRorc$RorcHCB
      # x = coralRorc$RorcHCO
      # x = coralRorc$RorcSD_SI
      # x = coralRorc$RorcRC
      # x = coralRorc$RorcCoralRichness
      # nbState = 5
      # zeroState = TRUE
      # plot = TRUE
      # varName = "VariableHCO"
      # eo test zone
      # print(x)
      
      # Process data
        # Convert to numeric if it isn't already
        x <- as.numeric(x)
        
        dataPos <- x[x > 0]
        
        dataPos_tr <- log(dataPos)

      # GETTING THRESHOLDS
        
        # GMM
          # Run model on positive log transformed data
          # safe_Mclust is a custom trycatch function from Mclust that yields NULL if there is an error (= model cannot run on the data)
          model <- safe_Mclust(dataPos_tr) #, G = nbState-1  , G = c(1:6)
          
          # if null model, then create empty res df
          if(is.null(model)){
            
            res <- t(as.data.frame(c("GMM",rep(NA, nbState))))
            colnames(res) <- c("Method", paste0("State", seq_len(nbState)))
            rownames(res) <- varName
            
          }else{
            
            # print(model)
            # Model checking
            # Number of modes found
            # model$G
            # Summary
            print(summary(model))
            # Uncertainty plot
            # plot(model, what = "uncertainty")
            # Uncertainty mean (closer to zero = good)
            # mean(model$uncertainty)
            # Model modes
            # model$parameters$mean
            # Transforming back to "raw values"
            print(exp(model$parameters$mean))
            # BIC
            # model$BIC
            
            # Extracting values
            params <- model$parameters
            G <- model$G                     # number of clusters
            varType <- params$variance$modelName #Variance type (equal or unequal)
            pro <- params$pro                 # mixing proportions
            mean_tr <- params$mean               # means
            sd_tr   <- sqrt(params$variance$sigmasq)  # standard deviations (univariate only)
            params
            
            # GMM Thresholds 
            
            # We've obtained breaking values from GMM, potentially expressing ecological regimes/grouping
            # We let the algorithm decide how many modes we get, however now we need to establish thresholds based on the number of States we chose
            # Now, the idea is that if the number of modes doesn't fit the number of states:
            # If modes > states, then we use quantiles to regroup closely related modes vector in a simple datadriven way (by adding 0 and max(data) in the vector)
            # If we count only zeros as a state, then we substract one quantile prob  # Deprecated : (shunt the first quantile value and start from zero instead of the first quantile value to determine the second state range)
            # if(G >= nbState) {
              if(zeroState == TRUE) {
                thres <- quantile(c(0,exp(mean_tr), max(dataPos)), probs = seq(0,1, 1/(nbState-1))) #max(dataPos or 100 doesnt change anything in computation)
                # thres <- quantile(c(0,exp(mean_tr), 100), probs = seq(0,1, 1/(nbState-1)))
                # thres <- thres[-length(thres)]
              }else{
                thres <- quantile(c(0,exp(mean_tr), max(dataPos)), probs = seq(0,1, 1/(nbState)))
                # thres <- thres[-c(1, length(thres))]
                
              }
              
            # }else{
            #   if(G < nbState){
            #     
            #   thres <- c(0,exp(mean_tr))
            #     
            #   }else{
            #     
            #   }
            # }
            
            
            # Now create a dataframe out of this:
            res <- t(as.data.frame(c("GMM",thres)))
            colnames(res) <- c("Method", paste0("State", seq_len(nbState)))
            rownames(res) <- varName
            
            
          }
          
          
        # QUANTILES 
          # Calculate quantiles to get both GMM and quantiles
          # Also, if GMM fails (not enough data/too many low values), we still have quantiles
          # Check if zero state, if yes then 1 less quantile as zero will be a single value for a state
          if(zeroState == TRUE) {
            thresQTL <- quantile(c(dataPos_tr), probs = seq(0,1, 1/(nbState-1))) 
            thresQTL <- exp(thresQTL)
            # Replace the first quantile by 0 or the zero state, because the 0% is irrelevant to get (for 5 states) :
            # State 0: zero
            # State 1: (0, Q25]
            # State 2: (Q25, Q50]
            # State 3: (Q50, Q75]
            # State 4: (Q75, Q100]
            thresQTL[1] <- 0
            
          }else{
            thresQTL <- quantile(c(dataPos_tr), probs = seq(0,1, 1/(nbState)))
            thresQTL <- exp(thresQTL)
            
          }
          
          # Convert back to "raw values"
          # Small twick : since we're on pos values, change the first value to zero
          # Now create a dataframe out of this:
          resQTL <- t(as.data.frame(c("Quantiles",thresQTL)))
          colnames(resQTL) <- c("method", paste0("State", seq_len(nbState)))
          rownames(resQTL) <- varName
          
          
        # Rbound GMM and quantiles 
        resAll <- as.data.frame(rbind(res, resQTL))
          
        # PLOT
          
          # Only if plot == true then printing the plot 
          if(plot == TRUE) {
            
            if(!is.null(model)){
              
              xgrid <- seq(min(dataPos, na.rm = TRUE),
                           max(dataPos, na.rm = TRUE),
                           length.out = 300)
              
              # Convert to log scale
              xgrid_tr <- log(xgrid)
              
              # Component densities (without zeros)
              # Compute component densities back in raw data space and not log
              dens_df <- lapply(1:G, function(g) {
                # Density(log scale) * Jacobian
                # g = 2
                # Check if equal variance or unequal variance
                if(varType == "V"){
                  density_raw <- pro[g] * dnorm(xgrid_tr, mean_tr[g], sd_tr[g]) * (1/xgrid)
                }else{
                  density_raw <- pro[g] * dnorm(xgrid_tr, mean_tr[g], sd_tr) * (1/xgrid)
                }
                
                tibble(
                  x = xgrid,
                  density = density_raw,
                  component = paste0("Comp ", g)
                )
              }) %>% bind_rows()
              
              # Total mixture density (no zeros)
              dens_total <- dens_df %>%
                group_by(x) %>%
                summarise(density = sum(density)) %>%
                mutate(component = "Mixture")
              
            }
            
            
            data <- as.data.frame(x)
            colnames(data) <- varName
            
            
            temp <- as.data.frame(t(resAll))
            temp <- temp[-1,]
            vlineGMM <- data.frame(ThresholdsGMM = temp[,1])
            vlineGMM$Type <- "Mixture boundaries"
            vlineGMM$ThresholdsGMM <- as.numeric(vlineGMM$ThresholdsGMM)
            
            vlineQTL <- data.frame(ThresholdsQTL = temp[,2])
            vlineQTL$Type <- "Quantiles"
            vlineQTL$ThresholdsQTL <- as.numeric(vlineQTL$ThresholdsQTL)
            
            # rownames(vline_df) <- vline_df$Method
            # vline_df$Method <- NULL
            # vline_df <- as.data.frame(t(vline_df))
            # vline_df$GMM <- as.numeric(vline_df$GMM)
            # vline_df$Quantiles <- as.numeric(vline_df$Quantiles)
            # vline_df$State <- c("Zero","Zero - Low","Low - Intermediate","Intermediate - Good","Good - Very Good")
            # vline_df$Method <- "Quantiles"
            
            thresQTL_Plot <- c(-max(thresQTL)/100, thresQTL)
            state_colors <- c("#d7191c", "orange", "yellow", "#a6d96a", "#1a9641") 
            
            rectDf <- data.frame(xmin = thresQTL_Plot[-length(thresQTL_Plot)],
                                 xmax = thresQTL_Plot[-1],
                                 fillC = state_colors,
                                 state = factor(c("Zero","Low","Intermediate","Good","Very Good"), levels = c("Zero","Low","Intermediate","Good","Very Good"))
                                 )
            
            # Basic plot
            p <- ggplot(data = data, aes(x = !!sym(varName))) +
              # State rects
              geom_rect(data = rectDf,
                        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = 0, fill = state), alpha = 0.8, inherit.aes = FALSE) +
              scale_fill_manual(values = setNames(rectDf$fillC, rectDf$state), name = "State") +
              
              geom_histogram(aes(y = after_stat(density)), fill = "darkgrey", color = "grey35", bins = 50) + #  binwidth = 1
              geom_density(fill = "#69b3a2", color = "lightgrey", alpha = 0.4)
              # geom_histogram(aes(y = ..density..), fill = "darkgrey", color = "grey35")+ # bins = 30
              # geom_density(aes(y = ..scaled..), fill = "#69b3a2", color = "lightgrey", alpha = 0.4)+
            
            # Model
              if(!is.null(model)) {
                components <- sort(unique(as.character(dens_df$component)))
                comp_colors <- setNames(scales::hue_pal()(length(components)), components)
                
                p <- p +
                  ggnewscale::new_scale_fill() +
                  geom_line(data = dens_total,
                          aes(x = x, y = density),
                          linewidth = 1.1, color = "black", alpha = 0.9) +
                  geom_area(data = dens_df, aes(x = x, y = density, fill = component), alpha = 0.2) +
                  geom_line(data = dens_df,
                            aes(x = x, y = density, color = component),
                            linewidth = 0.7, alpha = 0.9, show.legend = TRUE) 
                  # geom_segment(data = vlineGMM, aes(x = ThresholdsGMM, xend = ThresholdsGMM, y = 0, yend = max(density(x)$y), linetype = "Mixture boundaries", color = "Mixture boundaries"), linewidth = 0.7, inherit.aes = FALSE)
                  # geom_segment(data = vlineGMM, aes(x = ThresholdsGMM, xend = ThresholdsGMM, y = 0, yend = max(5)), linetype = "dashed", linewidth = 0.7, color = "lightblue")
                  # geom_vline(xintercept = as.numeric(resAll[1,2:nbState+1]), linetype = "dashed", linewidth = 0.7, color = "lightblue", alpha = 0.7)
                
                
              }else{
                comp_colors <- character(0)
              }
            # Quantiles
              p <- p + geom_segment(data = vlineQTL, aes(x = ThresholdsQTL, xend = ThresholdsQTL, y = 0, yend = max(density(x)$y), linetype = "Quantiles", color = "Quantiles"), linewidth = 0.7, inherit.aes = FALSE)
              # p <- p + geom_segment(data = vline_df, aes(x = Quantiles, xend = Quantiles, y = 0, yend = max(5)), linetype = "dashed", linewidth = 0.7, color = "red") +
              # p <- p + geom_vline(xintercept = as.numeric(resAll[2,2:nbState+1]), linetype = "dashed", linewidth = 0.7, color = "red", alpha = 0.7) +
              
            # Legend 
              # Base mapping: always include Quantiles
              color_values <- c("Quantiles" = "red",
                                # "Mixture boundaries" = "lightblue",
                                comp_colors)
              
              linetype_values <- c("Quantiles" = "dashed") #,
                                   # "Mixture boundaries" = "dotdash")
              
              # Component colors go in fill legend ONLY
              if (!is.null(model)) {
                p <- p +
                  # scale_fill_manual(values = comp_colors, name = "Component") +
                  scale_color_manual(values = color_values, guide = "none") +   # hide color legend
                  guides(
                    fill = guide_legend(
                      override.aes = list(color = comp_colors, linewidth = 1.2)  # show line + fill together
                    )
                  )+
                  # Unified line legend
                  scale_linetype_manual(values = linetype_values, name = NULL) +
                  guides(
                    # linetype = guide_legend(order = 1, override.aes = list(shape = NA, linewidth = 1.2, color = color_values[c("Mixture boundaries", "Quantiles")]))
                    linetype = guide_legend(order = 1, override.aes = list(shape = NA, linewidth = 1.2, color = color_values[c("Quantiles")]))
                  ) +
                  labs(x = varName, y = "Density") +
                  theme_classic()
                
              } else{
                p <- p + scale_color_manual(values = color_values, guide = "none") +
                  # Unified line legend
                  scale_linetype_manual(values = linetype_values, name = NULL) +
                  guides(
                    linetype = guide_legend(order = 1, override.aes = list(shape = NA, linewidth = 1.2, color = color_values[c("Quantiles")]))
                  ) +
                  labs(x = varName, y = "Density") +
                  theme_classic()
              }
              
             

              
              # p <- p +
              # scale_x_continuous(breaks = seq(0,100,10)) +
              # labs(x = varName,
              #      y = "Density",
              #      color = "Component", 
              #      fill = "Component") +
              # theme_classic()
            
            print(p)
            
            if(savePlot == TRUE){
              
              dir.create(savePath, showWarnings = FALSE)
              ggsave(filename = file.path(savePath, paste0("Thresholds_plot_",varName,".pdf")), plot = p, device = "pdf", width = 190, height = 150, units = "mm")
              
            }
            
          }
          
          # exp(mean_tr)
        
        
        
        
        resAllFinal <- cbind(data.frame(Variable = rownames(resAll)), resAll)
        rownames(resAllFinal) <- NULL
        resAllFinal$Variable <- gsub(".1","", resAllFinal$Variable)
        
        # Remove GMM quantiles
        resAllFinal <- resAllFinal[2,]
        
        # Add final warning column if duplicated values between states of a same variable
        resAllFinal <- resAllFinal %>%
          rowwise() %>% 
          mutate(
            warning = if(any(duplicated(c_across(State1:last_col())))) "Warning" else NA
          ) %>%
          ungroup()
        
        
        return(resAllFinal)
        
        
    } 
    
    
    # Function tests
    # makeDDThresholds(coralRorc$RorcHCB, varName = "RorcHCB", nbState = 5, zeroState = TRUE, plot = TRUE, savePlot = FALSE, savePath = file.path(pathPro, "Thresholds","General"))
    makeDDThresholds(coralRorcMean$RorcHCB, varName = "RorcHCB", nbState = 5, zeroState = TRUE, plot = TRUE, savePlot = FALSE, savePath = file.path(pathPro, "Thresholds","General"))
    
    # makeDDThresholds(coralRorc$RorcHCO, varName = "RorcHCO", nbState = 5, zeroState = TRUE, plot = TRUE)
    makeDDThresholds(coralRorcMean$RorcHCO, varName = "RorcHCO", nbState = 5, zeroState = TRUE, plot = TRUE)
    
    # makeDDThresholds(coralRorc$RorcCoralRichness, varName = "RorcSR", nbState = 5, zeroState = TRUE, plot = TRUE)
    makeDDThresholds(coralRorcMean$RorcCoralRichness, varName = "RorcSR", nbState = 5, zeroState = TRUE, plot = TRUE)
    
    # makeDDThresholds(coralRorc$RorcHCM, varName = "RorcHCM", nbState = 5, zeroState = TRUE, plot = TRUE)
    makeDDThresholds(coralRorcMean$RorcHCM, varName = "RorcHCM", nbState = 5, zeroState = TRUE, plot = TRUE)
    
    makeDDThresholds(coralRorcMean$RorcRC, varName = "RorcRC", nbState = 5, zeroState = TRUE, plot = TRUE)
    
    makeDDThresholds(coralRorcMean$RorcCoralCover, varName = "RorcCC", nbState = 5, zeroState = TRUE, plot = TRUE)
    
    makeDDThresholds(fishRorc$RorcCorallivoreAbund, varName = "RorcCorallivores", nbState = 5, zeroState = TRUE, plot = TRUE)
    makeDDThresholds(fishRorc$RorcFishRichness, varName = "RorcFishSR", nbState = 5, zeroState = TRUE, plot = TRUE)
    makeDDThresholds(fishRorc$RorcCarnivoreAbund, varName = "RorcCarnivores", nbState = 5, zeroState = TRUE, plot = TRUE)
    
    makeDDThresholds(invRorc$RorcInvAbund, varName = "RorcInvAbund", nbState = 5, zeroState = TRUE, plot = TRUE)
    makeDDThresholds(invRorc$RorcInvRichness, varName = "RorcInvSR", nbState = 5, zeroState = TRUE, plot = TRUE)
      
      
    # Compute Coral General Thresholds
    thresholdsCoral <- do.call(rbind, lapply(varCoral, function(x, data = coralRorc) {
      
      # x = varNames[1]
      
      res <- makeDDThresholds(data[,x], varName = x, nbState = 5, zeroState = TRUE, plot = TRUE, savePlot = TRUE, savePath = file.path(pathPro, "Thresholds","General"))
      
      return(res)
    }
    ))
      
    
    # Compute Fish General Thresholds
    thresholdsFish <- do.call(rbind, lapply(varFish, function(x, data = fishRorc) {
      
      # x = varNames[1]
      
      res <- makeDDThresholds(data[,x], varName = x, nbState = 5, zeroState = TRUE, plot = TRUE, savePlot = TRUE, savePath = file.path(pathPro, "Thresholds","General"))
      
      return(res)
    }
    ))  
      
      
    # Compute Invertebrates Thresholds 
    
    thresholdsInv <- do.call(rbind, lapply(varInv, function(x, data = invRorc) {
      
      # x = varNames[1]
      
      res <- makeDDThresholds(data[,x], varName = x, nbState = 5, zeroState = TRUE, plot = TRUE, savePlot = TRUE, savePath = file.path(pathPro, "Thresholds","General"))
      
      return(res)
    }
    ))  
    
    
    thresholdsGeneral <- rbind(thresholdsCoral, thresholdsFish, thresholdsInv)
    
    write.csv(thresholdsGeneral, file = file.path(pathPro, "Thresholds","Thresholds_General_DD.csv"), row.names = FALSE)
    
  # Assign Thresholds to data and export new csv to be read for hierarchical model construction
    
    thdGeneral <- read.csv(file.path(pathPro, "Thresholds","Thresholds_General_DD.csv"))
    thdGeneralStates <- c("Zero","Low","Medium","Good","Very_Good")
    colnames(thdGeneral)[3:(length(colnames(thdGeneral))-1)] <- thdGeneralStates
    
    # Uniform Function to assign between different datasets
    assignDDThresholds <- function(colname, data, thds){
      
      # data = coralRorcMean
      # colname = "RorcHCB"
      # thds = thdGeneral
      
      val <- data[[colname]] 
      thd <- thds[thds$Variable == colname,]
      
      res <- character(length = length(val))
      res[val == thd$Zero] <- "Zero"
      res[val > thd$Zero & val <= thd$Low] <- "Low"
      res[val > thd$Low & val <= thd$Medium] <- "medium"
      res[val > thd$Medium & val <= thd$Good] <- "Good"
      res[val > thd$Good] <- "Very_Good"
      
      return(res)
      
    }
    
    # create/duplicate dataset to be filled by states 
    coralRorcGeneralStates <- coralRorcMean
    fishRorcGeneralStates <- fishRorcMean
    invRorcGeneralStates <- invRorcMean
    
    
    coralRorcGeneralStates[varCoral] <- lapply(varCoral, assignDDThresholds, data = coralRorcMean, thds = thdGeneral)
    
    fishRorcGeneralStates[varFish] <- lapply(varFish, assignDDThresholds, data = fishRorcMean, thds = thdGeneral)
    
    invRorcGeneralStates[varInv] <- lapply(varInv, assignDDThresholds, data = invRorcMean, thds = thdGeneral)
    
    # Write to csv
    write.csv(coralRorcGeneralStates, file = file.path(pathPro,"RORC_Coral_States_hdbn.csv"), row.names = FALSE)
    write.csv(fishRorcGeneralStates, file = file.path(pathPro,"RORC_Fish_States_hdbn.csv"), row.names = FALSE)
    write.csv(invRorcGeneralStates, file = file.path(pathPro,"RORC_Inv_States_hdbn.csv"), row.names = FALSE)
    
    
    
    # eo general thresholds ----
    
    
    
    
    
    
  # LOCAL THRESHOLDS
    
    # States here are more "Degrading, "Stagnating, Recovering"
    # In fact its simple for this one : only check value of previous year and consider if it is going back up or stagnating or going down
    # But, to make it robust, we can do a few things :
    # 
    

    # station_df must contain columns: Site, Station, Year, value (numeric)
    # k_recent = number of last years to average (k_recent = 1 => last year only)
    # k_prev   = number of previous years to average (k_prev = 1 => previous year only)
    assign_trend_balanced <- function(st_year_df,
                                      st_year_df_sd,
                                      varName,
                                      k_recent = 1,
                                      k_prev = 1,
                                      pct_threshold = 0.20,      # 20% change
                                      log_threshold = 0.20,      # ~20% multiplicative
                                      se_multiplier = 1.1,       # 1.1 × SE
                                      min_tp = 3,
                                      small = 1e-6) {
      
      # st_year_df = coralRorcMean[coralRorcMean$Station == "casy",]
      # st_year_df_sd = coralRorcSD[coralRorcSD$Station == "casy",]
      # varName = "RorcHCB"
      # varName = "RorcCoralRichness"
      # k_recent = 1
      # k_prev = 1
      # pct_threshold = 0.50     # 10% change
      # log_threshold = 0.50      # ~10% multiplicative
      # se_multiplier = 1.1      # 1 × SE
      # min_tp = 3
      # small = 1e-6
      
      # Simple check just in case
      # must have columns: Site, Station, Year, mean_value, sd_value, n_samples
      st_year_df <- st_year_df %>% arrange(Year)
      st_year_df_sd <- st_year_df_sd %>% arrange(Year)
      
      if (nrow(st_year_df) < min_tp) {
        return(tibble(
          Site = unique(st_year_df$Site),
          Station = unique(st_year_df$Station),
          Variable = varName,
          Trend = NA_character_,
          pct_change = NA_real_,
          log_ratio = NA_real_,
          delta = NA_real_,
          SE_delta = NA_real_,
          reason = NA_character_
        ))
      }
      
      # extract vectors
      vals <- st_year_df[[varName]]
      sds  <- st_year_df_sd[[varName]]
      ns   <- st_year_df$n_samples
      yrs  <- st_year_df$Year
      n    <- length(vals)
      
      # define windows
      recent_idx <- seq(max(1, n - k_recent + 1), n)
      prev_idx   <- seq(max(1, n - k_recent - k_prev + 1), max(1, n - k_recent))
      
      # fallback if too short
      if (length(prev_idx) == 0) prev_idx <- max(1, n - 1)
      
      # compute means
      recent_mean <- mean(vals[recent_idx], na.rm = TRUE)
      prev_mean   <- mean(vals[prev_idx], na.rm = TRUE)
      
      # compute SE of difference
      recent_var <- mean((sds[recent_idx]^2) / pmax(ns[recent_idx], 1), na.rm = TRUE)
      prev_var   <- mean((sds[prev_idx]^2) / pmax(ns[prev_idx], 1),   na.rm = TRUE)
      
      SE_delta <- sqrt(recent_var + prev_var)
      
      # percent change
      denom <- ifelse(prev_mean == 0,
                      max(median(vals[vals>0], na.rm = TRUE) * small, small),
                      prev_mean)
      
      pct_change <- (recent_mean - prev_mean) / denom
      
      # log-ratio
      log_ratio <- log1p(recent_mean) - log1p(prev_mean)
      
      # raw difference
      delta <- recent_mean - prev_mean
      
      # --- decision gates -----------------------------------------------
      gate_pct <- abs(pct_change) >= pct_threshold
      
      gate_log <- abs(log_ratio) >= log_threshold
      
      gate_SE  <- abs(delta) >= se_multiplier * SE_delta
      
      is_up   <- (recent_mean > prev_mean) & (gate_pct | gate_log | gate_SE)
      is_down <- (recent_mean < prev_mean) & (gate_pct | gate_log | gate_SE)
      
      Trend <- dplyr::case_when(
        is_up   ~ "Increasing",
        is_down ~ "Decreasing",
        TRUE    ~ "Stable"
      )
      
      reason <- paste(
        c(
          if (gate_pct) "pct" else NULL,
          if (gate_log) "log" else NULL,
          if (gate_SE)  "SE"  else NULL
        ),
        collapse = "+"
      )
      
      tibble(
        Site = unique(st_year_df$Site),
        Station = unique(st_year_df$Station),
        Variable = varName,
        recent_mean = recent_mean,
        prev_mean = prev_mean,
        pct_change = pct_change,
        log_ratio = log_ratio,
        delta = delta,
        SE_delta = SE_delta,
        Trend = Trend,
        reason = ifelse(reason == "", "none", reason)
      )
    }
    
    
    # Function per variable for all station
    coralRorcStationTrends <- do.call(rbind, lapply())
    
    
    
    
    
    
    # Go see the code back in the trash
    
    # Take the GMM and force 3 states to identify 3 different groups as firstly thought (no quantiles on modes)
    
    # # Function to apply each variable to establish a state
    # 
    # 
    # # Gaussian mixture models :
    # # A data driven way to find thresholds based on each station temporal series.
    # # Gaussian because its the most instinctive distribution
    # # However it is a choice not to pick another distribution
    # # The zero inflated data may bias or deserve another distribution but we'll roll with this "standard" for now and see how it works
    # # Ecologically relevant, works.
    # # A bit quirky because it doesn't consider bound data (notably here cover 0-100) so zero inflated will go in negatives
    # # However, since we're looking at mode means, they cannot be negative and should stay relevant
    
    createStates <- function(x, data, by){

      # TEST ZONE
      x = "RorcHCB"
      data = coralRorcMean

      # x = "RorcCarnivoreAbund"
      # x = "RorcFishRichness"
      # data = fishRorc

      by = c("Year","Station")
      # eo test zone

      # Extracting the data
      df <- data[,c(by,x)]

      # Converting to numeric as a safe quality check (to convert integers to numeric)
      df[[x]] <- as.numeric(df[[x]])
      
      # Keep only positive data
      df <- df[df[[x]] > 0,]

      # Log transform
      # df[[x]] <- log(df[[x]])
      

      # lapply function to rbind all station chunks after being individually evaluated for state generation
      df1 <- do.call(rbind, lapply(unique(df$Station), function(y) {
        # Test zone
        # y = "mbere"
        # y = "mwaremwa"
        # y = "kanga_daa"
        y = "bancs_du_nord"
        # y = "koe"
        # y = "lekiny"
        # eo test zone

        print(y)

        # Extracting the station dataset chunk
        test <- df[df$Station == y,]
        # test <- df


        # Need to check for zeros or very low values
        if(sum(test[[x]]) %in% c(0:1)){
          # If only zeros, apply basic thresholds that will only go to positive values
          # 0% Low
          # 0-30% medium
          # > 30% = High
          # To do in the future : take the threshold of the closest station ?

          # 2 thresholds
          test$classification <- NA
          test$classification[test[[x]] == 0] <- "Low"
          test$classification[test[[x]] > 0 & test[[x]] < 30] <- "Intermediate"
          test$classification[test[[x]] >= 30] <- "High"
          test$classification <- factor(test$classification, levels = c("High","Intermediate","Low"))


        }else{

          # Mixture Model to identify multiple gaussian modes
          model <- Mclust(test[[x]])

          # Get state thresholds
          # If there is no mode resolution (too few data) we refer to the closest station ?
          # If there is only one mode, we chose the data quartiles
          # If there are more than two modes, keep the two extremes and leave the middle to "normal/medium/transitory" state
          if(is.null(model$G)){
            # To change to a warning with manual thresholds
            stop("Gaussian mixture model is null, no solution found. Potentially not enough data.")

          }else{
            if(model$G == 1) {
              # Take the quantile method
              # HOWEVER : make a small adaptation to the number of observation. Notably if there is few data, don't let outliers become the norm
              # Example with RorcHCB on koe (8 observations in total over 2 years (4 transects)). There are 4 zeros, 3 2.5 and one 5.
              # With 0.33-0.66, no values are considered "Normal". Although one could say 2.5 is normal, zero is not, 5 is high
              # So if we have few observation, let's say less than 5 years (20 observations over 4 transects), we use extremes 5% on each side and leave the 95% in "normal"
              # for consistent data, we use 1/3 and 2/3 quantiles OR 1/4 - 3/4
              # After some tests, values are too close one to another to exert the 1/3 or 1/4. Let's try and stick to extremes only
              # if(length(test[[x]]) > 20) {
              #   thres <- quantile(test[[x]], probs = c(0.25, 0.75))
              # } else {
              thres <- quantile(test[[x]], probs = c(0.025, 0.975))
              # }

            }

            if(model$G > 1) {
              # Take the first and last mode means as thresholds and keep all intermediate potential distributions in intermediate
              thres <- c(first(model$parameters$mean), last(model$parameters$mean))
            }
          }
          
          # Convert back to raw values
          # thres <- exp(thres)

          
          
          

          # 2 thresholds
          test$classification <- NA
          test$classification[test[[x]] <= thres[1]] <- "Degrading"
          test$classification[test[[x]] > thres[1] & test[[x]] < thres[2]] <- "Stable"
          test$classification[test[[x]] >= thres[2]] <- "Recovering"
          test$classification <- factor(test$classification, levels = c("Recovering","Stable","Degrading"))
          
          ggplot(data = test, aes(x = Year, y = !!sym(x), group = 1))+
            geom_path() +
            geom_point(aes(color = classification)) +
            scale_color_manual(values = c("#009E73","#56B4E9","#D55E00"), drop = FALSE) +
            facet_wrap(.~Station) +
            geom_hline(yintercept = thres[1], linetype = "dashed", linewidth = 0.1) +
            geom_hline(yintercept = thres[2], linetype = "dashed", linewidth = 0.1) +
            theme_classic()
          
          

        }

        return(test)

      }))


      return(df1)



    }

    res <- createStates(x = "RorcHCB", data = coralRorc, by = c("Year","Station", "Sample"))

  
  
    
    
    
    
    
  
  
  
# TRASH ----
    
    # coralRorcGeneralStates[varCoral] <- lapply(varCoral, function(x){
    #   
    #   # x = varCoral[1]
    #   
    #   val <- coralRorcMean[[x]] 
    #   thd <- thdGeneral[thdGeneral$Variable == x,]
    #   
    #   res <- character(length = length(val))
    #   res[val == thd$Zero] <- "Zero"
    #   res[val > thd$Zero & val <= thd$Low] <- "Low"
    #   res[val > thd$Low & val <= thd$Medium] <- "medium"
    #   res[val > thd$Medium & val <= thd$Good] <- "Good"
    #   res[val > thd$Good] <- "Very_Good"
    #   
    #   return(res)
    #   
    # })
    
  
  
  
  # DDThresholds <- do.call(rbind, lapply(varNames[1:6], function(x, data = coralRorc, nbState = 5, zeroState = TRUE){
  #   
  #   # Test Zone
  #   # x = "RorcSD_SI"
  #   # data = coralRorc
  #   # nbState = 5
  #   # zeroState = TRUE
  #   # eo test zone
  #   print(x)
  #   
  #   # Process data
  #   # Convert to numeric if it isn't already
  #   data[[x]] <- as.numeric(data[[x]])
  #   
  #   dataPos <- data[[x]][data[[x]] > 0]
  #   
  #   dataPos_tr <- log(dataPos)
  #   # Try to work on mean station values instead of directly observed numbers ?
  #   # dat <- data %>%
  #   #   select(-Sample) %>%           # drop the transect column
  #   #   group_by(Year, Country, Region, Sector, Site, Station) %>%   # your metadata columns
  #   #   summarise(across(all_of(varNames), mean), .groups = "drop")
  #   # 
  #   # 
  #   # datSD <- data %>%
  #   #   select(-Sample) %>%           # drop the transect column
  #   #   group_by(Year, Country, Region, Sector, Site, Station) %>%   # your metadata columns
  #   #   summarise(across(all_of(varNames), sd), .groups = "drop")
  #   
  #   
  #   # GAUSSIAN MIXTURE MODEL
  #   
  #   # Run model on positive log transformed data
  #   model <- Mclust(dataPos_tr) #, G = nbState-1  , G = c(1:6)
  #   
  #   # And Quantiles
  #   # thresQTL <- quantile(data[[x]], probs = seq(0,1, 1/nbState))
  #   
  #   # Model checking
  #   # Number of modes found
  #   model$G
  #   # Summary
  #   summary(model)
  #   # Uncertainty plot
  #   # plot(model, what = "uncertainty")
  #   # Uncertainty mean (closer to zero = good)
  #   # mean(model$uncertainty)
  #   # Model modes
  #   model$parameters$mean
  #   exp(model$parameters$mean)
  #   # BIC
  #   model$BIC
  #   
  #   # Try
  #   params <- model$parameters
  #   G <- model$G                     # number of clusters
  #   pro <- params$pro                 # mixing proportions
  #   mean_tr <- params$mean               # means
  #   sd_tr   <- sqrt(params$variance$sigmasq)  # standard deviations (univariate only)
  #   
  #   xgrid <- seq(min(dataPos, na.rm = TRUE),
  #                max(dataPos, na.rm = TRUE),
  #                length.out = 400)
  #   
  #   # Convert to log scale
  #   xgrid_tr <- log(xgrid)
  #   
  #   # Component densities (without zeros)
  #   # Compute component densities back in raw data space and not log
  #   dens_df <- lapply(1:G, function(g) {
  #     # Density(log scale) * Jacobian
  #     density_raw <- pro[g] * dnorm(xgrid_tr, mean_tr[g], sd_tr) * (1/xgrid)
  #     
  #     tibble(
  #       x = xgrid,
  #       density = density_raw,
  #       component = paste0("Comp ", g)
  #     )
  #   }) %>% bind_rows()
  #   
  #   # Total mixture density (no zeros)
  #   dens_total <- dens_df %>%
  #     group_by(x) %>%
  #     summarise(density = sum(density)) %>%
  #     mutate(component = "Mixture")
  #   
  #   ggplot(data = data, aes(x = !!sym(x))) +
  #     geom_histogram(aes(y = ..density..), fill = "darkgrey", color = "grey35", bins = 50)+ #  binwidth = 1
  #     geom_density(fill = "#69b3a2", color = "lightgrey", alpha = 0.4)+
  #     # geom_histogram(aes(y = ..density..), fill = "darkgrey", color = "grey35")+ # bins = 30
  #     # geom_density(aes(y = ..scaled..), fill = "#69b3a2", color = "lightgrey", alpha = 0.4)+
  #     geom_line(data = dens_total,
  #               aes(x = x, y = density),
  #               linewidth = 1.1, color = "black", alpha = 0.9) +
  #     geom_line(data = dens_df,
  #               aes(x = x, y = density, color = component),
  #               linewidth = 0.7, alpha = 0.6) +
  #     geom_vline(xintercept = exp(mean_tr), linetype = "dashed", linewidth = 0.7, color = "lightblue", alpha = 0.7) +
  #     scale_x_continuous(breaks = seq(0,100,10)) +
  #     labs(x = x,
  #          y = "Density",
  #          color = "Component") +
  #     theme_classic()
  #   
  #   exp(mean_tr)
  #   
  #   
  #   
  #   # THRESHOLDS 
  #   
  #   # We've obtained breaking values from GMM, potentially expressing ecological regimes/grouping
  #   # We let the algorithm decide how many modes we get, however now we need to establish thresholds based on the number of States we chose
  #   # Now, the idea is that if the number of modes doesn't fit the number of states:
  #   # If modes > states, then we use quantiles to regroup closely related modes vector in a simple datadriven way (by adding 0 and max(data) in the vector)
  #   # If we count only zeros as a state, then we substract one quantile prob  # Deprecated : (shunt the first quantile value and start from zero instead of the first quantile value to determine the second state range)
  #   if(zeroState == TRUE) {
  #     # thres <- quantile(c(0,exp(mean_tr), max(dataPos)), probs = seq(0,1, 1/(nbState-1))) #max(dataPos or 100 doesnt change anything in computation)
  #     thres <- quantile(c(0,exp(mean_tr), 100), probs = seq(0,1, 1/(nbState-1))) 
  #     # thres <- thres[-length(thres)]
  #   }else{
  #     thres <- quantile(c(0,exp(mean_tr), max(dataPos)), probs = seq(0,1, 1/(nbState)))
  #     # thres <- thres[-c(1, length(thres))]
  #     
  #   }
  #   
  #   # Now create a dataframe out of this:
  #   res <- t(as.data.frame(thres))
  #   colnames(res) <- paste0("State", seq_len(nbState))
  #   
  #   rownames(res) <- x
  #   
  #   return(res)
  #   
  #   
  # } ))
  # 
  # 
  # 
  # DDThresholds
  # 
  # 
  # # Density plots
  # ggplot(data = data, aes(x = !!sym(x))) +
  #   geom_histogram(aes(y = ..density..), bins = 50)+
  #   geom_density(fill = "#69b3a2", color = "grey", alpha = 0.4)+
  #   geom_vline(xintercept = model$parameters$mean, linewidth = 0.4, color = "lightblue") +
  #   # geom_vline(xintercept = thresQTL, linewidth = 0.4, color = "red") +
  #   theme_classic()
  # 
  # # density(data[[x]])
  # par(mfrow = c(1,1))
  # hist(data[[x]], prob = TRUE, breaks = 20, col = "grey80",
  #      main = "", #Coral cover with Gaussian mixture fit
  #      xlab = x)
  # # plot(data[[x]], to.draw.arg="d")
  # plot(model, what = "density", xlab = x, add = TRUE)
  # abline(v = model$parameters$mean, col = "lightblue")
  # 
  # 
  # # GO TO THRESHOLDS
  # # Generate threshold intervals depending on the number of inputed states number
  # # AND on the gaussian mixture model number
  # # Cannot do more than 2 states if only one mode
  # # Cannot do more than three states if 2 modes
  # # Can do less states with mode modes though
  # # Potential errors where :
  # # We want 5 states and there's only 2 modes
  # # We want 3 states but there's only 1 mode
  # # More states...
  # # 
  # 
  # 
  # 
  # # 2 thresholds
  # test$classification <- NA 
  # test$classification[test[[x]] <= thres[1]] <- "Low"
  # test$classification[test[[x]] > thres[1] & test[[x]] < thres[2]] <- "Intermediate"
  # test$classification[test[[x]] >= thres[2]] <- "High"
  # test$classification <- factor(test$classification, levels = c("High","Intermediate","Low"))
  # 
  # 
  # 
  # 
  # 
  # 
  # abline(v = thres, col = "lightblue")
  # 
  # 
  # 
  # 
  # 
  # 
  # # Function to apply each variable to establish a state
  # 
  # 
  # # Gaussian mixture models :
  # # A data driven way to find thresholds based on each station temporal series.
  # # Gaussian because its the most instinctive distribution
  # # However it is a choice not to pick another distribution
  # # The zero inflated data may bias or deserve another distribution but we'll roll with this "standard" for now and see how it works
  # # Ecologically relevant, works.
  # # A bit quirky because it doesn't consider bound data (notably here cover 0-100) so zero inflated will go in negatives
  # # However, since we're looking at mode means, they cannot be negative and should stay relevant
  # 
  # createStates <- function(x, data, by){
  #   
  #   # TEST ZONE 
  #   x = "RorcHCB"
  #   data = coralRorc
  #   
  #   # x = "RorcCarnivoreAbund"
  #   # x = "RorcFishRichness"
  #   # data = fishRorc
  #   
  #   by = c("Year","Station", "Sample")
  #   # eo test zone
  #   
  #   # Extracting the data
  #   df <- data[,c(by,x)]
  #   
  #   # Converting to numeric as a safe quality check (to convert integers to numeric)
  #   df[[x]] <- as.numeric(df[[x]])
  #   
  #   # Normalize data to reduce zero inflated (doesnt change anything because already bounded for covers between 0-100)
  #   # df[[x]] <- scale(df[[x]])
  #   
  #   
  #   # lapply function to rbind all station chunks after being individually evaluated for state generation
  #   df1 <- do.call(rbind, lapply(unique(df$Station), function(y) {
  #     # Test zone
  #     # x = "mbere"
  #     # x = "mwaremwa"
  #     # x = "kanga_daa"
  #     # x = "bancs_du_nord"
  #     # x = "koe"
  #     # y = "lekiny"
  #     # eo test zone
  #     
  #     print(y)
  #     
  #     # Extracting the station dataset chunk
  #     test <- df[df$Station == y,]
  #     test <- df
  #     
  #     
  #     # Need to check for zeros or very low values
  #     if(sum(test[[x]]) %in% c(0:1)){
  #       # If only zeros, apply basic thresholds that will only go to positive values
  #       # 0% Low
  #       # 0-30% medium
  #       # > 30% = High
  #       # To do in the future : take the threshold of the closest station ?
  #       
  #       # 2 thresholds
  #       test$classification <- NA 
  #       test$classification[test[[x]] == 0] <- "Low"
  #       test$classification[test[[x]] > 0 & test[[x]] < 30] <- "Intermediate"
  #       test$classification[test[[x]] >= 30] <- "High"
  #       test$classification <- factor(test$classification, levels = c("High","Intermediate","Low"))
  #       
  #       
  #     }else{
  #       
  #       # Mixture Model to identify multiple gaussian modes 
  #       model <- Mclust(test[[x]])
  #       
  #       # Get state thresholds
  #       # If there is no mode resolution (too few data) we refer to the closest station ?
  #       # If there is only one mode, we chose the data quartiles
  #       # If there are more than two modes, keep the two extremes and leave the middle to "normal/medium/transitory" state
  #       if(is.null(model$G)){
  #         # To change to a warning with manual thresholds
  #         stop("Gaussian mixture model is null, no solution found. Potentially not enough data.")
  #         
  #       }else{
  #         if(model$G == 1) {
  #           # Take the quantile method
  #           # HOWEVER : make a small adaptation to the number of observation. Notably if there is few data, don't let outliers become the norm
  #           # Example with RorcHCB on koe (8 observations in total over 2 years (4 transects)). There are 4 zeros, 3 2.5 and one 5.
  #           # With 0.33-0.66, no values are considered "Normal". Although one could say 2.5 is normal, zero is not, 5 is high
  #           # So if we have few observation, let's say less than 5 years (20 observations over 4 transects), we use extremes 5% on each side and leave the 95% in "normal"
  #           # for consistent data, we use 1/3 and 2/3 quantiles OR 1/4 - 3/4
  #           # After some tests, values are too close one to another to exert the 1/3 or 1/4. Let's try and stick to extremes only
  #           # if(length(test[[x]]) > 20) {
  #           #   thres <- quantile(test[[x]], probs = c(0.25, 0.75))
  #           # } else {
  #           thres <- quantile(test[[x]], probs = c(0.025, 0.975))
  #           # }
  #           
  #         }
  #         
  #         if(model$G > 1) {
  #           # Take the first and last mode means as thresholds and keep all intermediate potential distributions in intermediate
  #           thres <- c(first(model$parameters$mean), last(model$parameters$mean))
  #         }
  #       }
  #       
  #       
  #       # 2 thresholds
  #       test$classification <- NA 
  #       test$classification[test[[x]] <= thres[1]] <- "Low"
  #       test$classification[test[[x]] > thres[1] & test[[x]] < thres[2]] <- "Intermediate"
  #       test$classification[test[[x]] >= thres[2]] <- "High"
  #       test$classification <- factor(test$classification, levels = c("High","Intermediate","Low"))
  #       
  #     }
  #     
  #     return(test)
  #     
  #   }))
  #   
  #   
  #   return(df1)
  #   
  #   
  #   
  # }
  # 
  # res <- createStates(x = "RorcHCB", data = coralRorc, by = c("Year","Station", "Sample"))
  # 
  # 
  # 
  #   # CODE TO GENERATE STATES FROM MIXTURE MODEL THRESHOLDS 
  #   
  #   # Function to apply each variable to establish a state
  #   # Gaussian mixture models :
  #   # A data driven way to find thresholds based on each station temporal series.
  #   # Gaussian because its the most instinctive distribution
  #   # However it is a choice not to pick another distribution
  #   # The zero inflated data may bias or deserve another distribution but we'll roll with this "standard" for now and see how it works
  #   # Ecologically relevant, works.
  #   # A bit quirky because it doesn't consider bound data (notably here cover 0-100) so zero inflated will go in negatives
  #   # However, since we're looking at mode means, they cannot be negative and should stay relevant
  #   
  #   createStates <- function(x, data, by){
  #     
  #     # TEST ZONE 
  #     x = "RorcHCB"
  #     data = coralRorc
  #     
  #     # x = "RorcCarnivoreAbund"
  #     # x = "RorcFishRichness"
  #     # data = fishRorc
  #     
  #     by = c("Year","Station", "Sample")
  #     # eo test zone
  #     
  #     # Extracting the data
  #     df <- data[,c(by,x)]
  #     
  #     # Converting to numeric as a safe quality check (to convert integers to numeric)
  #     df[[x]] <- as.numeric(df[[x]])
  #     
  #     # Normalize data to reduce zero inflated (doesnt change anything because already bounded for covers between 0-100)
  #     # df[[x]] <- scale(df[[x]])
  #     
  #     
  #     # lapply function to rbind all station chunks after being individually evaluated for state generation
  #     df1 <- do.call(rbind, lapply(unique(df$Station), function(y) {
  #       # Test zone
  #       # x = "mbere"
  #       # x = "mwaremwa"
  #       # x = "kanga_daa"
  #       # x = "bancs_du_nord"
  #       # x = "koe"
  #       # y = "lekiny"
  #       # eo test zone
  #       
  #       print(y)
  #       
  #       # Extracting the station dataset chunk
  #       test <- df[df$Station == y,]
  #       test <- df
  #       
  #       
  #       # Need to check for zeros or very low values
  #       if(sum(test[[x]]) %in% c(0:1)){
  #         # If only zeros, apply basic thresholds that will only go to positive values
  #         # 0% Low
  #         # 0-30% medium
  #         # > 30% = High
  #         # To do in the future : take the threshold of the closest station ?
  #         
  #         # 2 thresholds
  #         test$classification <- NA 
  #         test$classification[test[[x]] == 0] <- "Low"
  #         test$classification[test[[x]] > 0 & test[[x]] < 30] <- "Intermediate"
  #         test$classification[test[[x]] >= 30] <- "High"
  #         test$classification <- factor(test$classification, levels = c("High","Intermediate","Low"))
  #         
  #         
  #       }else{
  #         
  #         # Mixture Model to identify multiple gaussian modes 
  #         model <- Mclust(test[[x]])
  #         
  #         # Get state thresholds
  #         # If there is no mode resolution (too few data) we refer to the closest station ?
  #         # If there is only one mode, we chose the data quartiles
  #         # If there are more than two modes, keep the two extremes and leave the middle to "normal/medium/transitory" state
  #         if(is.null(model$G)){
  #           # To change to a warning with manual thresholds
  #           stop("Gaussian mixture model is null, no solution found. Potentially not enough data.")
  #           
  #         }else{
  #           if(model$G == 1) {
  #             # Take the quantile method
  #             # HOWEVER : make a small adaptation to the number of observation. Notably if there is few data, don't let outliers become the norm
  #             # Example with RorcHCB on koe (8 observations in total over 2 years (4 transects)). There are 4 zeros, 3 2.5 and one 5.
  #             # With 0.33-0.66, no values are considered "Normal". Although one could say 2.5 is normal, zero is not, 5 is high
  #             # So if we have few observation, let's say less than 5 years (20 observations over 4 transects), we use extremes 5% on each side and leave the 95% in "normal"
  #             # for consistent data, we use 1/3 and 2/3 quantiles OR 1/4 - 3/4
  #             # After some tests, values are too close one to another to exert the 1/3 or 1/4. Let's try and stick to extremes only
  #             # if(length(test[[x]]) > 20) {
  #             #   thres <- quantile(test[[x]], probs = c(0.25, 0.75))
  #             # } else {
  #             thres <- quantile(test[[x]], probs = c(0.025, 0.975))
  #             # }
  #             
  #           }
  #           
  #           if(model$G > 1) {
  #             # Take the first and last mode means as thresholds and keep all intermediate potential distributions in intermediate
  #             thres <- c(first(model$parameters$mean), last(model$parameters$mean))
  #           }
  #         }
  #         
  #         
  #         # 2 thresholds
  #         test$classification <- NA 
  #         test$classification[test[[x]] <= thres[1]] <- "Low"
  #         test$classification[test[[x]] > thres[1] & test[[x]] < thres[2]] <- "Intermediate"
  #         test$classification[test[[x]] >= thres[2]] <- "High"
  #         test$classification <- factor(test$classification, levels = c("High","Intermediate","Low"))
  #         
  #       }
  #       
  #       return(test)
  #       
  #     }))
  #     
  #     
  #     return(df1)
  #     
  #     
  #     
  #   }
  #   
  #   res <- createStates(x = "RorcHCB", data = coralRorc, by = c("Year","Station", "Sample"))
  #   
  #   
  #   
  #   
  #   
  # # Need to find a way to deal with zero inflated data or more than that : mostly zero data for each category
  # # Usually there are two modes : close to zero and another one further away. These two modes can actually define our three states (before the first thres, between both, after the second thres)
  # # and, run mixture model with any G, but keep the two extremes and assign potential intermediate modes to "normal"
  # # Also need to check that there is at least two distinguished modes
  # # All of this leans toward a hybrid approach ;-)
  # # And convert always integers to numeric because for richness for example it can find 4 gaussians because 4 different richness values (1 - 2 - 3 - 4 for example)
  # # Quartiles for mono gaussian
  # test <- df[df$Station == "mbere",]
  # test <- df[df$Station == "mwaremwa",]
  # test <- df[df$Station == "kanga_daa",]
  # test <- df[df$Station == "bancs_du_nord",]
  # test <- df[df$Station == "koe",]
  # test <- df[df$Station == "lekiny",]
  # 
  # # Get temporal series plot per sample
  # par(mfrow =c(1,1))
  # ggplot(data = test, aes(x = as.factor(Year), y = !!sym(x), group = 1))+
  #   geom_point(aes(color = as.factor(Sample)))+
  #   geom_path(aes(color = as.factor(Sample)))+
  #   facet_wrap(.~Sample) +
  #   theme_classic()
  # 
  # 
  # # Mixture Model to identify multiple gaussian modes 
  # model <- Mclust(test[[x]])
  # model$G
  # # model <- Mclust(test[[x]], G = 2)
  # model$parameters$mean
  # # model$parameters$variance$sigmasq
  # model$BIC
  # summary(model)
  # plot(model, what = "uncertainty")
  # mean(model$uncertainty)
  # # model$z
  # # model$classification
  # # model$data
  # 
  # # Get state thresholds
  # # If there is no mode resolution (too few data) we refer to the closest station ?
  # # If there is only one mode, we chose the data quartiles
  # # If there are more than two modes, keep the two extremes and leave the middle to "normal/medium/transitory" state
  # if(is.null(model$G)){
  #   # To change to a warning with manual thresholds
  #   stop("Gaussian mixture model is null, no solution found. Potentially not enough data.")
  #   
  # }else{
  #   if(model$G == 1) {
  #     # Take the quantile method
  #     # HOWEVER : make a small adaptation to the number of observation. Notably if there is few data, don't let outliers become the norm
  #     # Example with RorcHCB on koe (8 observations in total over 2 years (4 transects)). There are 4 zeros, 3 2.5 and one 5.
  #     # With 0.33-0.66, no values are considered "Normal". Although one could say 2.5 is normal, zero is not, 5 is high
  #     # So if we have few observation, let's say less than 5 years (20 observations over 4 transects), we use extremes 5% on each side and leave the 95% in "normal"
  #     # for consistent data, we use 1/3 and 2/3 quantiles OR 1/4 - 3/4
  #     # After some tests, values are too close one to another to exert the 1/3 or 1/4. Let's try and stick to extremes only
  #     # if(length(test[[x]]) > 20) {
  #     #   thres <- quantile(test[[x]], probs = c(0.25, 0.75))
  #     # } else {
  #     thres <- quantile(test[[x]], probs = c(0.025, 0.975))
  #     # }
  #     
  #   }
  #   
  #   if(model$G > 1) {
  #     # Take the first and last mode means as thresholds and keep all intermediate potential distributions in intermediate
  #     thres <- c(first(model$parameters$mean), last(model$parameters$mean))
  #   }
  # }
  # 
  # thres
  # 
  # # Check histogram and model density distribution
  # par(mfrow = c(1,2))
  # hist(test[[x]], prob = TRUE, breaks = 20, col = "grey80",
  #      main = "", #Coral cover with Gaussian mixture fit
  #      xlab = x)    
  # plot(model, what = "density", add = TRUE, xlab = x)
  # abline(v = thres, col = "lightblue")
  # 
  # # 2 thresholds
  # test$classification <- NA 
  # test$classification[test[[x]] <= thres[1]] <- "Low"
  # test$classification[test[[x]] > thres[1] & test[[x]] < thres[2]] <- "Intermediate"
  # test$classification[test[[x]] >= thres[2]] <- "High"
  # test$classification <- factor(test$classification, levels = c("High","Intermediate","Low"))
  # 
  # # The palette with black:
  # # cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  # 
  # ggplot(data = test, aes(x = Year, y = !!sym(x), group = 1))+
  #   geom_path() +
  #   geom_point(aes(color = classification)) +
  #   scale_color_manual(values = c("#009E73","#56B4E9","#D55E00"), drop = FALSE) +
  #   facet_wrap(.~Station) +
  #   geom_hline(yintercept = thres[1], linetype = "dashed", linewidth = 0.1) +
  #   geom_hline(yintercept = thres[2], linetype = "dashed", linewidth = 0.1) +
  #   theme_classic()
  # 
  # 
  # 
  # 
  # 
  # 
  # # ggplot(data = test, aes(x = Year, y = RorcHCB, group = 1))+
  # #   geom_point(aes(color = as.factor(Sample)))+
  # #   geom_path(aes(color = as.factor(Sample)))+
  # #   theme_classic()
  # 
  # # hist(test$RorcHCB)
  # # temp <- mean(test$RorcHCB)
  # # test$meanInd <- ifelse(test$RorcHCB >= temp, "+","-")
  # # table(test$meanInd)
  # # 
  # # ggplot(data = test, aes(x = Year, y = RorcHCB, group = 1))+
  # #   geom_point(aes(color = as.factor(Sample), shape = meanInd))+
  # #   geom_path(aes(color = as.factor(Sample)))+
  # #   facet_wrap(.~Sample) +
  # #   theme_classic()
  # # 
  # # quantile(test$RorcHCB, probs = c(0, 0.33,0.66,1))
  # # kmeans(test$RorcHCB, 3)
  # # summary(test)
  # 
  # 
  # # Function to apply each variable to establish a state
  # # Gaussian mixture models :
  # # A data driven way to find thresholds based on each station temporal series.
  # # Gaussian because its the most instinctive distribution
  # # However it is a choice not to pick another distribution
  # # The zero inflated data may bias or deserve another distribution but we'll roll with this "standard" for now and see how it works
  # # Ecologically relevant, works.
  # # A bit quirky because it doesn't consider bound data (notably here cover 0-100) so zero inflated will go in negatives
  # # However, since we're looking at mode means, they cannot be negative and should stay relevant
  # 
  # createStates <- function(x, data, by){
  #   
  #   # TEST ZONE 
  #   x = "RorcHCB"
  #   data = coralRorc
  #   
  #   # x = "RorcCarnivoreAbund"
  #   # x = "RorcFishRichness"
  #   # data = fishRorc
  #   
  #   by = c("Year","Station", "Sample")
  #   # eo test zone
  #   
  #   # Extracting the data
  #   df <- data[,c(by,x)]
  #   
  #   # Converting to numeric as a safe quality check (to convert integers to numeric)
  #   df[[x]] <- as.numeric(df[[x]])
  #   
  #   # Normalize data to reduce zero inflated (doesnt change anything because already bounded for covers between 0-100)
  #   # df[[x]] <- scale(df[[x]])
  #   
  #   
  #   # lapply function to rbind all station chunks after being individually evaluated for state generation
  #   df1 <- do.call(rbind, lapply(unique(df$Station), function(y) {
  #     # Test zone
  #     # x = "mbere"
  #     # x = "mwaremwa"
  #     # x = "kanga_daa"
  #     # x = "bancs_du_nord"
  #     # x = "koe"
  #     # y = "lekiny"
  #     # eo test zone
  #     
  #     print(y)
  #     
  #     # Extracting the station dataset chunk
  #     test <- df[df$Station == y,]
  #     test <- df
  #     
  #     
  #     # Need to check for zeros or very low values
  #     if(sum(test[[x]]) %in% c(0:1)){
  #       # If only zeros, apply basic thresholds that will only go to positive values
  #       # 0% Low
  #       # 0-30% medium
  #       # > 30% = High
  #       # To do in the future : take the threshold of the closest station ?
  #       
  #       # 2 thresholds
  #       test$classification <- NA 
  #       test$classification[test[[x]] == 0] <- "Low"
  #       test$classification[test[[x]] > 0 & test[[x]] < 30] <- "Intermediate"
  #       test$classification[test[[x]] >= 30] <- "High"
  #       test$classification <- factor(test$classification, levels = c("High","Intermediate","Low"))
  #       
  #       
  #     }else{
  #       
  #       # Mixture Model to identify multiple gaussian modes 
  #       model <- Mclust(test[[x]])
  #       
  #       # Get state thresholds
  #       # If there is no mode resolution (too few data) we refer to the closest station ?
  #       # If there is only one mode, we chose the data quartiles
  #       # If there are more than two modes, keep the two extremes and leave the middle to "normal/medium/transitory" state
  #       if(is.null(model$G)){
  #         # To change to a warning with manual thresholds
  #         stop("Gaussian mixture model is null, no solution found. Potentially not enough data.")
  #         
  #       }else{
  #         if(model$G == 1) {
  #           # Take the quantile method
  #           # HOWEVER : make a small adaptation to the number of observation. Notably if there is few data, don't let outliers become the norm
  #           # Example with RorcHCB on koe (8 observations in total over 2 years (4 transects)). There are 4 zeros, 3 2.5 and one 5.
  #           # With 0.33-0.66, no values are considered "Normal". Although one could say 2.5 is normal, zero is not, 5 is high
  #           # So if we have few observation, let's say less than 5 years (20 observations over 4 transects), we use extremes 5% on each side and leave the 95% in "normal"
  #           # for consistent data, we use 1/3 and 2/3 quantiles OR 1/4 - 3/4
  #           # After some tests, values are too close one to another to exert the 1/3 or 1/4. Let's try and stick to extremes only
  #           # if(length(test[[x]]) > 20) {
  #           #   thres <- quantile(test[[x]], probs = c(0.25, 0.75))
  #           # } else {
  #           thres <- quantile(test[[x]], probs = c(0.025, 0.975))
  #           # }
  #           
  #         }
  #         
  #         if(model$G > 1) {
  #           # Take the first and last mode means as thresholds and keep all intermediate potential distributions in intermediate
  #           thres <- c(first(model$parameters$mean), last(model$parameters$mean))
  #         }
  #       }
  #       
  #       
  #       # 2 thresholds
  #       test$classification <- NA 
  #       test$classification[test[[x]] <= thres[1]] <- "Low"
  #       test$classification[test[[x]] > thres[1] & test[[x]] < thres[2]] <- "Intermediate"
  #       test$classification[test[[x]] >= thres[2]] <- "High"
  #       test$classification <- factor(test$classification, levels = c("High","Intermediate","Low"))
  #       
  #     }
  #     
  #     return(test)
  #     
  #   }))
  #   
  #   
  #   return(df1)
  #   
  #   
  #   
  # }
  # 
  # res <- createStates(x = "RorcHCB", data = coralRorc, by = c("Year","Station", "Sample"))
  # 



