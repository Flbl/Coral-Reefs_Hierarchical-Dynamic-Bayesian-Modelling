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

  dir.create(file.path(pathProSpe,"Thresholds"), showWarnings = FALSE)
  dir.create(file.path(pathProSpe,"Thresholds","General"), showWarnings = FALSE)
  

  # eo INIT ----



# READ DATA ----

# RORC ----

  # Coral
  coralRorc <- read.csv(file.path(pathSpe,"RORC_Coral_Data_hdbn.csv"))
  # coralRorc$Station <- gsub(" ","_", coralRorc$Station)
  # coralRorc$Site <- gsub(" ","_", coralRorc$Site)
  
  # Fish
  fishRorc <- read.csv(file.path(pathSpe,"RORC_Fish_Data_hdbn.csv"))
  
  # Inv
  invRorc <- read.csv(file.path(pathSpe,"RORC_Inv_Data_hdbn.csv"))


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
    group_by(Year, Country, Region, Sector, Site, Station) %>%   # your metadata columns
    summarise(
      n_samples = n_distinct(Sample),
      across(all_of(varCoral), sd),
      .groups = "drop"
    )

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
    # select(-Sample) %>%           # drop the transect column
    group_by(Year, Country, Region, Sector, Site, Station) %>%   # your metadata columns
    summarise(
      n_samples = n_distinct(Sample),
      across(all_of(varFish), sd),
      .groups = "drop"
    )
  
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
    # select(-Sample) %>%           # drop the transect column
    group_by(Year, Country, Region, Sector, Site, Station) %>%   # your metadata columns
    summarise(
      n_samples = n_distinct(Sample),
      across(all_of(varInv), sd),
      .groups = "drop"
    )
  
  
# GENERAL THRESHOLD ----
  
  # All data is used to calulate data driven thresholds using GMM and quantiles to establish global states of quality
  # This assumes that there is enough data in the temporal series to distinguish the "health" gradient from zero/bad to High/very High
  
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
      # x = coralRorc$RorcFS
      # nbState = 5
      # zeroState = TRUE
      # plot = TRUE
      # varName = "VariableFS"
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
            # Uncertainty mean (closer to zero = High)
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
            # vline_df$State <- c("Zero","Zero - Low","Low - Intermediate","Intermediate - High","High - Very High")
            # vline_df$Method <- "Quantiles"
            
            thresQTL_Plot <- c(-max(thresQTL)/100, thresQTL)
            state_colors <- c("#d7191c", "orange", "yellow", "#a6d96a", "#1a9641") 
            
            rectDf <- data.frame(xmin = thresQTL_Plot[-length(thresQTL_Plot)],
                                 xmax = thresQTL_Plot[-1],
                                 fillC = state_colors,
                                 state = factor(c("Zero","Low","Intermediate","High","Very High"), levels = c("Zero","Low","Intermediate","High","Very High"))
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
    makeDDThresholds(coralRorcMean$RorcHCB, varName = "RorcHCB", nbState = 5, zeroState = TRUE, plot = TRUE, savePlot = FALSE, savePath = file.path(pathProSpe, "Thresholds","General"))
    
    # makeDDThresholds(coralRorc$RorcHCO, varName = "RorcHCO", nbState = 5, zeroState = TRUE, plot = TRUE)
    makeDDThresholds(coralRorcMean$RorcHCO, varName = "RorcHCO", nbState = 5, zeroState = TRUE, plot = TRUE)
    
    makeDDThresholds(coralRorcMean$RorcFS, varName = "RorcFS", nbState = 5, zeroState = TRUE, plot = TRUE)
    
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
    makeDDThresholds(invRorc$RorcUrchinAbund, varName = "RorcUrchinAbund", nbState = 5, zeroState = TRUE, plot = TRUE)
    
      
    # Compute Coral General Thresholds
    thresholdsCoral <- do.call(rbind, lapply(varCoral, function(x, data = coralRorc) {
      
      # x = varNames[1]
      
      res <- makeDDThresholds(data[,x], varName = x, nbState = 5, zeroState = TRUE, plot = TRUE, savePlot = TRUE, savePath = file.path(pathProSpe, "Thresholds","General"))
      
      return(res)
    }
    ))
      
    
    # Compute Fish General Thresholds
    thresholdsFish <- do.call(rbind, lapply(varFish, function(x, data = fishRorc) {
      
      # x = varNames[1]
      
      res <- makeDDThresholds(data[,x], varName = x, nbState = 5, zeroState = TRUE, plot = TRUE, savePlot = TRUE, savePath = file.path(pathProSpe, "Thresholds","General"))
      
      return(res)
    }
    ))  
      
      
    # Compute Invertebrates Thresholds 
    
    thresholdsInv <- do.call(rbind, lapply(varInv, function(x, data = invRorc) {
      
      # x = varNames[1]
      
      res <- makeDDThresholds(data[,x], varName = x, nbState = 5, zeroState = TRUE, plot = TRUE, savePlot = TRUE, savePath = file.path(pathProSpe, "Thresholds","General"))
      
      return(res)
    }
    ))  
    
    
    thresholdsGeneral <- rbind(thresholdsCoral, thresholdsFish, thresholdsInv)
    
    # To be modified manually when there's warnings
    write.csv(thresholdsGeneral, file = file.path(pathProSpe, "Thresholds","Thresholds_General_DD.csv"), row.names = FALSE)
    
    
  # Assign Thresholds to data and export new csv to be read for hierarchical model construction
    
    thdGeneral <- read.csv(file.path(pathProSpe, "Thresholds","Thresholds_General_DD.csv"))
    thdGeneralStates <- c("Zero","Low","Medium","High","Very_High")
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
      res[val > thd$Low & val <= thd$Medium] <- "Medium"
      res[val > thd$Medium & val <= thd$High] <- "High"
      res[val > thd$High] <- "Very_High"
        
      
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
    write.csv(coralRorcGeneralStates, file = file.path(pathProSpe,"RORC_Coral_Station_General_States_hdbn.csv"), row.names = FALSE)
    write.csv(fishRorcGeneralStates, file = file.path(pathProSpe,"RORC_Fish_Station_General_States_hdbn.csv"), row.names = FALSE)
    write.csv(invRorcGeneralStates, file = file.path(pathProSpe,"RORC_Inv_Station_General_States_hdbn.csv"), row.names = FALSE)
    
    
    
    # eo general thresholds ----
    
    
    
    
    
    
  # LOCAL THRESHOLDS
    
    # States here are more "Degrading, "Stagnating, Recovering"
    # In fact its simple for this one : only check value of previous year and consider if it is going back up or stagnating or going down
    
    # /!\  ADD IF DIFF = 0 AT T0 FOR CUM2 THEN STABLE, DONT KEEP ONLY THE ADDITION OF ONE PART
    
    assign_trend_rules <- function(dfMean,
                                   dfSD,
                                   varName,
                                   hist_mult = 1.0,       # multiplier for historical sd gate (single year)
                                   se_mult = 1.0,         # multiplier for SE gate (single year)
                                   cum_mult = 1.0,        # multiplier for hist_sd for cumulative gate
                                   min_tp = 3,           # minimum years to attempt any classification
                                   cum_window = 2        # how many consecutive diffs to accumulate (2 = t + t-1)
    ) {
      
      # Test zone
      # dfMean = coralRorcMean
      # dfSD = coralRorcSD
      # varName = varCoral[1]
      # hist_mult = 1.0
      # se_mult = 1.0
      # cum_mult = 1.0
      # min_tp = 3
      # cum_window = 2
      
      
      # Join mean and sd info. Expect dfMean has Year, Country, Region, Sector, Site, Station, n_samples, and varName present.
      df <- dfMean %>%
        select(Year, Country, Region, Sector, Site, Station, n_samples, mean_value = !!rlang::sym(varName)) %>%
        left_join(
          dfSD %>% select(Year, Site, Station, sd_value = !!rlang::sym(varName)),
          by = c("Year", "Site", "Station")
        ) %>%
        arrange(Country, Region, Sector, Site, Station, Year)
      
      # group_modify: operate per (Site, Station)
      out <- df %>%
        group_by(Site, Station) %>%
        arrange(Year) %>%
        group_modify(~ {
          d <- .x
          # number of years available
          n <- nrow(d)
          
          # initialize columns
          d$dif <- NA_real_
          d$SE_delta <- NA_real_
          d$hist_sd <- NA_real_
          d$cum2 <- NA_real_
          d$Trend <- NA_character_
          
          # if not enough years, set Trend NA (or "Stable") and return the group's df
          if (n < 2) {
            d$Trend <- NA_character_
            return(d)
          }
          
          # Year-to-year diffs of means: dif[1] = NA, dif[i] = mean_value[i] - mean_value[i-1]
          difs <- c(NA_real_, diff(d$mean_value))
          d$dif <- difs
          
          # SE of mean for each year: SE = sd / sqrt(n_samples)
          # If sd_value or n_samples missing, produce NA
          SE_year <- rep(NA_real_, n)
          valid_idx <- which(!is.na(d$sd_value) & !is.na(d$n_samples) & d$n_samples > 0)
          if (length(valid_idx) > 0) {
            SE_year[valid_idx] <- d$sd_value[valid_idx] / sqrt(d$n_samples[valid_idx])
          }
          
          # SE of difference between year t and t-1 = sqrt(SE_t^2 + SE_t-1^2)
          SE_delta <- rep(NA_real_, n)
          for (i in 2:n) {
            if (!is.na(SE_year[i]) && !is.na(SE_year[i-1])) {
              SE_delta[i] <- sqrt(SE_year[i]^2 + SE_year[i-1]^2)
            } else {
              SE_delta[i] <- NA_real_
            }
          }
          d$SE_delta <- SE_delta
          
          # historical sd of diffs excluding the first NA
          hist_sd <- sd(difs[-1], na.rm = TRUE)
          if (is.na(hist_sd) || hist_sd == 0) hist_sd <- 0
          d$hist_sd <- hist_sd
          
          
          # cumulative window of diffs (sign-consistent only)
          cum_vals <- rep(NA_real_, n)
          
          if (cum_window <= 1) {
            cum_vals <- difs
          } else {
            for (i in 1:n) {
              idx <- seq(i - (cum_window - 1), i)
              idx <- idx[idx >= 2 & idx <= n]
              
              if (length(idx) == cum_window) {
                window_vals <- difs[idx]
                prev <- window_vals[1]
                curr <- window_vals[2]
                
                # --- NEW RULES ---
                if (curr == 0) {
                  cum_vals[i] <- 0                # momentum broken
                  next
                }
                
                # curr != 0 → cumulative = curr
                cum_vals[i] <- curr
                next
              }
              
              cum_vals[i] <- NA_real_
            }
          }
          
          # Old version
          # if (cum_window <= 1) {
          #   cum_vals <- difs
          # } else {
          #   for (i in 1:n) {
          #     idx <- seq(i - (cum_window - 1), i)
          #     idx <- idx[idx >= 2 & idx <= n]
          #     
          #     # only compute cumulative sum if we have exactly cum_window difs
          #     if (length(idx) == cum_window) {
          #       window_vals <- difs[idx]
          #       
          #       if (any(window_vals == 0)) {
          #         cum_vals[i] <- 0
          #         next
          #       }
          #       
          #       # check if all difs have same sign (ignoring zeros)
          #       signs <- sign(window_vals)
          #       signs <- signs[signs != 0]     # remove zeros; zero ruins sign logic but contains no direction
          #       
          #       if (length(signs) == 0) {
          #         # all zeros → cumulative = 0, but can never trip a gate, safe to keep
          #         cum_vals[i] <- 0
          #       } else if (all(signs == signs[1])) {
          #         # consistent signs → sum is meaningful
          #         cum_vals[i] <- sum(window_vals)
          #       } else {
          #         # inconsistent directional signals → ignore cumulative window
          #         cum_vals[i] <- NA_real_
          #       }
          #       
          #     } else {
          #       cum_vals[i] <- NA_real_
          #     }
          #   }
          # }
          
          
          d$cum2 <- cum_vals
          
          # Determine Trend per year (compare year i to year i-1)
          Trend <- rep("Stable", n)
          Trend[1] <- NA_character_
          
          for (i in 2:n) {
            dif_i <- difs[i]
            se_i <- SE_delta[i]
            cum_i <- cum_vals[i]
            
            # single-year gates: historical sd gate and SE gate
            gate_hist <- if (hist_sd == 0) FALSE else abs(dif_i) > hist_mult * hist_sd
            gate_se   <- if (is.na(se_i)) FALSE else abs(dif_i) > se_mult * se_i
            gate_single <- gate_hist & gate_se
            
            # cumulative gate (compare cumulative sum to cumulative multiplier * hist_sd)
            gate_cum <- if (is.na(cum_i) || hist_sd == 0) FALSE else abs(cum_i) > cum_mult * hist_sd
            
            # classification: require either single-year gate OR cumulative gate
            if (gate_single) {
              if (dif_i > 0) Trend[i] <- "Increasing" else if (dif_i < 0) Trend[i] <- "Decreasing" else Trend[i] <- "Stable"
            } else if (gate_cum) {
              if (cum_i > 0) Trend[i] <- "Increasing" else if (cum_i < 0) Trend[i] <- "Decreasing" else Trend[i] <- "Stable"
            } else {
              Trend[i] <- "Stable"
            }
          }
          
          d$Trend <- Trend
          
          # enforce min_tp: if group has fewer than min_tp years, set Trend to NA
          if (n < min_tp) d$Trend <- NA_character_
          
          d
        }) %>%
        ungroup()
      
      return(out)
    }
    
    # test <- assign_trend_rules(coralRorcMean, coralRorcSD, "RorcCoralCover")
    test <- assign_trend_rules(coralRorcMean, coralRorcSD, "RorcHCB")
    # test <- assign_trend_rules(coralRorcMean, coralRorcSD, "RorcHCO")
    # test$Trend
    
    unique(test$Station)
    test1 <- test[test$Station == "maitre",]
    test1 <- test[test$Station == "signal",]
    test1 <- test[test$Station == "casy",]
    test1 <- test[test$Station == "baie_des_citrons",]
    
    test1$Trend <- ifelse(is.na(test1$Trend), "NA", test1$Trend)
    
    ggplot(data = test1, aes(x = Year, y = mean_value)) +
      geom_hline(yintercept = c(0, 2.5,7.5,17.5), colour = c("red","orange","yellow","lightgreen"), linetype = "dashed") +
      geom_line()+
      geom_errorbar(aes(ymin=mean_value-sd_value, ymax=mean_value+sd_value), width=.2,
                    position=position_dodge(0.05), color = "grey", linetype = "dashed") +
      geom_point(aes(color = Trend))+
      scale_color_manual(
        values = c(
          "Decreasing"  = "#D55E00",
          "Increasing" = "#009E73",
          "Stable"     = "#56B4E9",
          "NA"         = "#999999"   # now a real category → legend appears
        ),
        limits = c("Decreasing","Increasing","Stable","NA"),
        drop = FALSE
      ) +
      theme_classic()
    
    
    # Compute for all variables
    pathTrends <- file.path(pathProSpe,"Thresholds","Trends")
    dir.create(pathTrends, showWarnings = FALSE)
    # Coral
    lapply(varCoral, function(x) {
      dataVar <- assign_trend_rules(coralRorcMean, coralRorcSD, x)
      # cb <- dataVar[,c("Site","Station","Year","Country","Region","Sector","n_samples","Trend")]
      # colnames(cb[colnames(cb) == "Trend"])
      write.csv(dataVar, file = file.path(pathTrends, paste0("RORC_Coral_", x, "_Station_Trends_States_hdbn.csv")), row.names = FALSE)
    })
    # Fish
    lapply(varFish, function(x) {
      dataVar <- assign_trend_rules(fishRorcMean, fishRorcSD, x)
      write.csv(dataVar, file = file.path(pathTrends, paste0("RORC_Fish_", x, "_Station_Trends_States_hdbn.csv")), row.names = FALSE)
    })
    
    # Inv
    lapply(varInv, function(x) {
      dataVar <- assign_trend_rules(invRorcMean, invRorcSD, x)
      write.csv(dataVar, file = file.path(pathTrends, paste0("RORC_Inv_", x, "_Station_Trends_States_hdbn.csv")), row.names = FALSE)
    })
    
    
    # Wrapper to concatenate all trends per variables
    # Coral
    coralMeanTrend <- coralRorcMean
    meta_cols <- c("Year","Country","Region","Sector","Site","Station","n_samples")
    var_list <- varCoral
    
    for(v in var_list){
      trendsdf <- assign_trend_rules(coralRorcMean, coralRorcSD, v)
      trendsdf <- trendsdf[,c(meta_cols, "Trend")]
      # rename Trend column to the variable name so join can overwrite cleanly
      trendsdf <- trendsdf %>% 
        rename(!!v := Trend)
      
      # merge by the meta columns, without depending on row order
      coralMeanTrend <- coralMeanTrend %>%
        select(-all_of(v)) %>%       # remove old values of v
        left_join(trendsdf, by = meta_cols)
      
    }    
    write.csv(coralMeanTrend, file = file.path(pathProSpe, paste0("RORC_Coral_Station_Trends_States_hdbn.csv")), row.names = FALSE)
    
    # Fish
    fishMeanTrend <- fishRorcMean
    meta_cols <- c("Year","Country","Region","Sector","Site","Station","n_samples")
    var_list <- varFish
    
    for(v in var_list){
      trendsdf <- assign_trend_rules(fishRorcMean, fishRorcSD, v)
      trendsdf <- trendsdf[,c(meta_cols, "Trend")]
      # rename Trend column to the variable name so join can overwrite cleanly
      trendsdf <- trendsdf %>% 
        rename(!!v := Trend)
      
      # merge by the meta columns, without depending on row order
      fishMeanTrend <- fishMeanTrend %>%
        select(-all_of(v)) %>%       # remove old values of v
        left_join(trendsdf, by = meta_cols)
      
    }    
    write.csv(fishMeanTrend, file = file.path(pathProSpe, paste0("RORC_Fish_Station_Trends_States_hdbn.csv")), row.names = FALSE)
    
    
    # Inv
    invMeanTrend <- invRorcMean
    meta_cols <- c("Year","Country","Region","Sector","Site","Station","n_samples")
    var_list <- varInv
    
    for(v in var_list){
      trendsdf <- assign_trend_rules(invRorcMean, invRorcSD, v)
      trendsdf <- trendsdf[,c(meta_cols, "Trend")]
      # rename Trend column to the variable name so join can overwrite cleanly
      trendsdf <- trendsdf %>% 
        rename(!!v := Trend)
      
      # merge by the meta columns, without depending on row order
      invMeanTrend <- invMeanTrend %>%
        select(-all_of(v)) %>%       # remove old values of v
        left_join(trendsdf, by = meta_cols)
      
    }    
    write.csv(invMeanTrend, file = file.path(pathProSpe, paste0("RORC_Inv_Station_Trends_States_hdbn.csv")), row.names = FALSE)
  
  
    
# some plots for presentations
    # Read general thresholds
    # coralMean <- coralRorcMean
    # coralSD <- coralRorcSD
    # coralStates <- read.csv(file.path(pathProSpe,"RORC_Coral_Station_General_States_hdbn.csv"))
    
    
    
    # Read evolution
  
    
    
    
    
    
    
    
# TRASH ----
    # # Go see the code back in the trash
    # 
    # # Take the GMM and force 3 states to identify 3 different groups as firstly thought (no quantiles on modes)
    # 
    # # # Function to apply each variable to establish a state
    # # 
    # # 
    # # # Gaussian mixture models :
    # # # A data driven way to find thresholds based on each station temporal series.
    # # # Gaussian because its the most instinctive distribution
    # # # However it is a choice not to pick another distribution
    # # # The zero inflated data may bias or deserve another distribution but we'll roll with this "standard" for now and see how it works
    # # # Ecologically relevant, works.
    # # # A bit quirky because it doesn't consider bound data (notably here cover 0-100) so zero inflated will go in negatives
    # # # However, since we're looking at mode means, they cannot be negative and should stay relevant
    # 
    # createStates <- function(x, data, by){
    # 
    #   # TEST ZONE
    #   x = "RorcHCB"
    #   data = coralRorcMean
    # 
    #   # x = "RorcCarnivoreAbund"
    #   # x = "RorcFishRichness"
    #   # data = fishRorc
    # 
    #   by = c("Year","Station")
    #   # eo test zone
    # 
    #   # Extracting the data
    #   df <- data[,c(by,x)]
    # 
    #   # Converting to numeric as a safe quality check (to convert integers to numeric)
    #   df[[x]] <- as.numeric(df[[x]])
    #   
    #   # Keep only positive data
    #   df <- df[df[[x]] > 0,]
    # 
    #   # Log transform
    #   # df[[x]] <- log(df[[x]])
    #   
    # 
    #   # lapply function to rbind all station chunks after being individually evaluated for state generation
    #   df1 <- do.call(rbind, lapply(unique(df$Station), function(y) {
    #     # Test zone
    #     # y = "mbere"
    #     # y = "mwaremwa"
    #     # y = "kanga_daa"
    #     y = "bancs_du_nord"
    #     # y = "koe"
    #     # y = "lekiny"
    #     # eo test zone
    # 
    #     print(y)
    # 
    #     # Extracting the station dataset chunk
    #     test <- df[df$Station == y,]
    #     # test <- df
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
    #       # Convert back to raw values
    #       # thres <- exp(thres)
    # 
    #       
    #       
    #       
    # 
    #       # 2 thresholds
    #       test$classification <- NA
    #       test$classification[test[[x]] <= thres[1]] <- "Degrading"
    #       test$classification[test[[x]] > thres[1] & test[[x]] < thres[2]] <- "Stable"
    #       test$classification[test[[x]] >= thres[2]] <- "Recovering"
    #       test$classification <- factor(test$classification, levels = c("Recovering","Stable","Degrading"))
    #       
    #       ggplot(data = test, aes(x = Year, y = !!sym(x), group = 1))+
    #         geom_path() +
    #         geom_point(aes(color = classification)) +
    #         scale_color_manual(values = c("#009E73","#56B4E9","#D55E00"), drop = FALSE) +
    #         facet_wrap(.~Station) +
    #         geom_hline(yintercept = thres[1], linetype = "dashed", linewidth = 0.1) +
    #         geom_hline(yintercept = thres[2], linetype = "dashed", linewidth = 0.1) +
    #         theme_classic()
    #       
    #       
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
    # 
    # 
    
    # # Permutational test of the mean between year x and year x-1 for each variable
    # # The function takes the values of a variable column from the dataframe at sample scale
    # # For each station, the function takes the whole temporal series and checks if there is enough data (2 series)
    # # It computes the sd of differences on the whole station's temporal series
    # # rbind Lapply in each station for each year, that runs the lmperm
    # # Check the significance, check the sd difference, if both are good
    # 
    # # Lmperm function
    # perm_p <- function(values_year_tm1, values_year_t) {
    #   # Test zone
    #   # values_year_tm1 = df2[df2$Year == yr-1,y]
    #   # values_year_t = df2[df2$Year == yr,y]
    #   values_year_tm1 = c(0,10,20,40)
    #   values_year_t = c(0,0,12.5,7.5)
    #   
    #   df <- tibble(
    #     value = c(values_year_tm1, values_year_t),
    #     year  = factor(c(rep("Tminus1", length(values_year_tm1)),
    #                      rep("T", length(values_year_t))))
    #   )
    #   
    #   fit <- lmp(value ~ year, data = df)
    #   summary(fit)$coefficients["year1", "Pr(Exact)"]
    # }
    # 
    # # Test perm_p
    # perm_p(c(0,1), c(0,0))
    # perm_p(c(0,0), c(0,0))
    # perm_p(c(0), c(0))
    # 
    # 
    # # Function
    # assignDDStationTrends <- function(st, data, varNames){ # x is the station name
    #   
    #   # Test zone
    #   st = "casy"
    #   data = coralRorc
    #   varNames = varCoral
    #   
    #   # Remember, the data should be structured at least with Year, Station and Samples outside the varNames
    #   
    #   # Extract the df
    #   df <- data[data$Station == st,]
    #   
    #   # Compute per-year means + sample counts
    #   dfMean <- df %>%
    #     group_by(Year, Station) %>%
    #     summarise(
    #       n_samples = n_distinct(Sample),
    #       across(all_of(varNames), mean),
    #       .groups = "drop"
    #     ) %>%
    #     arrange(Year)
    #   
    #   years <- dfMean$Year
    #   
    #   # ---- Historical SD-based gates per variable ----
    #   hist_sd_change <- sapply(varNames, function(v){
    #     x <- dfMean[[v]]
    #     
    #     dx <- diff(x)
    #     if (length(dx) < 2) {
    #       sdChange <- sd(x, na.rm = TRUE)
    #     } else {
    #       sdChange <- sd(dx, na.rm = TRUE)
    #     }
    #     
    #     if (is.na(sdChange)) sdChange <- 0
    #     return(sdChange)
    #   })
    #   
    #   
    #   # Run the function for each column of varNames
    #   # ---- For each variable, compute rolling year-to-year trend ----
    #   for (v in varNames) {
    #     # v = "RorcHCB"
    #     
    #     sdGate <- hist_sd_change[[v]]
    #     
    #     res <- sapply(seq_along(years), function(i){
    #       # i = 2
    #       
    #       if (i == 1) return(NA_character_)
    #       
    #       Tm1 <- df[df$Year == years[i-1], v]
    #       T0  <- df[df$Year == years[i],   v]
    #       
    #       dif <- mean(T0) - mean(Tm1)
    #       
    #       case_when(
    #         dif >  sdGate ~ "Recovering",
    #         dif < -sdGate ~ "Degrading",
    #         TRUE          ~ "Stable"
    #       )
    #     })
    #     
    #     # add the column to dfMean
    #     dfMean[[paste0(v, "_trend")]] <- res
    #   }
    #   
    #   
    #   
    #   return(dfMean)
    # }
    # 
    # 
    #   dfMean[varNames] <- lapply(dfMean[varNames], function(var, dat = df) {
    #     # Test
    #     # var = dfMean[varNames][2]
    # 
    #     # Define colname
    #     varCol <- colnames(var)
    #     
    #     print(varCol)
    #     
    #     # Recreate df with variable
    #     df1 <- dat[,c("Year", "Station", "Sample", varCol)]
    #     
    #     # isolate the sd threshold gate
    #     sdGate <- hist_sd_change[varCol]
    #     
    #     # Isolate first year
    #     firstYear <- df1$Year[1]
    #     
    #     # Compute yearly permutation tests of mean (no test for the first value), 
    #     # check the pvalue and if positive, compare to sd threshold to estimate state
    #     res <- unlist(lapply(unique(df1$Year), function(yr){
    #       # yr = 2014
    #       # yr = 2022
    #       
    #       if(yr == firstYear){trend <- NA
    #       }else{
    #         
    #         # df2 <- df1[df1$Year %in% c(yr-1, yr),]
    #         Tm1 <- df1[df1$Year == yr-1, varCol]
    #         T0 <- df1[df1$Year == yr, varCol]
    #         # re <- perm_p(Tm1, T0)
    #         dif <- mean(T0) - mean(Tm1)
    #         
    #         # 
    #         trend <- dplyr::case_when(
    #           # dif > sdGate && dif > 0 && re <= 0.05  ~ "Recovering",
    #           # dif < -sdGate && dif < 0 && re <= 0.05 ~ "Degrading",
    #           dif > sdGate && dif > 0  ~ "Recovering",
    #           dif < -sdGate && dif < 0 ~ "Degrading",
    #           TRUE ~ "Stable"
    #         )
    #         
    #         return(trend)
    #         
    #       }
    #       
    #     }))
    #     
    #     # ggplot(data = dfMean, aes(x = Year, y = !!sym(varCol))) +
    #     #   geom_line()+
    #     #   geom_point(aes(color = res))+
    #     #   theme_classic()
    #     
    #     return(res)
    #     
    #   })
    #   
    #   return(dfMean)
    #   
    # } 
    # 
    # # Then I do call rbind so I can replicate it simply for each dataset
    # assignDDStationTrends("casy", coralRorc, varCoral)
    # coralRorcStationTrends <- do.call(rbind, lapply(unique(coralRorc$Station), assignDDStationTrends, data = coralRorc, varNames = varCoral))
    # 
    # 
    # 
    # # permutation test for difference in means (two-sample)
    # # Brute force
    # # perm_test_mean <- function(x, y, B = 2000, seed = 1) {
    # #   set.seed(seed)
    # #   obs <- mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)
    # #   pooled <- c(x, y)
    # #   n_x <- length(x)
    # #   diffs <- replicate(B, {
    # #     s <- sample(pooled, length(pooled), replace = FALSE)
    # #     mean(s[1:n_x], na.rm = TRUE) - mean(s[(n_x+1):length(pooled)], na.rm = TRUE)
    # #   })
    # #   p_two <- mean(abs(diffs) >= abs(obs))
    # #   list(stat = obs, p_perm = p_two)
    # # }
    # 
    # # Lmperm
    # perm_p <- function(values_year_t, values_year_tm1) {
    #   df <- tibble(
    #     value = c(values_year_t, values_year_tm1),
    #     year  = factor(c(rep("T", length(values_year_t)),
    #                      rep("Tminus1", length(values_year_tm1))))
    #   )
    #   
    #   fit <- lmp(value ~ year, data = df)
    #   summary(fit)$coefficients["yearT", "Pr(Prob)"]
    # }
    # 
    # 
    # # ----------------------------------------------
    # # 2. Define rolling windows
    # # ----------------------------------------------
    # recent_idx <- seq(max(1, n - k_recent + 1), n)
    # prev_idx   <- seq(max(1, n - k_recent - k_prev + 1),
    #                   max(1, n - k_recent))
    # if (length(prev_idx) == 0) prev_idx <- max(1, n - 1)
    # 
    # recent_vals <- vals[recent_idx]
    # prev_vals   <- vals[prev_idx]
    # recent_sds  <- sds[recent_idx]
    # prev_sds    <- sds[prev_idx]
    # 
    # recent_mean <- mean(recent_vals, na.rm = TRUE)
    # prev_mean   <- mean(prev_vals,   na.rm = TRUE)
    # 
    # denom <- ifelse(prev_mean == 0,
    #                 max(median(vals[vals>0], na.rm=TRUE) * small, small),
    #                 prev_mean)
    # 
    # pct_change <- (recent_mean - prev_mean) / denom
    # log_ratio  <- log1p(recent_mean) - log1p(prev_mean)
    # 
    # # ----------------------------------------------
    # # 3. Permutation test via lmperm::lmp()
    # # ----------------------------------------------
    # perm_df <- tibble(
    #   value = c(prev_vals, recent_vals),
    #   year  = factor(c(rep("prev", length(prev_vals)),
    #                    rep("recent", length(recent_vals))))
    # )
    # 
    # perm_fit <- try(lmp(value ~ year, data = perm_df), silent = TRUE)
    # 
    # perm_p <- if (inherits(perm_fit, "try-error")) {
    #   NA_real_
    # } else {
    #   summary(perm_fit)$coefficients["yearrecent", "Pr(Prob)"]
    # }
    # 
    # 
    # # ----------------------------------------------
    # # 4. Historical relevance gate
    # # ----------------------------------------------
    # hist_sd_change <- sd(diff(vals), na.rm = TRUE)
    # if (is.na(hist_sd_change)) hist_sd_change <- 0
    # 
    # effect_size <- abs(recent_mean - prev_mean)
    # effect_gate <- effect_size > (hist_mult * hist_sd_change)
    # 
    # 
    # 
    # 
    # assign_trend_adaptive <- function(station_df,
    #                                   station_df_sd,
    #                                   varName,
    #                                   k_recent = 1,
    #                                   k_prev = 1,
    #                                   pct_q = 0.70,
    #                                   log_sd_mult = 0.5,
    #                                   min_tp = 3,
    #                                   p_thresh = 0.10,        # permutation p-value threshold
    #                                   hist_mult = 0.5,        # historical-effect-size threshold multiplier
    #                                   use_SE_gate = TRUE,     # optional
    #                                   small = 1e-6) {
    #   
    #   # Test Zone
    #   station_df = coralRorcMean[coralRorcMean$Station == "casy",]
    #   station_df_sd = coralRorcSD[coralRorcSD$Station == "casy",]
    #   varName = "RorcHCB"
    #   k_recent = 1
    #   k_prev = 1
    #   pct_q = 0.70
    #   log_sd_mult = 0.5
    #   min_tp = 2
    #   p_thresh = 0.10       # permutation p-value threshold
    #   hist_mult = 0.5        # historical-effect-size threshold multiplier
    #   use_SE_gate = TRUE     # optional
    #   small = 1e-6
    #   
    #   
    #   station_df <- station_df %>% arrange(Year)
    #   station_df_sd <- station_df_sd %>% arrange(Year)
    #   
    #   vals <- station_df[[varName]]
    #   sds  <- station_df_sd[[varName]]
    #   yrs  <- station_df$Year
    #   n    <- length(vals)
    #   
    #   if (n < min_tp) {
    #     return(tibble(
    #       Site = unique(station_df$Site),
    #       Station = unique(station_df$Station),
    #       n_tp = n,
    #       recent_mean = NA_real_,
    #       prev_mean = NA_real_,
    #       pct_change = NA_real_,
    #       log_ratio = NA_real_,
    #       hist_pct_thresh = NA_real_,
    #       hist_log_thresh = NA_real_,
    #       perm_p = NA_real_,
    #       effect_size = NA_real_,
    #       hist_sd_change = NA_real_,
    #       Trend = NA_character_
    #     ))
    #   }
    #   
    #   # ----------------------------------------------
    #   # 1. Historical percent and log adaptive thresholds
    #   # ----------------------------------------------
    #   pct_changes_hist <- numeric(n-1)
    #   for (i in 2:n) {
    #     prev <- vals[i-1]
    #     denom <- ifelse(prev == 0,
    #                     max(median(vals[vals>0], na.rm=TRUE) * small, small),
    #                     prev)
    #     pct_changes_hist[i-1] <- (vals[i] - prev) / denom
    #   }
    #   
    #   abs_pct_hist <- abs(pct_changes_hist)
    #   hist_pct_thresh <- quantile(abs_pct_hist, pct_q, na.rm = TRUE)
    #   if (is.na(hist_pct_thresh) || hist_pct_thresh == 0) hist_pct_thresh <- 0.05
    #   
    #   log_hist <- diff(log1p(vals))
    #   hist_log_sd <- sd(log_hist, na.rm = TRUE)
    #   if (is.na(hist_log_sd) || hist_log_sd == 0) hist_log_sd <- 0.01
    #   
    #   hist_log_thresh <- log_sd_mult * hist_log_sd
    #   
    #   # ----------------------------------------------
    #   # 2. Define rolling windows
    #   # ----------------------------------------------
    #   recent_idx <- seq(max(1, n - k_recent + 1), n)
    #   prev_idx   <- seq(max(1, n - k_recent - k_prev + 1),
    #                     max(1, n - k_recent))
    #   if (length(prev_idx) == 0) prev_idx <- max(1, n - 1)
    #   
    #   recent_vals <- vals[recent_idx]
    #   prev_vals   <- vals[prev_idx]
    #   recent_sds  <- sds[recent_idx]
    #   prev_sds    <- sds[prev_idx]
    #   
    #   recent_mean <- mean(recent_vals, na.rm = TRUE)
    #   prev_mean   <- mean(prev_vals,   na.rm = TRUE)
    #   
    #   denom <- ifelse(prev_mean == 0,
    #                   max(median(vals[vals>0], na.rm=TRUE) * small, small),
    #                   prev_mean)
    #   
    #   pct_change <- (recent_mean - prev_mean) / denom
    #   log_ratio  <- log1p(recent_mean) - log1p(prev_mean)
    #   
    #   # ----------------------------------------------
    #   # 3. Permutation test via lmperm::lmp()
    #   # ----------------------------------------------
    #   perm_df <- tibble(
    #     value = c(prev_vals, recent_vals),
    #     year  = factor(c(rep("prev", length(prev_vals)),
    #                      rep("recent", length(recent_vals))))
    #   )
    #   
    #   perm_fit <- try(lmp(value ~ year, data = perm_df), silent = TRUE)
    #   
    #   perm_p <- if (inherits(perm_fit, "try-error")) {
    #     NA_real_
    #   } else {
    #     summary(perm_fit)$coefficients["yearrecent", "Pr(Prob)"]
    #   }
    #   
    #   # ----------------------------------------------
    #   # 4. Historical relevance gate
    #   # ----------------------------------------------
    #   hist_sd_change <- sd(diff(vals), na.rm = TRUE)
    #   if (is.na(hist_sd_change)) hist_sd_change <- 0
    #   
    #   effect_size <- abs(recent_mean - prev_mean)
    #   effect_gate <- effect_size > (hist_mult * hist_sd_change)
    #   
    #   # ----------------------------------------------
    #   # 5. Optional SE gate (for within-year replication)
    #   # ----------------------------------------------
    #   if (use_SE_gate) {
    #     recent_var <- mean(recent_sds^2, na.rm = TRUE) / length(recent_sds)
    #     prev_var   <- mean(prev_sds^2,   na.rm = TRUE) / length(prev_sds)
    #     SE_delta <- sqrt(recent_var + prev_var)
    #     SE_gate <- effect_size > SE_delta
    #   } else {
    #     SE_gate <- TRUE
    #   }
    #   
    #   # ----------------------------------------------
    #   # 6. Combine gates + adaptive thresholds → Trend
    #   # ----------------------------------------------
    #   # Gate A: significance (low p-value)
    #   sig_gate <- !is.na(perm_p) && perm_p < p_thresh
    #   
    #   # Gate B: station historical relevance
    #   hist_gate <- effect_gate
    #   
    #   # Gate C: optional SE restraint
    #   final_gate <- sig_gate && hist_gate && SE_gate
    #   
    #   # Adaptive directionality
    #   is_up   <- pct_change >  hist_pct_thresh && log_ratio >  hist_log_thresh
    #   is_down <- pct_change < -hist_pct_thresh && log_ratio < -hist_log_thresh
    #   
    #   Trend <- dplyr::case_when(
    #     final_gate && is_up   ~ "Recovering",
    #     final_gate && is_down ~ "Degrading",
    #     TRUE ~ "Stable"
    #   )
    #   
    #   # ----------------------------------------------
    #   # 7. Return summary
    #   # ----------------------------------------------
    #   tibble(
    #     Site = unique(station_df$Site),
    #     Station = unique(station_df$Station),
    #     n_tp = n,
    #     recent_mean = recent_mean,
    #     prev_mean = prev_mean,
    #     pct_change = pct_change,
    #     log_ratio = log_ratio,
    #     hist_pct_thresh = hist_pct_thresh,
    #     hist_log_thresh = hist_log_thresh,
    #     perm_p = perm_p,
    #     effect_size = effect_size,
    #     hist_sd_change = hist_sd_change,
    #     SE_gate = SE_gate,
    #     sig_gate = sig_gate,
    #     hist_gate = hist_gate,
    #     Trend = Trend
    #   )
    # }
    # 
    # 
    # 
    # 
    # # station_df must contain columns: Site, Station, Year, value (numeric)
    # # k_recent = number of last years to average (k_recent = 1 => last year only)
    # # k_prev   = number of previous years to average (k_prev = 1 => previous year only)
    # assign_trend_adaptive <- function(station_df,
    #                                   station_df_sd,
    #                                   varName,
    #                                   # year,
    #                                   k_recent = 1,
    #                                   k_prev = 1,
    #                                   pct_q = 0.70,   # percentile for pct_change threshold
    #                                   log_sd_mult = 0.5, # multiplier for SD of historical log changes
    #                                   min_tp = 3,
    #                                   small = 1e-6) {
    #   
    #   
    #   # Test zone
    #   # station_df = coralRorcMean[coralRorcMean$Station == "casy",]
    #   # station_df_sd = coralRorcSD[coralRorcSD$Station == "casy",]
    #   # varName = "RorcHCB"
    #   # # varName = "RorcCoralRichness"
    #   # k_recent = 1
    #   # k_prev = 1
    #   # pct_q = 0.70
    #   # log_sd_mult = 0.5
    #   # min_tp = 3
    #   # small = 1e-6
    #   
    #   
    #   station_df <- station_df %>% arrange(Year)
    #   station_df_sd <- station_df_sd %>% arrange(Year)
    #   # Filter by the available year?
    #   # station_df <- station_df[station_df$Year <= year]
    #   # station_df_sd <- station_df_sd[station_df_sd$Year <= year]
    #   
    #   vals <- station_df[[varName]]
    #   sds <- station_df_sd[[varName]]
    #   yrs  <- station_df$Year
    #   n <- length(vals)
    #   
    #   
    #   # --- return NA results if too few time points ---
    #   if (n < min_tp) {
    #     return(tibble(
    #       Site = unique(station_df$Site),
    #       Station = unique(station_df$Station),
    #       n_tp = n,
    #       recent_mean = NA_real_,
    #       prev_mean = NA_real_,
    #       pct_change = NA_real_,
    #       log_ratio = NA_real_,
    #       hist_pct_thresh = NA_real_,
    #       hist_log_thresh = NA_real_,
    #       Trend = NA_character_
    #     ))
    #   }
    #   
    #   # ------------------------------
    #   # 1. Historical PERCENT CHANGES
    #   # ------------------------------
    #   pct_changes_hist <- c()
    #   # Take all possible differences between value and n-1 value (2:n) 
    #   for (i in 2:n) {
    #     # i = 2
    #     prev <- vals[i-1]
    #     denom <- ifelse(prev == 0,
    #                     max(median(vals[vals>0], na.rm = TRUE) * small, small), # check that
    #                     prev)
    #     pct_changes_hist[i-1] <- (vals[i] - prev) / denom
    #   }
    #   
    #   abs_pct_hist <- abs(pct_changes_hist)
    #   
    #   # Adaptive threshold: 70th percentile (default)
    #   hist_pct_thresh <- quantile(abs_pct_hist, pct_q, na.rm = TRUE)
    #   if (is.na(hist_pct_thresh) || hist_pct_thresh == 0) {
    #     hist_pct_thresh <- 0.05  # fallback small value of 5%
    #   }
    #   
    #   # ------------------------------
    #   # 2. Historical LOG CHANGES
    #   # ------------------------------
    #   log_hist <- diff(log1p(vals))
    #   hist_log_sd <- sd(log_hist, na.rm = TRUE)
    #   if (is.na(hist_log_sd) || hist_log_sd == 0) hist_log_sd <- 0.01 # small fall back expressing a very small variability for single values
    #   
    #   # Adaptive log threshold:
    #   hist_log_thresh <- log_sd_mult * hist_log_sd
    #   
    #   # ------------------------------
    #   # 3. Define recent & previous windows
    #   # ------------------------------
    #   recent_idx <- seq(max(1, n - k_recent + 1), n)
    #   prev_idx   <- seq(max(1, n - k_recent - k_prev + 1),
    #                     max(1, n - k_recent))
    #   
    #   if (length(prev_idx) == 0) prev_idx <- max(1, n - 1)
    #   
    #   recent_vals <- vals[recent_idx]
    #   prev_vals   <- vals[prev_idx]
    #   
    #   # Mean only if more than one year for each window
    #   recent_mean <- mean(recent_vals, na.rm = TRUE)
    #   prev_mean   <- mean(prev_vals, na.rm = TRUE)
    #   
    #   denom <- ifelse(prev_mean == 0,
    #                   max(median(vals[vals > 0], na.rm = TRUE) * small, small),
    #                   prev_mean)
    #   
    #   pct_change <- (recent_mean - prev_mean) / denom
    #   log_ratio  <- log1p(recent_mean) - log1p(prev_mean)
    #   
    #   # ------------------------------
    #   # 4. Classification using ADAPTIVE thresholds
    #   # ------------------------------
    #   is_up <- (pct_change >  hist_pct_thresh) &
    #     (log_ratio  >  hist_log_thresh)
    #   
    #   is_down <- (pct_change < -hist_pct_thresh) &
    #     (log_ratio  < -hist_log_thresh)
    #   
    #   Trend <- dplyr::case_when(
    #     is_up   ~ "Recovering",
    #     is_down ~ "Degrading",
    #     TRUE    ~ "Stable"
    #   )
    #   
    #   # ------------------------------
    #   # 5. Return results
    #   # ------------------------------
    #   tibble(
    #     Site = unique(station_df$Site),
    #     Station = unique(station_df$Station),
    #     n_tp = n,
    #     recent_mean = recent_mean,
    #     prev_mean = prev_mean,
    #     pct_change = pct_change,
    #     log_ratio = log_ratio,
    #     hist_pct_thresh = hist_pct_thresh,
    #     hist_log_thresh = hist_log_thresh,
    #     Trend = Trend
    #   )
    # }
    # 
    # 
    # assign_trend_adaptive(station_df = coralRorcMean[coralRorcMean$Station == "casy",], 
    #                       station_df_sd = coralRorcSD[coralRorcSD$Station == "casy",],
    #                       varName = "RorcHCB"
    #                       )
    # 
    # 
    # assign_trend_balanced <- function(st_year_df,
    #                                   st_year_df_sd,
    #                                   varName,
    #                                   k_recent = 1,
    #                                   k_prev = 1,
    #                                   pct_threshold = 0.20,      # 20% change
    #                                   log_threshold = 0.20,      # ~20% multiplicative
    #                                   se_multiplier = 1.1,       # 1.1 × SE
    #                                   min_tp = 3,
    #                                   small = 1e-6) {
    #   
    #   # st_year_df = coralRorcMean[coralRorcMean$Station == "casy",]
    #   # st_year_df_sd = coralRorcSD[coralRorcSD$Station == "casy",]
    #   # varName = "RorcHCB"
    #   # varName = "RorcCoralRichness"
    #   # k_recent = 1
    #   # k_prev = 1
    #   # pct_threshold = 0.50     # 10% change
    #   # log_threshold = 0.50      # ~10% multiplicative
    #   # se_multiplier = 1.1      # 1 × SE
    #   # min_tp = 3
    #   # small = 1e-6
    #   
    #   # Simple check just in case
    #   # must have columns: Site, Station, Year, mean_value, sd_value, n_samples
    #   st_year_df <- st_year_df %>% arrange(Year)
    #   st_year_df_sd <- st_year_df_sd %>% arrange(Year)
    #   
    #   if (nrow(st_year_df) < min_tp) {
    #     return(tibble(
    #       Site = unique(st_year_df$Site),
    #       Station = unique(st_year_df$Station),
    #       Variable = varName,
    #       Trend = NA_character_,
    #       pct_change = NA_real_,
    #       log_ratio = NA_real_,
    #       delta = NA_real_,
    #       SE_delta = NA_real_,
    #       reason = NA_character_
    #     ))
    #   }
    #   
    #   # extract vectors
    #   vals <- st_year_df[[varName]]
    #   sds  <- st_year_df_sd[[varName]]
    #   ns   <- st_year_df$n_samples
    #   yrs  <- st_year_df$Year
    #   n    <- length(vals)
    #   
    #   # define windows
    #   recent_idx <- seq(max(1, n - k_recent + 1), n)
    #   prev_idx   <- seq(max(1, n - k_recent - k_prev + 1), max(1, n - k_recent))
    #   
    #   # fallback if too short
    #   if (length(prev_idx) == 0) prev_idx <- max(1, n - 1)
    #   
    #   # compute means
    #   recent_mean <- mean(vals[recent_idx], na.rm = TRUE)
    #   prev_mean   <- mean(vals[prev_idx], na.rm = TRUE)
    #   
    #   # compute SE of difference
    #   recent_var <- mean((sds[recent_idx]^2) / pmax(ns[recent_idx], 1), na.rm = TRUE)
    #   prev_var   <- mean((sds[prev_idx]^2) / pmax(ns[prev_idx], 1),   na.rm = TRUE)
    #   
    #   SE_delta <- sqrt(recent_var + prev_var)
    #   
    #   # percent change
    #   denom <- ifelse(prev_mean == 0,
    #                   max(median(vals[vals>0], na.rm = TRUE) * small, small),
    #                   prev_mean)
    #   
    #   pct_change <- (recent_mean - prev_mean) / denom
    #   
    #   # log-ratio
    #   log_ratio <- log1p(recent_mean) - log1p(prev_mean)
    #   
    #   # raw difference
    #   delta <- recent_mean - prev_mean
    #   
    #   # --- decision gates -----------------------------------------------
    #   gate_pct <- abs(pct_change) >= pct_threshold
    #   
    #   gate_log <- abs(log_ratio) >= log_threshold
    #   
    #   gate_SE  <- abs(delta) >= se_multiplier * SE_delta
    #   
    #   is_up   <- (recent_mean > prev_mean) & (gate_pct | gate_log | gate_SE)
    #   is_down <- (recent_mean < prev_mean) & (gate_pct | gate_log | gate_SE)
    #   
    #   Trend <- dplyr::case_when(
    #     is_up   ~ "Increasing",
    #     is_down ~ "Decreasing",
    #     TRUE    ~ "Stable"
    #   )
    #   
    #   reason <- paste(
    #     c(
    #       if (gate_pct) "pct" else NULL,
    #       if (gate_log) "log" else NULL,
    #       if (gate_SE)  "SE"  else NULL
    #     ),
    #     collapse = "+"
    #   )
    #   
    #   tibble(
    #     Site = unique(st_year_df$Site),
    #     Station = unique(st_year_df$Station),
    #     Variable = varName,
    #     recent_mean = recent_mean,
    #     prev_mean = prev_mean,
    #     pct_change = pct_change,
    #     log_ratio = log_ratio,
    #     delta = delta,
    #     SE_delta = SE_delta,
    #     Trend = Trend,
    #     reason = ifelse(reason == "", "none", reason)
    #   )
    # }
    # 
    # 
    # # Function per variable for all station
    # coralRorcStationTrends <- do.call(rbind, lapply())
    # 
    # 
    # 
    # 
    
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



