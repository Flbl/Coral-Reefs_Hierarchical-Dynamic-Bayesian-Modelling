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

  # This script will read the data and convert each biological variable observation per year into a state (High, Normal, Low, Zero), depending on the available data 
  # For that it takes the whole temporal series of the station and checks for trends and mean value, checks the last value and coordinates a choice based on one or each results depending on data availability


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
  
# PROCESS ----
  
  # Create Thresholds from the available data
  # For each variable, make a gaussian Mixture Model and estimate the thresholds. Calculate also the quantiles 
  # Export thresholds to a CSV
  # 1st column is variable name
  # Next columns are threshold1, threshold2, threshold3...
  # Make it threshold number variable
  # (Export two files : "_GMM" for mixture models and "_QTLS" for quantiles)
  # 
  
    # Coral:
    varNames <- colnames(coralRorc)[grep("Rorc", colnames(coralRorc))]  
    # varNames <- colnames(fishRorc)[grep("Rorc", colnames(coralRorc))]  
    # varNames <- colnames(fishRorc)[grep("Rorc", colnames(fishRorc))]  
  
    # Function
     makeDDThresholds <- function(x, nbState, zeroState = TRUE, plot = TRUE, varName = "Variable", save = FALSE, savePath){
      
      # Test Zone
      # x = coralRorc$RorcHCB
      # x = coralRorc$RorcSD_SI
      # x = coralRorc$RorcRC
      # x = coralRorc$RorcCoralRichness
      # nbState = 5
      # zeroState = TRUE
      # plot = TRUE
      # varName = "VariableSR"
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
          }else{
            thresQTL <- quantile(c(dataPos_tr), probs = seq(0,1, 1/(nbState)))

          }
          
          # Convert back to "raw values"
          thresQTL <- exp(thresQTL)
          # Small twick : since we're on pos values, change the first value to zero
          thresQTL[1] <- 0
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
            
            # Basic plot
            p <- ggplot(data = data, aes(x = !!sym(varName))) +
              geom_histogram(aes(y = after_stat(density)), fill = "darkgrey", color = "grey35", bins = 50) + #  binwidth = 1
              geom_density(fill = "#69b3a2", color = "lightgrey", alpha = 0.4)
              # geom_histogram(aes(y = ..density..), fill = "darkgrey", color = "grey35")+ # bins = 30
              # geom_density(aes(y = ..scaled..), fill = "#69b3a2", color = "lightgrey", alpha = 0.4)+
            
            # Model
              if(!is.null(model)) {
                components <- sort(unique(as.character(dens_df$component)))
                comp_colors <- setNames(scales::hue_pal()(length(components)), components)
                
                p <- p +
                geom_line(data = dens_total,
                        aes(x = x, y = density),
                        linewidth = 1.1, color = "black", alpha = 0.9) +
                geom_area(data = dens_df, aes(x = x, y = density, fill = component), alpha = 0.3) +
                geom_line(data = dens_df,
                          aes(x = x, y = density, color = component),
                          linewidth = 0.7, alpha = 0.6, show.legend = TRUE) +
                geom_segment(data = vlineGMM, aes(x = ThresholdsGMM, xend = ThresholdsGMM, y = 0, yend = max(density(x)$y), linetype = "Mixture boundaries", color = "Mixture boundaries"), linewidth = 0.7, inherit.aes = FALSE)
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
                                "Mixture boundaries" = "lightblue",
                                comp_colors)
              
              linetype_values <- c("Quantiles" = "dashed",
                                   "Mixture boundaries" = "dotdash")
              
              # Component colors go in fill legend ONLY
              if (!is.null(model)) {
                p <- p +
                  scale_fill_manual(values = comp_colors, name = "Component") +
                  scale_color_manual(values = color_values, guide = "none") +   # hide color legend
                  guides(
                    fill = guide_legend(
                      override.aes = list(color = comp_colors, linewidth = 1.2)  # show line + fill together
                    )
                  )+
                  # Unified line legend
                  scale_linetype_manual(values = linetype_values, name = NULL) +
                  guides(
                    linetype = guide_legend(order = 1, override.aes = list(shape = NA, linewidth = 1.2, color = color_values[c("Mixture boundaries", "Quantiles")]))
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
            
            if(save == TRUE){
              
              dir.create(savePath, showWarnings = FALSE)
              ggsave(filename = file.path(savePath, paste0("Thresholds_plot_",varName,".pdf")), plot = p, device = "pdf", width = 190, height = 150, units = "mm")
              
            }
            
          }
          
          # exp(mean_tr)
        
        
        
        
        resAllFinal <- cbind(data.frame(Variable = rownames(resAll)), resAll)
        rownames(resAllFinal) <- NULL
        resAllFinal$Variable <- gsub(".1","", resAllFinal$Variable)
        
        return(resAllFinal)
        
        
    } 
  
     
  
  
  makeDDThresholds(coralRorc$RorcHCB, varName = "RorcHCB", nbState = 5, zeroState = TRUE, plot = TRUE, save = TRUE, savePath = file.path(pathPro, "Thresholds","General"))
  makeDDThresholds(coralRorc$RorcHCO, varName = "RorcHCO", nbState = 5, zeroState = TRUE, plot = TRUE)
  makeDDThresholds(coralRorc$RorcCoralRichness, varName = "RorcSR", nbState = 5, zeroState = TRUE, plot = TRUE)
  makeDDThresholds(coralRorc$RorcHCM, varName = "RorcHCM", nbState = 5, zeroState = TRUE, plot = TRUE)
  makeDDThresholds(coralRorc$RorcRC, varName = "RorcRC", nbState = 5, zeroState = TRUE, plot = TRUE)
  makeDDThresholds(coralRorc$RorcCoralCover, varName = "RorcCC", nbState = 5, zeroState = TRUE, plot = TRUE)
  
  makeDDThresholds(fishRorc$RorcCorallivoreAbund, varName = "RorcCorallivores", nbState = 5, zeroState = TRUE, plot = TRUE)
  makeDDThresholds(fishRorc$RorcFishRichness, varName = "RorcFishSR", nbState = 5, zeroState = TRUE, plot = TRUE)
  makeDDThresholds(fishRorc$RorcCarnivoreAbund, varName = "RorcCarnivores", nbState = 5, zeroState = TRUE, plot = TRUE)
  
  makeDDThresholds(invRorc$RorcInvAbund, varName = "RorcInvAbund", nbState = 5, zeroState = TRUE, plot = TRUE)
  makeDDThresholds(invRorc$RorcInvRichness, varName = "RorcInvSR", nbState = 5, zeroState = TRUE, plot = TRUE)
  
  
  # CORAL
  thresholdsCoral <- do.call(rbind, lapply(varNames, function(x, data = coralRorc) {
    
    # x = varNames[1]
    
    res <- makeDDThresholds(data[,x], varName = x, nbState = 5, zeroState = TRUE, plot = TRUE, save = TRUE, savePath = file.path(pathPro, "Thresholds","General"))
    
    return(res)
  }
  ))
    
    
  
  
  
  # FISH
  
  varNames <- colnames(fishRorc)[grep("Rorc", colnames(fishRorc))]
  # varNames <- colnames(fishRorc)[grep("Rorc", colnames(fishRorc))]  
  
  
  thresholdsFish <- do.call(rbind, lapply(varNames, function(x, data = fishRorc) {
    
    # x = varNames[1]
    
    res <- makeDDThresholds(data[,x], varName = x, nbState = 5, zeroState = TRUE, plot = TRUE, save = TRUE, savePath = file.path(pathPro, "Thresholds","General"))
    
    return(res)
  }
  ))  
    
    
    
    
  # Invertebrates
  
  varNames <- colnames(invRorc)[grep("Rorc", colnames(invRorc))]
  # varNames <- colnames(fishRorc)[grep("Rorc", colnames(fishRorc))]  
  
  
  thresholdsInv <- do.call(rbind, lapply(varNames, function(x, data = invRorc) {
    
    # x = varNames[1]
    
    res <- makeDDThresholds(data[,x], varName = x, nbState = 5, zeroState = TRUE, plot = TRUE, save = TRUE, savePath = file.path(pathPro, "Thresholds","General"))
    
    return(res)
  }
  ))  
  
  
  thresholdsGeneral <- rbind(thresholdsCoral, thresholdsFish, thresholdsInv)
  
  write.csv(thresholdsGeneral, file = file.path(pathPro, "Thresholds","Thresholds_General_DD.csv"), row.names = FALSE)
    
    
# TRASH ----
  
  
  
  
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
  #   facet_wrap(.~Sample) +
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



