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
  # c("Year","Country","Region","Sector","Site","Stationt","Sample")

  # It has biological data with columns starting with the dataset name
  # ex: RorcCoralRichness
  # This helps identifying meta from biological data

  # Later : maybe separate meta from biological obs and keep only one id column to cross with meta
  # (Slightly dangerous due to needing to always update both datasets accordingly)

  # This script will read the data and convert each biological variable observation per year into a state (High, Normal, Low), depending on the available data 
  # For that it takes the whole temporal serie of the station and checks for trends and mean value, checks the last value and coordinates a choice based on one or each results depending on data availability


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
  
  
# Function to apply each variable to establish a state
  # Gaussian mixture models :
  # A data driven way to find thresholds based on each station temporal series.
  # Gaussian because its the most instinctive distribution
  # However it is a choice not to pick another distribution
  # The zero inflated data may bias or deserve another distribution but we'll roll with this "standard" for now and see how it works
  # Ecologically relevant, works.
  # A bit quirky because it doesn't consider bound data (notably here cover 0-100) so zero inflated will go in negatives
  # However, since we're looking at mode means, they cannot be negative and should stay relevant
  
  createStates <- function(x, data, by){
    
    # TEST ZONE 
      # x = "RorcHCB"
      # data = coralRorc
      
      # x = "RorcCarnivoreAbund"
      # x = "RorcFishRichness"
      # data = fishRorc
      
      # by = c("Year","Station", "Sample")
    # eo test zone
    
    # Extracting the data
    df <- data[,c(by,x)]
    
    # Converting to numeric as a safe quality check (to convert integers to numeric)
    df[[x]] <- as.numeric(df[[x]])
    
    
    # lapply function to rbind all station chunks after being individually evaluated for state generation
    df1 <- do.call(rbind, lapply(unique(df$Station), function(y) {
      # Test zone
      # x = "mbere"
      # x = "mwaremwa"
      # x = "kanga_daa"
      # x = "bancs_du_nord"
      # x = "koe"
      # y = "lekiny"
      # eo test zone
      
      print(y)
      
      # Extracting the station dataset chunk
      test <- df[df$Station == y,]
      
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
        
        
        # 2 thresholds
        test$classification <- NA 
        test$classification[test[[x]] <= thres[1]] <- "Low"
        test$classification[test[[x]] > thres[1] & test[[x]] < thres[2]] <- "Intermediate"
        test$classification[test[[x]] >= thres[2]] <- "High"
        test$classification <- factor(test$classification, levels = c("High","Intermediate","Low"))
        
      }
      
      return(test)
      
    }))
    
    
    return(df1)
    
    
    
  }
  
  res <- createStates(x = "RorcHCB", data = coralRorc, by = c("Year","Station", "Sample"))
  
# TRASH ----
  
  # Need to find a way to deal with zero inflated data or more than that : mostly zero data for each category
  # Usually there are two modes : close to zero and another one further away. These two modes can actually define our three states (before the first thres, between both, after the second thres)
  # and, run mixture model with any G, but keep the two extremes and assign potential intermediate modes to "normal"
  # Also need to check that there is at least two distinguished modes
  # All of this leans toward a hybrid approach ;-)
  # And convert always integers to numeric because for richness for example it can find 4 gaussians because 4 different richness values (1 - 2 - 3 - 4 for example)
  # Quartiles for mono gaussian
  test <- df[df$Station == "mbere",]
  test <- df[df$Station == "mwaremwa",]
  test <- df[df$Station == "kanga_daa",]
  test <- df[df$Station == "bancs_du_nord",]
  test <- df[df$Station == "koe",]
  test <- df[df$Station == "lekiny",]
  
  # Get temporal series plot per sample
  par(mfrow =c(1,1))
  ggplot(data = test, aes(x = as.factor(Year), y = !!sym(x), group = 1))+
    geom_point(aes(color = as.factor(Sample)))+
    geom_path(aes(color = as.factor(Sample)))+
    facet_wrap(.~Sample) +
    theme_classic()
  
  
  # Mixture Model to identify multiple gaussian modes 
  model <- Mclust(test[[x]])
  model$G
  # model <- Mclust(test[[x]], G = 2)
  model$parameters$mean
  # model$parameters$variance$sigmasq
  # model$BIC
  summary(model)
  plot(model, what = "uncertainty")
  mean(model$uncertainty)
  # model$z
  # model$classification
  # model$data
  
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
  
  thres
  
  # Check histogram and model density distribution
  par(mfrow = c(1,2))
  hist(test[[x]], prob = TRUE, breaks = 20, col = "grey80",
       main = "", #Coral cover with Gaussian mixture fit
       xlab = x)    
  plot(model, what = "density", add = TRUE, xlab = x)
  abline(v = thres, col = "lightblue")
  
  # 2 thresholds
  test$classification <- NA 
  test$classification[test[[x]] <= thres[1]] <- "Low"
  test$classification[test[[x]] > thres[1] & test[[x]] < thres[2]] <- "Intermediate"
  test$classification[test[[x]] >= thres[2]] <- "High"
  test$classification <- factor(test$classification, levels = c("High","Intermediate","Low"))
  
  # The palette with black:
  # cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  ggplot(data = test, aes(x = Year, y = !!sym(x), group = 1))+
    geom_path() +
    geom_point(aes(color = classification)) +
    scale_color_manual(values = c("#009E73","#56B4E9","#D55E00"), drop = FALSE) +
    facet_wrap(.~Sample) +
    geom_hline(yintercept = thres[1], linetype = "dashed", linewidth = 0.1) +
    geom_hline(yintercept = thres[2], linetype = "dashed", linewidth = 0.1) +
    theme_classic()
  
  
  
  
  
  
  # ggplot(data = test, aes(x = Year, y = RorcHCB, group = 1))+
  #   geom_point(aes(color = as.factor(Sample)))+
  #   geom_path(aes(color = as.factor(Sample)))+
  #   theme_classic()
  
  # hist(test$RorcHCB)
  # temp <- mean(test$RorcHCB)
  # test$meanInd <- ifelse(test$RorcHCB >= temp, "+","-")
  # table(test$meanInd)
  # 
  # ggplot(data = test, aes(x = Year, y = RorcHCB, group = 1))+
  #   geom_point(aes(color = as.factor(Sample), shape = meanInd))+
  #   geom_path(aes(color = as.factor(Sample)))+
  #   facet_wrap(.~Sample) +
  #   theme_classic()
  # 
  # quantile(test$RorcHCB, probs = c(0, 0.33,0.66,1))
  # kmeans(test$RorcHCB, 3)
  # summary(test)
  
  




