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
  
  createStates <- function(x, data, by){
    
    # TEST ZONE 
      x = "RorcHCB"
      data = coralRorc
      x = "RorcCarnivoreAbund"
      x = "RorcFishRichness"
      data = fishRorc
      
      by = c("Year","Station", "Sample")
    # eo test zone
    
    
    df <- data[,c(by,x)]
    
    # Need to find a way to deal with zero inflated data or more than that : mostly zero data for each category
    
    # Usually there are two modes : close to zero and another one further away. These two modes can actually define our three states (before the first thres, between both, after the second thres)
    # and, run mixture model with any G, but keep the two extremes and assign potential intermediate modes to "normal"
    # Also need to check that there is at least two distinguished modes
    # All of this leans toward a hybrid approach ;-)
    test <- df[df$Station == "mbere",]
    test <- df[df$Station == "mwaremwa",]
    test <- df[df$Station == "kanga_daa",]
    test <- df[df$Station == "bancs_du_nord",]
    
    # Get temporal series plot per sample
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
    model$parameters$variance$sigmasq
    # model$z
    # model$classification
    # model$data
    
    # Get state thresholds
    # If there are more than two modes, keep the two extremes and leave the middle to "normal/medium/transitory" state
    thres <- c(first(model$parameters$mean), last(model$parameters$mean))
    
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
    test$classification[test[[x]] > thres[1] & test[[x]] < thres[2]] <- "Normal"
    test$classification[test[[x]] >= thres[2]] <- "High"
    test$classification <- factor(test$classification, levels = c("High","Normal","Low"))
    
    # The palette with black:
    # cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    
    ggplot(data = test, aes(x = Year, y = !!sym(x), group = 1))+
      geom_path() +
      geom_point(aes(color = classification)) +
      scale_color_manual(values = c("#009E73","#56B4E9","#D55E00")) +
      facet_wrap(.~Sample) +
      geom_hline(yintercept = thres[1], linetype = "dashed", linewidth = 0.1) +
      geom_hline(yintercept = thres[2], linetype = "dashed", linewidth = 0.1) +
      theme_classic()
    
    
  }
  
  
  
# TRASH ----
  
  
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
  
  




