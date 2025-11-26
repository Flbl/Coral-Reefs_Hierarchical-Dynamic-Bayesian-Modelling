#############################################################################################
#                                                                                           #
#                                                                                           #
########################      Booting all pipeline prerequisites          ###################
#                                                                                           #
#                                                                                           #
#############################################################################################



# Common paths

pathDat <- file.path("Data")
pathSpe <- file.path("Data", "Species")
pathEnv <- file.path("Data", "Environment")
pathPro <- file.path("Data", "Processed")
pathRes <- file.path("Results")

# Create dirs
dir.create("Data", showWarnings = FALSE)
dir.create(file.path("Data","Processed"), showWarnings = FALSE)
dir.create(file.path("Data","Environment"), showWarnings = FALSE)
dir.create(file.path("Data","Species"), showWarnings = FALSE)


# Load packages ----

# Common Packages

# if (!require("reshape2")) install.packages("reshape2")
# if (!require("plyr")) install.packages("plyr")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("ggnewscale")) install.packages("ggnewscale")
# if (!require("ggrepel")) install.packages("ggrepel")
# if (!require("ggsignif")) install.packages("ggsignif")


# Specific packages for specific scripts

# 01_State_data_conversion 

if(init1 == "01_Data_To_States.R") {
  
  
  # Travel time
  # This package will install pretty much all other needed packages for SIG processing
  if (!require("brms")) install.packages("brms")
  if (!require("mclust")) install.packages("mclust")
  if (!require("purrr")) install.packages("purrr")
  # if (!require("Kendall")) install.packages("Kendall")
  # if (!require("traveltime")) install.packages("traveltime", repos = c("https://idem-lab.r-universe.dev"))
  # if (!require("mapsapi")) install.packages("mapsapi")
  # if (!require("rSDM")) install.packages("rSDM", repos = c("https://pakillo.r-universe.dev", "https://cloud.r-project.org"))
  # if (!require("concaveman")) install.packages("concaveman")
  # if (!require("nngeo")) install.packages("nngeo")
  
  
  
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


