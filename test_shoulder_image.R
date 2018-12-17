# Test Shoulder Images

# Origin: Landsat vs Ground Truth Data Exploration

# 20181002
# Ryan Nagelkirk

# Description:
# 1. Plots the r-squared values of all bands from min, max and shoulder images from all parks to see 
#    which has the better relationship across bands 

########  Reminder #################################################
# This code doesn't attempt to merge point values from overlapping imagery. Instead,
# each image is treated independently and has it's own row in the final.df
# If you want to handle those points that do overlap, that code still has to be written


#######  Presentations with results from this code (not all ppts are listed)  ######################
# E:\Dropbox\Permanent\Grad School\Meetings\Kyla\20181130_Testing Transition Image n Search for consistent regression\20181130_Testing Transition Image & Search for consistent regression

# Libraries
library(raster)
library(rgdal)
library(stats)


###########################################################################################################
#  Variables and Directories 
###########################################################################################################
# Set working directory to ParkData folder on dropbox
# setwd("E:/Dropbox/Permanent/Grad School/Projects/EleTree/data/ParkData")
setwd("~/Dropbox/Permanent/Grad School/Projects/EleTree/data/ParkData")

# Set the folder for all the CSVs
csv.folder <- "E:/Dropbox/Permanent/Grad School/Projects/EleTree/analysis/Landsat_v_woody_cover_regressions/"

# Source the functions
# source("E:/Dropbox/Permanent/Grad School/src_functions/src_masterfunctions.R")
source("~/Dropbox/Permanent/Grad School/src_functions/src_masterfunctions.R")

# final.df.w.SpatialContext <- "E:/Dropbox/Permanent/Grad School/Projects/EleTree/analysis/Landsat_v_woody_cover_regressions/TwelvePark_Landsat_n_GT_values_w_SpatContext.csv"
final.df.w.SpatialContext <- "~/Dropbox/Permanent/Grad School/Projects/EleTree/analysis/Landsat_v_woody_cover_regressions/TwelvePark_Landsat_n_GT_values_w_SpatContext.csv"

# final.df.tran <- "E:/Dropbox/Permanent/Grad School/Projects/EleTree/analysis/Landsat_v_woody_cover_regressions/AllTwelveParks_Landsat_n_GT_values_TRANS.csv"
final.df.tran <- "~/Dropbox/Permanent/Grad School/Projects/EleTree/analysis/Landsat_v_woody_cover_regressions/AllTwelveParks_Landsat_n_GT_values_TRANS.csv"

# List the parks
park.names <- c("Chobe", "Kruger", "Limpopo", "Mpala", "Murchison", "North_Luangwa", "QWE", "Ruaha", "Selous", "Serengeti", "South_Luangwa")

# Identify the folder with the Landsat images
landsat.folder <- "X:/nagelki4/Projects/EleTree/data/Landsat/20km_Buffer"

# List only the images that are max or min, but not the minmax 
# (that's actually how I exported them from GEE, so we're good to go by just specifying "reflectance")
landsat.rasters <- list.files(landsat.folder, pattern = "minmax_NDVI_wo_indices_allBands.tif", full.names = T)



#########################################################################################################################
# Read in the files with the Landsat band values for the images (min, max and shoulder images)
#########################################################################################################################

# Read in the file with spatial context
final.df <- read.csv(final.df.w.SpatialContext)
# Read in the file with the Landsat band values from the shoulder image
tran.final.df <- read.csv(final.df.tran)

# The tran df has more rows. I looked at the rows carefully, and my guess is that the tran images overlapped
# more than the min and max images, so more points fell in both paths than did in the min and max. Though
# the first place to have a dropped number was Mpala, so some of it has to do with something else well. 
# I decided that because this is only a test, getting perfect matching between all the points was not necessary
# (I don't expect the difference in point counts to affect the r-squared significantly)


# Combine the two df into a list to be read through the following loop
df.list <- list(final.df, tran.final.df)




###########################################################################################################
###  Compute r-squared values
###########################################################################################################

# Loop through the df.list, computing the r2 as you go
for(r in 1:length(df.list)){
  final.df <- as.data.frame(df.list[r])

  # Create an empty table with the same column names and add PCA columns
  gtbl <- final.df
  gtbl <- gtbl[, -c(2:7)] # remove columns
  gtbl <- gtbl[-c(2:NROW(gtbl)), ] # remove all but first row
  
  # Make a list of the column names that will be regressed
  colNamesList <- colnames(gtbl)[-1] # remove "park" from the column list
  
  # Create an empty row for reuse in the loops below
  empty.row <- gtbl
  
  # Loop through the parks, plus one for the All Park analysis
  for(k in 1:(length(park.names)+1)){ # Add one for the inter-park analysis
    
    ### Get the values wanted for the regression 
    # Set the park name 
    if(k <= length(park.names)){
      park.name <- park.names[k]
    }else{
      park.name <- "AllParks"
    }
    
    # Subset to the particular park or all parks
    if(park.name == "AllParks"){
      f.df <- final.df
    }else{
      f.df <- final.df[final.df$park == park.name, ]
    }
    
    
    ####  Loop through the simple ones (the reflectance bands)  ####################################
    # Get the xvar (w.c.) first 
    xvar <- f.df$hand.tree
    
    # Get all the r2 values of the simple variables
    for(i in 1:length(colNamesList)){
      # Get the adjusted r2
      empty.row[, i+1] <-  r2(xvar, f.df[, colNamesList[i]])
    }
    
    # Add the park name as the front of the list
    empty.row$park <- park.name

    # Rbind (unless it's the first)
    if(k == 1){
      gtbl <- empty.row
    }else{
      gtbl <- rbind(gtbl, empty.row)
    }
  }
  
  # Move the gtbl to a tran variable if necessary (this is included so there 
  # are two gtbl versions for tran testing below)
  if(r == 1){
    gtbl.minmax <- gtbl
  }else{
    gtbl.tran <- gtbl
  }
}

# Move the minmax tbl back to gtbl for the section of code beyond the next chunk to still work
# (The sections that don't have to do with the shoulder image)
gtbl <- gtbl.minmax

#########################################################################################################
#######  Plot trans, min and max r-squared values 
#########################################################################################################

# Create individual min max tables for cleaner code later
gtbl.min <- gtbl.minmax[, c(1:10)]
gtbl.max <- gtbl.minmax[, c(1, 11:19)]

# Then do an eval, populating a new table with the name of the best image for each band
# In other words, the image that has the highest r-squared for each band/park combo. 
eval.tbl <- gtbl.tran
names(eval.tbl) <- c("park", "band1", "band2", "band3", "band4", "band5", "band6", "band7", "band8", "band9")

# For each row and column, eval and paste in the name of the image that had the highest r-squared value
for(i in 1:(ncol(eval.tbl)-1)){ # Had to take one off for the park column
  for(j in 1:nrow(eval.tbl)){
    # Create a list of the values
    val.list <- c(gtbl.min[j, i+1], gtbl.tran[j, i+1], gtbl.max[j, i+1])
    # Create list of image names
    img.list <- c("min", "tran", "max")
    # Pick the image name that corresponds to the largest number in the list
    img <- img.list[which(val.list == max(val.list))]
    
    # Rename it, inserting max value and the difference between the best r-squared value and the runner up
    marg <- round(max(val.list) - median(val.list), 3)
    img.renamed <- paste0(img," | ", marg, " | ", max(val.list))
    # Assign the name 
    eval.tbl[j, i+1] <- img.renamed
  }
}

##  Create a plotting function to show the r-squared values
# Combine the tables
gtbl.master <- cbind(gtbl.min, gtbl.tran, gtbl.max)
# Get rid of the two park columns that aren't needed
gtbl.master <- gtbl.master[, -c(11, 21)]

rSQRplotter <- function(parkname){
  # Get the max and have all the plots on the same scale
  ymax <- round(max(gtbl.master[, c(2:27)]), 1)
  
  # Get just one park
  gtbl.park <- gtbl.master[which(gtbl.master$park == parkname), ]
  gtbl.park <- t(gtbl.park[, names(gtbl.park) != "park"]) # Take of the park column 
  
  plot(gtbl.park[grep("min", row.names(gtbl.park))], pch = "M", col = "red", main = parkname, ylab = "r-squared", xlab = "Band", ylim = c(0, ymax))
  points(gtbl.park[grep("tran", row.names(gtbl.park))], pch = "T")
  points(gtbl.park[grep("max", row.names(gtbl.park))], pch = "X", col = "blue")
}


# Loop through and produce the plots
for(park in park.names){
  rSQRplotter(park)
}


