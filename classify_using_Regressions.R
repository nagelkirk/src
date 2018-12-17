# Classify Landsat Imagery Using Regressions 

# Origin: Landsat vs Ground Truth Data Exploration

# 20181002
# Ryan Nagelkirk

# Description:
# 1. Classify woody cover using linear regressions, PCAs and forward stepwise regression 
# 2. Output the results, with the 70 training sites reported as NA, and the testing sites having the
#    classification results 

# Specifically do all the bands and indices, along with the PCA results in a multiple regression

########  Reminder ###################################################################
# This code doesn't attempt to merge point values from overlapping imagery. Instead,
# each image is treated independently and has it's own row in the final.df
# If you want to handle those points that do overlap, that code still has to be written

#######  Presentations with results from this code (not all ppts are listed)  ######################
# E:\Dropbox\Permanent\Grad School\Meetings\Kyla\20181130_Testing Transition Image n Search for consistent regression\20181130_Testing Transition Image & Search for consistent regression

# Libraries
library(raster)
library(rgdal)
library(stats)
library(psych) # for PCA
library(nFactors) # for scree plot
library(QuantPsyc) # for lm.beta



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

# List the parks
park.names <- c("Chobe", "Kruger", "Limpopo", "Mpala", "Murchison", "North_Luangwa", "QWE", "Ruaha", "Selous", "Serengeti", "South_Luangwa")


#########################################################################################################################
# Read in the files with the Landsat, Spatial Context and Reference Data/Ground Truth values 
#########################################################################################################################
# Read in the file with spatial context
final.df <- read.csv(final.df.w.SpatialContext)



#########################################################################################################################
# Band Regressions 
#########################################################################################################################

# Add all the derived variables to the final.df. This way, all the correlations (except the multiband)
# can be ran in a loop
{
  final.df$band1Diff <- final.df$min1 - final.df$max1
  final.df$band2Diff <- final.df$min2 - final.df$max2
  final.df$band3Diff <- final.df$min3 - final.df$max3
  final.df$band4Diff <- final.df$min4 - final.df$max4
  final.df$band5Diff <- final.df$min5 - final.df$max5
  final.df$band6Diff <- final.df$min6 - final.df$max6
  final.df$band7Diff <- final.df$min7 - final.df$max7
  final.df$band8Diff <- final.df$min8 - final.df$max8
  final.df$band9Diff <- final.df$min9 - final.df$max9
  
  final.df$band1Sum <- final.df$min1 + final.df$max1
  final.df$band2Sum <- final.df$min2 + final.df$max2
  final.df$band3Sum <- final.df$min3 + final.df$max3
  final.df$band4Sum <- final.df$min4 + final.df$max4
  final.df$band5Sum <- final.df$min5 + final.df$max5
  final.df$band6Sum <- final.df$min6 + final.df$max6
  final.df$band7Sum <- final.df$min7 + final.df$max7
  final.df$band8Sum <- final.df$min8 + final.df$max8
  final.df$band9Sum <- final.df$min9 + final.df$max9
  
  # Also have the difference between key bands: red and NIR, NIR and SWIR2 (soil is high in SWIR and veg is low)
  final.df$min_NIR_minus_R <- final.df$min5 - final.df$min4
  final.df$max_NIR_minus_R <- final.df$max5 - final.df$max4
  final.df$min_NIR_minus_SWIR2 <- final.df$min5 - final.df$min7
  final.df$max_NIR_minus_SWIR2 <- final.df$max5 - final.df$max7
  final.df$NIRdiff_minus_Rdiff <- (final.df$max5 - final.df$min5)  - (final.df$min4 - final.df$max4) # This is just a wild card that might yield something. Computed so that each subtraction should yield a positive value.
  
  # Calc the index values
  final.df$NDVI_min <- (final.df$min5 - final.df$min4) / (final.df$min5 + final.df$min4)
  final.df$NDVI_max <- (final.df$max5 - final.df$max4) / (final.df$max5 + final.df$max4)
  final.df$NDVIdiff <- final.df$NDVI_max - final.df$NDVI_min
  
  # MSAVI2 needs reflectance values to be from 0-1, so redefine the bands
  minNIR <- final.df$min5 / 10000
  maxNIR <- final.df$max5 / 10000
  minRed <- final.df$min4 / 10000
  maxRed <- final.df$max4 / 10000
  minSWIR2 <- final.df$min7 / 10000
  maxSWIR2 <- final.df$max7 / 10000
  minSWIR1 <- final.df$min6 / 10000
  maxSWIR1 <- final.df$max6 / 10000
  
  # MSAVI2
  final.df$MSAVI2_min <- (2*minNIR + 1 - (((2*minNIR + 1)^2) - (8*(minNIR - minRed)))^.5) / 2
  final.df$MSAVI2_max <- (2*maxNIR + 1 - (((2*maxNIR + 1)^2) - (8*(maxNIR - maxRed)))^.5) / 2
  
  # SATVI (I set L in the equation to 0.5 based on Qi's MSAVI2 paper saying that's what the MSAVI developer used,
  # and the ppl that made SATVI reference the L used in MSAVI, so I assume they also use 0.5)
  final.df$SATVI_min <- (((minNIR - minRed) / (minNIR + minRed + 0.5)) * 1.5) - (minSWIR2 / 2)
  final.df$SATVI_max <- (((maxNIR - maxRed) / (maxNIR + maxRed + 0.5)) * 1.5) - (maxSWIR2 / 2)
  
  # SNDI 
  final.df$SNDI_min <- (minRed - (minNIR + minSWIR1)) / (minRed + (minNIR + minSWIR1))
  final.df$SNDI_max <- (maxRed - (maxNIR + maxSWIR1)) / (maxRed + (maxNIR + maxSWIR1))
}


  
##################################################################################################
#####  Classify WC using Regressions #############################################################
##################################################################################################
# Approach: Use 70% of the points for training and 30% for testing.


# Create an empty table with the same column names as the final.df and add PCA columns
VE.tbl <- final.df
VE.tbl <- VE.tbl[, -c(2:7)] # remove columns
VE.tbl <- VE.tbl[-c(2:NROW(VE.tbl)), ] # remove all but first row

# Add all of the other columns here
{
  
  # Add PCA columns
  VE.tbl$minAll3CompsWC <- 0
  VE.tbl$min_Comp1_var <- NA # these will be actual names of bands 
  VE.tbl$min_Comp2_var <- NA
  VE.tbl$min_Comp3_var <- NA
  
  VE.tbl$maxAll3Comps <- 0
  VE.tbl$max_Comp1_var <- NA
  VE.tbl$max_Comp2_var <- NA
  VE.tbl$max_Comp3_var <- NA
  
  
  VE.tbl$minmaxComp1 <- 0
  VE.tbl$minmaxComp2 <- 0
  VE.tbl$minmaxComp3 <- 0
  VE.tbl$minmaxAll3Comps <- 0
  VE.tbl$minmax_Comp1_var <- NA
  VE.tbl$minmax_Comp2_var <- NA
  VE.tbl$minmax_Comp3_var <- NA
  VE.tbl$minmax_Comp1_var_VE <- 0
  VE.tbl$minmax_Comp2_var_VE <- 0
  VE.tbl$minmax_Comp3_var_VE <- 0
  VE.tbl$minmax_All3Comps_var_VE <- 0
}

# Create an empty row for reuse in the loops below
empty.row <- VE.tbl


####################################################################################
###  Create the models and compute the VEcv values
####################################################################################

# First, write a function that will create the model and calc the VEcv
model_VEcv <- function(xval, yval){
  
  # Put the values together as a new table
  df <- cbind(xval, yval)
  df <- as.data.frame(df)
  
  # Create the training and testing data
  seventy_p_rows <- round(nrow(df) * 0.7) # this will return the number of rows == 70% of the rows
  
  # Get a random draw of 70
  random.draw <- sample(1:NROW(df), seventy_p_rows, replace = F)
  
  # Subset to the rows wanted for training and testing datasets 
  df.train <- df[random.draw, ]
  df.test <- df[-random.draw, ]
  
  # Create the lm
  mod <- lm(xval ~ yval, data = df.train)
  
  # Subset to just the data you want (the band)
  test_band <- subset(df.test, select = c(yval))
  
  # Predict WC using the model and the testing data
  wc.predicted <- predict(mod, test_band)
  
  # Return the VEcv of the predictions 
  return(round(VEcv(df.test$xval, wc.predicted), 3))
  
}



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
  
  for(i in 1:length(colNamesList)){
    # Get the name of the column, then the column number that corresponds with 
    cName <- colNamesList[i]
    cNumber <- which(colnames(f.df) == cName)
    # Find the complete cases (some of the radius values will be NA)
    comp.cases <- complete.cases(f.df[, cNumber])
    
    # Run the model and calc VEcv function 
    empty.row[, i+1] <- model_VEcv(xval = f.df$hand.tree[comp.cases], yval = f.df[comp.cases, cNumber]) 
  }
  
  
  {
    # Create the training and testing dfs
    seventy_p_rows <- round(nrow(f.df) * 0.7) # this will return the number of rows == 70% of the rows
    
    # Get a random draw of 70
    random.draw <- sample(1:NROW(f.df), seventy_p_rows, replace = F)
    
    # Subset to the rows wanted for training and testing datasets 
    df.train <- f.df[random.draw, ]
    df.test <- f.df[-random.draw, ]
    
    ####  Do the band differences etc. one at a time
    # Directly below are the 3 band regressions
    md <- lm(hand.tree ~ min3 + min4 + min8, data = df.train)
    # Subset to just the data you want (the band)
    test_band <- subset(df.test, select = c(min3, min4, min8))
    # Predict WC using the model and the testing data
    wc.predicted <- predict(md, test_band)
    # Return the VEcv of the predictions 
    empty.row$min348 <- VEcv(df.test$hand.tree, wc.predicted)
    
    # Do the next (max)
    md <- lm(hand.tree ~ max3 + max4 + max8, data = df.train)
    # Subset to just the data you want (the band)
    test_band <- subset(df.test, select = c(max3, max4, max8))
    # Predict WC using the model and the testing data
    wc.predicted <- predict(md, test_band)
    # Return the VEcv of the predictions 
    empty.row$max348 <- VEcv(df.test$hand.tree, wc.predicted)
    
    # Do the next (minmax)
    md <- lm(hand.tree ~ min3 + max3, data = df.train)
    # Subset to just the data you want (the band)
    test_band <- subset(df.test, select = c(min3, max3))
    # Predict WC using the model and the testing data
    wc.predicted <- predict(md, test_band)
    # Return the VEcv of the predictions 
    empty.row$minmax33 <- VEcv(df.test$hand.tree, wc.predicted)
    
  }
  
  
  ######  PCA  ###########################################################################################
  library(psych) # for PCA
  library(nFactors) # for scree plot
  library(QuantPsyc) # for lm.beta
  
  
  # Function to model woody cover using the PCA results and report the VEcv of the model
  PCA_VEcv <- function(df, xval){
    # Rotate with 3 dim extracted
    # Run the rotated PCA 
    fit.rot <- principal(df, nfactors = 3, rotate = "varimax") # c(-13:-20)
    
    # Put the scores into the df as a new column 
    df$Comp1 <- fit.rot$scores[, 1]
    df$Comp2 <- fit.rot$scores[, 2]
    df$Comp3 <- fit.rot$scores[, 3]
    
    # Model woody cover using the components and get the VEcv
    vals <- c(model_VEcv(xval, df$Comp1), model_VEcv(xval, df$Comp2), model_VEcv(xval, df$Comp3))
    
    
    # Do all three
    # Create the training and testing dfs
    seventy_p_rows <- round(nrow(df) * 0.7) # this will return the number of rows == 70% of the rows
    
    # Get a random draw of 70
    random.draw <- sample(1:NROW(df), seventy_p_rows, replace = F)
    
    # Subset to the rows wanted for training and testing datasets 
    df.train <- df[random.draw, ]
    df.test <- df[-random.draw, ]
    x.train <- xval[random.draw]
    x.test <- xval[-random.draw]
    
    # Combine into training and testing dfs
    df.train <- cbind(df.train, x.train)
    df.test <- cbind(df.test, x.test)
    
    # Create the model
    md <- lm(x.train ~ Comp1 + Comp2 + Comp3, data = df.train) 
    
    # Subset to just the data you want for the predictions 
    test_band <- subset(df.test, select = c(Comp1, Comp2, Comp3))
    
    # Predict WC using the model and the testing data
    wc.predicted <- predict(md, test_band)
    
    # Add the model VEcv to the list 
    vals <- c(vals, round(VEcv(df.test$x.test, wc.predicted), 3)) # x.test is woody cover (hand.tree)
    
    # Now get the highest loaded bands 
    high1 <- names(which(abs(fit.rot$loadings[, 1]) == max(abs(fit.rot$loadings[, 1]))))
    high2 <- names(which(abs(fit.rot$loadings[, 3]) == max(abs(fit.rot$loadings[, 3])))) # the numbering is right - for some reason the columns are switched btwn RC2 and 3
    high3 <- names(which(abs(fit.rot$loadings[, 2]) == max(abs(fit.rot$loadings[, 2]))))
    
    # Put those names in the val list
    vals <- c(vals, high1, high2, high3)
    
    # Use the band names to get the respective column's data and get the VEcv of the models
    high1_vals <- df[, which(names(df) == high1)]
    high2_vals <- df[, which(names(df) == high2)]
    high3_vals <- df[, which(names(df) == high3)]
    
    # VEcv
    vals <- c(vals, model_VEcv(xval, high1_vals), model_VEcv(xval, high2_vals), model_VEcv(xval, high3_vals))
    
    
    # Subset to the values
    threeVar.train <- subset(df.train, select = c(high1, high2, high3, 'x.train'))
    threeVar.test <- subset(df.test, select = c(high1, high2, high3, 'x.test'))
    
    # Change the names
    names(threeVar.train) <- c("b1", "b2", "b3", "x.train")
    names(threeVar.test) <- c("b1", "b2", "b3", "x.test")
    
    # Do all three
    md <- lm(x.train ~ b1 + b2 + b3, data = threeVar.train)
    
    # Subset to just the data you want for the predictions 
    test_band <- subset(threeVar.test, select = c(b1, b2, b3))
    
    # Predict WC using the model and the testing data
    wc.predicted <- predict(md, test_band)
    
    # Add the model VEcv to the list 
    vals <- c(vals, round(VEcv(threeVar.test$x.test, wc.predicted), 3))
    
    return(vals)
  }
  
  
  
  {
    # Enter the min image bands to get PCA regressions etc.
    minPCAresults <- PCA_VEcv(f.df[, c(8:16)], f.df$hand.tree)
    maxPCAresults <- PCA_VEcv(f.df[, c(17:25)], f.df$hand.tree)
    minmaxPCAresults <- PCA_VEcv(f.df[, c(8:25)], f.df$hand.tree)
    # Break the spatial context up into measures across images
    desired.cols.min <- which(substr(colnames(f.df), 1, 3) == "min" & nchar(colnames(f.df)) > 4)
    desired.cols.max <- which(substr(colnames(f.df), 1, 3) == "max" & nchar(colnames(f.df)) > 4)
    min.tbl <- f.df[, desired.cols.min]
    max.tbl <- f.df[, desired.cols.max]
    # Use only the rows with complete cases
    spatialContextPCAresults.min <- PCA_VEcv(min.tbl[complete.cases(min.tbl), ] , f.df$hand.tree[complete.cases(min.tbl)])   
    spatialContextPCAresults.max <- PCA_VEcv(max.tbl[complete.cases(max.tbl), ] , f.df$hand.tree[complete.cases(max.tbl)]) 
    
    # Plug the values into the empty row and bind it to the growing df
    # Min NDVI image
    empty.row$minComp1 <- as.numeric(minPCAresults[1])
    empty.row$minComp2 <- as.numeric(minPCAresults[2])
    empty.row$minComp3 <- as.numeric(minPCAresults[3])
    empty.row$minAll3Comps <- as.numeric(minPCAresults[4])
    empty.row$min_Comp1_var <- minPCAresults[5]
    empty.row$min_Comp2_var <- minPCAresults[6]
    empty.row$min_Comp3_var <- minPCAresults[7]
    empty.row$min_Comp1_var_VE <- as.numeric(minPCAresults[8])
    empty.row$min_Comp2_var_VE <- as.numeric(minPCAresults[9])
    empty.row$min_Comp3_var_VE <- as.numeric(minPCAresults[10])
    empty.row$min_All3Comps_var_VE <- as.numeric(minPCAresults[11])
    
    # Max NDVI image
    empty.row$maxComp1 <- as.numeric(maxPCAresults[1])
    empty.row$maxComp2 <- as.numeric(maxPCAresults[2])
    empty.row$maxComp3 <- as.numeric(maxPCAresults[3])
    empty.row$maxAll3Comps <- as.numeric(maxPCAresults[4])
    empty.row$max_Comp1_var <- maxPCAresults[5]
    empty.row$max_Comp2_var <- maxPCAresults[6]
    empty.row$max_Comp3_var <- maxPCAresults[7]
    empty.row$max_Comp1_var_VE <- as.numeric(maxPCAresults[8])
    empty.row$max_Comp2_var_VE <- as.numeric(maxPCAresults[9])
    empty.row$max_Comp3_var_VE <- as.numeric(maxPCAresults[10])
    empty.row$max_All3Comps_var_VE <- as.numeric(maxPCAresults[11])
    
    # Min Max image
    empty.row$minmaxComp1 <- as.numeric(minmaxPCAresults[1])
    empty.row$minmaxComp2 <- as.numeric(minmaxPCAresults[2])
    empty.row$minmaxComp3 <- as.numeric(minmaxPCAresults[3])
    empty.row$minmaxAll3Comps <- as.numeric(minmaxPCAresults[4])
    empty.row$minmax_Comp1_var <- minmaxPCAresults[5]
    empty.row$minmax_Comp2_var <- minmaxPCAresults[6]
    empty.row$minmax_Comp3_var <- minmaxPCAresults[7]
    empty.row$minmax_Comp1_var_VE <- as.numeric(minmaxPCAresults[8])
    empty.row$minmax_Comp2_var_VE <- as.numeric(minmaxPCAresults[9])
    empty.row$minmax_Comp3_var_VE <- as.numeric(minmaxPCAresults[10])
    empty.row$minmax_All3Comps_var_VE <- as.numeric(minmaxPCAresults[11])
    
    # Spatial Context PCA data: min image
    empty.row$min_SCcomp1 <- as.numeric(spatialContextPCAresults.min[1])
    empty.row$min_SCcomp2 <- as.numeric(spatialContextPCAresults.min[2])
    empty.row$min_SCcomp3 <- as.numeric(spatialContextPCAresults.min[3])
    empty.row$min_SCall3Comps <- as.numeric(spatialContextPCAresults.min[4])
    empty.row$min_SC_Comp1_var <- spatialContextPCAresults.min[5]
    empty.row$min_SC_Comp2_var <- spatialContextPCAresults.min[6]
    empty.row$min_SC_Comp3_var <- spatialContextPCAresults.min[7]
    empty.row$min_SC_Comp1_var_VE <- as.numeric(spatialContextPCAresults.min[8])
    empty.row$min_SC_Comp2_var_VE <- as.numeric(spatialContextPCAresults.min[9])
    empty.row$min_SC_Comp3_var_VE <- as.numeric(spatialContextPCAresults.min[10])
    empty.row$min_SC_All3Comps_var_VE <- as.numeric(spatialContextPCAresults.min[11])
    
    # Spatial Context PCA data: max image
    empty.row$max_SCcomp1 <- as.numeric(spatialContextPCAresults.max[1])
    empty.row$max_SCcomp2 <- as.numeric(spatialContextPCAresults.max[2])
    empty.row$max_SCcomp3 <- as.numeric(spatialContextPCAresults.max[3])
    empty.row$max_SCall3Comps <- as.numeric(spatialContextPCAresults.max[4])
    empty.row$max_SC_Comp1_var <- spatialContextPCAresults.max[5]
    empty.row$max_SC_Comp2_var <- spatialContextPCAresults.max[6]
    empty.row$max_SC_Comp3_var <- spatialContextPCAresults.max[7]
    empty.row$max_SC_Comp1_var_VE <- as.numeric(spatialContextPCAresults.max[8])
    empty.row$max_SC_Comp2_var_VE <- as.numeric(spatialContextPCAresults.max[9])
    empty.row$max_SC_Comp3_var_VE <- as.numeric(spatialContextPCAresults.max[10])
    empty.row$max_SC_All3Comps_var_VE <- as.numeric(spatialContextPCAresults.max[11])
  }
  
  # Add the park name as the front of the list
  empty.row$park <- park.name
  
  # Rbind (unless it's the first)
  if(k == 1){
    VE.tbl <- empty.row
  }else{
    VE.tbl <- rbind(VE.tbl, empty.row)
  }
  
}

# Print the table
# write.csv(VE.tbl, "~/Dropbox/Permanent/Grad School/Projects/EleTree/analysis/Landsat_v_woody_cover_regressions/Regression_Models_VEcv_ALLpark_SpatContxt.csv")
write.csv(VE.tbl, "E:/Dropbox/Permanent/Grad School/Projects/EleTree/analysis/Landsat_v_woody_cover_regressions/Regression_Models_VEcv_ALLpark_SpatContxt.csv")





#########################################################################################################
####  Forward Stepwise Regression  
#########################################################################################################

# I did the below for all the parks together. All of this could be written as a function
# that would automatically find the bands to include, but having two versus one band
# only made the regression slightly better, so I think it would be fair to leave out
# the stepwise regression from the paper. Unless I can come up with a theory for why
# the bands selected here, green (3) and NIR (5), interacted to be the best two variables
# The r-squared improved only from 0.16 to 0.17 when going from one band in the regression
# to two bands. 

# Load necessary libraries
library(HH)

## Build the table for no variables in the equation
## This will be an lm of each and pulling the t and p values
var.tab.1 <- matrix(ncol = 4, nrow = ncol(f.df)-7)
var.tab.1 <- as.data.frame(var.tab.1)
names(var.tab.1) <- c("Variable", "Tolerance", "t-value", "p-value")

# Fill in the variable names and tolerances, which are 1 for the first go
var.tab.1$Variable <- c(names(f.df)[8:ncol(f.df)])
var.tab.1$Tolerance <- 1

# Loop through and put in the correlation, tvalue and pvalue
i <- 1
for(i in 1:nrow(var.tab.1)){
  var <- which(names(f.df) == var.tab.1$Variable[i])
  mod <- lm(f.df$hand.tree ~ f.df[, var])
  mod.vals <- summary(mod)
  tval <- mod.vals$coefficients[2, 3]
  pval <- mod.vals$coefficients[2, 4]
  
  # Plug into table
  var.tab.1$`t-value`[i] <- tval
  var.tab.1$`p-value`[i] <- pval/2
}


##### Now add the most significant variable (max3)
{
  # First, look at the model
  first.mod <- lm(f.df$hand.tree ~ f.df$max3)
  first.sum <- summary(first.mod)
  
  # Make a table
  results.1 <- matrix(ncol = 6, nrow = 1)
  results.1 <- as.data.frame(results.1)
  names(results.1) <- c("Variable", "Coefficient", "Standard.Error", "Tolerance", "t-value", "p-value")
  
  # Plug in the values
  tval <- first.sum$coefficients[2, 3]
  pval <- first.sum$coefficients[2, 4]
  se <- first.sum$coefficients[2, 2]
  coeff <- first.sum$coefficients[2, 1]
  
  # Plug into table
  results.1$`t-value` <- tval
  results.1$`p-value` <- pval/2
  results.1$Tolerance <- 1
  results.1$Variable <- "max3"
  results.1$Coefficient <- coeff
  results.1$Standard.Error <- se
  
  
  #### Evaluate the new adds
  var.tab.2 <- matrix(ncol = 4, nrow = ncol(f.df)-8)
  var.tab.2 <- as.data.frame(var.tab.2)
  names(var.tab.2) <- c("Variable", "Tolerance", "t-value", "p-value")
  
  # Fill in the variable names and tolerances, which are 1 for the first go
  r.names <- c(names(f.df)[8:ncol(f.df)])
  L4.num <- which(r.names == "max3")
  var.tab.2$Variable <- r.names[-L4.num]
}

# Loop through and put in the correlation, tvalue and pvalue
i <- 1
for(i in 1:nrow(var.tab.2)){
  var <- which(names(f.df) == var.tab.2$Variable[i])
  
  # Make table 
  temp <- matrix(ncol = 3, nrow = nrow(f.df))
  temp <- as.data.frame(temp)
  temp$V1 <- f.df$hand.tree
  temp$V2 <- f.df$max3
  temp$V3 <- f.df[[var]]
  
  
  mod <- lm(temp$V1 ~ temp$V2 + temp$V3)
  mod.vals <- summary(mod)
  tval <- mod.vals$coefficients[3, 3]
  pval <- mod.vals$coefficients[3, 4]
  
  # Plug into table
  var.tab.2$`t-value`[i] <- tval
  var.tab.2$`p-value`[i] <- pval/2
  var.tab.2$Tolerance[i] <- 1/HH::vif(mod)[1]
}



### Add max5 and run again! (It had a really high tolerance, which made it better than some of the more signifcant
# variables, like max6 
# First, look at the model

first.mod <- lm(f.df$hand.tree ~ f.df$max3 + f.df$max5)
first.sum <- summary(first.mod)

# Make a table
results.1 <- matrix(ncol = 6, nrow = 2)
results.1 <- as.data.frame(results.1)
names(results.1) <- c("Variable", "Coefficient", "Standard.Error", "Tolerance", "t-value", "p-value")

# Plug in the values
for(t in 1:nrow(results.1)){
  var.list <- c("max3", "max5")
  tval <- first.sum$coefficients[t+1, 3]
  pval <- first.sum$coefficients[t+1, 4]
  se <- first.sum$coefficients[t+1, 2]
  coeff <- first.sum$coefficients[t+1, 1]
  
  # Plug into table
  results.1$`t-value`[t] <- tval
  results.1$`p-value`[t] <- pval/2
  results.1$Tolerance[t] <- 1/HH::vif(first.mod)[t]
  results.1$Variable[t] <- var.list[t]
  results.1$Coefficient[t] <- coeff
  results.1$Standard.Error[t] <- se
}


#### Evaluate the new adds
var.tab.3 <- matrix(ncol = 4, nrow = ncol(f.df)-9)
var.tab.3 <- as.data.frame(var.tab.3)
names(var.tab.3) <- c("Variable", "Tolerance", "t-value", "p-value")

# Fill in the variable names and tolerances
r.names <- c(names(f.df)[8:ncol(f.df)])
L4.num <- which(r.names == "max3")
L6.num <- which(r.names == "max5")
var.tab.3$Variable <- r.names[-c(L4.num, L6.num)]


# Loop through and put in the correlation, tvalue and pvalue
i <- 1
for(i in 1:nrow(var.tab.3)){
  var <- which(names(f.df) == var.tab.3$Variable[i])
  
  # Make table to use
  temp <- matrix(ncol = 4, nrow = nrow(f.df))
  temp <- as.data.frame(temp)
  temp$V1 <- f.df$hand.tree
  temp$V2 <- f.df$max3
  temp$V3 <- f.df$max5
  temp$V4 <- f.df[[var]]
  
  mod <- lm(temp$V1 ~ temp$V2 + temp$V3 + temp$V4)
  mod.vals <- summary(mod)
  tval <- mod.vals$coefficients[3, 3]
  pval <- mod.vals$coefficients[3, 4]
  
  # Plug into table
  var.tab.3$`t-value`[i] <- tval
  var.tab.3$`p-value`[i] <- pval/2
  var.tab.3$Tolerance[i] <- 1/HH::vif(mod)[1]
}



## Discussion
mod <- lm(f.df$hand.tree ~ f.df$max3 + f.df$max5)
# plot(mod)
lm.beta(mod)

# Plot against the two variables
plot(f.df$max3, resid(mod), main = "Residuals vs. Landsat.4")
abline(0,0)

plot(f.df$max6, resid(mod), main = "Residuals vs. L6.norm")
abline(0,0)


