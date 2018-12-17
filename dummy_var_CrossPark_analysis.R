# Regression Dummy Variable Analysis

# Origin: Landsat vs Ground Truth Data Exploration

# 20181002
# Ryan Nagelkirk

# Description: 
# Does a dummy variable analysis to determine whether relationships between landsat bands and woody cover are
# consistent across parks



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

# Create an empty table with the same column names and add PCA columns
gtbl <- final.df
gtbl <- gtbl[, -c(2:7)] # remove columns
gtbl <- gtbl[-c(2:NROW(gtbl)), ] # remove all but first row

# Make a list of the column names that will be regressed
colNamesList <- colnames(gtbl)[-1] # remove "park" from the column list

# Create an empty row for reuse in the loops below
empty.row <- gtbl

####################################################################################
###  Compute r-squared values
####################################################################################

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
    # Get the adjusted r2, slope and y-intercept
    r2.slope.y <- r2SlopeYint(xvar, f.df[, colNamesList[i]])
    
    empty.row[, i+1] <- r2.slope.y[[1]]
    s.y.empty[, i+i] <- r2.slope.y[[2]]
    s.y.empty[, i+i+1] <- r2.slope.y[[3]]
  }
  
  
  
  {
    ####  Do the combo band regressions
    md <- lm(xvar ~ f.df$min3 + f.df$min4 + f.df$min8)
    empty.row$min348 <- round(summary(md)$r.squared, 3)
    
    md <- lm(xvar ~ f.df$max3 + f.df$max4 + f.df$max8)
    empty.row$max348 <- round(summary(md)$r.squared, 3)
    
    md <- lm(xvar ~ f.df$min3 + f.df$max3)
    empty.row$minmax33 <- round(summary(md)$r.squared, 3)
  }
  
  
  ######  PCA  ###########################################################################################
  library(psych) # for PCA
  library(nFactors) # for scree plot
  library(QuantPsyc) # for lm.beta
  
  
  # Function to get the r2 of PCA components, name the highest loaded bands and regress those
  PCA_r2 <- function(df, xval){
    # Rotate with 3 dim extracted
    # Run the rotated PCA 
    fit.rot <- principal(df, nfactors = 3, rotate="varimax") # c(-13:-20)
    
    # Put the scores into the df as a new column 
    df$Comp1 <- fit.rot$scores[, 1]
    df$Comp2 <- fit.rot$scores[, 2]
    df$Comp3 <- fit.rot$scores[, 3]
    
    # Regress the components vs. woody cover
    vals <- c(r2(xval, df$Comp1), r2(xval, df$Comp2), r2(xval, df$Comp3))
    # Do all three
    md <- lm(xval ~ df$Comp1 + df$Comp2 + df$Comp3) 
    vals <- c(vals, round(summary(md)$r.squared, 3))
    
    # Now get the highest loaded bands 
    high1 <- names(which(abs(fit.rot$loadings[, 1]) == max(abs(fit.rot$loadings[, 1]))))
    high2 <- names(which(abs(fit.rot$loadings[, 3]) == max(abs(fit.rot$loadings[, 3])))) # the numbering is right - for some reason the columns are switched btwn RC2 and 3
    high3 <- names(which(abs(fit.rot$loadings[, 2]) == max(abs(fit.rot$loadings[, 2]))))
    
    # Put those names in the val list
    vals <- c(vals, high1, high2, high3)
    
    # Use the band names to get the respective column's data and get the r2
    high1_vals <- df[, which(names(df) == high1)]
    high2_vals <- df[, which(names(df) == high2)]
    high3_vals <- df[, which(names(df) == high3)]
    # r2
    vals <- c(vals, r2(xval, high1_vals), r2(xval, high2_vals), r2(xval, high3_vals))
    # Do all three
    md <- lm(xval ~ high1_vals + high2_vals + high3_vals)
    vals <- c(vals, round(summary(md)$r.squared, 3))
    
    return(vals)
  }
  
  # Enter the min and max image bands to get PCA regressions etc.
  minPCAresults <- PCA_r2(f.df[, c(8:16)], xvar)
  maxPCAresults <- PCA_r2(f.df[, c(17:25)], xvar)
  minmaxPCAresults <- PCA_r2(f.df[, c(8:25)], xvar)
  # Break the spatial context up into measures across images
  desired.cols.min <- which(substr(colnames(f.df), 1, 3) == "min" & nchar(colnames(f.df)) > 4)
  desired.cols.max <- which(substr(colnames(f.df), 1, 3) == "max" & nchar(colnames(f.df)) > 4)
  min.tbl <- f.df[, desired.cols.min]
  max.tbl <- f.df[, desired.cols.max]
  # Use only the rows with complete cases
  spatialContextPCAresults.min <- PCA_r2(min.tbl[complete.cases(min.tbl), ] , xvar[complete.cases(min.tbl)])
  spatialContextPCAresults.max <- PCA_r2(max.tbl[complete.cases(max.tbl), ] , xvar[complete.cases(max.tbl)])
  
  # Plug the results into the empty row and bind it to the growing df
  # Min NDVI image
  empty.row$minComp1 <- as.numeric(minPCAresults[1])
  empty.row$minComp2 <- as.numeric(minPCAresults[2])
  empty.row$minComp3 <- as.numeric(minPCAresults[3])
  empty.row$minAll3Comps <- as.numeric(minPCAresults[4])
  empty.row$min_Comp1_var <- minPCAresults[5]
  empty.row$min_Comp2_var <- minPCAresults[6]
  empty.row$min_Comp3_var <- minPCAresults[7]
  empty.row$min_Comp1_var_r2 <- as.numeric(minPCAresults[8])
  empty.row$min_Comp2_var_r2 <- as.numeric(minPCAresults[9])
  empty.row$min_Comp3_var_r2 <- as.numeric(minPCAresults[10])
  empty.row$min_All3Comps_var_r2 <- as.numeric(minPCAresults[11])
  
  # Max NDVI image
  empty.row$maxComp1 <- as.numeric(maxPCAresults[1])
  empty.row$maxComp2 <- as.numeric(maxPCAresults[2])
  empty.row$maxComp3 <- as.numeric(maxPCAresults[3])
  empty.row$maxAll3Comps <- as.numeric(maxPCAresults[4])
  empty.row$max_Comp1_var <- maxPCAresults[5]
  empty.row$max_Comp2_var <- maxPCAresults[6]
  empty.row$max_Comp3_var <- maxPCAresults[7]
  empty.row$max_Comp1_var_r2 <- as.numeric(maxPCAresults[8])
  empty.row$max_Comp2_var_r2 <- as.numeric(maxPCAresults[9])
  empty.row$max_Comp3_var_r2 <- as.numeric(maxPCAresults[10])
  empty.row$max_All3Comps_var_r2 <- as.numeric(maxPCAresults[11])
  
  # Min Max image
  empty.row$minmaxComp1 <- as.numeric(minmaxPCAresults[1])
  empty.row$minmaxComp2 <- as.numeric(minmaxPCAresults[2])
  empty.row$minmaxComp3 <- as.numeric(minmaxPCAresults[3])
  empty.row$minmaxAll3Comps <- as.numeric(minmaxPCAresults[4])
  empty.row$minmax_Comp1_var <- minmaxPCAresults[5]
  empty.row$minmax_Comp2_var <- minmaxPCAresults[6]
  empty.row$minmax_Comp3_var <- minmaxPCAresults[7]
  empty.row$minmax_Comp1_var_r2 <- as.numeric(minmaxPCAresults[8])
  empty.row$minmax_Comp2_var_r2 <- as.numeric(minmaxPCAresults[9])
  empty.row$minmax_Comp3_var_r2 <- as.numeric(minmaxPCAresults[10])
  empty.row$minmax_All3Comps_var_r2 <- as.numeric(minmaxPCAresults[11])
  
  # Spatial Context PCA data: min image
  empty.row$min_SCcomp1 <- as.numeric(spatialContextPCAresults.min[1])
  empty.row$min_SCcomp2 <- as.numeric(spatialContextPCAresults.min[2])
  empty.row$min_SCcomp3 <- as.numeric(spatialContextPCAresults.min[3])
  empty.row$min_SCall3Comps <- as.numeric(spatialContextPCAresults.min[4])
  empty.row$min_SC_Comp1_var <- spatialContextPCAresults.min[5]
  empty.row$min_SC_Comp2_var <- spatialContextPCAresults.min[6]
  empty.row$min_SC_Comp3_var <- spatialContextPCAresults.min[7]
  empty.row$min_SC_Comp1_var_r2 <- as.numeric(spatialContextPCAresults.min[8])
  empty.row$min_SC_Comp2_var_r2 <- as.numeric(spatialContextPCAresults.min[9])
  empty.row$min_SC_Comp3_var_r2 <- as.numeric(spatialContextPCAresults.min[10])
  empty.row$min_SC_All3Comps_var_r2 <- as.numeric(spatialContextPCAresults.min[11])
  
  # Spatial Context PCA data: max image
  empty.row$max_SCcomp1 <- as.numeric(spatialContextPCAresults.max[1])
  empty.row$max_SCcomp2 <- as.numeric(spatialContextPCAresults.max[2])
  empty.row$max_SCcomp3 <- as.numeric(spatialContextPCAresults.max[3])
  empty.row$max_SCall3Comps <- as.numeric(spatialContextPCAresults.max[4])
  empty.row$max_SC_Comp1_var <- spatialContextPCAresults.max[5]
  empty.row$max_SC_Comp2_var <- spatialContextPCAresults.max[6]
  empty.row$max_SC_Comp3_var <- spatialContextPCAresults.max[7]
  empty.row$max_SC_Comp1_var_r2 <- as.numeric(spatialContextPCAresults.max[8])
  empty.row$max_SC_Comp2_var_r2 <- as.numeric(spatialContextPCAresults.max[9])
  empty.row$max_SC_Comp3_var_r2 <- as.numeric(spatialContextPCAresults.max[10])
  empty.row$max_SC_All3Comps_var_r2 <- as.numeric(spatialContextPCAresults.max[11])
}

# Add the park name as the front of the list
empty.row$park <- park.name
s.y.empty$park <- park.name

# Rbind (unless it's the first)
if(k == 1){
  gtbl <- empty.row
  gtbl.slope.yint <- s.y.empty
}else{
  gtbl <- rbind(gtbl, empty.row)
  gtbl.slope.yint <- rbind(gtbl.slope.yint, s.y.empty)
}



# Get the standard deviations of the slope and y-int
# Don't take the row that corresponds to all parks (number 12)
x <- gtbl.slope.yint[-nrow(gtbl.slope.yint), ]
sd.slope.yint <- x[1, ]
# Loop through the columns and calc the sd
for(n in 1:ncol(sd.slope.yint)){
  if(n == 1){
    sd.slope.yint[, n] <- "stdev"
  }else{
    sd.slope.yint[, n] <- sd(x[, n])
  }
}

# Now transpose the table for easier viewing
sd.slope.yint <- t(sd.slope.yint)



###########################################################################################
####  Dummy Variable Analysis  
###########################################################################################

# Take the table from before with the slope and yint column names and
# run a loop on the relevant final.df columns, filling in the slope and yint
# table with significance values as it goes (remove the AllParks row)
dummy.table <- gtbl.slope.yint[-nrow(gtbl.slope.yint), ]

# Loop through the final.df.dums relevant columns 
ind.var.tbl <- final.df[, -c(2:7)] # remove columns without variable values 

for(i in 2:ncol(ind.var.tbl)){
  ## Run the regression
  # Get the variable name
  var.name <- names(ind.var.tbl)[i]
  
  # run the dummy regression (park names are used as the dummy)
  m <- lm(hand.tree ~ final.df[, var.name] * as.factor(park), data = final.df) # I checked to make sure this was doing what it was supposed to
  summary(m)
  
  # Get the pvalues
  # Chobe will have the slope and y int of the first row
  # Then the others are adjustment to that
  # the rows with ":" are the slope factors 
  tb <- as.data.frame(summary(m)$coefficients)
  
  # Plug the pvalues into the appropriate columns in the dummy.table 
  # Get the slope and yint p values (the second half of the tb)
  sl.vals <- tb[c(2, ((nrow(tb)/2)+2) : nrow(tb)), 4] # The second row is slope for Chobe
  y.vals <- tb[c(1, 3:(nrow(dummy.table)+1)), 4]  # dummy table get's me to 11 faster than trying to halve the tb - 1
  
  # Find the columns in the dummy.table
  cols <- which(substr(names(dummy.table), 1, nchar(var.name)) == var.name)
  names(dummy.table)[cols] # show the names
  
  # Plug in the values (the first column is the slope, the second is the y-int)
  dummy.table[ , cols[1]] <- round(sl.vals, 5)
  dummy.table[ , cols[2]] <- round(y.vals, 5)
}


### Add a row that says whether all the results above it were
### above 0.05. Have another row that is the total number
### above 0.05
### Remember that Chobe needs to be less than 0.05, so omit that row (first row)
#####  max6 was the first to have a high number of sig similar equations @ 8 out of 11  <- 

# When you require that all the other slopes and yints be above 0.05, no slope/yint combos satisfy that
# So let's set a variable that allows us to change how many out of 10 we want (Chobe is excluded, remember (=> 10 and not 11))
tot.over05 <- 8

{
  # Create the empty rows
  yes.no.row <- dummy.table[c(1:2), ]
  yes.no.row$park <- c("y.n", "num_above_0.05")
  
  # Loop through each column and determine if all are
  for(i in 2:ncol(dummy.table)){
    # Get the number above 0.05
    above05 <- sum(dummy.table[c(2:11), i] > 0.05)
    chobeSig <- sum(dummy.table[1, i] <= 0.05)
    sigSum <- above05 + chobeSig
    # Set whether all were above 0.05
    if(chobeSig == 1 & sigSum >= tot.over05){ #Chobe has to be significant no matter what, otherwise there is no sig trend
      yn.var <- "Yes"
    }else{
      yn.var <- "No"
    }
    # Set the values
    yes.no.row[1, i] <- yn.var
    yes.no.row[2, i] <- sigSum
  }
  
  # Rbind the two rows to the dummy.table - the 1:11 allows this to rewrite
  dummy.table <- rbind(dummy.table[c(1:11), ], yes.no.row)
  
  # Get the rows that had Yes
  yes.cols <- dummy.table[, which(dummy.table[12,] == "Yes")]
  
  ### Get only the columns that are slope and yint matches
  # Get the colnames
  y.col.names <- names(yes.cols)
  
  # Loop through the names and take off slope/yint
  for(i in 1:length(y.col.names)){
    if(grepl("slope", y.col.names[i])){
      y.col.names[i] <- substr(y.col.names[i], 1, nchar(y.col.names[i])-6)
    }else{
      y.col.names[i] <- substr(y.col.names[i], 1, nchar(y.col.names[i])-5)
    }
  }
  
  # Get the columns that have duplicates
  dup.cols <- duplicated(y.col.names) | duplicated(y.col.names, fromLast = T) # you have to reverse the dup function to get both dup values
  # Subset
  statistically.similar.equations <- yes.cols[, dup.cols]
  
  # Add the parknames back on 
  statistically.similar.equations$park <- c(park.names, ">= 8?", "# non-sig")
  statistically.similar.equations <- cbind(statistically.similar.equations, statistically.similar.equations)[, c(3:5)]
}

### OK, so max6 had the most parks with similar relationships to woody cover, but how good were those relationships?
# i.e. what were the r-squareds like?
gtbl.max6.only <- gtbl[, c("park", "max6")]

# Much of the above was reported in this ppt to Kyla
# E:\Dropbox\Permanent\Grad School\Meetings\Kyla\20181130_Testing Transition Image n Search for consistent regression\20181130_Testing Transition Image & Search for consistent regression





# # Test dummy
# dep <- c(2, 9, 15, 21, 1, 3, 5, 7)
# ind1 <- c(1, 2, 6, 8)
# ind2 <- ind1 + 5
# dum1 <- c(0, 0, 0, 0, 1, 1, 1, 1)
# inds <- c(ind1, ind2)
# df <- as.data.frame(cbind(dep, inds, dum1))
# mdd <- lm(inds ~ dep * as.factor(dum1), data = df)
# summary(mdd)
# # Test y int
# mdd <- lm(inds ~ dep * dum1)
# summary(mdd)
# m <- lm(hand.tree ~ min1 %in% park, data = final.df)
# summary(m)
# m <- lm(hand.tree ~ min1 * as.factor(park), data = final.df)
# summary(m)


# Print the table
# write.csv(gtbl, "~/Dropbox/Permanent/Grad School/Projects/EleTree/analysis/Landsat_v_woody_cover_regressions/CorrelationCoefficent_ALLpark_SpatContxt.csv")
write.csv(gtbl, "E:/Dropbox/Permanent/Grad School/Projects/EleTree/analysis/Landsat_v_woody_cover_regressions/CorrelationCoefficent_ALLpark_SpatContxt.csv")



#########################################################################################
####  Evaluate the r-squared results
#########################################################################################

# First, take out columns that can't be quantitatively evaluated
# Which columns have text instead of numbers?
no.text.tbl <- subset(gtbl, select = -c(min_Comp1_var, min_Comp2_var, min_Comp3_var, 
                                        max_Comp1_var, max_Comp2_var, max_Comp3_var, 
                                        minmax_Comp1_var, minmax_Comp2_var, minmax_Comp3_var,
                                        min_SC_Comp1_var, min_SC_Comp2_var, min_SC_Comp3_var,
                                        max_SC_Comp1_var, max_SC_Comp2_var, max_SC_Comp3_var))

# Which bands are best on average?
# Make table and fill
mx <- matrix(NA, ncol(no.text.tbl)-1, 1)
df <- as.data.frame(mx)
row.names(df) <- names(no.text.tbl)[c(2:ncol(no.text.tbl))]
band.performance.r2 <- df


# Take the mean r-squared of reflectance vs WC for each park (requires getting rid of the last row)
for(i in 1:nrow(band.performance.r2)){
  band.performance.r2[i, 1] <- mean(no.text.tbl[c(1:11), i+1]) 
}




##################################################################################################
#####  Test Models using VEcv  ###################################################################
##################################################################################################

# Approach: Use 70% of the points to make a simple linear regression. Then use
#           that regression equation to predict woody cover at 30% of the points.
#           Evaluate the results using VEcv 
#           - Do this for individual parks and for all parks together 


# Create an empty table with the same column names and add PCA columns
VE.tbl <- final.df
VE.tbl <- VE.tbl[, -c(2:7)] # remove columns
VE.tbl <- VE.tbl[-c(2:NROW(VE.tbl)), ] # remove all but first row

# Add all of the other columns here
{
  # Add the 3 band regression (bands 3, 4 and 8), which showed the highest r-squared values 
  # after running this code previously
  VE.tbl$min348 <- 0
  VE.tbl$max348 <- 0
  VE.tbl$minmax33 <- 0 # This one is just band three from both images 
  
  # Add PCA columns
  VE.tbl$minComp1 <- 0
  VE.tbl$minComp2 <- 0
  VE.tbl$minComp3 <- 0
  VE.tbl$minAll3Comps <- 0
  VE.tbl$min_Comp1_var <- NA # these will be actual names of bands 
  VE.tbl$min_Comp2_var <- NA
  VE.tbl$min_Comp3_var <- NA
  VE.tbl$min_Comp1_var_VE <- 0
  VE.tbl$min_Comp2_var_VE <- 0
  VE.tbl$min_Comp3_var_VE <- 0
  VE.tbl$min_All3Comps_var_VE <- 0
  
  
  VE.tbl$maxComp1 <- 0
  VE.tbl$maxComp2 <- 0
  VE.tbl$maxComp3 <- 0
  VE.tbl$maxAll3Comps <- 0
  VE.tbl$max_Comp1_var <- NA
  VE.tbl$max_Comp2_var <- NA
  VE.tbl$max_Comp3_var <- NA
  VE.tbl$max_Comp1_var_VE <- 0
  VE.tbl$max_Comp2_var_VE <- 0
  VE.tbl$max_Comp3_var_VE <- 0
  VE.tbl$max_All3Comps_var_VE <- 0
  
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




#########################################################################################
####  Evaluate the results
#########################################################################################

# First, take out columns that can't be quantitatively evaluated
# Which columns have text instead of numbers?
no.text.tbl.VE <- subset(VE.tbl, select = -c(min_Comp1_var, min_Comp2_var, min_Comp3_var,
                                             max_Comp1_var, max_Comp2_var, max_Comp3_var,
                                             minmax_Comp1_var, minmax_Comp2_var, minmax_Comp3_var,
                                             min_SC_Comp1_var, min_SC_Comp2_var, min_SC_Comp3_var,
                                             max_SC_Comp1_var, max_SC_Comp2_var, max_SC_Comp3_var))

# Which bands are best on average?
# Make table and fill
mx <- matrix(NA, ncol(no.text.tbl.VE)-1, 1)
df <- as.data.frame(mx)
row.names(df) <- names(no.text.tbl.VE)[c(2:ncol(no.text.tbl.VE))]
band.performance.VE <- df


# Take the mean VEcv of reflectance vs WC for each park (requires getting rid of the last row)
for(i in 1:nrow(band.performance.VE)){
  band.performance.VE[i, 1] <- mean(no.text.tbl.VE[c(1:11), i+1]) 
}



##################################################################################################
#####  Test Models using E1  ###################################################################
##################################################################################################

# Approach: Use 70% of the points to make a simple linear regression. Then use
#           that regression equation to predict woody cover at 30% of the points.
#           Evaluate the results using VEcv 
#           - Do this for individual parks and for all parks together 


# Create an empty table with the same column names and add PCA columns
E1.tbl <- final.df
E1.tbl <- E1.tbl[, -c(2:7)] # remove columns
E1.tbl <- E1.tbl[-c(2:NROW(E1.tbl)), ] # remove all but first row

# Add all of the other columns here
{
  # Add the 3 band regression (bands 3, 4 and 8), which showed the highest r-squared values 
  # after running this code previously
  E1.tbl$min348 <- 0
  E1.tbl$max348 <- 0
  E1.tbl$minmax33 <- 0 # This one is just band three from both images 
  
  # Add PCA columns
  E1.tbl$minComp1 <- 0
  E1.tbl$minComp2 <- 0
  E1.tbl$minComp3 <- 0
  E1.tbl$minAll3Comps <- 0
  E1.tbl$min_Comp1_var <- NA # these will be actual names of bands 
  E1.tbl$min_Comp2_var <- NA
  E1.tbl$min_Comp3_var <- NA
  E1.tbl$min_Comp1_var_E1 <- 0
  E1.tbl$min_Comp2_var_E1 <- 0
  E1.tbl$min_Comp3_var_E1 <- 0
  E1.tbl$min_All3Comps_var_E1 <- 0
  
  
  E1.tbl$maxComp1 <- 0
  E1.tbl$maxComp2 <- 0
  E1.tbl$maxComp3 <- 0
  E1.tbl$maxAll3Comps <- 0
  E1.tbl$max_Comp1_var <- NA
  E1.tbl$max_Comp2_var <- NA
  E1.tbl$max_Comp3_var <- NA
  E1.tbl$max_Comp1_var_E1 <- 0
  E1.tbl$max_Comp2_var_E1 <- 0
  E1.tbl$max_Comp3_var_E1 <- 0
  E1.tbl$max_All3Comps_var_E1 <- 0
  
  E1.tbl$minmaxComp1 <- 0
  E1.tbl$minmaxComp2 <- 0
  E1.tbl$minmaxComp3 <- 0
  E1.tbl$minmaxAll3Comps <- 0
  E1.tbl$minmax_Comp1_var <- NA
  E1.tbl$minmax_Comp2_var <- NA
  E1.tbl$minmax_Comp3_var <- NA
  E1.tbl$minmax_Comp1_var_E1 <- 0
  E1.tbl$minmax_Comp2_var_E1 <- 0
  E1.tbl$minmax_Comp3_var_E1 <- 0
  E1.tbl$minmax_All3Comps_var_E1 <- 0
}

# Create an empty row for reuse in the loops below
empty.row <- E1.tbl


####################################################################################
###  Create the models and compute the E1 values
####################################################################################

# First, write a function that will create the model and calc the VEcv
model_E1 <- function(xval, yval){
  
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
  return(round(E1(df.test$xval, wc.predicted), 5))
  
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
    empty.row[, i+1] <- model_E1(xval = f.df$hand.tree[comp.cases], yval = f.df[comp.cases, cNumber]) 
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
    empty.row$min348 <- E1(df.test$hand.tree, wc.predicted)
    
    # Do the next (max)
    md <- lm(hand.tree ~ max3 + max4 + max8, data = df.train)
    # Subset to just the data you want (the band)
    test_band <- subset(df.test, select = c(max3, max4, max8))
    # Predict WC using the model and the testing data
    wc.predicted <- predict(md, test_band)
    # Return the VEcv of the predictions 
    empty.row$max348 <- E1(df.test$hand.tree, wc.predicted)
    
    # Do the next (minmax)
    md <- lm(hand.tree ~ min3 + max3, data = df.train)
    # Subset to just the data you want (the band)
    test_band <- subset(df.test, select = c(min3, max3))
    # Predict WC using the model and the testing data
    wc.predicted <- predict(md, test_band)
    # Return the VEcv of the predictions 
    empty.row$minmax33 <- E1(df.test$hand.tree, wc.predicted)
    
  }
  
  
  ######  PCA  ###########################################################################################
  library(psych) # for PCA
  library(nFactors) # for scree plot
  library(QuantPsyc) # for lm.beta
  
  
  # Function to model woody cover using the PCA results and report the E1 of the model
  PCA_E1 <- function(df, xval){
    # Rotate with 3 dim extracted
    # Run the rotated PCA 
    fit.rot <- principal(df, nfactors = 3, rotate = "varimax") # c(-13:-20)
    
    # Put the scores into the df as a new column 
    df$Comp1 <- fit.rot$scores[, 1]
    df$Comp2 <- fit.rot$scores[, 2]
    df$Comp3 <- fit.rot$scores[, 3]
    
    # Model woody cover using the components and get the E1
    vals <- c(model_E1(xval, df$Comp1), model_E1(xval, df$Comp2), model_E1(xval, df$Comp3))
    
    
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
    
    # Add the model E1 to the list 
    vals <- c(vals, round(E1(df.test$x.test, wc.predicted), 3)) # x.test is woody cover (hand.tree)
    
    # Now get the highest loaded bands 
    high1 <- names(which(abs(fit.rot$loadings[, 1]) == max(abs(fit.rot$loadings[, 1]))))
    high2 <- names(which(abs(fit.rot$loadings[, 3]) == max(abs(fit.rot$loadings[, 3])))) # the numbering is right - for some reason the columns are switched btwn RC2 and 3
    high3 <- names(which(abs(fit.rot$loadings[, 2]) == max(abs(fit.rot$loadings[, 2]))))
    
    # Put those names in the val list
    vals <- c(vals, high1, high2, high3)
    
    # Use the band names to get the respective column's data and get the E1 of the models
    high1_vals <- df[, which(names(df) == high1)]
    high2_vals <- df[, which(names(df) == high2)]
    high3_vals <- df[, which(names(df) == high3)]
    
    # E1
    vals <- c(vals, model_E1(xval, high1_vals), model_E1(xval, high2_vals), model_E1(xval, high3_vals))
    
    
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
    
    # Add the model E1 to the list 
    vals <- c(vals, round(E1(threeVar.test$x.test, wc.predicted), 3))
    
    return(vals)
  }
  
  
  
  {
    # Enter the min image bands to get PCA regressions etc.
    minPCAresults <- PCA_E1(f.df[, c(8:16)], f.df$hand.tree)
    maxPCAresults <- PCA_E1(f.df[, c(17:25)], f.df$hand.tree)
    minmaxPCAresults <- PCA_E1(f.df[, c(8:25)], f.df$hand.tree)
    # Break the spatial context up into measures across images
    desired.cols.min <- which(substr(colnames(f.df), 1, 3) == "min" & nchar(colnames(f.df)) > 4)
    desired.cols.max <- which(substr(colnames(f.df), 1, 3) == "max" & nchar(colnames(f.df)) > 4)
    min.tbl <- f.df[, desired.cols.min]
    max.tbl <- f.df[, desired.cols.max]
    # Use only the rows with complete cases
    spatialContextPCAresults.min <- PCA_E1(min.tbl[complete.cases(min.tbl), ] , f.df$hand.tree[complete.cases(min.tbl)])   
    spatialContextPCAresults.max <- PCA_E1(max.tbl[complete.cases(max.tbl), ] , f.df$hand.tree[complete.cases(max.tbl)]) 
    
    # Plug the values into the empty row and bind it to the growing df
    # Min NDVI image
    empty.row$minComp1 <- as.numeric(minPCAresults[1])
    empty.row$minComp2 <- as.numeric(minPCAresults[2])
    empty.row$minComp3 <- as.numeric(minPCAresults[3])
    empty.row$minAll3Comps <- as.numeric(minPCAresults[4])
    empty.row$min_Comp1_var <- minPCAresults[5]
    empty.row$min_Comp2_var <- minPCAresults[6]
    empty.row$min_Comp3_var <- minPCAresults[7]
    empty.row$min_Comp1_var_E1 <- as.numeric(minPCAresults[8])
    empty.row$min_Comp2_var_E1 <- as.numeric(minPCAresults[9])
    empty.row$min_Comp3_var_E1 <- as.numeric(minPCAresults[10])
    empty.row$min_All3Comps_var_E1 <- as.numeric(minPCAresults[11])
    
    # Max NDVI image
    empty.row$maxComp1 <- as.numeric(maxPCAresults[1])
    empty.row$maxComp2 <- as.numeric(maxPCAresults[2])
    empty.row$maxComp3 <- as.numeric(maxPCAresults[3])
    empty.row$maxAll3Comps <- as.numeric(maxPCAresults[4])
    empty.row$max_Comp1_var <- maxPCAresults[5]
    empty.row$max_Comp2_var <- maxPCAresults[6]
    empty.row$max_Comp3_var <- maxPCAresults[7]
    empty.row$max_Comp1_var_E1 <- as.numeric(maxPCAresults[8])
    empty.row$max_Comp2_var_E1 <- as.numeric(maxPCAresults[9])
    empty.row$max_Comp3_var_E1 <- as.numeric(maxPCAresults[10])
    empty.row$max_All3Comps_var_E1 <- as.numeric(maxPCAresults[11])
    
    # Min Max image
    empty.row$minmaxComp1 <- as.numeric(minmaxPCAresults[1])
    empty.row$minmaxComp2 <- as.numeric(minmaxPCAresults[2])
    empty.row$minmaxComp3 <- as.numeric(minmaxPCAresults[3])
    empty.row$minmaxAll3Comps <- as.numeric(minmaxPCAresults[4])
    empty.row$minmax_Comp1_var <- minmaxPCAresults[5]
    empty.row$minmax_Comp2_var <- minmaxPCAresults[6]
    empty.row$minmax_Comp3_var <- minmaxPCAresults[7]
    empty.row$minmax_Comp1_var_E1 <- as.numeric(minmaxPCAresults[8])
    empty.row$minmax_Comp2_var_E1 <- as.numeric(minmaxPCAresults[9])
    empty.row$minmax_Comp3_var_E1 <- as.numeric(minmaxPCAresults[10])
    empty.row$minmax_All3Comps_var_E1 <- as.numeric(minmaxPCAresults[11])
    
    # Spatial Context PCA data: min image
    empty.row$min_SCcomp1 <- as.numeric(spatialContextPCAresults.min[1])
    empty.row$min_SCcomp2 <- as.numeric(spatialContextPCAresults.min[2])
    empty.row$min_SCcomp3 <- as.numeric(spatialContextPCAresults.min[3])
    empty.row$min_SCall3Comps <- as.numeric(spatialContextPCAresults.min[4])
    empty.row$min_SC_Comp1_var <- spatialContextPCAresults.min[5]
    empty.row$min_SC_Comp2_var <- spatialContextPCAresults.min[6]
    empty.row$min_SC_Comp3_var <- spatialContextPCAresults.min[7]
    empty.row$min_SC_Comp1_var_E1 <- as.numeric(spatialContextPCAresults.min[8])
    empty.row$min_SC_Comp2_var_E1 <- as.numeric(spatialContextPCAresults.min[9])
    empty.row$min_SC_Comp3_var_E1 <- as.numeric(spatialContextPCAresults.min[10])
    empty.row$min_SC_All3Comps_var_E1 <- as.numeric(spatialContextPCAresults.min[11])
    
    # Spatial Context PCA data: max image
    empty.row$max_SCcomp1 <- as.numeric(spatialContextPCAresults.max[1])
    empty.row$max_SCcomp2 <- as.numeric(spatialContextPCAresults.max[2])
    empty.row$max_SCcomp3 <- as.numeric(spatialContextPCAresults.max[3])
    empty.row$max_SCall3Comps <- as.numeric(spatialContextPCAresults.max[4])
    empty.row$max_SC_Comp1_var <- spatialContextPCAresults.max[5]
    empty.row$max_SC_Comp2_var <- spatialContextPCAresults.max[6]
    empty.row$max_SC_Comp3_var <- spatialContextPCAresults.max[7]
    empty.row$max_SC_Comp1_var_E1 <- as.numeric(spatialContextPCAresults.max[8])
    empty.row$max_SC_Comp2_var_E1 <- as.numeric(spatialContextPCAresults.max[9])
    empty.row$max_SC_Comp3_var_E1 <- as.numeric(spatialContextPCAresults.max[10])
    empty.row$max_SC_All3Comps_var_E1 <- as.numeric(spatialContextPCAresults.max[11])
  }
  
  # Add the park name as the front of the list
  empty.row$park <- park.name
  
  # Rbind (unless it's the first)
  if(k == 1){
    E1.tbl <- empty.row
  }else{
    E1.tbl <- rbind(E1.tbl, empty.row)
  }
  
}

# Print the table
# write.csv(E1.tbl, "~/Dropbox/Permanent/Grad School/Projects/EleTree/analysis/Landsat_v_woody_cover_regressions/Regression_Models_E1_ALLpark_SpatContxt.csv")
write.csv(E1.tbl, "E:/Dropbox/Permanent/Grad School/Projects/EleTree/analysis/Landsat_v_woody_cover_regressions/Regression_Models_E1_ALLpark_SpatContxt.csv")




#########################################################################################
####  Evaluate the results
#########################################################################################

# First, take out columns that can't be quantitatively evaluated
# Which columns have text instead of numbers?
no.text.tbl.E1 <- subset(E1.tbl, select = -c(min_Comp1_var, min_Comp2_var, min_Comp3_var,
                                             max_Comp1_var, max_Comp2_var, max_Comp3_var,
                                             minmax_Comp1_var, minmax_Comp2_var, minmax_Comp3_var,
                                             min_SC_Comp1_var, min_SC_Comp2_var, min_SC_Comp3_var,
                                             max_SC_Comp1_var, max_SC_Comp2_var, max_SC_Comp3_var))

# Which bands are best on average?
# Make table and fill
mx <- matrix(NA, ncol(no.text.tbl.E1)-1, 1)
df <- as.data.frame(mx)
row.names(df) <- names(no.text.tbl.E1)[c(2:ncol(no.text.tbl.E1))]
band.performance.E1 <- df


# Take the mean E1 of reflectance vs WC for each park (requires getting rid of the last row)
for(i in 1:nrow(band.performance.E1)){
  band.performance.E1[i, 1] <- mean(no.text.tbl.E1[c(1:11), i+1]) 
}




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


