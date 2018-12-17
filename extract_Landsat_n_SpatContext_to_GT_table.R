# Extract Landsat values and the Spatial Context Values at all Reference Data Locations 

# Origin: Landsat vs Ground Truth Data Exploration

# 20181002
# Ryan Nagelkirk

# Description:
# 1. Pull out the Landsat band values from all the ground truth points

########  Reminder #################################################
# This code doesn't attempt to merge point values from overlapping imagery. Instead,
# each image is treated independently and has it's own row in the final.df
# If you want to handle those points that do overlap, that code still has to be written


# Libraries
library(raster)
library(rgdal)
library(stats)

###########################################################################################################
#  Switches 
###########################################################################################################
run.spatial.context <- T # Set this to true if you want to extract the spatial context from each point


###########################################################################################################
#  Variables and Directories 
###########################################################################################################
# Set working directory to ParkData folder on dropbox
setwd("E:/Dropbox/Permanent/Grad School/Projects/EleTree/data/ParkData")
# setwd("~/Dropbox/Permanent/Grad School/Projects/EleTree/data/ParkData")

# Set the folder for all the CSVs
csv.folder <- "E:/Dropbox/Permanent/Grad School/Projects/EleTree/analysis/Landsat_v_woody_cover_regressions/"

# Source the functions
source("E:/Dropbox/Permanent/Grad School/src_functions/src_masterfunctions.R")
# source("~/Dropbox/Permanent/Grad School/src_functions/src_masterfunctions.R")

final.df.csv <- "E:/Dropbox/Permanent/Grad School/Projects/EleTree/analysis/Landsat_v_woody_cover_regressions/TwelvePark_Landsat_n_GT_values.csv"
# final.df.csv <- "~/Dropbox/Permanent/Grad School/Projects/EleTree/analysis/Landsat_v_woody_cover_regressions/TwelvePark_Landsat_n_GT_values.csv"

final.df.w.SpatialContext <- "E:/Dropbox/Permanent/Grad School/Projects/EleTree/analysis/Landsat_v_woody_cover_regressions/TwelvePark_Landsat_n_GT_values_w_SpatContext.csv"
# final.df.w.SpatialContext <- "~/Dropbox/Permanent/Grad School/Projects/EleTree/analysis/Landsat_v_woody_cover_regressions/TwelvePark_Landsat_n_GT_values_w_SpatContext.csv"

final.df.tran <- "E:/Dropbox/Permanent/Grad School/Projects/EleTree/analysis/Landsat_v_woody_cover_regressions/AllTwelveParks_Landsat_n_GT_values_TRANS.csv"
# final.df.tran <- "~/Dropbox/Permanent/Grad School/Projects/EleTree/analysis/Landsat_v_woody_cover_regressions/AllTwelveParks_Landsat_n_GT_values_TRANS.csv"

# List the parks
park.names <- c("Chobe", "Kruger", "Limpopo", "Mpala", "Murchison", "North_Luangwa", "QWE", "Ruaha", "Selous", "Serengeti", "South_Luangwa")

# Identify the folder with the Landsat images
landsat.folder <- "X:/nagelki4/Projects/EleTree/data/Landsat/20km_Buffer"

# List only the images that are max min or tran but not the minmaxtran
landsat.rasters <- list.files(landsat.folder, pattern = "_min_NDVI_wo_indices_allBands.tif", full.names = T)
landsat.rasters <- c(landsat.rasters, list.files(landsat.folder, pattern = "_max_NDVI_wo_indices_allBands.tif", full.names = T))
landsat.rasters <- c(landsat.rasters, list.files(paste0(landsat.folder, "/With_shoulder"), pattern = "_tran_", full.names = T))
landsat.rasters <- landsat.rasters[!grepl("00000000", landsat.rasters)]


#########################################################################################################################
# Loop through the parks and images, extracting the values and plugging them into the table
#########################################################################################################################
for(i in 1:length(park.names)){
  park <- park.names[i]
  
  # Construct point folder
  pt.folder <- paste0("./", park, "/ground_truth/wholepark")
  # Load the points
  pts <- readOGR(pt.folder, "gt_points")
  # Get the number of points
  num.o.pts <- length(pts)
  
  # Get the ground truth
  # List the files
  truth.files <- list.files(pt.folder, pattern = "groundtruth", full.names = T)
  # Get most recent
  truth.file <- sort(truth.files, decreasing = T)[1]
  # Read in the file
  truth.df <- read.csv(truth.file)
  # Remove rows with NA and water
  if(park == "Kruger" | park == "Mpala"){
    truth.df <- truth.df[!duplicated(truth.df$park.pxl.num), ] # this was unique to how the kruger pts file was made
  }else{
    truth.df <- truth.df[complete.cases(truth.df), ]
    truth.df <- truth.df[truth.df$water.flag != 1, ]
    
    # Eliminate ground truth that doesn't have a point in "pts"
    pts.cells <- pts@data$pxl
    truth.cells <- truth.df$park.pxl.num
    
    # Which truth has a value that isn't in pts.cells?
    TF.list <- truth.cells %in% pts.cells
    # Now remove that truth value
    truth.df <- truth.df[TF.list, ]
  }
  
  # The extracted values are going to be ordered by cell number, so do the same here to make the points match up
  truth.df <- truth.df[order(truth.df$park.pxl.num), ]
  
 
  
  # Subset the image list to only images from the current park.
  park.images <- grep(park, landsat.rasters, value = T)
  
  for(j in 1:length(park.images)){
    # Load the raster
    minmax.image <- brick(park.images[j])
    
    # If the data has already been extracted, skip
    if(gsub("_max_NDVI", "", basename(park.images[j]) %in% final.df$image)){
      print(paste0("Skipped: ", park.images[j]))
      next
    }
    
    if(tran.image.included){
      LS.tran.pts <- extract(minmax.image[[c(1:9)]], pts, method = 'simple', df = T, cellnumbers = T)
      
      # Create an empty df for the Landsat values
      mx <- matrix(0, nrow = num.o.pts, ncol = 16)
      c.names <- c("park", "image", "point", "cellnumber", "hand.tree", "hand.tree90", "hand.tree180",
                   "tran1", "tran2", "tran3", "tran4", "tran5", "tran6", "tran7", "tran8", "tran9")
      colnames(mx) <- c.names
      df <- as.data.frame(mx)
      
      # Plug the values into the table
      df$park <- park
      df$image <- gsub("_tran_NDVI", "", basename(park.images[j]))
      df$point <- LS.tran.pts$ID
      df$cellnumber <- LS.tran.pts$cells
      df[ , grep("tran", c.names)] <- LS.tran.pts[ , c(3:11)]
      
    }else{
      # Extract the min and max image points
      LS.min.pts <- extract(minmax.image[[c(10:18)]], pts, method = 'simple', df = T, cellnumbers = T)
      LS.max.pts <- extract(minmax.image[[c(1:9)]], pts, method = 'simple', df = T, cellnumbers = T)
      
      # Make sure the cell numbers agree between the two
      if(all((LS.max.pts$cells == LS.min.pts$cells) != TRUE)){
        print(paste0(park, ": Extracted points did not line up between images. Quit the extraction loop."))
        break
      }
      
      # Create an empty df for the Landsat values
      mx <- matrix(0, nrow = num.o.pts, ncol = 25)
      c.names <- c("park", "image", "point", "cellnumber", "hand.tree", "hand.tree90", "hand.tree180",
                   "min1", "min2", "min3", "min4", "min5", "min6", "min7", "min8", "min9",
                   "max1", "max2", "max3", "max4", "max5", "max6", "max7", "max8", "max9")
      colnames(mx) <- c.names
      df <- as.data.frame(mx)
      
      # Plug the values into the table
      df$park <- park
      df$image <- gsub("_max_NDVI", "", basename(park.images[j]))
      df$point <- LS.min.pts$ID
      df$cellnumber <- LS.min.pts$cells
      df[ , grep("min", c.names)] <- LS.min.pts[ , c(3:11)]
      df[ , grep("max", c.names)] <- LS.max.pts[ , c(3:11)]
    }
    
    if(park == "Kruger" | park == "Mpala"){
      df$hand.tree <- truth.df$hand.tree
      df$hand.tree90 <- NA
      df$hand.tree180 <- NA
    }else{
      df$hand.tree <- truth.df$hand.tree.30
      df$hand.tree90 <- truth.df$hand.tree.90
      df$hand.tree180 <- truth.df$hand.tree.180
    }
    
    if(tran.image.included){
      # Now eliminate incomplete cases
      df <- df[!is.na(df$tran1), ]
    }else{
      # Now eliminate incomplete cases
      df <- df[!is.na(df$min1), ]
      
      # Make sure the min and max images are actually the min and max
      if(mean((df$min5-df$min4)/(df$min5 + df$min4)) >  mean((df$max5-df$max4)/(df$max5 + df$max4))){
        print(paste0(park, ": The min and max NDVI images are reversed. Stopped here."))
        break
      }
    }
    
    # Now join the tables
    if(i == 1 & j == 1){
      final.df <- df
    }else{
      final.df <- rbind(final.df, df, make.row.names = F)
    }
  }
  # Write the file
  if(tran.image.included){
    write.csv(final.df, "E:/Dropbox/Permanent/Grad School/Projects/EleTree/analysis/Landsat_v_woody_cover_regressions/AllTwelveParks_Landsat_n_GT_values_TRANS.csv", row.names = FALSE)
  }else{
    write.csv(final.df, "E:/Dropbox/Permanent/Grad School/Projects/EleTree/analysis/Landsat_v_woody_cover_regressions/TwelvePark_Landsat_n_GT_values.csv", row.names = FALSE)
  }
}




#########################################################################################################################
# Loop through the parks and images, extracting the values at different ranges
#########################################################################################################################

# Only run this if elected, and then only for the min and max images (no tran image)
if(run.spatial.context && !tran.image.included){
  for(i in 1:length(park.names)){
    park <- park.names[i]
    
    # Construct point folder
    pt.folder <- paste0("./", park, "/ground_truth/wholepark")
    # Load the points
    pts <- readOGR(pt.folder, "gt_points")
    # Get the number of points
    num.o.pts <- length(pts)
    
    
    # Subset the image list to only images from the current park.
    park.images <- grep(park, landsat.rasters, value = T)
    
    # Make a list of radii
    rad.list <- c(50, 100, 200, 500, 1000)
    
    # Do the extractions for all the images of the park
    for(j in 1:length(park.images)){
      
      # Load the raster
      minmax.image <- brick(park.images[j])
      
      # # If the data has already been extracted, skip
      # if(gsub("_max_NDVI", "", basename(park.images[j]) %in% final.df$image)){
      #   print(paste0("Skipped: ", park.images[j]))
      #   next
      # }
      
      # Do the extraction and calcs for each radius
      for(x in rad.list){
        
        # If values already extracted, skip the extraction and load the csv that was saved
        if(file.exists(paste0(csv.folder, park, "_", substr(basename(park.images[j]), 1, nchar(basename(park.images[j])) - 4), "_minImage_", x, "_extractionRadius.csv"))){
          # Print the status
          print(paste0("Currently working on ", x, "-meter  ", substr(basename(park.images[j]), 1, nchar(basename(park.images[j])) - 4)))
          LS.min.pts <- read.csv(paste0(csv.folder, park, "_", substr(basename(park.images[j]), 1, nchar(basename(park.images[j])) - 4), "_minImage_", x, "_extractionRadius.csv"))
          LS.max.pts <- read.csv(paste0(csv.folder, park, "_", substr(basename(park.images[j]), 1, nchar(basename(park.images[j])) - 4), "_maxImage_", x, "_extractionRadius.csv"))
        }else{
          # Print the status
          print(paste0("Currently working on ", x, "-meter  ", substr(basename(park.images[j]), 1, nchar(basename(park.images[j])) - 4)))
          
          # Extract the 90 x 90
          LS.min.pts <- extract(minmax.image[[c(10:18)]], pts, method = 'simple', df = T, cellnumbers = T, buffer = x)
          LS.max.pts <- extract(minmax.image[[c(1:9)]], pts, method = 'simple', df = T, cellnumbers = T, buffer = x)
          
          # Write the files
          write.csv(LS.min.pts, paste0(csv.folder, park, "_",
                                       substr(basename(park.images[j]), 1, nchar(basename(park.images[j])) - 4), "_minImage_", x, "_extractionRadius.csv"), row.names = FALSE)
          write.csv(LS.max.pts, paste0(csv.folder, park, "_",
                                       substr(basename(park.images[j]), 1, nchar(basename(park.images[j])) - 4), "_maxImage_", x, "_extractionRadius.csv"), row.names = FALSE)
        }
        
        # Take both the final.df and the extracted points, computes the stat and plugs it into the final.df
        # First make into df
        LS.min.pts <- as.data.frame(LS.min.pts)
        LS.max.pts <- as.data.frame(LS.max.pts)
        
        
        # Get the unique IDs, which correspond to the points
        pt.list <- unique(LS.min.pts$ID)
        
        # Subset the final.df to the park
        park.df <- final.df[final.df$park == park & final.df$image == basename(park.images[j]), ]
        
        # Loop through the points, computing the stat and plugging it in
        for(k in pt.list){
          # Get the rows with that point
          min.pt.vals <- LS.min.pts[LS.min.pts$ID == k, ]
          max.pt.vals <- LS.max.pts[LS.max.pts$ID == k, ]
          
          # Get the center cell
          center.cell <- park.df$cellnumber[park.df$point == k]
          
          # Divide the center pixel by the mean of the green band (was the highest weighted band in the allpark PCA first component)
          # Also divide by the mean brightness of the pixels
          # Also divide by the mean NDVI
          
          
          #######  Green
          # # Get the single green value and the mean of the pixels around it
          min.green <- mean(min.pt.vals[, 5][min.pt.vals$cells == center.cell])
          max.green <- mean(max.pt.vals[, 5][max.pt.vals$cells == center.cell])
          
          # Get the mean of the rows without the center pixel (column 5 has the green band)
          min.meanG.vals <- mean(min.pt.vals[, 5][min.pt.vals$cells != center.cell])
          max.meanG.vals <- mean(max.pt.vals[, 5][max.pt.vals$cells != center.cell])
          
          
          ######  Mean brightness
          # Get the mean of the pixel
          min.meanB.val <- mean(as.matrix(min.pt.vals[, c(3:9)][min.pt.vals$cells == center.cell, ]))
          max.meanB.val <- mean(as.matrix(max.pt.vals[, c(3:9)][max.pt.vals$cells == center.cell, ]))
          
          # Get the mean of all the other pixels
          min.meanB.vals <- mean(as.matrix(min.pt.vals[, c(3:9)][min.pt.vals$cells != center.cell, ]))
          max.meanB.vals <- mean(as.matrix(max.pt.vals[, c(3:9)][max.pt.vals$cells != center.cell, ]))
          
          
          ######  NDVI ratio
          # Calc the NDVIs
          minRed <- min.pt.vals[, c(6)][min.pt.vals$cells == center.cell]
          minNIR <- min.pt.vals[, c(7)][min.pt.vals$cells == center.cell]
          maxRed <- max.pt.vals[, c(6)][max.pt.vals$cells == center.cell]
          maxNIR <- max.pt.vals[, c(7)][max.pt.vals$cells == center.cell]
          minReds <- min.pt.vals[, c(6)][min.pt.vals$cells != center.cell]
          minNIRs <- min.pt.vals[, c(7)][min.pt.vals$cells != center.cell]
          maxReds <- max.pt.vals[, c(6)][max.pt.vals$cells != center.cell]
          maxNIRs <- max.pt.vals[, c(7)][max.pt.vals$cells != center.cell]
          
          # Calc the NDVIs
          minNDVI <- (minNIR - minRed) / (minNIR + minRed)
          maxNDVI <- (maxNIR - maxRed) / (maxNIR + maxRed)
          minNDVIs <- (minNIRs - minReds) / (minNIRs + minReds)
          maxNDVIs <- (maxNIRs - maxReds) / (maxNIRs + maxReds)
          
          
          # Plug in the ratios
          if(x == rad.list[1]){
            final.df$minGreenContext50[final.df$cellnumber == center.cell] <- min.green / min.meanG.vals
            final.df$maxGreenContext50[final.df$cellnumber == center.cell] <- max.green / max.meanG.vals
            final.df$minBrightContext50[final.df$cellnumber == center.cell] <- min.meanB.val / min.meanB.vals
            final.df$maxBrightContext50[final.df$cellnumber == center.cell] <- max.meanB.val / max.meanB.vals
            final.df$minNDVIContext50[final.df$cellnumber == center.cell] <- minNDVI / mean(minNDVIs)
            final.df$maxNDVIContext50[final.df$cellnumber == center.cell] <- maxNDVI / mean(maxNDVIs)
          }else if(x == rad.list[2]){
            final.df$minGreenContext100[final.df$cellnumber == center.cell] <- min.green / min.meanG.vals
            final.df$maxGreenContext100[final.df$cellnumber == center.cell] <- max.green / max.meanG.vals
            final.df$minBrightContext100[final.df$cellnumber == center.cell] <- min.meanB.val / min.meanB.vals
            final.df$maxBrightContext100[final.df$cellnumber == center.cell] <- max.meanB.val / max.meanB.vals
            final.df$minNDVIContext100[final.df$cellnumber == center.cell] <- minNDVI / mean(minNDVIs)
            final.df$maxNDVIContext100[final.df$cellnumber == center.cell] <- maxNDVI / mean(maxNDVIs)
          }else if(x == rad.list[3]){
            final.df$minGreenContext200[final.df$cellnumber == center.cell] <- min.green / min.meanG.vals
            final.df$maxGreenContext200[final.df$cellnumber == center.cell] <- max.green / max.meanG.vals
            final.df$minBrightContext200[final.df$cellnumber == center.cell] <- min.meanB.val / min.meanB.vals
            final.df$maxBrightContext200[final.df$cellnumber == center.cell] <- max.meanB.val / max.meanB.vals
            final.df$minNDVIContext200[final.df$cellnumber == center.cell] <- minNDVI / mean(minNDVIs)
            final.df$maxNDVIContext200[final.df$cellnumber == center.cell] <- maxNDVI / mean(maxNDVIs)
          }else if(x == rad.list[4]){
            final.df$minGreenContext500[final.df$cellnumber == center.cell] <- min.green / min.meanG.vals
            final.df$maxGreenContext500[final.df$cellnumber == center.cell] <- max.green / max.meanG.vals
            final.df$minBrightContext500[final.df$cellnumber == center.cell] <- min.meanB.val / min.meanB.vals
            final.df$maxBrightContext500[final.df$cellnumber == center.cell] <- max.meanB.val / max.meanB.vals
            final.df$minNDVIContext500[final.df$cellnumber == center.cell] <- minNDVI / mean(minNDVIs)
            final.df$maxNDVIContext500[final.df$cellnumber == center.cell] <- maxNDVI / mean(maxNDVIs)
          }else{
            final.df$minGreenContext1000[final.df$cellnumber == center.cell] <- min.green / min.meanG.vals
            final.df$maxGreenContext1000[final.df$cellnumber == center.cell] <- max.green / max.meanG.vals
            final.df$minBrightContext1000[final.df$cellnumber == center.cell] <- min.meanB.val / min.meanB.vals
            final.df$maxBrightContext1000[final.df$cellnumber == center.cell] <- max.meanB.val / max.meanB.vals
            final.df$minNDVIContext1000[final.df$cellnumber == center.cell] <- minNDVI / mean(minNDVIs)
            final.df$maxNDVIContext1000[final.df$cellnumber == center.cell] <- maxNDVI / mean(maxNDVIs)
          }
        }
      }
    }
  }
  
  # Write the file
  write.csv(final.df, "E:/Dropbox/Permanent/Grad School/Projects/EleTree/analysis/Landsat_v_woody_cover_regressions/TwelvePark_Landsat_n_GT_values_w_SpatContext.csv", row.names = FALSE)
  
}
