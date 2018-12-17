
# Extract MESMA values to the GT table for all the parks 

# Ryan Nagelkirk
# 2018/08/20

# Description
# This code was written to go through the first round of the 2016 MESMA results and compare them to ground truth.
# Outputs table and plots.




#############################  LOAD LIBRARIES  ######################################################################
{
  library(raster)
  library(rgdal)
  library(jpeg)
  library(ggmap)
  library(ggplot2)
  library(RgoogleMaps)
  library(png)
  library(maps)
  library(fields)
  library(sp)
  library(dplyr)
}


#################################  Functions  ##################################################################
# Source the functions
source("E:/Dropbox/Permanent/Grad School/src_functions/src_masterfunctions.R")
# source("~/Dropbox/Permanent/Grad School/src_functions/src_masterfunctions.R")

# Change the memory limit to handle the larger files 
memory.limit(size = 50000) # default: 16308

#############################  SET WORKING DIRECTORY  ######################################################################
# Set working directory to where the results will be saved
setwd("E:/Dropbox/Permanent/Grad School/Projects/EleTree/pubs/MethodsPaper/analysis/accuracy_assessments")
# setwd("~/Dropbox/Permanent/Grad School/Projects/EleTree/z_pub_archives/MethodsPaper/analysis/accuracy_assessments")

#################################  Variables  #################################################################
park.names <- c("Chobe", "Kruger", "Limpopo", "Mpala", "Murchison", "North_Luangwa", "QWE", "Ruaha", "Selous", "Serengeti", "South_Luangwa")
                

                
##############  Create table of Ground truth vs MESMA for each MESMA image and each park  #########################
for(i in 1:length(park.names)){
  {
    # Start for loop here
    park.name <- park.names[i]
    
    ###############################  DIRECTORIES  ##############################################################
    # Folder with shape files of ground truth area
    shape.folder <- paste0("E:/Dropbox/Permanent/Grad School/Projects/EleTree/data/ParkData/", park.name, "/Boundary")
    # shape.folder <- paste0("~/Dropbox/Permanent/Grad School/Projects/EleTree/data/ParkData/", park.name, "/Boundary")
    
    # Folder with ground truth data
    ground.folder <- paste0("E:/Dropbox/Permanent/Grad School/Projects/EleTree/data/ParkData/", park.name, "/ground_truth")
    # ground.folder <- paste0("~/Dropbox/Permanent/Grad School/Projects/EleTree/data/ParkData/", park.name, "/ground_truth")
    
    # Folder with SMA results (Auto plot: just set experiment.folder)
    SMA.folder <- paste0("X:/nagelki4/Projects/EleTree/data/ParkData/", park.name, "/MESMA/hdr_results/Methods_Paper_Results_Year2016")
    # SMA.folder <- paste0("~/Dropbox/Permanent/Grad School/Projects/EleTree/data/ParkData/", park.name, "/MESMA/hdr_results")
    
    #######################  Create List of the MESMAs for the Park  ##############################################
    # Get list of MESMA files
    hdr.list <- list.files(SMA.folder, pattern = ".hdr", full.names = T)
    hdr.list <- hdr.list[!grepl(pattern = "class", hdr.list)] # exclude the hdr files of the classifications
    resid.list <- hdr.list[grepl(pattern = "resids", hdr.list)] # make a list of the files with the residuals to plot as well
    hdr.list <- hdr.list[!grepl(pattern = "resids", hdr.list)] # make a list of the actual MESMA results (i.e. the unmixed percent covers)
    
    #######################  VCF file  ########################################################################
    # Get the name of the VCF file and add it to the hdr.list
    vcf.file <- list.files(paste0("X:/nagelki4/Projects/EleTree/data/ParkData/", park.name, "/VCF"), pattern = "clipped", full.names = T)
    # Add it to hdr.list
    hdr.list <- c(hdr.list, vcf.file)
    
    ############  Park Boundary  ##############################################################################
    # Get the boundary
    if(park.name == "Serengeti"){
      park.boundary <- "Serengeti_Mara_dissolve"
    }else if(park.name == "Kruger"){
      park.boundary <- "Kruger_Boundary"
    }else if(park.name == "Mpala"){
      park.boundary <- "Mpala_Boundary"
    }else{
      park.boundary <- park.name
    }
    
    # Read in the file
    park.border <- readOGR(shape.folder, park.boundary)
    plot(park.border)
    
    ######  GROUND TRUTH TABLE  #############################################################################
    # List the groundtruth files
    gt.list <- list.files(paste0(ground.folder, "/wholepark"), pattern = "groundtruth_doneByHand", full.names = T)
    # Get the most recent one
    most.recent.gt <- sort(gt.list, decreasing = T)[1]
    # Read in the groundtruth table
    ground.truth <- read.csv(most.recent.gt)
    
    # Get rid of any rows that have tree soil and grass all at zero
    ground.truth <- ground.truth[!apply(ground.truth[, c(4:6)]==0, 1, all),]
    
    # Get rid of an NA rows
    ground.truth <- ground.truth[!is.na(ground.truth$park.pxl.num), ]
    
    # Get rid of water.flag rows
    ground.truth <- ground.truth[ground.truth$water.flag == 0, ]
  }
  #############  Create GT shapefile - or- Load GT shapefile if created already ########################################
  # Here, the ground truth data needs to be translated to coordinates. This is done differently
  # for Kruger versus all the other parks, because Kruger's ground truth was collected using
  # an earlier form of the code (Logan's work). 
  
  
  # If the points were already created, load that file. Else, create the points.
  if(file.exists(paste0(ground.folder, "/wholepark/gt_points.shp"))){
    gt.points <- readOGR(paste0(ground.folder, "/wholepark"), "gt_points")
  }else{
    # For Kruger, clip the background image and use the coordinates from the groundtruth file
    # to do raster to points
    if(park.name == "Kruger" | park.name == "Mpala"){
      
      # Get backdrop tif file
      if(park.name == "Kruger"){
        backdrop <- raster(paste0(ground.folder, "/kruger_backdrop_mosaic.tif")) # needs to be this (not backdropImage file) because it's what it indexed from
      }else{
        backdrop <- raster(paste0(ground.folder, "/mpala_backdrop_mosaic.tif")) # needs to be this (not backdropImage file) because it's what it indexed from
      }
      
      # backdrop <- raster(paste0(ground.folder, "/Kruger_LC08_20161029_20160420_backdropImage.tif"))
      new.border <- spTransform(park.border, crs(backdrop))
      
      # If the table was already created, pull it in
      if(file.exists(paste0(ground.folder, "/wholepark/gt_w_UTM_coords.csv"))){
        gt.points <- read.csv(paste0(ground.folder, "/wholepark/gt_w_UTM_coords.csv"))
      }else{
        
        # Get rid of duplicates
        ground.truth <- ground.truth[!duplicated(ground.truth$park.pxl.num), ]
        
        # Get the cell numbers
        cellnumbers <- ground.truth$park.pxl.num
        
        # Clip the backdrop to the park boundary
        blank.park <- clipTIF(tifname = backdrop, clipboundary = new.border)
        plot(blank.park)
        plot(new.border, add = T)
        
        # First have to extract the points from the raster via this somewhat clumsy manipulation
        # Assign the sample points 1000000 and then set rest to NA
        blank.park[cellnumbers] <- 1000000
        blank.park[blank.park < 1000000] <- NA
        
        # Convert cells to points and add them to the plot
        gt.points <- rasterToPoints(blank.park)
        gt.df <- as.data.frame(gt.points[, c(1:2)])
        
        # Add the coordinates to the ground truth
        # Luckily, both coordinate systems (the pixel indexing and UTM) are ordered to read top down, left to right
        # That is only true for this quadrant of the UTM coord system. 
        # So order the pixel.num column.
        new.gt.df <- ground.truth[order(ground.truth$park.pxl.num), ]
        # Add the coordinates and save for later use
        gt.points <- cbind(gt.df, new.gt.df)
        write.csv(gt.points, paste0(ground.folder, "/wholepark/gt_w_UTM_coords.csv"), row.names = F)
      }
      
      # Transform to points
      coordinates(gt.points) <- c("x", "y")
      
      # First, get the backdrop image
      landsat.image <- list.files(ground.folder, pattern = "backdropImage.tif", full.names = T) # this needs to be a landsat image that covers the entire sample area. It serves as a template for point creation
      backdrop <- raster(landsat.image)
      
      # Reproject to the raster (the raster is in UTM. The points have UTM coords while the border just needs to be transformed)
      new.border <- spTransform(park.border, crs(backdrop))
      crs(gt.points) <- crs(backdrop)
      
      # Now plot them
      plot(new.border)
      points(gt.points, pch = 20, cex = .1)
      
      # Save the points as shapefile for use outside R. 
      writeOGR(obj = gt.points, dsn = paste0(ground.folder, "/wholepark"), layer = "gt_points", driver = "ESRI Shapefile", overwrite_layer = T)
    }
    
    
    
    # For the rest, do the same thing as above (create the gt.points) but using the coordinates in the GT names rather
    # than using the old pixel indexing approach. This method is more robust. 
    if(park.name != "Kruger" & park.name != "Mpala"){
      # Get a list of all the GT image names
      image.names <- list.files(paste0(ground.folder, "/wholepark/handImages"), pattern = "3030_classified")
      # Loop through the names, filling a df as you go with the park.pxl.num, UTM coords 
      for(j in 1:length(image.names)){
        newrow <- c(
          pxl = as.numeric(gsub('^.*_pxlnum_\\s*|\\s*_sample_.*$', '', image.names[j])),
          x = as.numeric(gsub('^.*_x_\\s*|\\s*_y_.*$', '', image.names[j])),
          y = as.numeric(gsub('^.*_y_\\s*|\\s*_MATCHED_.*$', '', image.names[j])))
        if(j == 1){
          gt.points <- newrow
        }else{
          gt.points <- rbind(gt.points, newrow)
        }
      }
      
      # Make into a df
      gt.points <- as.data.frame(gt.points)
      
      # ### TEMP
      # # Bring in the temp files
      # ground.truth <- read.csv("C:/Users/nagelki-4/Desktop/groundTruth.csv")
      # gt.points <- read.csv("C:/Users/nagelki-4/Desktop/gtPoints.csv")
      
      # Rename the groundtruth's park.pxl column to match gt.points's
      names(ground.truth)[which(names(ground.truth) == "park.pxl.num")] <- "pxl"
      
      # Merge based on the park.pxl.num (now called pxl)
      gt.points <- merge(ground.truth, gt.points, by = "pxl")
      
      # Transform to points
      coordinates(gt.points) <- c("x", "y")
      
      #############  REPROJECT and PLOT POINTS  ####################################################################################
      # First, get the backdrop image
      landsat.image <- list.files(ground.folder, pattern = "backdropImage.tif", full.names = T) # this needs to be a landsat image that covers the entire sample area. It serves as a template for point creation
      backdrop <- raster(landsat.image)
      
      # Reproject to the raster (the raster is in UTM. The points have UTM coords while the border just needs to be transformed)
      new.border <- spTransform(park.border, crs(backdrop))
      crs(gt.points) <- crs(backdrop)
      
      # Now plot them
      plot(new.border)
      points(gt.points, pch = 20, cex = .1)
      
      # Eliminate the points that either had water or couldn't be classified (info in the ground.truth df)
      gt.points <- gt.points[gt.points$water.flag != 1, ]
      
      # Save the points as shapefile for use outside R. 
      writeOGR(obj = gt.points, dsn = paste0(ground.folder, "/wholepark"), layer = "gt_points", driver = "ESRI Shapefile", overwrite_layer = T)
    }
  }
  
  
  
  

  
  # Now we have points, regardless of whether the park was Kruger or any of the others.
  # So now we need to:
  # 1. For each MESMA result, extract the point values
  # 2. For each MESMA result, extract the nine point average
  # 3. The steps above will yield a table per MESMA result. Here, merge tables that have the same MESMA settings 
  # 4. Produce the final table of MESMA settings and the results as r2 values etc.
  
  #########  Extract the point values  ########################################################
  # h <- 1
  list.of.empty.files <- c()
  for(h in 1:length(hdr.list)){
    t1 <- tick()
    # Get file name
    mesma.base <- basename(hdr.list[h])
    current.SMA.file <- substr(mesma.base, 0, nchar(mesma.base)-4) # Take off the .hdr, which isn't needed
    full.current.SMA.file <- substr(hdr.list[h], 0, nchar(hdr.list[h])-4)
    print(paste(h, "of", length(hdr.list), ":", current.SMA.file))
    
    # If the points were already extracted and saved, skip the file
    if(file.exists(paste0(getwd(), "/MESMA_vs_GT_ptBYpt/", current.SMA.file, "_vs_GT.csv"))){
      print("Skipped")
      next 
    }
   
    
    # Check if the file name corresponds to a file that is already known to be empty
    # These are primary the shadeNorm files, that are based on other mesma rasters, and have values == Inf or NaN,
    # so they aren't caught in the other file checks, which only check for values == 0
    if(length(list.of.empty.files) > 0){
      binary <- FALSE
      for(j in 1:length(list.of.empty.files)){ # Have to loop, because grepl can only take a single argument
        if(length(current.SMA.file[grepl(pattern = list.of.empty.files[j], current.SMA.file)]) > 0 ){
          binary <- TRUE
        } 
      }
      if(binary == T){
        print("Based on empty raster")
        next
      } 
    }
    
    
    ###########  Determine if TS or TG and Set Band Numbers for indexing on image  ################################
    if(grepl(pattern = "TG", current.SMA.file)){
      TG_only <- TRUE
    }else{TG_only <- FALSE}
    
    # Set the band order for cover types in the SMA
    if(TG_only == T){ # the code is being funny, so running grass as a duplicate of soil. Won't be plotted
      tr.bandnum <- 1 
      gr.bandnum <- 2
      so.bandnum <- 2 
    }else{
      tr.bandnum <- 1 
      gr.bandnum <- 2  
      so.bandnum <- 3 
    }
    
    ############  Determine if it is a VCF file  ##################################################################
    if(grepl(pattern = "VCF", current.SMA.file)){
      VCF_file <- TRUE
    }else{VCF_file <- FALSE}
    
    ###########  Extract the points  ##############################################################################
    # Tifs need the full name with ".tif". MESMA loads without ".hdr", which is why there has to be an
    # if statement here
    if(VCF_file){
      current.mesma.raster <- raster(hdr.list[h])
      # Reported as 100*percent, so divide by 100
      current.mesma.raster <- current.mesma.raster/100
    }else{
      # Load the MESMA file
      current.mesma.raster <- brick(full.current.SMA.file)
    }
    
    
    # The points need to be transformed to the projection of the MESMA raster (not sure why, but had to)
    gt.points <- spTransform(gt.points, crs(current.mesma.raster))
    
    # # Make sure points overlap with MESMA raster
    # plot(current.mesma.raster[[tr.bandnum]])
    # points(gt.points, pch = 20, cex = .1)

    # Extract the point values
    gt.points.w.MESMA <- data.frame(coordinates(gt.points), gt.points, extract(current.mesma.raster[[tr.bandnum]], gt.points, method = 'simple'))
    names(gt.points.w.MESMA)[ncol(gt.points.w.MESMA)] <- 'MESMA.30'
    
    # Extract the 9-point mean by using a buffer that reaches from the center point out to the centers of the 8 surrounding pixels
    gt.points.w.MESMA.and.Nine.pt <- data.frame(coordinates(gt.points), gt.points.w.MESMA, 
                                                extract(current.mesma.raster[[tr.bandnum]], gt.points, method = 'simple', 
                                                        buffer = 50, fun = mean)) # I tested this, and it does just get the 9 pixels
    
    names(gt.points.w.MESMA.and.Nine.pt)[ncol(gt.points.w.MESMA.and.Nine.pt)] <- "MESMA.90"

    # Get rid of the excess columns 
    MESMA.gt.df <- subset(gt.points.w.MESMA.and.Nine.pt, select = -c(x.2, y.2, optional, x.1, y.1))
    
    # Write the table of MESMA and groundtruth values
    write.csv(MESMA.gt.df, paste0(getwd(), "/MESMA_vs_GT_ptBYpt/", current.SMA.file, "_vs_GT.csv"), row.names = F)
    tock(t1)
  }
  


  
  
  
  ############  Combine and save tables from same MESMA configurations but different Landsat Path  #################################
  # The steps above will yield a table of point values per MESMA result. Here, merge tables that have the same MESMA settings
  # taking the means when points overlap
  
  # Get the file names
  xls.files <- list.files(paste0(getwd(), "/MESMA_vs_GT_ptBYpt"), pattern = park.name)
  # Get rid of any that are "stitched" files (I didn't write this function)
  if(length(grep("stchd", xls.files, value = F)) > 0){
    xls.files <- xls.files[-grep("stchd", xls.files, value = F)]
  }
  
  # If some match methods but only differ in path, combine those
  methods <- gsub('^.*_\\d{8}_\\d{8}_*', '', xls.files) # the "d{}" calls out digits and they don't have to be defined. So it's just saying 8 digits
  # Get rid of the path name 
  methods.no.date <- gsub('*_\\d{8}_*', '_', methods) # This line specifically was added for Serengeti, but should work fine for the others 
  methods.no.path <- gsub('*_\\d{3}_*', '_', methods.no.date) # This doesn't work perfectly, but it should work (it deletes three digits from the date when it shouldn't)
  
  ##  If there are duplicated methods files, combine their MESMA point values taking the mean
  ##  Then save that new file to be read in in the next step
  if(max(duplicated(methods.no.path)) == 1){
    
    # Start an empty list of files to loop through later 
    complete.MESMAs <- c()
  
    # Create a list of matched files
    xls.files.dup <- unique(methods.no.path)
    
    # If the VCF file was listed in there, get it and take it out
    vcf.takenOut <- xls.files.dup[grepl("VCF", xls.files.dup)]
    xls.files.dup <- xls.files.dup[-grep("VCF", xls.files.dup)]
    
    # For each file, find the duplicates, load and cbind them
    for(v in 1:length(xls.files.dup)){
      # Get the names of the files with matching methods
      dups <- xls.files[which(methods.no.path == xls.files.dup[v])]
      # For each file, read it in and paste columns on as necessary
      for(b in 1:length(dups)){
        fl <- read.csv(paste0(getwd(), "/MESMA_vs_GT_ptBYpt/", dups[b]))
        if(b == 1){
          combo.file <- fl
        }else{
          new.cols <- subset(fl, select = c(MESMA.30, MESMA.90))
          names(new.cols) <- paste0(names(new.cols), ".", b)
          combo.file <- cbind(combo.file, new.cols)
        }
      }
      # Now take the averages and plug them into a new file
      # Get just the 30 m columns
      thirties <- combo.file[, grep("MESMA.30", names(combo.file))]
      nineties <- combo.file[, grep("MESMA.90", names(combo.file))]
      # Set zeros to NA to allow the mean calculations
      thirties[thirties == 0] <- NA
      nineties[nineties == 0] <- NA
      
      fl$MESMA.30 <- rowMeans(thirties, na.rm = T)
      fl$MESMA.90 <- rowMeans(nineties, na.rm = T)
      
      # Save the new file as an *stchd.xls
      combined.csv.name <- paste0(substr(dups[1], 1, nchar(dups[1]) - 4), "_stchd.csv")
      write.csv(fl, paste0(getwd(), "/MESMA_vs_GT_ptBYpt/", combined.csv.name), row.names = F)
      # Add the name of that file to a list of files to loop through
      complete.MESMAs <- c(complete.MESMAs, combined.csv.name)
    }
    # Add the VCF back
    complete.MESMAs <- c(complete.MESMAs, vcf.takenOut)
  }else{
    complete.MESMAs <- xls.files
  }
  
  
  
  ########################################################################################################
  #  Stitch all the tables together 
  ########################################################################################################
  
  # Now stitch them all into one big table of all 12 parks
  for(c in 1:length(complete.MESMAs)){
    
    # Read in the table
    tbl <- read.csv(paste0(getwd(), "/MESMA_vs_GT_ptBYpt/", complete.MESMAs[c]))
    # Insert the MESMA raster name to the front of the tbl(remove "_vs_GT.csv")
    tbl <- cbind(substr(complete.MESMAs[c], 1, nchar(complete.MESMAs[c]) - 10), tbl)
    names(tbl)[1] <- "Classification"
    
    # Kruger and Mpala don't have the pxl, 90 and 180 columns, nor "image.not.available, they also have extra columns, "X" and "mean"s
    # However, Mpala has some things that Kruger doesn't or vice versa, so have to have separate if statements 
    if(max(grepl("mean.", names(tbl)))) tbl <- tbl[, -grep("mean.", names(tbl))]
    if("X" %in% names(tbl) && "hand.tree" %in% names(tbl)){
      tbl <- tbl[, -grep("X", names(tbl))]
      names(tbl)[grep("hand.", names(tbl))] <- c("hand.tree.30", "hand.grass.30", "hand.soil.30")
      tbl <- cbind(tbl[, 1:3],
                  NA,
                  tbl[, 4:grep("truth.maj", names(tbl))],
                  NA, NA, NA, NA, NA, NA, NA,
                  tbl[, (grep("truth.maj", names(tbl))+1):ncol(tbl)])
      # Name the new columns 
      names(tbl)[grep("NA", names(tbl))] <- c("pxl", "hand.tree.90", "hand.grass.90", "hand.soil.90", "hand.tree.180", "hand.grass.180", "hand.soil.180", "image.not.available")
    }
    
    # The VCF files, for some reason, don't have a park.pxl.num column, so add it 
    if(!"park.pxl.num" %in% names(tbl)) tbl <- cbind(tbl[, 1:5], park.pxl.num = NA, tbl[, 6:ncol(tbl)])
    
    # Now append the table to the growing table 
    if(i == 1 && c == 1){
      complete.table <- tbl
    }else{
      complete.table <- rbind(complete.table, tbl)
    }
  }
}
  
# Write the file 
write.csv(complete.table, paste0(getwd(), "/all_MESMAs_appended_to_GT.csv"), row.names = F)

