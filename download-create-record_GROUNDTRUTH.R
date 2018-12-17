#################################  CHECKING GROUND TRUTH BY HAND  ##############################################
# Author: Ryan Nagelkirk
# Date: 5/15/2017

# Description:
# This script is specifically written for the manual classification of Google imagery. The classification serves as
# groundtruth for evaluating spectral unmixing results of Landsat imagery. For that reason, the Google imagery is 
# classified within the bounds of Landsat pixels (30x30 meter areas). The script helps the person classifying the imagery
# by automatically creating folders, selecting the sampling locations, downloading the Google imagery for those locations, 
# displaying classification results as you go, and then saving those results to predefined locations and file types. 
# The person classifying the imagery needs to define folder paths and variables at the beginning, and then sets the 
# classification thresholds in the "Classify Landcover" section. The classifier (you) also needs to manually change
# the sample number as they work through the sample points ("image.num" variable below "START of MANUAL CLASSIFICATION"
# section header). 
# NOTE: Unlike most scripts, this shouldn't be ran start to finish in one go. Instead, it aids the classifier in 
#       their work, accomplishing minor tasks at different points in the workflow as the classifier works through that workflow. 


#################################  Before You Begin  ###############################################################
# 1. Create a "ground_truth" folder within the park folder. 
# 2. Save a landsat image, or images, that cover the different areas that are going to be classified. For a small
#    park, this will generally be a single image that covers the entire park. For larger parks, several images might
#    be needed. It doesn't matter what band is used because this image is only used for indexing cells. 
#    What does matter is whether the image has the same shape and extent as all the other corresponding path/row images. 
#    This can be checked in ArcGIS. Ex: if the scene you're using is from path/row 168/078, then the image you chose 
#    should perfectly overlay other images of 168/078. Sometimes this is not the case and those images need to be discarded.
# 3. Once 1-2 are done, you're ready to go! The code itself is heavily annotated to allow you to work through it from here.
#    If there is not annotation, just run those lines - it isn't necessary for you to do anything else with them. 


#################################  Load Libraries  #############################################################
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
  library(R.utils)
}


{
  #################################  Functions  ##################################################################
  # Source the functions
  # source("E:/Dropbox/Permanent/Grad School/src_functions/src_masterfunctions.R")
  source("~/Dropbox/Permanent/Grad School/src_functions/src_masterfunctions.R")
  
  # Change the memory limit to handle the larger files 
  memory.limit(size = 50000) # default: 16308
  
  
  #################################  Variables  ##################################################################
  # Set the park
  park <- "Kruger"  # This is the folder name inside the working directory with all the files for the park
  
  
  #################################  Control Panel  ##############################################################
  # These control various parts of the code and whether they execute. 
  # If all are set to false, then nothing gets recorded or saved if the entire script is ran. That is good for when just trying to learn the code.
  # Once you have the hang of things, set them to true. In general, though, running the entire script in one go is not recommended. 
  # It is better to run only sections of this script at a time, since the classification process is manual and the script 
  # is just here to make classification streamlined, not entirely automated. 
  get.Google.image <- T # This controls whether the google images are downloaded. That process should only be done onces for each sample area (ref: "Get Google Imagery")
  save.rgb.pixel <- T  # Should the cropped 30x30 m RGB of the Google imagery be saved? Again, only needs to be done once. (ref: "Crop to 30x30 Meters")
  save.class.pixel <- T  # Should the classified image be saved to the image folder? Only save images once you have the classification you want. (ref: "Combine Covers")
  record.class.values <- TRUE # Should the classified veg cover values be plugged into the dataframe? (ref: "Record Values" section)
  save.final.class.df <- TRUE # Should the dataframe be saved? This should be done only once all the different images have been classified. (ref: "Save DF")
  
  
  
  # Classifier note: From here all the way down until the "START OF MANUAL CLASSIFICATION" section can be ran without any input.
  # Just make sure "get.Google.image" and "save.rgb.pixel" are set to false if the imagery is already downloaded. Otherwise it 
  # will download again and you'll have to wait for that. Also, the imagery might have changed and you'll have to do the entire 
  # classification all over if the old images are overwritten.
  
  # These variables will need to be changed as you move to new sample areas and/or parks.
  pixel.size.meters <- 30  # this shouldn't change. It is the length of one side of the area you want to classify. Typically will be set to landsat imagery resolution 
  larger.pxl.size <- 90  # this is the size of the image that I'll be using to classify the image
  truth.num <- 600  # the total number of points to take from the park.
  hand.truth.num <- 100  # the number of points to classify by hand. It can't go higher than truth.num
  min.hand.point.density <- 200 # This is the min density of sampling points (1 point every 200 km^2)
  hand.truth.image.folder <- "handImages" # this should stay the same. It is the folder within the park folder that will hold the hand-classified images
  ground.truth.image.folder <- "extraImages"  # this should also stay the same. Like above, except all sample images are saved here, not just the hand-classified
  
  
  # Get the working folder sorted out
  # working.dir <- "E:/Dropbox/Permanent/Grad School/Projects/EleTree/data/ParkData/" # set this to the folder path that has all the park folders in it. Make sure it has "/" at the end.
  working.dir <- "~/Dropbox/Permanent/Grad School/Projects/EleTree/data/ParkData/" # set this to the folder path that has all the park folders in it. Make sure it has "/" at the end.
  
  # List the park folder names
  folder.names <- list.dirs(working.dir, full.names = F, recursive = F)
  # Get the folder for this park
  folder.num <- which(substr(folder.names, 1, 3) %in% substr(park, 1, 3))
  # Set working directory
  setwd(paste0(working.dir, folder.names[folder.num], "/ground_truth"))
}


##  STOP HERE IF IMAGES HAVE ALREADY BEEN DOWNLOADED




# Get the tif file
landsat.image <- list.files(".", pattern = "backdropImage.tif") # this needs to be a landsat image that covers the entire sample area. It serves as a template for point creation

# Get the boundary
if(park == "Serengeti_Mara"){
  park.boundary <- "Serengeti_Mara_dissolve"
}else if(park == "Kruger"){
  park.boundary <- "Kruger_Boundary"
}else{
  park.boundary <- park
}

# par(mfrow = c(1,1)) # just in case needed for changing later
par(mar=c(1,1,2,1))



#############################  Create Points  #################################################################
###### Read in Landsat and shapefiles
# Landsat
landsat.image <- raster(landsat.image)
landsat.image[values(landsat.image) == 0] <- NA
# plot(landsat.image)
# Shapefiles
park.border <- readOGR("../Boundary", park.boundary)
# Convert shapefile to lat long, which everything else is in
park.border <- spTransform(park.border, crs(landsat.image))
plot(park.border)


##########################  Figure out the number of points needed  ###########################################
# Get the total area of the images in the park
tot.img.area <- area(park.border)/1e6
# Calculate the point density given the area 
hand.point.density <- tot.img.area / hand.truth.num

# If the density is lower than the minimum I've set, then calculate how many points are needed to meet that density requirement
if(hand.point.density < min.hand.point.density){
  hand.truth.num <- hand.truth.num
}else{
  hand.truth.num <- ceiling(tot.img.area / min.hand.point.density)
}


############################  Create DF  ######################################################################
# DF for manual classification entries
mx <- matrix(0, nrow = hand.truth.num, ncol = 13)
c.names <- c("Count", "park.pxl.num", "sample.pxl.num", "hand.tree.30", "hand.grass.30", "hand.soil.30", "truth.maj.cover", "hand.tree.90", "hand.grass.90", "hand.soil.90", "hand.tree.180", "hand.grass.180", "hand.soil.180")
colnames(mx) <- c.names # above, pxl.num is the pixel number within the sample area, park.pxl.num is from the entire park
hand.truth.df <- as.data.frame(mx)


#############################  Point Generation  ###############################################################


# Create background Park image
if(length(list.files(".", pattern = "blankpark")) < 1){
  blank.park <- clipTIF(tifname = landsat.image, clipboundary = park.border)
  blank.park[!is.na(blank.park)] <- 1
  
  # Save the raster
  writeRaster(blank.park, "blankpark.tif")
}else blank.park <- raster("blankpark.tif")

# Take a look
plot(blank.park, col = "lightgray", axes=F,box = F)
plot(park.border, add = T)


####  Generate List of Points 
# Create list of cells in Mpala for sampling
# First, how many cells have a value?
cell.vals <- length(na.omit(blank.park[])) # this used to be na.omit(mpala[])

# Use that to assign new values to the cells (numbers them)
blank.park[!is.na(blank.park)] <- 1:cell.vals
plot(blank.park, axes=F,box = F)

# Randomly select cells from the raster (the cells are represented as a list of numbers at this point)
set.seed(1) # set this for now. Can change later
sample.cells <- sample(1:cell.vals, truth.num)

####  Plot the Points
# First, get the cell numbers
# This is done because the entire raster is bigger than just the sample area. Outside the area,
# there are cells with NA values that are white in the plot. So have to find the cell number in the
# whole sample area raster, not just in the area that has values. 
sample.cell.numbers <- c()
for(i in sample.cells){
  sing.cell <- Which(blank.park == i, cells= TRUE) # the blank.park values aren't actual cells, but just cells values
  sample.cell.numbers <- c(sample.cell.numbers, sing.cell)
}

# Switch the sample area to NA and then give the sample points their pixel number value
trial <- blank.park
trial[] <- NA
trial[sample.cell.numbers] <- sample.cell.numbers # assign values for the cells based on their pixel number
plot(trial)

# Convert cells to points and add them to the plot
hand.points <- rasterToPoints(trial) # make the ground truth hand classified points - this is when they go out of order
par(mar=c(1,1,2,1))
plot(blank.park, axes=F,box = F, legend = FALSE, col = "lightgray")
plot(park.border, add = TRUE)
points(hand.points, pch = 0, cex = .3)



#############################  Park/Sample Area Cell Number Redundancy  ##############################################
#####  Get the cell numbers for the raster of the whole park (will be used in naming and added to hand.check.df)
#####  This is done for redundancy - so two different types of indexing are possible. So if park boundary shp file ever changes, 
#####  can just the sample area shp file to index. 
# Extract cell numbers
park.cell.nums <- extract(x = blank.park, y = hand.points[, 1:2], cellnumbers = TRUE) # this is what was wanted
park.cell.nums[, 2] <- hand.points[, 3] # reattach the pixel numbers from the sample area
# Set those points to high value and set all else to NA
park.ras <- blank.park
park.ras[] <- NA
park.ras[park.cell.nums[, 1]] <- park.cell.nums[, 1]

# Convert the remaining cells to points
park.points <- rasterToPoints(park.ras) # these points will match the hand.points above (they do)
# Compare the points
plot(hand.points, pch = 0, cex = .1)
points(park.points, pch = 0, cex = .1) # nothing should change when these points are added (nothing does)

## NOTE: If you plot using the cell numbers, they will be slightly offset. I checked this rigorously, and it's something
# to do with the plot function/display, not the actual data. Using the cell numbers yields the same coordinates whether it's
# from the whole park or from the sample area (using the correct cell numbers for each, of course) and using either 
# type of cell number (park- or sample area-derived) results in the very same image being pulled from Google. So there should 
# be no problem using either the park or sample area cell numbers when comparing to MESMA results, so long as, again, the 
# MESMA data has been clipped to each's value-derived areas (park boundary or sample area boundary, respectively)


#############################  Restore Cell Number Order  #########################################################

# First, get park.cell.nums back in the original order of sample.cell.numbers so we can use it and the first 50
# will still be distributed randomly
# First make into a df
cellnums.df <- as.data.frame(hand.points)
names(cellnums.df) <- c("x", "y", "pxl_num")
head(cellnums.df)


# Order to match sample.cell.numbers ordering
cellnums.df <- cellnums.df[match(sample.cell.numbers, cellnums.df$pxl_num),]
rownames(cellnums.df) <- 1:nrow(cellnums.df) # order the row names so they aren't confusing
# Add back the pxl numbers for just the park (when clipped, ie the number of the pixel inside the park when not counting NA pixels)
cellnums.df$sample_pxl <- sample.cells


#############################  Create Image Name Lists  ###########################################################
# Create lists of the files to crop (auto and hand)
auto.images <- c()
for(t in 1:truth.num){
  # Create the image name for that pixel
  image.name <- paste0("./wholepark/", ground.truth.image.folder, "/", park, "_site", t, "_", "pxlnum_", cellnums.df$pxl_num[t], 
                       "_sample_", cellnums.df$sample_pxl[t], "_x_", cellnums.df$x[t], "_y_", cellnums.df$y[t],  ".png")
  auto.images <- c(auto.images, image.name)
}

# Create hand image list
hand.images <- c()
for(t in 1:hand.truth.num){
  image.name.new <- paste0("./wholepark/", hand.truth.image.folder,  "/", park, "_site", t, "_", "pxlnum_", cellnums.df$pxl_num[t], 
                           "_sample_", cellnums.df$sample_pxl[t], "_x_", cellnums.df$x[t], "_y_", cellnums.df$y[t],".png")
  hand.images <- c(hand.images, image.name.new)
}


# Create the google image folders if not already done
if(file.exists(paste0("./wholepark/", ground.truth.image.folder)) == F){
  dir.create(paste0("./wholepark/", ground.truth.image.folder), recursive = T)
}

if(file.exists(paste0("./wholepark/", hand.truth.image.folder)) == F){
  dir.create(paste0("./wholepark/", hand.truth.image.folder), recursive = T)
}


###################### Plot and export point images  ####################################################################
# First, get a list of the cell numbers by pulling them from the names 
hand.cells <- c()
for(q in 1:length(hand.images)){
  cellNum <- as.numeric(gsub('^.*_pxlnum_\\s*|\\s*_sample_.*$', '', hand.images[q]))
  # Add to the list
  hand.cells <- c(hand.cells, cellNum)
}
  
# Switch the sample area to NA and then give the sample points their pixel number value
trial <- blank.park
trial[] <- NA
trial[hand.cells] <- hand.cells # assign values for the cells based on their pixel number
plot(trial)

# Convert cells to points and add them to the plot
hand.points.subset <- rasterToPoints(trial) # make the ground truth hand classified points - this is when they go out of order
par(mar=c(1,1,2,1))
# Plot the entire set as well, to make sure they match
plot(park.border)
points(hand.points, pch = 0, cex = .3)
# Now plot the subset
plot(park.border)
points(hand.points.subset, pch = 0, cex = .3)

# Now export them
{png("extraImages.png") # Save the image
plot(blank.park, axes=F,box = F, legend = FALSE, col = "lightgray")
plot(park.border, add = TRUE)
points(hand.points, pch = 0, cex = .3)
dev.off()}

{png("handImages.png") # Save the image
plot(blank.park, axes=F,box = F, legend = FALSE, col = "lightgray")
plot(park.border, add = TRUE)
points(hand.points.subset, pch = 0, cex = .3)
dev.off()}



#############################  Get Google Imagery  #################################################################
# This is the section of code that actually downloads and saves the google imagery
if(get.Google.image){
  # Create a blank list of the list of new autoimage names
  new_autoimages <- c()
  
  
  # Get the Google image for each pixel
  for(t in 1:length(auto.images)){
    # Create the image name for that pixel
    image.name <- auto.images[t] # this looks risky, but the naming above follows this same order
    
    # Check if the file already exists. If it does, go to the next
    google.dwnlds <- list.files("./wholepark/extraImages", full.names = T)
    google.dwnlds <- gsub('\\s*_x_.*$', '', google.dwnlds)
    short.img.name <- gsub('\\s*_x_.*$', '', image.name)
    if(short.img.name %in% google.dwnlds){
      next
    }
    
      
    # GET THE IMAGE FROM GOOGLE!! 
    # Function in src_masterfunctions.R (it returns the portion of the name added in the function itself)
    appendix <- getGoogleimage(rastername = blank.park, ras.pixel.num = cellnums.df$pxl_num[t],
                        destfilename = image.name, typeofmap = "satellite", zoomlevel = 19)

    # Function in src_masterfunctions.R (it returns the portion of the name added in the function itself)
    # appendix <- getGoogleimage(rastername = blank.park, ras.pixel.num = cellnums.df$pxl_num[t],
    #                            destfilename = paste0(substr(image.name, 1, nchar(image.name) - 15), "_out.png"), typeofmap = "satellite", zoomlevel = 16)

    
    
    # j has periods in it that need to be removed, otherwise they confuse the readpng function farther below
    appendix_new <- gsub("\\.", "pnt", substr(appendix, 1, nchar(appendix) - 4))
    appendix_new <- paste0(appendix_new, ".png")
    
    # It won't work in the function (I think the name becomes too long), so append to the name out here
    # List and delete the rda file
    file.remove(list.files("./wholepark/extraImages", pattern = ".rda", full.names = T))
    
    # List the file that was just created
    # flz <- list.files("./wholepark/extraImages", pattern = substr(basename(image.name), 1, nchar(basename(image.name)) - 4), full.names = T)
    
    # Rename with the appendage - you have to subset to the first name because the code sometimes picks up the rda file as well
    new.name <- paste0(substr(image.name, 1, nchar(image.name) - 4), appendix_new)
    file.rename(image.name, new.name)
    
    # # Add the name to a list of new autoimages
    # if(t <= hand.truth.num){
    new_autoimages <- c(new_autoimages, new.name)
    # }
  }
  
  # Move the first n (hand.truth.num) images to the hand classifying folder (they are now in random order (a good thing: it's the original order from random number generator))
  for(t in 1:length(hand.images)){
    image.name <- new_autoimages[t]
    image.name.new <- gsub("extraImages", "handImages", new_autoimages[t])
    file.copy(image.name, image.name.new, overwrite = T)
    print(paste0(t, ". ", image.name.new, " moved to hand folder: ", file.exists(image.name.new)))
  }
}

# Below, I want to classify and record the cover for a 90 by 90 and 30x30 area. 
# Will classify using the whole image, then crop it to get the 90 and 30 meter values
# Then both sets of values will be recorded

#############################  Create 30x30 & 90x90 Shapefile  ##########################################################
# This isn't currently used to crop the images because it was noticeably slower than the crop.google line
{
  cell.size <- pixel.size.meters/190 # 190m is the length of one side of the Google image
  larger.cell.size <- larger.pxl.size/190
  landscape.size <- 180/190 # I'm doing this to remove the text at the border in the images, which will register as trees
  
  half.pix <- cell.size / 2
  large.half.pix <- larger.cell.size / 2
  land.half.pix <- landscape.size / 2
  
  center <- .5 # Middle of the image. The reference point for indexing. Will always be the same.
  
  # Create the boundaries of the upper left cell. The 9-pixel loop will modify 
  # these as it counts through the cells (extent.list)
  bottom <- center - half.pix
  right <- center + half.pix
  top <- center + half.pix 
  left <- center - half.pix 
  l.bottom <- center - large.half.pix
  l.right <- center + large.half.pix
  l.top <- center + large.half.pix 
  l.left <- center - large.half.pix
  la.bottom <- center - land.half.pix
  la.right <- center + land.half.pix
  la.top <- center + land.half.pix 
  la.left <- center - land.half.pix
  
  #Create extent object 
  cr.ext <- extent(c(left, right, bottom, top))
  l.cr.ext <- extent(c(l.left, l.right, l.bottom, l.top))
  la.cr.ext <- extent(c(la.left, la.right, la.bottom, la.top))
  
  # Coerce extent to a SpatialPolygons object 
  cell.poly <- as(cr.ext, 'SpatialPolygons')
  l.cell.poly <- as(l.cr.ext, 'SpatialPolygons')
  landscape.cell.poly <- as(la.cr.ext, 'SpatialPolygons')
}

#############################  Crop to 30x30 Meters  ##########################################################
# List the hand image
hand.list <- gsub("extraImages", "handImages", new_autoimages[1:hand.truth.num])
# hand.list <- list.files("./wholepark/handImages", pattern = ".png", full.names = T)

# Create the hand image 30x30 rgb name list
rgb.3030.names <- c()
for(x in 1:length(hand.list)){
  rgb.3030.name <- paste0(substr(hand.list[x], 1, nchar(hand.list[x]) - 4), "_3030_rgb.png")
  rgb.3030.names <- c(rgb.3030.names, rgb.3030.name)
}

# Create the hand image 30x30 rgb name list
rgb.9090.names <- c()
for(x in 1:length(hand.list)){
  rgb.9090.name <- paste0(substr(hand.list[x], 1, nchar(hand.list[x]) - 4), "_9090_rgb.png")
  rgb.9090.names <- c(rgb.9090.names, rgb.9090.name)
}


## Loop through and crop the hand images. Later code can do the same for the auto.images
if(save.rgb.pixel){
  for(x in 1:length(hand.list)){
    # Crop the image (function in src_masterfunctions.R)
    img <- crop.google(image.name = hand.list[x], pixel.res = pixel.size.meters)
    img2 <- crop.google(image.name = hand.list[x], pixel.res = larger.pxl.size)
    
    # Plot and save the RGB image of the cell - not needed if they are already saved
    png(rgb.3030.names[x]) # Save the image
    plotRGB(img, scale = 1)
    dev.off()
    
    png(rgb.9090.names[x]) # Save the image
    plotRGB(img2, scale = 1)
    dev.off()
  }
}

#############################  Save DF  ##########################################################################
# Plug in the pxl numbers 
hand.truth.df$park.pxl.num <- cellnums.df$pxl_num[1:hand.truth.num]
hand.truth.df$sample.pxl.num <- cellnums.df$sample_pxl[1:hand.truth.num]
hand.truth.df$Count <- c(1:hand.truth.num)

if(save.final.class.df){
  write.csv(hand.truth.df, paste0("./wholepark/groundtruth_doneByHand_", YMD(), ".csv"), row.names = F)
}







#############################  START of MANUAL CLASSIFICATION  ###################################################
# Load the table
if(length(list.files("./wholepark", pattern = "groundtruth_doneByHand")) >= 1){
  num.saved.files <- length(list.files("./wholepark", pattern = "groundtruth_doneByHand"))
  saved.file.name <- list.files("./wholepark", pattern = "groundtruth_doneByHand")[num.saved.files]
  hand.truth.df <- read.csv(paste0("./wholepark/", saved.file.name))
  if(length(hand.truth.df$tree_equation) > 0){
    hand.truth.df$tree_equation <- as.character(hand.truth.df$tree_equation)
    hand.truth.df$grass_equation <- as.character(hand.truth.df$grass_equation)
    hand.truth.df$soil_equation <- as.character(hand.truth.df$soil_equation)
    hand.truth.df$truth.maj.cover <- as.character(hand.truth.df$truth.maj.cover)
  }
}


# If images already exist, create list names from those images, else create them from the variables already existing
if(length(list.files("./wholepark/handImages", pattern = "rgb")) >= 1){
  rgb90 <- list.files("./wholepark/handImages", pattern = "9090_rgb", full.names = T)
  rgb30 <- list.files("./wholepark/handImages", pattern = "3030_rgb", full.names = T)
  handp <- gsub("_9090_rgb", "", rgb90)
  # Get the pixel numbers, which will be in the correct order here
  sample.pxl.num <- as.numeric(hand.truth.df$sample.pxl.num)
  park.pxl.num <- as.numeric(hand.truth.df$park.pxl.num)
  
  # Get the unordered list
  unordered.sample <- as.numeric(gsub('^.*_sample_\\s*|\\s*_x_.*$', '', handp))
  unordered.park <- as.numeric(gsub('^.*_pxlnum_\\s*|\\s*_sample_.*$', '', handp))
  
  # Put the name lists in the correct order
  # Start the empty list
  hand.list <- c()
  rgb.9090.names <- c()
  rgb.3030.names <- c()
  
  for(x in 1:length(handp)){
    # Get the image name that matches botht he park and sample numbers
    position <- which(unordered.park == park.pxl.num[x] & unordered.sample == sample.pxl.num[x])
    hand.list <- c(hand.list, handp[position])
    rgb.9090.names <- c(rgb.9090.names, rgb90[position])
    rgb.3030.names <- c(rgb.3030.names, rgb30[position])
  }
} 




#############################  PLUG in Image Number  ##########################################################
image.num <- image.num + 1 # Manually set this as you go. It represents the point you are classifying. So will be in range of 1-50.
# After the first point is classified, set this to 2, then 3 after that, and so on, until you hit 50. After 50, you're done!

# 435 was the old number 

# This isn't currently used to crop the images because it was noticeably slower than the crop.google line
{
  cell.size <- pixel.size.meters/190 # 190m is the length of one side of the Google image
  larger.cell.size <- larger.pxl.size/190
  landscape.size <- 180/190 # I'm doing this to remove the text at the border in the images, which will register as trees
  
  half.pix <- cell.size / 2
  large.half.pix <- larger.cell.size / 2
  land.half.pix <- landscape.size / 2
  
  center <- .5 # Middle of the image. The reference point for indexing. Will always be the same.
  
  # Create the boundaries of the upper left cell. The 9-pixel loop will modify 
  # these as it counts through the cells (extent.list)
  bottom <- center - half.pix
  right <- center + half.pix
  top <- center + half.pix 
  left <- center - half.pix 
  l.bottom <- center - large.half.pix
  l.right <- center + large.half.pix
  l.top <- center + large.half.pix 
  l.left <- center - large.half.pix
  la.bottom <- center - land.half.pix
  la.right <- center + land.half.pix
  la.top <- center + land.half.pix 
  la.left <- center - land.half.pix
  
  #Create extent object 
  cr.ext <- extent(c(left, right, bottom, top))
  l.cr.ext <- extent(c(l.left, l.right, l.bottom, l.top))
  la.cr.ext <- extent(c(la.left, la.right, la.bottom, la.top))
  
  # Coerce extent to a SpatialPolygons object 
  cell.poly <- as(cr.ext, 'SpatialPolygons')
  l.cell.poly <- as(l.cr.ext, 'SpatialPolygons')
  landscape.cell.poly <- as(la.cr.ext, 'SpatialPolygons')
  
  
  
  #############################  Index Calculations  ###############################################################
  # Read in the single cell and landscape images
  current.png <- readPNG(hand.list[image.num])
  ninetyM.png <- readPNG(rgb.9090.names[image.num])
  thirtyM.png <- readPNG(rgb.3030.names[image.num])
  
  # Pull out the red, green and blue bands from both images
  red.c <- raster(current.png[,,1]) # the "c" is just a reminder that it's the color value for the cropped image
  green.c <- raster(current.png[,,2])
  blue.c <- raster(current.png[,,3])
  
  red90 <- raster(ninetyM.png[,,1]) 
  green90 <- raster(ninetyM.png[,,2])
  blue90 <- raster(ninetyM.png[,,3])
  
  red30 <- raster(thirtyM.png[,,1]) 
  green30 <- raster(thirtyM.png[,,2])
  blue30 <- raster(thirtyM.png[,,3])
  
  # Stack so can create images later
  rgb.stack <- stack(red.c, green.c, blue.c)
  ninety.stack <- stack(red90, green90, blue90)
  thirty.stack <- stack(red30, green30, blue30)
  
  # Sum the RGB values. Low values are dark areas which correspond to trees. 
  sum.im.pre <- green.c + red.c + blue.c
  
  # Calculate the Green-Red Vegetation Index
  grvd.pre <- (green.c - red.c) / (green.c + red.c)
  
  # Calculate summary stats for images to determine if it is mostly grass or soil
  #  This is done because images with a lot of grass or soil require different thresholds for classifying tree cover
  
  # Find the percent area of the cell that has a brightness sum over 2
  bright <- red30 + green30 + blue30
  bright.area <- bright > 2
  bright.frac <- length(bright.area[bright.area == 1]) / length(bright.area)
  
  # Look at pixels with higher red values than green
  red.area <- red30 > green30
  red.frac <- length(red.area[red.area == 1]) / length(red.area) 
  
  
  
  # Below, classifier is able to set thresholds as they like. At end, all thresholds along with every index and measure 
  # are recorded in the df.That way, can analyze the values (simple regression) to see if there are relationships between
  # the thresholds chosen and the different brightness values etc.
  
  
  
  #############################  Classify Landcover  ###############################################################
  # This is the section where you will have to iteratively change the thresholds to create the classification you like.
  # First get tree cover to where you want it, then move on to grass cover. Soil cover will be done automatically as the
  # remaining area. 
  
  # First, plot the image of the overall area with the cell overlaid
  {
    plotRGB(rgb.stack, scale = 1)
    plot(cell.poly, add = TRUE)
    plot(l.cell.poly, add = TRUE)
  }
  
  # Is their an image?
  image.not.available <- F
  # IMPORTANT: Was there water in the classified area? If yes, set this variable to "TRUE"
  water.was.present <- F
  
  # If there is no image, this will flag that. Then just manually move to the next image
  if(image.not.available == TRUE){
    hand.truth.df[image.num, 4:22] <- NA # put NA in all the columns that would have had entries
    hand.truth.df[image.num, 14] <- 1
  }else{hand.truth.df$image.not.available[image.num] <- 0}
}


# Do trees, then grass. Soil will automatically be calculated after grass (soil = remaining pixels)
# TREES 
tree.bright.thresh <- .62
{ # I did this to make it all run at once
  trees <- sum.im.pre < tree.bright.thresh
  tree.equation <- paste0("sum.im.pre < ", tree.bright.thresh)
  # Compare classification to image
  par(mar=c(1,1,2,1))
  plot(trees) # this is done twice to overlay the rgb image in the second plot
  plotRGB(rgb.stack, scale = 1, add = TRUE)
  plot(trees) # tree cover will be green
}

# Now do the same for the other two images 
{
  ninety <- crop.class.image(trees, 90) # try crop.google
  thirty <- crop.class.image(trees, 30)
  par(mar=c(1,1,2,1))
  
  plot(trees) # this is done twice to overlay the rgb image in the second plot
  plotRGB(ninety.stack, scale = 1, add = T)
  plot(ninety) # tree cover will be green
  plot(trees) # this is done twice to overlay the rgb image in the second plot
  plotRGB(thirty.stack, scale = 1, add = TRUE)
  plot(thirty) # tree cover will be green
}
 
  
# Now remove the tree pixels from the image that will go into identifying soil and grass
# create tree mask
{
  tree.mask <- trees
  tree.mask[tree.mask == 1] <- 50
  holy.image <- grvd.pre + tree.mask
  holy.image.3 <- sum.im.pre + tree.mask
  holy.image[holy.image > 5] <- NA
  holy.image.3[holy.image.3 > 5] <- NA
  grvd <- holy.image
  sum.im <- holy.image.3
}



##### GRASS
# Split soil and grass using GRVD and brightness. Soil is nonveg and bright
grass.bright.thresh <- 2.1
grass.grvd.thresh <- -.2
{
  grass <-  sum.im < grass.bright.thresh & grvd > grass.grvd.thresh
  grass.equation <- paste0("sum.im < ", grass.bright.thresh, " & grvd > ", grass.grvd.thresh)
  
  # Compare classification to image
  plot(grass) 
  plotRGB(rgb.stack, scale = 1, add = T)
  plot(grass) # grass will now be green
}



 { 
  # Now do the same for the other two images 
  ninety <- crop.class.image(grass, 90)
  thirty <- crop.class.image(grass, 30)
  
  par(mar=c(1,1,2,1))
  plot(grass) # this is done twice to overlay the rgb image in the second plot
  plotRGB(ninety.stack, scale = 1, add = T)
  plot(ninety) # tree cover will be green
  plot(grass) # this is done twice to overlay the rgb image in the second plot
  plotRGB(thirty.stack, scale = 1, add = TRUE)
  plot(thirty) # tree cover will be green
  
  
  # SOIL 
  # Soil is the rest
  soil <- grass == 0
  soil.equation <- "remainder"
  
 
  
  #############################   Combine Covers  #################################################################
  # Give the different cover types different values
  
  grass.ad <- grass
  grass.ad[grass.ad == 1] <- 2
  grass.ad[is.na(grass.ad)] <- 0
  
  soil.ad <- soil # don't need to set soil to a different value. Just let it be 1
  soil.ad[is.na(soil.ad)] <- 0
  
  tree.ad <- trees
  tree.ad[tree.ad == 1] <- 3
  tree.ad[is.na(tree.ad)] <- 0
  
  # Combine
  tot.img <- tree.ad + soil.ad + grass.ad
  tot.img[tot.img == 0] <- NA
  
  
  
  # One last comparison
  par(mar=c(1,1,2,1)) #plotRGB makes the margins small, so have to set it back to this
  {plot(tot.img)
    plotRGB(rgb.stack, scale = 1, add = T)}
  plot(tot.img, axes=F,box = F, main = "3 = Trees, 2 = Grass, 1 = Soil")
  plot(cell.poly, add = TRUE)
  plot(l.cell.poly, add = TRUE) 
  
  
  
  
  #########  Crop for exporting  ################################################################################################

  ninety <- crop.class.image(tot.img, 90)
  thirty <- crop.class.image(tot.img, 30)
  one80 <- crop.class.image(tot.img, 180)
  
  
  
  # Plot and save final classification (if wanted)
  if(save.class.pixel){ # This will overwrite already existing files (the creation date and time doesn't change, though)
    class.Landscape.name <- paste0(substr(hand.list[image.num], 1, nchar(hand.list[image.num]) - 4), "_classified.png")
    # Save the image
    png(class.Landscape.name)
    plot(tot.img, axes=F,box = F, main = "3 = Trees, 2 = Grass, 1 = Soil")
    plot(cell.poly, add = TRUE)
    plot(l.cell.poly, add = TRUE)
    dev.off()
    
    class.9090.name <- paste0(substr(hand.list[image.num], 1, nchar(hand.list[image.num]) - 4), "_9090_classified.png")
    # Save the image
    png(class.9090.name)
    plot(ninety, axes=F,box = F, main = "3 = Trees, 2 = Grass, 1 = Soil")
    plot(cell.poly, add = TRUE)
    plot(l.cell.poly, add = TRUE)
    dev.off()
    
    class.3030.name <- paste0(substr(hand.list[image.num], 1, nchar(hand.list[image.num]) - 4), "_3030_classified.png")
    # Save the image
    png(class.3030.name)
    plot(thirty, axes=F,box = F, main = "3 = Trees, 2 = Grass, 1 = Soil")
    plot(cell.poly, add = TRUE)
    dev.off()
  }
  
  
  
  
  #############################   Calculate Cover %  ###################################################################
  # Calculate the % covers
  # Get the number of cells in each
  
  tr.cells <- length(thirty[thirty == 3])
  gr.cells <- length(thirty[thirty == 2])
  sl.cells <- length(thirty[thirty == 1])
  tot.cells <- ncell(thirty)
  
  # Compute percentages and what the majority cover is
  p.tree <- tr.cells / tot.cells
  p.grass <- gr.cells / tot.cells
  p.soil <- sl.cells / tot.cells
  
  tr.cells <- length(ninety[ninety == 3])
  gr.cells <- length(ninety[ninety == 2])
  sl.cells <- length(ninety[ninety == 1])
  tot.cells <- ncell(ninety)
  
  # Compute percentages and what the majority cover is
  p.tree.90 <- tr.cells / tot.cells
  p.grass.90 <- gr.cells / tot.cells
  p.soil.90 <- sl.cells / tot.cells
  
  tr.cells <- length(one80[one80 == 3])
  gr.cells <- length(one80[one80 == 2])
  sl.cells <- length(one80[one80 == 1])
  tot.cells <- ncell(one80)
  
  # Compute percentages and what the majority cover is
  p.tree.180 <- tr.cells / tot.cells
  p.grass.180 <- gr.cells / tot.cells
  p.soil.180 <- sl.cells / tot.cells
  
  
  # Define majority cover for Google classification
  if(p.tree > p.grass & p.tree > p.soil){
    maj.cov <- "Trees"
  }else if(p.grass > p.tree & p.grass > p.soil){
    maj.cov <- "Grass"
  }else{maj.cov <- "Soil"}
  
  
  
  #############################  Record Values  ###################################################################
  p.tree
  p.grass
  p.soil
  
  
  
  if(record.class.values){
    hand.truth.df$Count[image.num] <- image.num
    hand.truth.df$hand.tree.30[image.num] <- round(p.tree, 5)
    hand.truth.df$hand.soil.30[image.num] <- round(p.soil, 5)
    hand.truth.df$hand.grass.30[image.num] <- round(p.grass, 5)
    hand.truth.df$truth.maj.cover[image.num] <- maj.cov
    hand.truth.df$hand.tree.90[image.num] <- round(p.tree.90, 5)
    hand.truth.df$hand.soil.90[image.num] <- round(p.soil.90, 5)
    hand.truth.df$hand.grass.90[image.num] <- round(p.grass.90, 5)
    hand.truth.df$hand.tree.180[image.num] <- round(p.tree.180, 5)
    hand.truth.df$hand.soil.180[image.num] <- round(p.soil.180, 5)
    hand.truth.df$hand.grass.180[image.num] <- round(p.grass.180, 5)
    hand.truth.df$tree_equation[image.num] <- tree.equation
    hand.truth.df$grass_equation[image.num] <- grass.equation
    hand.truth.df$soil_equation[image.num] <- soil.equation
    hand.truth.df$tree.bright.thresh[image.num] <- tree.bright.thresh
    hand.truth.df$grass.bright.thresh[image.num] <- grass.bright.thresh
    hand.truth.df$grass.grvd.thresh[image.num] <- grass.grvd.thresh
    hand.truth.df$bright.frac[image.num] <- bright.frac
    hand.truth.df$red.frac[image.num] <- red.frac
    # Add the a flag if the pixel had water in it (1), otherwise set to zero
    if(water.was.present){
      hand.truth.df$water.flag[image.num] <- 1
    }else{hand.truth.df$water.flag[image.num] <- 0}
  }
}


#############################  Save DF  ##########################################################################
if(save.final.class.df){
  write.csv(hand.truth.df, paste0("./wholepark/groundtruth_doneByHand_", YMD(), ".csv"), row.names = F)
}





########################  GET EXTRA IMAGES AS NEEDED  ############################################################

# Get the number of pictures that either weren't available, or had water
water <- sum(na.omit(hand.truth.df$water.flag))
noImage <- sum(hand.truth.df$image.not.available)
# water <- sum(na.omit(hand.truth.df$water.flag)[c(455:nrow(hand.truth.df))])
# noImage <- sum(hand.truth.df$image.not.available[c(455:nrow(hand.truth.df))])
total.missing <- water + noImage

# List the images in the wholepark folder
img.list <- list.files("./wholepark/extraImages", pattern = "MATCHED", full.names = T)

hand.truth.num <- nrow(hand.truth.df)
# Create the site numbers you'll want
site.num <- c(hand.truth.num : (hand.truth.num + total.missing)) + 1
# site.num <- c(nrow(hand.truth.df): (nrow(hand.truth.df) + (total.missing-1))) + 1

# Get an empty row to fill in the hand.df
empty.row <- hand.truth.df[1, ]
empty.row[] <- 0

# Now loop through and get the images that match the site.numbers
for(w in 1:length(site.num)){
 
  # What is the location of the specific image 
  image.index <- which(as.numeric(gsub('^.*_site\\s*|\\s*_pxlnum_.*$', '', img.list)) == site.num[w])
  # Get the specific file name
  image.txt <- img.list[image.index]
  # Contrust the new name and move the file
  image.name.new <- gsub("extraImages", "handImages", image.txt)
  file.copy(image.txt, image.name.new, overwrite = T)
  
  
  ###########  Crop to 30 and 90  #################
  # Create the hand image 30x30 rgb name list
  rgb.3030.name <- paste0(substr(image.name.new, 1, nchar(image.name.new) - 4), "_3030_rgb.png")
  rgb.9090.name <- paste0(substr(image.name.new, 1, nchar(image.name.new) - 4), "_9090_rgb.png")
    

  ## Loop through and crop the hand images. Later code can do the same for the auto.images
  if(save.rgb.pixel){
    
    # Crop the image (function in src_masterfunctions.R)
    img <- crop.google(image.name = image.name.new, pixel.res = pixel.size.meters)
    img2 <- crop.google(image.name = image.name.new, pixel.res = larger.pxl.size)
    
    # Plot and save the RGB image of the cell - not needed if they are already saved
    png(rgb.3030.name) # Save the image
    plotRGB(img, scale = 1)
    dev.off()
    
    png(rgb.9090.name) # Save the image
    plotRGB(img2, scale = 1)
    dev.off()
  }
  
  # Add the names to the image lists used in the classification section above
  hand.list <- c(hand.list, image.name.new)
  rgb.3030.names <- c(rgb.3030.names, rgb.3030.name)
  rgb.9090.names <- c(rgb.9090.names, rgb.9090.name)
  
  # hand.list <- insert(hand.list, 437, NA)
  # rgb.3030.names <- insert(rgb.3030.names, 437, NA)
  # rgb.9090.names <- insert(rgb.9090.names, 437, NA)
  
  # Find and plug in the park and section numbers
  samp.num <- as.numeric(gsub('^.*_sample_\\s*|\\s*_x_.*$', '', image.name.new))
  park.num <- as.numeric(gsub('^.*_pxlnum_\\s*|\\s*_sample_.*$', '', image.name.new))
  
  # Fill the empty row
  empty.row$Count <- site.num[w]
  empty.row$park.pxl.num <- park.num
  empty.row$sample.pxl.num <- samp.num
  
  # Append the row to the whole df
  hand.truth.df <- rbind(hand.truth.df, empty.row)
}

if(save.final.class.df){
  write.csv(hand.truth.df, paste0("./wholepark/groundtruth_doneByHand_", YMD(), ".csv"), row.names = F)
}



# 
# for(k in 1:length(hand.truth.df$image.not.available)){
#   if(is.na(hand.truth.df$image.not.available[k])){
#     hand.truth.df$image.not.available[k] <- 1
#   }
# }
# hand.truth.df$image.not.available



nrow(hand.truth.df) - total.missing


