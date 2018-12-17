#################################  CHECKING GROUND TRUTH BY HAND  ##############################################
# Author: Ryan Nagelkirk
# Date: 5/15/2017

# Description:

# This was edited on 20170904 to get images from the entire park as a back up dataset in case changes are made to the Google Earth API

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
# 2. Within the ground_truth folder, create a "gis" folder. Save the park boundary shapefile here. The boundary
#    will need to be the file with the boundaries of the google imagery.
# 3. Save a landsat image, or images, that cover the different areas that are going to be classified. For a small
#    park, this will generally be a single image that covers the entire park. For larger parks, several images might
#    be needed. It doesn't matter what band is used because this image is only used for indexing cells. 
#    What does matter is whether the image has the same shape and extent as all the other corresponding path/row images. 
#    This can be checked in ArcGIS. Ex: if the scene you're using is from path/row 168/078, then the image you chose 
#    should perfectly overlay other images of 168/078. Sometimes this is not the case and those images need to be discarded.
# 4. Once 1-3 are done, you're ready to go! The code itself is heavily annotated to allow you to work through it from here.
#    If there is not annotation, just run those lines - it isn't necessary for you to do anything else with them. 


#################################  Load Libraries  #############################################################
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


#################################  Functions  ##################################################################
# Source the functions
function.folder <- "E:/Dropbox/Permanent/Grad School/src_functions/" # set this to the path containing the src_masterfunction.R file. Be sure to end the path with "/"
source(paste0(function.folder, "src_masterfunctions.R")) # this reads the functions


#################################  Variables  ##################################################################
# These variables, especially the first half of the set, will need to be changed as you move to new sample areas and/or parks.
park.section <- 1  # this is the section of the boundary file that you want your ground truth points to be created within
park <- "Tuli"  # This is the folder name inside the working directory with all the files for the park
landsat.image <- "tuli_backdrop_mosaic.tif"  # this needs to be a landsat image that covers the entire sample area. It serves as a template for point creation
working.dir <- "E:/Dropbox/Permanent/Grad School/Projects/EleTree/data/ParkData/" # set this to the folder path that has all the park folders in it. Make sure it has "/" at the end.
truth.num <- 500  # the total number of points to take from the park.
hand.truth.num <- 0  # the number of points to classify by hand. It can't go higher than truth.num
hand.truth.image.folder <- "manual_ground_truth_images" # this should stay the same. It is the folder within the park folder that will hold the hand-classified images
ground.truth.image.folder <- "automated_ground_truth_images"  # this should also stay the same. Like above, except all sample images are saved here, not just the hand-classified
pixel.size.meters <- 90  # this shouldn't change. It is the length of one side of the area you want to classify. Typically will be set to landsat imagery resolution 
park.boundary <- "Tuli_Entire_Boundary" # this is the simple outline of the park
image.boundaries <- park.boundary

#################################  Control Panel  ##############################################################
# These control various parts of the code and whether they execute. 
# If all are set to false, then nothing gets recorded or saved if the entire script is ran. That is good for when just trying to learn the code.
# Once you have the hang of things, set them to true. In general, though, running the entire script in one go is not recommended. 
# It is better to run only sections of this script at a time, since the classification process is manual and the script 
# is just here to make classification streamlined, not entirely automated. 
get.Google.image <- TRUE # This controls whether the google images are downloaded. That process should only be done onces for each sample area (ref: "Get Google Imagery")
save.rgb.pixel <- TRUE  # Should the cropped 30x30 m RGB of the Google imagery be saved? Again, only needs to be done once. (ref: "Crop to 30x30 Meters")
save.class.pixel <- FALSE  # Should the classified image be saved to the image folder? Only save images once you have the classification you want. (ref: "Combine Covers")
record.class.values <- FALSE # Should the classified veg cover values be plugged into the dataframe? (ref: "Record Values" section)
save.final.class.df <- FALSE # Should the dataframe be saved? This should be done only once all the different images have been classified. (ref: "Save DF")
using.area.dependent.point.generation <- FALSE # this enables the code to automatically determine how many points should be hand classified from each image/sample area

# Classifier note: From here all the way down until the "START OF MANUAL CLASSIFICATION" section can be ran without any input.
# Just make sure "get.Google.image" and "save.rgb.pixel" are set to false if the imagery is already downloaded. Otherwise it 
# will download again and you'll have to wait for that. Also, the imagery might have changed and you'll have to do the entire 
# classification all over if the old images are overwritten. 

#############################  Set Working Directory  ##########################################################
# Set working directory to where the imagery will be saved, along with the final table
setwd(paste0(working.dir, park, "/ground_truth"))
# par(mfrow = c(1,1)) # just in case needed for changing later
par(mar=c(1,1,2,1))

############################  Create DFs  ######################################################################
# DF for manual classification entries
mx <- matrix(0, nrow = hand.truth.num, ncol = 7)
c.names <- c("Count", "park.pxl.num", "sample.pxl.num", "hand.tree", "hand.grass", "hand.soil", "truth.maj.cover")
colnames(mx) <- c.names # above, pxl.num is the pixel number within the sample area, park.pxl.num is from the entire park
hand.truth.df <- as.data.frame(mx)

# DF for automated classification
mx <- matrix(0, nrow = truth.num, ncol = 7)
c.names <- c("Count", "park.pxl.num", "pxl.num", "hand.tree", "hand.grass", "hand.soil", "maj.truth")
colnames(mx) <- c.names # above, pxl.num is the pixel number within the sample area, park.pxl.num is from the entire park
auto.truth.df <- as.data.frame(mx)



#############################  Create Points  #################################################################
###### Read in Landsat and shapefiles
# Landsat
landsat.image <- raster(landsat.image)

# Shapefiles
park.border <- readOGR("../Boundary", park.boundary)
google.boundary <- park.border
# Convert Mpala shapefile to lat long, which everything else is in
google.boundary <- spTransform(google.boundary, crs(landsat.image))
park.border <- spTransform(park.border, crs(landsat.image))
plot(google.boundary)


# Get the total area of the images in the park
tot.img.area <- sum(google.boundary$Shape_Area)


# Create background Park image
blank.park <- clipTIF(tifname = landsat.image, clipboundary = park.border)
blank.park[!is.na(blank.park)] <- 0
plot(blank.park, col = "lightgray")



#############################  Point Generation  ###############################################################

# Select the section of the park to extract points from
ground.truth.frame <- as.SpatialPolygons.PolygonsList(google.boundary@polygons[park.section], crs(landsat.image))
plot(ground.truth.frame)


# Clip out the raster and set the values to 1
sample.area.raster <- blank.park 
sample.area.raster[!is.na(sample.area.raster)] <- 1
plot(sample.area.raster, axes=F,box = F)


####  Generate List of Points 
# Create list of cells in Mpala for sampling
# First, how many cells have a value?
cell.vals <- length(na.omit(sample.area.raster[])) # this used to be na.omit(mpala[])

# Use that to assign new values to the cells (numbers them)
sample.area.raster[!is.na(sample.area.raster)] <- 1:cell.vals
plot(sample.area.raster, axes=F,box = F)

#### Redefine what truth.num and hand.truth.num are at this point if these numbers are being created based on the image area
if(using.area.dependent.point.generation == TRUE){
  truth.num <- google.boundary$tot_pts[park.section]
  hand.truth.num <- google.boundary$hand_pts[park.section]
}

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
  sing.cell <- Which(sample.area.raster == i, cells= TRUE) # the sample.area.raster values aren't actual cells, but just cells values
  sample.cell.numbers <- c(sample.cell.numbers, sing.cell)
}

# Switch the sample area to NA and then give the sample points their pixel number value
trial <- sample.area.raster
trial[] <- NA
trial[sample.cell.numbers] <- sample.cell.numbers # assign values for the cells based on their pixel number
plot(trial)

# Convert cells to points and add them to the plot
hand.points <- rasterToPoints(trial) # make the ground truth hand classified points
par(mar=c(1,1,2,1)) 
plot(blank.park, axes=F,box = F, legend = FALSE, col = "lightgray")
plot(ground.truth.frame, add = TRUE)
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
cellnums.df <- as.data.frame(park.cell.nums)
names(cellnums.df) <- c("park_cell", "sample_cell")

# Order to match sample.cell.numbers ordering
cellnums.df <- cellnums.df[match(sample.cell.numbers, cellnums.df$sample_cell),]
rownames(cellnums.df) <- 1:nrow(cellnums.df) # order the row names so they aren't confusing


#############################  Create Image Name Lists  ###########################################################
# Create lists of the files to crop (auto and hand)
auto.images <- c()
for(t in 1:truth.num){
  # Create the image name for that pixel
  image.name <- paste0("./section_", park.section, "/", ground.truth.image.folder, "/pixel_park", cellnums.df$park_cell[t], 
                       "_sample", cellnums.df$sample_cell[t], "_", park, ".png")
  auto.images <- c(auto.images, image.name)
}

# Create hand image list
hand.images <- c()
for(t in 1:hand.truth.num){
  image.name.new <- paste0("./section_", park.section, "/", hand.truth.image.folder, "/pixel_park", cellnums.df$park_cell[t], 
                           "_sample", cellnums.df$sample_cell[t], "_", park, ".png")
  hand.images <- c(hand.images, image.name.new)
}

# Create the hand image 30x30 rgb name list
rgb.3030.names <- c()
for(x in 1:length(auto.images)){
  rgb.3030.name <- paste0(substr(hand.images[x], 1, nchar(hand.images[x]) - 4), "_9090_rgb.png")
  rgb.3030.names <- c(rgb.3030.names, rgb.3030.name)
}

# Create the google image folders if not already done
if(file.exists(paste0("./section_", park.section, "/", ground.truth.image.folder)) == FALSE){
  dir.create(paste0("./section_", park.section, "/", ground.truth.image.folder), recursive = T)
}

if(file.exists(paste0("./section_", park.section, "/", hand.truth.image.folder)) == FALSE){
  dir.create(paste0("./section_", park.section, "/", hand.truth.image.folder), recursive = T)
}


##
# Extract that pixel
pxl <- rasterFromCells(sample.area.raster, cellnums.df$sample_cell[t])


cellnums.df$sample_cell[t]

#############################  Get Google Imagery  #################################################################
# This is the section of code that actually downloads and saves the google imagery
if(get.Google.image){
  
  # Get the Google image for each pixel
  for(t in 1:truth.num){
    
    # Create the image name for that pixel
    image.name <- auto.images[t] # this looks risky, but the naming above follows this same order
    
    # Check if the file already exists first
    if(!file.exists(image.name)){
      # GET THE IMAGE FROM GOOGLE!! 
      # Function in src_masterfunctions.R
      j <- getGoogleimage(rastername = sample.area.raster, ras.pixel.num = cellnums.df$sample_cell[t], 
                          destfilename = image.name, typeofmap = "satellite", zoomlevel = 19)
    }
    gc() #frees up memory and allows function to run with fewer interruptions
  }
}







