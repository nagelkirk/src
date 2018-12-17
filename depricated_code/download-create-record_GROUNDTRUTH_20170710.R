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
function.folder <- "F:/Dropbox/Permanent/Grad School/src_functions/" # set this to the path containing the src_masterfunction.R file. Be sure to end the path with "/"
source(paste0(function.folder, "src_masterfunctions.R")) # this reads the functions


#################################  Variables  ##################################################################
# These variables, especially the first half of the set, will need to be changed as you move to new sample areas and/or parks.
park.section <- 1  # this is the section of the boundary file that you want your ground truth points to be created within
park <- "Kruger"  # This is the folder name inside the working directory with all the files for the park
park.boundary <- "Kruger_Boundary" # this is the simple outline of the park
image.boundaries <- "Kruger_Boundary_GE_Buffers_excluded" # this is the boundary file that also has the boundaries of the google imagery
landsat.image <- "Kruger_backdrop_mosaic.tif"  # this needs to be a landsat image that covers the entire sample area. It serves as a template for point creation
working.dir <- "F:/Dropbox/Permanent/Grad School/Projects/EleTree/data/ParkData/" # set this to the folder path that has all the park folders in it. Make sure it has "/" at the end.
truth.num <- 1000  # the total number of points to take from the park.
hand.truth.num <- 200  # the number of points to classify by hand. It can't go higher than truth.num
hand.truth.image.folder <- "manual_ground_truth_images" # this should stay the same. It is the folder within the park folder that will hold the hand-classified images
ground.truth.image.folder <- "automated_ground_truth_images"  # this should also stay the same. Like above, except all sample images are saved here, not just the hand-classified
pixel.size.meters <- 30  # this shouldn't change. It is the length of one side of the area you want to classify. Typically will be set to landsat imagery resolution 


#################################  Control Panel  ##############################################################
# These control various parts of the code and whether they execute. 
# If all are set to false, then nothing gets recorded or saved if the entire script is ran. That is good for when just trying to learn the code.
# Once you have the hang of things, set them to true. In general, though, running the entire script in one go is not recommended. 
# It is better to run only sections of this script at a time, since the classification process is manual and the script 
# is just here to make classification streamlined, not entirely automated. 
get.Google.image <- FALSE # This controls whether the google images are downloaded. That process should only be done onces for each sample area (ref: "Get Google Imagery")
save.rgb.pixel <- FALSE  # Should the cropped 30x30 m RGB of the Google imagery be saved? Again, only needs to be done once. (ref: "Crop to 30x30 Meters")
save.class.pixel <- FALSE  # Should the classified image be saved to the image folder? Only save images once you have the classification you want. (ref: "Combine Covers")
record.class.values <- FALSE # Should the classified veg cover values be plugged into the dataframe? (ref: "Record Values" section)
save.final.class.df <- FALSE # Should the dataframe be saved? This should be done only once all the different images have been classified. (ref: "Save DF")
using.area.dependent.point.generation <- TRUE # this enables the code to automatically determine how many points should be hand classified from each image/sample area

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
google.boundary <- readOGR("./gis", image.boundaries)
park.border <- readOGR("../Boundary", park.boundary)
# Convert Mpala shapefile to lat long, which everything else is in
google.boundary <- spTransform(google.boundary, crs(landsat.image))
park.border <- spTransform(park.border, crs(landsat.image))
plot(google.boundary)


# Get the total area of the images in the park
tot.img.area <- sum(google.boundary$Shape_Area)
# Calculate the point density given the area 
hand.point.density <- hand.truth.num / tot.img.area
tot.point.density <- truth.num / tot.img.area
# Create a new column with number of points for both hand checking and the total drawn from the park
google.boundary$hand_pts <- round(google.boundary$Shape_Area * hand.point.density)
google.boundary$tot_pts <- round(google.boundary$Shape_Area * tot.point.density)
hand.missing.pts <- hand.truth.num - sum(google.boundary$hand_pts) # this is to figure out how many we'll need to add yet
tot.missing.pts <- truth.num - sum(google.boundary$tot_pts) # this is to figure out how many we'll need to add yet
# Find largest image index number and add the remaining points to that image
max.img <- which(google.boundary$Shape_Area == max(google.boundary$Shape_Area))
google.boundary$hand_pts[max.img] <- google.boundary$hand_pts[max.img] + hand.missing.pts 
google.boundary$tot_pts[max.img] <- google.boundary$tot_pts[max.img] + tot.missing.pts
sum(google.boundary$hand_pts)
sum(google.boundary$tot_pts)

# Create background Park image
blank.park <- clipTIF(tifname = landsat.image, clipboundary = park.border)
blank.park[!is.na(blank.park)] <- 0
plot(blank.park, col = "lightgray")



#############################  Point Generation  ###############################################################

# Select the section of the park to extract points from
ground.truth.frame <- as.SpatialPolygons.PolygonsList(google.boundary@polygons[park.section], crs(landsat.image))
plot(ground.truth.frame)


# Clip out the raster and set the values to 1
sample.area.raster <- clipTIF(tifname = blank.park, clipboundary =  ground.truth.frame) 
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
for(x in 1:length(hand.images)){
  rgb.3030.name <- paste0(substr(hand.images[x], 1, nchar(hand.images[x]) - 4), "_3030_rgb.png")
  rgb.3030.names <- c(rgb.3030.names, rgb.3030.name)
}

# Create the google image folders if not already done
if(file.exists(paste0("./section_", park.section, "/", ground.truth.image.folder)) == FALSE){
  dir.create(paste0("./section_", park.section, "/", ground.truth.image.folder))
}

if(file.exists(paste0("./section_", park.section, "/", hand.truth.image.folder)) == FALSE){
  dir.create(paste0("./section_", park.section, "/", hand.truth.image.folder))
}



#############################  Get Google Imagery  #################################################################
# This is the section of code that actually downloads and saves the google imagery
if(get.Google.image){
  
  # Get the Google image for each pixel
  for(t in 1:truth.num){
    # Create the image name for that pixel
    image.name <- auto.images[t] # this looks risky, but the naming above follows this same order
    # GET THE IMAGE FROM GOOGLE!! 
    # Function in src_masterfunctions.R
    j <- getGoogleimage(rastername = sample.area.raster, ras.pixel.num = cellnums.df$sample_cell[t], 
                        destfilename = image.name, typeofmap = "satellite", zoomlevel = 19)
  }
  
  # Move the first n (hand.truth.num) images to the hand classifying folder (they are now in random order (a good thing: it's the original order from random number generator))
  for(t in 1:hand.truth.num){
    image.name <- auto.images[t]
    image.name.new <- hand.images[t]
    file.copy(image.name, image.name.new)
  }
}



#############################  Create 30x30 Shapefile  ##########################################################
# This isn't currently used to crop the images because it was noticeably slower than the crop.google line
cell.size <- pixel.size.meters/190 # 190m is the length of one side of the Google image
half.pix <- cell.size / 2
center <- .5 # Middle of the image. The reference point for indexing

# Create the boundaries of the upper left cell. The 9-pixel loop will modify 
# these as it counts through the cells (extent.list)
bottom <- center - half.pix
right <- center + half.pix
top <- center + half.pix 
left <- center - half.pix # might be able to simplify these. Write whole loop first

#Create extent object 
cr.ext <- extent(c(left, right, bottom, top))

# Coerce extent to a SpatialPolygons object and add to plot
cell.poly <- as(cr.ext, 'SpatialPolygons')



#############################  Crop to 30x30 Meters  ##########################################################
## Loop through and crop the hand images. Later code can do the same for the auto.images
for(x in 1:length(hand.images)){
  # Crop the image (function in src_masterfunctions.R)
  img <- crop.google(image.name = hand.images[x], pixel.res = pixel.size.meters)
  
  # Plot and save the RGB image of the cell - not needed if they are already saved
  if(save.rgb.pixel){
    png(rgb.3030.names[x]) # Save the image
    plotRGB(img, scale = 1)
    dev.off()
  }
}





#############################  START of MANUAL CLASSIFICATION  ###################################################
image.num <- 1 # Manually set this as you go. It represents the point you are classifying. So will be in range of 1-50.
# After the first point is classified, set this to 2, then 3 after that, and so on, until you hit 50. After 50, you're done!

#############################  Index Calculations  ###############################################################
# Read in the single cell and landscape images
current.png <- readPNG(rgb.3030.names[image.num])
landscape.png <- readPNG(hand.images[image.num])

# Pull out the red, green and blue bands from both images
red.c <- raster(current.png[,,1]) # the "c" is just a reminder that it's the color value for the cropped image
green.c <- raster(current.png[,,2])
blue.c <- raster(current.png[,,3])

red <- raster(landscape.png[,,1]) 
green <- raster(landscape.png[,,2])
blue <- raster(landscape.png[,,3])

# Stack so can create images later
rgb.stack <- stack(red.c, green.c, blue.c)
landscape.stack <- stack(red, green, blue)

# Sum the RGB values. Low values are dark areas which correspond to trees. 
sum.im.pre <- green.c + red.c + blue.c

# Calculate the Green-Red Vegetation Index
grvd.pre <- (green.c - red.c) / (green.c + red.c)

# Calculate summary stats for images to determine if it is mostly grass or soil
#  This is done because images with a lot of grass or soil require different thresholds for classifying tree cover

# Find the percent area of the cell that has a brightness sum over 2
bright <- red.c + green.c + blue.c
bright.area <- bright > 2
bright.frac <- length(bright.area[bright.area == 1]) / length(bright.area)


# Look at pixels with higher red values than green
red.area <- red.c > green.c
red.frac <- length(red.area[red.area == 1]) / length(red.area) 



# Below, classifier is able to set thresholds as they like. At end, all thresholds along with every index and measure 
# are recorded in the df.That way, can analyze the values (simple regression) to see if there are relationships between
# the thresholds chosen and the different brightness values etc.



#############################  Classify Landcover  ###############################################################
# This is the section where you will have to iteratively change the thresholds to create the classification you like.
# First get tree cover to where you want it, then move on to grass cover. Soil cover will be done automatically as the
# remaining area. 

# First, plot the image of the overall area with the cell overlaid
plotRGB(landscape.stack, scale = 1)
plot(cell.poly, add = TRUE)

# Do trees, then grass. Soil will automatically be calculated after grass (soil = remaining pixels)
# TREES 
tree.bright.thresh <- 1.4  # Change this setting until you like the tree cover 
trees <- sum.im.pre < tree.bright.thresh
tree.equation <- paste0("sum.im.pre < ", tree.bright.thresh)
# Compare classification to image
par(mar=c(1,1,2,1))
plot(trees) # this is done twice to overlay the rgb image in the second plot
plotRGB(rgb.stack, scale = 1, add = TRUE)
plot(trees) # tree cover will be green

# Now remove the tree pixels from the image that will go into identifying soil and grass
# create tree mask
tree.mask <- trees
tree.mask[tree.mask == 1] <- 50
holy.image <- grvd.pre + tree.mask
holy.image.3 <- sum.im.pre + tree.mask
holy.image[holy.image > 5] <- NA
holy.image.3[holy.image.3 > 5] <- NA
grvd <- holy.image
sum.im <- holy.image.3


# GRASS
# Split soil and grass using GRVD and brightness. Soil is nonveg and bright
grass.bright.thresh <- 2.1
grass.grvd.thresh <- -0.05
grass <-  sum.im < grass.bright.thresh & grvd > grass.grvd.thresh
grass.equation <- paste0("sum.im < ", grass.bright.thresh, " & grvd > ", grass.grvd.thresh)

# Compare classification to image
plot(grass) 
plotRGB(rgb.stack, scale = 1, add = T)
plot(grass) # grass will now be green

# SOIL 
# Soil is the rest
soil <- grass == 0
soil.equation <- "remainder"

# IMPORTANT: Was there water in the classified area? If yes, set this variable to "TRUE"
water.was.present <- "FALSE"

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
plotRGB(landscape.stack, scale = 1)
plot(cell.poly, add = TRUE)
par(mar=c(1,1,2,1)) #plotRGB makes the margins small, so have to set it back to this
plot(tot.img)
plotRGB(rgb.stack, scale = 1, add = T)
plot(tot.img, axes=F,box = F, main = "3 = Trees, 2 = Grass, 1 = Soil")

# Plot and save final classification (if wanted)
if(save.class.pixel){ # This will overwrite already existing files (the creation date and time doesn't change, though)
  class.3030.name <- paste0(substr(hand.images[x], 1, nchar(hand.images[x]) - 4), "_3030_classified.png")
  # Save the image
  png(class.3030.name)
  plot(tot.img, axes=F,box = F, main = "3 = Trees, 2 = Grass, 1 = Soil")
  dev.off()
}


#############################   Calculate Cover %  ###################################################################
# Calculate the % covers
# Get the number of cells in each
tr.cells <- length(tot.img[tot.img == 3])
gr.cells <- length(tot.img[tot.img == 2])
sl.cells <- length(tot.img[tot.img == 1])
tot.cells <- ncell(tot.img)

# Compute percentages and what the majority cover is
p.tree <- tr.cells / tot.cells
p.grass <- gr.cells / tot.cells
p.soil <- sl.cells / tot.cells


# Define majority cover for Google classification
if(p.tree > p.grass & p.tree > p.soil){
  maj.cov <- "Trees"
} else if(p.grass > p.tree & p.grass > p.soil){
  maj.cov <- "Grass"
} else{maj.cov <- "Soil"}



#############################  Record Values  ###################################################################
if(record.class.values){
  hand.truth.df$Count[image.num] <- image.num
  hand.truth.df$park.pxl.num[image.num] <- cellnums.df$park_cell[image.num]
  hand.truth.df$sample.pxl.num[image.num] <- cellnums.df$sample_cell[image.num]
  hand.truth.df$hand.tree[image.num] <- p.tree
  hand.truth.df$hand.soil[image.num] <- p.soil
  hand.truth.df$hand.grass[image.num] <- p.grass
  hand.truth.df$truth.maj.cover[image.num] <- maj.cov
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


#############################  Save DF  ##########################################################################
if(save.final.class.df){
  write.csv(hand.truth.df, paste0("./section_", park.section, "/groundtruth_doneByHand_", YMD(), ".csv"))
}











#### THIS was specific to Mpala. Was moving over old imagery that had been collected from Google so that wouldn't have to do again
# # Move the images over to new folder from the old 2017 folder
# old.folder.path <- "C:/Users/nagelki-4/Desktop/nagelki4/Grad School/Projects/EleTree Analysis/SMA/Rasters/test_images_20170102"
# image.list <- list.files(old.folder.path, full.names = T)
# 
# # Only take the big ones
# final.list <- c() # this will be filled
# for(i in 1:length(image.list)){
#   size <- file.size(image.list[i])
#   if(size > 350000){
#     final.list <- c(final.list, image.list[i])
#   }
# }
# 
# for(t in 1:500){
#   image.name <- final.list[t]
#   image.name.new <- paste0("./", ground.truth.image.folder, "/pixel_park", park.cell.nums[t, 1], 
#                            "_sample", park.cell.nums[t, 2], "_", park, ".png")
#   file.copy(image.name, image.name.new)
# }
