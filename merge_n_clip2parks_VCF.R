# VCF merger

# Ryan Nagelkirk
# 20180913

# Description
# Reads in and merges all the Landsat VCF (technically called TC) for each park, merges them and clips the park out

# Purpose:
# Create VCF layers for each park

# Load libraries
library(raster)
library(rgdal)

# Set working directory
setwd("X:/nagelki4/Projects/EleTree/data/ParkData")

# Change the memory limit to handle the larger files 
memory.limit(size = 50000) # default: 16308

# Source the functions
source("E:/Dropbox/Permanent/Grad School/src_functions/src_masterfunctions.R")

#################################  Variables  #################################################################
park.names <- c("Chobe", "Kruger", "Limpopo", "Mpala", "Murchison", "North_Luangwa", "QWE", "Ruaha", "Selous", "Serengeti", "South_Luangwa")

# Load africa
africa <- readOGR("E:/Dropbox/Permanent/Grad School/Projects/EleTree/data/GIS/AfricanCountries", "AfricanCountires")

for(i in 1:length(park.names)){
  # Start for loop here
  park.name <- park.names[i]
  # If the clipped raster already exists, go to the next park
  if(file.exists(paste0(getwd(), "/", park.name, "/VCF/", park.name, "_merged_n_clipped_VCF.tif"))){
    print(paste0(park.name, "_merged_n_clipped_VCF.tif", " already exists. Went to next park."))
    next
  }
  
  
  ###############################  DIRECTORIES  ##############################################################
  # Folder with shape files of ground truth area
  shape.folder <- paste0("E:/Dropbox/Permanent/Grad School/Projects/EleTree/data/ParkData/", park.name, "/Boundary")
  # shape.folder <- paste0("~/Dropbox/Permanent/Grad School/Projects/EleTree/data/ParkData/", park.name, "/Boundary")
  
  ############  Park Boundary  ##############################################################################
  # Get the boundary
  if(park.name == "Serengeti_Mara"){
    park.boundary <- "Serengeti_Mara_dissolve"
  }else if(park.name == "Kruger"){
    park.boundary <- "Kruger_Boundary"
  }else if(park.name == "Mpala"){
    park.boundary <- "Mpala_Boundary"
  }else{
    park.boundary <- park.name
  }
  
  ########################  LOAD VCF FILES  #########################################################
  # Read in the file
  park.border <- readOGR(shape.folder, park.boundary)
  plot(park.border)

  # Get a list of the files
  vcf.files <- list.files(paste0("./", park.name, "/VCF"), full.names = T, pattern = "TC_2015")
  # Get rid of the .gz files
  vcf.files <- vcf.files[-grep(pattern = ".gz", vcf.files)]
  # # Get rid of any merged file that might exist
  # if(length(list.files(paste0("./", park.name, "/VCF"), full.names = T, pattern = "merged")) > 0){
  #   vcf.files <- vcf.files[-grep(pattern = "merged", vcf.files)]
  # }
  
  # Load each raster
  for(n in 1:length(vcf.files)){
    new.file <- raster(vcf.files[n])
    if(n == 1){
      raster.list <- list(new.file)
    }else{
      raster.list <- c(raster.list, new.file)
    }
  }
  
  ###########  MERGE  ###############################################################################
  if(length(raster.list) > 1){
    # Prepare for merge call
    raster.list$filename <- paste0(getwd(), "/", park.name, "/VCF/", park.name, "_merged_VCF.tif")
    raster.list$overwrite <- TRUE
    
    # Merge the rasters 
    m <- do.call(merge, raster.list)
  }else{
    m <- new.file
  }
  
  ############  CLIP  ###############################################################################
  new.border <- spTransform(park.border, crs(new.file))
  plot(m)
  plot(new.border, add = T)
  
  # Clip the raster
  m.clip <- clipTIF(tifname = m, clipboundary = new.border)
  
  ###############  Save the file  ###################################################################
  writeRaster(m.clip, paste0(getwd(), "/", park.name, "/VCF/", park.name, "_merged_n_clipped_VCF.tif"))
  print(paste(park.name, "complete"))
}



