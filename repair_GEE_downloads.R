# GEE Image Output Merge Code

# PURPOSE
# Take the parsed images and merge them into their complete path as an hdr file. 
# Images that don't require merging are converted to an hdr file as well. 
# All files are placed back in the original folder. 
# The original, unmerged tif files are kept the first time the code runs. Run it once more to delete them - this allows
# the user to check whether the merge was effective before deleting the files. Tifs that don't require merging are kept. 


# Library
library(raster)

# Set working directory
setwd("X:/nagelki4/Projects/EleTree/data")
#setwd("E:/Dropbox/Permanent/Grad School/Projects/EleTree/data")
#setwd("~/Dropbox/Permanent/Grad School/Projects/EleTree")

# Source the functions
source("E:/Dropbox/Permanent/Grad School/src_functions/src_masterfunctions.R")
#source("~/Dropbox/Permanent/Grad School/src_functions/src_masterfunctions.R")

# Change the memory limit to handle the larger files 
memory.limit(size = 50000) # default: 16308

#########  CREATE VARIABLES  #############################################################################################
# Read in a table with the park names 
image_download_df <- read.csv("E:/Dropbox/Permanent/Grad School/Projects/EleTree/data/LandsatDownload_inputs/Image_Date_tables/Landsat_image_4year_span_dates_singleParks_twentykmbuffer_shp_n_parkCloudCover.csv")
# Create lists of names
park_names <- unique(image_download_df$GEE_park)

# Set a destination for the final hdr files (what ENVI reads)
image_folder <- "./Landsat/20km_Buffer/With_shoulder"

# Start a timer for the entire code
process_start <- tick()

# Start a timer to time how long it takes to merge the files that need it
all_merge_start <- tick()

#k <- 2
# Loop through parks
for(k in 1:length(park_names)){
  # Find the files names that should be stitched based on park name and "0000000" identifier
  fls <- list.files(image_folder, pattern = "0000000000", full.names = T)
  park_fls <- fls[grepl(park_names[k], fls)]
  park_fls <- park_fls[grepl("reflectanceBands", park_fls)]
  
  # If there are no files like that for the park, go to next park
  if(length(park_fls) == 0){
    print(paste0("No images to merge for ", park_names[k], ", went to next park."))
    next
  }
  # Create a list of the dates from that list
  # Remove the 0000000000000
  fls_short <- gsub("-.*", "", park_fls)
  
  # Get unique paths and names (will need unique years, too, once going back in time)
  b.names <- basename(fls_short)
  paths.long <- sapply(strsplit(b.names, "LC08"), "[", 2) # This removed the text before the path number, so trim off the rest now
  paths <- substr(paths.long, 1, 4)  # keep the underscore - it helps later on
  landsat_file_names <- unique(fls_short)
  unique.paths <- unique(paths)
  
  
  
  #########  LOOP THROUGH DATES #############################################################################################
  # For every date, find the files that correspond and stitch them if needed 
  
  #p <- 1
  for(p in 1:length(unique.paths)){
    path.images <- park_fls[grepl(unique.paths[p], park_fls)]
    path_prefixes <- unique(gsub("-.*", "", path.images))
    
    #i <- 1
    # Loop through the files from the same path, loading the rasters, stitching them, and exporting
    for(i in 1:length(path_prefixes)){
      # Start timer
      one.date.time <- tick()
      # Get the images for that date
      images_to_merge <- park_fls[grepl(path_prefixes[i], park_fls)]
      file_prefix <- unique(gsub("-.*", "", images_to_merge))
      
      # Get the satellite number
      sat.long <- sapply(strsplit(file_prefix, "LC0"), "[", 2) # This removed the text before the path number, so trim off the rest now
      sat_num <- substr(sat.long, 1, 1)
      
      
      # Define break between normal bands and indices based on the satellite number (LC8 has an extra band)
      # The indices are the last 4 bands (2 for each stacked image)
      if(sat_num == 8){
        band_breaker <- 9
      } else {
        band_breaker <- 8
      }
      
      ###########  CHECK IF FILE ALREADY EXISTS, NEXT IF IT DOES  ###################################################################
      # Create the final name
      final.tif.name <- file_prefix
    
      # If the merged hdr and tif files exist, delete the unmerged files and go to next pair in loop
      if(file.exists(paste0(substr(final.tif.name, 1, nchar(final.tif.name)-25), "maxmintran_NDVI_wo_indices_allBands.tif"))){
        print(paste0("Files already merged and converted to tif. Parent files being deleted: ", images_to_merge))
        #file.remove(images_to_merge)
        next
      }
    
      ###########  LOOP THROUGH IMAGES FROM SAME DATE AND MERGE  ###################################################################
      # Loop through the two images
      #j <- 1
      for(j in 1:length(images_to_merge)){
        
        # Merge if needed, if not, move to the next
        if(length(images_to_merge) == 1){
          print(paste0("Single image (no merging): ", images_to_merge[j]))
          next
        }else{
          print(paste0(i, " of ", length(path_prefixes), ". File being created: ", final.tif.name))
        }
        
        # Load the image as a brick
        brick1 <- brick(images_to_merge[j])
        #plot(brick1[[1]])
        
        ########  WORKING TO FIX STITCHING  ##############
        # Do a subset
        fixstitch <- F
        if(fixstitch){
          if(j == 1){
            # Want the bottom rows
            # Get the dims
            brick1 <- crop(brick1, extent(brick1, dim(brick1)[1]-99, dim(brick1)[1], 1, dim(brick1)[2]))
            plot(brick1)
          }else{
            brick1 <- crop(brick1, extent(brick1, 1, 100, 1, dim(brick1)[2]))
            plot(brick1)
          }
        }
        
        # Build the list of images (each brick is a landsat image with all its bands)
        if(j == 1){
          brick_list <- c(brick1)
        }else{
          brick_list <- c(brick_list, brick1)
        }
      }
      
      ##########  MERGE  ####################################################################################################
      # Just in case it doesn't break out of the function above with images_to_merge == 1, this is written as an if statement
      # This is a hold over from another code. But it doesn't hurt. Images_to_merge should always be greater than 1
      # Also loads any prior merged file
      if(file.exists(paste0(substr(final.tif.name, 1, nchar(final.tif.name)-17), "_wo_indices_allBands.tif"))){
        brick_merged <- brick(paste0(substr(final.tif.name, 1, nchar(final.tif.name)-17), "_wo_indices_allBands.tif"))
      }else if(length(images_to_merge) > 1){
        merge.start <- tick()
        print(paste0("Merging file: ", final.tif.name))
        brick_merged <- do.call(merge, brick_list) 
        print(paste("Merge time:", tock(merge.start)))
      }
      
      
      #########  CHANGE ZEROS, HIGH and LOW VALUES to NA  ################################################################################
      # For every band, set invalid values to NA
      #t <- 1
      for(t in 1:dim(brick_merged)[3]){
        # Set 0 to NA
        NAvalue(brick_merged[[t]]) <- 0
        
        # Trim the top and bottom values based on the band (NDVI and MSAVI2 have diff value ranges)
        if(t <= band_breaker){ 
          # values below zero or >10000 are no good
          values(brick_merged[[t]])[values(brick_merged[[t]]) > 10000] <- NA # Set high and low values to NA
          values(brick_merged[[t]])[values(brick_merged[[t]]) <= 0] <- NA # Set anything below zero to NA. This wasn't done for Kruger and Mpala, but those points will be outliers in the TGS maps, and eliminated there. 
        }else{
          # Change the values outside the range to NA (these are the index bands, though they should no longer exist)
          values(brick_merged[[t]])[values(brick_merged[[t]]) > 1000] <- 1000
          values(brick_merged[[t]])[values(brick_merged[[t]]) < -1000] <- -1000 # Could also increase index values to 10,000 scale by *10 here
        }
      }
      
      
      
      ##########  EXPORT  ####################################################################################################
      # Save the raster stack as a tif and hdr
      export_start <- tick()
      
      # Define band numbers
      # Because these can be different for the min and max images, the images are controlled and stacked in the if statements below this section
      if(sat_num == 8){
        red <- 4
        nir <- 5
        swir <- 6
      }else{
        red <- 3
        ndvi <- 4
        swir <- 5
      }
      
      # Write the merged rasters (but only if they haven't been already)
      if(!file.exists(paste0(substr(final.tif.name, 1, nchar(final.tif.name)-17), "_wo_indices_allBands.tif"))){
        writeRaster(brick_merged, paste0(substr(final.tif.name, 1, nchar(final.tif.name)-17), "_wo_indices_allBands"), format = "GTiff", overwrite = TRUE, datatype = 'INT2S')
      }
      if(!file.exists(paste0(substr(final.tif.name, 1, nchar(final.tif.name)-17), "_wo_indices_allBands.hdr"))){
        writeRaster(brick_merged, paste0(substr(final.tif.name, 1, nchar(final.tif.name)-17), "_wo_indices_allBands"), format = "ENVI", overwrite = TRUE, datatype = 'INT2S')
        writeRaster(brick_merged[[c(red, nir, swir)]], paste0(substr(final.tif.name, 1, nchar(final.tif.name) - 17), "_red_nir_swir1"), format = "ENVI", overwrite = TRUE, datatype = 'INT2S') 
      }
      print(paste("Write time:", tock(export_start)))
      print(paste("Total time for image:", tock(one.date.time))) 
      #########  Construct and write the maxmintrans images  
      # Stack the min and max images here
      if(i == 1){
        maxmintran_brick <- brick_merged
        rednirswir_brick <- brick_merged[[c(red, nir, swir)]]
      }else if(i > 1){
        # Stack them
        maxmintran_brick <- stack(c(maxmintran_brick, brick_merged))
        rednirswir_brick <- stack(c(rednirswir_brick, brick_merged[[c(red, nir, swir)]]))
      }
    }
    
    # Export
    # Write the maxmintran images 
    writeRaster(rednirswir_brick, paste0(substr(final.tif.name, 1, nchar(final.tif.name) - 25), "maxmintran_NDVI_red_nir_swir1"), format = "ENVI", overwrite = TRUE, datatype = 'INT2S') 
    writeRaster(maxmintran_brick, paste0(substr(final.tif.name, 1, nchar(final.tif.name)-25), "maxmintran_NDVI_wo_indices_allBands"), format = "ENVI", overwrite = TRUE, datatype = 'INT2S')
    
    # Then remove the temp files that have been stored on the C drive. Hopefully everything but the rasters was stored on RAM
    removeTmpFiles()
  }
}
# Report the total time it took to merge all the images  
print(paste("Total time to merge all images:", tock(all_merge_start)))






#########  CONVERT COMPLETE TIFS TO HDR (THE ONES SKIPPED ABOVE) #############################################################################################

tif_2_hdr_files_start <- tick()

# Make a tif of the minmax, then hdfs of the min, max and minmax for allbands and rednirswir
# k <- 4
for(k in 1:length(park_names)){
  tif_2_hdr_files_start <- tick()
  
  fls2 <- list.files(image_folder, pattern = ".tif$", full.names = T) # the dollar sign designates .tif as the end of the string
  fls2 <- fls2[grepl(park_names[k], fls2)]
  fls2 <- fls2[grepl(pattern = "reflectanceBands", fls2)]
  # fls2 <- fls2[grepl(pattern = "LC08", fls2)]
  
  # Get only the files that are complete (don't have 000000)
  complete_fls <- fls2[!grepl("00000000", fls2)]
  
  # Remove the tif
  complete_fls <- gsub(".tif.*", "", complete_fls)
  
  # Get unique paths and names (will need unique years, too, once going back in time)
  b.names <- basename(complete_fls)
  paths <- substr(b.names, nchar(b.names)-47, nchar(b.names)-44) # keep the underscore - it helps later on
  landsat_file_names <- unique(complete_fls)
  unique.paths <- unique(paths)
  
  # If there were no files, skip to the next park (the if statement below won't catch this, it just tests whether the envi image was created already)
  if(length(unique.paths) < 1){
    print(paste0("No reflectanceBands images for ", park_names[k], ". Going to next park."))
    next
  }
  
  
  # p <- 1
  for(p in 1:length(unique.paths)){
    # Get the images that are in the specific path
    path.images <- landsat_file_names[grepl(unique.paths[p], landsat_file_names)]
    
    
    # i <- 1
    # Loop through the files from the same path, loading the rasters and exporting.
    # At the end, stack them into the minmax and export that
    for(i in 1:length(path.images)){
      # Start timer
      one.date.time <- tick()

      # Get the images for that date
      images_to_merge <- park_fls[grepl(path_prefixes[i], park_fls)]
      file_prefix <- unique(gsub("-.*", "", images_to_merge))
      
      # Get the satellite number
      sat.long <- sapply(strsplit(file_prefix, "LC0"), "[", 2) # This removed the text before the path number, so trim off the rest now
      sat_num <- substr(sat.long, 1, 1)
      
      # Define break between normal bands and indices based on the satellite number (LC8 has an extra band)
      # The indices are the last 4 bands (2 for each stacked image)
      if(sat_num == 8){
        band_breaker <- 9
      } else {
        band_breaker <- 8
      }
      
      
      ###########  CHECK IF FILE ALREADY EXISTS, NEXT IF IT DOES  ###################################################################
      
      # If the file exists, go to next in loop. This is just checking one of the many hdrs that might be created by this loop. 
      if(file.exists(paste0(substr(file_prefix, 1, nchar(file_prefix)-25), "maxmintran_NDVI_wo_indices_allBands", ".tif"))){
        print(paste0("Already created ", substr(file_prefix, 1, nchar(file_prefix)-25), "maxmintran_NDVI_wo_indices_allBands", ".tif"))
        next
      }
      
      # If it passes the test above, announce that it is being created 
      print(paste0( i, " of ", length(path.images), ". File being created: ", substr(file_prefix, 1, nchar(file_prefix)-22), "_NDVI_wo_indices_allBands", ".tif"))
      
      # Load the image as a brick
      brick2 <- brick(images_to_merge)
      #plot(brick2[[1]])
      
      # Get rid of the cloud and qa bands to reduce the number of layers that need processing
      num_of_bands <- dim(brick2)[3]
      
      
      #########  CHANGE ZEROS, HIGH and LOW VALUES to NA  ################################################################################
      # For every band, set invalid values to NA
      for(t in 1:dim(brick2)[3]){
        # Set zero to NA
        NAvalue(brick2[[t]]) <- 0
        
        # have to replace <= with %in% in the if statement below if doing a multi-image image
        # Trim the top and bottom values based on the band (NDVI and MSAVI2 have diff value ranges)
        if(t <= band_breaker){ # This is specific to LC8, the others have one less band. Subtracting 3 is to account for removing the three bands earlier
          # values below zero or >10000 are no good
          values(brick2[[t]])[values(brick2[[t]]) > 10000] <- NA
          values(brick2[[t]])[values(brick2[[t]]) < 0] <- NA # do this to avoid the zero becoming an NA
        }else{
          # The last two bands are NDVI and MSAVI2, which need to be adjusted to a diff range of values than the landsat bands 
          # The indices were multiplied by 1000 in GEE, so the range for NDVI would be 1000 to -1000. The same goes for MSAVI2 
          # Change the values outside the range to NA
          values(brick2[[t]])[values(brick2[[t]]) > 1000] <- 1000 
          values(brick2[[t]])[values(brick2[[t]]) < -1000] <- -1000
        }
      }
      
      ##########  EXPORT  ####################################################################################################
      # Save the raster stack as an ENVI hdr file
      # Select specific bands wanted for classification later in MESMA
      
      export_start <- tick()
      
      # Define band numbers
      # Because these can be different for the min and max images, the images are controlled and stack in the if statements below this section
      if(sat_num == 8){
        red <- c(4, 13)
        nir <- c(5, 14)
        swir <- c(6, 15)
      }else{
        red <- c(3, 11)
        ndvi <- c(4, 12)
        swir <- c(5, 13)
      }
      
      # Write the min or max hdr for all bands and rns
      writeRaster(brick2, paste0(substr(file_prefix, 1, nchar(file_prefix)-17), "_wo_indices_allBands"), format = "ENVI", overwrite = TRUE, datatype = 'INT2S')
      writeRaster(brick2[[c(red[1], nir[1], swir[1])]], paste0(substr(file_prefix, 1, nchar(file_prefix)-17), "_red_nir_swir1"), format = "ENVI", overwrite = TRUE, datatype = 'INT2S') 
      
      #########  Construct and write the MINMAX images  
      # Stack the min and max images here
      if(i == 1){
        maxmintran_brick <- brick2
      }else if(i > 1){
        # Stack them
        maxmintran_brick <- stack(c(maxmintran_brick, brick2))
      }
      
      # Export
      # Write the minmax allband tif and the two hdrs 
      writeRaster(maxmintran_brick[[c(red[1], nir[1], swir[1], red[2], nir[2], swir[2])]], paste0(substr(file_prefix, 1, nchar(file_prefix) - 25), "maxmintran_NDVI_red_nir_swir1"), format = "ENVI", overwrite = TRUE, datatype = 'INT2S') 
      writeRaster(maxmintran_brick, paste0(substr(file_prefix, 1, nchar(file_prefix)-25), "maxmintran_NDVI_wo_indices_allBands"), format = "ENVI", overwrite = TRUE, datatype = 'INT2S')
      gc() # Garbage collection.
      
      print(paste("Write time:", tock(export_start)))
      print(paste("Total time for image:", tock(one.date.time))) 
    }
    # Then remove any temp files that have been stored on the C drive.The functions specifically targets rasters 
    removeTmpFiles()
    
    # Report time it took to convert all the files
    print(paste("Time to convert tifs to hdr:", tock(tif_2_hdr_files_start)))
  }
}
# Report time for entire code to run 
print(paste("Time for entire code to run:", tock(process_start)))


