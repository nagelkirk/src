# Helper for GEE Image Download Code
# 20180131
# Ryan Nagelkirk

# Description
# Code takes elephant census data and GEE sourced NDVI values for all possible images and subsets the image dates
# based on ele census years and an NDVI criteria (primary annual min NDVI)

# Purpose
# My GEE code crashes if all possible images are downloaded. This selects only the images needed for woody cover classification.
# That way, only a few images will need to be downloaded, hopefully allowing the GEE code to run.

# Libraries

# Set working directory
setwd("E:/Dropbox/Permanent/Grad School/Projects/EleTree/data")
#setwd("~/Dropbox/Permanent/Grad School/Projects/EleTree/data")


#############  FUNCTIONS  #############################################################################################
# Source the functions
source("E:/Dropbox/Permanent/Grad School/src_functions/src_masterfunctions.R")
#source("~/Dropbox/Permanent/Grad School/src_functions/src_masterfunctions.R")

# Function takes the system.index and finds the LT identifier 
getSatName <- function(x){
  # Get the index
  idx <- x["system.index"]
  
  # Index through the string until you hit an "L", then take that and the next two characters
  suffix_idx <- sub("^[^L]*", "", idx)
  cl_idx <- substr(suffix_idx, 1, 3) 
  
  return(cl_idx)
}


# Function that takes the park name, the df, and adds years 
addYears <- function(df, parkName, years_to_add){
  # Make a blank row for plugging values into
  temp_row <- df[1, ]
  temp_row[] <- NA
  # For each year, add a row with the GEE name and Year
  for(i in 1:length(years_to_add)){
    # Set values
    temp_row$Year <- years_to_add[i]
    temp_row$GEE_park <- parkName
    # bind the row to the end
    df <- rbind(df, temp_row)
  }
  # Remove the first row ()
  df <- df[-1, ]
  return(df)
}

#################  VARIABLES  #####################################################################
summary_stat_boundary <- "twentykmbuffer_shp"  # What boundary do you want to use for calculating NDVI mean? 

# What if did this set of years for each park?: 1984, 1987, 1990, 1994, 1998, 2002, 2006, 2010, 2015 (9) or 1984, 1988, 1992, 1996, 2000, 2004, 2008, 2012, 2016 (9)
image_years <- c(1984, 1988, 1992, 1996, 2000, 2004, 2008, 2012, 2016)


################# LOAD DATA  ######################################################################
# ELE DATA
ele_table <- read.csv("./Census Numbers/MASTER_Compiled_Elephant_Census_Data_raw.csv")

# Compile all the NDVI data into a single table with all parks
NDVI_files <- list.files("./LandsatDownload_inputs/NDVI_MSAVI2_tables/20kmBuffer_Individual_Parks_using_Park_cloud_cover", full.names = T, pattern = ".csv")

# Loop through and combine
for(i in 1:length(NDVI_files)){
  cv <- read.csv(NDVI_files[i], header = T, stringsAsFactors = F)
  if(i == 1){
    df <- cv
  }else{
    df <- rbind(df, cv)
  }
}


########## SUBSET THE DATA  ################################################################################  
# Subset the tables and reassign
NDVI_all <- df[c("system.index", "Boundary", "DATE_ACQUIRED", "NDVI_mean", "Park", "WRS_PATH", "Year", "CLOUD_COVER", "cfmask")]
# Get only the values from the boundary selected in the Variables section up above
NDVI_all_buf <- subset(NDVI_all, Boundary == summary_stat_boundary)
# Landsat 7's scan line corrector (SLC) failed on May 31, 2003, so remove any L7 data past that date
NDVI_all_buf$Satellite <- apply(NDVI_all_buf, 1, getSatName) # create the column and assign satellite names
NDVI_table <- NDVI_all_buf[!(NDVI_all_buf$Satellite == "LE7" & NDVI_all_buf$DATE_ACQUIRED >= "2003-05-31"), ]
NDVI_table <- NDVI_table[!is.na(NDVI_table$NDVI_mean), ]

# Get a list of the protected areas as used in GEE
GEE_pa_list <- unique(NDVI_table$Park)


##########  CREATE ROWS WITH THE WANTED IMAGE YEARS  ################################################################################
# Create a blank row
mx <- matrix(ncol = 2)
df <- as.data.frame(mx)
colnames(df) <- c("GEE_park", "Year")
temp_row <- df

# Used blank row in a loop to build the table with park names 
for(i in 1:length(GEE_pa_list)){
  # For each year, add a row with the GEE name and Year
  for(j in 1:length(image_years)){
    # Set values
    temp_row$Year <- image_years[j]
    temp_row$GEE_park <- GEE_pa_list[i]
    # bind the row to the end
    df <- rbind(df, temp_row)
  }
}

# Remove the first row (it's empty)
ele_table_w_no_extraYears <- df[-1, ]


## GET DATE OF NDVI MIN FOR EACH ROW (Year)
# Loop through the new ele data, subsetting the NDVI data by path and minNDVI for each ele data row and plugging the
# image path, year, park, and date into a new table (so if a park covers three paths, then there will be three rows for each year)

# Duplicate table for final output
dates_table <- ele_table_w_no_extraYears

i <- 48
for(i in 1:nrow(ele_table_w_no_extraYears)){
  # Get one all the NDVI values for one year's census at a park
  sub <- subset(NDVI_table , Park == ele_table_w_no_extraYears$GEE_park[i] & Year == ele_table_w_no_extraYears$Year[i])
  
  # How many paths?
  paths <- unique(NDVI_table[NDVI_table$Park == ele_table_w_no_extraYears$GEE_park[i], ]$WRS_PATH)
  
  # Check whether all paths have imagery
  path_check <- length(unique(sub$WRS_PATH)) == length(paths)
  
  # If there is no imagery for the year of the census, go to the next year. Add years until a year with imagery for all paths is found.
  # The first year change is to add a year, then subtract. So it goes for a later date at first 
  yr_addition <- 1 # This will change until it reaches a level that finds a year with imagery
  while(path_check == FALSE){
    sub <- subset(NDVI_table , Park == ele_table_w_no_extraYears$GEE_park[i] & Year == ele_table_w_no_extraYears$Year[i] - yr_addition)
    # Check whether all paths have imagery
    path_check <- length(unique(sub$WRS_PATH)) == length(paths)
    if(path_check == FALSE){
      sub <- subset(NDVI_table , Park == ele_table_w_no_extraYears$GEE_park[i] & Year == ele_table_w_no_extraYears$Year[i] + yr_addition)
      # Check whether all paths have imagery
      path_check <- length(unique(sub$WRS_PATH)) == length(paths)
    } 
    yr_addition <- yr_addition + 1
  }
  
  # Loop and subset by each path, get the min NDVI, assign it to a row specific to that path. Additional rows are added for subsequent paths on the same year.
  #j <- 1
  for(j in 1:length(paths)){
    # Subset by path
    path_df <- subset(sub, WRS_PATH == paths[j])
    
    # Get the image date with the min NDVI
    minNDVI <- min(path_df$NDVI_mean, na.rm = T)
    img_date <- path_df[which(path_df$NDVI_mean == minNDVI), "DATE_ACQUIRED"]
    img.index <- path_df[which(path_df$NDVI_mean == minNDVI), "Satellite"]
    
    # Assign the date to dates_table
    if(j == 1){
      dates_table$WRS_PATH[i] <- paths[j]
      dates_table$DATE_ACQUIRED[i] <- img_date
      dates_table$Satellite[i] <- img.index
    }else{
      new_row <- dates_table[i, ] # create a new row to be added and assign the path and date values
      new_row$WRS_PATH <- paths[j]
      new_row$DATE_ACQUIRED <- img_date
      new_row$Satellite <- img.index
      dates_table <- rbind(dates_table, new_row) # Put in the new row
    }
  }
}

# Order the table based on PA
LS_dates <- dates_table[order(dates_table$GEE_park), ] 
unique(LS_dates$Satellite)

# Now get only the unique GEE_park, path, and dates. This is so imagery isn't downloaded twice
names(LS_dates)
LS_dates_noReps <- unique(LS_dates[c("GEE_park", "WRS_PATH", "DATE_ACQUIRED", "Satellite")]) # Any more columns than this will start causing replicates (e.g., "Year")

# Export the new table to be loaded in to fusion tables and imported to GEE
write.csv(LS_dates_noReps, paste0("./LandsatDownload_inputs/Image_Date_tables/Landsat_image_4year_span_dates_singleParks_", summary_stat_boundary, "_n_parkCloudCover.csv"), row.names = F)







