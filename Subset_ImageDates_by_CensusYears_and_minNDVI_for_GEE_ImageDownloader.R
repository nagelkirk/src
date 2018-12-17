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

# Function that takes the park name, the df, and years wanted and removes unwanted years 
subsetYears <- function(df, parkName, years_to_keep){
  # Create empty list of rows to delete
  rows_to_remove <- c()
  # If GEE name == park name, see if it has one of the years. If not, record that row number (can't let the original get smaller)
  for(i in 1:nrow(df)){
    if(df$GEE_park[i] == parkName){
      # If the year is NOT in the list of years to keep, record the row
      if(!df$Year[i] %in% years_to_keep){
        rows_to_remove <- c(rows_to_remove, i)
      }
    }
  }
  # With the loop complete, subset the table, getting rid of the rows identified in rows_to_remove
  return(df[-rows_to_remove, ])
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
  return(df)
}

#################  VARIABLES  #####################################################################
summary_stat_boundary <- "twentykmbuffer_shp"  # What boundary do you want to use for calculating NDVI mean? 


# What if did this set of years for each park?: 1984, 1987, 1990, 1994, 1998, 2002, 2006, 2010, 2015 (9) or 1984, 1988, 1992, 1996, 2000, 2004, 2008, 2012, 2016 (9)
kruger_years_kept <- c(1977, 1987, 1995, 1998, 2002, 2006, 2010, 2015) # Kruger and some of the others have so many years of data, have to pick years
chobe_years_kept <- c(2014, 2010, 2006, 2002, 1999, 1994)
limpopo_years_kept <- c(2002, 2006, 2010, 2014)
mara_years_kept <- c(1977, 1984, 1987, 1992, 1997, 2002, 2006, 2010, 2014)
ruaha_years_kept <- c(1984, 1993, 1996, 2002, 2006, 2009, 2014)
selous_years_kept <- c(1977, 1986, 1994, 1998, 2002, 2006, 2009, 2014)
serengeti_years_kept <- c(1977, 1987, 1994, 1998, 2003, 2008, 2014)  
qwe_years_kept <- c(1977, 1987, 1992, 1995, 1998, 2002, 2006, 2012, 2014) 
murchison_years_kept <- c(1977, 1987, 1992, 1999, 2002, 2005, 2010, 2014)
Lu_S_years_kept <- c(1987, 1996, 2002, 2009, 2012, 2015)
Lu_N_years_kept <- c(1987, 1993, 1996, 2001, 2003, 2007, 2011, 2015)
mpala_years_kept <- c(1987, 1992, 1996, 2002, 2012, 2014) 

# Add years for Limpopo, which only goes back to 2002 for censuses
kruger_years_add <- c(1991)
chobe_years_add <- c(1992, 1988, 1984)
limpopo_years_add <- c(1999, 1996, 1992, 1988, 1984)
ruaha_years_add <- c(1988, 1999)
selous_years_add <- c(1990)
serengeti_years_add <- c(2011, 1990)
murchison_years_add <- c(1995)
Lu_S_years_add <- c(1984, 1992, 2005)
Lu_N_years_add <- c(1984)
mpala_years_added <- c(1984, 1999, 2007) 





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
Stratum_list <- unique(ele_table$Stratum)
Stratum_list[order(Stratum_list)]
Ecosystem_list <- unique(ele_table$Ecosystem)

# WHO ACTUALLY HAS AN ENTIRE ECOSYSTEM VALUE?
Eco_Strat <- unique(ele_table[c("Ecosystem", "Stratum")])
Eco_Strat[Eco_Strat$Stratum == "ENTIRE ECOSYSTEM", ]

# Go through and assign the GEE names to the Census names 
# Since some rows will satisfy two parks, do an if that adds a row if a park name already exists in that GEE_park cell 
for(i in 1:nrow(ele_table)){
  if(ele_table$Ecosystem[i] == "Northern Botswana" & any(ele_table$Stratum[i] == "Chobe district" | ele_table$Stratum[i] == "Chobe NP" | ele_table$Stratum[i] == "ENTIRE ECOSYSTEM")){
    ele_table$GEE_park[i] <- "Chobe"
    ele_table$eco_strat[i] <- paste0(ele_table$Ecosystem[i], "_", ele_table$Stratum[i])
  }else if(ele_table$Ecosystem[i] == "Kruger NP" & any(ele_table$Stratum[i] == "Kruger NP" | ele_table$Stratum[i] == "ENTIRE ECOSYSTEM")){
    ele_table$GEE_park[i] <- "Kruger"
    ele_table$eco_strat[i] <- paste0(ele_table$Ecosystem[i], "_", ele_table$Stratum[i])
  }else if(ele_table$Ecosystem[i] == "Limpopo NP" & any(ele_table$Stratum[i] == "Limpopo NP" | ele_table$Stratum[i] == "ENTIRE ECOSYSTEM")){  
    ele_table$GEE_park[i] <- "Limpopo"
    ele_table$eco_strat[i] <- paste0(ele_table$Ecosystem[i], "_", ele_table$Stratum[i])
  }else if(ele_table$Ecosystem[i] == "Laikipia-Samburu" & any(ele_table$Stratum[i] == "Laikipia-Samburu" | ele_table$Stratum[i] == "Laikipia District" | ele_table$Stratum[i] == "ENTIRE ECOSYSTEM")){
    ele_table$GEE_park[i] <- "Mpala"
    ele_table$eco_strat[i] <- paste0(ele_table$Ecosystem[i], "_", ele_table$Stratum[i])
  }else if(ele_table$Ecosystem[i] == "Murchison Falls Protected Area" & any(ele_table$Stratum[i] == "Murchison Falls PA" | ele_table$Stratum[i] == "Murchison South NP" | ele_table$Stratum[i] == "ENTIRE ECOSYSTEM")){
    ele_table$GEE_park[i] <- "Murchison"
    ele_table$eco_strat[i] <- paste0(ele_table$Ecosystem[i], "_", ele_table$Stratum[i])
  }else if(ele_table$Ecosystem[i] == "Luangwa" & any(ele_table$Stratum[i] == "Luangwa North NP" | ele_table$Stratum[i] == "ENTIRE ECOSYSTEM" | ele_table$Stratum[i] == "Luangwa North and South NP")){
    ele_table$GEE_park[i] <- "North_Luangwa"
    ele_table$eco_strat[i] <- paste0(ele_table$Ecosystem[i], "_", ele_table$Stratum[i])
  }else if(ele_table$Ecosystem[i] == "Luangwa" & any(ele_table$Stratum[i] == "Luangwa South NP" | ele_table$Stratum[i] == "ENTIRE ECOSYSTEM" | ele_table$Stratum[i] == "Luangwa North and South NP")){
    ele_table$GEE_park[i] <- "South_Luangwa"
    ele_table$eco_strat[i] <- paste0(ele_table$Ecosystem[i], "_", ele_table$Stratum[i])
  }else if(ele_table$Ecosystem[i] == "Greater Virunga Landscape" & any(ele_table$Stratum[i] == "Queen Elizabeth NP" | ele_table$Stratum[i] == "ENTIRE ECOSYSTEM")){
    ele_table$GEE_park[i] <- "QWE"
    ele_table$eco_strat[i] <- paste0(ele_table$Ecosystem[i], "_", ele_table$Stratum[i])
  }else if(ele_table$Ecosystem[i] == "Selous-Mikumi" & any(ele_table$Stratum[i] == "Selous-Mikumi" | ele_table$Stratum[i] == "Selous GR" | ele_table$Stratum[i] == "ENTIRE ECOSYSTEM")){
    ele_table$GEE_park[i] <- "Selous"
    ele_table$eco_strat[i] <- paste0(ele_table$Ecosystem[i], "_", ele_table$Stratum[i])
  }else if(ele_table$Ecosystem[i] == "Serengeti" & any(ele_table$Stratum[i] == "Serengeti" | ele_table$Stratum[i] == "ENTIRE ECOSYSTEM")){ # Is there a Serengeti-Mara stratum?
    ele_table$GEE_park[i] <- "Serengeti"
    ele_table$eco_strat[i] <- paste0(ele_table$Ecosystem[i], "_", ele_table$Stratum[i])
  }else if(ele_table$Ecosystem[i] == "Masai Mara" & any(ele_table$Stratum[i] == "Masai Mara R" | ele_table$Stratum[i] == "ENTIRE ECOSYSTEM")){
    ele_table$GEE_park[i] <- "Mara"
    ele_table$eco_strat[i] <- paste0(ele_table$Ecosystem[i], "_", ele_table$Stratum[i])
  }else if(ele_table$Ecosystem[i] == "Ruaha-Rungwa" & any(ele_table$Stratum[i] == "Ruaha-Rungwa" | ele_table$Stratum[i] == "Ruaha NP" | ele_table$Stratum[i] == "ENTIRE ECOSYSTEM")) {
    ele_table$GEE_park[i] <- "Ruaha"
    ele_table$eco_strat[i] <- paste0(ele_table$Ecosystem[i], "_", ele_table$Stratum[i])
  }else{
    ele_table$GEE_park[i] <- NA
    ele_table$eco_strat[i] <- NA
  }
  # Below covers the contingencies of two parks being within a stratum (N and S Luangwa in Luangwa ENTIRE ECOSYSTEM)
  # Creates another row for the other park
  if(ele_table$Ecosystem[i] == "Luangwa" & any(ele_table$Stratum[i] == "ENTIRE ECOSYSTEM" | ele_table$Stratum[i] == "Luangwa North and South NP")){
    temp_row <- ele_table[i, ]
    temp_row$GEE_park <- "South_Luangwa"
    temp_row$eco_strat <- paste0(temp_row$Ecosystem, "_", temp_row$Stratum)
    ele_table <- rbind(ele_table, temp_row)
  }
}


## Subset the ele data to only the parks and columns you're interested in 
ele_table_w_GEEnames <- ele_table[!is.na(ele_table$GEE_park), c("Ecosystem", "Stratum", "Year", "Elephants", "GEE_park", "eco_strat", "Stratum.area", "Total.area.surveyed", "My.source")]


#########  REMOVE EXTRA UNWANTED YEARS  #############################################################################################################
# Remove years in parks that have so much data that they will generate too many scenes for downloading (e.g., Kruger)
ele_table_w_no_extraYears <- subsetYears(df = ele_table_w_GEEnames, parkName = "Kruger", years_to_keep = kruger_years_kept)
ele_table_w_no_extraYears <- subsetYears(df = ele_table_w_no_extraYears, parkName = "Mara", years_to_keep = mara_years_kept)
ele_table_w_no_extraYears <- subsetYears(df = ele_table_w_no_extraYears, parkName = "Chobe", years_to_keep = chobe_years_kept)
ele_table_w_no_extraYears <- subsetYears(df = ele_table_w_no_extraYears, parkName = "QWE", years_to_keep = qwe_years_kept)
ele_table_w_no_extraYears <- subsetYears(df = ele_table_w_no_extraYears, parkName = "Ruaha", years_to_keep = ruaha_years_kept)
ele_table_w_no_extraYears <- subsetYears(df = ele_table_w_no_extraYears, parkName = "Selous", years_to_keep = selous_years_kept)
ele_table_w_no_extraYears <- subsetYears(df = ele_table_w_no_extraYears, parkName = "Serengeti", years_to_keep = serengeti_years_kept)
ele_table_w_no_extraYears <- subsetYears(df = ele_table_w_no_extraYears, parkName = "Murchison", years_to_keep = murchison_years_kept)
ele_table_w_no_extraYears <- subsetYears(df = ele_table_w_no_extraYears, parkName = "North_Luangwa", years_to_keep = Lu_N_years_kept)
ele_table_w_no_extraYears <- subsetYears(df = ele_table_w_no_extraYears, parkName = "South_Luangwa", years_to_keep = Lu_S_years_kept)


##########  ADD YEARS WITHOUT ELE DATA BUT THAT ADD TO WC DATA  ######################################################################################
ele_table_w_no_extraYears <- addYears(df = ele_table_w_no_extraYears, parkName = "Kruger", years_to_add = kruger_years_add)
ele_table_w_no_extraYears <- addYears(df = ele_table_w_no_extraYears, parkName = "Chobe", years_to_add = chobe_years_add)
ele_table_w_no_extraYears <- addYears(df = ele_table_w_no_extraYears, parkName = "Mpala", years_to_add = mpala_years_added)
ele_table_w_no_extraYears <- addYears(df = ele_table_w_no_extraYears, parkName = "Limpopo", years_to_add = limpopo_years_add)
ele_table_w_no_extraYears <- addYears(df = ele_table_w_no_extraYears, parkName = "Ruaha", years_to_add = ruaha_years_add)
ele_table_w_no_extraYears <- addYears(df = ele_table_w_no_extraYears, parkName = "Selous", years_to_add = selous_years_add)
ele_table_w_no_extraYears <- addYears(df = ele_table_w_no_extraYears, parkName = "Serengeti", years_to_add = serengeti_years_add)
ele_table_w_no_extraYears <- addYears(df = ele_table_w_no_extraYears, parkName = "Murchison", years_to_add = murchison_years_add)
ele_table_w_no_extraYears <- addYears(df = ele_table_w_no_extraYears, parkName = "North_Luangwa", years_to_add = Lu_N_years_add)
ele_table_w_no_extraYears <- addYears(df = ele_table_w_no_extraYears, parkName = "South_Luangwa", years_to_add = Lu_S_years_add)






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
    sub <- subset(NDVI_table , Park == ele_table_w_no_extraYears$GEE_park[i] & Year == ele_table_w_no_extraYears$Year[i] + yr_addition)
    # Check whether all paths have imagery
    path_check <- length(unique(sub$WRS_PATH)) == length(paths)
    if(path_check == FALSE){
      sub <- subset(NDVI_table , Park == ele_table_w_no_extraYears$GEE_park[i] & Year == ele_table_w_no_extraYears$Year[i] - yr_addition)
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
LS_dates <- dates_table[order(dates_table$Stratum), ] 
unique(LS_dates$Satellite)

# Now get only the unique GEE_park, path, and dates. This is so imagery isn't downloaded twice
names(LS_dates)
LS_dates_noReps <- unique(LS_dates[c("GEE_park", "WRS_PATH", "DATE_ACQUIRED", "Satellite")]) # Any more columns than this will start causing replicates (e.g., "Year")


# Export the new table to be loaded in to fusion tables and imported to GEE
write.csv(LS_dates_noReps, paste0("./LandsatDownload_inputs/Image_Date_tables/Landsat_image_dates_toExtract_using_singleParks_", summary_stat_boundary, "_n_parkCloudCover.csv"), row.names = F)


length(LS_dates$DATE_ACQUIRED[LS_dates_noReps$Satellite == "LC8"])
nrow(LS_dates_noReps)





