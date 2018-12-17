# Combine all SMA and MCU results from the GEE classification

# Ryan Nagelkirk
# 20181213

# Description:
# Takes all the csvs from the MCU and SMA GEE code and combines them into one large csv

# Set the working directory
setwd("~/Dropbox/Permanent/Grad School/Projects/EleTree/pubs/MethodsPaper/data/MCU_n_SMA/Extracted_Values_from_GEE/20181205")

# We'll need the park list
park.names <- c("Chobe", "Kruger", "Limpopo", "Mpala", "Murchison", "North_Luangwa", "QWE", "Ruaha", "Selous", "Serengeti", "South_Luangwa")

# Get the list of files
f.list <- list.files(pattern = "GEE_AutoMCUnSMA_pointVals")

# Loop through the files and combine
# Just get the necessary columns and the rest can be added back later
# Doing that because of Mpala and Kruger once again having different numbers of columns
# It's just easier to export what we need and move on
for(i in 1:length(f.list)){
  # Get the file
  f.iteration <- read.csv(f.list[i])
  
  # Get just the columns needed 
  f.short <- f.iteration[, c("system.index", "Count", "first", "name")]
  # Get the "pxl" column as well. However, that column has different names, so use "px" as the pattern
  f.short$park.pxl <- f.iteration[, grep("px", names(f.iteration))]
  
  # Then add the park names by lookin at the first three letters and finding the name in the park list
  # that matches
  park.name3 <- substr(f.list[i], 1, 3)
  park.name <- park.names[grep(park.name3, substr(park.names, 1, 3))]
  f.short$park <- park.name
  
  #Add it to a growing tbl
  if(i == 1){
    growing.tbl <- f.short
  }else{
    growing.tbl <- rbind(growing.tbl, f.short)
  }
}

# Write the file
write.csv(growing.tbl, file = "./all_MCUnSMA_results_combined.csv")

# Ruaha, Selous, 
names(f.iteration)
