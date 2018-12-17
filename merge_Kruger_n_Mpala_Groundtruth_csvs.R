# Code to combine all the Ground truth from Logan's work into one table

#################  KRUGER  ############################################################################################

# Set working directory
setwd("E:/Dropbox/Permanent/Grad School/Projects/EleTree/data/ParkData/Kruger/ground_truth/From_Logan/ground_truth")

# Get the folder names that have the ground truth. Will loop through these
dir.list <- list.files(".", recursive = F, pattern = "section")

# Loop through each folder, pulling the most recent file and adding it to an accumulating table
for(i in 1:length(dir.list)){
  
  # Construct the name
  folder.name <- paste0("./", dir.list[i])
  
  # List the groundtruth files 
  gt.list <- list.files(folder.name, pattern = "groundtruth", full.names = T)
  # Get the most recent one
  most.recent.gt <- sort(gt.list, decreasing = T)[1]
  # Read in the groundtruth table 
  ground.truth <- read.csv(most.recent.gt)
  
  ## Move the cannot.classify values over to the water flag (some tables don't have the column)
  # List the column names
  cNames <- colnames(ground.truth)
  
  # If the column exists, move the 1's over by adding the two columns
  if("cannot.classify" %in% cNames){
    ground.truth$water.flag <- ground.truth$water.flag + ground.truth$cannot.classify
    # Take the cannot.classify column off
    ground.truth <- ground.truth[, -ncol(ground.truth)]
  }
  
  # Make sure the 18 columns of data are the only thing listed (take off the variable number of header columns)
  ground.truth <- ground.truth[, c((ncol(ground.truth)-17) : ncol(ground.truth))]
  print(i)
  print(NROW(ground.truth))
  
  if(i == 1){
    growing.tbl <- ground.truth
  }else{
    growing.tbl <- rbind(growing.tbl, ground.truth)
  }
}

# Renumber the points
growing.tbl$Count <- 1:NROW(growing.tbl)


# write the file
write.csv(growing.tbl, paste0("E:/Dropbox/Permanent/Grad School/Projects/EleTree/data/ParkData/Kruger/ground_truth/wholepark/groundtruth_doneByHand_combineFromLogan.csv" ))



###########  MPALA  ##############################################################################################################

# Set working directory
setwd("E:/Dropbox/Permanent/Grad School/Projects/EleTree/data/ParkData/Mpala/ground_truth")

# Get the folder names that have the ground truth. Will loop through these
dir.list <- list.files(".", recursive = F, pattern = "Section")

# Loop through each folder, pulling the most recent file and adding it to an accumulating table
for(i in 1:length(dir.list)){
  
  # Construct the name
  folder.name <- paste0("./", dir.list[i])
  
  # List the groundtruth files 
  gt.list <- list.files(paste0(folder.name, "/manual_ground_truth_images"), pattern = "groundtruth", full.names = T)
  # Get the most recent one
  most.recent.gt <- sort(gt.list, decreasing = T)[1]
  # Read in the groundtruth table 
  ground.truth <- read.csv(most.recent.gt)
  
  ## Move the cannot.classify values over to the water flag (some tables don't have the column)
  # List the column names
  cNames <- colnames(ground.truth)
  
  # If the column exists, move the 1's over by adding the two columns
  if("cannot.classify" %in% cNames){
    ground.truth$water.flag <- ground.truth$water.flag + ground.truth$cannot.classify
    # Take the cannot.classify column off
    ground.truth <- ground.truth[, -ncol(ground.truth)]
  }
  
  # Make sure the 18 columns of data are the only thing listed (take off the variable number of header columns)
  ground.truth <- ground.truth[, c((ncol(ground.truth)-17) : ncol(ground.truth))]
  
  if(i == 1){
    growing.tbl <- ground.truth
  }else{
    growing.tbl <- rbind(growing.tbl, ground.truth)
  }
}

# Renumber the points
growing.tbl$Count <- 1:NROW(growing.tbl)

# Get rid of any rows that have tree soil and grass all at zero
growing.tbl <- growing.tbl[!apply(growing.tbl[, c(5:7)]==0, 1, all),]

# Get rid of anything with a water flag
growing.tbl <- growing.tbl[!growing.tbl$water.flag == 1, ]


# write the file
write.csv(growing.tbl, paste0("E:/Dropbox/Permanent/Grad School/Projects/EleTree/data/ParkData/Mpala/ground_truth/wholepark/groundtruth_doneByHand_combined.csv" ))





