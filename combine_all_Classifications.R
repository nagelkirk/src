# Combine all Classification Results

# Description: Takes all the result from the different classifications and puts them into a single table for
#              accuracy assessment.

# Ryan Nagelkirk
# 20181217

#########################################################################################################
#  Working Directory and Functions
#########################################################################################################
# Set working directory
# setwd("E:/Dropbox/Permanent/Grad School/Projects/EleTree/pubs/MethodsPaper")
setwd("~/Dropbox/Permanent/Grad School/Projects/EleTree/pubs/MethodsPaper")

# Source the functions
# source("E:/Dropbox/Permanent/Grad School/src_functions/src_masterfunctions.R")
source("~/Dropbox/Permanent/Grad School/src_functions/src_masterfunctions.R")

#########################################################################################################
#  Load the data 
#########################################################################################################
mcu.sma.raw <- read.csv("./data/MCU_n_SMA/results/20181205/all_MCUnSMA_results_combined.csv")
mesma.raw <- read.csv("./data/MESMA/results/all_MESMAs_appended_to_GT.csv")
regress.raw <- "not done yet"
rf.raw <- "not done yet"

#########################################################################################################
# Combine the tables
#########################################################################################################
# Instead of columns for each classification, I'm rbinding the tables to prevent creating
# too many columns. This way, the table is long (tall) and can be subset in the
# accuracy assessment code by unique classification names. It also allows me to write a simple
# accuracy assessment function (columns won't change)

# Put the tables into identical format
names(mesma.raw)
names(mcu.sma.raw)


#########################################################################################################
#  MESMA  

# Rename the MESMA column names to be more generic (just the last two)
names(mesma.raw)[(ncol(mesma.raw)-1):ncol(mesma.raw)] <- c("class.30", "class.90")

# Subset the columns to only those needed
mesma.short <- mesma.raw[, c("Classification", "park.pxl.num", "Count", "hand.tree.30",
                            "hand.tree.90", "class.30", "class.90")]


#########################################################################################################
#  MCU & SMA

# Take the MCU results and put them in a table identical to the MESMA
mx <- as.data.frame(matrix(NA, nrow(mcu.sma.raw), ncol(mesma.short)))
names(mx) <- names(mesma.short)

# Fill in the values
mx$Classification <- paste0(mcu.sma.raw$park, "_", mcu.sma.raw$name)
mx$Count <- mcu.sma.raw$Count
mx$class.30 <- mcu.sma.raw$first
mx$park.pxl.num <- mcu.sma.raw$park.pxl
mx$hand.tree.30 <- mcu.sma.raw$hand.tree.30
mx$hand.tree.90 <- mcu.sma.raw$hand.tree.90
mcu.short <- mx


#########################################################################################################
#  Regressions



#########################################################################################################
#  RF



#########################################################################################################
#  Combine and Save
#########################################################################################################
all.classifications <- rbind(mesma.short, mcu.short)

# Save
write.csv(all.classifications, file = paste0("./analysis/accuracy_assessments/all_class_results_combined_", YMD(), ".csv"))



