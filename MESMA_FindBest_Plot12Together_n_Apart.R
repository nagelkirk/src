# MESMA Accuracy Sorter

# Ryan Nagelkirk
# 2018/09/12



#############
### THIS CODE HAS BEEN MODIFIED TO MAKE PLOTS FOR FORESTSAT CONFERENCE POSTER  ###############
#########




# Description
# 1. Finds the MESMAs that had the highest accuracy and best slope for each of the parks
# 2. Evaluates TG, TGS, min, max, minmax, etc. to see what did better overall in each park
# 3. Evaluates TG, TGS, min, max, minmax, etc. to see what did better overall in all parks



# NOTES:

# So far, the data looks best when all the values outside of 0-100% are discarded. Then the RMSE, MAE and the VEcv are the best






#############################  LOAD LIBRARIES  ######################################################################
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
  library(dplyr)
}



#################################  Functions  ##################################################################
# Source the functions
source("E:/Dropbox/Permanent/Grad School/src_functions/src_masterfunctions.R")
# source("~/Dropbox/Permanent/Grad School/src_functions/src_masterfunctions.R")



###############################  DIRECTORIES  ##############################################################
# Set working directory
setwd("E:/Dropbox/Permanent/Grad School/Projects/EleTree/pubs/MethodsPaper/analysis/accuracy_assessments")
# setwd("~/Dropbox/Permanent/Grad School/Projects/EleTree/pubs/MethodsPaper/analysis/accuracy_assessments")

# Folder with MESMA accuracy assessment results 
hdr.eval.folder <- "./hdr_eval_dfs"

##############################  Switches  #################################################################
save.plot <- F


#########################  Get best MESMAS  ###############################################################
# List the parks using the file names
csv.names <- list.files(hdr.eval.folder, pattern = "HDR_eval", full.names = T)
# Pull out the names
park.names <- gsub('^.*_eval_\\s*|\\s*_201.*$', '', csv.names)
park.names <- unique(park.names)

# For each park, pull in the file and determine the best MESMAs within it
for(i in 1:length(park.names)){
  # Get the files for the park
  hdr.evals <- grep(pattern = park.names[i], csv.names, value = T)
  # There might be multiple, so get the most recent
  hdr.eval <- sort(hdr.evals, decreasing = T)[1]
  
  # Read in the file
  hdr.eval.df <- read.csv(hdr.eval)
  hdr.eval.df <- subset(hdr.eval.df, select = c(hdr.name, RMSE, adj.r.sq, slope, y.int, p.value, num.points.plotted, VEcv))
  
  
  # One way to rank these would be to normalize/standardize each by a quasi legitimate method
  # Create scores that are: 1. how far the slope was from 1 and 2. how far the r squared is from 1
  # Then weight those values if you want and add them together into a single score, then sort by that score 
  hdr.eval.df$slope.from.one <- abs(1 - hdr.eval.df$slope)
  hdr.eval.df$r2.from.one <- abs(0.5 - hdr.eval.df$adj.r.sq)
  
  # Add the values to make a new score
  hdr.eval.df$accuracy.score <- hdr.eval.df$slope.from.one + hdr.eval.df$r2.from.one
  
  # Sort based on the new score
  sorted.eval <- hdr.eval.df[order(hdr.eval.df$accuracy.score), ] # I changed this to slope for looking at all the rows manually. Was accuracy.score before
  
  # 20181017 - I added VEcv to this code and the MESMA_vs_Groundtruth_all12Parks_2080820 code
  # That's a better value to sort by
  sorted.eval <- sorted.eval[order(sorted.eval$VEcv, decreasing = T), ]
  
  
  # Put the best row in a new df that will be growing
  if(i == 1){
    all.park.bests <- sorted.eval[c(1:5), ]
  }else{
    all.park.bests <- rbind(all.park.bests, sorted.eval[c(1:5), ], make.row.names = F)
  }
  
}

# Save the table
write.csv(all.park.bests, "E:/Dropbox/Permanent/Grad School/Projects/EleTree/pubs/MethodsPaper/analysis/accuracy_assessments/best_MESMA_methods.csv")


# Record the names of the best MESMAs for each park here
best.mesmas <- c("Chobe_2_LC08_174_20161016_20160423_minmax_NDVI_red_nir_swir1_20161016_174_2image_roi_HML_TGS_vs_GT",
                 "Kruger_1_LC08_168_20161006_20160328_min_NDVI_wo_indices_allBands_20161006_168_2image_roi_spectLib_HML_all_TGS_vs_GT_stchd",
                 "Limpopo_1_LC08_168_20161022_20160413_min_NDVI_red_nir_swir1_20161022_168_2image_roi_spectLib_HML_TG_vs_GT",
                 "Mpala_LC08_168_20161022_20160531_minmax_NDVI_wo_indices_allBands_20161022_2image_roi_spectLib_EMC_TG_vs_GT", 
                 "Murchison_LC08_172_20160120_20160527_minmax_NDVI_red_nir_swir1_20160120_172_2image_roi_HML_TG_vs_GT",
                 "North_Luangwa_2_LC08_170_20161020_20160310_max_NDVI_wo_indices_allBands_20161020_170_2image_roi_HML_TG_vs_GT",
                 "QWE_LC08_173_20160127_20161126_max_NDVI_wo_indices_allBands_20160127_173_2image_roi_EMC_TG_vs_GT",
                 "Ruaha_1_LC08_168_20161022_20160312_max_NDVI_wo_indices_allBands_20161022_168_2image_roi_EMC_TG_vs_GT_stchd",
                 "Selous_1_LC08_166_20161008_20160704_minmax_NDVI_wo_indices_allBands_20161008_166_2image_roi_spectLib_EMC_TG_vs_GT_stchd",
                 "Serengeti_Mara_1_LC08_169_20160911_20160216_max_NDVI_wo_indices_allBands_20160911_169_2image_roi_EMC_TG_vs_GT_stchd",
                 "South_Luangwa_LC08_170_20161020_20160326_min_NDVI_red_nir_swir1_20161020_170_2image_roi_EMC_TGS_vs_GT") # there were serengeti and Murchison ones that were close seconds and might look better (or maybe mpala?)

close.seconds <- c("South_Luangwa_LC08_170_20161020_20160326_minmax_NDVI_wo_indices_allBands_20161020_170_2image_roi_EMC_TG_vs_GT",
                   "Serengeti_Mara_1_LC08_169_20160911_20160216_min_NDVI_red_nir_swir1_20160911_169_2image_roi_HML_TG_vs_GT_stchd",
                   "Mpala_LC08_168_20161022_20160531_max_NDVI_red_nir_swir1_20161022_2image_roi_spectLib_HML_all_TG_vs_GT",
                   "Murchison_LC08_172_20160120_20160527_max_NDVI_wo_indices_allBands_20160120_172_2image_roi_EMC_TG_vs_GT")


best.VEcv.mesmas <- c("Chobe_2_LC08_174_20161016_20160423_max_NDVI_wo_indices_allBands_20161016_174_2image_roi_HML_TGS_vs_GT",
                      "Kruger_1_LC08_168_20161006_20160328_max_NDVI_red_nir_swir1_20161006_168_2image_roi_spectLib_EMC_TGS_vs_GT_stchd",
                      "Limpopo_1_LC08_168_20161022_20160413_minmax_NDVI_wo_indices_allBands_20161022_168_2image_roi_spectLib_HML_TG_vs_GT",
                      "Mpala_LC08_168_20161022_20160531_min_NDVI_wo_indices_allBands_20161022_2image_roi_spectLib_EMC_TGS_vs_GT", 
                      "Murchison_LC08_172_20160120_20160527_min_NDVI_red_nir_swir1_20160120_172_2image_roi_EMC_TG_vs_GT",
                      "North_Luangwa_2_LC08_170_20161020_20160310_max_NDVI_red_nir_swir1_20161020_170_2image_roi_EMC_TG_vs_GT",
                      "QWE_LC08_173_20160127_20161126_minmax_NDVI_wo_indices_allBands_20160127_173_2image_roi_EMC_TG_vs_GT",
                      "Ruaha_1_LC08_168_20161022_20160312_minmax_NDVI_red_nir_swir1_20161022_168_2image_roi_EMC_TG_vs_GT_stchd",
                      "Selous_1_LC08_166_20161008_20160704_min_NDVI_wo_indices_allBands_20161008_166_2image_roi_spectLib_EMC_TGS_vs_GT_stchd",
                      "Serengeti_Mara_1_LC08_169_20160911_20160216_minmax_NDVI_red_nir_swir1_20160911_169_2image_roi_EMC_TGS_vs_GT_stchd",
                      "South_Luangwa_LC08_170_20161020_20160326_min_NDVI_red_nir_swir1_20161020_170_2image_roi_EMC_TGS_vs_GT")

# How many of the mesmas match? Only 1??
best.mesmas == best.VEcv.mesmas


##############################  Join all the Best MESMAS into one table  ###############################################
# Loop through and combine the best.mesmas
# First, create a table to plug hdr.name, RMSE, slope, 1-slope, Adj r-sq, yint into
mx <- matrix(0, nrow = 1, ncol = 7) # set number of rows to zero so the plot1to1 function can add rows as it goes
c.names <- c("hdr.name", "RMSE", "adj.r.sq", "slope", "y-int", "p.value", "num.points.plotted")
colnames(mx) <- c.names
hdr.eval.df <- as.data.frame(mx)

# Start a vE and RMSE table
mx <- matrix(0, nrow = 22, ncol = 3) 
c.names <- c("name", "RMSE", "VEcv")
colnames(mx) <- c.names
VE.RMSE.tbl <- as.data.frame(mx)


# Now loop through, getting the y-intercepts as you go
for(j in 1:length(best.mesmas)){
  # read in the table
  current.table <- read.csv(paste0("./MESMA_vs_GT_ptBYpt/", best.mesmas[j], ".csv"))
  
  # Get just the hand.tree and 30 meter data
  hand.tree.col <- grep("hand.tree", names(current.table))[1]
  MESMA.col <- grep("MESMA.30", names(current.table))
  # Get those suckas
  current.table <- current.table[ , c(hand.tree.col, MESMA.col)]
  # Rename columns so they can bind
  names(current.table) <- c("hand.tree", "MESMA")
  
  # plot
  hdr.eval.df <- plot1to1(x = current.table$hand.tree, y = current.table$MESMA, xlab = "Observed % Tree Cover", ylab = "MESMA % Tree Cover",
                          main = paste0("Percent Tree Cover: Observed (30m) vs. MESMA (30m) ", j), ylim = c(-100, 300), add_one2one_line = TRUE, add.reg.line = TRUE,
                          save.plot = save.plot, save.plot.name = paste0("./images/", best.mesmas[j], "best.png"), hdr.table = hdr.eval.df)

  # Record the y intercept in a new column
  current.table$y_int <- hdr.eval.df$`y-int`[grep(best.mesmas[j], hdr.eval.df$hdr.name)]

  # Adjust for y intercept
  current.table$y_adjusted_MESMA <- as.numeric(current.table$MESMA) - as.numeric(current.table$y_int)

  # Plot that and add to the table
  hdr.eval.df <- plot1to1(x = current.table$hand.tree, y = current.table$y_adjusted_MESMA, xlab = "Observed % Tree Cover", ylab = "MESMA % Tree Cover",
                          main = "Percent Tree Cover: Observed (30m) vs. MESMA (30m)", ylim = c(-100, 300), add_one2one_line = TRUE, add.reg.line = TRUE,
                          save.plot = save.plot, save.plot.name = paste0("./images/", best.mesmas[j], "best_y_Adjusted.png"), hdr.table = hdr.eval.df)


  
  # Make a new table just for plotting
  current.plot.tbl <- current.table 
  
  # Adjust all the values so that if a value is over 1 or under 0, it goes to those values
  current.plot.tbl$MESMA[current.plot.tbl$MESMA > 1] <- 1
  current.plot.tbl$MESMA[current.plot.tbl$MESMA < 0] <- 0
  
  # Filter out any NAs or other incomplete cases 
  current.plot.tbl <- current.plot.tbl[complete.cases(current.plot.tbl), ]
  
  # Plot for the ForestSAT poster
  # Define x and y first
  xvar <- current.plot.tbl$hand.tree*100
  yvar <- current.plot.tbl$MESMA*100
  
  # Get the RMSE
  rmse <- sqrt(mean((yvar - xvar)^2 , na.rm = TRUE))
  
  # Get the VEcv
  VE <- VEcv(observed = xvar, predicted = yvar)
  
  # Add to table
  VE.RMSE.tbl[j, ] <- c(as.character(best.mesmas[j]), round(rmse, 0), round(VE, 0))
  
  
  {
    # png(paste0("E:/Dropbox/Permanent/Grad School/Meetings/ForestSAT/2018/Poster/Images/", gsub(".csv", "", best.mesmas[j]), "_ForestSAT.png"), 
    #     width = 1200, height = 1200)
    # 
    # plot(x = yvar,
    #      y = xvar,
    #      main = paste0(best.mesmas[j], " \nRMSE = ", round(rmse, 0), "% VEcv = ", round(VE, 0), "%"),
    #      xlab = "Reference % Woody Cover", ylab = "MESMA % Woody Cover",
    #      xlim = c(0, 100), ylim = c(0, 100), cex = 4, lwd = 3, pch = 19)
    # 
    # # Add 1:1 line
    # abline(0,1, lwd = 2) # Add the 1:1 line
    # 
    # # Add a regression line
    # modl <- lm(xvar ~ yvar) # calc regression
    # abline(modl, lty = 2, lwd = 2)
    # # # Add another line, shifted up to 1:1 line
    # # mdl2 <- modl
    # # mdl2$coefficients[[1]] <- 0
    # # abline(mdl2, lty = 2, col = "red")
    # 
    # dev.off()
    
    
  }
  
  
  
  # Combine
  if(j == 1){
    growing.table <- current.table
  }else{
    growing.table <- rbind(growing.table, current.table)
  }
}


##########################  Combine the VCFs from each Park  ##################################################
# List the vcf files
vcf.files <- list.files("./MESMA_vs_GT_ptBYpt", pattern = "VCF")

# Loop through, load, y adjust and plot
# Now loop through, getting the y-intercepts as you go
for(k in 1:length(vcf.files)){
  # read in the table
  current.table.vcf <- read.csv(paste0("./MESMA_vs_GT_ptBYpt/", vcf.files[k]))
  
  # Get just the hand.tree and 30 meter data
  hand.tree.col <- grep("hand.tree", names(current.table.vcf))[1]
  MESMA.col <- grep("MESMA.30", names(current.table.vcf))
  # Get those suckas
  current.table.vcf <- current.table.vcf[ , c(hand.tree.col, MESMA.col)]
  # Rename columns so they can bind
  names(current.table.vcf) <- c("hand.tree", "VCF")
  
  # plot
  hdr.eval.df <- plot1to1(x = current.table.vcf$hand.tree, y = current.table.vcf$VCF, xlab = "Observed % Tree Cover", ylab = "VCF % Tree Cover",
                          main = "Percent Tree Cover: Observed (30m) vs. VCF (30m)", ylim = c(-100, 300), add_one2one_line = TRUE, add.reg.line = TRUE,
                          save.plot = save.plot, save.plot.name = paste0("./images/", gsub(".csv", "", vcf.files[k]), "best.png"), hdr.table = hdr.eval.df)


  # Record the y intercept in a new column
  current.table.vcf$y_int <- hdr.eval.df$`y-int`[grep(paste0(gsub(".csv", "", vcf.files[k]), "best.png"), hdr.eval.df$hdr.name)]

  # Adjust for y intercept
  current.table.vcf$y_adjusted_VCF <- as.numeric(current.table.vcf$VCF) - as.numeric(current.table.vcf$y_int)

  # Plot that and add to the table
  hdr.eval.df <- plot1to1(x = current.table.vcf$hand.tree, y = current.table.vcf$y_adjusted_VCF, xlab = "Observed % Tree Cover", ylab = "VCF % Tree Cover",
                          main = "Percent Tree Cover: Observed (30m) vs. VCF (30m)", ylim = c(-100, 300), add_one2one_line = TRUE, add.reg.line = TRUE,
                          save.plot = save.plot, save.plot.name = paste0("./images/", gsub(".csv", "", vcf.files[k]), "best_y_Adjusted.png"), hdr.table = hdr.eval.df)

  
  # Filter out any NAs or other incomplete cases 
  current.table.vcf <- current.table.vcf[complete.cases(current.table.vcf), ]
  

  # Plot for the ForestSAT poster
  # Define x and y first
  xvar <- current.table.vcf$hand.tree*100
  yvar <- current.table.vcf$VCF*100

  # Get the RMSE
  rmse <- sqrt(mean((yvar - xvar)^2 , na.rm = TRUE))

  # Get the VEcv
  VE <- VEcv(observed = xvar, predicted = yvar)

  # Add to table
  VE.RMSE.tbl[k + 11, ] <- c(as.character(vcf.files[k]), round(rmse, 0), round(VE, 0))
  
  {
    # png(paste0("E:/Dropbox/Permanent/Grad School/Meetings/ForestSAT/2018/Poster/Images/", gsub(".csv", "", vcf.files[k]), "_ForestSAT.png"), 
    #     width = 1200, height = 1200)
    # 
    # 
    # plot(x = yvar,
    #      y = xvar,
    #      main = paste0(vcf.files[k], " \nRMSE = ", round(rmse, 0), "% VEcv = ", round(VE, 0), "%"),
    #      xlab = "Reference % Woody Cover", ylab = "VCF % Tree Cover",
    #      xlim = c(0, 100), ylim = c(0, 100), cex = 4, lwd = 3, col = "#FF5722", pch = 19)
    # 
    # # Add 1:1 line
    # abline(0,1, lwd = 2) # Add the 1:1 line
    # 
    # # Add a regression line
    # modl <- lm(xvar ~ yvar) # calc regression
    # abline(modl, lty = 2, lwd = 2)
    # # # Add another line, shifted up to 1:1 line
    # # mdl2 <- modl
    # # mdl2$coefficients[[1]] <- 0
    # # abline(mdl2, lty = 2, col = "red")
    # 
    # dev.off()
    
  }
  
  
  
  
  # Combine
  if(k == 1){
    growing.vcf.table <- current.table.vcf
  }else{
    growing.vcf.table <- rbind(growing.vcf.table, current.table.vcf)
  }
  
}


# Check the correlation 
mesma.cor <- cor(current.table$hand.tree, current.table$MESMA, method = "spearman", use = "complete.obs")
vcf.cor <- cor(current.table.vcf$hand.tree, current.table.vcf$VCF, method = "spearman", use = "complete.obs")



# Plot VE and RMSE from MESMA and VCF against each other (MESMA on x-axis)
plot(x = VE.RMSE.tbl$VEcv[c(1:11)], y = VE.RMSE.tbl$VEcv[c(12:22)], xlab = "MESMA", ylab = "VCF")
abline(0, 1)
plot(x = VE.RMSE.tbl$RMSE[c(1:11)], y = VE.RMSE.tbl$RMSE[c(12:22)], xlab = "MESMA", ylab = "VCF")
abline(0, 1)

# Clean the data to have just the complete cases
growing.table <- growing.table[complete.cases(growing.table), ]
growing.vcf.table <- growing.vcf.table[complete.cases(growing.vcf.table), ]

# Create a quality assurance column that will have data that has been screened
growing.table$QAQC_MESMA <- growing.table$MESMA
growing.vcf.table$QAQC_VCF <- growing.vcf.table$VCF

new.high.value <- 1
new.low.value <- 0

{
  # Adjust all the values so that if a value is over 1 or under 0, it goes to those values
  growing.table$QAQC_MESMA[growing.table$MESMA > 1] <- new.high.value
  growing.table$QAQC_MESMA[growing.table$MESMA < 0] <- new.low.value
  growing.table$y_adjusted_MESMA[growing.table$y_adjusted_MESMA > 1] <- new.high.value
  growing.table$y_adjusted_MESMA[growing.table$y_adjusted_MESMA < 0] <- new.low.value
  growing.vcf.table$QAQC_VCF[growing.vcf.table$VCF > 1] <- new.high.value
  growing.vcf.table$QAQC_VCF[growing.vcf.table$VCF < 0] <- new.low.value
  growing.vcf.table$y_adjusted_VCF[growing.vcf.table$y_adjusted_VCF > 1] <- new.high.value
  growing.vcf.table$y_adjusted_VCF[growing.vcf.table$y_adjusted_VCF < 0] <- new.low.value
}

# Remove NA rows
growing.table.new <- growing.table[!is.na(growing.table$QAQC_MESMA), ]
growing.vcf.table.new <- growing.vcf.table[!is.na(growing.vcf.table$QAQC_VCF), ]

nrow(growing.table.new)
nrow(growing.vcf.table.new)


# Plot for the ForestSAT poster
# Define x and y first
xvar <- growing.table.new$hand.tree*100
yvar <- growing.table.new$QAQC_MESMA*100

xvar2 <- growing.vcf.table.new$hand.tree * 100
yvar2 <- growing.vcf.table.new$QAQC_VCF * 100

# Get the means
mean.mesma <- mean(yvar)
mean.vcf <- mean(yvar2)
mean.vcf.hand <- mean(xvar2)

# Get the RMSE
rmse <- sqrt(mean((yvar - xvar)^2 , na.rm = TRUE))
rmse_VCF <- sqrt(mean((yvar2 - xvar2)^2 , na.rm = TRUE))

# Get the VEcv
VE <- VEcv(observed = xvar, predicted = yvar)
VE_VCF <- VEcv(observed = xvar2, predicted = yvar2)


# Add to table
VE.RMSE.tbl[k + 11, ] <- c(as.character(vcf.files[k]), round(rmse, 0), round(VE, 0))

{
  png(paste0("E:/Dropbox/Permanent/Grad School/Meetings/ForestSAT/2018/Poster/Images/AllParks_MESMA_VCF_v_GT_ForestSAT.png"), 
      width = 1200, height = 1200)
  
  
  plot(x = yvar,
       y = xvar,
       main = paste0(vcf.files[k], " \nRMSE = ", round(rmse, 0), "% VEcv = ", round(VE, 0), "%"),
       ylab = "Reference % Woody Cover", xlab = "VCF % Tree Cover",
       xlim = c(0, 100), ylim = c(0, 100), cex = 2, pch = 19)
  
  points(x = yvar2,
       y = xvar2,
       cex = 2, col = "#FF5722", pch = 19)
  
  
  # Add 1:1 line
  abline(0,1, lwd = 4) # Add the 1:1 line
  
  # Add a regression line
  modl <- lm(xvar ~ yvar) # calc regression
  abline(modl, lty = 2, lwd = 5)
  
  
  modl <- lm(xvar2 ~ yvar2) # calc regression
  abline(modl, lty = 2, lwd = 6, col = "#FF5722")
  # # Add another line, shifted up to 1:1 line
  # mdl2 <- modl
  # mdl2$coefficients[[1]] <- 0
  # abline(mdl2, lty = 2, col = "red")
  
  dev.off()
}


# Plot the combined MESMA plot for normal and y adjusted
hdr.eval.df <- plot1to1(x = growing.table.new$hand.tree, y = growing.table.new$QAQC_MESMA, xlab = "Observed % Tree Cover", ylab = "MESMA % Tree Cover", 
                        main = "Percent Tree Cover: Observed (30m) vs. MESMA (30m)", ylim = c(0, 100), add_one2one_line = TRUE, add.reg.line = TRUE,
                        save.plot = save.plot, save.plot.name = "./images/All_Parks_MESMA_vs_GT.png", hdr.table = hdr.eval.df)

# hdr.eval.df <- plot1to1(x = growing.table.new$hand.tree, y = growing.table.new$y_adjusted_MESMA, xlab = "Observed % Tree Cover", ylab = "MESMA % Tree Cover", 
#                         main = "Percent Tree Cover: Observed (30m) vs. MESMA (30m)", ylim = c(-100, 300), add_one2one_line = TRUE, add.reg.line = TRUE,
#                         save.plot = save.plot, save.plot.name = "./images/All_Parks_MESMA_vs_GT_y_Adjusted.png", hdr.table = hdr.eval.df)
# Do the same for VCF
hdr.eval.df <- plot1to1(x = growing.vcf.table.new$hand.tree, y = growing.vcf.table.new$QAQC_VCF, xlab = "Observed % Tree Cover", ylab = "VCF % Tree Cover", 
                        main = "Percent Tree Cover: Observed (30m) vs. VCF (30m)", ylim = c(0, 100), add_one2one_line = TRUE, add.reg.line = TRUE,
                        save.plot = save.plot, save.plot.name = "./images/All_Parks_VCF_vs_GT.png", hdr.table = hdr.eval.df)

# hdr.eval.df <- plot1to1(x = growing.vcf.table.new$hand.tree, y = growing.vcf.table.new$y_adjusted_VCF, xlab = "Observed % Tree Cover", ylab = "VCF % Tree Cover", 
#                         main = "Percent Tree Cover: Observed (30m) vs. VCF (30m)", ylim = c(-100, 300), add_one2one_line = TRUE, add.reg.line = TRUE,
#                         save.plot = save.plot, save.plot.name = "./images/All_Parks_VCF_vs_GT_y_Adjusted.png", hdr.table = hdr.eval.df)




# Can use the plots above or write up new plots here.
# The last plot in the loops looked really good - slope of 1.04
# Maybe just use that, but give the stats from all the parks?
# Otherwise, the plot looks pretty sloppy. 
# Could also adjust all the data to the y-int for each park and see if that total plot looks better










################  Accuracy Assessment  #########################################################################
# Write a quick bin function
bin10 <- function(x){
  for(t in 1:length(x)){
    if(is.na(x[t])){
      val <- 0
    }else{
      if(x[t] < 0){
        val <- 0
      }else if(x[t] <= .1){
        val <- 1
      }else if(x[t] <= .2){
        val <- 2
      }else if(x[t] <= .3){
        val <- 3
      }else if(x[t] <= .4){
        val <- 4
      }else if(x[t] <= .5){
        val <- 5
      }else if(x[t] <= .6){
        val <- 6
      }else if(x[t] <= .7){
        val <- 7
      }else if(x[t] <= .8){
        val <- 8
      }else if(x[t] <= .9){
        val <- 9
      }else if(x[t] <= 1){
        val <- 10
      }else if(x[t] > 1){
        val <- 11
      }
    }
    
    if(t == 1){
      val.list <- c(val)
    }else{
      val.list <- c(val.list, val)
    }
  }
  return(val.list)
}


bin25 <- function(x){
  for(t in 1:length(x)){
    if(is.na(x[t])){
      val <- 0
    }else{
      if(x[t] < 0){
        val <- 0
      }else if(x[t] <= .25){
        val <- 1
      }else if(x[t] <= .50){
        val <- 2
      }else if(x[t] <= .75){
        val <- 3
      }else if(x[t] <= 1.0){
        val <- 4
      }else if(x[t] > 1){
        val <- 5
      }
    }
    
    if(t == 1){
      val.list <- c(val)
    }else{
      val.list <- c(val.list, val)
    }
  }
  return(val.list)
}



# Variance explained of original data
VEcv(observed = growing.table$hand.tree, predicted = growing.table$MESMA)
VEcv(observed = growing.vcf.table$hand.tree, predicted = growing.vcf.table$VCF)

VEcv(observed = growing.table.new$hand.tree, predicted = growing.table.new$QAQC_MESMA)
VEcv(observed = growing.vcf.table.new$hand.tree, predicted = growing.vcf.table.new$QAQC_VCF)




# Bin the hand.tree, VCF and MESMA values in 0-10, 11-20% etc
# Do for the raw values only. No y-intercept adjusting
growing.table$hand.binned <- bin10(growing.table$hand.tree)
growing.table$MESMA.binned <- bin10(growing.table$MESMA)
growing.vcf.table$hand.binned <- bin10(growing.vcf.table$hand.tree)
growing.vcf.table$VCF.binned <- bin10(growing.vcf.table$VCF)

# Bin the hand.tree, VCF and MESMA values in 0-10, 11-20% etc
# Do for the raw values only. No y-intercept adjusting
growing.table$hand.binned <- bin25(growing.table$hand.tree)
growing.table$MESMA.binned <- bin25(growing.table$MESMA)
growing.vcf.table$hand.binned <- bin25(growing.vcf.table$hand.tree)
growing.vcf.table$VCF.binned <- bin25(growing.vcf.table$VCF)

percent.correct.MESMA <- sum(growing.table$hand.binned == growing.table$MESMA.binned) / nrow(growing.table)
percent.correct.VCF <- sum(growing.vcf.table$hand.binned == growing.vcf.table$VCF.binned) / nrow(growing.vcf.table)



# The data is skewed, so need to use Spearman's correlation, not Pearson's
hist(growing.table$hand.tree)
hist(growing.table$MESMA)
hist(growing.vcf.table$hand.tree)
hist(growing.vcf.table$VCF)

# Check the correlation 
mesma.cor <- cor(growing.table$hand.tree, growing.table$MESMA, method = "spearman", use = "complete.obs")
vcf.cor <- cor(growing.vcf.table$hand.tree, growing.vcf.table$VCF, method = "spearman", use = "complete.obs")

# Check RMSE
rmse <- function(x, y){
  error <- x - y
  return(sqrt(mean(error^2)))
}

mae <- function(x, y){
  error <- x - y
  return(mean(abs(error)))
}

# RMSE - VCF does better when there is no data cleaning 
rmse(growing.table$hand.tree, growing.table$MESMA)
rmse(growing.vcf.table$hand.tree, growing.vcf.table$VCF)
mae(growing.table$hand.tree, growing.table$MESMA)
mae(growing.vcf.table$hand.tree, growing.vcf.table$VCF)

# RMSE filtered data
rmse(growing.table.new$hand.tree, growing.table.new$QAQC_MESMA)
rmse(growing.vcf.table.new$hand.tree, growing.vcf.table.new$QAQC_VCF)
mae(growing.table.new$hand.tree, growing.table.new$QAQC_MESMA)
mae(growing.vcf.table.new$hand.tree, growing.vcf.table.new$QAQC_VCF)

# Look at the linear models
mesma.lm <- lm(growing.table$MESMA ~ growing.table$hand.tree)
vcf.lm <- lm(growing.vcf.table$VCF ~ growing.vcf.table$hand.tree)
summary(mesma.lm)
summary(vcf.lm)
