##################################################################
### This script can be used to run the ClustMod model, a model to 
### predict snow distribution patterns based from observations. 
### This script is devided into:
###
### 1) Model Preparation: Preparation of training and test sites.
###    This includes the determination of model features and clusters
###    using the ClustSnow workflow.
### 2) Merging training data.
### 3) Training of the ClustMod model
### 4) Evaluation of the ClustMod model
### 5) Predicting clusters to testsites 
### 6) Application of predicted clusters.
###
##################################################################
### Version 4 
### Last Update 2025-01-05 by Joschka Geissler
##################################################################

rm(list = ls())

# Read required Packages
library(rstudioapi)
library(raster)
library(gstat)
library(tidyverse)
library(sf)
library(whitebox)
library(caret)
library(assertthat)
library(ForestTools)
library(tmap) 
library(R.utils)
library(reshape2)
library(windninjr) # devtools::install_github("SCAR-sandpit/windninjr")
#-------------------------------------------------------
# Load Functions and define working directory
# NOTE: This Rscript must be placed in the same directory as the ClustMod_functions and ClustSnow_functions scripts.
#-------------------------------------------------------

# Load functions
working_directory <- dirname(rstudioapi::getSourceEditorContext()$path) # Get Location of Rscript
source(paste0(working_directory,'\\ClustMod_Functions_v4.R'))
source(paste0(working_directory,'\\ClustSnow_Functions_v4.R'))
source(paste0(working_directory,'\\windninjr.R'))

#-------------------------------------------------------
# Prepare Training Data and Features
# NOTE: Refer to folders of Input and Output. Note the folder structure and minimum requirements explicitly described in the readme file.
#-------------------------------------------------------
setwd(working_directory)

# Read Model Configuration file
path_config <- list.files(working_directory,pattern='config_ClustMod',full.names = T,recursive = T)
readConfigFile(path_config)

path_in  <- paste0(working_directory,'\\',path_in) # Path to Input Data Folders
path_out <- paste0(working_directory,'\\',path_out) # Path to write Results to

# Create Output Folder if Non-Existent
if (!dir.exists(path_out)) {
  dir.create(path_out, recursive = TRUE)
  cat("Output folder created at:", path_out, "\n")
} else {
  cat("Output folder already exists at:", path_out, "\n")
}

# Logfile handling
logfile <- paste0(path_out,"\\ClustMod_logfile.txt")
if (file.exists(logfile)) {
  file.remove(logfile)
}
  write(paste0(Sys.time(), " - Reading configuration file: ", path_config), file = logfile, append = TRUE)
  write("\n", file = logfile, append = TRUE)
  config_lines <- apply(read.delim(path_config), 1, paste, collapse = "\t")  # Combine columns with tabs
  write(config_lines, file = logfile,append = TRUE)
  write("\n", file = logfile, append = TRUE)

# -------------------------------------------------------
# 1) Model Preparation
# -------------------------------------------------------
# Loop over all domains and run 'prepareData' Function
dirs <- normalizePath(list.dirs(path_in,recursive = F,full.names = T))

for (dir in dirs) {

setwd(dir) 
  
  # Read ConfigFiles
  path_config <- list.files(dir,pattern = 'config_data_preparation',full.names = T)
  readConfigFile(path_config)

  #Logfile Handling
  write(paste0(Sys.time()," - Start Data Preparation of", dir), file = logfile, append = TRUE)
  write("\n", file = logfile, append = TRUE)
  write(paste0(Sys.time(), " - Reading configuration file: ", path_config), file = logfile, append = TRUE)
  write("\n", file = logfile, append = TRUE)
  config_lines <- apply(read.delim(path_config), 1, paste, collapse = "\t")  # Combine columns with tabs
  write(config_lines, file = logfile,append = TRUE)
  write("\n", file = logfile, append = TRUE)
  

# Read Data with informations provided by configuration file
dtm_1 <- raster(dtm_path,crs=crs_dtm)
chm_1 <- raster(chm_path,crs=crs_chm)
SRTM <- raster(dtm_largescale_path,crs=crs_dtm_largescale) #Instead the getSRTM(country,dir_out) function could be used.

# Handling different input resolutions.
res_chm <- round(res(chm_1)[1],0)
res_dtm <- round(res(dtm_1)[1],0)

if (res_chm!=1 | res_dtm!=1) {
  write(paste0(Sys.time(),' - Resampling of Input Data (CHM or DTM) to Resolution of 1m!'), file = logfile, append = TRUE)
  warning('Resampling of Input Data to Resolution of 1m!')
  
  chm_1 <- raster::aggregate(chm_1,1/res_chm)
  res(chm_1)<-c(1,1)
  is.na(chm_1[])<-0
  dtm_1 <- raster::resample(dtm_1,chm_1)
}

# Align chm and dtm and fill data gaps
dtm_1 <- raster::resample(dtm_1,chm_1)
chm_1 <-  focal(chm_1, w = matrix(1, 3, 3), fun = fill_na)
chm_1 <-  focal(chm_1, w = matrix(1, 3, 3), fun = fill_na)
dtm_1 <-  focal(dtm_1, w = matrix(1, 5, 5), fun = fill_na)

# Compute Features when FINALIZE is set to FALSE.
if (!finalize) {

write(paste0(Sys.time(),' - Start prepareData Function for ', dir), file = logfile, append = TRUE)

  prepareData(chm_1 = chm_1,
             dtm_1 = dtm_1,
             dtm_largescale = SRTM,
             path_out=dir,
             DCE_step_nr = DCE_step_nr,
             wind_direction = wind_direction,
             wind_velocity = wind_velocity,
             path_to_wn_exe = path_to_wn_exe)

}

if (finalize) {
  
  write(paste0(Sys.time(),'- ',dir,' is finalized. No Data is being computed.'), file = logfile, append = TRUE)
  warning(paste(dir,'is finalized. No Data is being computed.'))
}

# Compute Clusters
if (isTrain) {
  warning(paste0(Sys.time(),': Start Computing Clusters'))
  write(paste0(Sys.time(),'- ',dir,' is finalized. No Data is being computed.'), file = logfile, append = TRUE)
  
  # Read HS maps
  files <- list.files(snow_depth_path,full.names = T)
  files <- files[stringr::str_detect(files,pattern='tif')]
  snow_stack <- stack(files)
  
  # Calculate Cluster
  cali_data <- data.frame(sample_length,n_class,kmeans_maxiter,kmeans_nstart,mtry,n_trees)
  
  set.seed(123)
  cluster <- getCluster(calibration_data = cali_data, data_for_prediction=snow_stack)
  cluster <- orderCluster(cluster,snow_stack)

  if (crs_snowdata != crs_chm) {
    crs(cluster)<-crs_snowdata
    cluster <- raster::projectRaster(cluster,chm_1,method='ngb')
    
  }
  
  # Save Cluster
  cluster <- raster::resample(cluster,chm_1,method='ngb')
  abs_cluster <- cluster[[nlayers(cluster)]]
  prob_cluster <- cluster[[-nlayers(cluster)]]
  
  writeRaster(abs_cluster,paste0(dir,'\\cluster.asc'),overwrite=T)
  writeRaster(prob_cluster,paste0(dir,'\\prob_cluster.tif'),overwrite=T)
  
  }

}

# -------------------------------------------------------
# 2) Putting together training data
# -------------------------------------------------------
# Read Train Data and merge to dataframe
setwd(path_in)

study_sites_path <- list.dirs(path_in,recursive=F)[list.dirs(path_in,recursive=F,full.names = F)%in%train_sites]

data <- FeaturesToDataFrame(features,
                            'prob_cluster.tif',
                            study_sites_path)
data <- na.omit(data)
data$id <- 1:nrow(data)

# -------------------------------------------------------
# 3) Train Model
# -------------------------------------------------------
write(paste0(Sys.time(),'- ','Start Training of ClustMod'), file = logfile, append = TRUE)

# Split Data in Training and Test Data
set.seed(789)
id_test <- sample(data$id,(1-train_test_split)*nrow(data))

train_data <- data[!(data$id%in%id_test),]
test_data <- data[(data$id%in%id_test),]

# Delete ID Column
train_data <- train_data[,-ncol(train_data)]
test_data <- test_data[,-ncol(test_data)]

x = train_data[,!colnames(train_data)%in%paste0('X',1:4)]
y1 = as.numeric(train_data[,colnames(train_data)%in%'X1'])
y2 = as.numeric(train_data[,colnames(train_data)%in%'X2'])
y3 = as.numeric(train_data[,colnames(train_data)%in%'X3'])
y4 = as.numeric(train_data[,colnames(train_data)%in%'X4'])


model_x1 <- randomForest::randomForest(x,y1)
model_x2 <- randomForest::randomForest(x,y2)
model_x3 <- randomForest::randomForest(x,y3)
model_x4 <- randomForest::randomForest(x,y4)

write(paste0(Sys.time(),'- ','Training of ClustMod completed'), file = logfile, append = TRUE)
# -------------------------------------------------------
# Overwrite/Use Pre-Saved Data
# -------------------------------------------------------
save.image(paste0(path_out,'\\ClustMod_POST_TRAINING.RData'))
#load(paste0(path_out,'\\ClustMod_POST_TRAINING.RData'))

# -------------------------------------------------------
# 4) Model Evaluation
# -------------------------------------------------------
x_test  <- test_data[,!colnames(test_data)%in%paste0('X',1:4)]

y1_test <- as.numeric(test_data[,colnames(test_data)%in%'X1'])
y2_test <- as.numeric(test_data[,colnames(test_data)%in%'X2'])
y3_test <- as.numeric(test_data[,colnames(test_data)%in%'X3'])
y4_test <- as.numeric(test_data[,colnames(test_data)%in%'X4'])

y1_test_pred <- predict(model_x1,x_test)
y2_test_pred <- predict(model_x2,x_test)
y3_test_pred <- predict(model_x3,x_test)
y4_test_pred <- predict(model_x4,x_test)

predicted_probs <- data.frame(
  cl1 = y1_test_pred,
  cl2 = y2_test_pred,
  cl3 = y3_test_pred,
  cl4 = y4_test_pred
)

obs_probs <- data.frame(
  cl1 = y1_test,
  cl2 = y2_test,
  cl3 = y3_test,
  cl4 = y4_test
)

# Determine the winner class by finding the class with the highest probability
winner_pred <- apply(predicted_probs, 1, which.max)
winner_obs <- apply(obs_probs, 1, which.max)

# Create confusion matrix
confusion_matrix <- table(Actual = winner_obs, Predicted = winner_pred)

# Calculate Overall Accuracy
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
write.csv2(confusion_matrix,paste0(path_out,'\\confusionmatrix.csv'))

# Logfile Handling
write('\n', file = logfile, append = TRUE)
write(paste0(Sys.time(),' - ',"Accuracy: ", round(accuracy,2)), file = logfile, append = TRUE)
write('\n', file = logfile, append = TRUE)
write('confusion_matrix:', file = logfile, append = TRUE)
write(confusion_matrix, file = logfile, append = TRUE)


# -------------------------------------------------------
# 5) Predictions
# -------------------------------------------------------

for (site in test_sites) {

  write('\n', file = logfile, append = TRUE)
  write(paste0(Sys.time(),"- Start Cluster Prediction for site ",site), file = logfile, append = TRUE)
  
  dir_in <- paste0(path_in,'\\',site)

  testsite_data <- stack(paste0(dir_in,'\\',features,'.asc'))
  testsite_data <- raster::aggregate(testsite_data,1)
  testsite_data$TWI_5[is.na(testsite_data$TWI_5)] <- 10
  
  testsite_data <- as.data.frame(testsite_data,xy=T)
  testsite_data <- na.omit(testsite_data)
  
  pred.x1 <- predict(model_x1,testsite_data[,!(colnames(testsite_data)%in%c('x','y'))])
  pred.x2 <- predict(model_x2,testsite_data[,!(colnames(testsite_data)%in%c('x','y'))])
  pred.x3 <- predict(model_x3,testsite_data[,!(colnames(testsite_data)%in%c('x','y'))])
  pred.x4 <- predict(model_x4,testsite_data[,!(colnames(testsite_data)%in%c('x','y'))])
  
  result <- data.frame(x=testsite_data$x,
                       y=testsite_data$y,
                       X1=pred.x1,
                       X2=pred.x2,
                       X3=pred.x3,
                       X4=pred.x4)
  
  result$cluster <- max.col(result[, c("X1", "X2", "X3", "X4")], ties.method = "first")
  
  r_result <- rasterFromXYZ(result)
  writeRaster(r_result,paste0(path_out,'/',site,'_Predictions.tif'),overwrite=T)
  
  write(paste0(Sys.time(),"- End Prediction for site ",site), file = logfile, append = TRUE)
  
}

# -------------------------------------------------------
# Overwrite/Use Pre-Saved Data
# -------------------------------------------------------
save.image(paste0(path_out,'\\ClustMod_POST_Prediction.RData'))
#load(paste0(path_out,'\\ClustMod_POST_Prediction.RData'))

# Logfile handling
closeAllConnections()

# -------------------------------------------------------
# 6) Application
# -------------------------------------------------------
#load(paste0(path_out,'\\ClustMod_POST_TRAINING.RData'))

 # percentage_samples       = 0.01
 # 
 # snow_depth <- raster('PATH TO YOUR REFERENCE SNOW DEPTH MAP YOU WANT TO REPRODUCE.')
 # cluster <- stack('PATH TO YOUR PREDICTED CLUSTERS.')
 # aoi <- shapefile('PATH TO SHAPEFILE OF YOUR AOI,')
 # chm <- raster('PATH TO YOUR OF THE STUDY SITE THE CLUSTERS WERE PREDICTED FOR.')
 # 
 # cluster <- crop(cluster,aoi)
 # snow_depth <- raster::resample(snow_depth,cluster)
 # chm <- raster::resample(chm,cluster)
 # 
 # sensorlocations <- sample_sensorlocations(prob_cluster=cluster,
 #                                           n_points=round(ncell(cluster)*percentage_samples/4))
 # HS_cluster_Map <- interpolateSensorMeasurements(prob_cluster=cluster,
 #                                                 snow_depth_raster=snow_depth,
 #                                                 sensorlocations)
 # HS_predicted <- HS_cluster_Map$X1*cluster[[1]]+HS_cluster_Map$X2*cluster[[2]]+HS_cluster_Map$X3*cluster[[3]]+HS_cluster_Map$X4*cluster[[4]]
 # err <- err_assessment(HS_predicted,snow_depth)
 # err$percentage_samples <- percentage_samples
 # print(err)

