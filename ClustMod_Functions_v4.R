#####################################################################
### Skript containing all functions needed for the ClusMod Model. ###
#####################################################################

# ------------------------------------------
# 1) Functions to determine model features
# ------------------------------------------

calc_DCE <- function(chm_input,step_nr,wind_direction=NA){
  
  require(raster)
  require(stringr)
  
  message(paste(Sys.time(),': Compute DCE with step_nr ', step_nr))

  # Binarize CHM
  chm_input[chm_input <= 2] = 0; 
  chm_input[chm_input > 2] = 1;   
  
  # CALCULATIONS 
  # 1. Convert to 1m resolution 
  
  chm_data = chm_input; 
  agg.factor = 1/res(chm_data)[1]
  
  if(agg.factor!=1){
    chm_data = aggregate(chm_data,agg.factor)
  }
  
  chm_data[is.na(chm_data)]<-0

  chm_data.data = as.data.frame(chm_data,xy=T)
  colnames(chm_data.data)[3]<-'chm_data'
  chm_data.data$chm_data[chm_data.data$chm_data >= 0.5] = 1
  chm_data.data$chm_data[chm_data.data$chm_data < 0.5] = 0    
  # 2. Calculate non-directional DCE 
  
  # initialize DCE of open pixels: matrix with -1 = canopy, 0 = open 
  chm_data.data$opnclasses = -chm_data.data$chm_data;   
  # initialize DCE of canopy pixels: matrix with -1 = open, 0 = canopy 
  chm_data.data$canclasses = chm_data.data$chm_data-1;  
  # input to first step: binary chm grid with 1 = canopy pixels, 0 = open
  # pixels
  chm_data.data$ingrid_opn = chm_data.data$chm_data;
  chm_data.data$ingrid_can = -(chm_data.data$chm_data-1);
  
  
  kernel_disk       <- matrix(c(0,1,0,1,1,1,0,1,0),nrow=3,ncol=3,byrow = T)
  kernel_disk_south <- create_directional_kernel(0)
  kernel_disk_north <- create_directional_kernel(180)
  
  x <- chm_data.data$x
  y <- chm_data.data$y
  
  # iterative edge detection for open and canopy pixel DCE separately 
  for(ssx in 1:step_nr) { 

    # generate smoothed chm
    chm_data.data$chmsm_opn = as.data.frame(focal(rasterFromXYZ(data.frame(x,y,chm_data.data$ingrid_opn)),kernel_disk,pad=T))       
    chm_data.data$chmsm_can = as.data.frame(focal(rasterFromXYZ(data.frame(x,y,chm_data.data$ingrid_can)),kernel_disk,pad=T))      
    
    # detect edges and attribute value = step nr. 
    intersection_can <- intersect(which(chm_data.data$canclasses == 0), which(chm_data.data$chmsm_can > 0))
    intersection_opn <- intersect(which(chm_data.data$opnclasses == 0), which(chm_data.data$chmsm_opn > 0))
    
    # Assign the value 0.001 to the elements in opnclasses based on the intersection
    chm_data.data$opnclasses[intersection_opn] <- ssx
    chm_data.data$canclasses[intersection_can] <- ssx
    
    # smoothed and re-binarized grids create input for next iteration;
    # use max filter for open and min filter for canopy pixels
    chm_data.data$ingrid_opn <- chm_data.data$chmsm_opn
    chm_data.data$ingrid_opn[chm_data.data$ingrid_opn > 0] <-  1
    
    chm_data.data$ingrid_can <- chm_data.data$chmsm_can
    chm_data.data$ingrid_can[chm_data.data$ingrid_can > 0] <-  1
    
  } 
  
  # set all canopy pixels (open pixels) to 0 in opn (can) DCE matrices 
  # (to allow merging of DCE of canopy and open pixels later on)
  chm_data.data$opnclasses[chm_data.data$opnclasses < 0] = 0; 
  chm_data.data$canclasses[chm_data.data$canclasses < 0] = 0;
  
  # merge DCE of canopy and open pixels 
  chm_data.data$dceall_grid = chm_data.data$canclasses-chm_data.data$opnclasses;
  # set values that have not been defined to NaN
  chm_data.data$dceall_grid[chm_data.data$dceall_grid == 0] = NA;
  
  # initialize output struct
  dce_output = raster::resample(rasterFromXYZ(data.frame(x,y,-chm_data.data$dceall_grid)),chm_data)

  # 3. Calculate directional DCE 
  # DCE-north (for north-exposed edges) 
  # initialize
  chm_data.data$opnclasses = -chm_data.data$chm_data;
  chm_data.data$canclasses = chm_data.data$chm_data-1;     
  chm_data.data$ingrid_opn = chm_data.data$chm_data;
  chm_data.data$ingrid_can = -(chm_data.data$chm_data-1);
  
  for (ssx in (1:step_nr)){ 
    # smoothing with asymmetric kernel
    chm_data.data$chmsm_opn = as.data.frame(focal(rasterFromXYZ(data.frame(x,y,chm_data.data$ingrid_opn)),kernel_disk_north,pad=T))       
    chm_data.data$chmsm_can = as.data.frame(focal(rasterFromXYZ(data.frame(x,y,chm_data.data$ingrid_can)),kernel_disk_south,pad=T))      
    
    
    # detect edges and attribute value = step nr. 
    intersection_can <- intersect(which(chm_data.data$canclasses == 0), which(chm_data.data$chmsm_can > 0))
    intersection_opn <- intersect(which(chm_data.data$opnclasses == 0), which(chm_data.data$chmsm_opn > 0))
    
    # Assign ssx to the elements in opnclasses based on the intersection
    chm_data.data$opnclasses[intersection_opn] <- ssx
    chm_data.data$canclasses[intersection_can] <- ssx
    
    # compute input to next smoothing iteration
    chm_data.data$ingrid_opn <- chm_data.data$chmsm_opn
    chm_data.data$ingrid_opn[chm_data.data$ingrid_opn > 0] <-  1
    
    chm_data.data$ingrid_can <- chm_data.data$chmsm_can
    chm_data.data$ingrid_can[chm_data.data$ingrid_can > 0] <-  1   
    
  }
  
  # merge DCE-north of open and canopy pixels
  chm_data.data$ndceall_grid = NA 
  chm_data.data$ndceall_grid[chm_data.data$opnclasses > 0 & !is.na(chm_data.data$opnclasses)] = -chm_data.data$opnclasses[chm_data.data$opnclasses > 0 & !is.na(chm_data.data$opnclasses)] 
  chm_data.data$ndceall_grid[chm_data.data$canclasses > 0 & !is.na(chm_data.data$canclasses)] =  chm_data.data$canclasses[chm_data.data$canclasses > 0 & !is.na(chm_data.data$canclasses)] 
  
  chm_data.data$ndceall_grid[chm_data.data$canclasses == 0] = NA
  chm_data.data$ndceall_grid[chm_data.data$opnclasses == 0] = NA
  
  ndce_output = raster::resample(rasterFromXYZ(data.frame(x,y,-chm_data.data$ndceall_grid)),chm_data)
  
  # DCE-south (for south-exposed edges) 
  chm_data.data$opnclasses = -chm_data.data$chm_data;
  chm_data.data$canclasses = chm_data.data$chm_data-1; 
  chm_data.data$ingrid_opn = chm_data.data$chm_data;
  chm_data.data$ingrid_can = -(chm_data.data$chm_data-1);
  
  for (ssx in (1:step_nr)){ 
    # smoothing with asymmetric kernel
    chm_data.data$chmsm_opn = as.data.frame(focal(rasterFromXYZ(data.frame(x,y,chm_data.data$ingrid_opn)),kernel_disk_south,pad=T))       
    chm_data.data$chmsm_can = as.data.frame(focal(rasterFromXYZ(data.frame(x,y,chm_data.data$ingrid_can)),kernel_disk_north,pad=T))      
    
    # detect edges and attribute value = step nr. 
    intersection_can <- intersect(which(chm_data.data$canclasses == 0), which(chm_data.data$chmsm_can > 0))
    intersection_opn <- intersect(which(chm_data.data$opnclasses == 0), which(chm_data.data$chmsm_opn > 0))
    
    # Assign the value 0.001 to the elements in opnclasses based on the intersection
    chm_data.data$opnclasses[intersection_opn] <- ssx
    chm_data.data$canclasses[intersection_can] <- ssx
    
    # compute input to next smoothing iteration
    chm_data.data$ingrid_opn <- chm_data.data$chmsm_opn
    chm_data.data$ingrid_opn[chm_data.data$ingrid_opn > 0] <-  1
    
    chm_data.data$ingrid_can <- chm_data.data$chmsm_can
    chm_data.data$ingrid_can[chm_data.data$ingrid_can > 0] <-  1   
    
  }
  
  chm_data.data$sdceall_grid = NA 
  chm_data.data$sdceall_grid[chm_data.data$opnclasses > 0 & !is.na(chm_data.data$opnclasses)] = -chm_data.data$opnclasses[chm_data.data$opnclasses > 0 & !is.na(chm_data.data$opnclasses)] 
  chm_data.data$sdceall_grid[chm_data.data$canclasses > 0 & !is.na(chm_data.data$canclasses)] =  chm_data.data$canclasses[chm_data.data$canclasses > 0 & !is.na(chm_data.data$canclasses)] 
  
  chm_data.data$sdceall_grid[chm_data.data$canclasses == 0] = NA
  chm_data.data$sdceall_grid[chm_data.data$opnclasses == 0] = NA
  
  sdce_output = raster::resample(rasterFromXYZ(data.frame(x,y,-chm_data.data$sdceall_grid)),chm_data)
  
  dce_stack <- stack(dce_output,ndce_output,sdce_output)
  names(dce_stack)<-c('DCE_1','NDCE_1','SDCE_1')
  
  if (!is.na(wind_direction)) {
    
    for (i in 1:2) {
      direction <- c(wind_direction,wind_direction-180)
      
      if (wind_direction<180) {
        direction[2]<-direction[2]+360
      }
      
      direction <- direction[i]
      
      
      # Wind_DCE (for wind-facing canopy edges)
      # Initialize
      chm_data.data$opnclasses = -chm_data.data$chm_data;
      chm_data.data$canclasses = chm_data.data$chm_data-1; 
      chm_data.data$ingrid_opn = chm_data.data$chm_data;
      chm_data.data$ingrid_can = -(chm_data.data$chm_data-1);
      
      # Create Kernel
      kernel_disk_windfacing <- create_directional_kernel(direction)
      kernel_disk_windshaded <- kernel_disk_windfacing[nrow(kernel_disk_windfacing):1,ncol(kernel_disk_windfacing):1]
      
      # And GO!
      for (ssx in (1:step_nr)){ 
        #print(ssx)
        # smoothing with asymmetric kernel
        chm_data.data$chmsm_opn = as.data.frame(focal(rasterFromXYZ(data.frame(x,y,chm_data.data$ingrid_opn)),kernel_disk_windfacing,pad=T))       
        chm_data.data$chmsm_can = as.data.frame(focal(rasterFromXYZ(data.frame(x,y,chm_data.data$ingrid_can)),kernel_disk_windshaded,pad=T))      
        
        # detect edges and attribute value = step nr. 
        intersection_can <- intersect(which(chm_data.data$canclasses == 0), which(chm_data.data$chmsm_can > 0))
        intersection_opn <- intersect(which(chm_data.data$opnclasses == 0), which(chm_data.data$chmsm_opn > 0))
        
        # Assign the value 0.001 to the elements in opnclasses based on the intersection
        chm_data.data$opnclasses[intersection_opn] <- ssx
        chm_data.data$canclasses[intersection_can] <- ssx
        
        # compute input to next smoothing iteration
        chm_data.data$ingrid_opn <- chm_data.data$chmsm_opn
        chm_data.data$ingrid_opn[chm_data.data$ingrid_opn > 0] <-  1
        
        chm_data.data$ingrid_can <- chm_data.data$chmsm_can
        chm_data.data$ingrid_can[chm_data.data$ingrid_can > 0] <-  1   
        
      }
      
      chm_data.data$windceall_grid = NA 
      chm_data.data$windceall_grid[chm_data.data$opnclasses > 0 & !is.na(chm_data.data$opnclasses)] = -chm_data.data$opnclasses[chm_data.data$opnclasses > 0 & !is.na(chm_data.data$opnclasses)] 
      chm_data.data$windceall_grid[chm_data.data$canclasses > 0 & !is.na(chm_data.data$canclasses)] =  chm_data.data$canclasses[chm_data.data$canclasses > 0 & !is.na(chm_data.data$canclasses)] 
      
      chm_data.data$windceall_grid[chm_data.data$canclasses == 0] = NA
      chm_data.data$windceall_grid[chm_data.data$opnclasses == 0] = NA
      
      WIND_dceall_output = raster::resample(rasterFromXYZ(data.frame(x,y,-chm_data.data$windceall_grid)),chm_data)
      
      if(nlayers(dce_stack)>0){
        names_stack <- names(dce_stack)
        dce_stack <- stack(dce_stack,WIND_dceall_output)
        names(dce_stack)<-c(names_stack,paste0('WIND_DCE_',i))
      }
      
      if(nlayers(dce_stack)==0){
        
        dce_stack <- stack(dce_stack,WIND_dceall_output)
        names(dce_stack)<-c(paste0('WIND_DCE_',i))
        
      }
      
      
    }
    
  }
  
  dce_final_data <- (as.data.frame(dce_stack,xy=T))
  dce_final_data$DCE_1[is.na(dce_final_data$DCE_1)]<-max(dce_final_data$DCE_1,na.rm=T)
  dce_final_data$NDCE_1[is.na(dce_final_data$NDCE_1)]<-max(dce_final_data$NDCE_1,na.rm=T)
  dce_final_data$SDCE_1[is.na(dce_final_data$SDCE_1)]<-max(dce_final_data$SDCE_1,na.rm=T)
  dce_final_data$WIND_DCE_1[is.na(dce_final_data$WIND_DCE_1)]<-max(dce_final_data$WIND_DCE_1,na.rm=T)
  dce_final_data$WIND_DCE_2[is.na(dce_final_data$WIND_DCE_2)]<-max(dce_final_data$WIND_DCE_2,na.rm=T)
  
  dce_stack_final <- rasterFromXYZ(dce_final_data)
  plot(dce_stack_final)
  return(dce_stack_final)
  
}

calcTWI <- function(dtm_1,tmp_folder='D:\\Test'){
  require(tidyverse)
  require(raster)
  require(sf)
  require(whitebox)
  require(tmap) 
  
  # Set WD to tempfolder
  old_wd <- getwd()
  setwd(tmp_folder)
  
  writeRaster(dtm_1,'dtm.tif',overwrite=T)
  
  wbt_breach_depressions_least_cost(
    dem = ".\\dtm.tif",
    output = ".\\dtm_breached.tif",
    dist = 5,
    fill = TRUE)
  
  wbt_fill_depressions_wang_and_liu(
    dem = ".\\dtm_breached.tif",
    output = ".\\dtm_filled_breached.tif")
  
  wbt_d_inf_flow_accumulation(input = ".\\dtm_filled_breached.tif",
                              output = ".\\DinfFAsca.tif",
                              out_type = "Specific Contributing Area")
  
  wbt_slope(dem = ".\\dtm_filled_breached.tif",
            output = ".\\demslope.tif",
            units = "degrees")
  
  wbt_wetness_index(sca = ".\\DinfFAsca.tif",
                    slope = ".\\demslope.tif",
                    output = ".\\TWI.tif")
  
  TWI <- raster("TWI.tif")
  
  file.remove(c('DinfFAsca.tif','demslope.tif','dtm_filled_breached.tif','dtm_breached.tif','dtm.tif'))
  
  setwd(old_wd)
  
  return(TWI)
}

applyWindNinja <- function(dtm_1, SRTM, wind_direction, wind_velocity, path_out, path_to_wn_exe){
  
  require(raster)
  require(windninjr)
  require(assertthat)
  
  w <- matrix(1, 3, 3)
  
  # Making sure that dtm_1 does not contain any NA values.
  repeat{
    dtm_1 <- focal(dtm_1, w, mean, na.rm=TRUE, NAonly=TRUE, pad=TRUE)
    
    # If no NA-Values are present in dtm, break loop. 
    if(is.na((table(as.data.frame(is.na(dtm_1))))['TRUE'])){
      break()
    }
    
    
  }
  
  closeAllConnections()
  writeRaster(dtm_1,paste0(path_out,'\\dem_wn.tif'),overwrite=T)
  
  ## set up a domain-average model run with 10m/s winds blowing towards the west
  config <- wn_config_domain_average(elevation = paste0(path_out,'\\dem_wn.tif'),
                                     input_speed = wind_velocity,
                                     input_direction = wind_direction,
                                     input_speed_units = "mps",
                                     output_speed_units = "mps",
                                     input_wind_height = 10,
                                     units_input_wind_height = "m",
                                     number_of_iterations = 300L,
                                     output_wind_height = 2,
                                     units_output_wind_height = "m")
  
  wn_find_exe(path_to_wn_exe)
  res <- wn_run(config)
  
  x <- wn_read(res$output_dir)
  names(x) <- c('WNDIR','WNVEL')
  
  x <- raster::resample(x,dtm_1)
  plot(x)
  
  WNDIR=x$WNDIR
  WNVEL=x$WNVEL
  
  writeRaster(WNDIR,paste0(path_out,'\\WNDIR_1.asc'),overwrite=T)
  writeRaster(WNVEL,paste0(path_out,'\\WNVEL_1.asc'),overwrite=T)
  
  wind <- stack(WNDIR,WNVEL)
  
  file.remove(paste0(path_out,'\\dem_wn.tif'))
  return(wind)
}

distanceToTree <- function(chm_1, path_out){
  
  # Load required libraries
  require(raster)
  require(sf)
  require('ForestTools')
  
  # Load CHM raster
  chm <- chm_1
  
  # Detect Trees using ForestTools
  # Function for defining dynamic window size
  lin <- function(x){x * 0.05 + 0.6}
  
  # Detect treetops
  ttops <- vwf(chm, winFun = lin, minHeight = 2)
  
  # Convert to SpatialPointsDataFrame and save as shapefile
  tree_spdf <- SpatialPointsDataFrame(ttops, data = data.frame(ID = 1:length(ttops)))
  tree_points <- st_as_sf(tree_spdf)
  
  distance_raster <- chm
  distance_raster[] <- NA
  
  # Identify cells corresponding to tree locations and set them as NA in CHM raster
  for (i in 1:length(ttops)) {
    cell_number <- cellFromXY(chm, ttops[i,]@coords)
    distance_raster[cell_number] <- 0
  }
  
  # Calculate distance to the closest tree
  dist_to_tree <- distance(distance_raster)
  
  dist_to_tree <- raster::resample(dist_to_tree,chm)
  # Save the distance raster
  writeRaster(dist_to_tree, filename = paste0(path_out,"\\DIST_1.asc"),overwrite=T)
  
}

# ------------------------------------------
# 2) Helper Functions

# ------------------------------------------

# Function to create a directional kernel for the DCE Algorithm
create_directional_kernel <- function(direction_degrees) {
  closestPossible <- seq(0,360,45)[abs(seq(0,360,45)-direction_degrees)==min(abs(seq(0,360,45)-direction_degrees))]
  
  # Initialize a 3x3 matrix filled with zeros
  kernel <- matrix(0, nrow = 3, ncol = 3)
  
  # Convert the direction to radians
  direction_radians <- (closestPossible-90) * pi / 180
  
  # Calculate the x and y coordinates for the direction
  x_coord <- round(cos(direction_radians))
  y_coord <- round(sin(direction_radians))
  
  # Set the center cell of the kernel to 1
  kernel[2, 2] <- 1
  
  # Set the cell in the specified direction
  if (any(x_coord == 0) && any(y_coord == -1)) {
    kernel[1, 2] <- 1
  } 
  if (any(x_coord == 1) && any(y_coord == -1)) {
    kernel[1, 3] <- 1
  } 
  if (any(x_coord == 1) && any(y_coord == 0)) {
    kernel[2, 3] <- 1
  } 
  if (any(x_coord == 1) && any(y_coord == 1)) {
    kernel[3, 3] <- 1
  } 
  if (any(x_coord == 0) && any(y_coord == 1)) {
    kernel[3, 2] <- 1
  } 
  if (any(x_coord == -1) && any(y_coord == 1)) {
    kernel[3, 1] <- 1
  } 
  if (any(x_coord == -1) && any(y_coord == 0)) {
    kernel[2, 1] <- 1
  } 
  if (any(x_coord == -1) && any(y_coord == -1)) {
    kernel[1, 1] <- 1
  }
  
  
  return(kernel)
}

# Function to calculate combined score.
combined_score <- function(RMSE,R){
  v_RMSE <- abs(min(RMSE)-RMSE)/(max(RMSE)-min(RMSE))
  v_R    <- (max(R)-R)/(max(R)-min(R))
  
  combined_score <- v_RMSE+v_R
  
  return(combined_score)
  
}

# Function to set na values to mean of surronding cells.
fill_na <- function(values) {
  if (is.na(values[5])) {
    return(mean(values, na.rm = TRUE))
  } else {
    return(values[5])
  }
}

# ------------------------------------------
# 3) Overarching Functions for the ClustMod Workflow
# ------------------------------------------
# Function to create all model features from chm and dtm
prepareData <- function(chm_1,dtm_1,dtm_largescale,path_out,DCE_step_nr,wind_direction,wind_velocity,path_to_wn_exe){
  
  require(raster)
  require(stringr)
  
  message(paste0(Sys.time(),': Start Data Preparation \n \n'))
  
  ### Derive all variables
  
  focal.5  <- focalWeight(chm_1, 5, type='circle') 
  focal.20  <- focalWeight(chm_1, 20, type='circle') 
  focal.50 <- focalWeight(chm_1, 50, type='circle')
  
  focal.5_nonwweight <- focal.5
  focal.5_nonwweight[focal.5_nonwweight>0] <- 1
  
  chm_binary <- chm_1 > 2
  
  # Canopy Closure 5 and 50
  CC_5 <- focal(chm_binary,focal.5,fun='sum',na.rm=F)
  CC_5 <- raster::resample(CC_5,chm_1)
  CC_50 <- focal(chm_binary,focal.50,fun='sum',na.rm=F)
  CC_50 <- raster::resample(CC_50,chm_1)
  
  writeRaster(CC_5,paste0(path_out,'\\CC_5.asc'),overwrite=T)
  writeRaster(CC_50,paste0(path_out,'\\CC_50.asc'),overwrite=T)
  
  rm(CC_5,CC_50)
  
  # CHM 5 and 50
  rcl <- matrix(c(1,1,0,NA),ncol=2,byrow=T)
  chm_bool <- reclassify(chm_binary,rcl)
  
  CHM_5 <- focal(chm_1*chm_bool,focal.5,fun='sum',na.rm=T)
  CHM_5 <- raster::resample(CHM_5,chm_1)
  CHM_5[is.na(CHM_5)] <- 0
  
  CHM_20 <- focal(chm_1*chm_bool,focal.20,fun='sum',na.rm=T)
  CHM_20 <- raster::resample(CHM_20,chm_1)
  CHM_20[is.na(CHM_20)] <- 0
  
  CHM_50 <- focal(chm_1*chm_bool,focal.50,fun='sum',na.rm=T)
  CHM_50 <- raster::resample(CHM_50,chm_1)
  CHM_50[is.na(CHM_50)] <- 0
  
  writeRaster(chm_1,paste0(path_out,'\\CHM_1.asc'),overwrite=T)
  writeRaster(CHM_5,paste0(path_out,'\\CHM_5.asc'),overwrite=T)
  writeRaster(CHM_50,paste0(path_out,'\\CHM_50.asc'),overwrite=T)
  
  # Median of Canopy Height
  CHMmed <-   focal(chm_1*chm_bool,focal.5_nonwweight,fun=median,na.rm=T)
  CHMmed <- raster::resample(CHMmed,chm_1)
  CHMmed[is.na(CHMmed)] <- 0
  
  writeRaster(CHMmed,paste0(path_out,'\\CHMmed_5.asc'),overwrite=T)
  
  rm(CHMmed)
  
  # TPI
  TPI_1 <- terrain(dtm_1,'TPI')
  TPI_5 <- focal(TPI_1,focal.5,fun='sum',na.rm=F)
  
  TPI_90 <- terrain(SRTM,'TPI')
  TPI_90 <- projectRaster(TPI_90,chm_1)
  TPI_90 <- raster::resample(TPI_90,chm_1)
  
  TPI_1 <- raster::resample(TPI_1,chm_1)
  TPI_5 <- raster::resample(TPI_5,chm_1)
  TPI_90 <- raster::resample(TPI_90,chm_1)
  
  writeRaster(TPI_1,paste0(path_out,'\\TPI_1.asc'),overwrite=T)
  writeRaster(TPI_5,paste0(path_out,'\\TPI_5.asc'),overwrite=T)
  writeRaster(TPI_90,paste0(path_out,'\\TPI_90.asc'),overwrite=T)
  
  rm(TPI_1,TPI_5,TPI_90)
  
  # TWI
  TWI_1 <- calcTWI(dtm_1)
  TWI_5 <- focal(TWI_1,focal.5,fun='sum',na.rm=F)
  
  TWI_1 <- raster::resample(TWI_1,chm_1)
  TWI_5 <- raster::resample(TWI_5,chm_1)
  
  writeRaster(TWI_1,paste0(path_out,'\\TWI_1.asc'),overwrite=T)
  writeRaster(TWI_5,paste0(path_out,'\\TWI_5.asc'),overwrite=T)
  
  rm(TWI_1,TWI_5)

  # Northness of DTM  
  ASPECT <- terrain(aggregate(dtm_1,5),'aspect',unit='radians')
  ASPECT <- raster::resample(ASPECT,chm_1)
  SLOPE <- terrain(aggregate(dtm_1,5),'slope',unit='radians')
  SLOPE <- raster::resample(SLOPE,chm_1)
  
  NORTHNESS_1 <- cos(ASPECT)*sin(SLOPE)
  NORTHNESS_5 <- focal(NORTHNESS_1,focal.5,fun='sum',na.rm=F)
  NORTHNESS_20 <- focal(NORTHNESS_1,focal.20,fun='sum',na.rm=F)
  NORTHNESS_50 <- focal(NORTHNESS_1,focal.50,fun='sum',na.rm=F)
  
  writeRaster(NORTHNESS_5,paste0(path_out,'\\NN_DTM_5.asc'),overwrite=T)
  writeRaster(NORTHNESS_20,paste0(path_out,'\\NN_DTM_20.asc'),overwrite=T)
  writeRaster(NORTHNESS_50,paste0(path_out,'\\NN_DTM_50.asc'),overwrite=T)
  
  rm(ASPECT,SLOPE,NORTHNESS_5,NORTHNESS_20,NORTHNESS_50,NORTHNESS_1)
  
  # NORTHNESS OF CHM
  ASPECT_CHM <- terrain(aggregate(CHM_5,1),'aspect',unit='radians')
  SLOPE_CHM <- terrain(aggregate(CHM_5,1),'slope',unit='radians')
  NORTHNESS_CHM_5 <- cos(ASPECT_CHM)*sin(SLOPE_CHM)
  NORTHNESS_CHM_5 <- raster::resample(NORTHNESS_CHM_5,chm_1)
  
  ASPECT_CHM <- terrain(aggregate(CHM_20,1),'aspect',unit='radians')
  SLOPE_CHM <- terrain(aggregate(CHM_20,1),'slope',unit='radians')
  NORTHNESS_CHM_20 <- cos(ASPECT_CHM)*sin(SLOPE_CHM)
  NORTHNESS_CHM_20 <- raster::resample(NORTHNESS_CHM_20,chm_1)
  
  ASPECT_CHM <- terrain(aggregate(CHM_50,1),'aspect',unit='radians')
  SLOPE_CHM <- terrain(aggregate(CHM_50,1),'slope',unit='radians')
  NORTHNESS_CHM_50 <- cos(ASPECT_CHM)*sin(SLOPE_CHM)
  NORTHNESS_CHM_50 <- raster::resample(NORTHNESS_CHM_50,chm_1)
  
  writeRaster(NORTHNESS_CHM_5,paste0(path_out,'\\NN_CHM_5.asc'),overwrite=T)
  writeRaster(NORTHNESS_CHM_20,paste0(path_out,'\\NN_CHM_20.asc'),overwrite=T)
  writeRaster(NORTHNESS_CHM_50,paste0(path_out,'\\NN_CHM_50.asc'),overwrite=T)
  
  rm(ASPECT_CHM,SLOPE_CHM,NORTHNESS_CHM_5,NORTHNESS_CHM_20,NORTHNESS_CHM_50)
  
  # DIST
  distanceToTree(chm_1,path_out)
  
  # DCE
  dce_stack_final <- calc_DCE(chm_1,step_nr = DCE_step_nr,wind_direction = wind_direction)
  dce_stack_final <- raster::resample(dce_stack_final,chm_1)
  
  writeRaster(dce_stack_final$DCE_1,paste0(path_out,'\\DCE_1.asc'),overwrite=T)
  writeRaster(dce_stack_final$NDCE_1,paste0(path_out,'\\NDCE_1.asc'),overwrite=T)
  writeRaster(dce_stack_final$SDCE_1,paste0(path_out,'\\SDCE_1.asc'),overwrite=T)
  writeRaster(dce_stack_final$WIND_DCE_1,paste0(path_out,'\\LWDCE_1.asc'),overwrite=T)
  writeRaster(dce_stack_final$WIND_DCE_2,paste0(path_out,'\\WFDCE_1.asc'),overwrite=T)
  
  # Windninja
  wind_stack <- applyWindNinja(dtm_1,SRTM,wind_direction,wind_velocity,path_out = path_out,path_to_wn_exe)
  wind_stack <- raster::resample(wind_stack,chm_1)
  
  writeRaster(wind_stack$WNDIR,paste0(path_out,'\\WNDIR_1.asc'),overwrite=T)
  writeRaster(wind_stack$WNVEL,paste0(path_out,'\\WNVEL_1.asc'),overwrite=T)
  
  
}

# Function to evaluate HS maps with reference HS maps (e.g. from LiDAR observations)
err_assessment <- function(stack_hs_test,stack_hs_reference){
  stack_hs_test<-raster::resample(stack_hs_test,stack_hs_reference)
  err_map <- stack_hs_test-stack_hs_reference
  err_map_df <- na.omit(as.data.frame(err_map))
  err_map_df <- reshape2::melt(err_map_df)
  
  test.df <-reshape2::melt(as.data.frame(stack_hs_test))$value
  ref.df <- reshape2::melt(as.data.frame(stack_hs_reference))$value
  
  df <- data.frame(test.df,ref.df)
  df <- na.omit(df)
  
  R <- cor(df$test.df,df$ref.df)
  MEA <- mean(abs(df$test.df-df$ref.df))
  RMSE <- sqrt(mean(((df$test.df-df$ref.df))^2))
  NRMSE <- RMSE/mean(df$ref.df) 
  NMEA <- MEA/mean(df$ref.df)
  
  result <-data.frame(data='hs_daily',
                      validata='Reference',
                      mean= mean(df$ref.df),
                      n=nrow(df),
                      RMSE,
                      NRMSE,
                      MEA,
                      NMEA,
                      R)
  
  
  return(result[1,]) 
}

# Function to randomly select n_points Sensor locations within each cluster!
sample_sensorlocations <- function(prob_cluster,n_points){
  require(dplyr)
  require(sp)
  
  names(prob_cluster) <- c('X1','X2','X3','X4','cluster')
  
  # Prepare prob_cluster
  prob_cluster.data <- na.omit(as.data.frame(prob_cluster,xy=T))
  
  if(nrow(prob_cluster.data)>0){
    
    # Choose n reference locations based on probabilities as weights
    for (cl in unique(prob_cluster.data$cluster)) {
      names_cl <- names(prob_cluster.data)
      names_cl <- names_cl[!(names_cl%in%c('x','y','cluster'))]
      
      cl_data <- prob_cluster.data %>%
        arrange(desc(eval(parse(text=names_cl))))  %>%
        filter(cluster==cl)
      
      cl_sample <- sample(1:nrow(cl_data),n_points)
      cl_data <- cl_data[cl_sample,]
      
      if (cl==unique(prob_cluster.data$cluster)[1]) {
        sensor_locations <- cl_data
      }else{
        sensor_locations <- rbind(sensor_locations,cl_data)
      }
      
    }
    
    sensor_locations_data <- sensor_locations
    coordinates(sensor_locations)=~x+y

    sensor_locations$ID <- paste0('cl',sensor_locations$cluster,'_nr',1:n_points)
    
    return(sensor_locations)
  }
}

# Function to interpolate local measurements based on clusters.
interpolateSensorMeasurements <- function(prob_cluster,snow_depth_raster,sensor_locations){
  require(gstat)
  
  sensor_locations$hs <- raster::extract(snow_depth_raster,sensor_locations,method='bilinear')
  sensor_locations <- sensor_locations[!is.na(sensor_locations$hs),]
  # Use the extent of the existing raster to define the model domain
  model_domain <- raster::extent(prob_cluster[[nlayers(prob_cluster)]])
  
  # Create an empty raster based on the extent for interpolation
  interpolation_grid <- raster(model_domain, resolution = res(prob_cluster[[nlayers(prob_cluster)]]))  # Define resolution as needed

  hs_cluster_map <- stack()
  # Perform IDW interpolation within the defined model domain
  for (c in 1:n_class) {

        # Perform Interpolation
    interpolation_model <- gstat(id = "hs",
                                 formula = hs ~ 1,
                                 data = sensor_locations[sensor_locations$cluster==c,],
                                 nmax=5,
                                 maxdist=500)
    
    interpolation <- interpolate(interpolation_grid, interpolation_model)
    hs_cluster_map <- stack(hs_cluster_map,interpolation)
    
  }
  
  names(hs_cluster_map)<-1:n_class
  return(hs_cluster_map)
  
}

# Function to read configuration files.
readConfigFile <- function(path){
  require(stringr)
  
  config_data <- read.delim(path,quote = "")  
  
  for (row in 1:nrow(config_data)) {
    if (stringr::str_starts(config_data[row,],"#",negate = T)){
      eval(parse(text=config_data[row,]),envir = .GlobalEnv)
      
    }
  }
  message(paste0(Sys.time(),': Successfully read configfile \n \n'))
}

# Function to read and combine feature and cluster data from all sites to one data.frame.
FeaturesToDataFrame <-function(features,target_variable,study_sites_path){
  require(raster)
  
  final_data <- data.frame()
  
  for (site in study_sites_path) {
    
    print(site)
    
    # Path-handling
    files <- list.files(site,full.names = T)
    feature_raster <- grep(paste0("(", paste(paste0(features,'.asc'), collapse = "|"), ")$"), files, value = TRUE)
    target_raster <- files[str_detect(list.files(site,full.names = T),pattern=paste0(target_variable,'$'))]
    
    # Read features and aggregate to 3 m spatial resolution
    feature_stack <- stack(feature_raster)
    feature_stack <- raster::aggregate(feature_stack,3)
    
    # Read cluster probabilities and align with feature_stack
    target_raster <- stack(target_raster)
    names(target_raster) <- c(paste0('X',1:nlayers(target_raster)))
    target_raster <- raster::resample(target_raster,feature_stack)
    
    # Merge to data frame
    data <- stack(target_raster,feature_stack)
    data <- na.omit(as.data.frame(data))
    
    final_data <- rbind(final_data,data)
    
  }
  
  return(final_data)
}
