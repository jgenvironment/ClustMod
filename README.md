# ClustMod 
**Version 4.0**  
Author: Joschka Geissler  
Last modified: 25 January 2025

## Background
This documentation introduces the required steps for applying the ClustMod model. ClustMod is a model that allows the prediction of snow distribution patterns based on Canopy Height Models (CHM) and Digital Terrain Models (DTM). ClustMod is implemented in R (Version 4.1.0) and requires the following libraries to be installed:

### Table 1: R-packages required for running the ClustMod model

| Name             | Version | Literature                                               |
|------------------|---------|----------------------------------------------------------|
| Rstudioapi       | 0.13    | Ushey et al. (2014)                                       |
| raster           | 3.4-13  | Hijmans (2021)                                           |
| gstat            | 2.0-8   | Pebesma and Graeler (2003)                                |
| tidyverse        | 1.3.2   | Wickham et al. (2019)                                     |
| sf               | 1.0-9   | Pebesma (2018)                                           |
| whitebox         | 2.3.0   | Wu and Brown (2022) and Lindsay (2016)                   |
| assertthat       | 0.2.1   | Wickham (2013)                                           |
| ForestTools      | 0.2.5   | Plowright (2017)                                         |
| tmap             | 3.3-4   | Tennekes (2018)                                          |
| R.utils          | 2.12.2  | Bengtsson (2005)                                         |
| reshape2         | 1.4.4   | Gräler et al. (2016)                                      |
| windninjr        | 0.0.3   | Raymond (2020)                                           |
| caret            | 6.0-88  | Kuhn (2008)                                              |

Snow distribution patterns are based on the ClustSnow workflow, a workflow that performs an unsupervised classification (clustering using k-Means and random forest algorithms) of a stack of multitemporal HS maps. ClustSnow was first introduced by Geissler et al. (2023) and Geissler, Mazzotti, et al. (2024) and is available for download from Geissler and Weiler (2024).

## Workflow
The ClustMod model first determines all model features for each study site using the `PrepareData` function. This function involves the application of several third party models and algorithms: 
- The WindNinja (Raymond, 2020) model is applied to spatialize values of windspeed and direction across the study domains. WindNinja must be installed and the path to the executable must be indicated in the configuration files. Note that functions from the windninjr package are available as an independent script as an update to the windninjr package lead to problems in executing ClustMod.
- The DCE algorithm, first presented by Mazzotti et al. (2019), to determine forest structure metrics based on the distance to the canopy edges (DCE).
- A workflow to derive the Topographic Wetness Index as presented by Lindsay (2016).

For all study sites where the `isTrain` argument of the `config_data_preparation` configuration file is TRUE, clusters are derived from available HS maps using the `getCluster` function. Subsequently, one random forest model is trained for each of the four (`n_class`) clusters, based on the features (listed in `features`) as independent variables and cluster probabilities as dependent variable for all sites that are listed as `train_sites`. Note that for this version of ClustMod, only the number of clusters and therefore random forest models must be four (`n_class`=4).

The available code further allows for an evaluation of the ClustMod model. Therefore, the Overall Accuracy (OA) is derived for a small, randomly selected subset of the training dataset (size determined by 1 -  `train_test_split`). Therefore, for each grid cell, the winner clusters (cluster number with the highest probability) are compared between the predicted clusters and the cluster derived from the observations using the `getCluster` function. 

Subsequently, the uploaded code allows for a prediction of clusters for the sites that are listed in the configurations file (file `config_ClustMod`; variable `test_sites`).

Predicted clusters can serve as a basis for extrapolating point measurements or model simulations. For more information we refer to Geissler (2025).

## Setting Up the Folder and Data Structure
The steps required for running ClustMod involve setting up i) the right folder and data structure, ii) creating the individual configuration files and finally iii) run the R-scripts. This documentation will explain each of these steps in the following.

This documentation presents the version 4.0 of ClustMod, as presented in Geissler (2025). To run this specific version of ClustMod, different external datasets are required that are listed in Table 2. From these data sets, the CHM, DTM and available HS maps are required together with a coarser-level DTM containing each of the study sites. DTM, CHM and HS maps must be resampled to 1m spatial resolution. Large-scale DTMs are available worldwide, among others, from the Shuttle Radar Topography Mission (SRTM).

### Table 2: Data sets required for reproduction ClustMod as presented by Geissler (2025)

| Study Site         | Country       | Last Year of Data Acquisition | Citation                                                         |
|--------------------|---------------|-------------------------------|------------------------------------------------------------------|
| Alptal             | Switzerland   | 2023                          | Geissler, Rathmann, and Weiler (2024)                            |
| Schauinsland       | Germany       | 2023                          | Geissler, Rathmann, and Weiler (2024)                            |
| Fluela_North       | Switzerland   | 2020                          | Koutantou et al. (2022)                                          |
| Fluela_South       | Switzerland   | 2020                          | Koutantou et al. (2022)                                          |
| Fluela             | Switzerland   | 2017                          | Mazzotti et al. (2023) and Mazzotti and Jonas (2022)             |

For a successful application of ClustMod, a specific folder structure is required. The working directory must contain the ClustMod script (`ClustMod_v4.R`) as well as ClustMod’s and ClustSnow’s functions (`ClustMod_Functions_v4.R`, `ClustSnow_Functions_v4.R` and `windninjr.R`). Data sets must be stored in individual folders located within the `Data_ClustMod` folder. Besides these data sets, the `Data_ClustMod` folder contains the configuration file `config_ClustMod.txt`. 

For all data sets, the CHM and DTM must be placed in the respective folders and, if they should be considered as a training site, a stack with all HS observations is required additionally. File paths to these raster data must be specified in the `config_data_preparation` configuration file that is located within each individual data folder. Figure 1 illustrates this data and folder structure for the ClustMod model as presented by Geissler (2025).

![grafik](https://github.com/user-attachments/assets/1dc4e91b-8088-4de2-8f7d-03fb959892ba)
### Figure 1: Folder Structure for setting up the ClustMod Model

## Configuration file “config_ClustMod”
The configuration file `config_ClustMod` must be stored in the `Data_ClustMod` folder. This file contains all parameters required to define the ClustMod model. Each line in `config_ClustMod` must be written in valid R syntax to ensure it can be interpreted by R.

Key parameters in the configuration file include:
- `path_in`: Specifies the directory path where the input data is stored.
- `path_out`: Defines the directory path where the model output will be written.
- `experiment` and `model`: Specify the name of the individual model experiment.
- `train_sites`: A vector containing the names of the study sites to be used for training the ClustMod model. The names must exactly match the corresponding folder names in Data_ClustMod.
- `test_sites`: A vector containing the names of the study sites for which clusters should be predicted. These names must also match their respective folder names.
- `features`: A vector listing all the independent variables to be used during training of the ClustMod model.
- `n_class`: The number of clusters and thus random forest models the ClustMod model should be built upon (Only value 4 possible in the published version of ClustMod).
- `train_test_split`: Defines how many data points of the training data set should be used for training of the ClustMod model, and how many for the subsequent determination of the OA (testing).

The parameters `sample_length`, `kmeans_maxiter`, `kmeans_nstart`, `mtry` and `n_trees` are parameters related to the ClustSnow workflow. More details can be found in Geissler, Mazzotti, et al. (2024) and Geissler et al. (2023) including a sensitivity analysis and calibration results.

The following lines show the `config_ClustMod` configuration file for the example presented by Geissler (2025).

```r
# Parameters for ClustMod
path_in 		= “Data_ClustMod”
path_out 		= “Output”
experiment               	= "20250101"
model			= "ClustMod_v4"
train_sites		= c("Fluela_North", "Fluela_South")
test_sites		= c("Alptal", "Schauinsland")
features               	= c("CHM", "DTM", "WindSpeed", "Topographic Wetness Index")
n_class                	= 4
train_test_split      	= 0.8
sample_length         	= 2000
kmeans_maxiter        	= 100
kmeans_nstart         	= 25
mtry                    = 3
n_trees                 = 300
```

## Configuration file “config_data_preparation”

As for `config_ClustMod`, each line of the `config_data_preparation` configuration file must be written in valid R syntax. Each data set must contain such a configuration file. It contains all relevant data paths and information needed for the individual study sites.

- `finalize`: Using the parameter `finalize` can allow the applicant of the ClustMod model to select whether (TRUE) or not (FALSE) the `PrepareData` function should be applied to this study site and thus if features should be (re)calculated for this study site. This parameter therefore allows the applicant to avoid the recalculation of features after small adjustments have been made to the model.
- `dtm_path`, `dtm_largescale_path` and `chm_path`: Defines the absolute or relative (within the folder of the data set) paths of the DTM, larger-scale DTM (e.g., SRTM) and CHM.
- `crs_dtm`, `crs_dtm_largescale` and `crs_chm`: Define the coordinate systems of provided DTMs and CHMs using EPSG codes.
- `dir_out`: Defines the folder where data and logfiles should be written to.
- `buffer_size`: Defines the buffer size, in meters, by which the SRTM data should extend beyond the provided DTM. This buffer is necessary for calculating large-scale topographic parameters, such as the topographic position index (TPI), or for the application of the WindNinja model.
- `isTrain`: This parameter defines whether or not clusters should be calculated from existing HS data. Note that for all data sets listed as `train_sites`, this parameter must be TRUE. For all data sets listed in `test_sites`, this parameter should be FALSE.
- `snow_depth_path`: This parameter must indicate the path to the HS data. This is required when `isTrain` is TRUE or when the application will be performed. Otherwise, this parameter can be left out.
- `crs_snowdata`: This parameter must indicate the coordinate system of the HS data using EPSG codes. This is required when `isTrain` is TRUE or when the application will be performed. Otherwise, this parameter can be left out.
- `wind_direction`: The dominant wind direction during winter months within the domain of the data set in degrees. In Geissler (2025), this value was obtained from the ERA5 meteorological reanalysis model (Hersbach et al., 2020).
- `wind_speed`: The average wind speed during winter months within the domain of the data set in meters per second, measured at 10 m above ground. In Geissler (2025), this value was obtained from the ERA5 meteorological reanalysis model (Hersbach et al., 2020).
- `path_to_wn_exe`: For running the WindNinja model, the path to the WindNinja executable must be specified.

The following lines show the `config_data_preparation` configuration file for the example presented by Geissler (2025) and the Alptal data set.

```r
# ClustMod Parameters
dtm_path                  = ".\\RAW\\dtm_1.asc"
dtm_largescale_path       = ".\\SRTM.asc"
chm_path                  = ".\\RAW\\chm_1.asc"
crs_dtm                   = "EPSG:32632"
crs_chm                   = "EPSG:32632"
crs_dtm_largescale        = "EPSG:4326"

finalize = FALSE
buffer_size               = 90
DCE_step_nr               = 300

isTrain                   = TRUE
snow_depth_path           = ".\\snow-depth"
crs_snowdata              = "EPSG:2056"

# WindNinja Parameters
wind_direction            = 120 
wind_velocity             = 1.3
path_to_wn_exe            = "C:\\WindNinja\\WindNinja-3.6.0\\bin\\WindNinja_cli.exe"
```
## Data Output

All output of ClustMod is stored in the defined output path (`dir_out`).

For all sites listed in `test_sites`, the clusters are predicted using the trained ClustMod model. The application of ClustMod produces logfiles to store essential information regarding the model execution. These logfiles include the following details:

- **Model parameters**: Key settings and configurations used during the model's operation.
- **Runtime information**: Details about the execution process, including timestamps and processing durations.
- **Messages, warnings, and errors**: Any messages generated by R during execution, including warnings and errors, providing insights for debugging and ensuring reproducibility.

## References

- Bengtsson, H. (2005). R.utils: CRAN Contributed Packages [Software]. [https://doi.org/10.32614/CRAN.package.R.utils](https://doi.org/10.32614/CRAN.package.R.utils)
- Geissler, J. (2025). Determining and Predicting Spatiotemporal Snow Patterns in Forested Environments [In Press]. Professur für Hydrologie, University of Freiburg, Germany.
- Geissler, J., Mazzotti, G., Rathmann, L., Webster, C., & Weiler, M. (2024). ClustSnow: Utilizing temporally persistent forest snow patterns under variable environmental conditions. [https://doi.org/10.22541/essoar.172222597.78203131/v1](https://doi.org/10.22541/essoar.172222597.78203131/v1)
- Geissler, J., Rathmann, L., & Weiler, M. (2023). Combining Daily Sensor Observations and Spatial LiDAR Data for Mapping Snow Water Equivalent in a Sub‐Alpine Forest. *Water Resources Research*, 59(9), Article e2023WR034460. [https://doi.org/10.1029/2023WR034460](https://doi.org/10.1029/2023WR034460)
- Geissler, J., Rathmann, L., & Weiler, M. (2024). Spatiotemporal observations on the distribution of snow in forests for two study sites and three winter seasons [Research Data]. [https://doi.org/10.6094/UNIFR/255332](https://doi.org/10.6094/UNIFR/255332)
- Geissler, J., & Weiler, M. (2024). ClustSnow (Version v2.0) [Computer software]. GitHub. [https://github.com/jgenvironment/ClustSnow](https://github.com/jgenvironment/ClustSnow)
- Gräler, B., Pebesma, E., & Heuvelink, G. (2016). Spatio-Temporal Interpolation using gstat. *The R Journal*, 8(1), 204–218. [https://journal.r-project.org/archive/2016/RJ-2016-014/index.html](https://journal.r-project.org/archive/2016/RJ-2016-014/index.html)
- Hersbach, H., Bell, B., Berrisford, P., Hirahara, S., Horányi, A., Muñoz‐Sabater, J., Nicolas, J., Peubey, C., Radu, R., Schepers, D., Simmons, A., Soci, C., Abdalla, S., Abellan, X., Balsamo, G., Bechtold, P., Biavati, G., Bidlot, J., Bonavita, M., . . . Thépaut, J.‑N. (2020). The ERA5 global reanalysis. *Quarterly Journal of the Royal Meteorological Society*, 146(730), 1999–2049. [https://doi.org/10.1002/qj.3803](https://doi.org/10.1002/qj.3803)
- Hijmans, R. J. (2021). raster: Geographic Data Analysis and Modeling [Software]. [https://CRAN.R-project.org/package=raster](https://CRAN.R-project.org/package=raster)
- Koutantou, K., Mazzotti, G., Brunner, P., Webster, C., & Jonas, T. (2022). Exploring snow distribution dynamics in steep forested slopes with UAV-borne LiDAR. *Cold Regions Science and Technology*, 200(5), 103587. [https://doi.org/10.1016/j.coldregions.2022.103587](https://doi.org/10.1016/j.coldregions.2022.103587)
- Kuhn, M. (2008). Building Predictive Models in R Using the caret Package. *Journal of Statistical Software*, 28(5). [https://doi.org/10.18637/jss.v028.i05](https://doi.org/10.18637/jss.v028.i05)
- Lindsay, J. B. (2016). Whitebox GAT: A case study in geomorphometric analysis. *Computers & Geosciences*, 95, 75–84. [https://doi.org/10.1016/j.cageo.2016.07.003](https://doi.org/10.1016/j.cageo.2016.07.003)
- Mazzotti, G., Currier, W. R., Deems, J. S., Pflug, J. M., Lundquist, J. D., & Jonas, T. (2019). Revisiting Snow Cover Variability and Canopy Structure Within Forest Stands: Insights From Airborne Lidar Data. *Water Resources Research*, 55(7), 6198–6216. [https://doi.org/10.1029/2019WR024898](https://doi.org/10.1029/2019WR024898)
- Mazzotti, G., & Jonas, T. (2022). Input datasets for forest snow modelling in Fluela valley, WY 2016-21. [https://doi.org/10.16904/envidat.338](https://doi.org/10.16904/envidat.338)
- Mazzotti, G., Webster, C., Quéno, L., Cluzet, B., & Jonas, T. (2023). Canopy structure, topography and weather are equally important drivers of small-scale snow cover dynamics in sub-alpine forests. *Hydrology and Earth System Sciences*, 27, 2099–2121. [https://doi.org/10.5194/hess-2022-273](https://doi.org/10.5194/hess-2022-273)
- Pebesma, E. (2018). Simple Features for R: Standardized Support for Spatial Vector Data. *The R Journal*, 10(1), 439–446. [https://doi.org/10.32614/RJ-2018-009](https://doi.org/10.32614/RJ-2018-009)
- Pebesma, E., & Graeler, B. (2003). gstat: CRAN Contributed Packages [Software]. [https://doi.org/10.32614/CRAN.package.gstat](https://doi.org/10.32614/CRAN.package.gstat)
- Plowright, A. (2017). ForestTools: CRAN Contributed Packages [Software]. [https://doi.org/10.32614/CRAN.package.ForestTools](https://doi.org/10.32614/CRAN.package.ForestTools)
- Raymond, B. (2020). windninjr (Version 0.0.3) [R-Package]. [https://github.com/SCAR-sandpit/windninjr/tree/main](https://github.com/SCAR-sandpit/windninjr/tree/main)
- Tennekes, M. (2018). tmap: Thematic Maps in R. *Journal of Statistical Software*, 84(6), 1–39. [https://doi.org/10.18637/jss.v084.i06](https://doi.org/10.18637/jss.v084.i06)
- Ushey, K., Allaire, J. J., Wickham, H., & Ritchie, G. (2014). rstudioapi: CRAN Contributed Packages [Software]. [https://doi.org/10.32614/CRAN.package.rstudioapi](https://doi.org/10.32614/CRAN.package.rstudioapi)
- Wickham, H. (2013). assertthat: CRAN Contributed Packages [Software]. [https://doi.org/10.32614/CRAN.package.assertthat](https://doi.org/10.32614/CRAN.package.assertthat)
- Wickham, H., Averick, M., Bryan, J., Chang, W., McGowan, L. D., François, R., Grolemund, G., Hayes, A., Henry, L., Hester, J., Kuhn, M., Pedersen, T. L., Miller, E., Bache, S. M., Müller, K., Ooms, J., Robinson, D., Seidel, D. P., Spinu, V., . . . Yutani, H. (2019). Welcome to the tidyverse. *Journal of Open Source Software*, 4(43), 1686. [https://doi.org/10.21105/joss.01686](https://doi.org/10.21105/joss.01686)
- Wu, Q., & Brown, A. (2022). whitebox (Version 2.3) [R-Package]. [https://CRAN.R-project.org/package=whitebox](https://CRAN.R-project.org/package=whitebox)
