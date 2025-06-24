# STDM
This GitHub repository contains codes and data to reproduce simulation results, figures, and real data analysis results from the paper "Structural Identification for Spatio-Temporal Dynamic Models".

## Package Environments
The codes in this repository requires R version 4.5.0 and the following R packages:
1. ```splines 4.5.0```
2. ```glmnet 4.1.8```
3. ```KRLS 1.0.0```
4. ```mosum 1.2.7```
5. ```MASS 7.3.65```
6. ```cluster 2.1.8.1```
7. ```fossil 0.4.0```
8. ```arrow 20.0.0.2```
9. ```dplyr 2.5.0```
10. ```ggplot2 3.5.2```
11. ```doParallel 1.0.17```
12. ```foreach 1.5.2```

The codes also requires Python version 3.12.7 and the following Python packages:
1. ```geopandas 1.01```
2. ```pandas 2.2.2```
3. ```matplotlib 3.9.2```

## simulation-code
This folder contains the codes to reproduce the full simulation results presented in Section 5 and Appendix B.1 in the paper. As the full simulation results invovles multiple competing methods, various settings, and 200 replications for each setting, the codes in Section 5.2 utilized parallel computing in order to reduce the running time. We would recommend to reproduce the full simulation results with a computing system of at least 16 cores and 500GB of storage. The following list provides the guidance for each file in the folder.
1. Under ```Part 1``` folder, the R code will generate results for part of Section 5.2 and Appendix B.1
   - ```Figure 1 & Table 1 & Table B1.R``` will plot Figure 1 and generate the tables for Table 1 and Table B1. It will also generate intermediate temporary data which is stored under the folder ```Results/Simulation-code/Part 1```.
2. Under ```Part 2``` folder, the R code will generate results for part of Section 5.2
   - ```Table 2.R``` will Table 2. It will also generate intermediate temporary data which is stored under the folder ```Results/Simulation-code/Part 2```.
3. Under ```Part 3``` folder, the R code will generate results for part of Section 5.3 and Appendix B.1
   - ```Table 3 & Figure 2 & Figure 3 & Table B2.R``` will plot Figure 2 and Figure 3, and generate the tables for Table 3 and Table B2. It will also generate intermediate temporary data which is stored under the folder ```Results/Simulation-code/Part 3```.

## real-data-code
This folder contains the dataset and the codes to reproduce the real data analysis results presented in Section 6 and Appendix B.2 in the paper.
1. Under ```data``` folder, it contains the orginal dataset and generated dataset for the real data analysis.
   - ```taxi_zones``` folder contains necessary location information. The dataset is publicly avaliable at [TLC Trip Record Data](https://www.nyc.gov/site/tlc/about/tlc-trip-record-data.page). Please download the data file "Taxi Zone Shapefile (PARQUET)".
   - ```combined_data.csv``` is created by ```data_combine.R```. We combined all the data documents into one single file.
   - 
