# Rotated Quadtree-Based Spatial Adaptive Partitioning for Soil Temperature Modeling

This repository contains the R implementation of the **rotated quadtree spatial partitioning framework** developed for the study:

**“Spatially adaptive estimation of multi-layer soil temperature at a daily time-step across China during 2010-2020”**

The code constructs **six rotated quadtree structures** at different orientation angles and exports the resulting polygons as shapefiles. These partitions are used to train localized machine-learning models for high-resolution soil temperature (Ts) prediction.

---

## 🔍 Overview

A quadtree is a hierarchical spatial data structure that recursively subdivides space into four quadrants.  
However, a single quadtree partition may suffer from:
  
- edge effects at grid boundaries  
- instability for boundary-area samples  

To address this, we generate **six rotated quadtree partitions** using rotation angles:
0°, 15°, 30°, 45°, 60°, 75°


---

## 🧰 Dependencies

The code is written in R and uses the following packages:

- data.table  
- AQuadtree  
- raster  
- terra  
- sf  
- ggspatial  


## 📦 Input File Format (CSV)

Each CSV file must contain the following columns:

| Column Name     | Description                                  | Type     |
|-----------------|----------------------------------------------|----------|
| Station_Id_C    | Unique station identifier                    | Integer / String |
| Lon             | Longitude (WGS84)                            | Numeric  |
| Lat             | Latitude (WGS84)                             | Numeric  |

Example CSV structure:
Station_Id_C, Lon, Lat, 
50136, 122.5167, 52.9667
50246, 124.7167, 52.3500
50349, 124.4000, 51.6667

## 📦 Output Data and Formats
All outputs are exported as ESRI **Shapefiles (.shp)**
- Six rotated-quadtree shapefiles for spatially adaptive modeling  
- CRS: WGS84  
- Each file contains hierarchical spatial partitions used for localized predictions
data_s_quadtree_rotation_1.shp
data_s_quadtree_rotation_2.shp
data_s_quadtree_rotation_3.shp
data_s_quadtree_rotation_4.shp
data_s_quadtree_rotation_5.shp
data_s_quadtree_rotation_6.shp
