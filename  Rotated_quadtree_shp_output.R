# Clean workspace
rm(list = ls())

library(data.table)
library(AQuadtree)
library(raster)
library(terra)
library(ggspatial)
library(sf)

# 获取文件夹中所有 CSV 文件的路径
file_list <- list.files(path = "E:/Table", 
                        pattern = "\\.csv$", full.names = TRUE)

# 读取所有 CSV 文件并合并
df_train <- rbindlist(lapply(file_list, fread))

# 查看合并后的数据
head(df_train)
unique_station_count <- length(unique(df_train$Station_Id_C))
print(unique_station_count)

# aggregate到每个点
d01 <- df_train[, .(Lon = Lon[1], Lat = Lat[1]), by = Station_Id_C][, Station_Id_C := NULL]
colnames(d01)<-c("x","y")
d01 <- d01[order(x,y)]



# 设定旋转角度，共旋转6次
xuanzhuan <- c(0,pi/12,2*pi/12,3*pi/12,4*pi/12,5*pi/12)
d01 <- as.data.frame(d01)


################################################
# 建立六次旋转的四叉树分区
data_scs <- c()
for(i in c(1:6)){
  xz <- xuanzhuan[i]
  # 数据点进行旋转
  d01$x1 <- sqrt((d01$x)^2+(d01$y)^2)*
    ((d01$x/sqrt((d01$x)^2+(d01$y)^2))*cos(xz)-
       (d01$y/sqrt((d01$x)^2+(d01$y)^2))*sin(xz))
  d01$y1 <- sqrt((d01$x)^2+(d01$y)^2)*
    ((d01$y/sqrt((d01$x)^2+(d01$y)^2))*cos(xz)+
       (d01$x/sqrt((d01$x)^2+(d01$y)^2))*sin(xz))
  
  d01$x1 <- (d01$x1+360)*10
  d01$y1 <- (d01$y1+360)*10
  
  # 转化为spdataframe数据类型
  xy <- as.matrix(subset(d01, select = c(x1, y1))) 
  x <- SpatialPointsDataFrame(xy,d01[,3:4],match.ID=TRUE,
                              proj4string=CRS(as.character("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")), bbox = NULL)
  
  # 建立四叉树
  bcn.QT2 <- AQuadtree(x,
                       # colnames=c('pre'), 
                       # funs=c('mean'),
                       # thresholdField=c('pre'),
                       #threshold=5,
                       colnames = NULL,
                       threshold = 30,
                       thresholdField = NULL,
                       funs = NULL,
                       dim=1300,
                       layers=5)
  
  b <- bcn.QT2@data
  b <- order(b$level)
  # 打印当前四叉树分区的数量
  cat("旋转角度", xuanzhuan[i], "对应的四叉树分区数为:", length(bcn.QT2@polygons), "\n")
  # 保存每个四叉树四个顶点坐标
  for(k in 1:length(b)){
    j <- b[k]
    d <- c()
    d$n <- i
    # 四个顶点坐标
    d$r <- sqrt(bcn.QT2@polygons[[j]]@area)
    d$x <- (bcn.QT2@polygons[[j]]@labpt[1]-d$r/2)/10-360
    d$y <- (bcn.QT2@polygons[[j]]@labpt[2]-d$r/2)/10-360
    d$x1 <- (bcn.QT2@polygons[[j]]@labpt[1]+d$r/2)/10-360
    d$y1 <- (bcn.QT2@polygons[[j]]@labpt[2]+d$r/2)/10-360
    d$x2 <- d$x
    d$y2 <- d$y1
    d$x3 <- d$x1
    d$y3 <- d$y
    # 四个顶点坐标旋转回去
    d$x12 <- sqrt((d$x1)^2+(d$y1)^2)*
      ((d$x1/sqrt((d$x1)^2+(d$y1)^2))*cos(2*pi-xz)-
         (d$y1/sqrt((d$x1)^2+(d$y1)^2))*sin(2*pi-xz))
    d$y12 <- sqrt((d$x1)^2+(d$y1)^2)*
      ((d$y1/sqrt((d$x1)^2+(d$y1)^2))*cos(2*pi-xz)+
         (d$x1/sqrt((d$x1)^2+(d$y1)^2))*sin(2*pi-xz))
    d$x11 <- sqrt((d$x)^2+(d$y)^2)*
      ((d$x/sqrt((d$x)^2+(d$y)^2))*cos(2*pi-xz)-
         (d$y/sqrt((d$x)^2+(d$y)^2))*sin(2*pi-xz))
    d$y11 <- sqrt((d$x)^2+(d$y)^2)*
      ((d$y/sqrt((d$x)^2+(d$y)^2))*cos(2*pi-xz)+
         (d$x/sqrt((d$x)^2+(d$y)^2))*sin(2*pi-xz))
    d$x13 <- sqrt((d$x2)^2+(d$y2)^2)*
      ((d$x2/sqrt((d$x2)^2+(d$y2)^2))*cos(2*pi-xz)-
         (d$y2/sqrt((d$x2)^2+(d$y2)^2))*sin(2*pi-xz))
    d$y13 <- sqrt((d$x2)^2+(d$y2)^2)*
      ((d$y2/sqrt((d$x2)^2+(d$y2)^2))*cos(2*pi-xz)+
         (d$x2/sqrt((d$x2)^2+(d$y2)^2))*sin(2*pi-xz))
    d$x14 <- sqrt((d$x3)^2+(d$y3)^2)*
      ((d$x3/sqrt((d$x3)^2+(d$y3)^2))*cos(2*pi-xz)-
         (d$y3/sqrt((d$x3)^2+(d$y3)^2))*sin(2*pi-xz))
    d$y14 <- sqrt((d$x3)^2+(d$y3)^2)*
      ((d$y3/sqrt((d$x3)^2+(d$y3)^2))*cos(2*pi-xz)+
         (d$x3/sqrt((d$x3)^2+(d$y3)^2))*sin(2*pi-xz))
    # 保存顶点数据
    data_scs <- as.data.frame(rbind(data_scs,d))
  }
  print(i)
}
data_scs <- as.data.frame(data_scs)
# 将指定的列转换为数值类型
cols_to_convert <- c("n", "r", "x", "y", "x1", "y1", "x2", "y2", "x3", "y3", 
                     "x11", "y11", "x12", "y12", "x13", "y13", "x14", "y14")

data_scs[cols_to_convert] <- lapply(data_scs[cols_to_convert], as.numeric)


# 创建矩形框的多边形对象
polygons_list <- list()
# 批量导出六次旋转的矩形框
for (i in 1:6) {
  # 选择当前旋转角度对应的矩形框数据
  data_s <- data_scs[which(data_scs$n == i), ]
  
  # 创建矩形框的多边形对象
  polygons_list <- list()
  for (j in 1:nrow(data_s)) {
    # 创建矩形框的顶点顺序（左下 -> 左上 -> 右上 -> 右下 -> 回到左下）
    coords <- matrix(c(data_s$x12[j], data_s$y12[j],  # 左下角
                       data_s$x14[j], data_s$y14[j],  # 左上角
                       data_s$x11[j], data_s$y11[j],  # 右上角
                       data_s$x13[j], data_s$y13[j],  # 右下角
                       data_s$x12[j], data_s$y12[j]), # 回到左下角（闭合矩形）
                     ncol = 2, byrow = TRUE)
    
    # 将这些顶点转为一个多边形
    polygons_list[[j]] <- st_polygon(list(coords))
  }
  
  # 将所有矩形框转换为 sf 对象，并指定坐标系
  data_s_sf <- st_sf(geometry = st_sfc(polygons_list), crs = 4326)  # 假设坐标系为 WGS 84
  
  # 添加其他需要的信息列
  data_s_sf <- cbind(data_s_sf, data_s[, c("n", "r")])
  
  # 创建输出文件名（每次旋转对应一个不同的文件名）
  output_filename <- paste0("D:/AQuadtree_out/data_s_quadtree_rotation_", i, ".shp")
  
  # 将矩形框导出为 shapefile
  st_write(data_s_sf, output_filename, delete_dsn = TRUE)
  
  cat("已导出旋转角度", i, "的矩形框到文件:", output_filename, "\n")
}








