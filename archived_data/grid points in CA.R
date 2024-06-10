
library(sf)
library(raster)
library(ggplot2)

shp <- getData('GADM', country = 'USA', level = 2) %>%
  subset(NAME_1 == 'California') %>% # 選擇加州（California）
  st_as_sf()

grid <- st_make_grid(shp, n = c(40, 40),
                     what = "centers",
                     square = TRUE) %>% st_intersection(shp)

ggplot() +
  geom_sf(data = shp) +
  geom_sf(data = grid, col = "navy")

# 提取經緯度
grid_coords <- st_coordinates(grid)

# 創建經緯度資料框
coords_df <- data.frame(lon = grid_coords[, "X"],
                        lat = grid_coords[, "Y"])

# 匯出經緯度資料框到CSV檔案
write.csv(coords_df, "D:\\研究所Meeting\\COVID19\\CA grid points\\grid_coordinates.csv", row.names = TRUE)
plot(coords_df)
