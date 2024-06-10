library(sf)
library(raster)
library(ggplot2)
library(maps)
library(grid)

# load some spatial data. Administrative Boundary
aut <- getData('GADM', country = 'aut', level = 0)%>% st_as_sf()

# ggplot() + 
# geom_sf(data = aut)

grid <-st_make_grid(aut, n=c(20, 20), 
                    what = "centers", 
                    square=TRUE) %>% st_intersection(aut)                            # only within the polygon

 ggplot() + 
   geom_sf(data = aut) + 
   geom_sf(data = grid, col ="red")

 #====================================================================
shpp<-map('state', region = 'california', fill=TRUE, col="gray95",interior = T)%>% st_as_sf()
 grid <-st_make_grid(shpp, n=c(20, 20), 
                     what = "centers", 
                     square=TRUE) %>% st_intersection(shpp)   

 ggplot() + 
   geom_sf(data = shpp) + 
   geom_sf(data = grid, col ="red")

 
 states    <- c('California')
 us <- getData("GADM",country="USA",level=1)
 us.states <- us[us$NAME_1 %in% states,]%>% st_as_sf()
 us.bbox <- bbox(us.states)
 xlim <- c(min(us.bbox[1,1],ca.bbox[1,1]),max(us.bbox[1,2],ca.bbox[1,2]))
 ylim <- c(min(us.bbox[2,1],ca.bbox[2,1]),max(us.bbox[2,2],ca.bbox[2,2]))
 plot(us.states, xlim=xlim, ylim=ylim)
 plot(ca.provinces, xlim=xlim, ylim=ylim, add=T)