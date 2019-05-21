#TSA R code -------
library(ggplot2)
library(ggmap)
library(maps)
library(mapproj)
library(sp)  # classes for spatial data
library(raster)  # grids, rasters
library(rasterVis)  # raster visualisation
library(maptools)
library(rgeos)
library(dismo)
library(gstat)
attach(TRMMENSO)
library(lubridate)
byyear <- aggregate(list(rainanom, alt.x), by=list(Year,lon,lat), data=STFDF,FUN=mean)
names(byyear)
library(plyr)
byyear=rename(byyear, c("Group.1"="year","Group.2" ="lon","Group.3"="lat", "c.71.2159751279424...14.5643359200576...8.9899486560576..6.1612728239424.."="anom","c.382..483..1049..959..67..452..309..336..539..1007..1075..464.."="alt"))
#Year 1998-------
year1998 <-subset(byyear, year==1998)
attach(year1998)
year1998$x1=lat.standard=(2*lat-(max(lat)+min(lat)))/(max(lat)-min(lat))
year1998$x2=lon.standard=(2*lon-(max(lon)+min(lon)))/(max(lon)-min(lon))

year1998$X2=year1998$x2^2
year1998$X1=year1998$x1^2
year1998$X13=year1998$x1^3
year1998$X23=year1998$x2^3
year1998$X2x1=year1998$X2*year1998$x1
year1998$X1x2=year1998$X1*year1998$x2

rain1998mean=lm(anom~1, year1998)
anova(rain1998mean)
summary(rain1998mean)

rain1998linear=lm(anom~x1+x2, year1998)
anova(rain1998linear)
summary(rain1998linear)

rain1998quad=lm(anom~x1*x2+X1+X2, year1998)
anova(rain1998quad)
summary(rain1998quad)

rain1998quadalt=lm(anom~x1*x2+X1+X2+alt, year1998)
anova(rain1998quadalt)
summary(rain1998quadalt)

rain1998cubic=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2, year1998)
anova(rain1998cubic)
summary(rain1998cubic)

rain1998cubicalt=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2+alt, year1998)
anova(rain1998cubicalt)
summary(rain1998cubicalt)


year1998$fitted1=fitted(rain1998cubic)
year1998$resid1=resid(rain1998cubic) #List of residuals
plot(density(resid(rain1998cubic))) #A density plot
qqnorm(resid(rain1998cubic)) # A quantile normal plot - good for checking normality
qqline(resid(rain1998cubic))
year1998fit=data.frame(cbind(year1998$lon, year1998$lat,year1998$fitted1))
year1998res=data.frame(cbind(year1998$lon, year1998$lat,year1998$resid1))

anova(rain1998mean,rain1998linear, rain1998quad, rain1998cubic, rain1998cubicalt)

library(lattice)
w98=wireframe(X3~X1*X2,data=year1998fit,
          xlab = "Longitude", ylab = "Latitude", zlab="yhat", 
          main = "1998", scales = list(arrows = FALSE, distance=2),
          drape = TRUE,aspect = c(70/87, 0.41),
          colorkey = TRUE, screen = list(z = 30, x = -25, y = 10))
w98
c98= contourplot(X3~X1*X2,data=year1998fit,
            cuts = 20, region = TRUE,
            xlab = "Longitude",
            ylab = "Latitude",
            main = "1998")
c98 

attach(year1998fit)
X1=X2=seq(-1,1,.1)
model=function(X1,X2){13.393-7.824*X1 -138.389*X2-60.968 *X1^2-28.688*X2^2
+28.689*X1^3+66.610*X2^3-64.299*X2^2*X1+19.485*X1^2*X2-46.510*X1*X2}
model(X1,X2)

z=outer(X1,X2,model)
contour(X1,X2,z,nlevels=20, main="1998")



wireframe(X3~X1*X2,data=year1998res,
          xlab = "Lon", ylab = "Lat", zlab="ehat", 
          main = "Precip.residuals year 1998", scales = list(arrows = FALSE, distance=2),
          drape = TRUE,aspect = c(70/87, 0.22),
          colorkey = TRUE, screen = list(z = 50, x = -55, y = 2))

## Kriging 1998---- 
coordinates(year1998res)=~X1 +X2
proj4string(year1998res) <- "+proj=utm +zone=37 +datum=WGS84 +units=m"
class(year1998res)
summary(year1998res)

# For working with spatial (and spatio-temporal) data, we use the gstat package, which includes functionality for kriging, among other many things.
library(sp)
library(gstat)
# packages for manipulation & visualization
suppressPackageStartupMessages({
  library(dplyr) # for "glimpse"
  library(ggplot2)
  library(scales) # for "comma"
  library(magrittr)
})

# load up area shape file:
library(maptools)
library(rgdal)
library(sp) 
### Read spatial shapefile (Ethiopian adm0 from Global Admintrative area)
area <-readRDS("gadm36_ETH_0_sp.rds")

##function summary provides us with some useful information
summary(area)
##manually define the coordinate system by setting the coordinate system information 
plot(area)
library(raster)
library(rgdal) # for spTransform
area<-spTransform(area, CRS("+proj=longlat +zone=37 +datum=WGS84 +units=m"))

grd=spsample(area, n = 1272, "regular")
gridded(grd) = TRUE # Make it a grid
class(grd)

grd<-spTransform(grd, CRS("+proj=utm +zone=37 +datum=WGS84 +units=m"))
summary(grd)
plot(grd)

##a variogram could be fit as simply as the following code:
library(gstat)
library(sp)

vgm <- variogram(X3~1, year1998res) # calculates sample variogram values 
plot(vgm)
library(gstat)
vgm.fit98<-fit.variogram(vgm,vgm(180, "Sph", 2.1,0)) # fit model
plot(vgm, vgm.fit98)

# Arrange the data for the ggplot2 plot
# add the semivariance values of v2 to v1
Fitted <- data.frame(dist = seq(0.12, max(vgm$dist), length = 1272))
Fitted$Spherical <- variogramLine(vgm.fit98, dist_vector = Fitted$dist)$gamma
#convert the dataframes to a long format
library(reshape2)
Empirical <- melt(vgm, id.vars = "dist", measure.vars = "gamma")
Modeled <- melt(Fitted, id.vars = "dist", measure.vars = "Spherical")
library(ggplot2)
#both variogram on the same plot with different colours
v98=ggplot(Modeled, aes(x = dist, y = value, colour = variable )) + 
  geom_line() +  
  xlab("Distance") +
  ylab("Semivariance") +   ggtitle("1998")+ 
  geom_point(data = Empirical)
v98

library(magrittr)
library(gstat)
library(sp)
rain.oK98<- krige(X3 ~ 1, year1998res, grd,model = vgm.fit98)
x <- krige.cv(X3 ~ 1, year1998res, grd,model = vgm.fit98, nfold=100)
#RMSE (root mean squared error)
sqrt(mean(x$residual^2)) 

library(ggplot2)
library(scales)
k98ps=rain.oK98 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKpred"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k98ps=k98ps+xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("1998")
k98ps 

k98sd=rain.oK98 %>% as.data.frame %>%
    ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.var)) + coord_equal() +
    scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKerror"))+
    scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
    theme_bw() 
k98sd=k98sd+  xlab("Longitude") +
    ylab("Latitude") +
    ggtitle("1998")
  k98sd

# Some parameters for plotting
par(font.main=2, cex.main=1.5,cex.lab=0.4, cex.sub=1)
# Use the ColorBrewer library for color ramps
library(RColorBrewer)
precip.pal = colorRampPalette(c("yellow", "white", "blue"))

#### Countor map with values : Cokriging
ok98ps=spplot(rain.oK98,  zcol='var1.pred',col.regions=precip.pal, contour=TRUE, scales=list(draw=F),
  col='black', pretty=TRUE,main="1998", 
  labels = list(cex = 0.7), label.style = 'align', margin = TRUE,   width = 2, cex = 1.5,  xlab="Longitude",
    ylab="Latitude")
ok98ps


#Year 1999-------
year1999 <-subset(byyear, year==1999)
attach(year1999)
year1999$x1=lat.standard=(2*lat-(max(lat)+min(lat)))/(max(lat)-min(lat))
year1999$x2=lon.standard=(2*lon-(max(lon)+min(lon)))/(max(lon)-min(lon))

year1999$X2=year1999$x2^2
year1999$X1=year1999$x1^2
year1999$X13=year1999$x1^3
year1999$X23=year1999$x2^3
year1999$X2x1=year1999$X2*year1999$x1
year1999$X1x2=year1999$X1*year1999$x2

rain1999mean=lm(anom~1, year1999)
anova(rain1999mean)
summary(rain1999mean)

rain1999linear=lm(anom~x1+x2, year1999)
anova(rain1999linear)
summary(rain1999linear)

rain1999quad=lm(anom~x1*x2+X1+X2, year1999)
anova(rain1999quad)
summary(rain1999quad)

rain1999quadalt=lm(anom~x1*x2+X1+X2+alt, year1999)
anova(rain1999quadalt)
summary(rain1999quadalt)

rain1999cubic=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2, year1999)
anova(rain1999cubic)
summary(rain1999cubic)

rain1999cubicalt=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2+alt, year1999)
anova(rain1999cubicalt)
summary(rain1999cubicalt)


year1999$fitted1=fitted(rain1999cubic)
year1999$resid1=resid(rain1999cubic) #List of residuals
plot(density(resid(rain1999cubic))) #A density plot
qqnorm(resid(rain1999cubic)) # A quantile normal plot - good for checking normality
qqline(resid(rain1999cubic))
year1999fit=data.frame(cbind(year1999$lon, year1999$lat,year1999$fitted1))
year1999res=data.frame(cbind(year1999$lon, year1999$lat,year1999$resid1))

anova(rain1999mean,rain1999linear, rain1999quad, rain1999cubic, rain1999cubicalt)

library(lattice)
w99=wireframe(X3~X1*X2,data=year1999fit,
              xlab = "Longitude", ylab = "Latitude", zlab="yhat", 
              main = "1999", scales = list(arrows = FALSE, distance=2),
              drape = TRUE,aspect = c(70/87, 0.41),
              colorkey = TRUE, screen = list(z = 30, x = -25, y = 10))
w99
c99= contourplot(X3~X1*X2,data=year1999fit,
                 cuts = 20, region = TRUE,
                 xlab = "Longitude",
                 ylab = "Latitude",
                 main = "1999")
c99 

attach(year1999fit)
X1=X2=seq(-1,1,.1)
model=function(X1,X2){13.393-7.824*X1 -138.389*X2-60.968 *X1^2-28.688*X2^2
  +28.689*X1^3+66.610*X2^3-64.299*X2^2*X1+19.485*X1^2*X2-46.510*X1*X2}
model(X1,X2)

z=outer(X1,X2,model)
contour(X1,X2,z,nlevels=20, main="1999")

wireframe(X3~X1*X2,data=year1999res,
          xlab = "Lon", ylab = "Lat", zlab="ehat", 
          main = "Precip.residuals year 1999", scales = list(arrows = FALSE, distance=2),
          drape = TRUE,aspect = c(70/87, 0.22),
          colorkey = TRUE, screen = list(z = 50, x = -55, y = 2))

## Kriging 1999---- 
coordinates(year1999res)=~X1 +X2
proj4string(year1999res) <- "+proj=utm +zone=37 +datum=WGS84 +units=m"
class(year1999res)
summary(year1999res)

# For working with spatial (and spatio-temporal) data, we use the gstat package, which includes functionality for kriging, among other many things.
library(sp)
library(gstat)
# packages for manipulation & visualization
suppressPackageStartupMessages({
  library(dplyr) # for "glimpse"
  library(ggplot2)
  library(scales) # for "comma"
  library(magrittr)
})


##a variogram could be fit as simply as the following code:
vgm <- variogram(X3~1, year1999res) # calculates sample variogram values 
plot(vgm)
library(gstat)
vgm.fit99<- fit.variogram(vgm, model=vgm(187, "Sph", 2.1,0)) # fit model
plot(vgm, vgm.fit99)

# Arrange the data for the ggplot2 plot
# add the semivariance values of v2 to v1
Fitted <- data.frame(dist = seq(0.12, max(vgm$dist), length = 1272))
Fitted$Spherical <- variogramLine(vgm.fit99, dist_vector = Fitted$dist)$gamma
#convert the dataframes to a long format
library(reshape2)
Empirical <- melt(vgm, id.vars = "dist", measure.vars = "gamma")
Modeled <- melt(Fitted, id.vars = "dist", measure.vars = "Spherical")
library(ggplot2)
#both variogram on the same plot with different colours
v99=ggplot(Modeled, aes(x = dist, y = value, colour = variable )) + 
  geom_line() +  
  xlab("Distance") +
  ylab("Semivariance") +   ggtitle("1999")+ 
  geom_point(data = Empirical)
v99

library(magrittr)
library(gstat)
library(sp)
rain.oK99<- krige(X3 ~ 1, year1999res, grd,model = vgm.fit99)
x <- krige.cv(X3 ~ 1, year1999res, grd,model = vgm.fit99, nfold=100)
#RMSE (root mean squared error)
sqrt(mean(x$residual^2)) 

library(ggplot2)
library(scales)
k99ps=rain.oK99 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKpred"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k99ps=k99ps+xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("1999")
k99ps 

k99sd=rain.oK99 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.var)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKerror"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k99sd=k99sd+  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("1999")
k99sd

# Some parameters for plotting
par(font.main=2, cex.main=1.5,cex.lab=0.4, cex.sub=1)
# Use the ColorBrewer library for color ramps
library(RColorBrewer)
precip.pal = colorRampPalette(c("yellow", "white", "blue"))

#### Countor map with values : Cokriging
ok99ps=spplot(rain.oK99,  zcol='var1.pred',col.regions=precip.pal, contour=TRUE, scales=list(draw=T),
              col='black', pretty=TRUE,main="1999", 
              labels = list(cex = 0.7), label.style = 'align', margin = TRUE,   width = 2, cex = 1.5,  xlab="Longitude",
              ylab="Latitude")
ok99ps
#Year 2000-------
year2000 <-subset(byyear, year==2000)
attach(year2000)
year2000$x1=lat.standard=(2*lat-(max(lat)+min(lat)))/(max(lat)-min(lat))
year2000$x2=lon.standard=(2*lon-(max(lon)+min(lon)))/(max(lon)-min(lon))

year2000$X2=year2000$x2^2
year2000$X1=year2000$x1^2
year2000$X13=year2000$x1^3
year2000$X23=year2000$x2^3
year2000$X2x1=year2000$X2*year2000$x1
year2000$X1x2=year2000$X1*year2000$x2

rain2000mean=lm(anom~1, year2000)
anova(rain2000mean)
summary(rain2000mean)

rain2000linear=lm(anom~x1+x2, year2000)
anova(rain2000linear)
summary(rain2000linear)

rain2000quad=lm(anom~x1*x2+X1+X2, year2000)
anova(rain2000quad)
summary(rain2000quad)

rain2000quadalt=lm(anom~x1*x2+X1+X2+alt, year2000)
anova(rain2000quadalt)
summary(rain2000quadalt)

rain2000cubic=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2, year2000)
anova(rain2000cubic)
summary(rain2000cubic)

rain2000cubicalt=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2+alt, year2000)
anova(rain2000cubicalt)
summary(rain2000cubicalt)


year2000$fitted1=fitted(rain2000cubic)
year2000$resid1=resid(rain2000cubic) #List of residuals
plot(density(resid(rain2000cubic))) #A density plot
qqnorm(resid(rain2000cubic)) # A quantile normal plot - good for checking normality
qqline(resid(rain2000cubic))
year2000fit=data.frame(cbind(year2000$lon, year2000$lat,year2000$fitted1))
year2000res=data.frame(cbind(year2000$lon, year2000$lat,year2000$resid1))

anova(rain2000mean,rain2000linear, rain2000quad, rain2000cubic, rain2000cubicalt)

library(lattice)
w00=wireframe(X3~X1*X2,data=year2000fit,
              xlab = "Longitude", ylab = "Latitude", zlab="yhat", 
              main = "2000", scales = list(arrows = FALSE, distance=2),
              drape = TRUE,aspect = c(70/87, 0.41),
              colorkey = TRUE, screen = list(z = 30, x = -25, y = 10))
w00
c00= contourplot(X3~X1*X2,data=year2000fit,
                 cuts = 20, region = TRUE,
                 xlab = "Longitude",
                 ylab = "Latitude",
                 main = "2000")
c00 

attach(year2000fit)
X1=X2=seq(-1,1,.1)
model=function(X1,X2){13.393-7.824*X1 -138.389*X2-60.968 *X1^2-28.688*X2^2
  +28.689*X1^3+66.610*X2^3-64.200*X2^2*X1+19.485*X1^2*X2-46.510*X1*X2}
model(X1,X2)

z=outer(X1,X2,model)
contour(X1,X2,z,nlevels=20, main="2000")



wireframe(X3~X1*X2,data=year2000res,
          xlab = "Lon", ylab = "Lat", zlab="ehat", 
          main = "Precip.residuals year 2000", scales = list(arrows = FALSE, distance=2),
          drape = TRUE,aspect = c(70/87, 0.22),
          colorkey = TRUE, screen = list(z = 50, x = -55, y = 2))

## Kriging 2000---- 
coordinates(year2000res)=~X1 +X2
proj4string(year2000res) <- "+proj=utm +zone=37 +datum=WGS84 +units=m"
class(year2000res)
summary(year2000res)

# For working with spatial (and spatio-temporal) data, we use the gstat package, which includes functionality for kriging, among other many things.
library(sp)
library(gstat)
# packages for manipulation & visualization
suppressPackageStartupMessages({
  library(dplyr) # for "glimpse"
  library(ggplot2)
  library(scales) # for "comma"
  library(magrittr)
})


##a variogram could be fit as simply as the following code:
vgm <- variogram(X3~1, year2000res) # calculates sample variogram values 
plot(vgm)
library(gstat)
vgm.fit00<- fit.variogram(vgm, model=vgm(187, "Sph", 2.1,0)) # fit model
plot(vgm, vgm.fit00)

# Arrange the data for the ggplot2 plot
# add the semivariance values of v2 to v1
Fitted <- data.frame(dist = seq(0.12, max(vgm$dist), length = 1272))
Fitted$Spherical <- variogramLine(vgm.fit00, dist_vector = Fitted$dist)$gamma
#convert the dataframes to a long format
library(reshape2)
Empirical <- melt(vgm, id.vars = "dist", measure.vars = "gamma")
Modeled <- melt(Fitted, id.vars = "dist", measure.vars = "Spherical")
library(ggplot2)
#both variogram on the same plot with different colours
v00=ggplot(Modeled, aes(x = dist, y = value, colour = variable )) + 
  geom_line() +  
  xlab("Distance") +
  ylab("Semivariance") +   ggtitle("2000")+ 
  geom_point(data = Empirical)
v00

library(magrittr)
library(gstat)
library(sp)
rain.oK00<- krige(X3 ~ 1, year2000res, grd,model = vgm.fit00)
x <- krige.cv(X3 ~ 1, year2000res, grd,model = vgm.fit00, nfold=100)
#RMSE (root mean squared error)
sqrt(mean(x$residual^2)) 

library(ggplot2)
library(scales)
k00ps=rain.oK00 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKpred"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k00ps=k00ps+xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2000")
k00ps 

k00sd=rain.oK00 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.var)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKerror"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k00sd=k00sd+  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2000")
k00sd

# Some parameters for plotting
par(font.main=2, cex.main=1.5,cex.lab=0.4, cex.sub=1)
# Use the ColorBrewer library for color ramps
library(RColorBrewer)
precip.pal = colorRampPalette(c("yellow", "white", "blue"))

#### Countor map with values : Cokriging
ok00ps=spplot(rain.oK00,  zcol='var1.pred',col.regions=precip.pal, contour=TRUE, scales=list(draw=T),
              col='black', pretty=TRUE,main="2000", 
              labels = list(cex = 0.7), label.style = 'align', margin = TRUE,   width = 2, cex = 1.5,  xlab="Longitude",
              ylab="Latitude")
ok00ps
#Year 2001-------
year2001 <-subset(byyear, year==2001)
attach(year2001)
year2001$x1=lat.standard=(2*lat-(max(lat)+min(lat)))/(max(lat)-min(lat))
year2001$x2=lon.standard=(2*lon-(max(lon)+min(lon)))/(max(lon)-min(lon))

year2001$X2=year2001$x2^2
year2001$X1=year2001$x1^2
year2001$X13=year2001$x1^3
year2001$X23=year2001$x2^3
year2001$X2x1=year2001$X2*year2001$x1
year2001$X1x2=year2001$X1*year2001$x2

rain2001mean=lm(anom~1, year2001)
anova(rain2001mean)
summary(rain2001mean)

rain2001linear=lm(anom~x1+x2, year2001)
anova(rain2001linear)
summary(rain2001linear)

rain2001quad=lm(anom~x1*x2+X1+X2, year2001)
anova(rain2001quad)
summary(rain2001quad)

rain2001quadalt=lm(anom~x1*x2+X1+X2+alt, year2001)
anova(rain2001quadalt)
summary(rain2001quadalt)

rain2001cubic=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2, year2001)
anova(rain2001cubic)
summary(rain2001cubic)

rain2001cubicalt=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2+alt, year2001)
anova(rain2001cubicalt)
summary(rain2001cubicalt)


year2001$fitted1=fitted(rain2001cubic)
year2001$resid1=resid(rain2001cubic) #List of residuals
plot(density(resid(rain2001cubic))) #A density plot
qqnorm(resid(rain2001cubic)) # A quantile normal plot - good for checking normality
qqline(resid(rain2001cubic))
year2001fit=data.frame(cbind(year2001$lon, year2001$lat,year2001$fitted1))
year2001res=data.frame(cbind(year2001$lon, year2001$lat,year2001$resid1))

anova(rain2001mean,rain2001linear, rain2001quad, rain2001cubic, rain2001cubicalt)

library(lattice)
w01=wireframe(X3~X1*X2,data=year2001fit,
              xlab = "Longitude", ylab = "Latitude", zlab="yhat", 
              main = "2001", scales = list(arrows = FALSE, distance=2),
              drape = TRUE,aspect = c(70/87, 0.41),
              colorkey = TRUE, screen = list(z = 30, x = -25, y = 10))
w01
c01= contourplot(X3~X1*X2,data=year2001fit,
                 cuts = 20, region = TRUE,
                 xlab = "Longitude",
                 ylab = "Latitude",
                 main = "2001")
c01 

attach(year2001fit)
X1=X2=seq(-1,1,.1)
model=function(X1,X2){13.393-7.824*X1 -138.389*X2-60.968 *X1^2-28.688*X2^2
  +28.689*X1^3+66.610*X2^3-64.201*X2^2*X1+19.485*X1^2*X2-46.510*X1*X2}
model(X1,X2)

z=outer(X1,X2,model)
contour(X1,X2,z,nlevels=20, main="2001")



wireframe(X3~X1*X2,data=year2001res,
          xlab = "Lon", ylab = "Lat", zlab="ehat", 
          main = "Precip.residuals year 2001", scales = list(arrows = FALSE, distance=2),
          drape = TRUE,aspect = c(70/87, 0.22),
          colorkey = TRUE, screen = list(z = 50, x = -55, y = 2))

## Kriging 2001---- 
coordinates(year2001res)=~X1 +X2
proj4string(year2001res) <- "+proj=utm +zone=37 +datum=WGS84 +units=m"
class(year2001res)
summary(year2001res)

# For working with spatial (and spatio-temporal) data, we use the gstat package, which includes functionality for kriging, among other many things.
library(sp)
library(gstat)
# packages for manipulation & visualization
suppressPackageStartupMessages({
  library(dplyr) # for "glimpse"
  library(ggplot2)
  library(scales) # for "comma"
  library(magrittr)
})


##a variogram could be fit as simply as the following code:
vgm <- variogram(X3~1, year2001res) # calculates sample variogram values 
plot(vgm)
library(gstat)
vgm.fit01<- fit.variogram(vgm, model=vgm(200, "Sph", 2.5,0)) # fit model
plot(vgm, vgm.fit01)

# Arrange the data for the ggplot2 plot
# add the semivariance values of v2 to v1
Fitted <- data.frame(dist = seq(0.12, max(vgm$dist), length = 1272))
Fitted$Spherical <- variogramLine(vgm.fit01, dist_vector = Fitted$dist)$gamma
#convert the dataframes to a long format
library(reshape2)
Empirical <- melt(vgm, id.vars = "dist", measure.vars = "gamma")
Modeled <- melt(Fitted, id.vars = "dist", measure.vars = "Spherical")
library(ggplot2)
#both variogram on the same plot with different colours
v01=ggplot(Modeled, aes(x = dist, y = value, colour = variable )) + 
  geom_line() +  
  xlab("Distance") +
  ylab("Semivariance") +   ggtitle("2001")+ 
  geom_point(data = Empirical)
v01

library(magrittr)
library(gstat)
library(sp)
rain.oK01<- krige(X3 ~ 1, year2001res, grd,model = vgm.fit01)
x <- krige.cv(X3 ~ 1, year2001res, grd,model = vgm.fit01, nfold=100)
#RMSE (root mean squared error)
sqrt(mean(x$residual^2)) 

library(ggplot2)
library(scales)
k01ps=rain.oK01 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKpred"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k01ps=k01ps+xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2001")
k01ps 

k01sd=rain.oK01 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.var)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKerror"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k01sd=k01sd+  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2001")
k01sd

# Some parameters for plotting
par(font.main=2, cex.main=1.5,cex.lab=0.4, cex.sub=1)
# Use the ColorBrewer library for color ramps
library(RColorBrewer)
precip.pal = colorRampPalette(c("yellow", "white", "blue"))

#### Countor map with values : Cokriging
ok01ps=spplot(rain.oK01,  zcol='var1.pred',col.regions=precip.pal, contour=TRUE, scales=list(draw=T),
              col='black', pretty=TRUE,main="2001", 
              labels = list(cex = 0.7), label.style = 'align', margin = TRUE,   width = 2, cex = 1.5,  xlab="Longitude",
              ylab="Latitude")
ok01ps

#Year 2002-------
year2002 <-subset(byyear, year==2002)
attach(year2002)
year2002$x1=lat.standard=(2*lat-(max(lat)+min(lat)))/(max(lat)-min(lat))
year2002$x2=lon.standard=(2*lon-(max(lon)+min(lon)))/(max(lon)-min(lon))

year2002$X2=year2002$x2^2
year2002$X1=year2002$x1^2
year2002$X13=year2002$x1^3
year2002$X23=year2002$x2^3
year2002$X2x1=year2002$X2*year2002$x1
year2002$X1x2=year2002$X1*year2002$x2

rain2002mean=lm(anom~1, year2002)
anova(rain2002mean)
summary(rain2002mean)

rain2002linear=lm(anom~x1+x2, year2002)
anova(rain2002linear)
summary(rain2002linear)

rain2002quad=lm(anom~x1*x2+X1+X2, year2002)
anova(rain2002quad)
summary(rain2002quad)

rain2002quadalt=lm(anom~x1*x2+X1+X2+alt, year2002)
anova(rain2002quadalt)
summary(rain2002quadalt)

rain2002cubic=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2, year2002)
anova(rain2002cubic)
summary(rain2002cubic)

rain2002cubicalt=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2+alt, year2002)
anova(rain2002cubicalt)
summary(rain2002cubicalt)


year2002$fitted1=fitted(rain2002cubic)
year2002$resid1=resid(rain2002cubic) #List of residuals
plot(density(resid(rain2002cubic))) #A density plot
qqnorm(resid(rain2002cubic)) # A quantile normal plot - good for checking normality
qqline(resid(rain2002cubic))
year2002fit=data.frame(cbind(year2002$lon, year2002$lat,year2002$fitted1))
year2002res=data.frame(cbind(year2002$lon, year2002$lat,year2002$resid1))

anova(rain2002mean,rain2002linear, rain2002quad, rain2002cubic, rain2002cubicalt)

library(lattice)
w02=wireframe(X3~X1*X2,data=year2002fit,
              xlab = "Longitude", ylab = "Latitude", zlab="yhat", 
              main = "2002", scales = list(arrows = FALSE, distance=2),
              drape = TRUE,aspect = c(70/87, 0.41),
              colorkey = TRUE, screen = list(z = 30, x = -25, y = 10))
w02
c02= contourplot(X3~X1*X2,data=year2002fit,
                 cuts = 20, region = TRUE,
                 xlab = "Longitude",
                 ylab = "Latitude",
                 main = "2002")
c02 

attach(year2002fit)
X1=X2=seq(-1,1,.1)
model=function(X1,X2){13.393-7.824*X1 -138.389*X2-60.968 *X1^2-28.688*X2^2
  +28.689*X1^3+66.610*X2^3-64.202*X2^2*X1+19.485*X1^2*X2-46.510*X1*X2}
model(X1,X2)

z=outer(X1,X2,model)
contour(X1,X2,z,nlevels=20, main="2002")



wireframe(X3~X1*X2,data=year2002res,
          xlab = "Lon", ylab = "Lat", zlab="ehat", 
          main = "Precip.residuals year 2002", scales = list(arrows = FALSE, distance=2),
          drape = TRUE,aspect = c(70/87, 0.22),
          colorkey = TRUE, screen = list(z = 50, x = -55, y = 2))

## Kriging 2002---- 
coordinates(year2002res)=~X1 +X2
proj4string(year2002res) <- "+proj=utm +zone=37 +datum=WGS84 +units=m"
class(year2002res)
summary(year2002res)

# For working with spatial (and spatio-temporal) data, we use the gstat package, which includes functionality for kriging, among other many things.
library(sp)
library(gstat)
# packages for manipulation & visualization
suppressPackageStartupMessages({
  library(dplyr) # for "glimpse"
  library(ggplot2)
  library(scales) # for "comma"
  library(magrittr)
})


##a variogram could be fit as simply as the following code:
vgm <- variogram(X3~1, year2002res) # calculates sample variogram values 
plot(vgm)
library(gstat)
vgm.fit02<- fit.variogram(vgm, model=vgm(180, "Sph", 2.45,0)) # fit model
plot(vgm, vgm.fit02)

# Arrange the data for the ggplot2 plot
# add the semivariance values of v2 to v1
Fitted <- data.frame(dist = seq(0.12, max(vgm$dist), length = 1272))
Fitted$Spherical <- variogramLine(vgm.fit02, dist_vector = Fitted$dist)$gamma
#convert the dataframes to a long format
library(reshape2)
Empirical <- melt(vgm, id.vars = "dist", measure.vars = "gamma")
Modeled <- melt(Fitted, id.vars = "dist", measure.vars = "Spherical")
library(ggplot2)
#both variogram on the same plot with different colours
v02=ggplot(Modeled, aes(x = dist, y = value, colour = variable )) + 
  geom_line() +  
  xlab("Distance") +
  ylab("Semivariance") +   ggtitle("2002")+ 
  geom_point(data = Empirical)
v02

library(magrittr)
library(gstat)
library(sp)
rain.oK02<- krige(X3 ~ 1, year2002res, grd,model = vgm.fit02)
x <- krige.cv(X3 ~ 1, year2002res, grd,model = vgm.fit02, nfold=102)
#RMSE (root mean squared error)
sqrt(mean(x$residual^2)) 

library(ggplot2)
library(scales)
k02ps=rain.oK02 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKpred"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k02ps=k02ps+xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2002")
k02ps 

k02sd=rain.oK02 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.var)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKerror"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k02sd=k02sd+  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2002")
k02sd

# Some parameters for plotting
par(font.main=2, cex.main=1.5,cex.lab=0.4, cex.sub=1)
# Use the ColorBrewer library for color ramps
library(RColorBrewer)
precip.pal = colorRampPalette(c("yellow", "white", "blue"))

#### Countor map with values : Cokriging
ok02ps=spplot(rain.oK02,  zcol='var1.pred',col.regions=precip.pal, contour=TRUE, scales=list(draw=T),
              col='black', pretty=TRUE,main="2002", 
              labels = list(cex = 0.7), label.style = 'align', margin = TRUE,   width = 2, cex = 1.5,  xlab="Longitude",
              ylab="Latitude")
ok02ps


#Year 2003-------
year2003 <-subset(byyear, year==2003)
attach(year2003)
year2003$x1=lat.standard=(2*lat-(max(lat)+min(lat)))/(max(lat)-min(lat))
year2003$x2=lon.standard=(2*lon-(max(lon)+min(lon)))/(max(lon)-min(lon))

year2003$X2=year2003$x2^2
year2003$X1=year2003$x1^2
year2003$X13=year2003$x1^3
year2003$X23=year2003$x2^3
year2003$X2x1=year2003$X2*year2003$x1
year2003$X1x2=year2003$X1*year2003$x2

rain2003mean=lm(anom~1, year2003)
anova(rain2003mean)
summary(rain2003mean)

rain2003linear=lm(anom~x1+x2, year2003)
anova(rain2003linear)
summary(rain2003linear)

rain2003quad=lm(anom~x1*x2+X1+X2, year2003)
anova(rain2003quad)
summary(rain2003quad)

rain2003quadalt=lm(anom~x1*x2+X1+X2+alt, year2003)
anova(rain2003quadalt)
summary(rain2003quadalt)

rain2003cubic=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2, year2003)
anova(rain2003cubic)
summary(rain2003cubic)

rain2003cubicalt=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2+alt, year2003)
anova(rain2003cubicalt)
summary(rain2003cubicalt)


year2003$fitted1=fitted(rain2003cubic)
year2003$resid1=resid(rain2003cubic) #List of residuals
plot(density(resid(rain2003cubic))) #A density plot
qqnorm(resid(rain2003cubic)) # A quantile normal plot - good for checking normality
qqline(resid(rain2003cubic))
year2003fit=data.frame(cbind(year2003$lon, year2003$lat,year2003$fitted1))
year2003res=data.frame(cbind(year2003$lon, year2003$lat,year2003$resid1))

anova(rain2003mean,rain2003linear, rain2003quad, rain2003cubic, rain2003cubicalt)

library(lattice)
w03=wireframe(X3~X1*X2,data=year2003fit,
              xlab = "Longitude", ylab = "Latitude", zlab="yhat", 
              main = "2003", scales = list(arrows = FALSE, distance=2),
              drape = TRUE,aspect = c(70/87, 0.41),
              colorkey = TRUE, screen = list(z = 30, x = -25, y = 10))
w03
c03= contourplot(X3~X1*X2,data=year2003fit,
                 cuts = 20, region = TRUE,
                 xlab = "Longitude",
                 ylab = "Latitude",
                 main = "2003")
c03 

attach(year2003fit)
X1=X2=seq(-1,1,.1)
model=function(X1,X2){13.393-7.824*X1 -138.389*X2-60.968 *X1^2-28.688*X2^2
  +28.689*X1^3+66.610*X2^3-64.203*X2^2*X1+19.485*X1^2*X2-46.510*X1*X2}
model(X1,X2)

z=outer(X1,X2,model)
contour(X1,X2,z,nlevels=20, main="2003")



wireframe(X3~X1*X2,data=year2003res,
          xlab = "Lon", ylab = "Lat", zlab="ehat", 
          main = "Precip.residuals year 2003", scales = list(arrows = FALSE, distance=2),
          drape = TRUE,aspect = c(70/87, 0.22),
          colorkey = TRUE, screen = list(z = 50, x = -55, y = 2))

## Kriging 2003---- 
coordinates(year2003res)=~X1 +X2
proj4string(year2003res) <- "+proj=utm +zone=37 +datum=WGS84 +units=m"
class(year2003res)
summary(year2003res)

# For working with spatial (and spatio-temporal) data, we use the gstat package, which includes functionality for kriging, among other many things.
library(sp)
library(gstat)
# packages for manipulation & visualization
suppressPackageStartupMessages({
  library(dplyr) # for "glimpse"
  library(ggplot2)
  library(scales) # for "comma"
  library(magrittr)
})


##a variogram could be fit as simply as the following code:
vgm <- variogram(X3~1, year2003res) # calculates sample variogram values 
plot(vgm)
library(gstat)
vgm.fit03<- fit.variogram(vgm, model=vgm(155, "Sph", 2.1,0)) # fit model
plot(vgm, vgm.fit03)

# Arrange the data for the ggplot2 plot
# add the semivariance values of v2 to v1
Fitted <- data.frame(dist = seq(0.12, max(vgm$dist), length = 1272))
Fitted$Spherical <- variogramLine(vgm.fit03, dist_vector = Fitted$dist)$gamma
#convert the dataframes to a long format
library(reshape2)
Empirical <- melt(vgm, id.vars = "dist", measure.vars = "gamma")
Modeled <- melt(Fitted, id.vars = "dist", measure.vars = "Spherical")
library(ggplot2)
#both variogram on the same plot with different colours
v03=ggplot(Modeled, aes(x = dist, y = value, colour = variable )) + 
  geom_line() +  
  xlab("Distance") +
  ylab("Semivariance") +   ggtitle("2003")+ 
  geom_point(data = Empirical)
v03

library(magrittr)
library(gstat)
library(sp)
rain.oK03<- krige(X3 ~ 1, year2003res, grd,model = vgm.fit03)
x <- krige.cv(X3 ~ 1, year2003res, grd,model = vgm.fit03, nfold=103)
#RMSE (root mean squared error)
sqrt(mean(x$residual^2)) 

library(ggplot2)
library(scales)
k03ps=rain.oK03 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKpred"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k03ps=k03ps+xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2003")
k03ps 

k03sd=rain.oK03 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.var)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKerror"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k03sd=k03sd+  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2003")
k03sd

# Some parameters for plotting
par(font.main=2, cex.main=1.5,cex.lab=0.4, cex.sub=1)
# Use the ColorBrewer library for color ramps
library(RColorBrewer)
precip.pal = colorRampPalette(c("yellow", "white", "blue"))

#### Countor map with values : Cokriging
ok03ps=spplot(rain.oK03,  zcol='var1.pred',col.regions=precip.pal, contour=TRUE, scales=list(draw=T),
              col='black', pretty=TRUE,main="2003", 
              labels = list(cex = 0.7), label.style = 'align', margin = TRUE,   width = 2, cex = 1.5,  xlab="Longitude",
              ylab="Latitude")
ok03ps
#Year 2004-------
year2004 <-subset(byyear, year==2004)
attach(year2004)
year2004$x1=lat.standard=(2*lat-(max(lat)+min(lat)))/(max(lat)-min(lat))
year2004$x2=lon.standard=(2*lon-(max(lon)+min(lon)))/(max(lon)-min(lon))

year2004$X2=year2004$x2^2
year2004$X1=year2004$x1^2
year2004$X13=year2004$x1^3
year2004$X23=year2004$x2^3
year2004$X2x1=year2004$X2*year2004$x1
year2004$X1x2=year2004$X1*year2004$x2

rain2004mean=lm(anom~1, year2004)
anova(rain2004mean)
summary(rain2004mean)

rain2004linear=lm(anom~x1+x2, year2004)
anova(rain2004linear)
summary(rain2004linear)

rain2004quad=lm(anom~x1*x2+X1+X2, year2004)
anova(rain2004quad)
summary(rain2004quad)

rain2004quadalt=lm(anom~x1*x2+X1+X2+alt, year2004)
anova(rain2004quadalt)
summary(rain2004quadalt)

rain2004cubic=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2, year2004)
anova(rain2004cubic)
summary(rain2004cubic)

rain2004cubicalt=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2+alt, year2004)
anova(rain2004cubicalt)
summary(rain2004cubicalt)


year2004$fitted1=fitted(rain2004cubic)
year2004$resid1=resid(rain2004cubic) #List of residuals
plot(density(resid(rain2004cubic))) #A density plot
qqnorm(resid(rain2004cubic)) # A quantile normal plot - good for checking normality
qqline(resid(rain2004cubic))
year2004fit=data.frame(cbind(year2004$lon, year2004$lat,year2004$fitted1))
year2004res=data.frame(cbind(year2004$lon, year2004$lat,year2004$resid1))

anova(rain2004mean,rain2004linear, rain2004quad, rain2004cubic, rain2004cubicalt)

library(lattice)
w04=wireframe(X3~X1*X2,data=year2004fit,
              xlab = "Longitude", ylab = "Latitude", zlab="yhat", 
              main = "2004", scales = list(arrows = FALSE, distance=2),
              drape = TRUE,aspect = c(70/87, 0.41),
              colorkey = TRUE, screen = list(z = 30, x = -25, y = 10))
w04
c04= contourplot(X3~X1*X2,data=year2004fit,
                 cuts = 20, region = TRUE,
                 xlab = "Longitude",
                 ylab = "Latitude",
                 main = "2004")
c04 

attach(year2004fit)
X1=X2=seq(-1,1,.1)
model=function(X1,X2){13.393-7.824*X1 -138.389*X2-60.968 *X1^2-28.688*X2^2
  +28.689*X1^3+66.610*X2^3-64.204*X2^2*X1+19.485*X1^2*X2-46.510*X1*X2}
model(X1,X2)

z=outer(X1,X2,model)
contour(X1,X2,z,nlevels=20, main="2004")



wireframe(X3~X1*X2,data=year2004res,
          xlab = "Lon", ylab = "Lat", zlab="ehat", 
          main = "Precip.residuals year 2004", scales = list(arrows = FALSE, distance=2),
          drape = TRUE,aspect = c(70/87, 0.22),
          colorkey = TRUE, screen = list(z = 50, x = -55, y = 2))

## Kriging 2004---- 
coordinates(year2004res)=~X1 +X2
proj4string(year2004res) <- "+proj=utm +zone=37 +datum=WGS84 +units=m"
class(year2004res)
summary(year2004res)

# For working with spatial (and spatio-temporal) data, we use the gstat package, which includes functionality for kriging, among other many things.
library(sp)
library(gstat)
# packages for manipulation & visualization
suppressPackageStartupMessages({
  library(dplyr) # for "glimpse"
  library(ggplot2)
  library(scales) # for "comma"
  library(magrittr)
})


##a variogram could be fit as simply as the following code:
vgm <- variogram(X3~1, year2004res) # calculates sample variogram values 
plot(vgm)
library(gstat)
vgm.fit04<- fit.variogram(vgm, model=vgm(220, "Sph", 2.22,0)) # fit model
plot(vgm, vgm.fit04)

# Arrange the data for the ggplot2 plot
# add the semivariance values of v2 to v1
Fitted <- data.frame(dist = seq(0.12, max(vgm$dist), length = 1272))
Fitted$Spherical <- variogramLine(vgm.fit04, dist_vector = Fitted$dist)$gamma
#convert the dataframes to a long format
library(reshape2)
Empirical <- melt(vgm, id.vars = "dist", measure.vars = "gamma")
Modeled <- melt(Fitted, id.vars = "dist", measure.vars = "Spherical")
library(ggplot2)
#both variogram on the same plot with different colours
v04=ggplot(Modeled, aes(x = dist, y = value, colour = variable )) + 
  geom_line() +  
  xlab("Distance") +
  ylab("Semivariance") +   ggtitle("2004")+ 
  geom_point(data = Empirical)
v04

library(magrittr)
library(gstat)
library(sp)
rain.oK04<- krige(X3 ~ 1, year2004res, grd,model = vgm.fit04)
x <- krige.cv(X3 ~ 1, year2004res, grd,model = vgm.fit04, nfold=100)
#RMSE (root mean squared error)
sqrt(mean(x$residual^2)) 

library(ggplot2)
library(scales)
k04ps=rain.oK04 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKpred"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k04ps=k04ps+xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2004")
k04ps 

k04sd=rain.oK04 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.var)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKerror"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k04sd=k04sd+  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2004")
k04sd

# Some parameters for plotting
par(font.main=2, cex.main=1.5,cex.lab=0.4, cex.sub=1)
# Use the ColorBrewer library for color ramps
library(RColorBrewer)
precip.pal = colorRampPalette(c("yellow", "white", "blue"))

#### Countor map with values : Cokriging
ok04ps=spplot(rain.oK04,  zcol='var1.pred',col.regions=precip.pal, contour=TRUE, scales=list(draw=T),
              col='black', pretty=TRUE,main="2004", 
              labels = list(cex = 0.7), label.style = 'align', margin = TRUE,   width = 2, cex = 1.5,  xlab="Longitude",
              ylab="Latitude")
ok04ps
#Year 2005-------
year2005 <-subset(byyear, year==2005)
attach(year2005)
year2005$x1=lat.standard=(2*lat-(max(lat)+min(lat)))/(max(lat)-min(lat))
year2005$x2=lon.standard=(2*lon-(max(lon)+min(lon)))/(max(lon)-min(lon))

year2005$X2=year2005$x2^2
year2005$X1=year2005$x1^2
year2005$X13=year2005$x1^3
year2005$X23=year2005$x2^3
year2005$X2x1=year2005$X2*year2005$x1
year2005$X1x2=year2005$X1*year2005$x2

rain2005mean=lm(anom~1, year2005)
anova(rain2005mean)
summary(rain2005mean)

rain2005linear=lm(anom~x1+x2, year2005)
anova(rain2005linear)
summary(rain2005linear)

rain2005quad=lm(anom~x1*x2+X1+X2, year2005)
anova(rain2005quad)
summary(rain2005quad)

rain2005quadalt=lm(anom~x1*x2+X1+X2+alt, year2005)
anova(rain2005quadalt)
summary(rain2005quadalt)

rain2005cubic=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2, year2005)
anova(rain2005cubic)
summary(rain2005cubic)

rain2005cubicalt=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2+alt, year2005)
anova(rain2005cubicalt)
summary(rain2005cubicalt)


year2005$fitted1=fitted(rain2005cubic)
year2005$resid1=resid(rain2005cubic) #List of residuals
plot(density(resid(rain2005cubic))) #A density plot
qqnorm(resid(rain2005cubic)) # A quantile normal plot - good for checking normality
qqline(resid(rain2005cubic))
year2005fit=data.frame(cbind(year2005$lon, year2005$lat,year2005$fitted1))
year2005res=data.frame(cbind(year2005$lon, year2005$lat,year2005$resid1))

anova(rain2005mean,rain2005linear, rain2005quad, rain2005cubic, rain2005cubicalt)

library(lattice)
w05=wireframe(X3~X1*X2,data=year2005fit,
              xlab = "Longitude", ylab = "Latitude", zlab="yhat", 
              main = "2005", scales = list(arrows = FALSE, distance=2),
              drape = TRUE,aspect = c(70/87, 0.41),
              colorkey = TRUE, screen = list(z = 30, x = -25, y = 10))
w05
c05= contourplot(X3~X1*X2,data=year2005fit,
                 cuts = 20, region = TRUE,
                 xlab = "Longitude",
                 ylab = "Latitude",
                 main = "2005")
c05 

attach(year2005fit)
X1=X2=seq(-1,1,.1)
model=function(X1,X2){13.393-7.824*X1 -138.389*X2-60.968 *X1^2-28.688*X2^2
  +28.689*X1^3+66.610*X2^3-64.205*X2^2*X1+19.485*X1^2*X2-46.510*X1*X2}
model(X1,X2)

z=outer(X1,X2,model)
contour(X1,X2,z,nlevels=20, main="2005")



wireframe(X3~X1*X2,data=year2005res,
          xlab = "Lon", ylab = "Lat", zlab="ehat", 
          main = "Precip.residuals year 2005", scales = list(arrows = FALSE, distance=2),
          drape = TRUE,aspect = c(70/87, 0.22),
          colorkey = TRUE, screen = list(z = 50, x = -55, y = 2))

## Kriging 2005---- 
coordinates(year2005res)=~X1 +X2
proj4string(year2005res) <- "+proj=utm +zone=37 +datum=WGS84 +units=m"
class(year2005res)
summary(year2005res)

# For working with spatial (and spatio-temporal) data, we use the gstat package, which includes functionality for kriging, among other many things.
library(sp)
library(gstat)
# packages for manipulation & visualization
suppressPackageStartupMessages({
  library(dplyr) # for "glimpse"
  library(ggplot2)
  library(scales) # for "comma"
  library(magrittr)
})


##a variogram could be fit as simply as the following code:
vgm <- variogram(X3~1, year2005res) # calculates sample variogram values 
plot(vgm)
library(gstat)
vgm.fit05<- fit.variogram(vgm, model=vgm(300, "Sph", 2,30)) # fit model
plot(vgm, vgm.fit05)

# Arrange the data for the ggplot2 plot
# add the semivariance values of v2 to v1
Fitted <- data.frame(dist = seq(0.12, max(vgm$dist), length = 1272))
Fitted$Spherical <- variogramLine(vgm.fit05, dist_vector = Fitted$dist)$gamma
#convert the dataframes to a long format
library(reshape2)
Empirical <- melt(vgm, id.vars = "dist", measure.vars = "gamma")
Modeled <- melt(Fitted, id.vars = "dist", measure.vars = "Spherical")
library(ggplot2)
#both variogram on the same plot with different colours
v05=ggplot(Modeled, aes(x = dist, y = value, colour = variable )) + 
  geom_line() +  
  xlab("Distance") +
  ylab("Semivariance") +   ggtitle("2005")+ 
  geom_point(data = Empirical)
v05

library(magrittr)
library(gstat)
library(sp)
rain.oK05<- krige(X3 ~ 1, year2005res, grd,model = vgm.fit05)
x <- krige.cv(X3 ~ 1, year2005res, grd,model = vgm.fit05, nfold=100)
#RMSE (root mean squared error)
sqrt(mean(x$residual^2)) 

library(ggplot2)
library(scales)
k05ps=rain.oK05 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKpred"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k05ps=k05ps+xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2005")
k05ps 

k05sd=rain.oK05 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.var)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKerror"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k05sd=k05sd+  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2005")
k05sd

# Some parameters for plotting
par(font.main=2, cex.main=1.5,cex.lab=0.4, cex.sub=1)
# Use the ColorBrewer library for color ramps
library(RColorBrewer)
precip.pal = colorRampPalette(c("yellow", "white", "blue"))

#### Countor map with values : Cokriging
ok05ps=spplot(rain.oK05,  zcol='var1.pred',col.regions=precip.pal, contour=TRUE, scales=list(draw=T),
              col='black', pretty=TRUE,main="2005", 
              labels = list(cex = 0.7), label.style = 'align', margin = TRUE,   width = 2, cex = 1.5,  xlab="Longitude",
              ylab="Latitude")
ok05ps
#Year 2006-------
year2006 <-subset(byyear, year==2006)
attach(year2006)
year2006$x1=lat.standard=(2*lat-(max(lat)+min(lat)))/(max(lat)-min(lat))
year2006$x2=lon.standard=(2*lon-(max(lon)+min(lon)))/(max(lon)-min(lon))

year2006$X2=year2006$x2^2
year2006$X1=year2006$x1^2
year2006$X13=year2006$x1^3
year2006$X23=year2006$x2^3
year2006$X2x1=year2006$X2*year2006$x1
year2006$X1x2=year2006$X1*year2006$x2

rain2006mean=lm(anom~1, year2006)
anova(rain2006mean)
summary(rain2006mean)

rain2006linear=lm(anom~x1+x2, year2006)
anova(rain2006linear)
summary(rain2006linear)

rain2006quad=lm(anom~x1*x2+X1+X2, year2006)
anova(rain2006quad)
summary(rain2006quad)

rain2006quadalt=lm(anom~x1*x2+X1+X2+alt, year2006)
anova(rain2006quadalt)
summary(rain2006quadalt)

rain2006cubic=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2, year2006)
anova(rain2006cubic)
summary(rain2006cubic)

rain2006cubicalt=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2+alt, year2006)
anova(rain2006cubicalt)
summary(rain2006cubicalt)


year2006$fitted1=fitted(rain2006cubic)
year2006$resid1=resid(rain2006cubic) #List of residuals
plot(density(resid(rain2006cubic))) #A density plot
qqnorm(resid(rain2006cubic)) # A quantile normal plot - good for checking normality
qqline(resid(rain2006cubic))
year2006fit=data.frame(cbind(year2006$lon, year2006$lat,year2006$fitted1))
year2006res=data.frame(cbind(year2006$lon, year2006$lat,year2006$resid1))

anova(rain2006mean,rain2006linear, rain2006quad, rain2006cubic, rain2006cubicalt)

library(lattice)
w06=wireframe(X3~X1*X2,data=year2006fit,
              xlab = "Longitude", ylab = "Latitude", zlab="yhat", 
              main = "2006", scales = list(arrows = FALSE, distance=2),
              drape = TRUE,aspect = c(70/87, 0.41),
              colorkey = TRUE, screen = list(z = 30, x = -25, y = 10))
w06
c06= contourplot(X3~X1*X2,data=year2006fit,
                 cuts = 20, region = TRUE,
                 xlab = "Longitude",
                 ylab = "Latitude",
                 main = "2006")
c06 

attach(year2006fit)
X1=X2=seq(-1,1,.1)
model=function(X1,X2){13.393-7.824*X1 -138.389*X2-60.968 *X1^2-28.688*X2^2
  +28.689*X1^3+66.610*X2^3-64.206*X2^2*X1+19.485*X1^2*X2-46.510*X1*X2}
model(X1,X2)

z=outer(X1,X2,model)
contour(X1,X2,z,nlevels=20, main="2006")



wireframe(X3~X1*X2,data=year2006res,
          xlab = "Lon", ylab = "Lat", zlab="ehat", 
          main = "Precip.residuals year 2006", scales = list(arrows = FALSE, distance=2),
          drape = TRUE,aspect = c(70/87, 0.22),
          colorkey = TRUE, screen = list(z = 50, x = -55, y = 2))

## Kriging 2006---- 
coordinates(year2006res)=~X1 +X2
proj4string(year2006res) <- "+proj=utm +zone=37 +datum=WGS84 +units=m"
class(year2006res)
summary(year2006res)

# For working with spatial (and spatio-temporal) data, we use the gstat package, which includes functionality for kriging, among other many things.
library(sp)
library(gstat)
# packages for manipulation & visualization
suppressPackageStartupMessages({
  library(dplyr) # for "glimpse"
  library(ggplot2)
  library(scales) # for "comma"
  library(magrittr)
})


##a variogram could be fit as simply as the following code:
vgm <- variogram(X3~1, year2006res) # calculates sample variogram values 
plot(vgm)
library(gstat)
vgm.fit06<- fit.variogram(vgm, model=vgm(370, "Sph", 2.3,70)) # fit model
plot(vgm, vgm.fit06)

# Arrange the data for the ggplot2 plot
# add the semivariance values of v2 to v1
Fitted <- data.frame(dist = seq(0.12, max(vgm$dist), length = 1272))
Fitted$Spherical <- variogramLine(vgm.fit06, dist_vector = Fitted$dist)$gamma
#convert the dataframes to a long format
library(reshape2)
Empirical <- melt(vgm, id.vars = "dist", measure.vars = "gamma")
Modeled <- melt(Fitted, id.vars = "dist", measure.vars = "Spherical")
library(ggplot2)
#both variogram on the same plot with different colours
v06=ggplot(Modeled, aes(x = dist, y = value, colour = variable )) + 
  geom_line() +  
  xlab("Distance") +
  ylab("Semivariance") +   ggtitle("2006")+ 
  geom_point(data = Empirical)
v06

library(magrittr)
library(gstat)
library(sp)
rain.oK06<- krige(X3 ~ 1, year2006res, grd,model = vgm.fit06)
x <- krige.cv(X3 ~ 1, year2006res, grd,model = vgm.fit06, nfold=106)
#RMSE (root mean squared error)
sqrt(mean(x$residual^2)) 

library(ggplot2)
library(scales)
k06ps=rain.oK06 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKpred"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k06ps=k06ps+xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2006")
k06ps 

k06sd=rain.oK06 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.var)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKerror"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k06sd=k06sd+  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2006")
k06sd

# Some parameters for plotting
par(font.main=2, cex.main=1.5,cex.lab=0.4, cex.sub=1)
# Use the ColorBrewer library for color ramps
library(RColorBrewer)
precip.pal = colorRampPalette(c("yellow", "white", "blue"))

#### Countor map with values : Cokriging
ok06ps=spplot(rain.oK06,  zcol='var1.pred',col.regions=precip.pal, contour=TRUE, scales=list(draw=T),
              col='black', pretty=TRUE,main="2006", 
              labels = list(cex = 0.7), label.style = 'align', margin = TRUE,   width = 2, cex = 1.5,  xlab="Longitude",
              ylab="Latitude")
ok06ps
#Year 2007-------
year2007 <-subset(byyear, year==2007)
attach(year2007)
year2007$x1=lat.standard=(2*lat-(max(lat)+min(lat)))/(max(lat)-min(lat))
year2007$x2=lon.standard=(2*lon-(max(lon)+min(lon)))/(max(lon)-min(lon))

year2007$X2=year2007$x2^2
year2007$X1=year2007$x1^2
year2007$X13=year2007$x1^3
year2007$X23=year2007$x2^3
year2007$X2x1=year2007$X2*year2007$x1
year2007$X1x2=year2007$X1*year2007$x2

rain2007mean=lm(anom~1, year2007)
anova(rain2007mean)
summary(rain2007mean)

rain2007linear=lm(anom~x1+x2, year2007)
anova(rain2007linear)
summary(rain2007linear)

rain2007quad=lm(anom~x1*x2+X1+X2, year2007)
anova(rain2007quad)
summary(rain2007quad)

rain2007quadalt=lm(anom~x1*x2+X1+X2+alt, year2007)
anova(rain2007quadalt)
summary(rain2007quadalt)

rain2007cubic=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2, year2007)
anova(rain2007cubic)
summary(rain2007cubic)

rain2007cubicalt=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2+alt, year2007)
anova(rain2007cubicalt)
summary(rain2007cubicalt)


year2007$fitted1=fitted(rain2007cubic)
year2007$resid1=resid(rain2007cubic) #List of residuals
plot(density(resid(rain2007cubic))) #A density plot
qqnorm(resid(rain2007cubic)) # A quantile normal plot - good for checking normality
qqline(resid(rain2007cubic))
year2007fit=data.frame(cbind(year2007$lon, year2007$lat,year2007$fitted1))
year2007res=data.frame(cbind(year2007$lon, year2007$lat,year2007$resid1))

anova(rain2007mean,rain2007linear, rain2007quad, rain2007cubic, rain2007cubicalt)

library(lattice)
w07=wireframe(X3~X1*X2,data=year2007fit,
              xlab = "Longitude", ylab = "Latitude", zlab="yhat", 
              main = "2007", scales = list(arrows = FALSE, distance=2),
              drape = TRUE,aspect = c(70/87, 0.41),
              colorkey = TRUE, screen = list(z = 30, x = -25, y = 10))
w07
c07= contourplot(X3~X1*X2,data=year2007fit,
                 cuts = 20, region = TRUE,
                 xlab = "Longitude",
                 ylab = "Latitude",
                 main = "2007")
c07 

attach(year2007fit)
X1=X2=seq(-1,1,.1)
model=function(X1,X2){13.393-7.824*X1 -138.389*X2-60.968 *X1^2-28.688*X2^2
  +28.689*X1^3+66.610*X2^3-64.207*X2^2*X1+19.485*X1^2*X2-46.510*X1*X2}
model(X1,X2)

z=outer(X1,X2,model)
contour(X1,X2,z,nlevels=20, main="2007")



wireframe(X3~X1*X2,data=year2007res,
          xlab = "Lon", ylab = "Lat", zlab="ehat", 
          main = "Precip.residuals year 2007", scales = list(arrows = FALSE, distance=2),
          drape = TRUE,aspect = c(70/87, 0.22),
          colorkey = TRUE, screen = list(z = 50, x = -55, y = 2))

## Kriging 2007---- 
coordinates(year2007res)=~X1 +X2
proj4string(year2007res) <- "+proj=utm +zone=37 +datum=WGS84 +units=m"
class(year2007res)
summary(year2007res)

# For working with spatial (and spatio-temporal) data, we use the gstat package, which includes functionality for kriging, among other many things.
library(sp)
library(gstat)
# packages for manipulation & visualization
suppressPackageStartupMessages({
  library(dplyr) # for "glimpse"
  library(ggplot2)
  library(scales) # for "comma"
  library(magrittr)
})


##a variogram could be fit as simply as the following code:
vgm <- variogram(X3~1, year2007res) # calculates sample variogram values 
plot(vgm)
library(gstat)
vgm.fit07<- fit.variogram(vgm, model=vgm(235, "Sph", 2,60)) # fit model
plot(vgm, vgm.fit07)

# Arrange the data for the ggplot2 plot
# add the semivariance values of v2 to v1
Fitted <- data.frame(dist = seq(0.12, max(vgm$dist), length = 1272))
Fitted$Spherical <- variogramLine(vgm.fit07, dist_vector = Fitted$dist)$gamma
#convert the dataframes to a long format
library(reshape2)
Empirical <- melt(vgm, id.vars = "dist", measure.vars = "gamma")
Modeled <- melt(Fitted, id.vars = "dist", measure.vars = "Spherical")
library(ggplot2)
#both variogram on the same plot with different colours
v07=ggplot(Modeled, aes(x = dist, y = value, colour = variable )) + 
  geom_line() +  
  xlab("Distance") +
  ylab("Semivariance") +   ggtitle("2007")+ 
  geom_point(data = Empirical)
v07

library(magrittr)
library(gstat)
library(sp)
rain.oK07<- krige(X3 ~ 1, year2007res, grd,model = vgm.fit07)
x <- krige.cv(X3 ~ 1, year2007res, grd,model = vgm.fit07, nfold=107)
#RMSE (root mean squared error)
sqrt(mean(x$residual^2)) 

library(ggplot2)
library(scales)
k07ps=rain.oK07 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKpred"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k07ps=k07ps+xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2007")
k07ps 

k07sd=rain.oK07 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.var)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKerror"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k07sd=k07sd+  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2007")
k07sd

# Some parameters for plotting
par(font.main=2, cex.main=1.5,cex.lab=0.4, cex.sub=1)
# Use the ColorBrewer library for color ramps
library(RColorBrewer)
precip.pal = colorRampPalette(c("yellow", "white", "blue"))

#### Countor map with values : Cokriging
ok07ps=spplot(rain.oK07,  zcol='var1.pred',col.regions=precip.pal, contour=TRUE, scales=list(draw=T),
              col='black', pretty=TRUE,main="2007", 
              labels = list(cex = 0.7), label.style = 'align', margin = TRUE,   width = 2, cex = 1.5,  xlab="Longitude",
              ylab="Latitude")
ok07ps
#Year 2008-------
year2008 <-subset(byyear, year==2008)
attach(year2008)
year2008$x1=lat.standard=(2*lat-(max(lat)+min(lat)))/(max(lat)-min(lat))
year2008$x2=lon.standard=(2*lon-(max(lon)+min(lon)))/(max(lon)-min(lon))

year2008$X2=year2008$x2^2
year2008$X1=year2008$x1^2
year2008$X13=year2008$x1^3
year2008$X23=year2008$x2^3
year2008$X2x1=year2008$X2*year2008$x1
year2008$X1x2=year2008$X1*year2008$x2

rain2008mean=lm(anom~1, year2008)
anova(rain2008mean)
summary(rain2008mean)

rain2008linear=lm(anom~x1+x2, year2008)
anova(rain2008linear)
summary(rain2008linear)

rain2008quad=lm(anom~x1*x2+X1+X2, year2008)
anova(rain2008quad)
summary(rain2008quad)

rain2008quadalt=lm(anom~x1*x2+X1+X2+alt, year2008)
anova(rain2008quadalt)
summary(rain2008quadalt)

rain2008cubic=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2, year2008)
anova(rain2008cubic)
summary(rain2008cubic)

rain2008cubicalt=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2+alt, year2008)
anova(rain2008cubicalt)
summary(rain2008cubicalt)


year2008$fitted1=fitted(rain2008cubic)
year2008$resid1=resid(rain2008cubic) #List of residuals
plot(density(resid(rain2008cubic))) #A density plot
qqnorm(resid(rain2008cubic)) # A quantile normal plot - good for checking normality
qqline(resid(rain2008cubic))
year2008fit=data.frame(cbind(year2008$lon, year2008$lat,year2008$fitted1))
year2008res=data.frame(cbind(year2008$lon, year2008$lat,year2008$resid1))

anova(rain2008mean,rain2008linear, rain2008quad, rain2008cubic, rain2008cubicalt)

library(lattice)
w08=wireframe(X3~X1*X2,data=year2008fit,
              xlab = "Longitude", ylab = "Latitude", zlab="yhat", 
              main = "2008", scales = list(arrows = FALSE, distance=2),
              drape = TRUE,aspect = c(70/87, 0.41),
              colorkey = TRUE, screen = list(z = 30, x = -25, y = 10))
w08
c08= contourplot(X3~X1*X2,data=year2008fit,
                 cuts = 20, region = TRUE,
                 xlab = "Longitude",
                 ylab = "Latitude",
                 main = "2008")
c08 

attach(year2008fit)
X1=X2=seq(-1,1,.1)
model=function(X1,X2){13.393-7.824*X1 -138.389*X2-60.968 *X1^2-28.688*X2^2
  +28.689*X1^3+66.610*X2^3-64.208*X2^2*X1+19.485*X1^2*X2-46.510*X1*X2}
model(X1,X2)

z=outer(X1,X2,model)
contour(X1,X2,z,nlevels=20, main="2008")



wireframe(X3~X1*X2,data=year2008res,
          xlab = "Lon", ylab = "Lat", zlab="ehat", 
          main = "Precip.residuals year 2008", scales = list(arrows = FALSE, distance=2),
          drape = TRUE,aspect = c(70/87, 0.22),
          colorkey = TRUE, screen = list(z = 50, x = -55, y = 2))

## Kriging 2008---- 
coordinates(year2008res)=~X1 +X2
proj4string(year2008res) <- "+proj=utm +zone=37 +datum=WGS84 +units=m"
class(year2008res)
summary(year2008res)

# For working with spatial (and spatio-temporal) data, we use the gstat package, which includes functionality for kriging, among other many things.
library(sp)
library(gstat)
# packages for manipulation & visualization
suppressPackageStartupMessages({
  library(dplyr) # for "glimpse"
  library(ggplot2)
  library(scales) # for "comma"
  library(magrittr)
})


##a variogram could be fit as simply as the following code:
vgm <- variogram(X3~1, year2008res) # calculates sample variogram values 
plot(vgm)
library(gstat)
vgm.fit08<- fit.variogram(vgm, model=vgm(325, "Sph",2,60)) # fit model
plot(vgm, vgm.fit08)

# Arrange the data for the ggplot2 plot
# add the semivariance values of v2 to v1
Fitted <- data.frame(dist = seq(0.12, max(vgm$dist), length = 1272))
Fitted$Spherical <- variogramLine(vgm.fit08, dist_vector = Fitted$dist)$gamma
#convert the dataframes to a long format
library(reshape2)
Empirical <- melt(vgm, id.vars = "dist", measure.vars = "gamma")
Modeled <- melt(Fitted, id.vars = "dist", measure.vars = "Spherical")
library(ggplot2)
#both variogram on the same plot with different colours
v08=ggplot(Modeled, aes(x = dist, y = value, colour = variable )) + 
  geom_line() +  
  xlab("Distance") +
  ylab("Semivariance") +   ggtitle("2008")+ 
  geom_point(data = Empirical)
v08

library(magrittr)
library(gstat)
library(sp)
rain.oK08<- krige(X3 ~ 1, year2008res, grd,model = vgm.fit08)
x <- krige.cv(X3 ~ 1, year2008res, grd,model = vgm.fit08, nfold=100)
#RMSE (root mean squared error)
sqrt(mean(x$residual^2)) 

library(ggplot2)
library(scales)
k08ps=rain.oK08 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKpred"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k08ps=k08ps+xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2008")
k08ps 

k08sd=rain.oK08 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.var)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKerror"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k08sd=k08sd+  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2008")
k08sd

# Some parameters for plotting
par(font.main=2, cex.main=1.5,cex.lab=0.4, cex.sub=1)
# Use the ColorBrewer library for color ramps
library(RColorBrewer)
precip.pal = colorRampPalette(c("yellow", "white", "blue"))

#### Countor map with values : Cokriging
ok08ps=spplot(rain.oK08,  zcol='var1.pred',col.regions=precip.pal, contour=TRUE, scales=list(draw=T),
              col='black', pretty=TRUE,main="2008", 
              labels = list(cex = 0.7), label.style = 'align', margin = TRUE,   width = 2, cex = 1.5,  xlab="Longitude",
              ylab="Latitude")
ok08ps
#Year 2009-------
year2009 <-subset(byyear, year==2009)
attach(year2009)
year2009$x1=lat.standard=(2*lat-(max(lat)+min(lat)))/(max(lat)-min(lat))
year2009$x2=lon.standard=(2*lon-(max(lon)+min(lon)))/(max(lon)-min(lon))

year2009$X2=year2009$x2^2
year2009$X1=year2009$x1^2
year2009$X13=year2009$x1^3
year2009$X23=year2009$x2^3
year2009$X2x1=year2009$X2*year2009$x1
year2009$X1x2=year2009$X1*year2009$x2

rain2009mean=lm(anom~1, year2009)
anova(rain2009mean)
summary(rain2009mean)

rain2009linear=lm(anom~x1+x2, year2009)
anova(rain2009linear)
summary(rain2009linear)

rain2009quad=lm(anom~x1*x2+X1+X2, year2009)
anova(rain2009quad)
summary(rain2009quad)

rain2009quadalt=lm(anom~x1*x2+X1+X2+alt, year2009)
anova(rain2009quadalt)
summary(rain2009quadalt)

rain2009cubic=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2, year2009)
anova(rain2009cubic)
summary(rain2009cubic)

rain2009cubicalt=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2+alt, year2009)
anova(rain2009cubicalt)
summary(rain2009cubicalt)


year2009$fitted1=fitted(rain2009cubic)
year2009$resid1=resid(rain2009cubic) #List of residuals
plot(density(resid(rain2009cubic))) #A density plot
qqnorm(resid(rain2009cubic)) # A quantile normal plot - good for checking normality
qqline(resid(rain2009cubic))
year2009fit=data.frame(cbind(year2009$lon, year2009$lat,year2009$fitted1))
year2009res=data.frame(cbind(year2009$lon, year2009$lat,year2009$resid1))

anova(rain2009mean,rain2009linear, rain2009quad, rain2009cubic, rain2009cubicalt)

library(lattice)
w09=wireframe(X3~X1*X2,data=year2009fit,
              xlab = "Longitude", ylab = "Latitude", zlab="yhat", 
              main = "2009", scales = list(arrows = FALSE, distance=2),
              drape = TRUE,aspect = c(70/87, 0.41),
              colorkey = TRUE, screen = list(z = 30, x = -25, y = 10))
w09
c09= contourplot(X3~X1*X2,data=year2009fit,
                 cuts = 20, region = TRUE,
                 xlab = "Longitude",
                 ylab = "Latitude",
                 main = "2009")
c09 

attach(year2009fit)
X1=X2=seq(-1,1,.1)
model=function(X1,X2){13.393-7.824*X1 -138.389*X2-60.968 *X1^2-28.688*X2^2
  +28.689*X1^3+66.610*X2^3-64.209*X2^2*X1+19.485*X1^2*X2-46.510*X1*X2}
model(X1,X2)

z=outer(X1,X2,model)
contour(X1,X2,z,nlevels=20, main="2009")



wireframe(X3~X1*X2,data=year2009res,
          xlab = "Lon", ylab = "Lat", zlab="ehat", 
          main = "Precip.residuals year 2009", scales = list(arrows = FALSE, distance=2),
          drape = TRUE,aspect = c(70/87, 0.22),
          colorkey = TRUE, screen = list(z = 50, x = -55, y = 2))

## Kriging 2009---- 
coordinates(year2009res)=~X1 +X2
proj4string(year2009res) <- "+proj=utm +zone=37 +datum=WGS84 +units=m"
class(year2009res)
summary(year2009res)

# For working with spatial (and spatio-temporal) data, we use the gstat package, which includes functionality for kriging, among other many things.
library(sp)
library(gstat)
# packages for manipulation & visualization
suppressPackageStartupMessages({
  library(dplyr) # for "glimpse"
  library(ggplot2)
  library(scales) # for "comma"
  library(magrittr)
})


##a variogram could be fit as simply as the following code:
vgm <- variogram(X3~1, year2009res) # calculates sample variogram values 
plot(vgm)
library(gstat)
vgm.fit09<- fit.variogram(vgm, model=vgm(270, "Sph",2,50)) # fit model
plot(vgm, vgm.fit09)

# Arrange the data for the ggplot2 plot
# add the semivariance values of v2 to v1
Fitted <- data.frame(dist = seq(0.12, max(vgm$dist), length = 1272))
Fitted$Spherical <- variogramLine(vgm.fit09, dist_vector = Fitted$dist)$gamma
#convert the dataframes to a long format
library(reshape2)
Empirical <- melt(vgm, id.vars = "dist", measure.vars = "gamma")
Modeled <- melt(Fitted, id.vars = "dist", measure.vars = "Spherical")
library(ggplot2)
#both variogram on the same plot with different colours
v09=ggplot(Modeled, aes(x = dist, y = value, colour = variable )) + 
  geom_line() +  
  xlab("Distance") +
  ylab("Semivariance") +   ggtitle("2009")+ 
  geom_point(data = Empirical)
v09

library(magrittr)
library(gstat)
library(sp)
rain.oK09<- krige(X3 ~ 1, year2009res, grd,model = vgm.fit09)
x <- krige.cv(X3 ~ 1, year2009res, grd,model = vgm.fitg09, nfold=100)
#RMSE (root mean squared error)
sqrt(mean(x$residual^2)) 

library(ggplot2)
library(scales)
k09ps=rain.oK09 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKpred"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k09ps=k09ps+xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2009")
k09ps 

k09sd=rain.oK09 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.var)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKerror"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k09sd=k09sd+  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2009")
k09sd

# Some parameters for plotting
par(font.main=2, cex.main=1.5,cex.lab=0.4, cex.sub=1)
# Use the ColorBrewer library for color ramps
library(RColorBrewer)
precip.pal = colorRampPalette(c("yellow", "white", "blue"))

#### Countor map with values : Cokriging
ok09ps=spplot(rain.oK09,  zcol='var1.pred',col.regions=precip.pal, contour=TRUE, scales=list(draw=T),
              col='black', pretty=TRUE,main="2009", 
              labels = list(cex = 0.7), label.style = 'align', margin = TRUE,   width = 2, cex = 1.5,  xlab="Longitude",
              ylab="Latitude")
ok09ps
#Year 2010-------
year2010 <-subset(byyear, year==2010)
attach(year2010)
year2010$x1=lat.standard=(2*lat-(max(lat)+min(lat)))/(max(lat)-min(lat))
year2010$x2=lon.standard=(2*lon-(max(lon)+min(lon)))/(max(lon)-min(lon))

year2010$X2=year2010$x2^2
year2010$X1=year2010$x1^2
year2010$X13=year2010$x1^3
year2010$X23=year2010$x2^3
year2010$X2x1=year2010$X2*year2010$x1
year2010$X1x2=year2010$X1*year2010$x2

rain2010mean=lm(anom~1, year2010)
anova(rain2010mean)
summary(rain2010mean)

rain2010linear=lm(anom~x1+x2, year2010)
anova(rain2010linear)
summary(rain2010linear)

rain2010quad=lm(anom~x1*x2+X1+X2, year2010)
anova(rain2010quad)
summary(rain2010quad)

rain2010quadalt=lm(anom~x1*x2+X1+X2+alt, year2010)
anova(rain2010quadalt)
summary(rain2010quadalt)

rain2010cubic=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2, year2010)
anova(rain2010cubic)
summary(rain2010cubic)

rain2010cubicalt=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2+alt, year2010)
anova(rain2010cubicalt)
summary(rain2010cubicalt)


year2010$fitted1=fitted(rain2010cubic)
year2010$resid1=resid(rain2010cubic) #List of residuals
plot(density(resid(rain2010cubic))) #A density plot
qqnorm(resid(rain2010cubic)) # A quantile normal plot - good for checking normality
qqline(resid(rain2010cubic))
year2010fit=data.frame(cbind(year2010$lon, year2010$lat,year2010$fitted1))
year2010res=data.frame(cbind(year2010$lon, year2010$lat,year2010$resid1))

anova(rain2010mean,rain2010linear, rain2010quad, rain2010cubic, rain2010cubicalt)

library(lattice)
w10=wireframe(X3~X1*X2,data=year2010fit,
              xlab = "Longitude", ylab = "Latitude", zlab="yhat", 
              main = "2010", scales = list(arrows = FALSE, distance=2),
              drape = TRUE,aspect = c(70/87, 0.41),
              colorkey = TRUE, screen = list(z = 30, x = -25, y = 10))
w10
c10= contourplot(X3~X1*X2,data=year2010fit,
                 cuts = 20, region = TRUE,
                 xlab = "Longitude",
                 ylab = "Latitude",
                 main = "2010")
c10 

attach(year2010fit)
X1=X2=seq(-1,1,.1)
model=function(X1,X2){13.393-7.824*X1 -138.389*X2-60.968 *X1^2-28.688*X2^2
  +28.689*X1^3+66.610*X2^3-64.210*X2^2*X1+19.485*X1^2*X2-46.510*X1*X2}
model(X1,X2)

z=outer(X1,X2,model)
contour(X1,X2,z,nlevels=20, main="2010")



wireframe(X3~X1*X2,data=year2010res,
          xlab = "Lon", ylab = "Lat", zlab="ehat", 
          main = "Precip.residuals year 2010", scales = list(arrows = FALSE, distance=2),
          drape = TRUE,aspect = c(70/87, 0.22),
          colorkey = TRUE, screen = list(z = 50, x = -55, y = 2))

## Kriging 2010---- 
coordinates(year2010res)=~X1 +X2
proj4string(year2010res) <- "+proj=utm +zone=37 +datum=WGS84 +units=m"
class(year2010res)
summary(year2010res)

# For working with spatial (and spatio-temporal) data, we use the gstat package, which includes functionality for kriging, among other many things.
library(sp)
library(gstat)
# packages for manipulation & visualization
suppressPackageStartupMessages({
  library(dplyr) # for "glimpse"
  library(ggplot2)
  library(scales) # for "comma"
  library(magrittr)
})


##a variogram could be fit as simply as the following code:
vgm <- variogram(X3~1, year2010res) # calculates sample variogram values 
plot(vgm)
library(gstat)
vgm.fit10<- fit.variogram(vgm, model=vgm(255, "Sph",1.7,45)) # fit model
plot(vgm, vgm.fit10)

# Arrange the data for the ggplot2 plot
# add the semivariance values of v2 to v1
Fitted <- data.frame(dist = seq(0.12, max(vgm$dist), length = 1272))
Fitted$Spherical <- variogramLine(vgm.fit10, dist_vector = Fitted$dist)$gamma
#convert the dataframes to a long format
library(reshape2)
Empirical <- melt(vgm, id.vars = "dist", measure.vars = "gamma")
Modeled <- melt(Fitted, id.vars = "dist", measure.vars = "Spherical")
library(ggplot2)
#both variogram on the same plot with different colours
v10=ggplot(Modeled, aes(x = dist, y = value, colour = variable )) + 
  geom_line() +  
  xlab("Distance") +
  ylab("Semivariance") +   ggtitle("2010")+ 
  geom_point(data = Empirical)
v10

library(magrittr)
library(gstat)
library(sp)
rain.oK10<- krige(X3 ~ 1, year2010res, grd,model = vgm.fit10)
x <- krige.cv(X3 ~ 1, year2010res, grd,model = vgm.fit10, nfold=110)
#RMSE (root mean squared error)
sqrt(mean(x$residual^2)) 

library(ggplot2)
library(scales)
k10ps=rain.oK10 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKpred"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k10ps=k10ps+xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2010")
k10ps 

k10sd=rain.oK10 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.var)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKerror"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k10sd=k10sd+  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2010")
k10sd

# Some parameters for plotting
par(font.main=2, cex.main=1.5,cex.lab=0.4, cex.sub=1)
# Use the ColorBrewer library for color ramps
library(RColorBrewer)
precip.pal = colorRampPalette(c("yellow", "white", "blue"))

#### Countor map with values : Cokriging
ok10ps=spplot(rain.oK10,  zcol='var1.pred',col.regions=precip.pal, contour=TRUE, scales=list(draw=T),
              col='black', pretty=TRUE,main="2010", 
              labels = list(cex = 0.7), label.style = 'align', margin = TRUE,   width = 2, cex = 1.5,  xlab="Longitude",
              ylab="Latitude")
ok10ps
#Year 2011-------
year2011 <-subset(byyear, year==2011)
attach(year2011)
year2011$x1=lat.standard=(2*lat-(max(lat)+min(lat)))/(max(lat)-min(lat))
year2011$x2=lon.standard=(2*lon-(max(lon)+min(lon)))/(max(lon)-min(lon))

year2011$X2=year2011$x2^2
year2011$X1=year2011$x1^2
year2011$X13=year2011$x1^3
year2011$X23=year2011$x2^3
year2011$X2x1=year2011$X2*year2011$x1
year2011$X1x2=year2011$X1*year2011$x2

rain2011mean=lm(anom~1, year2011)
anova(rain2011mean)
summary(rain2011mean)

rain2011linear=lm(anom~x1+x2, year2011)
anova(rain2011linear)
summary(rain2011linear)

rain2011quad=lm(anom~x1*x2+X1+X2, year2011)
anova(rain2011quad)
summary(rain2011quad)

rain2011quadalt=lm(anom~x1*x2+X1+X2+alt, year2011)
anova(rain2011quadalt)
summary(rain2011quadalt)

rain2011cubic=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2, year2011)
anova(rain2011cubic)
summary(rain2011cubic)

rain2011cubicalt=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2+alt, year2011)
anova(rain2011cubicalt)
summary(rain2011cubicalt)


year2011$fitted1=fitted(rain2011cubic)
year2011$resid1=resid(rain2011cubic) #List of residuals
plot(density(resid(rain2011cubic))) #A density plot
qqnorm(resid(rain2011cubic)) # A quantile normal plot - good for checking normality
qqline(resid(rain2011cubic))
year2011fit=data.frame(cbind(year2011$lon, year2011$lat,year2011$fitted1))
year2011res=data.frame(cbind(year2011$lon, year2011$lat,year2011$resid1))

anova(rain2011mean,rain2011linear, rain2011quad, rain2011cubic, rain2011cubicalt)

library(lattice)
w11=wireframe(X3~X1*X2,data=year2011fit,
              xlab = "Longitude", ylab = "Latitude", zlab="yhat", 
              main = "2011", scales = list(arrows = FALSE, distance=2),
              drape = TRUE,aspect = c(70/87, 0.41),
              colorkey = TRUE, screen = list(z = 30, x = -25, y = 11))
w11
c11= contourplot(X3~X1*X2,data=year2011fit,
                 cuts = 20, region = TRUE,
                 xlab = "Longitude",
                 ylab = "Latitude",
                 main = "2011")
c11 

attach(year2011fit)
X1=X2=seq(-1,1,.1)
model=function(X1,X2){13.393-7.824*X1 -138.389*X2-60.968 *X1^2-28.688*X2^2
  +28.689*X1^3+66.611*X2^3-64.211*X2^2*X1+19.485*X1^2*X2-46.511*X1*X2}
model(X1,X2)

z=outer(X1,X2,model)
contour(X1,X2,z,nlevels=20, main="2011")



wireframe(X3~X1*X2,data=year2011res,
          xlab = "Lon", ylab = "Lat", zlab="ehat", 
          main = "Precip.residuals year 2011", scales = list(arrows = FALSE, distance=2),
          drape = TRUE,aspect = c(70/87, 0.22),
          colorkey = TRUE, screen = list(z = 50, x = -55, y = 2))

## Kriging 2011---- 
coordinates(year2011res)=~X1 +X2
proj4string(year2011res) <- "+proj=utm +zone=37 +datum=WGS84 +units=m"
class(year2011res)
summary(year2011res)

# For working with spatial (and spatio-temporal) data, we use the gstat package, which includes functionality for kriging, among other many things.
library(sp)
library(gstat)
# packages for manipulation & visualization
suppressPackageStartupMessages({
  library(dplyr) # for "glimpse"
  library(ggplot2)
  library(scales) # for "comma"
  library(magrittr)
})


##a variogram could be fit as simply as the following code:
vgm <- variogram(X3~1, year2011res) # calculates sample variogram values 
plot(vgm)
library(gstat)
vgm.fit11<- fit.variogram(vgm, model=vgm(337, "Sph",1.6,100)) # fit model
plot(vgm, vgm.fit11)

# Arrange the data for the ggplot2 plot
# add the semivariance values of v2 to v1
Fitted <- data.frame(dist = seq(0.12, max(vgm$dist), length = 1272))
Fitted$Spherical <- variogramLine(vgm.fit11, dist_vector = Fitted$dist)$gamma
#convert the dataframes to a long format
library(reshape2)
Empirical <- melt(vgm, id.vars = "dist", measure.vars = "gamma")
Modeled <- melt(Fitted, id.vars = "dist", measure.vars = "Spherical")
library(ggplot2)
#both variogram on the same plot with different colours
v11=ggplot(Modeled, aes(x = dist, y = value, colour = variable )) + 
  geom_line() +  
  xlab("Distance") +
  ylab("Semivariance") +   ggtitle("2011")+ 
  geom_point(data = Empirical)
v11

library(magrittr)
library(gstat)
library(sp)
rain.oK11<- krige(X3 ~ 1, year2011res, grd,model = vgm.fit11)
x <- krige.cv(X3 ~ 1, year2011res, grd,model = vgm.fit11, nfold=111)
#RMSE (root mean squared error)
sqrt(mean(x$residual^2)) 

library(ggplot2)
library(scales)
k11ps=rain.oK11 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKpred"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k11ps=k11ps+xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2011")
k11ps 

k11sd=rain.oK11 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.var)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKerror"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k11sd=k11sd+  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2011")
k11sd

# Some parameters for plotting
par(font.main=2, cex.main=1.5,cex.lab=0.4, cex.sub=1)
# Use the ColorBrewer library for color ramps
library(RColorBrewer)
precip.pal = colorRampPalette(c("yellow", "white", "blue"))

#### Countor map with values : Cokriging
ok11ps=spplot(rain.oK11,  zcol='var1.pred',col.regions=precip.pal, contour=TRUE, scales=list(draw=T),
              col='black', pretty=TRUE,main="2011", 
              labels = list(cex = 0.7), label.style = 'align', margin = TRUE,   width = 2, cex = 1.5,  xlab="Longitude",
              ylab="Latitude")
ok11ps
#Year 2012-------
year2012 <-subset(byyear, year==2012)
attach(year2012)
year2012$x1=lat.standard=(2*lat-(max(lat)+min(lat)))/(max(lat)-min(lat))
year2012$x2=lon.standard=(2*lon-(max(lon)+min(lon)))/(max(lon)-min(lon))

year2012$X2=year2012$x2^2
year2012$X1=year2012$x1^2
year2012$X13=year2012$x1^3
year2012$X23=year2012$x2^3
year2012$X2x1=year2012$X2*year2012$x1
year2012$X1x2=year2012$X1*year2012$x2

rain2012mean=lm(anom~1, year2012)
anova(rain2012mean)
summary(rain2012mean)

rain2012linear=lm(anom~x1+x2, year2012)
anova(rain2012linear)
summary(rain2012linear)

rain2012quad=lm(anom~x1*x2+X1+X2, year2012)
anova(rain2012quad)
summary(rain2012quad)

rain2012quadalt=lm(anom~x1*x2+X1+X2+alt, year2012)
anova(rain2012quadalt)
summary(rain2012quadalt)

rain2012cubic=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2, year2012)
anova(rain2012cubic)
summary(rain2012cubic)

rain2012cubicalt=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2+alt, year2012)
anova(rain2012cubicalt)
summary(rain2012cubicalt)


year2012$fitted1=fitted(rain2012cubic)
year2012$resid1=resid(rain2012cubic) #List of residuals
plot(density(resid(rain2012cubic))) #A density plot
qqnorm(resid(rain2012cubic)) # A quantile normal plot - good for checking normality
qqline(resid(rain2012cubic))
year2012fit=data.frame(cbind(year2012$lon, year2012$lat,year2012$fitted1))
year2012res=data.frame(cbind(year2012$lon, year2012$lat,year2012$resid1))

anova(rain2012mean,rain2012linear, rain2012quad, rain2012cubic, rain2012cubicalt)

library(lattice)
w12=wireframe(X3~X1*X2,data=year2012fit,
              xlab = "Longitude", ylab = "Latitude", zlab="yhat", 
              main = "2012", scales = list(arrows = FALSE, distance=2),
              drape = TRUE,aspect = c(70/87, 0.41),
              colorkey = TRUE, screen = list(z = 30, x = -25, y = 12))
w12
c12= contourplot(X3~X1*X2,data=year2012fit,
                 cuts = 20, region = TRUE,
                 xlab = "Longitude",
                 ylab = "Latitude",
                 main = "2012")
c12 

attach(year2012fit)
X1=X2=seq(-1,1,.1)
model=function(X1,X2){13.393-7.824*X1 -138.389*X2-60.968 *X1^2-28.688*X2^2
  +28.689*X1^3+66.612*X2^3-64.212*X2^2*X1+19.485*X1^2*X2-46.512*X1*X2}
model(X1,X2)

z=outer(X1,X2,model)
contour(X1,X2,z,nlevels=20, main="2012")



wireframe(X3~X1*X2,data=year2012res,
          xlab = "Lon", ylab = "Lat", zlab="ehat", 
          main = "Precip.residuals year 2012", scales = list(arrows = FALSE, distance=2),
          drape = TRUE,aspect = c(70/87, 0.22),
          colorkey = TRUE, screen = list(z = 50, x = -55, y = 2))

## Kriging 2012---- 
coordinates(year2012res)=~X1 +X2
proj4string(year2012res) <- "+proj=utm +zone=37 +datum=WGS84 +units=m"
class(year2012res)
summary(year2012res)

# For working with spatial (and spatio-temporal) data, we use the gstat package, which includes functionality for kriging, among other many things.
library(sp)
library(gstat)
# packages for manipulation & visualization
suppressPackageStartupMessages({
  library(dplyr) # for "glimpse"
  library(ggplot2)
  library(scales) # for "comma"
  library(magrittr)
})


##a variogram could be fit as simply as the following code:
vgm <- variogram(X3~1, year2012res) # calculates sample variogram values 
plot(vgm)
library(gstat)
vgm.fit12<- fit.variogram(vgm, model=vgm(380, "Sph",1.5,200)) # fit model
plot(vgm, vgm.fit12)

# Arrange the data for the ggplot2 plot
# add the semivariance values of v2 to v1
Fitted <- data.frame(dist = seq(0.12, max(vgm$dist), length = 1272))
Fitted$Spherical <- variogramLine(vgm.fit12, dist_vector = Fitted$dist)$gamma
#convert the dataframes to a long format
library(reshape2)
Empirical <- melt(vgm, id.vars = "dist", measure.vars = "gamma")
Modeled <- melt(Fitted, id.vars = "dist", measure.vars = "Spherical")
library(ggplot2)
#both variogram on the same plot with different colours
v12=ggplot(Modeled, aes(x = dist, y = value, colour = variable )) + 
  geom_line() +  
  xlab("Distance") +
  ylab("Semivariance") +   ggtitle("2012")+ 
  geom_point(data = Empirical)
v12

library(magrittr)
library(gstat)
library(sp)
rain.oK12<- krige(X3 ~ 1, year2012res, grd,model = vgm.fit12)
x <- krige.cv(X3 ~ 1, year2012res, grd,model = vgm.fit12, nfold=100)
#RMSE (root mean squared error)
sqrt(mean(x$residual^2)) 

library(ggplot2)
library(scales)
k12ps=rain.oK12 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKpred"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k12ps=k12ps+xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2012")
k12ps 

k12sd=rain.oK12 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.var)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKerror"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k12sd=k12sd+  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2012")
k12sd

# Some parameters for plotting
par(font.main=2, cex.main=1.5,cex.lab=0.4, cex.sub=1)
# Use the ColorBrewer library for color ramps
library(RColorBrewer)
precip.pal = colorRampPalette(c("yellow", "white", "blue"))

#### Countor map with values : Cokriging
ok12ps=spplot(rain.oK12,  zcol='var1.pred',col.regions=precip.pal, contour=TRUE, scales=list(draw=T),
              col='black', pretty=TRUE,main="2012", 
              labels = list(cex = 0.7), label.style = 'align', margin = TRUE,   width = 2, cex = 1.5,  xlab="Longitude",
              ylab="Latitude")
ok12ps
#Year 2013-------
year2013 <-subset(byyear, year==2013)
attach(year2013)
year2013$x1=lat.standard=(2*lat-(max(lat)+min(lat)))/(max(lat)-min(lat))
year2013$x2=lon.standard=(2*lon-(max(lon)+min(lon)))/(max(lon)-min(lon))

year2013$X2=year2013$x2^2
year2013$X1=year2013$x1^2
year2013$X13=year2013$x1^3
year2013$X23=year2013$x2^3
year2013$X2x1=year2013$X2*year2013$x1
year2013$X1x2=year2013$X1*year2013$x2

rain2013mean=lm(anom~1, year2013)
anova(rain2013mean)
summary(rain2013mean)

rain2013linear=lm(anom~x1+x2, year2013)
anova(rain2013linear)
summary(rain2013linear)

rain2013quad=lm(anom~x1*x2+X1+X2, year2013)
anova(rain2013quad)
summary(rain2013quad)

rain2013quadalt=lm(anom~x1*x2+X1+X2+alt, year2013)
anova(rain2013quadalt)
summary(rain2013quadalt)

rain2013cubic=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2, year2013)
anova(rain2013cubic)
summary(rain2013cubic)

rain2013cubicalt=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2+alt, year2013)
anova(rain2013cubicalt)
summary(rain2013cubicalt)


year2013$fitted1=fitted(rain2013cubic)
year2013$resid1=resid(rain2013cubic) #List of residuals
plot(density(resid(rain2013cubic))) #A density plot
qqnorm(resid(rain2013cubic)) # A quantile normal plot - good for checking normality
qqline(resid(rain2013cubic))
year2013fit=data.frame(cbind(year2013$lon, year2013$lat,year2013$fitted1))
year2013res=data.frame(cbind(year2013$lon, year2013$lat,year2013$resid1))

anova(rain2013mean,rain2013linear, rain2013quad, rain2013cubic, rain2013cubicalt)

library(lattice)
w13=wireframe(X3~X1*X2,data=year2013fit,
              xlab = "Longitude", ylab = "Latitude", zlab="yhat", 
              main = "2013", scales = list(arrows = FALSE, distance=2),
              drape = TRUE,aspect = c(70/87, 0.41),
              colorkey = TRUE, screen = list(z = 30, x = -25, y = 13))
w13
c13= contourplot(X3~X1*X2,data=year2013fit,
                 cuts = 20, region = TRUE,
                 xlab = "Longitude",
                 ylab = "Latitude",
                 main = "2013")
c13 

attach(year2013fit)
X1=X2=seq(-1,1,.1)
model=function(X1,X2){13.393-7.824*X1 -138.389*X2-60.968 *X1^2-28.688*X2^2
  +28.689*X1^3+66.613*X2^3-64.213*X2^2*X1+19.485*X1^2*X2-46.513*X1*X2}
model(X1,X2)

z=outer(X1,X2,model)
contour(X1,X2,z,nlevels=20, main="2013")



wireframe(X3~X1*X2,data=year2013res,
          xlab = "Lon", ylab = "Lat", zlab="ehat", 
          main = "Precip.residuals year 2013", scales = list(arrows = FALSE, distance=2),
          drape = TRUE,aspect = c(70/87, 0.22),
          colorkey = TRUE, screen = list(z = 50, x = -55, y = 2))

## Kriging 2013---- 
coordinates(year2013res)=~X1 +X2
proj4string(year2013res) <- "+proj=utm +zone=37 +datum=WGS84 +units=m"
class(year2013res)
summary(year2013res)

# For working with spatial (and spatio-temporal) data, we use the gstat package, which includes functionality for kriging, among other many things.
library(sp)
library(gstat)
# packages for manipulation & visualization
suppressPackageStartupMessages({
  library(dplyr) # for "glimpse"
  library(ggplot2)
  library(scales) # for "comma"
  library(magrittr)
})


##a variogram could be fit as simply as the following code:
vgm <- variogram(X3~1, year2013res) # calculates sample variogram values 
plot(vgm)
library(gstat)
vgm.fit13<- fit.variogram(vgm, model=vgm(253, "Sph",2,70)) # fit model
plot(vgm, vgm.fit13)

# Arrange the data for the ggplot2 plot
# add the semivariance values of v2 to v1
Fitted <- data.frame(dist = seq(0.13, max(vgm$dist), length = 1372))
Fitted$Spherical <- variogramLine(vgm.fit13, dist_vector = Fitted$dist)$gamma
#convert the dataframes to a long format
library(reshape2)
Empirical <- melt(vgm, id.vars = "dist", measure.vars = "gamma")
Modeled <- melt(Fitted, id.vars = "dist", measure.vars = "Spherical")
library(ggplot2)
#both variogram on the same plot with different colours
v13=ggplot(Modeled, aes(x = dist, y = value, colour = variable )) + 
  geom_line() +  
  xlab("Distance") +
  ylab("Semivariance") +   ggtitle("2013")+ 
  geom_point(data = Empirical)
v13

library(magrittr)
library(gstat)
library(sp)
rain.oK13<- krige(X3 ~ 1, year2013res, grd,model = vgm.fit13)
x <- krige.cv(X3 ~ 1, year2013res, grd,model = vgm.fit13, nfold=100)
#RMSE (root mean squared error)
sqrt(mean(x$residual^2)) 

library(ggplot2)
library(scales)
k13ps=rain.oK13 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKpred"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k13ps=k13ps+xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2013")
k13ps 

k13sd=rain.oK13 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.var)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKerror"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k13sd=k13sd+  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2013")
k13sd

# Some parameters for plotting
par(font.main=2, cex.main=1.5,cex.lab=0.4, cex.sub=1)
# Use the ColorBrewer library for color ramps
library(RColorBrewer)
precip.pal = colorRampPalette(c("yellow", "white", "blue"))

#### Countor map with values : Cokriging
ok13ps=spplot(rain.oK13,  zcol='var1.pred',col.regions=precip.pal, contour=TRUE, scales=list(draw=T),
              col='black', pretty=TRUE,main="2013", 
              labels = list(cex = 0.7), label.style = 'align', margin = TRUE,   width = 2, cex = 1.5,  xlab="Longitude",
              ylab="Latitude")
ok13ps
#Year 2014-------
year2014 <-subset(byyear, year==2014)
attach(year2014)
year2014$x1=lat.standard=(2*lat-(max(lat)+min(lat)))/(max(lat)-min(lat))
year2014$x2=lon.standard=(2*lon-(max(lon)+min(lon)))/(max(lon)-min(lon))

year2014$X2=year2014$x2^2
year2014$X1=year2014$x1^2
year2014$X13=year2014$x1^3
year2014$X23=year2014$x2^3
year2014$X2x1=year2014$X2*year2014$x1
year2014$X1x2=year2014$X1*year2014$x2

rain2014mean=lm(anom~1, year2014)
anova(rain2014mean)
summary(rain2014mean)

rain2014linear=lm(anom~x1+x2, year2014)
anova(rain2014linear)
summary(rain2014linear)

rain2014quad=lm(anom~x1*x2+X1+X2, year2014)
anova(rain2014quad)
summary(rain2014quad)

rain2014quadalt=lm(anom~x1*x2+X1+X2+alt, year2014)
anova(rain2014quadalt)
summary(rain2014quadalt)

rain2014cubic=lm(anom~x1 * x2 + X1 + X2 + X13 + X23 + X2x1 + X1x2, year2014)
anova(rain2014cubic)
summary(rain2014cubic)

rain2014cubicalt=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2+alt, year2014)
anova(rain2014cubicalt)
summary(rain2014cubicalt)


year2014$fitted1=fitted(rain2014cubic)
year2014$resid1=resid(rain2014cubic) #List of residuals
plot(density(resid(rain2014cubic))) #A density plot
qqnorm(resid(rain2014cubic)) # A quantile normal plot - good for checking normality
qqline(resid(rain2014cubic))
year2014fit=data.frame(cbind(year2014$lon, year2014$lat,year2014$fitted1))
year2014res=data.frame(cbind(year2014$lon, year2014$lat,year2014$resid1))

anova(rain2014mean,rain2014linear, rain2014quad, rain2014cubic, rain2014cubicalt)

library(lattice)
w14=wireframe(X3~X1*X2,data=year2014fit,
              xlab = "Longitude", ylab = "Latitude", zlab="yhat", 
              main = "2014", scales = list(arrows = FALSE, distance=2),
              drape = TRUE,aspect = c(70/87, 0.41),
              colorkey = TRUE, screen = list(z = 30, x = -25, y = 14))
w14
c14= contourplot(X3~X1*X2,data=year2014fit,
                 cuts = 20, region = TRUE,
                 xlab = "Longitude",
                 ylab = "Latitude",
                 main = "2014")
c14 

attach(year2014fit)
X1=X2=seq(-1,1,.1)
model=function(X1,X2){14.393-7.824*X1 -148.389*X2-60.968 *X1^2-28.688*X2^2
  +28.689*X1^3+66.614*X2^3-64.214*X2^2*X1+19.485*X1^2*X2-46.514*X1*X2}
model(X1,X2)

z=outer(X1,X2,model)
contour(X1,X2,z,nlevels=20, main="2014")



wireframe(X3~X1*X2,data=year2014res,
          xlab = "Lon", ylab = "Lat", zlab="ehat", 
          main = "Precip.residuals year 2014", scales = list(arrows = FALSE, distance=2),
          drape = TRUE,aspect = c(70/87, 0.22),
          colorkey = TRUE, screen = list(z = 50, x = -55, y = 2))

## Kriging 2014---- 
coordinates(year2014res)=~X1 +X2
proj4string(year2014res) <- "+proj=utm +zone=37 +datum=WGS84 +units=m"
class(year2014res)
summary(year2014res)

# For working with spatial (and spatio-temporal) data, we use the gstat package, which includes functionality for kriging, among other many things.
library(sp)
library(gstat)
# packages for manipulation & visualization
suppressPackageStartupMessages({
  library(dplyr) # for "glimpse"
  library(ggplot2)
  library(scales) # for "comma"
  library(magrittr)
})


##a variogram could be fit as simply as the following code:
vgm <- variogram(X3~1, year2014res) # calculates sample variogram values 
plot(vgm)
library(gstat)
vgm.fit14<- fit.variogram(vgm, model=vgm(422, "Sph",2.2,100)) # fit model
plot(vgm, vgm.fit14)

# Arrange the data for the ggplot2 plot
# add the semivariance values of v2 to v1
Fitted <- data.frame(dist = seq(0.14, max(vgm$dist), length = 1472))
Fitted$Spherical <- variogramLine(vgm.fit14, dist_vector = Fitted$dist)$gamma
#convert the dataframes to a long format
library(reshape2)
Empirical <- melt(vgm, id.vars = "dist", measure.vars = "gamma")
Modeled <- melt(Fitted, id.vars = "dist", measure.vars = "Spherical")
library(ggplot2)
#both variogram on the same plot with different colours
v14=ggplot(Modeled, aes(x = dist, y = value, colour = variable )) + 
  geom_line() +  
  xlab("Distance") +
  ylab("Semivariance") +   ggtitle("2014")+ 
  geom_point(data = Empirical)
v14

library(magrittr)
library(gstat)
library(sp)
rain.oK14<- krige(X3 ~ 1, year2014res, grd,model = vgm.fit14)
x <- krige.cv(X3 ~ 1, year2014res, grd,model = vgm.fit14, nfold=100)
#RMSE (root mean squared error)
sqrt(mean(x$residual^2)) 

library(ggplot2)
library(scales)
k14ps=rain.oK14 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKpred"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k14ps=k14ps+xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2014")
k14ps 

k14sd=rain.oK14 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.var)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKerror"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k14sd=k14sd+  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2014")
k14sd

# Some parameters for plotting
par(font.main=2, cex.main=1.5,cex.lab=0.4, cex.sub=1)
# Use the ColorBrewer library for color ramps
library(RColorBrewer)
precip.pal = colorRampPalette(c("yellow", "white", "blue"))

#### Countor map with values : Cokriging
ok14ps=spplot(rain.oK14,  zcol='var1.pred',col.regions=precip.pal, contour=TRUE, scales=list(draw=T),
              col='black', pretty=TRUE,main="2014", 
              labels = list(cex = 0.7), label.style = 'align', margin = TRUE,   width = 2, cex = 1.5,  xlab="Longitude",
              ylab="Latitude")
ok14ps
#Year 2015-------
year2015 <-subset(byyear, year==2015)
attach(year2015)
year2015$x1=lat.standard=(2*lat-(max(lat)+min(lat)))/(max(lat)-min(lat))
year2015$x2=lon.standard=(2*lon-(max(lon)+min(lon)))/(max(lon)-min(lon))

year2015$X2=year2015$x2^2
year2015$X1=year2015$x1^2
year2015$X13=year2015$x1^3
year2015$X23=year2015$x2^3
year2015$X2x1=year2015$X2*year2015$x1
year2015$X1x2=year2015$X1*year2015$x2

rain2015mean=lm(anom~1, year2015)
anova(rain2015mean)
summary(rain2015mean)

rain2015linear=lm(anom~x1+x2, year2015)
anova(rain2015linear)
summary(rain2015linear)

rain2015quad=lm(anom~x1*x2+X1+X2, year2015)
anova(rain2015quad)
summary(rain2015quad)

rain2015quadalt=lm(anom~x1*x2+X1+X2+alt, year2015)
anova(rain2015quadalt)
summary(rain2015quadalt)

rain2015cubic=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2, year2015)
anova(rain2015cubic)
summary(rain2015cubic)

rain2015cubicalt=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2+alt, year2015)
anova(rain2015cubicalt)
summary(rain2015cubicalt)


year2015$fitted1=fitted(rain2015cubic)
year2015$resid1=resid(rain2015cubic) #List of residuals
plot(density(resid(rain2015cubic))) #A density plot
qqnorm(resid(rain2015cubic)) # A quantile normal plot - good for checking normality
qqline(resid(rain2015cubic))
year2015fit=data.frame(cbind(year2015$lon, year2015$lat,year2015$fitted1))
year2015res=data.frame(cbind(year2015$lon, year2015$lat,year2015$resid1))

anova(rain2015mean,rain2015linear, rain2015quad, rain2015cubic, rain2015cubicalt)

library(lattice)
w15=wireframe(X3~X1*X2,data=year2015fit,
              xlab = "Longitude", ylab = "Latitude", zlab="yhat", 
              main = "2015", scales = list(arrows = FALSE, distance=2),
              drape = TRUE,aspect = c(70/87, 0.41),
              colorkey = TRUE, screen = list(z = 30, x = -25, y = 15))
w15
c15= contourplot(X3~X1*X2,data=year2015fit,
                 cuts = 20, region = TRUE,
                 xlab = "Longitude",
                 ylab = "Latitude",
                 main = "2015")
c15 

attach(year2015fit)
X1=X2=seq(-1,1,.1)
model=function(X1,X2){15.393-7.824*X1 -158.389*X2-60.968 *X1^2-28.688*X2^2
  +28.689*X1^3+66.615*X2^3-64.215*X2^2*X1+19.485*X1^2*X2-46.515*X1*X2}
model(X1,X2)

z=outer(X1,X2,model)
contour(X1,X2,z,nlevels=20, main="2015")



wireframe(X3~X1*X2,data=year2015res,
          xlab = "Lon", ylab = "Lat", zlab="ehat", 
          main = "Precip.residuals year 2015", scales = list(arrows = FALSE, distance=2),
          drape = TRUE,aspect = c(70/87, 0.22),
          colorkey = TRUE, screen = list(z = 50, x = -55, y = 2))

## Kriging 2015---- 
coordinates(year2015res)=~X1 +X2
proj4string(year2015res) <- "+proj=utm +zone=37 +datum=WGS84 +units=m"
class(year2015res)
summary(year2015res)

# For working with spatial (and spatio-temporal) data, we use the gstat package, which includes functionality for kriging, among other many things.
library(sp)
library(gstat)
# packages for manipulation & visualization
suppressPackageStartupMessages({
  library(dplyr) # for "glimpse"
  library(ggplot2)
  library(scales) # for "comma"
  library(magrittr)
})


##a variogram could be fit as simply as the following code:
vgm <- variogram(X3~1, year2015res) # calculates sample variogram values 
plot(vgm)
library(gstat)
vgm.fit15<- fit.variogram(vgm, model=vgm(370, "Sph",2,120)) # fit model
plot(vgm, vgm.fit15)

# Arrange the data for the ggplot2 plot
# add the semivariance values of v2 to v1
Fitted <- data.frame(dist = seq(0.15, max(vgm$dist), length = 1572))
Fitted$Spherical <- variogramLine(vgm.fit15, dist_vector = Fitted$dist)$gamma
#convert the dataframes to a long format
library(reshape2)
Empirical <- melt(vgm, id.vars = "dist", measure.vars = "gamma")
Modeled <- melt(Fitted, id.vars = "dist", measure.vars = "Spherical")
library(ggplot2)
#both variogram on the same plot with different colours
v15=ggplot(Modeled, aes(x = dist, y = value, colour = variable )) + 
  geom_line() +  
  xlab("Distance") +
  ylab("Semivariance") +   ggtitle("2015")+ 
  geom_point(data = Empirical)
v15

library(magrittr)
library(gstat)
library(sp)
rain.oK15<- krige(X3 ~ 1, year2015res, grd,model = vgm.fit15)
x <- krige.cv(X3 ~ 1, year2015res, grd,model = vgm.fit15, nfold=151)
#RMSE (root mean squared error)
sqrt(mean(x$residual^2)) 

library(ggplot2)
library(scales)
k15ps=rain.oK15 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKpred"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k15ps=k15ps+xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2015")
k15ps 

k15sd=rain.oK15 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.var)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKerror"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k15sd=k15sd+  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2015")
k15sd

# Some parameters for plotting
par(font.main=2, cex.main=1.5,cex.lab=0.4, cex.sub=1)
# Use the ColorBrewer library for color ramps
library(RColorBrewer)
precip.pal = colorRampPalette(c("yellow", "white", "blue"))

#### Countor map with values : Cokriging
ok15ps=spplot(rain.oK15,  zcol='var1.pred',col.regions=precip.pal, contour=TRUE, scales=list(draw=T),
              col='black', pretty=TRUE,main="2015", 
              labels = list(cex = 0.7), label.style = 'align', margin = TRUE,   width = 2, cex = 1.5,  xlab="Longitude",
              ylab="Latitude")
ok15ps

year2015 <-subset(byyear, year==2015)
attach(year2015)
year2015$x1=lat.standard=(2*lat-(max(lat)+min(lat)))/(max(lat)-min(lat))
year2015$x2=lon.standard=(2*lon-(max(lon)+min(lon)))/(max(lon)-min(lon))

year2015$X2=year2015$x2^2
year2015$X1=year2015$x1^2
year2015$X13=year2015$x1^3
year2015$X23=year2015$x2^3
year2015$X2x1=year2015$X2*year2015$x1
year2015$X1x2=year2015$X1*year2015$x2

rain2015mean=lm(anom~1, year2015)
anova(rain2015mean)
summary(rain2015mean)

rain2015linear=lm(anom~x1+x2, year2015)
anova(rain2015linear)
summary(rain2015linear)

rain2015quad=lm(anom~x1*x2+X1+X2, year2015)
anova(rain2015quad)
summary(rain2015quad)

rain2015quadalt=lm(anom~x1*x2+X1+X2+alt, year2015)
anova(rain2015quadalt)
summary(rain2015quadalt)

rain2015cubic=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2, year2015)
anova(rain2015cubic)
summary(rain2015cubic)

rain2015cubicalt=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2+alt, year2015)
anova(rain2015cubicalt)
summary(rain2015cubicalt)


year2015$fitted1=fitted(rain2015cubic)
year2015$resid1=resid(rain2015cubic) #List of residuals
plot(density(resid(rain2015cubic))) #A density plot
qqnorm(resid(rain2015cubic)) # A quantile normal plot - good for checking normality
qqline(resid(rain2015cubic))
year2015fit=data.frame(cbind(year2015$lon, year2015$lat,year2015$fitted1))
year2015res=data.frame(cbind(year2015$lon, year2015$lat,year2015$resid1))

anova(rain2015mean,rain2015linear, rain2015quad, rain2015cubic, rain2015cubicalt)

library(lattice)
w15=wireframe(X3~X1*X2,data=year2015fit,
              xlab = "Longitude", ylab = "Latitude", zlab="yhat", 
              main = "2015", scales = list(arrows = FALSE, distance=2),
              drape = TRUE,aspect = c(70/87, 0.41),
              colorkey = TRUE, screen = list(z = 30, x = -25, y = 15))
w15
c15= contourplot(X3~X1*X2,data=year2015fit,
                 cuts = 20, region = TRUE,
                 xlab = "Longitude",
                 ylab = "Latitude",
                 main = "2015")
c15 

attach(year2015fit)
X1=X2=seq(-1,1,.1)
model=function(X1,X2){15.393-7.824*X1 -158.389*X2-60.968 *X1^2-28.688*X2^2
  +28.689*X1^3+66.615*X2^3-64.215*X2^2*X1+19.485*X1^2*X2-46.515*X1*X2}
model(X1,X2)

z=outer(X1,X2,model)
contour(X1,X2,z,nlevels=20, main="2015")



wireframe(X3~X1*X2,data=year2015res,
          xlab = "Lon", ylab = "Lat", zlab="ehat", 
          main = "Precip.residuals year 2015", scales = list(arrows = FALSE, distance=2),
          drape = TRUE,aspect = c(70/87, 0.22),
          colorkey = TRUE, screen = list(z = 50, x = -55, y = 2))

#Year 2016-------
year2016 <-subset(byyear, year==2016)
attach(year2016)
year2016$x1=lat.standard=(2*lat-(max(lat)+min(lat)))/(max(lat)-min(lat))
year2016$x2=lon.standard=(2*lon-(max(lon)+min(lon)))/(max(lon)-min(lon))

year2016$X2=year2016$x2^2
year2016$X1=year2016$x1^2
year2016$X13=year2016$x1^3
year2016$X23=year2016$x2^3
year2016$X2x1=year2016$X2*year2016$x1
year2016$X1x2=year2016$X1*year2016$x2

rain2016mean=lm(anom~1, year2016)
anova(rain2016mean)
summary(rain2016mean)

rain2016linear=lm(anom~x1+x2, year2016)
anova(rain2016linear)
summary(rain2016linear)

rain2016quad=lm(anom~x1*x2+X1+X2, year2016)
anova(rain2016quad)
summary(rain2016quad)

rain2016quadalt=lm(anom~x1*x2+X1+X2+alt, year2016)
anova(rain2016quadalt)
summary(rain2016quadalt)

rain2016cubic=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2, year2016)
anova(rain2016cubic)
summary(rain2016cubic)

rain2016cubicalt=lm(anom~x1*x2+X1+X2+X13+X23+X2x1+X1x2+alt, year2016)
anova(rain2016cubicalt)
summary(rain2016cubicalt)


year2016$fitted1=fitted(rain2016cubic)
year2016$resid1=resid(rain2016cubic) #List of residuals
plot(density(resid(rain2016cubic))) #A density plot
qqnorm(resid(rain2016cubic)) # A quantile normal plot - good for checking normality
qqline(resid(rain2016cubic))
year2016fit=data.frame(cbind(year2016$lon, year2016$lat,year2016$fitted1))
year2016res=data.frame(cbind(year2016$lon, year2016$lat,year2016$resid1))

anova(rain2016mean,rain2016linear, rain2016quad, rain2016cubic, rain2016cubicalt)

library(lattice)
w16=wireframe(X3~X1*X2,data=year2016fit,
              xlab = "Longitude", ylab = "Latitude", zlab="yhat", 
              main = "2016", scales = list(arrows = FALSE, distance=2),
              drape = TRUE,aspect = c(70/87, 0.41),
              colorkey = TRUE, screen = list(z = 30, x = -25, y = 16))
w16
c16= contourplot(X3~X1*X2,data=year2016fit,
                 cuts = 20, region = TRUE,
                 xlab = "Longitude",
                 ylab = "Latitude",
                 main = "2016")
c16 

attach(year2016fit)
X1=X2=seq(-1,1,.1)
model=function(X1,X2){16.393-7.824*X1 -168.389*X2-60.968 *X1^2-28.688*X2^2
  +28.689*X1^3+66.616*X2^3-64.216*X2^2*X1+19.485*X1^2*X2-46.516*X1*X2}
model(X1,X2)

z=outer(X1,X2,model)
contour(X1,X2,z,nlevels=20, main="2016")



wireframe(X3~X1*X2,data=year2016res,
          xlab = "Lon", ylab = "Lat", zlab="ehat", 
          main = "Precip.residuals year 2016", scales = list(arrows = FALSE, distance=2),
          drape = TRUE,aspect = c(70/87, 0.22),
          colorkey = TRUE, screen = list(z = 50, x = -55, y = 2))

## Kriging 2016---- 
coordinates(year2016res)=~X1 +X2
proj4string(year2016res) <- "+proj=utm +zone=37 +datum=WGS84 +units=m"
class(year2016res)
summary(year2016res)

# For working with spatial (and spatio-temporal) data, we use the gstat package, which includes functionality for kriging, among other many things.
library(sp)
library(gstat)
# packages for manipulation & visualization
suppressPackageStartupMessages({
  library(dplyr) # for "glimpse"
  library(ggplot2)
  library(scales) # for "comma"
  library(magrittr)
})


##a variogram could be fit as simply as the following code:
vgm <- variogram(X3~1, year2016res) # calculates sample variogram values 
plot(vgm)
library(gstat)
vgm.fit16<- fit.variogram(vgm, model=vgm(230, "Sph",2, 110)) # fit model
plot(vgm, vgm.fit16)

# Arrange the data for the ggplot2 plot
# add the semivariance values of v2 to v1
Fitted <- data.frame(dist = seq(0.16, max(vgm$dist), length = 1672))
Fitted$Spherical <- variogramLine(vgm.fit16, dist_vector = Fitted$dist)$gamma
#convert the dataframes to a long format
library(reshape2)
Empirical <- melt(vgm, id.vars = "dist", measure.vars = "gamma")
Modeled <- melt(Fitted, id.vars = "dist", measure.vars = "Spherical")
library(ggplot2)
#both variogram on the same plot with different colours
v16=ggplot(Modeled, aes(x = dist, y = value, colour = variable )) + 
  geom_line() +  
  xlab("Distance") +
  ylab("Semivariance") +   ggtitle("2016")+ 
  geom_point(data = Empirical)
v16

library(magrittr)
library(gstat)
library(sp)
rain.oK16<- krige(X3 ~ 1, year2016res, grd,model = vgm.fit16)
x <- krige.cv(X3 ~ 1, year2016res, grd,model = vgm.fit16, nfold=100)
#RMSE (root mean squared error)
sqrt(mean(x$residual^2)) 

library(ggplot2)
library(scales)
k16ps=rain.oK16 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKpred"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k16ps=k16ps+xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2016")
k16ps 

k16sd=rain.oK16 %>% as.data.frame %>%
  ggplot(aes(x=x1, y=x2)) + geom_tile(aes(fill=var1.var)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="blue") +guides(fill=guide_legend(title="RKerror"))+
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() 
k16sd=k16sd+  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("2016")
k16sd

# Some parameters for plotting
par(font.main=2, cex.main=1.5,cex.lab=0.4, cex.sub=1)
# Use the ColorBrewer library for color ramps
library(RColorBrewer)
precip.pal = colorRampPalette(c("yellow", "white", "blue"))

#### Countor map with values : Cokriging
ok16ps=spplot(rain.oK16,  zcol='var1.pred',col.regions=precip.pal, contour=TRUE, scales=list(draw=T),
              col='black', pretty=TRUE,main="2016", 
              labels = list(cex = 0.7), label.style = 'align', margin = TRUE,   width = 2, cex = 1.5,  xlab="Longitude",
              ylab="Latitude")
ok16ps

# Export results to latex----------------------
library(stargazer)
stargazer(rain1998cubic, rain1999cubic, rain2000cubic, rain2001cubic,rain2002cubic)
stargazer(rain2003cubic,rain2004cubic,rain2005cubic,rain2006cubic,rain2007cubic)
stargazer(rain2008cubic,rain2009cubic, rain2010cubic, rain2012cubic, rain2012cubic)
stargazer(rain2013cubic, rain2014cubic, rain2015cubic,rain2016cubic)

# Put ggplot in one page -------
library(gridExtra)
grid.arrange(w98,w99,w00,w01, ncol=2)
grid.arrange(w02, w03, w04,w05, ncol=2)
grid.arrange(w06, w07, w08, w09, ncol=2)
grid.arrange(w10, w11, w12, w13, ncol=2)
grid.arrange(w14, w15, w16, ncol=2)


grid.arrange(c98,c99,c00,c01,ncol=2)
grid.arrange(c02, c03, c04, c05,ncol=2)
grid.arrange(c06, c07, c08, c09,ncol=2)
grid.arrange(c10, c11, c12, c13,ncol=2)
grid.arrange(c14, c15, c16,ncol=2)

grid.arrange(v98,v99,v00,v01, v02, v03, v04, v05, v06, ncol=3)
grid.arrange(v07, v08, v09, v10, v11, v12, v13, v14, v15, v16, ncol=3)

grid.arrange(ok98ps,ok99ps,ok00ps,ok01ps,ncol=2)
grid.arrange(ok02ps,ok03ps,ok04ps,ok05ps,ncol=2)
grid.arrange(ok06ps,ok07ps,ok08ps,ok09ps,ncol=2)
grid.arrange(ok10ps,ok11ps,ok12ps,ok13ps,ncol=2)
grid.arrange(ok14ps,ok15ps,ok16ps, ncol=2)


