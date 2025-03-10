# Visualizing a density population map according to geometrical distributions
# www.overfitting.net
# https://www.overfitting.net/


library(terra)  # read GeoTIFF, reprojection, crop and resample
library(tiff)  # save 16-bit TIFF's


###########################################################

# 1. PROCESS GEOTIFF DATA TO GET POPULATION DENSITY AS A MATRIX

# https://data.humdata.org/dataset/worldpop-population-density-for-spain
spainpop=rast("esp_pd_2020_1km.tif")
spainpop
plot(spainpop)
RESOLUTION=res(spainpop)[1]  # 0.008333 degrees resolution

# REPROJECT raster from Longitude Latitude (+proj=longlat)/WGS84
# to Lambert Conic Conformal (+proj=lcc)/WGS84
# by default crs="+proj=lcc +ellps=WGS84 +lat_1=33 +lat_2=45 +lon_0=0"
CRS="+proj=lcc +ellps=WGS84 +lat_1=33 +lat_2=45 +lon_0=0 +units=km"
spainpop=project(x=spainpop, y=CRS, threads=TRUE)
spainpop
plot(spainpop)
abline(v=0)  # Greenwich meridian
RESOLUTION=res(spainpop)[1]  # 0.8104684m grid resolution

# RESAMPLE raster to Full HD
# DIMY=1080
DIMX=1920
DIMY=DIMX*nrow(spainpop)/ncol(spainpop)
spainpoprs=rast(nrows=DIMY, ncols=DIMX, extent=ext(spainpop))
spainpoprs=resample(x=spainpop, y=spainpoprs, method='bilinear', threads=TRUE)
plot(spainpoprs)
abline(v=0)  # Greenwich meridian

# Convert to matrix and save as TIFF
DEM=matrix(as.array(spainpoprs), nrow=nrow(spainpoprs))
hist(DEM, breaks=1000)
DEM[is.na(DEM)]=0
DEM[DEM<0]=0  # let's preserve the ocean basin

# Solid map for masking in Photoshop
DEMsolid=DEM
DEMsolid[DEMsolid>0]=1
writeTIFF(DEMsolid, "spainpopsolid.tif", compression='LZW', bits.per.sample=16)

# Grayscale density map
writeTIFF((DEM/max(DEM))^(1/2.2), "spainpop.tif", compression='LZW', bits.per.sample=16)


###########################################################

# 2. FLAG DISTRIBUTION

population=sum(DEM)  # total population to be splitted into 1/3, 2/3
poprow=rowSums(DEM)
poprowacum=cumsum(poprow)

i=1
while (poprowacum[i]<population/3) i=i+1
isup=i
print(paste0("Row for 1/3 of population: ", isup))  # 355

while (poprowacum[i]<population*2/3) i=i+1
iinf=i
print(paste0("Row for 2/3 of population: ", iinf))  # 581

# Draw flag
img=DEM*0
img[isup:iinf,]=1
writeTIFF(img, "flag.tif", compression='LZW')


###########################################################

# 3. MEDIAN & CIRCLE DISTRIBUTION

population=sum(DEM)  # total population to be splitted into 3 circles
poprow=rowSums(DEM)
poprowacum=cumsum(poprow)
popcol=colSums(DEM)
popcolacum=cumsum(popcol)

# Calculate median
i=1
while (poprowacum[i]<population/2) i=i+1
irowcentre=i
print(paste0("Row for 1/2 of population: ", irowcentre))  # 447

i=1
while (popcolacum[i]<population/2) i=i+1
icolcentre=i
print(paste0("Column for 1/2 of population: ", icolcentre))  # 1292

# Draw medians
img=DEM*0
img[1:irowcentre,1:icolcentre]=1
img[(irowcentre+1):DIMY,1:icolcentre]=2/3
img[(irowcentre+1):DIMY,(icolcentre+1):DIMX]=1/3
writeTIFF(img, "quarters.tif", compression='LZW')


# Calculate 1/3 and 2/3 radius
r=0
sumpop=0
while (sumpop<population/3) {
    r=r+1
    sumpop=sum(DEM[which( ((row(tmp)-irowcentre))^2 + ((col(tmp)-icolcentre))^2 < r^2 )])
}
rinner=r
print(paste0("Radius for 1/3 of population: ", rinner))  # 251

while (sumpop<population*2/3) {
    r=r+1
    sumpop=sum(DEM[which( ((row(tmp)-irowcentre))^2 + ((col(tmp)-icolcentre))^2 < r^2 )])
}
router=r
print(paste0("Radius for 2/3 of population: ", router))  # 358

# Draw circles
img=DEM*0
img[which( ((row(img)-irowcentre))^2 + ((col(img)-icolcentre))^2 < router^2 )]=0.5
img[which( ((row(img)-irowcentre))^2 + ((col(img)-icolcentre))^2 < rinner^2 )]=1
writeTIFF(img, "circles.tif", compression='LZW')





