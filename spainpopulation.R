# Visualizing a density population map according to geometrical distributions
# www.overfitting.net
# https://www.overfitting.net/2025/03/visualizando-densidades-de-poblacion.html


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
CRS="+proj=lcc +ellps=WGS84 +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=0 +units=km"  # default
# CRS="+proj=lcc +ellps=WGS84 +lat_1=37 +lat_2=43 +lat_0=40 +lon_0=-3 +units=km"  # better for Peninsula
# CRS="+proj=lcc +ellps=WGS84 +lat_1=30 +lat_2=44 +lat_0=37 +lon_0=-3 +units=km"  # better for Spain incl. Canarias
# CRS="+proj=merc +ellps=WGS84 +datum=WGS84 +units=km"  # Mercator for Spain incl. Canarias

spainpop=project(x=spainpop, y=CRS, threads=TRUE)
spainpop
plot(spainpop)
abline(h=0, v=0)  # centre of reference
RESOLUTION=res(spainpop)[1]  # 0.8104684 km resolution

# RESAMPLE raster to Full HD
# DIMY=1080
DIMX=1920
DIMY=round(DIMX*nrow(spainpop)/ncol(spainpop))  # 1662 px
spainpoprs=rast(nrows=DIMY, ncols=DIMX, extent=ext(spainpop))
spainpoprs=resample(x=spainpop, y=spainpoprs, method='bilinear', threads=TRUE)
plot(spainpoprs)
abline(h=0, v=0)  # centre of reference

# Convert to matrix and save as TIFF
DEM=matrix(as.array(spainpoprs), nrow=nrow(spainpoprs))
hist(DEM, breaks=1000)
DEM[is.na(DEM)]=0
writeTIFF((DEM/max(DEM))^(1/2.2), "spainpop.tif", compression='LZW', bits.per.sample=16)

# Solid map for masking in Photoshop
DEMsolid=DEM
DEMsolid[DEMsolid>0]=1
writeTIFF(DEMsolid, "spainpopsolid.tif", compression='LZW')


###########################################################

# 2. FLAG DISTRIBUTION

population=sum(DEM)  # total population to be split

poprow=rowSums(DEM)
poprowacum=cumsum(poprow)

i=1
while (poprowacum[i] < population/3) i=i+1
isup=i
print(paste0("Row for 1/3 of population: ", isup))  # 355

while (poprowacum[i] < population*2/3) i=i+1
iinf=i
print(paste0("Row for 2/3 of population: ", iinf))  # 581

# Draw flag
img=DEM*0
img[isup:iinf,]=1
writeTIFF(img, "flag.tif", compression='LZW')


###########################################################

# 3. MEDIAN DISTRIBUTION

population=sum(DEM)  # total population to be split

poprow=rowSums(DEM)
poprowacum=cumsum(poprow)
popcol=colSums(DEM)
popcolacum=cumsum(popcol)

# Calculate median
i=1
while (poprowacum[i] < population/2) i=i+1
irowcentre=i
print(paste0("Row for 1/2 of population (median): ", irowcentre))  # 447

i=1
while (popcolacum[i] < population/2) i=i+1
icolcentre=i
print(paste0("Column for 1/2 of population (median): ", icolcentre))  # 1292

# Draw medians
img=DEM*0
img[1:irowcentre,1:icolcentre]=1
img[(irowcentre+1):DIMY,1:icolcentre]=2/3
img[(irowcentre+1):DIMY,(icolcentre+1):DIMX]=1/3
writeTIFF(img, "quadrants.tif", compression='LZW')

# Check that quadrants cross population match: NW=SE and NE=SW
quad=matrix(0, nrow=2, ncol=2)
quad[1,1]=sum(DEM[1:irowcentre, 1:icolcentre])/population  # ~22% (NW)
quad[2,1]=sum(DEM[(irowcentre+1):DIMY, 1:icolcentre])/population  # ~28% (SW)
quad[1,2]=sum(DEM[1:irowcentre, (icolcentre+1):DIMX])/population  # ~28% (NE)
quad[2,2]=sum(DEM[(irowcentre+1):DIMY, (icolcentre+1):DIMX])/population  # ~22% (SE)
write.csv2(quad, "quadrants.csv", quote=FALSE, row.names=FALSE)


###########################################################

# 4. CIRCLE DISTRIBUTION

# Calculate 1/3 and 2/3 radius
r=0
sumpop=0
while (sumpop < population/3) {
    r=r+1
    sumpop=sum(DEM[which( ((row(DEM)-irowcentre))^2 + ((col(DEM)-icolcentre))^2 < r^2 )])
}
rinner=r
print(paste0("Radius for 1/3 of population: ", rinner))  # 251

while (sumpop < population*2/3) {
    r=r+1
    sumpop=sum(DEM[which( ((row(DEM)-irowcentre))^2 + ((col(DEM)-icolcentre))^2 < r^2 )])
}
router=r
print(paste0("Radius for 2/3 of population: ", router))  # 358

# Draw circles
img=DEM*0
img[which( ((row(img)-irowcentre))^2 + ((col(img)-icolcentre))^2 < router^2 )]=0.5
img[which( ((row(img)-irowcentre))^2 + ((col(img)-icolcentre))^2 < rinner^2 )]=1
writeTIFF(img, "circles.tif", compression='LZW')


###########################################################

# 5. LOWEST VS HIGHEST POPULATION DENSITY

# Rank each value in DEM from max to min population density
rankedindices=order(DEM, decreasing=TRUE)
rankedvalues=DEM[rankedindices]
rankedvaluesacum=cumsum(rankedvalues)

PERC=0.8  # where do 80% of population live?
i=1
while (rankedvaluesacum[i] < population*PERC) i=i+1
icut=i
print(paste0("Count for ", PERC*100, "% of population: ", icut))  # 16358 points
print(paste0(PERC*100, "% of the population occupies ",
             round(icut/length(DEMsolid[DEMsolid==1])*100,1),"% of the land"))

# Draw half zones
img=DEMsolid/2
img[rankedindices[1:icut]]=1
writeTIFF(img, "halfs.tif", compression='LZW')



