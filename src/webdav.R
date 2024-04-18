# WebDAV: direct access to the maps in VRT format.
# https://files.isric.org/soilgrids/latest/data/
# SoilGrids files in VRT format are hosted on the WebDAV. They can be easily accessed using your file browser and different programming languages.
# https://www.isric.org/explore/soilgrids/soilgrids-access


# Load libraries, select boundary box and cell size
# Initially you need to load the following libraries and select the bounding box of the area you are interested in:

# Check if required packages geobr, sf, terra and gdalUtilities exist. If not, install them.
if (!require(geobr)) {
  install.packages("geobr")
}
if (!require(sf)) {
  install.packages("sf")
}
if (!require(terra)) {
  install.packages("terra")
}
if (!require(gdalUtilities)) {
  install.packages("gdalUtilities")
}

# To deal with an increasing number of inputs and computation demands, SoilGrids has since 2019 been computed on an equal-area projection. After a thorough comparison, the Homolosine projection was identified as the most efficient in an open source software framework (de Sousa et al. 2019). This projection is fully supported by the PROJ and GDAL libraries; therefore, it can be used with any GIS software. The Homolosine projection is now included in the PROJ database with the code ESRI:54052. The actual Spatial Reference System (SRS) of the SoilGrids maps is composed by the Homolosine projection applied to the WGS84 datum. The European Petroleum Survey Group (EPSG) never issued a code for this projection. However, some programmes like MapServer require any SRS to be associated with such a code. For that reason a pseudo EPSG code was created to refer to the SoilGrids SRS: EPSG:152160. The PROJ string for the Homolosine projection is:
homolosine <- "+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs"

# Load the vector of the target area of interest. As we will work with the entire Brazilian territory, we will use the geobr package to load the country's shapefile. The geobr package provides a series of spatial data sets for Brazil, including the shapefile of the country's boundaries and the boundaries of its states and municipalities. 
target <- geobr::read_country()
# target <- geobr::read_state(code_state = "AL")

# Next we transform the CRS of the area of interest to Homolosine projection and get the bounding box, the region of interest specified by georeference coordinates. The bounding box must be upper <ulx> <uly> <lrx> <lry>, where <ulx> is the X value of the upper left corner, <uly> is the Y value of the upper left corner, <lrx> is the X value of the lower right corner and <lry> is the Y value of the lower right corner. This is how GDAL requires the region of interest to be specified. This is different from the way R-package sf handles bounding boxes, where the order is c(xmin, ymin, xmax, ymax).
# bbox <- sf::st_bbox(sf::st_transform(target, homolosine))[c(1, 4, 3, 2)]
bbox <- sf::st_bbox(target)[c(1, 4, 3, 2)]

# Next we will download data from SoilGrids WebDAV service. But first we need to specify the URL of the WebDAV service with some additional parameters that are used to configure a Virtual File System (VFS) for accessing resources via web protocols, in this case, via cURL. Here's a breakdown of the parameters:
# - `max_retry=3`: This sets the maximum number of retries if the connection fails.
# - `retry_delay=1`: This sets the delay between retries in seconds.
# - `list_dir=no`: This indicates that directory listing is not required.
# - `url=https://files.isric.org/soilgrids/latest/data/`: This is the actual URL of the resource to be accessed.
soilgrids_url <- "https://files.isric.org/soilgrids/latest/data/"
webdav_call <- paste0("/vsicurl?max_retry=3&retry_delay=1&list_dir=no&url=", soilgrids_url)

# It gives direct access to the maps in VRT format. There are two options. First, download  a geotiff in Homolosine projection. The following GDAL command will create a local geotiff in the Homolosine projection. We will use the gdalUtilities::gdal_translate(), but this could be done using GDAL command line tools as well.
# In this example, we will work with the entire Brazilian territory and download data for the Acrisols soil class. The following command will download the VRT for the Acrisols class in the Homolosine projection. The VRT file is a lightweight XML file that contains metadata and pointers to the actual data files. The VRT file can be used as input for other GDAL commands. The gdal_translate() function will download the VRT file and save it as a local file named acrisols.vrt. The function will also clip the data to the bounding box of the area of interest.
# More information about the gdal_translate() function can be found at https://gdal.org/programs/gdal_translate.html
gdalUtilities::gdal_translate(
  src_dataset = paste0(webdav_call, "wrb/Acrisols.vrt"),
  dst_dataset = "data/acrisols.tif",
  # tr = c(250, 250), # target resolution <xres> <yres>
  projwin = bbox, # -projwin <ulx> <uly> <lrx> <lry>
  projwin_srs = "EPSG:4326"
  # projwin_srs = homolosine # -projwin_srs <srs_def>
)

# The second option is to download a geotiff in a different projection. In this case, we will download the VRT for the Acrisols class in the Homolosine projection and then reproject it to the WGS84 projection. The following command will download the VRT for the Acrisols class in the Homolosine projection. The VRT file is a lightweight XML file that contains metadata and pointers to the actual data files. The VRT file can be used as input for other GDAL commands. The gdal_translate() function will download the VRT file and save it as a local file named acrisols.vrt. The function will also clip the data to the bounding box of the area of interest. The following commands describe a workflow to obtain a VRT or a GeoTiff for an area of interest in a projection of your choice. In this example we will use EPSG=4326. To local VRT in homolosine (directly from the webdav connection). The first step is to obtain a VRT for the area of interest in the Homolosine projection. We suggest to use VRT for the intermediate steps to save space and computation times.

gdalUtilities::gdal_translate(
  src_dataset = paste0(soilgrids_url, "wrb/Acrisols.vrt"),
  dst_dataset = "data/acrisols.vrt",
  of = "VRT",
  # tr = c(250, 250),
  projwin = bbox,
  projwin_srs = homolosine
)

# To a VRT in, for example, LatLong
# The following command will generate a VRT in the projection of your choice:
gdalUtilities::gdalwarp("data/acrisols.vrt",
    "data/acrisols_wgs84.vrt", 
    s_srs = homolosine, 
    t_srs = "EPSG:4326", 
    of = "VRT",
    overwrite = TRUE)

# To a final Geotiff
# The following command will generate a Geotiff in the projection of your choice for the area of interest defined above
gdalUtilities::gdal_translate("data/acrisols_wgs84.vrt",
    "data/acrisols_wgs84.tif", 
    co = c("TILED=YES","COMPRESS=DEFLATE","PREDICTOR=2","BIGTIFF=YES"))

prob <- 45
cmd <- paste0(
  'gdal_calc.py -A data/acrisols.tif --outfile=data/acrisols_', prob, '.tif --calc="where(A>',
  prob, ', 1, 0)"'
)
system(cmd)

system('gdal_sieve.py data/acrisols_bin.tif -st 100 -8 data/acrisols_bin_sieve.tif')

system("gdal_polygonize.py data/acrisols_bin_sieve.tif data/acrisols_poly.shp")
