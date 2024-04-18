# Downloading and Handling SoilGrids Data with R and GDAL
# Alessandro Samuel-Rosa
# 2024-04-18 

# Original tutorial
# https://git.wur.nl/isric/soilgrids/soilgrids.notebooks/-/blob/master/markdown/webdav_from_R.md

# Load libraries
if (!require(geobr)) {
  install.packages("geobr")
}
if (!require(sf)) {
  install.packages("sf")
}
if (!require(gdalUtilities)) {
  install.packages("gdalUtilities")
}
if (!require(raster)) {
  install.packages("raster")
}

# Load vector of the area of interest
target <- geobr::read_state(code_state = "AC")

# Compute the bounding box of the area of interest
bbox <- sf::st_bbox(target)
bbox <- bbox[c("xmin", "ymax", "xmax", "ymin")]

# Configure WebDAV access
# - max_retry: This sets the maximum number of retries if the connection fails.
# - retry_delay: This sets the delay between retries in seconds.
# - list_dir: This indicates that directory listing is not required.
# - url: This is the actual URL of the resource to be accessed.
max_retry <- 3
retry_delay <- 1
list_dir <- "no"
url <- "https://files.isric.org/soilgrids/latest/data/"
webdav_call <- paste0(
  "/vsicurl?max_retry=", max_retry,
  "&retry_delay=", retry_delay,
  "&list_dir=", list_dir,
  "&url=", url
)

# Download WRB maps of soil classification
map_name <- "wrb/Acrisols.vrt"
output_file <- "data/acrisols.tif"
gdalUtilities::gdal_translate(
  src_dataset = paste0(webdav_call, map_name),
  dst_dataset = output_file,
  projwin = bbox,
  projwin_srs = "EPSG:4326"
)

# Read file and plot map
acrisols <- raster::raster(output_file)
png("res/acrisols.png")
image(acrisols, main = "Acrisols")
plot(target, add = TRUE, col = "transparent", lwd = 2)
dev.off()

# Create a binary map using gdal_calc.py
# The input file is the output file from the previous step.
input_file <- output_file
prob <- 40
calc <- paste0("where(A>", prob, ", 1, 0)")
output_file <- paste0("data/acrisols_bin", prob, ".tif")
cmd <- paste0(
  "gdal_calc.py -A ",input_file, " --outfile=", output_file, " --calc=\"", calc, "\""
)
system(cmd)

# Compute minimum legible area
mla <- ((0.6 * 250000) / 100000)^2

# Compute threshold for sieve filter
threshold <- round(((sqrt(mla) * 1000) / 250)^2)

# Filter binary map
# The input file is the output file from the previous step.
# -st: the threshold to remove isolated pixels
# -8: the 8-connectedness of the sieve filter
input_file <- output_file
threshold <- ceiling(((sqrt(2.5) * 1000) / 250)^2)
connectedness <- 8
output_file <- paste0(
  gsub(".tif", "", input_file), "_sieve", threshold, ".tif"
)
cmd <- paste0(
  "gdal_sieve.py ", input_file, " -st ", threshold, " -", connectedness, " ",
  output_file
)
system(cmd)

# Read file and plot map
acrisols <- raster::raster(output_file)
png("res/acrisols_bin40_sieve41.png")
image(acrisols, main = "Acrisols", col = c("white", "orange"))
plot(target, add = TRUE, col = "transparent", lwd = 2)
dev.off()

# Transform the raster data to polygons
# The input file is the output file from the previous step.
input_file <- output_file
output_file <- gsub(".tif", ".shp", input_file)
cmd <- paste0("gdal_polygonize.py ", input_file, " ", output_file)
system(cmd)








# To deal with an increasing number of inputs and computation demands, SoilGrids has since 2019 been computed on an equal-area projection. After a thorough comparison, the Homolosine projection was identified as the most efficient in an open source software framework (de Sousa et al. 2019). This projection is fully supported by the PROJ and GDAL libraries; therefore, it can be used with any GIS software. The Homolosine projection is now included in the PROJ database with the code ESRI:54052. The actual Spatial Reference System (SRS) of the SoilGrids maps is composed by the Homolosine projection applied to the WGS84 datum. The European Petroleum Survey Group (EPSG) never issued a code for this projection. However, some programmes like MapServer require any SRS to be associated with such a code. For that reason a pseudo EPSG code was created to refer to the SoilGrids SRS: EPSG:152160. The PROJ string for the Homolosine projection is:
homolosine <- "+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs"

# Next we transform the CRS of the area of interest to Homolosine projection and get the bounding box, the region of interest specified by georeference coordinates. 
bbox <- sf::st_bbox(sf::st_transform(target, homolosine))[c(1, 4, 3, 2)]

# Download  a geotiff in Homolosine projection. The following GDAL command will create a local geotiff in the Homolosine projection. We will use the gdalUtilities::gdal_translate(), but this could be done using GDAL command line tools as well.
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
