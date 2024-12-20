# Downloading and Handling MapBiomas Data with R and GDAL
# Alessandro Samuel-Rosa
# 2024-12-20

# Create necessary directories to store data and results
# If the directories already exist, the function will return a warning.
# - data: to store downloaded data
# - res: to store results
dir.create("data", showWarnings = TRUE)
dir.create("res", showWarnings = TRUE)

# Load necessary libraries
if (!require(geobr)) {
  install.packages("geobr")
}
if (!require(sf)) {
  install.packages("sf")
}
if (!require(gdalUtilities)) {
  install.packages("gdalUtilities")
}

# We will be downloading large files, so we need to increase the timeout
options(timeout = 3000)  # Increase timeout to 300 seconds

# Brazilian municipalities #########################################################################

# This is a vector of the Brazilian municipalities. The data is available at a scale of 1:250,000
# and uses the Geodetic reference system "SIRGAS2000" and CRS(4674).

# We will use the geobr package to download the data and save it in the GeoPackage format.
# The geobr package is a wrapper for the IBGE API, which provides access to spatial data from
# the Brazilian Institute of Geography and Statistics (IBGE). 
#   * We set the simplified argument to FALSE to keep the original spatial resolution.

# The GeoPackage format is an open, standards-based, platform-independent, portable,
# self-describing, compact format for transferring geospatial information.

# This is a large dataset (~260 MB). So we first check if it is already available in the data
# directory. If not, we download it and save it in the data directory.
output_file <- "data/municipalities.gpkg"
if (!file.exists(output_file)) {
  muni <- geobr::read_municipality(code_muni = "all", year = 2020, simplified = FALSE)
  muni <- sf::st_transform(muni, crs = 4326)
  sf::st_write(muni, output_file)
}

# Soil Map of Brazil ##############################################################################

# This is a vector of the Brazilian pedology map. The data is available at a scale of 1:250,000
# and uses the Geodetic reference system "WGS84" and CRS(4326).

# We will download the data from the IBGE GeoServer. The data is available in the GeoJSON format.
# https://geoservicos.ibge.gov.br/geoserver/ows

# This is a large dataset (~1.5 GB). So we first check if it is already available in the data
# directory. If not, we download it and save it in the data directory.
input_url <- paste0(
  "https://geoservicos.ibge.gov.br/geoserver/ows?service=WFS&version=1.0.0",
  "&request=GetFeature&typeName=BDIA:pedo_area&outputFormat=application/json"
)
output_file <- "data/pedology.json"
if (!file.exists(output_file)) {
  download.file(url = input_url, destfile = output_file)
}

# get file info using gdal/ogr
system("ogrinfo data/pedology.json")


# Land Use and Land Cover Map of Brazil ##########################################################

# This is a raster of the Land Use and Land Cover (LULC) map of Brazil. The data is available at a
# scale of 30 meters and uses the Geodetic reference system "WGS84" and CRS(4326).

# We will download the data from the MapBiomas API. The data is available in the GeoTIFF format.
# https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_9/lclu/coverage/

# This is a large dataset (~1.1 GB). So we first check if it is already available in the data
# directory. If not, we download it and save it in the data directory.
input_url <- paste0(
  "https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/",
  "collection_9/lclu/coverage/brasil_coverage_2023.tif"
)
output_file <- "data/mapbiomas2023.tif"
if (!file.exists(output_file)) {
  download.file(url = input_url, destfile = output_file)
}

# Data processing ##################################################################################

# Now we will create a soybean mask, that is a binary map where the value 1 represents the presence
# of soybean crops and the value 0 represents the absence of soybean crops.














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
