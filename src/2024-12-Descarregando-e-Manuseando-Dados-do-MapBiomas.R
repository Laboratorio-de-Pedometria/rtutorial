# Downloading and Handling MapBiomas Data with R and GDAL/OGR
# Alessandro Samuel-Rosa & Taciara Zborowski Horst
# 2024-12-20

# Processing municipality: Santana - AP
# muni_geom size: 0 Mb
# pedology_muni size: 2.9 Mb
# pedology_muni size: 0.2 Mb
# mapbiomas_muni size: 0 Mb
# Error : [] nrow dataframe does not match nrow geometry
# In addition: Warning messages:
# 1: [SpatVector from sf] not all geometries were transferred, use svc for a geometry collection 
# 2: [SpatVector from sf] not all geometries were transferred, use svc for a geometry collection 
# 3: [SpatVector from sf] not all geometries were transferred, use svc for a geometry collection 
# Error: [as,sf] coercion failed. You can try coercing via a Spatial* (sp) class
# > Sys.time() - t0
# Time difference of 58.59917 mins

# Summary
# In this exercise, we will identify the dominant soil class in areas with agriculture for each
# Brazilian municipality.

# Create necessary directories to store data and results
# If the directories already exist, the function will return a warning.
# - data: to store downloaded data
# - res: to store results
# - tmp: to store temporary files
dir.create("data", showWarnings = TRUE)
dir.create("res", showWarnings = TRUE)
dir.create("tmp", showWarnings = TRUE)

# Load necessary libraries, installing them if they are not available
if (!require(sf)) {
  install.packages("sf")
  library(sf)
}
if (!require(geobr)) {
  install.packages("geobr")
  library(geobr)
}
if (!require(terra)) {
  install.packages("terra")
  library(terra)
}

# We will be downloading large files, so we need to increase the timeout for the download
# The default timeout is 60 seconds. We will increase it to 3600 seconds (60 minutes).
options(timeout = 3600)

# Helper function to print object name and size
print_object_size <- function(object) {
  object_name <- deparse(substitute(object))
  size <- format(object.size(object), units = "MB")
  cat(paste0(object_name, " size: ", size, "\n"))
}

# Brazilian municipalities #########################################################################
# We start by downloading the vector of the Brazilian municipalities. The data is available at a
# scale of 1:250,000 and uses the Geodetic reference system "SIRGAS2000" and CRS(4674).
# We will use the geobr package to download the data and save it in the GeoPackage format.
# The geobr package is a wrapper for the IBGE API, which provides access to spatial data from
# the Brazilian Institute of Geography and Statistics (IBGE). While downloading, we will set the
# argument simplified = FALSE to keep the original spatial resolution.
# The GeoPackage format is an open, standards-based, platform-independent, portable,
# self-describing, compact format for transferring geospatial information. It is based on
# SQLite, which makes it easy to read and write using R.

# Attention: This is a somewhat large dataset (~260 MB)!
# We first check if the data is already available in the data directory. If not, we download it and
# save it in the data directory.
# We transform the data to the WGS84 CRS (EPSG:4326) to match the other datasets we will download.
output_file <- "data/municipalities.gpkg"
if (!file.exists(output_file)) {
  t0 <- Sys.time()
  muni <- geobr::read_municipality(code_muni = "all", year = 2020, simplified = FALSE)
  muni <- sf::st_transform(muni, crs = 4326)
  sf::st_write(muni, output_file)
  Sys.time() - t0
  rm(muni, t0)
} else {
  cat("The municipalities data is already available.\n")
}

# Soil Map of Brazil ##############################################################################
# The second dataset we will download is the vector of the Brazilian pedology map. The data is 
# available at a scale of 1:250,000 and uses the Geodetic reference system "SIRGAS2000" and 
# CRS(4674).
# We will download the data from the IBGE GeoServer. The data is available in the GeoJSON format.
# https://geoservicos.ibge.gov.br/geoserver/ows

# Attention: This is a large dataset (~1.2 GB)! Your internet connection may be slow, and the
# download may take a long time. Besides, the memory of your computer may be insufficient to handle
# this dataset in R.
# We first check if it is already available in the data directory. If not, we download it and save
# it in the data directory.
input_url <- paste0(
  "https://geoservicos.ibge.gov.br/geoserver/ows?service=WFS&version=1.0.0",
  "&request=GetFeature&typeName=BDIA:pedo_area&outputFormat=application/json"
)
output_file <- "data/pedology.json"
if (!file.exists(gsub("json", "gpkg", output_file))) {
  t0 <- Sys.time()
  download.file(url = input_url, destfile = output_file)
  cmd <- paste("ogr2ogr -t_srs EPSG:4326 -makevalid -of GPKG data/pedology.gpkg", output_file)
  system(cmd)
  Sys.time() - t0
} else {
  cat("The pedology data is already available.\n")
}

# Land Use and Land Cover Map of Brazil ##########################################################
# The third and last dataset we will download is the raster of the Land Use and Land Cover (LULC)
# map of Brazil. The data is available at a scale of 30 meters and uses the Geodetic reference
# system "WGS84" and CRS(4326). The data is from the MapBiomas project, which provides annual
# land use and land cover maps for Brazil.
# We will download the data from the MapBiomas API. The data is available in the GeoTIFF format.
# https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_9/lclu/coverage/

# Attention: This is a large dataset (~1.1 GB). Your internet connection may be slow, and the
# download may take a long time. Besides, the memory of your computer may be insufficient to handle
# this dataset in R.
# We first check if it is already available in the data() directory. If not, we download it and
# save it in the data directory.
input_url <- paste0(
  "https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/",
  "collection_9/lclu/coverage/brasil_coverage_2023.tif"
)
output_file <- "data/mapbiomas2023.tif"
if (!file.exists(output_file)) {
  t0 <- Sys.time()
  download.file(url = input_url, destfile = output_file)
  Sys.time() - t0
} else {
  cat("The MapBiomas data is already available.\n")
}

# Data analysis ####################################################################################
# We will identify, for each municipality, the dominant soil class in areas with agriculture.
# By agriculture, we mean the classes related to agriculture in the MapBiomas data.

# Read the municipalities
muni <- sf::st_read("data/municipalities.gpkg")
print_object_size(muni)
# Create columns to store the results
# - soil_class_first: the dominant soil class in the municipality
# - soil_prop_first: the proportion of the dominant soil class in the municipality
# - soil_class_second: the second dominant soil class in the municipality
# - soil_prop_second: the proportion of the second dominant soil class in the municipality
# - agri_area_ha: the area of agriculture in hectares
muni["soil_class_first"] <- NA_character_
muni["soil_prop_first"] <- NA_real_
muni["soil_class_second"] <- NA_character_
muni["soil_prop_second"] <- NA_real_
muni["agri_area_ha"] <- NA_real_
print(muni)

# Read the MapBiomas data
# When we read a raster with the terra package, it does not load the cell (pixel) values into
# memory (RAM). It only reads the parameters that describe the geometry of the raster, such as
# the number of rows and columns and the coordinate reference system. The actual values will be
# read when needed. This saves memory and allows us to work with large rasters.
mapbiomas <- terra::rast("data/mapbiomas2023.tif")
print(mapbiomas)
agri_classes <- c(
  18, # Agriculture
  19, # Temporary Crop
  39, # Soybean
  20, # Sugar cane
  40, # Rice
  62, # Cotton (beta)
  41, # Other Temporary Crops
  36, # Perennial Crop
  46, # Coffee
  47, # Citrus
  35, # Palm Oil
  48  # Other Perennial Crops
)
reclass_matrix <- matrix(cbind(agri_classes, rep(1, length(agri_classes))), ncol = 2)

# Loop over municipalities
# There are 5,570 municipalities in Brazil. This loop may take a long time to complete.
t0 <- Sys.time()
for (i in 1:nrow(muni)) {
  # Print municipality name
  name_muni <- muni[i, "name_muni"][[1]]
  abbrev_state <- muni[i, "abbrev_state"][[1]]
  cat(paste0("\nProcessing municipality: ", name_muni, " - ", abbrev_state, "\n"))
  
  # Get the geometry of the municipality
  muni_geom <- muni[i, "geom"]
  print_object_size(muni_geom)

  # Read the pedology data to the bounding box of the municipality
  wkt <- sf::st_as_text(sf::st_geometry(muni_geom))
  pedology_muni <- sf::st_read("data/pedology.gpkg", wkt_filter = wkt, quiet = TRUE)
  print_object_size(pedology_muni)
  # Set the aggregation function to "constant" to avoid aggregation of the data
  # A constant attribute means that the attribute values do not change within the spatial extent
  # of each geometry, in this case, the soil class.
  sf::st_agr(pedology_muni) <- "constant"  
  
  # Clip the pedology data to the boundaries of the selected municipality
  pedology_muni <- sf::st_intersection(pedology_muni, muni_geom)
  print_object_size(pedology_muni)
  # x11()
  # plot(pedology_muni["id"], main = name_muni)

  # Crop the MapBiomas raster data to the boundaries of the pedology data
  # This operation will require up to 12 GB of memory to process the municipality of Altamira (PA).
  mapbiomas_muni <- terra::crop(mapbiomas, pedology_muni)
  print_object_size(mapbiomas_muni)
  # plot(mapbiomas_muni, main = name_muni)
  # plot(pedology_muni["legenda"], add = TRUE)

  # Clip (mask) the MapBiomas data to the boundary of the pedology data
  # This operation will require up to 16 GB of memory to process the municipality of Altamira (PA).
  mapbiomas_muni <- terra::mask(mapbiomas_muni, pedology_muni)
  # plot(mapbiomas_muni, main = name_muni)
  # plot(pedology_muni["legenda"], add = TRUE)
  # print_object_size(mapbiomas_muni)

  # Create a mask for agriculture classes using terra::classify
  # This operation will require up to 14 GB of memory to process the municipality of Altamira (PA).
  mapbiomas_muni <- terra::classify(mapbiomas_muni, reclass_matrix, others = NA)
  print_object_size(mapbiomas_muni)
  # If all values are NA, skip the municipality
  # This operation will require up to 14 GB of memory to process the municipality of Altamira (PA).
  muni_all_na <- terra::global(mapbiomas_muni, "isNA") == terra::ncell(mapbiomas_muni)
  if (muni_all_na) {
    muni[i, "soil_class_first"] <- NA_character_
    muni[i, "soil_prop_first"] <- NA_real_
    muni[i, "soil_class_second"] <- NA_character_
    muni[i, "soil_prop_second"] <- NA_real_
    muni[i, "agri_area_ha"] <- NA_real_
  } else {
    # Rasterize the pedology data, taking the mapbiomas_muni as a template. This will ensure that
    # the pedology data has the same extent, resolution, and alignment as the mapbiomas data.
    # We will use the "legenda" field to rasterize the data.
    # This operation will require up to 10 GB of memory to process the municipality of Altamira (PA)
    pedology_muni <- terra::rasterize(pedology_muni, mapbiomas_muni, field = "legenda")
    print_object_size(pedology_muni)
    # Mask the pedology raster to the mapbiomas data. This will ensure that the pedology data
    # is available only for the agriculture areas.
    # This operation will require up to 16 GB of memory to process the municipality of Altamira (PA)
    pedology_muni <- terra::mask(pedology_muni, mapbiomas_muni)
    print_object_size(pedology_muni)
    # plot(pedology_muni)
    # Get the frequency of the pedology data and sort by 'count'  
    # This operation will require up to 14 GB of memory to process the municipality of Altamira (PA).
    pedology_table <- terra::freq(pedology_muni)
    pedology_table <- pedology_table[order(-pedology_table[, "count"]), ]
    # Record the two dominant soil classes in the municipality data
    muni[i, "soil_class_first"] <- pedology_table[1, "value"]
    muni[i, "soil_class_second"] <- pedology_table[2, "value"]
    # Record the proportion of the dominant soil class in the municipality data. The proportion is
    # calculated as the percentage of the pixels of the dominant soil class in relation
    muni[i, "soil_prop_first"] <- round(
      pedology_table[1, "count"] / sum(pedology_table[, "count"]) * 100)
    muni[i, "soil_prop_second"] <- round(
      pedology_table[2, "count"] / sum(pedology_table[, "count"]) * 100)
    # Record the area of agriculture in hectares
    # This operation will require up to 21 GB of memory to process the municipality of Altamira (PA).
    muni[i, "agri_area_ha"] <- round(
      sum(terra::cellSize(mapbiomas_muni, unit = "ha")[] * mapbiomas_muni[], na.rm = TRUE)
    )
    # Print the dominant soil class
    soil_class <- muni[i, "soil_class_first"][[1]]
    soil_prop <- muni[i, "soil_prop_first"][[1]]
    cat(paste0("Dominant soil class: ", soil_class, " - ", soil_prop, "%\n"))
  }
}
Sys.time() - t0

# Save the results
output_file <- "~/ownCloud/geoBioAS/res/2024-12-21/municipalities_soil_agri.gpkg"
sf::st_write(muni, output_file, driver = "GPKG")
