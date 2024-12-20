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
if (!require(terra)) {
  install.packages("terra")
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
  t0 <- Sys.time()
  muni <- geobr::read_municipality(code_muni = "all", year = 2020, simplified = FALSE)
  muni <- sf::st_transform(muni, crs = 4326)
  sf::st_write(muni, output_file)
  Sys.time() - t0
} else {
  cat("The municipalities data is already available.\n")
}

# Soil Map of Brazil ##############################################################################

# This is a vector of the Brazilian pedology map. The data is available at a scale of 1:250,000
# and uses the Geodetic reference system "WGS84" and CRS(4326).

# We will download the data from the IBGE GeoServer. The data is available in the GeoJSON format.
# https://geoservicos.ibge.gov.br/geoserver/ows

# This is a large dataset (~1.2 GB). So we first check if it is already available in the data
# directory. If not, we download it and save it in the data directory.
input_url <- paste0(
  "https://geoservicos.ibge.gov.br/geoserver/ows?service=WFS&version=1.0.0",
  "&request=GetFeature&typeName=BDIA:pedo_area&outputFormat=application/json"
)
output_file <- "data/pedology.json"
if (!file.exists(output_file)) {
  t0 <- Sys.time()
  download.file(url = input_url, destfile = output_file)
  pedology <- sf::st_read("data/pedology.json")
  pedology <- sf::st_transform(pedology, crs = 4326)
  pedology <- sf::st_make_valid(pedology)
  sf::st_write(pedology, "data/pedology.gpkg")
  Sys.time() - t0
} else {
  cat("The pedology data is already available.\n")
}

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
  t0 <- Sys.time()
  download.file(url = input_url, destfile = output_file)
  Sys.time() - t0
} else {
  cat("The MapBiomas data is already available.\n")
}

# Data analysis ####################################################################################
# We need to identify, for each municipality, the dominant soil class occurring in areas with
# agriculture.
# The areas of agriculture are defined by the MapBiomas data using the following codes:
# Agriculture 18
# Temporary Crop 19
# Soybean 39
# Sugar cane 20
# Rice 40
# Cotton (beta) 62
# Other Temporary Crops 41
# Perennial Crop 36
# Coffee 46
# Citrus 47
# Palm Oil 35
# Other Perennial Crops 48

# Read the municipalities
muni <- sf::st_read("data/municipalities.gpkg")
muni["soil_class"] <- NA_character_
muni["soil_prop"] <- NA_real_
muni["agri_area_ha"] <- NA_real_
print(muni)

# Read the MapBiomas data
mapbiomas <- terra::rast("data/mapbiomas2023.tif")
print(mapbiomas)
agriculture_classes <- c(18, 19, 39, 20, 40, 62, 41, 36, 46, 47, 35, 48)

# Loop over municipalities
for (i in 1:2) {
  # Print municipality name
  name_muni <- muni[i, "name_muni"][[1]]
  abbrev_state <- muni[i, "abbrev_state"][[1]]
  cat(paste0("Processing municipality: ", name_muni, " - ", abbrev_state, "\n"))
  # Get the geometry of the municipality
  muni_geom <- muni[i, "geom"]

  # Clip the pedology data to the municipality
  wkt <- sf::st_as_text(sf::st_geometry(muni_geom[1, ]))
  pedology <- sf::st_read("data/pedology.gpkg", wkt_filter = wkt, quiet = TRUE)
  pedology <- sf::st_intersection(pedology, muni_geom)

  # Crop the MapBiomas data (SpatRaster) to the pedology
  mapbiomas_clip <- terra::crop(mapbiomas, pedology)

  # Mask the MapBiomas data to the pedology
  mapbiomas_clip <- terra::mask(mapbiomas_clip, pedology)

  # Create a mask for agriculture classes: if the value is in the agriculture classes, set it
  # to 1, otherwise set it to 0
  mapbiomas_clip[] <- ifelse(mapbiomas_clip[] %in% agriculture_classes, 1, NA_real_)
  # If is all NA, skip the municipality
  if (all(is.na(mapbiomas_clip[]))) {
    muni[i, "soil_class"] <- NA_character_
  } else {
    # rasterize the pedology data, taking the mapbiomas_clip as a template
    pedology_raster <- terra::rasterize(pedology, mapbiomas_clip, field = "legenda")
    pedology_raster <- terra::mask(pedology_raster, mapbiomas_clip)
    pedology_table <- table(pedology_raster[])
    pedology_levels <- levels(pedology_raster)[[1]]
    pedology_levels_idx <- pedology_levels$ID == names(which.max(pedology_table))
    muni[i, "soil_class"] <- pedology_levels[pedology_levels_idx, "legenda"]
    muni[i, "soil_prop"] <- round(max(pedology_table) / sum(pedology_table) * 100)
    muni[i, "agri_area_ha"] <- round(
      sum(terra::cellSize(mapbiomas_clip, unit = "ha")[] * mapbiomas_clip[], na.rm = TRUE)
    )
    # Print the dominant soil class
    soil_class <- muni[i, "soil_class"][[1]]
    soil_prop <- muni[i, "soil_prop"][[1]]
    cat(paste0("Dominant soil class: ", soil_class, " - ", soil_prop, "%\n"))
  }
}

# Save the results
muni






length(terra::cellSize(mapbiomas_clip)[])
length(na.exclude(mapbiomas_clip[]))
