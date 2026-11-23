# Relieve sombreado para los mapas (figuras 1 y 6): mosaico del Copernicus GLO-90 (data/dem_glo90,
# descargado con data-raw/descargar_suelo_relieve.py) agregado a unos 1,8 km y sombreado con
# iluminación desde el noroeste (azimut 315°, altura 45°), recortado al Perú.
# Salidas: data/relieve_peru.tif (sombreado, 0-1) y data/elevacion_peru.tif (m, misma rejilla). Uso: Rscript data-raw/relieve_sombreado.R
suppressPackageStartupMessages({library(terra); library(sf)})
teselas <- list.files("data/dem_glo90", pattern = "_DEM\\.tif$", full.names = TRUE)
stopifnot(length(teselas) > 100)
gruesas <- lapply(teselas, function(f) aggregate(rast(f), fact = 20, fun = "mean", na.rm = TRUE))
dem <- do.call(merge, unname(gruesas))
peru <- st_union(st_read("data/limites/distritos_2012.gpkg", quiet = TRUE))
peru <- vect(st_transform(st_sf(geometry = peru), crs(dem)))
dem <- mask(crop(dem, peru), peru)
dem <- project(dem, "EPSG:32718", method = "bilinear")      # métrico, para que la pendiente sea correcta
pend <- terrain(dem, "slope", unit = "radians"); asp <- terrain(dem, "aspect", unit = "radians")
sombra <- shade(pend, asp, angle = 45, direction = 315)
names(sombra) <- "sombreado"
writeRaster(sombra, "data/relieve_peru.tif", overwrite = TRUE, gdal = "COMPRESS=DEFLATE")
names(dem) <- "elevacion"
writeRaster(dem, "data/elevacion_peru.tif", overwrite = TRUE, gdal = "COMPRESS=DEFLATE")
cat("relieve_peru.tif:", dim(sombra), "celdas de", round(res(sombra)[1]), "m\n")
