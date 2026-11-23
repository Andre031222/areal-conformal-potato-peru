# P1 · Definiciones comunes (plan de análisis aprobado el 16-nov-2026).
`%||%` <- function(a, b) if (is.null(a)) b else a

# Rutas relativas: todos los scripts se corren desde la raíz del repositorio.
DATOS <- "data"
RESULTADOS <- "04-resultados"
dir.create(RESULTADOS, showWarnings = FALSE, recursive = TRUE)

SEMILLA <- 20261116
K_BLOQUES <- 20
CRS_METRICO <- 32718          # UTM 18S
SIERRA_M <- 2300

# Variables de M1 (sección 2 del plan): clima con término cuadrático, tierra y controles.
VARIABLES_M1 <- c("tmin", "tmin2", "ppt", "ppt2", "helada", "spei_min",
                  "pendiente", "carbono", "arcilla", "riego", "viaje")
VARIABLES_FENO <- c("ndvi_max", "ndvi_amp", "verde_meses")
VARIABLES_M3_EXTRA <- c("altitud", "ph", "deficit")   # solo en el bosque aleatorio (D2, D3)

leer_distritos <- function() {
  t <- read.csv(file.path(DATOS, "tabla_distrital_2012.csv"), colClasses = c(ubigeo = "character"))
  g <- sf::st_read(file.path(DATOS, "limites", "distritos_2012.gpkg"), quiet = TRUE)
  d <- merge(g, t, by = "ubigeo")
  d <- d[order(d$ubigeo), ]
  d$k <- d$registros_papa
  d$n <- d$registros_cultivo
  d$y <- d$y_tau1
  d
}

# Variables crudas del distrito, con nombres cortos.
crudas <- function(d) {
  data.frame(
    tmin = d$tmin_c_ajuste, ppt = d$ppt_mm_ajuste, helada = d$meses_helada_ajuste,
    spei_min = d$spei3_min_ajuste, pendiente = d$pendiente_grados,
    carbono = d$carbono_org_gkg, arcilla = d$arcilla_pct,
    riego = d$fraccion_riego, viaje = log1p(d$viaje_min),
    ndvi_max = d$feno_ndvi_max, ndvi_amp = d$feno_ndvi_amplitud, verde_meses = d$feno_duracion_verde_meses,
    altitud = d$altitud_media_m, ph = d$ph, deficit = d$deficit_mm_ajuste,
    provincia = d$provincia
  )
}

# Imputa con la mediana de la provincia (o nacional) y estandariza, ambas cosas con los
# datos de `entrena` (dentro de cada pliegue, para no filtrar información del bloque de prueba).
preparar_predictores <- function(d, entrena) {
  x <- crudas(d); xe <- crudas(entrena)
  num <- setdiff(names(x), "provincia")
  for (v in num) {
    med_prov <- tapply(xe[[v]], xe$provincia, median, na.rm = TRUE)
    falta <- is.na(x[[v]])
    x[[v]][falta] <- med_prov[x$provincia[falta]]
    x[[v]][is.na(x[[v]])] <- median(xe[[v]], na.rm = TRUE)
    xe[[v]][is.na(xe[[v]])] <- median(xe[[v]], na.rm = TRUE)
    x[[v]] <- (x[[v]] - mean(xe[[v]])) / sd(xe[[v]])
  }
  x$tmin2 <- x$tmin^2
  x$ppt2 <- x$ppt^2
  x$provincia <- NULL
  x
}
