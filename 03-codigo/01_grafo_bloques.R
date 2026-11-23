# P1 · Paso 1: grafo de vecindad y bloques espaciales (plan de análisis, secciones 3 y 4).
#
# - Vecindad: distritos que comparten un tramo de límite (no basta un punto: queen = FALSE).
#   Los distritos sin vecinos (islas) se unen al distrito más cercano.
# - Tamaño de bloque: rango de autocorrelación de los residuos de deviance de M2
#   (binomial sin efecto espacial), estimado con blockCV::cv_spatial_autocor.
# - k = 20 bloques con blockCV::cv_spatial; semilla 20261116.
#
# Salidas en 04-resultados/: grafo.adj (formato INLA), pliegues.csv, bloques.txt (resumen).
suppressPackageStartupMessages({library(sf); library(spdep); library(blockCV)})
source("03-codigo/00_comun.R")   # correr desde la carpeta del paper

d <- leer_distritos()
cat("distritos:", nrow(d), "\n")

# --- Grafo -----------------------------------------------------------------
vecinos <- poly2nb(d, queen = FALSE, snap = 1e-6)
islas <- which(card(vecinos) == 0)
if (length(islas)) {
  centros <- st_coordinates(st_point_on_surface(st_geometry(st_transform(d, CRS_METRICO))))
  for (i in islas) {
    dist <- sqrt(colSums((t(centros) - centros[i, ])^2)); dist[i] <- Inf
    j <- which.min(dist)
    vecinos[[i]] <- sort(unique(c(setdiff(vecinos[[i]], 0L), j)))
    vecinos[[j]] <- sort(unique(c(setdiff(vecinos[[j]], 0L), i)))
  }
}
comp <- n.comp.nb(vecinos)
cat("islas unidas al más cercano:", length(islas), "| componentes conexos:", comp$nc, "\n")
stopifnot(comp$nc == 1)
nb2INLA(file.path(RESULTADOS, "grafo.adj"), vecinos)
saveRDS(vecinos, file.path(RESULTADOS, "vecinos.rds"))

# --- Residuos de M2 y rango de autocorrelación -----------------------------
X <- preparar_predictores(d, d)          # estandarización con todos los datos (solo para el rango)
m2 <- glm(cbind(k, n - k) ~ ., family = binomial, data = cbind(k = d$k, n = d$n, X[, VARIABLES_M1]))
d$residuo <- residuals(m2, type = "deviance")
puntos <- st_sf(residuo = d$residuo,
                geometry = st_point_on_surface(st_geometry(st_transform(d, CRS_METRICO))))
auto <- cv_spatial_autocor(x = puntos, column = "residuo", plot = FALSE, progress = FALSE)
rango_blockcv <- auto$range

# Desviación 1 (16-nov-2026): el ajuste automático de blockCV dio 25,5 km, pero el
# correlograma empírico de los residuos sigue en r > 0,07 hasta 100 km. Se usa el rango
# práctico (3a) de un modelo exponencial r(h) = c exp(-h/a) ajustado al correlograma,
# ponderando cada intervalo por su número de pares.
z <- (d$residuo - mean(d$residuo)) / sd(d$residuo)
D <- as.matrix(dist(st_coordinates(puntos))) / 1000
ut <- upper.tri(D)
cortes <- seq(0, 500, by = 10)
banda <- cut(D[ut], cortes)
correlograma <- data.frame(
  h = head(cortes, -1) + 5,
  r = as.numeric(tapply(outer(z, z)[ut], banda, mean)),
  pares = as.numeric(table(banda))
)
ajuste <- nls(r ~ c0 * exp(-h / a), data = correlograma, weights = pares,
              start = list(c0 = 0.3, a = 50))
rango <- 3 * coef(ajuste)[["a"]] * 1000
write.csv(correlograma, file.path(RESULTADOS, "correlograma_residuos_m2.csv"), row.names = FALSE)
cat("rango blockCV (km):", round(rango_blockcv / 1000, 1),
    "| rango práctico exponencial (km):", round(rango / 1000, 1), "\n")

# --- Bloques ---------------------------------------------------------------
set.seed(SEMILLA)
cv <- cv_spatial(x = puntos, size = rango, k = K_BLOQUES, selection = "random",
                 iteration = 200, seed = SEMILLA, plot = FALSE, progress = FALSE, report = FALSE)
pliegues <- data.frame(ubigeo = d$ubigeo, pliegue = cv$folds_ids)
write.csv(pliegues, file.path(RESULTADOS, "pliegues.csv"), row.names = FALSE)

tab <- table(pliegues$pliegue)
resumen <- c(
  sprintf("fecha: %s", Sys.Date()),
  sprintf("islas unidas: %d; componentes: %d", length(islas), comp$nc),
  sprintf("rango blockCV (km): %.1f; rango práctico exponencial usado (km): %.1f", rango_blockcv / 1000, rango / 1000),
  sprintf("k = %d; distritos por pliegue: mín %d, mediana %d, máx %d", K_BLOQUES, as.integer(min(tab)), as.integer(median(tab)), as.integer(max(tab))),
  sprintf("pliegues con al menos 5 distritos de sierra: %d", sum(tapply(d$sierra, pliegues$pliegue, sum) >= 5)),
  sprintf("negativos (y_tau1 = 0) por pliegue: %s", paste(tapply(1 - d$y, pliegues$pliegue, sum), collapse = " "))
)
writeLines(resumen, file.path(RESULTADOS, "bloques.txt"))
cat(resumen, sep = "\n")
