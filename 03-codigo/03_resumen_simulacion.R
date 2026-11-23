# P1 · Paso 3: resumen de la simulación (manuscrito, sección 5.1 y 5.2).
#
# Para cada réplica, diseño y pliegue de prueba:
#   - nominal:            cubierto si |u - 1/2| <= 0,45 (conjunto predictivo posterior al 90 %);
#   - conformal por bloques: q_b = cuantil conformal de los puntajes |u - 1/2| de los OTROS
#                          pliegues del mismo diseño de bloques;
#   - conformal aleatorio: q = cuantil conformal de los puntajes del diseño con pliegues
#                          aleatorios (interpolación), aplicado a los distritos de prueba de bloques.
# El ancho se interpola en la rejilla de q guardada por 02_simulacion.R.
# Salidas: 04-resultados/simulacion/cobertura_*.csv
source("03-codigo/00_comun.R")
ALFA <- 0.10
Q_REJILLA <- seq(0.30, 0.50, by = 0.01)
CORTES_DELTA <- c(0, 10, 25, 50, 100, Inf)

# Uso: Rscript 03-codigo/03_resumen_simulacion.R        (escenario correcto, partes_v2)
#      Rscript 03-codigo/03_resumen_simulacion.R mal    (escenarios mal especificados, partes_mal;
#      la columna "phi" de las salidas lleva entonces el nombre del escenario)
MAL <- "mal" %in% commandArgs(TRUE)
SUFIJO <- if (MAL) "_mal" else ""
partes <- list.files(file.path(RESULTADOS, "simulacion", if (MAL) "partes_mal" else "partes_v2"), full.names = TRUE)
x <- do.call(rbind, lapply(partes, function(a) { r <- readRDS(a); if (!"intento" %in% names(r)) r$intento <- NA; r }))
x$s <- abs(x$u - 0.5)
if (MAL) x$phi <- x$escenario
cat("réplicas:", length(partes), "| filas:", nrow(x), "\n")
if ("intento" %in% names(x)) { cat("intento de INLA usado (filas):\n"); print(table(x$intento, useNA = "ifany")) }

cuantil_conformal <- function(s) {
  m <- length(s); sort(s)[min(m, ceiling((1 - ALFA) * (m + 1)))]
}
# Umbral ponderado por bloque de la Proposición 1: inf{t : (1/B) sum_{c != b} G_c(t) >= 1 - alfa}.
cuantil_bloques <- function(s, grupo, B, alfa = 0.10) {
  t <- sort(unique(s)); g <- split(s, grupo)
  H <- Reduce(`+`, lapply(g, function(sc) findInterval(t, sort(sc)) / length(sc))) / B
  i <- which(H >= 1 - alfa - 1e-12)[1]; if (is.na(i)) Inf else t[i]
}
ancho_en <- function(fila_anchos, q) {
  q <- min(max(q, min(Q_REJILLA)), max(Q_REJILLA))
  approx(Q_REJILLA, fila_anchos, xout = q)$y
}
cols_ancho <- grep("^ancho", names(x))

SIERRA <- leer_distritos()$sierra == 1
# Mitad oriental (la partición de varianza del escenario no_estacionario_este), con los mismos puntos que 02_simulacion.R.
ESTE <- local({ d <- leer_distritos(); x <- sf::st_coordinates(sf::st_point_on_surface(sf::st_geometry(sf::st_transform(d, CRS_METRICO))))[, 1]; x > median(x) })
res <- list()
for (clave in split(seq_len(nrow(x)), list(x$phi, x$rep), drop = TRUE)) {
  r <- x[clave, ]
  aleat <- r[r$diseno == "aleatorio", ]
  q_aleat <- cuantil_conformal(aleat$s)
  for (dis in setdiff(unique(r$diseno), "aleatorio")) {
    b <- r[r$diseno == dis, ]
    for (f in unique(b$pliegue)) {
      prueba <- b$pliegue == f
      q_b <- cuantil_conformal(b$s[!prueba])
      q_w <- cuantil_bloques(b$s[!prueba], b$pliegue[!prueba], length(unique(b$pliegue)), ALFA)
      t <- b[prueba, ]
      # Mondrian por región: cuantil de los otros pliegues dentro de la misma región (exploratorio)
      reg_b <- SIERRA[b$area]
      q_m <- c(sierra = cuantil_conformal(b$s[!prueba & reg_b]), resto = cuantil_conformal(b$s[!prueba & !reg_b]))
      q_t <- q_m[ifelse(SIERRA[t$area], "sierra", "resto")]
      A <- as.matrix(t[, cols_ancho])
      res[[length(res) + 1]] <- data.frame(
        phi = t$phi, rep = t$rep, diseno = dis, area = t$area, delta = t$delta,
        cub_nominal = t$s <= 0.45, cub_bloque = t$s <= q_b, cub_aleatorio = t$s <= q_aleat, cub_mondrian = t$s <= q_t, cub_bloque_pesos = t$s <= q_w,
        ancho_nominal = A[, which.min(abs(Q_REJILLA - 0.45))],
        ancho_bloque = apply(A, 1, ancho_en, q = q_b),
        ancho_aleatorio = apply(A, 1, ancho_en, q = q_aleat))
    }
  }
}
res <- do.call(rbind, res)
res$tramo_delta <- cut(res$delta, CORTES_DELTA, right = FALSE)

media_ic <- function(v) { m <- mean(v); se <- sd(v) / sqrt(length(v)); c(m = m, lo = m - 1.96 * se, hi = m + 1.96 * se) }
# Cobertura por réplica primero (unidad de replicación), luego promedio e IC entre réplicas.
por_rep <- aggregate(cbind(cub_nominal, cub_bloque, cub_aleatorio, ancho_nominal, ancho_bloque, ancho_aleatorio, cub_bloque_pesos) ~
                       phi + diseno + rep, data = res, FUN = mean)
metricas <- c("cub_nominal", "cub_bloque", "cub_aleatorio", "ancho_nominal", "ancho_bloque", "ancho_aleatorio", "cub_bloque_pesos")
global <- do.call(rbind, lapply(split(por_rep, list(por_rep$phi, por_rep$diseno), drop = TRUE), function(g) {
  fila <- data.frame(phi = g$phi[1], diseno = g$diseno[1], replicas = nrow(g))
  for (m in metricas) {
    ic <- media_ic(g[[m]])
    fila[[paste0(m, "_media")]] <- ic[["m"]]; fila[[paste0(m, "_lo")]] <- ic[["lo"]]; fila[[paste0(m, "_hi")]] <- ic[["hi"]]
  }
  fila
}))
write.csv(global, file.path(RESULTADOS, "simulacion", paste0("cobertura_global", SUFIJO, ".csv")), row.names = FALSE)
write.csv(por_rep, file.path(RESULTADOS, "simulacion", paste0("cobertura_por_replica", SUFIJO, ".csv")), row.names = FALSE)

por_delta_rep <- aggregate(cbind(cub_nominal, cub_bloque, cub_aleatorio) ~ phi + diseno + tramo_delta + rep,
                           data = res, FUN = mean)
por_delta <- aggregate(cbind(cub_nominal, cub_bloque, cub_aleatorio) ~ phi + diseno + tramo_delta,
                       data = por_delta_rep, FUN = mean)
# IC 95 % entre réplicas por tramo (para las barras de la figura 2)
ic_tramo <- aggregate(cbind(cub_nominal, cub_bloque) ~ phi + diseno + tramo_delta, data = por_delta_rep,
                      FUN = function(v) 1.96 * sd(v) / sqrt(length(v)))
names(ic_tramo)[4:5] <- c("se_nominal", "se_bloque")
por_delta <- merge(por_delta, ic_tramo)
n_tramo <- aggregate(delta ~ phi + diseno + tramo_delta, data = res, FUN = length)
names(n_tramo)[4] <- "n_distritos_total"
por_delta <- merge(por_delta, n_tramo)
write.csv(por_delta, file.path(RESULTADOS, "simulacion", paste0("cobertura_por_distancia", SUFIJO, ".csv")), row.names = FALSE)

# Cobertura por región (sierra: altitud agrícola media >= 2300 m), media entre réplicas con IC 95 %.
# En el escenario no estacionario la desviación estándar del efecto es 1,5 en la sierra y 0,5 fuera.
res$region <- ifelse(SIERRA[res$area], "sierra", "resto")
por_region_rep <- aggregate(cbind(cub_nominal, cub_bloque, cub_mondrian) ~ phi + diseno + region + rep, data = res, FUN = mean)
por_region <- do.call(rbind, lapply(split(por_region_rep, list(por_region_rep$phi, por_region_rep$diseno, por_region_rep$region), drop = TRUE),
  function(g) { a <- media_ic(g$cub_nominal); b <- media_ic(g$cub_bloque); m <- media_ic(g$cub_mondrian)
    data.frame(phi = g$phi[1], diseno = g$diseno[1], region = g$region[1], replicas = nrow(g),
               cub_nominal = a[["m"]], cub_nominal_lo = a[["lo"]], cub_nominal_hi = a[["hi"]],
               cub_bloque = b[["m"]], cub_bloque_lo = b[["lo"]], cub_bloque_hi = b[["hi"]],
               cub_mondrian = m[["m"]], cub_mondrian_lo = m[["lo"]], cub_mondrian_hi = m[["hi"]]) }))
write.csv(por_region, file.path(RESULTADOS, "simulacion", paste0("cobertura_por_region", SUFIJO, ".csv")), row.names = FALSE)
# Cobertura por mitad este/oeste: muestra qué pasa cuando la partición de la calibración (sierra/resto)
# no coincide con la de la varianza (escenario no_estacionario_este).
res$mitad <- ifelse(ESTE[res$area], "este", "oeste")
por_mitad <- aggregate(cbind(cub_nominal, cub_bloque, cub_mondrian) ~ phi + diseno + mitad + rep, data = res, FUN = mean)
por_mitad <- do.call(rbind, lapply(split(por_mitad, list(por_mitad$phi, por_mitad$diseno, por_mitad$mitad), drop = TRUE), function(g) {
  f <- data.frame(phi = g$phi[1], diseno = g$diseno[1], mitad = g$mitad[1])
  for (m in c("cub_nominal", "cub_bloque", "cub_mondrian")) { ic <- media_ic(g[[m]]); f[[m]] <- ic[["m"]]; f[[paste0(m, "_lo")]] <- ic[["lo"]]; f[[paste0(m, "_hi")]] <- ic[["hi"]] }
  f }))
write.csv(por_mitad, file.path(RESULTADOS, "simulacion", paste0("cobertura_por_mitad", SUFIJO, ".csv")), row.names = FALSE)

options(width = 200)
cat("\nCobertura global (media entre réplicas; objetivo 0,90):\n")
print(global[order(global$diseno, global$phi), c("phi", "diseno", "replicas", "cub_nominal_media", "cub_bloque_media",
                                                  "cub_aleatorio_media", "ancho_nominal_media", "ancho_bloque_media", "cub_bloque_pesos_media")], digits = 3)
cat("\nCobertura por distancia al entrenamiento (bloques de 145 km):\n")
print(por_delta[por_delta$diseno == "bloque145", ], digits = 3)
cat("\nCobertura por mitad este/oeste (bloques de 145 km):\n"); print(por_mitad[por_mitad$diseno == "bloque145", ], digits = 3)
cat("\nCobertura por región (bloques de 145 km):\n")
print(por_region[por_region$diseno == "bloque145", ], digits = 3)
