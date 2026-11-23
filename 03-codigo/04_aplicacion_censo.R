# P1 · Paso 4: aplicación al censo (plan de análisis, secciones 3 a 6; manuscrito, sección 4).
#
# Validación por los 20 bloques de 145 km (04-resultados/pliegues.csv). En cada pliegue, con
# imputación y estandarización ajustadas solo con los distritos de entrenamiento:
#   M1  BYM2 binomial (INLA)                M1F  M1 + fenología NDVI (para H2)
#   M2  M1 sin efectos aleatorios (INLA)    M3   bosque aleatorio (ranger), afinado por bloques
#   M4  MaxEnt (maxnet, rasgos l y q)       M5   multicriterio ECOCROP + TOPSIS (no se ajusta)
# Para M1, M1F y M2 se guarda la predictiva de k_d (mezcla de binomiales sobre 1000 draws de p_d):
# media, cuantiles 5 % y 95 % de p_d, PIT aleatorizado y log score. Para la abstención se guarda
# el índice de disimilitud (DI) de la AOA (CAST), con pesos de la importancia de M3.
#
# Uso: sh 03-codigo/04_aplicacion_lanzar.sh  (los 20 pliegues, cada uno en su propio proceso)
#      o un pliegue:  Rscript 03-codigo/04_aplicacion_censo.R pliegue=7   o la prueba:  ... prueba
suppressPackageStartupMessages({
  library(sf); library(INLA); library(ranger); library(maxnet); library(CAST)
})
inla.setOption(num.threads = "1:1", blas.num.threads = 1)
options(mc.cores = 1)
source("03-codigo/00_comun.R")

args <- commandArgs(TRUE)
PRUEBA <- "prueba" %in% args
UNO <- suppressWarnings(as.integer(sub("^pliegue=", "", args[startsWith(args, "pliegue=")])))
DRAWS <- 1000
Q_REJILLA <- seq(0.02, 0.50, by = 0.01)  # niveles del conjunto {k : |F(k) - 1/2| <= q}
SALIDA <- file.path(RESULTADOS, "aplicacion"); dir.create(file.path(SALIDA, "pliegues"), recursive = TRUE, showWarnings = FALSE)

d <- leer_distritos()
D <- nrow(d)
pl <- read.csv(file.path(RESULTADOS, "pliegues.csv"), colClasses = c(ubigeo = "character"))
stopifnot(identical(pl$ubigeo, d$ubigeo))
pliegue <- pl$pliegue
grafo <- inla.read.graph(file.path(RESULTADOS, "grafo.adj"))
xy <- st_coordinates(st_point_on_surface(st_geometry(st_transform(d, CRS_METRICO)))) / 1000
dist_km <- as.matrix(dist(xy))

# --- M5: multicriterio ECOCROP + TOPSIS (Tabla M5 del plan) ----------------------------
trapecio <- function(x, a, b, c, e) {
  y <- ifelse(x <= a | x >= e, 0, ifelse(x < b, (x - a) / (b - a), ifelse(x <= c, 1, (e - x) / (e - c))))
  y[is.na(y)] <- NA; y
}
topsis_m5 <- function(d) {
  ph <- d$ph; ph[is.na(ph)] <- median(ph, na.rm = TRUE)
  m <- cbind(
    temp = trapecio((d$tmax_c_ajuste + d$tmin_c_ajuste) / 2, 7, 15, 25, 30),
    ppt = trapecio(d$ppt_mm_ajuste, 250, 500, 800, 2000),
    ph = trapecio(ph, 4.2, 5.0, 6.2, 8.5),
    helada = pmax(0, 1 - d$meses_helada_ajuste / 9))
  r <- sweep(m, 2, sqrt(colSums(m^2)), "/")
  v <- r * 0.25
  mejor <- apply(v, 2, max); peor <- apply(v, 2, min)
  s_mas <- sqrt(rowSums(sweep(v, 2, mejor)^2)); s_menos <- sqrt(rowSums(sweep(v, 2, peor)^2))
  ifelse(s_mas + s_menos > 0, s_menos / (s_mas + s_menos), 0)
}
p_m5 <- topsis_m5(d)

# --- Predictiva binomial a partir de los draws de p_d ------------------------------------
# Además de la media y el 5-95 % de p_d, se guardan los límites del conjunto predictivo de k_d / n_d
# para cada q de la rejilla (li_q, ls_q), con los que después se recalibra y se decide la abstención.
predictiva <- function(p, k, n) {
  p <- pmin(pmax(p[is.finite(p)], 1e-12), 1 - 1e-12)  # draws fuera de (0, 1) por interpolación
  kd <- rbinom(length(p), n, p)
  li <- quantile(kd, 0.5 - Q_REJILLA, type = 1, names = FALSE) / n
  ls <- quantile(kd, 0.5 + Q_REJILLA, type = 1, names = FALSE) / n
  c(p_media = mean(p), p_lo = unname(quantile(p, 0.05)), p_hi = unname(quantile(p, 0.95)),
    u = mean(pbinom(k - 1, n, p)) + runif(1) * mean(dbinom(k, n, p)),
    log_score = { l <- dbinom(k, n, p, log = TRUE); max(l) + log(mean(exp(l - max(l)))) },  # log-sum-exp
    setNames(li, sprintf("li%02d", round(100 * Q_REJILLA))),
    setNames(ls, sprintf("ls%02d", round(100 * Q_REJILLA))))
}

# faltan: filas con la respuesta en NA (por defecto, las de prueba; en 07_sensibilidad.R, también
# los vecinos excluidos). Se devuelven las predictivas de las filas de prueba.
ajustar_inla <- function(X, prueba, variables, espacial, faltan = prueba) {
  dd <- cbind(data.frame(k = d$k, n = d$n, id = seq_len(D)), X[, variables])
  dd$k[faltan] <- NA
  rhs <- paste(variables, collapse = " + ")
  if (espacial) rhs <- paste(rhs, "+ f(id, model = 'bym2', graph = grafo, scale.model = TRUE,",
                             "hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01)),",
                             "phi = list(prior = 'pc', param = c(0.5, 2/3))))")
  # Desviaciones 6 y 8: si el optimizador de INLA no converge, se reintenta con otros ajustes numéricos.
  intentos <- list(list(ci = list(), modo = "compact"), list(ci = list(cmin = 0), modo = "compact"),
                   list(ci = list(cmin = 0), modo = "classic"),
                   list(ci = list(cmin = 0, h = 0.001, tolerance = 1e-6), modo = "classic"),
                   list(ci = list(cmin = 0, int.strategy = "eb"), modo = "compact"))
  for (j in seq_along(intentos)) {
    conf <- intentos[[j]]
    fit <- tryCatch(inla(as.formula(paste("k ~", rhs)), family = "binomial", Ntrials = n, data = dd,
                         control.predictor = list(compute = TRUE, link = 1),
                         control.compute = list(return.marginals.predictor = TRUE),
                         control.inla = conf$ci, inla.mode = conf$modo, safe = FALSE, verbose = FALSE,
                         num.threads = "1:1", working.directory = tempfile("inla_"), keep = FALSE),
                    error = function(e) e)
    if (!inherits(fit, "error") && all(vapply(fit$marginals.fitted.values[which(prueba)],
                                              function(m) !is.null(m) && all(is.finite(m)), logical(1)))) break
  }
  if (inherits(fit, "error")) stop("INLA falló con todas las configuraciones")
  cbind(t(sapply(which(prueba), function(i)
    predictiva(inla.rmarginal(DRAWS, fit$marginals.fitted.values[[i]]), d$k[i], d$n[i]))), intento = j)
}

# --- Bosque aleatorio afinado con validación interna por bloques -------------------------
VAR_M3 <- c(setdiff(VARIABLES_M1, c("tmin2", "ppt2")), VARIABLES_M3_EXTRA)
ajustar_rf <- function(X, prueba) {
  tr <- which(!prueba)
  frac <- d$k / d$n
  grupos <- pliegue[tr]
  rejilla <- expand.grid(mtry = c(3, 5, 7), min.node.size = c(5, 20))
  internos <- split(seq_along(tr), (match(grupos, unique(grupos)) %% 5) + 1)
  error <- apply(rejilla, 1, function(h) {
    mean(sapply(internos, function(j) {
      m <- ranger(x = X[tr[-j], VAR_M3], y = frac[tr[-j]], case.weights = d$n[tr[-j]],
                  num.trees = 300, mtry = h[["mtry"]], min.node.size = h[["min.node.size"]],
                  seed = SEMILLA, num.threads = 1)
      p <- predict(m, X[tr[j], VAR_M3], num.threads = 1)$predictions
      weighted.mean((frac[tr[j]] - p)^2, d$n[tr[j]])
    }))
  })
  h <- rejilla[which.min(error), ]
  m <- ranger(x = X[tr, VAR_M3], y = frac[tr], case.weights = d$n[tr], num.trees = 500,
              mtry = h$mtry, min.node.size = h$min.node.size, importance = "permutation",
              seed = SEMILLA, num.threads = 1)
  list(p = predict(m, X[prueba, VAR_M3], num.threads = 1)$predictions,
       importancia = pmax(m$variable.importance, 0))
}

VAR_M4 <- setdiff(VARIABLES_M1, c("tmin2", "ppt2"))
ajustar_maxent <- function(X, prueba) {
  tr <- which(!prueba)
  m <- maxnet(p = d$y[tr], data = X[tr, VAR_M4], f = maxnet.formula(d$y[tr], X[tr, VAR_M4], classes = "lq"))
  as.numeric(predict(m, X[prueba, VAR_M4], type = "cloglog"))
}

# --- Un pliegue -------------------------------------------------------------------------------
correr_pliegue <- function(f) {
  archivo <- file.path(SALIDA, "pliegues", sprintf("pliegue%02d.rds", f))
  if (file.exists(archivo)) return(archivo)
  set.seed(SEMILLA + f)
  prueba <- pliegue == f
  X <- preparar_predictores(d, d[!prueba, ])
  m1 <- ajustar_inla(X, prueba, VARIABLES_M1, TRUE)
  m1f <- ajustar_inla(X, prueba, c(VARIABLES_M1, VARIABLES_FENO), TRUE)
  m2 <- ajustar_inla(X, prueba, VARIABLES_M1, FALSE)
  rf <- ajustar_rf(X, prueba)
  w <- as.data.frame(t(rf$importancia))
  ao <- aoa(newdata = X[prueba, VAR_M3], train = X[!prueba, VAR_M3], variables = VAR_M3,
            weight = w, verbose = FALSE)
  res <- data.frame(
    ubigeo = d$ubigeo[prueba], pliegue = f, sierra = d$sierra[prueba], y = d$y[prueba],
    k = d$k[prueba], n = d$n[prueba],
    delta = apply(dist_km[prueba, !prueba, drop = FALSE], 1, min),
    setNames(as.data.frame(m1), paste0(colnames(m1), "_m1")),
    setNames(as.data.frame(m1f), paste0(colnames(m1f), "_m1f")),
    setNames(as.data.frame(m2), paste0(colnames(m2), "_m2")),
    p_m3 = rf$p, p_m4 = ajustar_maxent(X, prueba), p_m5 = p_m5[prueba],
    di = as.numeric(ao$DI), aoa_umbral = as.numeric(ao$parameters$threshold), dentro_aoa = as.logical(ao$AOA))
  saveRDS(res, archivo)
  cat(format(Sys.time(), "%H:%M"), "pliegue", f, "listo\n", file = file.path(SALIDA, "avance.log"), append = TRUE)
  archivo
}

pliegues_a_correr <- if (PRUEBA) 1 else if (length(UNO) && !is.na(UNO)) UNO else seq_len(K_BLOQUES)
inicio <- Sys.time()
hechos <- lapply(pliegues_a_correr, function(f) tryCatch(correr_pliegue(f), error = function(e) {
  cat(format(Sys.time(), "%H:%M"), "ERROR pliegue", f, conditionMessage(e), "\n",
      file = file.path(SALIDA, "avance.log"), append = TRUE); e }))
errores <- Filter(function(h) inherits(h, "error"), hechos)
cat(sprintf("pliegues: %d | errores: %d | minutos: %.1f\n", length(hechos), length(errores),
            as.numeric(difftime(Sys.time(), inicio, units = "mins"))))
if (length(errores)) print(errores[[1]])
if (PRUEBA) { x <- readRDS(hechos[[1]]); print(dim(x)); print(head(x[, c("ubigeo", "y", "p_media_m1", "u_m1", "p_m3", "p_m4", "p_m5", "di")], 4)) }
