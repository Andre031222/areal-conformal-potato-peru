# P1 · Paso 6: análisis de la aplicación al censo (plan de análisis, secciones 5, 6 y 10.1-10.2).
#
# Entrada: 04-resultados/aplicacion/pliegues/pliegueNN.rds (predicciones fuera de pliegue, paso 4).
#   1. H1, H1b, H2 y métricas secundarias de M1-M5.
#   2. Calibración: cobertura del conjunto predictivo nominal al 90 % frente al conformal por bloques
#      (puntaje PIT binomial), al conformal sobre la fracción, a GeoConformal y al conformal por
#      tramos del índice de disimilitud (DI), en total y por distancia al entrenamiento.
#   3. Abstención con control de riesgo: umbral lambda del DI elegido con validación cruzada anidada
#      (los pliegues de prueba nunca intervienen en su propio umbral, cuantil ni lambda).
# Salidas: 04-resultados/aplicacion/*.csv
suppressPackageStartupMessages({library(sf); library(spdep); library(pROC)})
source("03-codigo/00_comun.R")
set.seed(SEMILLA)

ALFA <- 0.10            # conjuntos al 90 %
ALFA_R <- 0.10          # riesgo objetivo de la abstención
B <- 2000               # réplicas bootstrap por bloques
Q_REJILLA <- seq(0.02, 0.50, by = 0.01)  # igual que en 04_aplicacion_censo.R
BANDA_KM <- 145         # GeoConformal: ancho de banda = rango práctico de los residuos (desviación 9)
TRAMOS_DI <- 5          # conformal por tramos del DI (desviación 4)
CORTES_DELTA <- c(0, 10, 25, 50, 100, Inf)
SALIDA <- file.path(RESULTADOS, "aplicacion")

archivos <- list.files(file.path(SALIDA, "pliegues"), pattern = "\\.rds$", full.names = TRUE)
x <- do.call(rbind, lapply(archivos, readRDS))
x <- x[order(x$ubigeo), ]
FOLDS <- sort(unique(x$pliegue))
cat("pliegues:", length(FOLDS), "| distritos:", nrow(x), "\n")
x$frac <- x$k / x$n
for (m in c("m1", "m1f", "m2")) x[[paste0("s_", m)]] <- abs(x[[paste0("u_", m)]] - 0.5)

d <- leer_distritos()
d <- d[match(x$ubigeo, d$ubigeo), ]
xy <- st_coordinates(st_point_on_surface(st_geometry(st_transform(d, CRS_METRICO)))) / 1000
dist_km <- as.matrix(dist(xy))

# --- Utilidades ------------------------------------------------------------------------------------
auc <- function(y, p) {
  r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
cuantil_conformal <- function(s, alfa = ALFA) { m <- length(s); sort(s)[min(m, ceiling((1 - alfa) * (m + 1)))] }
cuantil_ponderado <- function(s, w, alfa = ALFA) {   # el punto de prueba aporta peso 1 en +infinito
  o <- order(s); s <- s[o]; w <- w[o]
  acum <- cumsum(w) / (sum(w) + 1)
  j <- which(acum >= 1 - alfa)[1]
  if (is.na(j)) Inf else s[j]
}
youden <- function(y, p) {
  o <- order(p, decreasing = TRUE); y <- y[o]; p <- p[o]
  sens <- cumsum(y == 1) / sum(y == 1); fpr <- cumsum(y == 0) / sum(y == 0)
  p[which.max(sens - fpr)]
}
cp_superior <- function(errores, n) if (n == 0) 0 else qbeta(0.95, errores + 1, n - errores)
# Límites del conjunto predictivo en q (interpolados en la rejilla; q >= 0,5 es [0, 1]).
limites <- function(filas, q, m) {
  li <- as.matrix(x[filas, sprintf("li%02d_%s", round(100 * Q_REJILLA), m)])
  ls <- as.matrix(x[filas, sprintf("ls%02d_%s", round(100 * Q_REJILLA), m)])
  q <- rep_len(q, length(filas))
  inf <- sup <- numeric(length(filas))
  for (i in seq_along(filas)) {
    if (!is.finite(q[i]) || q[i] >= 0.5) { inf[i] <- 0; sup[i] <- 1; next }
    qq <- max(q[i], min(Q_REJILLA))
    inf[i] <- approx(Q_REJILLA, li[i, ], xout = qq)$y; sup[i] <- approx(Q_REJILLA, ls[i, ], xout = qq)$y
  }
  cbind(inf = inf, sup = sup)
}
media_cp <- function(v) { b <- binom.test(sum(v), length(v)); c(est = mean(v), lo = b$conf.int[1], hi = b$conf.int[2]) }

# --- 1. Hipótesis H1, H1b y H2 -------------------------------------------------------------------
puntajes <- c(M1 = "p_media_m1", M1F = "p_media_m1f", M2 = "p_media_m2", M3 = "p_m3", M4 = "p_m4", M5 = "p_m5")
auc_obs <- sapply(puntajes, function(v) auc(x$y, x[[v]]))
por_bloque <- split(seq_len(nrow(x)), x$pliegue)
boot <- t(replicate(B, {
  i <- unlist(por_bloque[sample(length(por_bloque), replace = TRUE)], use.names = FALSE)
  sapply(puntajes, function(v) auc(x$y[i], x[[v]][i]))
}))
delong <- function(a, b) roc.test(roc(x$y, x[[puntajes[a]]], quiet = TRUE), roc(x$y, x[[puntajes[b]]], quiet = TRUE),
                                   method = "delong", paired = TRUE)$p.value
contraste <- function(a, b, umbral) {
  dd <- boot[, a] - boot[, b]
  ic <- quantile(dd, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  data.frame(contraste = paste(a, "-", b), dif_auc = unname(auc_obs[a] - auc_obs[b]),
             ic_lo = ic[1], ic_hi = ic[2], p_delong = delong(a, b), umbral = umbral, se_sostiene = ic[1] > umbral)
}
h1 <- rbind(contraste("M1", "M5", 0.05), contraste("M3", "M5", 0.05))
h1b <- contraste("M1", "M2", 0)
cob_m1 <- media_cp(x$s_m1 <= 0.45); cob_m2 <- media_cp(x$s_m2 <= 0.45)
h1b$cob_m1 <- cob_m1[["est"]]; h1b$cob_m2 <- cob_m2[["est"]]
h1b$se_sostiene <- h1b$se_sostiene && abs(cob_m1[["est"]] - 0.9) < abs(cob_m2[["est"]] - 0.9)

sie <- x[x$sierra == 1, ]
ls_bloque <- aggregate(cbind(log_score_m1, log_score_m1f) ~ pliegue, data = sie, FUN = mean)
ls_bloque <- ls_bloque[ls_bloque$pliegue %in% names(which(table(sie$pliegue) >= 5)), ]
w <- wilcox.test(ls_bloque$log_score_m1f, ls_bloque$log_score_m1, paired = TRUE, alternative = "greater")
h2 <- data.frame(contraste = "log score sierra M1F - M1", bloques = nrow(ls_bloque),
                 dif_mediana = median(ls_bloque$log_score_m1f - ls_bloque$log_score_m1),
                 V = unname(w$statistic), p = w$p.value, se_sostiene = w$p.value < 0.05)
hip <- rbind(
  data.frame(hipotesis = "H1", h1[, c("contraste", "dif_auc", "ic_lo", "ic_hi", "p_delong", "se_sostiene")]),
  data.frame(hipotesis = "H1b", h1b[, c("contraste", "dif_auc", "ic_lo", "ic_hi", "p_delong", "se_sostiene")]))
write.csv(hip, file.path(SALIDA, "hipotesis_auc.csv"), row.names = FALSE)
write.csv(h2, file.path(SALIDA, "hipotesis_h2.csv"), row.names = FALSE)

# Métricas secundarias. Brier ponderado y Moran sobre la fracción k/n (modelos que la predicen);
# calibración logística de y sobre logit(puntaje) para todos.
# El grafo sigue el orden por ubigeo de leer_distritos(), el mismo que x con los 20 pliegues.
lineas <- readLines(file.path(RESULTADOS, "grafo.adj"))[-1]
vec <- lapply(strsplit(trimws(lineas), "\\s+"), function(v) as.integer(v[-(1:2)]))
nb <- structure(lapply(vec, function(v) if (length(v)) v else 0L), class = "nb")
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
metricas <- do.call(rbind, lapply(names(puntajes), function(mod) {
  p <- x[[puntajes[mod]]]
  pc <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  cal <- coef(glm(x$y ~ qlogis(pc), family = binomial))
  de_fraccion <- mod %in% c("M1", "M1F", "M2", "M3")
  m <- tolower(mod)
  data.frame(modelo = mod, auc = auc_obs[[mod]],
             auc_ic_lo = quantile(boot[, mod], 0.025, names = FALSE), auc_ic_hi = quantile(boot[, mod], 0.975, names = FALSE),
             brier_fraccion = if (de_fraccion) weighted.mean((x$frac - p)^2, x$n) else NA,
             log_score = if (mod %in% c("M1", "M1F", "M2")) mean(x[[paste0("log_score_", m)]]) else NA,
             cobertura_90 = if (mod %in% c("M1", "M1F", "M2")) mean(x[[paste0("s_", m)]] <= 0.45) else NA,
             cal_intercepto = cal[[1]], cal_pendiente = cal[[2]],
             moran_residuos = if (de_fraccion && nrow(x) == length(nb)) unname(moran.test(x$frac - p, lw, zero.policy = TRUE)$estimate[1]) else NA)
}))
write.csv(metricas, file.path(SALIDA, "metricas_modelos.csv"), row.names = FALSE)

# --- 2. Calibración: cinco reglas para el conjunto al 90 % -----------------------------------------
# Para cada pliegue de prueba f, los cuantiles salen solo de los distritos de los otros pliegues.
tramo_di <- function(di, cortes) findInterval(di, cortes, all.inside = TRUE)
# Umbral ponderado por bloque de la Proposición 1: inf{t : (1/B) sum_{c != b} G_c(t) >= 1 - alfa}.
cuantil_bloques <- function(s, grupo, B, alfa = 0.10) {
  t <- sort(unique(s)); g <- split(s, grupo)
  H <- Reduce(`+`, lapply(g, function(sc) findInterval(t, sort(sc)) / length(sc))) / B
  i <- which(H >= 1 - alfa - 1e-12)[1]; if (is.na(i)) Inf else t[i]
}
calibrar <- function(m, excluir = integer(0)) {
  s <- x[[paste0("s_", m)]]
  q <- rep(NA_real_, nrow(x)); q_geo <- q; q_di <- q; q_frac <- q; q_mon <- q; q_pes <- q
  for (f in setdiff(FOLDS, excluir)) {
    te <- which(x$pliegue == f); ca <- which(!x$pliegue %in% c(f, excluir))
    q[te] <- cuantil_conformal(s[ca])
    q_pes[te] <- cuantil_bloques(s[ca], x$pliegue[ca], length(setdiff(FOLDS, excluir)))
    q_frac[te] <- cuantil_conformal(abs(x$frac[ca] - x[[paste0("p_media_", m)]][ca]))
    W <- exp(-0.5 * (dist_km[te, ca, drop = FALSE] / BANDA_KM)^2)
    q_geo[te] <- apply(W, 1, function(wi) cuantil_ponderado(s[ca], wi))
    cortes <- quantile(x$di[ca], seq(0, 1, length.out = TRAMOS_DI + 1), names = FALSE)
    tc <- tramo_di(x$di[ca], cortes); tt <- tramo_di(x$di[te], cortes)
    qt <- sapply(seq_len(TRAMOS_DI), function(j) cuantil_conformal(s[ca][tc == j]))
    q_di[te] <- qt[tt]
    # Mondrian por región (exploratorio, desviación 12): cuantil de los otros pliegues en la misma región
    for (r in 0:1) { tr <- te[x$sierra[te] == r]; if (length(tr)) q_mon[tr] <- cuantil_conformal(s[ca][x$sierra[ca] == r]) }
  }
  list(s = s, q = q, q_geo = q_geo, q_di = q_di, q_frac = q_frac, q_mon = q_mon, q_pes = q_pes)
}
cal <- list(m1 = calibrar("m1"), m2 = calibrar("m2"))
filas_cal <- do.call(rbind, lapply(names(cal), function(m) {
  k <- cal[[m]]; todas <- seq_len(nrow(x))
  reglas <- list(nominal = rep(0.45, nrow(x)), bloque = k$q, geoconformal = k$q_geo, tramos_di = k$q_di, mondrian = k$q_mon, bloque_pesos = k$q_pes)
  out <- do.call(rbind, lapply(names(reglas), function(r) {
    lim <- limites(todas, reglas[[r]], m)
    data.frame(modelo = toupper(m), regla = r, ubigeo = x$ubigeo, pliegue = x$pliegue, delta = x$delta, sierra = x$sierra,
               cubre = k$s <= reglas[[r]], ancho = lim[, "sup"] - lim[, "inf"])
  }))
  pm <- x[[paste0("p_media_", m)]]
  rbind(out, data.frame(modelo = toupper(m), regla = "fraccion", ubigeo = x$ubigeo, pliegue = x$pliegue,
                        delta = x$delta, sierra = x$sierra, cubre = abs(x$frac - pm) <= k$q_frac,
                        ancho = pmin(1, pm + k$q_frac) - pmax(0, pm - k$q_frac)))
}))
filas_cal$tramo_delta <- cut(filas_cal$delta, CORTES_DELTA, right = FALSE)
cob <- do.call(rbind, lapply(split(filas_cal, list(filas_cal$modelo, filas_cal$regla), drop = TRUE), function(g) {
  ic <- media_cp(g$cubre)
  data.frame(modelo = g$modelo[1], regla = g$regla[1], n = nrow(g), cobertura = ic[["est"]],
             ic_lo = ic[["lo"]], ic_hi = ic[["hi"]], ancho_medio = mean(g$ancho),
             cob_min_bloque = min(tapply(g$cubre, g$pliegue, mean)), cob_max_bloque = max(tapply(g$cubre, g$pliegue, mean)))
}))
cob_delta <- do.call(rbind, lapply(split(filas_cal, list(filas_cal$modelo, filas_cal$regla, filas_cal$tramo_delta), drop = TRUE),
  function(g) { ic <- media_cp(g$cubre)
    data.frame(modelo = g$modelo[1], regla = g$regla[1], tramo_delta = g$tramo_delta[1], n = nrow(g),
               cobertura = ic[["est"]], ic_lo = ic[["lo"]], ic_hi = ic[["hi"]], ancho_medio = mean(g$ancho)) }))
write.csv(cob, file.path(SALIDA, "cobertura_censo.csv"), row.names = FALSE)
write.csv(cob_delta, file.path(SALIDA, "cobertura_censo_distancia.csv"), row.names = FALSE)
# Ancho del conjunto al 90 % de M1 por distrito (figura 6b)
an <- filas_cal[filas_cal$modelo == "M1" & filas_cal$regla %in% c("nominal", "bloque", "mondrian"), c("ubigeo", "regla", "ancho", "cubre")]
write.csv(reshape(an, idvar = "ubigeo", timevar = "regla", direction = "wide"), file.path(SALIDA, "ancho_por_distrito.csv"), row.names = FALSE)
cob_region <- do.call(rbind, lapply(split(filas_cal, list(filas_cal$modelo, filas_cal$regla, filas_cal$sierra), drop = TRUE),
  function(g) { ic <- media_cp(g$cubre)
    data.frame(modelo = g$modelo[1], regla = g$regla[1], region = if (g$sierra[1] == 1) "sierra" else "resto", n = nrow(g),
               cobertura = ic[["est"]], ic_lo = ic[["lo"]], ic_hi = ic[["hi"]], ancho_medio = mean(g$ancho)) }))
write.csv(cob_region, file.path(SALIDA, "cobertura_censo_region.csv"), row.names = FALSE)

# --- 3. Abstención con control de riesgo -------------------------------------------------------
# Umbral de clasificación (Youden) y cuantil conformal de cada pliegue, sin usar ese pliegue.
# Abstención propuesta: el conjunto conformal por bloques contiene el umbral, o DI > lambda;
# lambda es el más permisivo cuya cota superior de Clopper-Pearson (95 %) del error en los retenidos
# no supera ALFA_R, calculada en los otros pliegues con sus propios umbrales y cuantiles (anidado).
REJILLA_LAMBDA <- c(quantile(x$di, seq(0.05, 1, by = 0.01), names = FALSE), Inf)
pm <- x$p_media_m1
umbral_sin <- function(excluir) youden(x$y[!x$pliegue %in% excluir], pm[!x$pliegue %in% excluir])
decision <- function(filas, umbral, q) {   # TRUE = el conjunto contiene el umbral (no concluyente)
  lim <- limites(filas, q, "m1"); lim[, "inf"] <= umbral & umbral <= lim[, "sup"]
}
t_f <- setNames(sapply(FOLDS, umbral_sin), FOLDS)
x$umbral <- t_f[as.character(x$pliegue)]
x$yhat <- as.integer(pm >= x$umbral)
x$error <- as.integer(x$yhat != x$y)
x$nc_bloque <- decision(seq_len(nrow(x)), x$umbral, cal$m1$q)
x$nc_geo <- decision(seq_len(nrow(x)), x$umbral, cal$m1$q_geo)
x$nc_tramos <- decision(seq_len(nrow(x)), x$umbral, cal$m1$q_di)
x$nc_bayes <- x$p_lo_m1 <= x$umbral & x$umbral <= x$p_hi_m1
x$fuera_aoa <- x$di > x$aoa_umbral

elegir_lambda <- function(nc, error, di) {
  ok <- sapply(REJILLA_LAMBDA, function(l) { r <- !nc & di <= l; cp_superior(sum(error[r]), sum(r)) <= ALFA_R })
  if (any(ok)) max(REJILLA_LAMBDA[ok]) else -Inf
}
x$lambda <- NA_real_
for (f in FOLDS) {
  otras <- which(x$pliegue != f)
  q_in <- calibrar("m1", excluir = f)$q                      # cuantiles sin f ni el pliegue propio
  t_in <- sapply(FOLDS, function(g) if (g == f) NA else umbral_sin(c(f, g)))
  u_in <- t_in[match(x$pliegue[otras], FOLDS)]
  nc_in <- decision(otras, u_in, q_in[otras])
  err_in <- as.integer(as.integer(pm[otras] >= u_in) != x$y[otras])
  x$lambda[x$pliegue == f] <- elegir_lambda(nc_in, err_in, x$di[otras])
}
x$abst_propuesta <- x$nc_bloque | x$di > x$lambda

# --- 3b. Abstención calibrada al riesgo (exploratoria; desviación 10, decidida tras ver la 3) -------
# En vez de fijar el conjunto conformal al 90 %, se elige su nivel 1 - a (a en ALFAS_CONJ) con la misma
# validación anidada: el a más alto (conjunto más estrecho, menos abstención) cuya cota de
# Clopper-Pearson del error en los retenidos de los otros pliegues no supera ALFA_R. El nivel mínimo
# es 10 % (la rejilla de límites guardada empieza en q = 0,02).
ALFAS_CONJ <- seq(0.10, 0.90, by = 0.05)
q_bloque <- function(alfa, excluir = integer(0)) {
  q <- rep(NA_real_, nrow(x))
  for (f in setdiff(FOLDS, excluir)) {
    te <- x$pliegue == f; ca <- !x$pliegue %in% c(f, excluir)
    q[te] <- cuantil_conformal(x$s_m1[ca], alfa)
  }
  q
}
x$alfa_conj <- NA_real_; x$nc_riesgo <- NA
for (f in FOLDS) {
  otras <- which(x$pliegue != f); te <- which(x$pliegue == f)
  t_in <- sapply(FOLDS, function(g) if (g == f) NA else umbral_sin(c(f, g)))
  u_in <- t_in[match(x$pliegue[otras], FOLDS)]
  err_in <- as.integer(as.integer(pm[otras] >= u_in) != x$y[otras])
  ok <- sapply(ALFAS_CONJ, function(a) {
    r <- !decision(otras, u_in, q_bloque(a, excluir = f)[otras])
    cp_superior(sum(err_in[r]), sum(r)) <= ALFA_R
  })
  a <- if (any(ok)) max(ALFAS_CONJ[ok]) else NA
  x$alfa_conj[te] <- a
  x$nc_riesgo[te] <- if (is.na(a)) TRUE else decision(te, x$umbral[te], q_bloque(a)[te])
}

reglas_abst <- list(
  "a. intervalo bayesiano 90 %" = x$nc_bayes,
  "b. solo AOA" = x$fuera_aoa,
  "c. GeoConformal" = x$nc_geo,
  "d. conformal por tramos del DI" = x$nc_tramos,
  "e. conformal por bloques" = x$nc_bloque,
  "f. propuesta: conformal por bloques + DI con control de riesgo" = x$abst_propuesta,
  "g. exploratoria: nivel del conjunto conformal calibrado al riesgo" = x$nc_riesgo,
  "sin abstención" = rep(FALSE, nrow(x)))
reglas_abst <- lapply(reglas_abst, function(a) a | is.na(a))   # sin información = abstenerse
abst <- do.call(rbind, lapply(names(reglas_abst), function(r) {
  a <- reglas_abst[[r]]; ret <- !a
  e <- sum(x$error[ret]); n <- sum(ret)
  ci <- if (n) binom.test(e, n)$conf.int else c(NA, NA)
  data.frame(regla = r, retenidos = n, fraccion_retenida = n / nrow(x), riesgo = if (n) e / n else NA,
             riesgo_ic_lo = ci[1], riesgo_ic_hi = ci[2], cp_superior = cp_superior(e, n),
             auc_retenidos = if (n && length(unique(x$y[ret])) == 2) auc(x$y[ret], pm[ret]) else NA,
             riesgo_max_bloque = max(tapply(x$error[ret], x$pliegue[ret], mean)),
             delta_mediana_abstenidos = if (any(a)) median(x$delta[a]) else NA)
}))
write.csv(abst, file.path(SALIDA, "abstencion.csv"), row.names = FALSE)

# Curva riesgo-cobertura: se barre lambda (con los conjuntos conformales ya cruzados) y, para
# comparar, se retiene por DI creciente sin la regla del intervalo.
curva <- do.call(rbind, lapply(REJILLA_LAMBDA, function(l) {
  r1 <- !x$nc_bloque & x$di <= l; r2 <- x$di <= l
  data.frame(lambda = l, metodo = c("conformal por bloques + DI", "solo DI"),
             fraccion_retenida = c(mean(r1), mean(r2)),
             riesgo = c(if (any(r1)) mean(x$error[r1]) else NA, if (any(r2)) mean(x$error[r2]) else NA))
}))
write.csv(curva, file.path(SALIDA, "riesgo_cobertura.csv"), row.names = FALSE)
write.csv(x[, c("ubigeo", "pliegue", "sierra", "y", "k", "n", "delta", "di", "p_media_m1", "umbral", "yhat",
                "error", "nc_bloque", "nc_bayes", "fuera_aoa", "lambda", "abst_propuesta", "alfa_conj", "nc_riesgo")],
          file.path(SALIDA, "decision_por_distrito.csv"), row.names = FALSE)

# Salidas para las figuras: cobertura por bloque, curva riesgo-cobertura barriendo el nivel del
# conjunto conformal (cruzado, sin anidar) y riesgo de los retenidos por bloque.
write.csv(aggregate(cbind(cobertura = cubre, ancho) ~ modelo + regla + pliegue, data = filas_cal, FUN = mean),
          file.path(SALIDA, "cobertura_por_bloque.csv"), row.names = FALSE)
curva_alfa <- do.call(rbind, lapply(ALFAS_CONJ, function(a) {
  r <- !decision(seq_len(nrow(x)), x$umbral, q_bloque(a))
  data.frame(alfa = a, metodo = "nivel del conjunto conformal", fraccion_retenida = mean(r), riesgo = mean(x$error[r]))
}))
write.csv(curva_alfa, file.path(SALIDA, "riesgo_cobertura_alfa.csv"), row.names = FALSE)
riesgo_bloque <- do.call(rbind, lapply(c("sin abstención", names(reglas_abst)[grepl("^[efg]\\.", names(reglas_abst))]), function(r) {
  ret <- !reglas_abst[[r]]
  g <- aggregate(cbind(error = x$error[ret]), by = list(pliegue = x$pliegue[ret]), FUN = mean)
  g$retenidos <- as.vector(table(factor(x$pliegue[ret], levels = g$pliegue)))
  g$distritos <- as.vector(table(x$pliegue)[as.character(g$pliegue)]); g$regla <- r; g
}))
write.csv(riesgo_bloque, file.path(SALIDA, "riesgo_por_bloque.csv"), row.names = FALSE)

# --- Resumen en pantalla -----------------------------------------------------------------------
options(width = 200)
cat("\nHipótesis (AUC; IC bootstrap por bloques, B =", B, "):\n"); print(hip, digits = 3)
cat("\nH2:\n"); print(h2, digits = 3)
cat("\nMétricas de los modelos:\n"); print(metricas, digits = 3)
cat("\nCobertura del conjunto al 90 %:\n"); print(cob, digits = 3)
cat("\nCobertura por región (M1):\n"); print(cob_region[cob_region$modelo == "M1", ], digits = 3)
cat("\nCobertura por distancia (M1):\n"); print(cob_delta[cob_delta$modelo == "M1", ], digits = 3)
cat("\nAbstención (riesgo objetivo", ALFA_R, "):\n"); print(abst, digits = 3)
cat("\nlambda por pliegue:\n"); print(tapply(x$lambda, x$pliegue, unique))
cat("\nnivel alfa del conjunto elegido por pliegue (exploratoria):\n"); print(tapply(x$alfa_conj, x$pliegue, unique))
if ("intento_m1" %in% names(x)) { cat("\nintento de INLA (M1, M1F, M2):\n"); print(table(x$intento_m1)); print(table(x$intento_m1f)); print(table(x$intento_m2)) }
