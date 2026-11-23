# P1 · Paso 2: estudio de simulación (plan de análisis, sección 10.3; manuscrito, sección 3).
#
# Grafo real de los 1 783 distritos y sus n_d del censo. Para cada phi y réplica:
#   - dos covariables suaves (campo gaussiano exponencial, rango 200 km), estandarizadas;
#   - efecto BYM2 con sigma_b = 1 y el ICAR escalado (varianza generalizada 1);
#   - logit p = -1 + 0.8 x1 - 0.5 x2 + b;  k ~ Binomial(n_d, p).
# Diseños de validación (20 pliegues cada uno): bloques de 25, 75 y 145 km, y pliegues
# aleatorios (para el conformal "aleatorio"). En cada pliegue se ajusta BYM2 con INLA con la
# respuesta del pliegue en NA y se guarda, por distrito de prueba:
#   - u: PIT aleatorizado bajo la predictiva (mezcla de binomiales sobre 1000 draws de p_d);
#   - delta: distancia (km) al distrito de entrenamiento más cercano;
#   - ancho del conjunto {k : |F(k) - 1/2| <= q} / n_d para una rejilla de q.
# La cobertura nominal y la conformal se calculan después, en 03_resumen_simulacion.R.
#
# Uso: sh 03-codigo/02_simulacion_lanzar.sh  (corre las 60 tareas, 8 a la vez, cada una en su
#      propio proceso de R)  o una sola tarea:  Rscript 03-codigo/02_simulacion.R phi=0.9 rep=3
suppressPackageStartupMessages({library(sf); library(INLA); library(blockCV)})
inla.setOption(num.threads = "1:1", blas.num.threads = 1)
options(mc.cores = 1)
source("03-codigo/00_comun.R")

ARGS <- commandArgs(TRUE)
argumento <- function(nombre, defecto) {
  a <- sub(paste0("^", nombre, "="), "", ARGS[startsWith(ARGS, paste0(nombre, "="))])
  if (length(a)) as.numeric(a) else defecto
}
PILOTO <- "piloto" %in% ARGS
# Desviación 11: escenarios mal especificados (el modelo ajustado sigue siendo el BYM2 gaussiano).
#   correcto          efecto BYM2 gaussiano, como en el modelo ajustado
#   colas             parte no estructurada t de Student con 3 g. l. (varianza 1): distritos atípicos
#   no_estacionario   desviación estándar del efecto 1,5 en la sierra y 0,5 fuera de ella
#   no_estacionario_este  igual, pero 1,5 en la mitad oriental de los distritos (coordenada x mayor
#                     que la mediana): regiones de varianza que no coinciden con la sierra (desviación 13)
ESCENARIO <- { a <- sub("^escenario=", "", ARGS[startsWith(ARGS, "escenario=")]); if (length(a)) a else "correcto" }
stopifnot(ESCENARIO %in% c("correcto", "colas", "no_estacionario", "no_estacionario_este"))
PHI <- argumento("phi", c(0.2, 0.5, 0.9))          # p. ej. phi=0.9 corre solo ese escenario
REPLICAS <- argumento("rep", if (PILOTO) 1 else 1:20)   # desviación 5: 30 -> 20 (16-nov-2026)
LADOS_KM <- c(25, 75, 145)
DRAWS <- 1000
Q_REJILLA <- seq(0.30, 0.50, by = 0.01)
BETA <- c(-1, 0.8, -0.5)
# Desviación 7: cada tarea corre en su propio proceso de R (sin mclapply) y cada ajuste de INLA
# en su propia carpeta de trabajo; los resultados van a partes_v2 (corrida completa nueva).
SALIDA <- file.path(RESULTADOS, "simulacion", if (ESCENARIO == "correcto") "partes_v2" else "partes_mal")
dir.create(SALIDA, recursive = TRUE, showWarnings = FALSE)

d <- leer_distritos()
n <- d$n
D <- nrow(d)
pts <- st_point_on_surface(st_geometry(st_transform(d, CRS_METRICO)))
xy <- st_coordinates(pts) / 1000
dist_km <- as.matrix(dist(xy))
grafo <- inla.read.graph(file.path(RESULTADOS, "grafo.adj"))

# Precisión ICAR escalada y su descomposición para simular.
Q <- -inla.graph2matrix(grafo); diag(Q) <- 0; diag(Q) <- -rowSums(Q)
Q <- inla.scale.model(Q, constr = list(A = matrix(1, 1, D), e = 0))
eig <- eigen(as.matrix(Q), symmetric = TRUE)
pos <- eig$values > 1e-8
simular_icar <- function() {
  z <- rnorm(sum(pos))
  u <- eig$vectors[, pos] %*% (z / sqrt(eig$values[pos]))
  as.numeric(u - mean(u))
}
L_cov <- chol(exp(-dist_km / (200 / 3)) + diag(1e-6, D))
campo <- function() { x <- as.numeric(crossprod(L_cov, rnorm(D))); (x - mean(x)) / sd(x) }

# Pliegues: bloques por lado y aleatorios. Fijos para todas las réplicas. Si ya existen, se
# reutilizan tal cual (así otra computadora usa exactamente los mismos pliegues).
archivo_pliegues <- file.path(RESULTADOS, "simulacion", "pliegues_simulacion.rds")
if (file.exists(archivo_pliegues)) {
  pliegues <- readRDS(archivo_pliegues)
} else {
  pliegues <- list()
  for (lado in LADOS_KM) {
    cv <- cv_spatial(x = st_sf(id = seq_len(D), geometry = pts), size = lado * 1000, k = K_BLOQUES,
                     selection = "random", iteration = 100, seed = SEMILLA,
                     plot = FALSE, progress = FALSE, report = FALSE)
    pliegues[[paste0("bloque", lado)]] <- cv$folds_ids
  }
  set.seed(SEMILLA); pliegues[["aleatorio"]] <- sample(rep_len(seq_len(K_BLOQUES), D))
  saveRDS(pliegues, archivo_pliegues)
}

formula <- k ~ x1 + x2 + f(id, model = "bym2", graph = grafo, scale.model = TRUE,
                           hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                                        phi = list(prior = "pc", param = c(0.5, 2 / 3))))

registrar <- function(texto) cat(format(Sys.time(), "%H:%M"), texto, "\n",
                                  file = file.path(RESULTADOS, "simulacion", "avance.log"), append = TRUE)

# En algunos pliegues el optimizador de Newton-Raphson de INLA no converge (error 12).
# Desviación 6: si falla, se reintenta con cmin = 0 y, si vuelve a fallar, además en modo
# clásico. El modelo y los priors no cambian; solo el método numérico de búsqueda de la moda.
# Desviación 8 (17-nov-2026): en unos pocos pliegues de 25 km con phi alto fallan los tres; se
# añaden un paso de diferenciación menor (h = 0.001) y, como último recurso, la integración
# empírica de Bayes sobre los hiperparámetros (int.strategy = "eb"). Se guarda qué intento sirvió.
CONFIG_INTENTOS <- list(
  list(control.inla = list(), inla.mode = "compact"),
  list(control.inla = list(cmin = 0), inla.mode = "compact"),
  list(control.inla = list(cmin = 0), inla.mode = "classic"),
  list(control.inla = list(cmin = 0, h = 0.001, tolerance = 1e-6), inla.mode = "classic"),
  list(control.inla = list(cmin = 0, int.strategy = "eb"), inla.mode = "compact"))
contador_inla <- 0
ajustar_inla <- function(dd, filas) {
  for (intento in seq_along(CONFIG_INTENTOS)) {
    contador_inla <<- contador_inla + 1
    carpeta <- file.path(tempdir(), sprintf("inla_%d_%04d", Sys.getpid(), contador_inla))
    conf <- CONFIG_INTENTOS[[intento]]
    fit <- tryCatch(inla(formula, family = "binomial", Ntrials = n, data = dd,
                         control.predictor = list(compute = TRUE, link = 1),
                         control.compute = list(config = FALSE, return.marginals.predictor = TRUE),
                         control.inla = conf$control.inla, inla.mode = conf$inla.mode,
                         num.threads = "1:1", verbose = FALSE, safe = FALSE,
                         working.directory = carpeta, keep = FALSE),
                    error = function(e) e)
    unlink(carpeta, recursive = TRUE)
    ok <- !inherits(fit, "error") && length(fit$marginals.fitted.values) >= max(filas) &&
      all(vapply(fit$marginals.fitted.values[filas], function(m) !is.null(m) && all(is.finite(m)), logical(1)))
    if (ok) { fit$intento <- intento; return(fit) }
    registrar(sprintf("aviso: intento %d de INLA fallido (%s)", intento,
                      if (inherits(fit, "error")) conditionMessage(fit) else "marginales no finitas"))
  }
  stop("INLA falló con las ", length(CONFIG_INTENTOS), " configuraciones")
}

ajustar_pliegue <- function(datos, prueba) {
  dd <- datos; dd$k[prueba] <- NA
  filas <- which(prueba)
  fit <- ajustar_inla(dd, filas)
  out <- lapply(filas, function(i) {
    # La interpolación de la marginal puede dar unos pocos draws fuera de (0, 1); se acotan.
    p <- inla.rmarginal(DRAWS, fit$marginals.fitted.values[[i]])
    p <- pmin(pmax(p[is.finite(p)], 1e-12), 1 - 1e-12)
    k <- datos$k[i]
    Fm1 <- mean(pbinom(k - 1, n[i], p)); f <- mean(dbinom(k, n[i], p))
    u <- Fm1 + runif(1) * f
    kd <- rbinom(DRAWS, n[i], p)
    ancho <- (quantile(kd, 0.5 + Q_REJILLA, type = 1, names = FALSE) -
              quantile(kd, 0.5 - Q_REJILLA, type = 1, names = FALSE)) / n[i]
    c(u = u, ancho = ancho)
  })
  m <- do.call(rbind, out)
  delta <- apply(dist_km[filas, !prueba, drop = FALSE], 1, min)
  data.frame(area = filas, delta = delta, m, intento = fit$intento, check.names = FALSE)
}

correr <- function(tarea) {
  phi <- tarea$phi; r <- tarea$rep
  archivo <- file.path(SALIDA, sprintf("%sphi%.1f_rep%02d.rds", if (ESCENARIO == "correcto") "" else paste0(ESCENARIO, "_"), phi, r))
  if (file.exists(archivo)) return(archivo)
  set.seed(SEMILLA + round(phi * 1000) + r)
  x1 <- campo(); x2 <- campo()
  v <- rnorm(D); u_icar <- simular_icar()
  if (ESCENARIO == "colas") v <- rt(D, df = 3) / sqrt(3)
  b <- sqrt(1 - phi) * v + sqrt(phi) * u_icar
  if (ESCENARIO == "no_estacionario") b <- b * ifelse(d$sierra == 1, 1.5, 0.5)
  if (ESCENARIO == "no_estacionario_este") b <- b * ifelse(xy[, 1] > median(xy[, 1]), 1.5, 0.5)
  eta <- BETA[1] + BETA[2] * x1 + BETA[3] * x2 + b
  datos <- data.frame(k = rbinom(D, n, plogis(eta)), n = n, x1 = x1, x2 = x2, id = seq_len(D))
  res <- list()
  disenos <- if (ESCENARIO == "correcto") names(pliegues) else c("bloque145", "aleatorio")
  for (dis in disenos) {
    for (f in seq_len(K_BLOQUES)) {
      prueba <- pliegues[[dis]] == f
      if (!any(prueba)) next
      out <- tryCatch(ajustar_pliegue(datos, prueba), error = function(e) {
        registrar(sprintf("ERROR phi=%.1f rep=%02d %s pliegue %d: %s", phi, r, dis, f, conditionMessage(e)))
        stop(e)
      })
      out$diseno <- dis; out$pliegue <- f
      res[[length(res) + 1]] <- out
      if (f == K_BLOQUES) registrar(sprintf("phi=%.1f rep=%02d %s listo", phi, r, dis))
      if (PILOTO && length(res) >= 3) break
    }
    if (PILOTO) break
  }
  res <- do.call(rbind, res)
  res$phi <- phi; res$rep <- r; res$escenario <- ESCENARIO
  saveRDS(res, archivo)
  archivo
}

tareas <- unlist(lapply(PHI, function(phi) lapply(REPLICAS, function(r) list(phi = phi, rep = r))),
                 recursive = FALSE)
if (PILOTO) tareas <- tareas[1]
for (t in tareas) {
  inicio <- Sys.time()
  correr(t)
  registrar(sprintf("TAREA phi=%.1f rep=%02d terminada en %.1f min", t$phi, t$rep,
                    as.numeric(difftime(Sys.time(), inicio, units = "mins"))))
}
