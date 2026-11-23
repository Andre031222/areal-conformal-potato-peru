# P1 · Paso 7: análisis de sensibilidad (plan de análisis, sección 7) y validación secundaria
# (sección 4: dejar fuera un distrito y sus vecinos). Reutiliza las funciones de 04_aplicacion_censo.R.
#
# Casos (mismos 20 pliegues de 145 km, salvo "vecinos"):
#   superficie   M1 con k = hectáreas de papa y n = hectáreas sembradas (redondeadas)       (7.2)
#   altitud_ph   M1 con altitud media y pH además de sus variables                           (7.3)
#   iid          M2 con un efecto aleatorio iid (sobredispersión sin estructura espacial)    (desv. 13)
#   pisco        M1 con el clima de PISCO en vez de TerraClimate                               (7.4)
#   tau0.5, tau5 M4 (MaxEnt) reajustado con la etiqueta y_tau0.5 o y_tau5; M1, M3 y M5 no
#                cambian porque no usan la etiqueta                                          (7.1)
#   vecinos      M1 dejando fuera cada distrito y sus vecinos de contigüidad; desde= y hasta=
#                eligen el tramo de distritos                                                 (4)
# La sensibilidad 7.6 (sin los 17 distritos urbanos) solo excluye distritos al evaluar y se hace en
# 08_analisis_sensibilidad.R.
#
# Uso: Rscript 03-codigo/07_sensibilidad.R caso=superficie pliegue=3
#      Rscript 03-codigo/07_sensibilidad.R caso=vecinos desde=1 hasta=50
lineas <- readLines("03-codigo/04_aplicacion_censo.R")
if (any(commandArgs(TRUE) == "caso=iid")) {   # M2 con efecto aleatorio no espacial (iid) en vez de BYM2
  lineas <- sub("f(id, model = 'bym2', graph = grafo, scale.model = TRUE,", "f(id, model = 'iid',", lineas, fixed = TRUE)
  lineas <- sub("\"hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01)),\",", "\"hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))\")", lineas, fixed = TRUE)
  lineas <- lineas[!grepl("phi = list(prior = 'pc', param = c(0.5, 2/3))))", lineas, fixed = TRUE)]
}
eval(parse(text = lineas[1:(grep("^# --- Un pliegue", lineas) - 1)]))

argumento <- function(nombre) {
  a <- sub(paste0("^", nombre, "="), "", args[startsWith(args, paste0(nombre, "="))])
  if (length(a)) a else NA
}
CASO <- argumento("caso")
stopifnot(CASO %in% c("superficie", "altitud_ph", "tau0.5", "tau5", "pisco", "iid", "vecinos"))
SAL <- file.path(RESULTADOS, "sensibilidad", CASO)
dir.create(SAL, recursive = TRUE, showWarnings = FALSE)
registrar <- function(texto) cat(format(Sys.time(), "%H:%M"), CASO, texto, "\n",
                                  file = file.path(RESULTADOS, "sensibilidad", "avance.log"), append = TRUE)

if (CASO == "superficie") {
  d$k <- round(d$ha_papa); d$n <- pmax(round(d$ha_sembrada), d$k, 1)
}
VARS <- if (CASO == "altitud_ph") c(VARIABLES_M1, "altitud", "ph") else VARIABLES_M1
if (CASO %in% c("tau0.5", "tau5")) d$y <- d[[paste0("y_", CASO)]]
if (CASO == "pisco") {   # clima de M1 con PISCO (data-raw/clima_pisco.py), campañas 1991-2012
  p <- read.csv(file.path(DATOS, "clima_pisco_distrito_2012.csv"), colClasses = c(ubigeo = "character"))
  p <- p[match(d$ubigeo, p$ubigeo), ]; stopifnot(!anyNA(p$ubigeo))
  d$tmin_c_ajuste <- p$tmin_c_pisco; d$ppt_mm_ajuste <- p$ppt_mm_pisco
  d$meses_helada_ajuste <- p$meses_helada_pisco; d$spei3_min_ajuste <- p$spei3_min_pisco
}

base <- function(prueba) data.frame(ubigeo = d$ubigeo[prueba], y = d$y[prueba], k = d$k[prueba], n = d$n[prueba],
                                    delta = apply(dist_km[prueba, !prueba, drop = FALSE], 1, min))

if (CASO != "vecinos") {
  f <- as.integer(argumento("pliegue"))
  archivo <- file.path(SAL, sprintf("pliegue%02d.rds", f))
  if (!file.exists(archivo)) {
    set.seed(SEMILLA + f)
    prueba <- pliegue == f
    X <- preparar_predictores(d, d[!prueba, ])
    res <- cbind(base(prueba), pliegue = f)
    if (CASO %in% c("tau0.5", "tau5")) {
      res$p_m4 <- ajustar_maxent(X, prueba)
    } else {
      m1 <- ajustar_inla(X, prueba, VARS, TRUE)
      res <- cbind(res, setNames(as.data.frame(m1), paste0(colnames(m1), "_m1")))
    }
    saveRDS(res, archivo)
    registrar(sprintf("pliegue %d listo", f))
  }
} else {
  desde <- as.integer(argumento("desde")); hasta <- min(as.integer(argumento("hasta")), D)
  archivo <- file.path(SAL, sprintf("distritos_%04d_%04d.rds", desde, hasta))
  if (!file.exists(archivo)) {
    res <- lapply(desde:hasta, function(i) {
      set.seed(SEMILLA + i)
      fuera <- rep(FALSE, D); fuera[c(i, grafo$nbs[[i]])] <- TRUE   # distrito y zona de exclusión
      prueba <- seq_len(D) == i
      X <- preparar_predictores(d, d[!fuera, ])
      m1 <- tryCatch(ajustar_inla(X, prueba, VARS, TRUE, faltan = fuera), error = function(e) { registrar(sprintf("ERROR distrito %d: %s", i, conditionMessage(e))); NULL })
      if (is.null(m1)) return(NULL)
      cbind(base(prueba), vecinos_excluidos = length(grafo$nbs[[i]]), setNames(as.data.frame(m1), paste0(colnames(m1), "_m1")))
    })
    saveRDS(do.call(rbind, res), archivo)
    registrar(sprintf("distritos %d-%d listos", desde, hasta))
  }
}
