# P1 · Paso 8: resumen de las sensibilidades (plan, sección 7) y de la validación con vecinos
# excluidos (sección 4). Cada fila compara un caso con el análisis principal en las mismas métricas:
# AUC de M1, M3 y M5 (y M4 en los casos de tau), diferencias de H1 con IC bootstrap por bloques,
# y cobertura del conjunto al 90 % de M1, nominal y conformal por bloques.
# Salida: 04-resultados/sensibilidad/resumen_sensibilidad.csv
source("03-codigo/00_comun.R")
set.seed(SEMILLA)
B <- 2000
SENS <- file.path(RESULTADOS, "sensibilidad")

auc <- function(y, p) { r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
cuantil_conformal <- function(s, alfa = 0.10) { m <- length(s); sort(s)[min(m, ceiling((1 - alfa) * (m + 1)))] }
leer <- function(caso) { a <- list.files(file.path(SENS, caso), pattern = "\\.rds$", full.names = TRUE)
  if (!length(a)) return(NULL); x <- do.call(rbind, lapply(a, readRDS)); x[order(x$ubigeo), ] }
dif_boot <- function(y, pa, pb, grupo) {
  obs <- auc(y, pa) - auc(y, pb); g <- split(seq_along(y), grupo)
  d <- replicate(B, { i <- unlist(g[sample(length(g), replace = TRUE)], use.names = FALSE); auc(y[i], pa[i]) - auc(y[i], pb[i]) })
  sprintf("%.3f (%.3f to %.3f)", obs, quantile(d, 0.025), quantile(d, 0.975))
}
cobertura <- function(x) {
  s <- abs(x$u_m1 - 0.5)
  q <- sapply(seq_len(nrow(x)), function(i) cuantil_conformal(s[x$pliegue != x$pliegue[i]]))
  c(nominal = mean(s <= 0.45), bloque = mean(s <= q))
}

base <- do.call(rbind, lapply(list.files(file.path(RESULTADOS, "aplicacion", "pliegues"), full.names = TRUE), readRDS))
base <- base[order(base$ubigeo), ]
tabla <- read.csv(file.path(DATOS, "tabla_distrital_2012.csv"), colClasses = c(ubigeo = "character"))
tabla <- tabla[match(base$ubigeo, tabla$ubigeo), ]

fila <- function(caso, x, y, p1, p3, p5, p4 = NULL, cob = NULL, nota = "") {
  data.frame(caso = caso, distritos = length(y), auc_m1 = auc(y, p1), auc_m3 = auc(y, p3), auc_m5 = auc(y, p5),
             auc_m4 = if (is.null(p4)) NA else auc(y, p4),
             h1_m1_m5 = dif_boot(y, p1, p5, x$pliegue), h1_m3_m5 = dif_boot(y, p3, p5, x$pliegue),
             cob_nominal = if (is.null(cob)) NA else cob[["nominal"]], cob_bloque = if (is.null(cob)) NA else cob[["bloque"]],
             nota = nota)
}

res <- list(fila("principal (tau = 1 %)", base, base$y, base$p_media_m1, base$p_m3, base$p_m5, base$p_m4, cobertura(base)))

for (caso in c("tau0.5", "tau5")) {
  s <- leer(caso); if (is.null(s) || nrow(s) < nrow(base)) next
  y <- tabla[[paste0("y_", caso)]]
  res[[caso]] <- fila(sprintf("umbral tau = %s %%", sub("tau", "", caso)), base, y, base$p_media_m1, base$p_m3, base$p_m5,
                      s$p_m4[match(base$ubigeo, s$ubigeo)], nota = "M1, M3 y M5 no dependen de la etiqueta; M4 reajustado")
}
for (caso in c("altitud_ph", "superficie", "pisco", "iid")) {
  s <- leer(caso); if (is.null(s) || nrow(s) < nrow(base)) next
  s <- s[match(base$ubigeo, s$ubigeo), ]
  nombre <- c(superficie = "respuesta en hectáreas", altitud_ph = "M1 con altitud y pH", pisco = "M1 con clima PISCO", iid = "M2 con efecto iid (en lugar de M1)")[[caso]]
  nota <- c(superficie = "k y n en hectáreas redondeadas", altitud_ph = "", pisco = "Tmín, precipitación, heladas y SPEI-3 mínimo de PISCO, campañas 1991-2012",
            iid = if (caso == "iid") sprintf("columnas _m1 = M2 + iid; M1 - (M2 + iid): %s", dif_boot(base$y, base$p_media_m1, s$p_media_m1, base$pliegue)) else "")[[caso]]
  res[[caso]] <- fila(nombre, s, base$y, s$p_media_m1, base$p_m3, base$p_m5, cob = cobertura(s), nota = nota)
}
urb <- tabla$base_agregacion != "distrito"
b2 <- base[urb, ]
res$urbano <- fila(sprintf("sin los %d distritos urbanos", sum(!urb)), b2, b2$y, b2$p_media_m1, b2$p_m3, b2$p_m5, b2$p_m4, cobertura(b2))

v <- leer("vecinos")
if (!is.null(v)) {
  v$pliegue <- base$pliegue[match(v$ubigeo, base$ubigeo)]
  bv <- base[match(v$ubigeo, base$ubigeo), ]
  sv <- abs(v$u_m1 - 0.5)
  res$vecinos <- data.frame(caso = "dejar fuera el distrito y sus vecinos", distritos = nrow(v), auc_m1 = auc(v$y, v$p_media_m1),
                            auc_m3 = NA, auc_m5 = auc(bv$y, bv$p_m5), auc_m4 = NA, h1_m1_m5 = dif_boot(v$y, v$p_media_m1, bv$p_m5, v$pliegue),
                            h1_m3_m5 = NA, cob_nominal = mean(sv <= 0.45), cob_bloque = NA,
                            nota = sprintf("%d de %d distritos; comparar con bloques de 145 km en los mismos distritos: AUC %.3f, cobertura nominal %.3f",
                                           nrow(v), nrow(base), auc(bv$y, bv$p_media_m1), mean(abs(bv$u_m1 - 0.5) <= 0.45)))
}
res <- do.call(rbind, res)
write.csv(res, file.path(SENS, "resumen_sensibilidad.csv"), row.names = FALSE)
options(width = 220); print(res[, -ncol(res)], digits = 3, row.names = FALSE); cat("\n"); print(res$nota)
