# P1 · Paso 9: confusión espacial (plan de análisis, sección 7, punto 5).
# Se ajustan M1 (BYM2) y M2 (sin efectos aleatorios) con los 1 783 distritos y se comparan los
# coeficientes agroclimáticos. Si alguno difiere en más de 2 desviaciones estándar posteriores de M1,
# se ajusta la regresión espacial restringida (Hanks et al., 2015): el efecto BYM2 se restringe a
# ser ortogonal a las covariables (restricciones lineales de INLA sobre la parte combinada del efecto).
# Salida: 04-resultados/aplicacion/confusion_espacial.csv
lineas <- readLines("03-codigo/04_aplicacion_censo.R")
eval(parse(text = lineas[1:(grep("^# --- Un pliegue", lineas) - 1)]))
set.seed(SEMILLA)

X <- preparar_predictores(d, d)
AGRO <- c("tmin", "tmin2", "ppt", "ppt2", "helada", "spei_min")
dd <- cbind(data.frame(k = d$k, n = d$n, id = seq_len(D)), X[, VARIABLES_M1])
rhs <- paste(VARIABLES_M1, collapse = " + ")
bym2 <- "f(id, model = 'bym2', graph = grafo, scale.model = TRUE, hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01)), phi = list(prior = 'pc', param = c(0.5, 2/3)))%s)"
ajustar <- function(formula) {
  for (ci in list(list(), list(cmin = 0))) {
    fit <- tryCatch(inla(as.formula(formula), family = "binomial", Ntrials = n, data = dd, control.inla = ci,
                         num.threads = "1:1", working.directory = tempfile("inla_"), safe = FALSE),
                    error = function(e) e)
    if (!inherits(fit, "error")) return(fit)
  }
  stop(conditionMessage(fit))
}
coef <- function(fit, modelo) {
  s <- fit$summary.fixed[c("(Intercept)", VARIABLES_M1), c("mean", "sd", "0.025quant", "0.975quant")]
  data.frame(modelo = modelo, variable = rownames(s), media = s$mean, de = s$sd, lo = s[, 3], hi = s[, 4])
}
m1 <- ajustar(paste("k ~", rhs, "+", sprintf(bym2, "")))
m2 <- ajustar(paste("k ~", rhs))
c1 <- coef(m1, "M1"); c2 <- coef(m2, "M2")
comp <- data.frame(variable = c1$variable, m1 = c1$media, de_m1 = c1$de, m2 = c2$media, de_m2 = c2$de)
comp$dif_en_de_m1 <- (comp$m1 - comp$m2) / comp$de_m1
comp$dif_en_de_conjunta <- (comp$m1 - comp$m2) / sqrt(comp$de_m1^2 + comp$de_m2^2)
comp$agroclimatica <- comp$variable %in% AGRO
dispara <- any(abs(comp$dif_en_de_m1[comp$agroclimatica]) > 2)
cat("¿Algún coeficiente agroclimático difiere en más de 2 DE de M1?", dispara, "\n")

res <- rbind(c1, c2)
if (dispara) {
  # Restricción A u = 0 con A = [1, X]^T sobre las primeras D componentes del efecto bym2 (u + v combinado).
  A <- cbind(t(cbind(1, as.matrix(X[, VARIABLES_M1]))), matrix(0, length(VARIABLES_M1) + 1, D))
  m1r <- ajustar(paste("k ~", rhs, "+", sprintf(bym2, ", extraconstr = list(A = A, e = rep(0, nrow(A)))")))
  c3 <- coef(m1r, "M1 restringido"); res <- rbind(res, c3)
  comp$m1_rsr <- c3$media; comp$de_m1_rsr <- c3$de
  comp$rsr_menos_m2_en_de <- (comp$m1_rsr - comp$m2) / comp$de_m1_rsr
}
write.csv(comp, file.path(RESULTADOS, "aplicacion", "confusion_espacial.csv"), row.names = FALSE)
write.csv(res, file.path(RESULTADOS, "aplicacion", "coeficientes_m1_m2.csv"), row.names = FALSE)
options(width = 200); print(comp, digits = 3, row.names = FALSE)
