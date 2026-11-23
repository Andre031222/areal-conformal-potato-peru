# Plan de análisis (fijado antes de ajustar los modelos al censo)

Versión 1, 16 de noviembre de 2026. **Aprobado por los autores el 16-nov-2026**, antes de ajustar los modelos competidores al censo.

Este documento fija todas las decisiones **antes de ver un solo resultado**. Cuando el director lo apruebe, se le pone la fecha de aprobación y ya no se cambia; si después algo tiene que cambiar, se agrega como «desviación» con su fecha y su motivo. Así la métrica y las variables no se eligen después de ver qué sale mejor (*p-hacking*).

Datos: `data/tabla_distrital_2012.csv`, 1 783 distritos.

## Decisiones que necesitan el visto bueno del director

Las tres primeras **se apartan del proyecto aprobado** y hay que justificarlas.

| # | Decisión | Proyecto aprobado | Propuesta | Por qué |
| --- | --- | --- | --- | --- |
| D1 | Métrica de H2 | AUC por bloque en la sierra, con Wilcoxon | **Log score binomial por bloque** (devianza predictiva fuera de pliegue), con la misma Wilcoxon pareada | En la sierra solo 37 de 1 138 distritos no tienen papa (τ = 1 %): el AUC por bloque no se puede calcular. Lo que varía en la sierra es **cuánta** papa hay |
| D2 | Altitud en el modelo | La altitud media es variable de «tierra» | **Fuera de BYM2**; queda en el bosque aleatorio y como sensibilidad | Altitud y Tmín correlacionan −0,96 (VIF ≈ 50). Con las dos, los coeficientes no se interpretan. Además, para H3 solo cambia el clima: si la altitud se queda fija, absorbe el efecto de la temperatura |
| D3 | pH en el modelo | pH entre las variables de suelo | **Fuera de BYM2**; queda en el bosque aleatorio y como sensibilidad | El pH de SoilGrids se predice en buena parte con el clima: r = −0,88 con la precipitación y 0,91 con el déficit. Sin pH ni déficit, todos los VIF bajan de 4 |
| D4 | Forma de la temperatura | Lineal | **Término cuadrático en Tmín y en precipitación** | La papa tiene un rango óptimo: no crece ni en la puna más fría ni en la selva caliente. Un logit lineal en Tmín diría «más calor, siempre mejor», que es falso |
| D5 | Umbral τ para el AUC | 1 % | 1 %; sensibilidad con 0,5 % y 5 % | Sin cambio |
| D6 | Abstención | Intervalo creíble al 90 % de BYM2 | **El intervalo de BYM2 se recalibra con conformal por bloques espaciales** (manuscrito, sección 2.5) antes de decidir la abstención, y se reporta la curva riesgo–cobertura | *Spatial Statistics* rechaza los «estudios de caso puros sin desarrollo metodológico» (guía para autores). La recalibración conformal por bloques de un modelo areal bayesiano es el aporte de método; la comparación con la zonificación experta queda como aplicación |

## 1. Respuesta

- Modelo: $k_d \sim \text{Binomial}(n_d, p_d)$, con $k_d$ = `registros_papa` y $n_d$ = `registros_cultivo`.
- AUC (H1 y H1b): contra `y_tau1`, a escala nacional (69 % de distritos con papa: clases equilibradas).

## 2. Predictores (periodo de ajuste: campañas 1991–2012)

| Grupo | Variables | Transformación |
| --- | --- | --- |
| Clima ($\mathbf{x}^A$) | `tmin_c_ajuste` + su cuadrado, `ppt_mm_ajuste` + su cuadrado, `meses_helada_ajuste`, `spei3_min_ajuste` | Estandarizadas (media 0, DE 1); el cuadrado se calcula después de estandarizar |
| Tierra ($\mathbf{x}^A$) | `pendiente_grados`, `carbono_org_gkg`, `arcilla_pct` | Estandarizadas |
| Controles ($\mathbf{x}^C$) | `fraccion_riego`, log(1 + `viaje_min`) | Estandarizadas |
| Fenología (solo H2) | `feno_ndvi_max`, `feno_ndvi_amplitud`, `feno_duracion_verde_meses` | Estandarizadas |

VIF del conjunto sin la fenología: todos < 4 (verificado el 16-nov-2026).
Faltantes: 2 distritos sin suelo y 19 sin fenología. Se imputan con la mediana de la provincia, y la imputación se ajusta **dentro de cada pliegue** con los datos de entrenamiento.
No se usan coordenadas como predictores (Meyer et al., 2019).

## 3. Modelos

| Modelo | Especificación | Software |
| --- | --- | --- |
| **M1, BYM2** (principal) | Binomial, logit, efecto BYM2 sobre la contigüidad (límite compartido) con priors PC: $P(\sigma > 1) = 0{,}01$ y $P(\phi < 0{,}5) = 2/3$ (Simpson et al., 2017) | R-INLA |
| **M2, sin efecto espacial** | M1 sin el efecto BYM2 | R-INLA |
| **M3, bosque aleatorio** | Regresión sobre `fraccion_registros` con pesos $n_d$; todas las variables de M1 más altitud, pH y déficit (D2, D3); 500 árboles; `mtry` y `min.node.size` por validación interna por bloques | ranger |
| **M4, MaxEnt** | Presencias = distritos con `y_tau1 = 1`; fondo = todos los distritos; mismas variables que M1, rasgos lineales y cuadráticos, regularización 1 | maxnet |
| **M5, multicriterio TOPSIS** | Funciones de pertenencia trapezoidales con los rangos de la ficha ECOCROP de *Solanum tuberosum* (FAO) para temperatura, precipitación, pH y pendiente; pesos iguales por criterio | Python o R, sin ajustar con datos |

Los rangos de M5 se tomaron el 16-nov-2026 de la base ECOCROP de la FAO que trae el paquete `Recocrop` 0.4-2 (Hijmans, 2021, doi:10.32614/CRAN.package.Recocrop), registro «Potato, *Solanum tuberosum* L.», **antes de abrir cualquier bloque de prueba**.

**Tabla M5.** Funciones de pertenencia trapezoidales (0 fuera del rango absoluto, 1 dentro del óptimo, lineal entre ambos).

| Criterio | Variable de la tabla distrital | Mínimo absoluto | Óptimo mínimo | Óptimo máximo | Máximo absoluto | Unidad |
| --- | --- | --- | --- | --- | --- | --- |
| Temperatura media | (`tmax_c_ajuste` + `tmin_c_ajuste`) / 2 | 7 | 15 | 25 | 30 | °C |
| Precipitación | `ppt_mm_ajuste` (campaña septiembre–mayo) | 250 | 500 | 800 | 2 000 | mm |
| pH del suelo | `ph` | 4,2 | 5,0 | 6,2 | 8,5 | — |
| Helada | `meses_helada_ajuste` (temperatura letal de ECOCROP: −1 °C) | — | 0 meses | 0 meses | 9 meses | meses |

Pesos iguales (0,25 cada uno), normalización vectorial y TOPSIS; la puntuación de cercanía relativa es la predicción de M5. Notas: ECOCROP da la precipitación anual y aquí se usa la de la campaña, que concentra la mayor parte de la lluvia en la sierra; ECOCROP no tiene pendiente, así que la pendiente sale de M5 (desviación 3).

Islas sin vecinos en la contigüidad: se unen al distrito más cercano para que el grafo sea conexo, y se reporta cuántas fueron.

## 4. Validación

- **Bloques espaciales:** k = 20 pliegues con blockCV. El tamaño del bloque es el rango de la autocorrelación de los residuos de deviance de M2 ajustado con todos los datos, estimado por variograma empírico. Semilla 20261116.
- **BYM2 en el bloque de prueba:** la respuesta se declara faltante y se predice con la distribución predictiva posterior.
- **Validación secundaria** (para mostrar dónde ayuda el efecto espacial): *leave-one-district-out* con zona de exclusión igual a 1 vecino.

## 5. Hipótesis y pruebas (α = 0,05)

| Hipótesis | Estadístico | Prueba | Se sostiene si |
| --- | --- | --- | --- |
| **H1** | ΔAUC = AUC(M1) − AUC(M5) y AUC(M3) − AUC(M5), con las predicciones fuera de pliegue | Intervalo bootstrap remuestreando bloques, B = 2 000; DeLong como complemento | El límite inferior del IC 95 % es > 0,05 |
| **H1b** | AUC(M1) − AUC(M2) y la cobertura empírica del intervalo predictivo al 90 % | Igual que H1, con límite > 0 para el AUC; cobertura reportada con su IC binomial | AUC(M1) ≥ AUC(M2) **y** cobertura de M1 más cerca de 0,90 que la de M2 |
| **H2** (D1) | Log score por bloque en la sierra, M1 con fenología frente a M1 sin ella | Wilcoxon pareado de rangos con signo, unilateral, sobre los bloques con al menos 5 distritos de sierra | p < 0,05 |

Métricas secundarias de todos los modelos: Brier, log score, pendiente e intercepto de calibración e I de Moran de los residuos.

## 6. Incertidumbre y abstención

- **Umbral de clasificación:** índice de Youden calculado en los pliegues de entrenamiento.
- **M1 «no concluyente»:** el intervalo creíble al 90 % de $p_d$ contiene el umbral.
- **M3:** predicción conformal por bloques al 90 %, con la regla equivalente.
- **Fuera del área de aplicabilidad** (índice de disimilitud de CAST; Meyer y Pebesma, 2021): se abstiene.
- **Se reporta** el porcentaje de distritos en abstención y el AUC de los demás.

## 7. Análisis de sensibilidad (se reportan en el suplemento, no cambian la conclusión principal)

1. τ = 0,5 % y 5 %.
2. Respuesta con la fracción de **superficie** en vez de registros.
3. M1 con altitud y con pH (D2, D3).
4. Clima con PISCO en vez de TerraClimate (campañas 1991–2012; PISCOt llega hasta 2016).
5. Regresión espacial restringida (Hanks et al., 2015), si los coeficientes de M1 y M2 difieren en más de 2 errores estándar.
6. Sin los 17 distritos urbanos con `base_agregacion = "distrito"`.

## 8. Lo que se reporta aunque no convenga

- Resultados de las tres hipótesis, se sostengan o no.
- Todas las métricas de todos los modelos (no solo la mejor).
- Las desviaciones de este plan, con su fecha.

## Aprobación

| | Nombre | Fecha | Observaciones |
| --- | --- | --- | --- |
| Director | Dr. Vladimiro Ibañez-Quispe | 16-nov-2026 | Aprobado |
| Coasesor | Ing. Fred Torres-Cruz | | |

## 10. Ampliación metodológica (decidida el 16-nov-2026, antes de ajustar cualquier modelo)

Motivo: una búsqueda de literatura mostró que la recalibración conformal por bloques, sola, ya está casi publicada, así que el trabajo necesitaba un aporte metodológico propio. Se agregan tres piezas. La aplicación al censo y las hipótesis H1, H1b y H2 no cambian.

**10.1 Puntaje de conformidad para conteos binomiales.** Para cada distrito de prueba se calcula el PIT aleatorizado de $k_d$ bajo la predictiva posterior fuera de pliegue de M1: $u_d = F_d(k_d - 1) + V\,f_d(k_d)$, con $V \sim U(0,1)$, y el puntaje $s_d = |u_d - 0{,}5|$. El cuantil conformal $\hat q_b$ se calcula con los distritos de los **otros bloques**. El conjunto de predicción es $\{k : |F_d(k) - 0{,}5| \le \hat q_b\}$, que respeta que $n_d$ varía de 1 a 79 419 registros (mediana 1 635). El puntaje de la versión anterior (sobre la fracción $k_d/n_d$) queda como comparación.

**10.2 Abstención con control de riesgo.** Índice de incertidumbre de cada distrito: el intervalo recalibrado de $p_d$ contiene el umbral, o el índice de disimilitud de la AOA supera $\lambda_{AOA}$. El umbral $\lambda$ se elige como el más permisivo cuyo límite superior de Clopper–Pearson al 95 % del error entre los distritos retenidos, estimado por bloques, no supera $\alpha_r = 0{,}10$ (enfoque de control de riesgo; Bai y Jin, 2026). Se reporta la curva riesgo–cobertura. Comparaciones: (a) solo intervalo bayesiano sin recalibrar; (b) solo máscara AOA (Meyer y Pebesma, 2021); (c) GeoConformal (Lou et al., 2025); (d) conformal por tramos del índice de disimilitud (cinco tramos por cuantiles del DI, al estilo de DI-CQR de Mustafa, 2026, doi:10.2139/ssrn.7563918).

**10.3 Estudio de simulación con verdad conocida.** Grafo real de los 1 783 distritos y los $n_d$ reales del censo. Respuesta simulada con el modelo (1): dos covariables con estructura espacial, $\sigma_b = 1$ y $\phi \in \{0{,}2;\ 0{,}5;\ 0{,}9\}$. Tamaño de bloque $\in \{25;\ 75;\ 145\}$ km. 20 réplicas por escenario y 20 pliegues por réplica (desviación 5). Métricas: cobertura empírica del 90 % del intervalo bayesiano, del conformal por bloques y del conformal aleatorio, en función de la **distancia del distrito al distrito de entrenamiento más cercano**; y el ancho medio de los intervalos. Semilla 20261116.

## Desviaciones

Toda desviación posterior al 16-nov-2026 se anota aquí con fecha y motivo.

| # | Fecha | Qué decía el plan | Qué se hizo | Motivo |
| --- | --- | --- | --- | --- |
| 1 | 16-nov-2026 | Tamaño de bloque = rango del variograma de los residuos de M2 estimado con blockCV | Rango práctico (3a) de un modelo exponencial ajustado por mínimos cuadrados ponderados al correlograma empírico de los mismos residuos: **144,9 km** | blockCV dio 25,5 km, pero el correlograma empírico sigue en r = 0,30 a 10–25 km, 0,18 a 25–50 km y 0,08 a 50–100 km, y se anula hacia los 200 km (`04-resultados/correlograma_residuos_m2.csv`). Con 25 km los bloques serían casi del tamaño de un distrito (lado típico de 14,5 km) y la validación no mediría la extrapolación. El cambio se decidió **antes** de ajustar cualquier modelo de comparación |
| 2 | 16-nov-2026 | Abstención con el intervalo recalibrado (D6) | Se agregan el puntaje binomial, el control de riesgo y la simulación (sección 10) | La revisión de novedad mostró que D6 sola ya está casi publicada. Se decidió **antes** de ajustar cualquier modelo |
| 3 | 16-nov-2026 | M5 con temperatura, precipitación, pH y pendiente | M5 con temperatura, precipitación, pH y helada (Tabla M5) | ECOCROP no trae umbrales de pendiente; la helada sí tiene umbral en ECOCROP (temperatura letal). Fijado antes de abrir los bloques |
| 4 | 16-nov-2026 | Tres comparaciones para la abstención | Se agrega (d): conformal por tramos del DI (Mustafa, 2026) | Es el trabajo más cercano encontrado (leído completo el 16-nov-2026); un revisor pediría la comparación. Decidido antes de ver resultados |
| 5 | 16-nov-2026 | 30 réplicas por escenario | 20 réplicas | Un ajuste de INLA con un hilo tarda 30 s, no 8 s como en el piloto (que usaba varios hilos); 30 réplicas daban unas 7 h. Con 20 réplicas × 1 783 distritos, el error Monte Carlo de la cobertura es de unos ±0,5 puntos porcentuales. Decidido antes de ver resultados |
| 6 | 16-nov-2026 | Ajuste de INLA con su configuración por defecto | Si el optimizador de Newton-Raphson no converge (error 12 de INLA, 29 de 60 réplicas afectadas en la primera corrida), se reintenta con `cmin = 0` y luego en modo clásico | Solo cambia el método numérico de búsqueda de la moda en los pliegues que fallan; el modelo y los priors no cambian. Se registra en `avance.log` cada reintento |
| 7 | 16-nov-2026 | Simulación en el Mac (macOS, M-series) con 9 procesos hijos de `mclapply` | Simulación completa rehecha en una laptop Linux (Arch, 16 núcleos, R 4.6.2, INLA 26.8.7), cada tarea en un proceso de R independiente y cada ajuste en su propia carpeta de trabajo; resultados en `partes_v2/` | El binario de INLA para macOS abortaba en algunos pliegues bajo carga («Abort trap» y error 12), incluso con reintentos; los mismos 20 pliegues que fallaban corrieron todos sin error en Linux con la configuración por defecto. La corrida parcial anterior (31 réplicas, `partes/`) se conserva solo para comparar |
| 8 | 17-nov-2026 | Tres configuraciones numéricas de INLA (desviación 6) | Cinco: se agregan `h = 0.001` con `tolerance = 1e-6` (modo clásico) y, como último recurso, `int.strategy = "eb"`. Los draws de $p_d$ que la interpolación de la marginal deja fuera de (0, 1) se acotan a $[10^{-12}, 1 - 10^{-12}]$. Lo mismo en la aplicación al censo | 22 de 60 réplicas fallaban siempre en el mismo pliegue (determinista, no por carga): 11 por no convergencia en bloques de 25 km y 11 por draws fuera de (0, 1). Las 38 réplicas ya terminadas no cambian. Cada fila guarda qué configuración sirvió (`intento`), para reportarlo y repetir el resumen sin las filas con `eb` |
| 9 | 17-nov-2026 | Sin detalle de implementación para GeoConformal, el umbral de Youden y la elección de $\lambda$ | GeoConformal con núcleo gaussiano de ancho de banda 145 km (el rango práctico de la sección 4) y cuantil ponderado con el punto de prueba en $+\infty$. Umbral de Youden, cuantil conformal y $\lambda$ elegidos con validación cruzada anidada: para el pliegue $f$, el umbral y el cuantil de cada pliegue $g \neq f$ salen de los pliegues distintos de $f$ y $g$; $\lambda$ se busca en la rejilla de cuantiles 5–100 % del DI más $+\infty$. Si ninguna $\lambda$ cumple, el pliegue se abstiene por completo. El AUC (H1, H1b) usa la media posterior de $p_d$ como puntaje | El plan no lo fijaba. Se decidió antes de ver resultados del censo, con solo los pliegues de prueba técnica |
| 10 | 17-nov-2026 | Abstención: el conjunto conformal al 90 % contiene el umbral, o DI > $\lambda$ (sección 10.2) | Se mantiene y se reporta como análisis principal. Se agrega un análisis **exploratorio**: el nivel del conjunto conformal ($1 - a$, con $a \in \{0{,}10; 0{,}15; \dots; 0{,}90\}$; para eso la aplicación al censo guarda los límites del conjunto desde $q = 0{,}02$ (antes desde 0,30), con las mismas semillas; las demás salidas cambian solo en la cuarta cifra decimal (por ejemplo, ΔAUC de M1 − M5: 0,2901 frente a 0,2900), por el muestreo Monte Carlo de INLA) se elige con la misma validación anidada, como el más estrecho cuya cota de Clopper–Pearson del error en los retenidos no supera $\alpha_r = 0{,}10$ | **Decidida después de ver los resultados del censo**: con la regla del plan, $\lambda$ salió $+\infty$ en los 20 pliegues y el error de los retenidos (4,2 %) quedó muy por debajo del 10 %, con 58 % de abstención. La variante usa el control de riesgo para fijar cuánto abstenerse. Se presenta como exploratoria, separada de las hipótesis |
| 11 | 17-nov-2026 | Simulación solo con el modelo bien especificado (sección 10.3) | Se agregan dos escenarios mal especificados, ajustando siempre el BYM2 gaussiano: (a) «colas», con la parte no estructurada del efecto t de Student de 3 g. l. escalada a varianza 1; (b) «no estacionario», con desviación estándar del efecto 1,5 en la sierra y 0,5 fuera de ella. $\phi = 0{,}5$, 20 réplicas por escenario, bloques de 145 km y pliegues aleatorios, mismas semillas que el escenario correcto. Resultados en `partes_mal/` | Con el modelo bien especificado la cobertura nominal ya es 0,90 (sección 10.3), así que la simulación no puede mostrar cuándo hace falta recalibrar. Los dos escenarios imitan fallas plausibles del censo (distritos atípicos y variabilidad mayor en la sierra). Diseño fijado antes de correrlos, pero **decidido después de ver** los resultados del escenario correcto; solo se corrió un piloto técnico de 3 pliegues por escenario |
| 12 | 17-nov-2026 | Conformal por bloques con un solo cuantil para todos los distritos (sección 10.1) | Se agrega, como análisis **exploratorio**, el conformal por bloques **Mondrian por región**: el cuantil del pliegue $f$ se calcula por separado con los distritos de los otros pliegues de la sierra (altitud agrícola media $\ge$ 2300 m) y con los del resto. Se aplica a la simulación (todos los escenarios) y al censo, con las salidas ya guardadas; no se reajusta ningún modelo. Nuevas salidas: `cobertura_por_region*.csv` y `cobertura_censo_region.csv`; la regla `mondrian` se agrega a `cobertura_censo*.csv` sin cambiar las demás filas | **Decidida después de ver los resultados**: en el escenario no estacionario el conformal por bloques dio 0,899 en total pero 0,842 en la sierra y 0,999 fuera de ella, y en el censo 0,951 en la sierra y 0,805 fuera. La garantía conformal es marginal; la partición por región es la forma estándar de obtener cobertura dentro de grupos fijados de antemano (Vovk et al., 2005). La región es la misma definición de sierra de H2, fijada en el plan |
| 13 | 17-nov-2026 | Sin más análisis | Tres análisis pedidos por una revisión interna, todos **después de ver los resultados**: (a) cobertura con el umbral ponderado por bloque de la Proposición 1 (`cub_bloque_pesos`, regla `bloque_pesos`); (b) escenario `no_estacionario_este`, con varianza 1,5 en la mitad oriental de los distritos, para ver la calibración por región cuando las regiones de varianza no coinciden con la sierra; (c) M2 con efecto aleatorio iid (`caso=iid`), porque M2 sin efectos aleatorios no admite sobredispersión y H1b mezclaba «espacial» con «sobredisperso» | La revisión señaló que la proposición cubre un umbral distinto del usado, que la calibración por región en la simulación usaba la misma partición que generó los datos y que M2 no era una comparación justa |
