# P1 · Paso 5: figuras del manuscrito (Elsevier: 190 mm a página completa; PNG y TIFF a 600 dpi y PDF).
# Uso: Rscript 03-codigo/05_figuras.R [fig1 ...]   (sin argumentos hace todas las disponibles)
suppressPackageStartupMessages({library(sf); library(ggplot2); library(patchwork)})
source("03-codigo/00_comun.R")

FIGURAS <- file.path(dirname(RESULTADOS), "figuras")
dir.create(FIGURAS, showWarnings = FALSE, recursive = TRUE)
AZUL <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")  # secuencial, un tono
TINTA <- "#0b0b0b"; TINTA2 <- "#52514e"; GRIS <- "#d9d8d3"
tema <- theme_void(base_family = "Latin Modern Roman", base_size = 8.4) +
  theme(legend.position = "bottom", legend.title = element_text(size = 8.4, colour = TINTA),
        legend.text = element_text(size = 7.8, colour = TINTA2),
        plot.tag = element_text(size = 9.6, face = "bold", colour = TINTA),
        plot.margin = margin(2, 2, 2, 2))

guardar <- function(g, nombre, ancho_mm, alto_mm, g_revista = g) {
  base <- file.path(FIGURAS, nombre)
  # PNG (README y presentaciones) con títulos interpretativos; PDF y TIFF (revista) sin ellos:
  # en el artículo el mensaje va en la leyenda de la figura.
  ggsave(paste0(base, ".png"), g, width = ancho_mm, height = alto_mm, units = "mm", dpi = 600, bg = "white")
  sin <- theme(plot.title = element_blank(), plot.subtitle = element_blank())
  g <- g_revista
  g <- if (inherits(g, "patchwork")) g & sin else g + sin
  ggsave(paste0(base, ".pdf"), g, width = ancho_mm, height = alto_mm, units = "mm", device = cairo_pdf)
  ggsave(paste0(base, ".tiff"), g, width = ancho_mm, height = alto_mm, units = "mm", dpi = 600,
         compression = "lzw", bg = "white")
  cat("guardada:", base, "(.pdf .png .tiff)\n")
}

# Relieve sombreado (data-raw/relieve_sombreado.R) como capa semitransparente sobre los colores:
# solo oscurece las laderas en sombra, de modo que los colores de la leyenda se siguen leyendo.
capa_relieve <- local({
  r <- NULL
  function() {
    archivo <- file.path(DATOS, "relieve_peru.tif")
    if (!file.exists(archivo)) return(NULL)
    if (is.null(r)) {
      t <- terra::rast(archivo)
      r <<- setNames(as.data.frame(t, xy = TRUE, na.rm = TRUE), c("x", "y", "s"))
      r$a <<- pmin(pmax((0.80 - r$s) * 0.75, 0), 0.42)
    }
    list(geom_raster(data = r, aes(x, y, alpha = a), fill = "black", inherit.aes = FALSE, show.legend = FALSE),
         scale_alpha_identity())
  }
})

# Globo de ubicación: proyección ortográfica centrada en el Perú; solo el hemisferio visible.
globo <- function() {
  orto <- "+proj=ortho +lat_0=-10 +lon_0=-68 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m"
  casquete <- st_buffer(st_sfc(st_point(c(-68, -10)), crs = 4326), 9.6e6)   # hemisferio visible (s2, metros)
  visible <- function(x) { x <- suppressWarnings(st_intersection(st_make_valid(x), casquete)); st_transform(x[!st_is_empty(x), ], orto) }
  mundo <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf")[, "admin"]
  mundo <- mundo[!mundo$admin %in% c("Antarctica", "Russia", "Fiji"), ]   # cruzan la línea de cambio de fecha; no se ven
  mundo <- mundo[s2::s2_is_valid(mundo), ]
  m <- visible(mundo); peru <- m[m$admin == "Peru", ]
  rejilla <- visible(st_graticule(lon = seq(-180, 150, 30), lat = seq(-60, 60, 30), crs = 4326, ndiscr = 200)[, 0])
  disco <- st_buffer(st_sfc(st_point(c(0, 0)), crs = orto), 6371000 * sin(9.6e6 / 6371000))
  ggplot() + geom_sf(data = disco, fill = "#eef4fb", colour = TINTA2, linewidth = 0.25) +
    geom_sf(data = rejilla, colour = "white", linewidth = 0.15) +
    geom_sf(data = m, fill = "#d9d8d3", colour = "white", linewidth = 0.08) +
    geom_sf(data = peru, fill = "#184f95", colour = NA) +
    coord_sf(crs = orto, expand = FALSE) + theme_void()
}

fig1 <- function() {
  d <- leer_distritos()
  pl <- read.csv(file.path(RESULTADOS, "pliegues.csv"), colClasses = c(ubigeo = "character"))
  d$pliegue <- pl$pliegue[match(d$ubigeo, pl$ubigeo)]
  d$fraccion <- d$k / d$n
  d <- st_transform(d, CRS_METRICO)
  bloques <- aggregate(d["pliegue"], by = list(pliegue = d$pliegue), FUN = function(x) x[1])
  bloques <- suppressWarnings(st_cast(st_make_valid(bloques), "POLYGON"))
  pais <- st_union(d)
  departamentos <- aggregate(d["departamento"], by = list(dep = d$departamento), FUN = function(x) x[1])

  # Ciudades de referencia (coordenadas de su plaza principal) y el lago Titicaca.
  ciudades <- st_transform(st_as_sf(data.frame(
    nombre = c("Lima", "Cajamarca", "Huancayo", "Cusco", "Arequipa", "Puno", "Iquitos"),
    lon = c(-77.043, -78.517, -75.205, -71.968, -71.537, -70.020, -73.254),
    lat = c(-12.046, -7.162, -12.065, -13.532, -16.409, -15.840, -3.744),
    dx = c(-1, 1, -1, 1, -1, 1, 1)), coords = c("lon", "lat"), crs = 4326), CRS_METRICO)
  ciudades$x <- st_coordinates(ciudades)[, 1]; ciudades$y <- st_coordinates(ciudades)[, 2]
  titicaca <- st_transform(st_sfc(st_point(c(-69.35, -15.75)), crs = 4326), CRS_METRICO)

  # Clases (niveles) de la fracción de registros con papa.
  cortes <- c(0, 0.01, 0.05, 0.10, 0.20, 0.30, 0.45, 1)
  etiquetas <- c("Below 0.01", "0.01 to 0.05", "0.05 to 0.10", "0.10 to 0.20", "0.20 to 0.30", "0.30 to 0.45", "0.45 or more")
  d$clase <- cut(d$fraccion, cortes, labels = etiquetas, right = FALSE, include.lowest = TRUE)

  # Sin rótulo para el lago Titicaca: se superponía con Puno (el lago se ve en el mapa).
  lugares <- data.frame(nombre = ciudades$nombre, x = ciudades$x, y = ciudades$y, ciudad = TRUE)

  capas_lugar <- function(g, con_departamentos = TRUE) {
    g <- g + capa_relieve()
    if (con_departamentos) g <- g + geom_sf(data = departamentos, fill = NA, colour = "white", linewidth = 0.25)
    g + geom_sf(data = pais, fill = NA, colour = TINTA2, linewidth = 0.3) +
      geom_point(data = lugares[lugares$ciudad, ], aes(x, y), shape = 21, size = 1.4, stroke = 0.45,
                 fill = "white", colour = TINTA) +
      ggrepel::geom_text_repel(data = lugares, aes(x, y, label = nombre, fontface = ifelse(ciudad, "plain", "italic")),
                               size = 2.52, colour = TINTA, family = "Latin Modern Roman", bg.colour = "white", bg.r = 0.12,
                               min.segment.length = 0.3, segment.size = 0.2, segment.colour = TINTA2,
                               box.padding = 0.35, point.padding = 0.2, seed = SEMILLA, max.overlaps = Inf) +
      ggspatial::annotation_scale(location = "bl", width_hint = 0.25, height = unit(1.2, "mm"),
                                  text_cex = 0.55, line_width = 0.3, style = "ticks") +
      ggspatial::annotation_north_arrow(location = "tr", height = unit(7, "mm"), width = unit(5, "mm"),
                                        style = ggspatial::north_arrow_orienteering(text_size = 5)) +
      coord_sf(crs = CRS_METRICO, datum = 4326, expand = FALSE) +
      scale_x_continuous(breaks = c(-80, -76, -72)) + scale_y_continuous(breaks = c(-4, -8, -12, -16)) +
      tema + theme(axis.text = element_text(size = 6.6, colour = TINTA2),
                   axis.ticks = element_line(colour = TINTA2, linewidth = 0.2), axis.ticks.length = unit(0.8, "mm"),
                   panel.border = element_rect(fill = NA, colour = TINTA2, linewidth = 0.3))
  }

  a <- capas_lugar(ggplot() + geom_sf(data = d, aes(fill = clase), colour = NA) +
    scale_fill_manual(values = AZUL, drop = FALSE, name = "Share of crop records with potato, 2012") +
    guides(fill = guide_legend(title.position = "top", nrow = 2, byrow = TRUE,
                               keywidth = unit(4, "mm"), keyheight = unit(2.5, "mm")))) + labs(tag = "a")

  b <- capas_lugar(ggplot() +
    geom_sf(data = d, aes(fill = factor(sierra, levels = c(0, 1),
                                        labels = c("Below 2300 m", "2300 m or more (highlands)"))), colour = NA) +
    geom_sf(data = bloques, fill = NA, colour = TINTA, linewidth = 0.35) +
    scale_fill_manual(values = c(GRIS, "#86b6ef"), name = "Mean elevation of agricultural land\n(black lines: 145 km validation blocks)") +
    guides(fill = guide_legend(title.position = "top", nrow = 1, keywidth = unit(4, "mm"),
                               keyheight = unit(2.5, "mm"))), con_departamentos = FALSE) + labs(tag = "b")

  # El globo de ubicación va solo en el PNG (README): Elsevier no admite en el artículo mapas de
  # ubicación más grandes que el área de estudio.
  con_globo <- a + inset_element(globo(), left = 0.01, bottom = 0.075, right = 0.27, top = 0.325, align_to = "panel")
  guardar(con_globo + b, "Fig1_area_estudio", 190, 135, g_revista = a + b)
}

# --- Figuras 2 a 5: resultados ----------------------------------------------------------------
# Colores categóricos de Okabe e Ito (seguros para daltonismo), en orden fijo por regla.
COL_REGLA <- c("Bayesian 90% set" = "#D55E00", "Block conformal (PIT score)" = "#0072B2",
               "Conformal on the fraction" = "#E69F00", "GeoConformal" = "#009E73",
               "Conformal by DI bins" = "#CC79A7", "Block conformal by region" = "#56B4E9")
NOMBRE_REGLA <- c(nominal = "Bayesian 90% set", bloque = "Block conformal (PIT score)",
                  fraccion = "Conformal on the fraction", geoconformal = "GeoConformal",
                  tramos_di = "Conformal by DI bins", mondrian = "Block conformal by region")
TRAMOS <- c("[0,10)" = "0-10", "[10,25)" = "10-25", "[25,50)" = "25-50", "[50,100)" = "50-100", "[100,Inf)" = "100 or more")
tema_ejes <- theme_minimal(base_family = "Latin Modern Roman", base_size = 8.4) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#ebeae6", linewidth = 0.25),
        axis.line.x = element_line(colour = TINTA2, linewidth = 0.3), axis.ticks.x = element_line(colour = TINTA2, linewidth = 0.25),
        axis.text = element_text(size = 7.2, colour = TINTA2), axis.title = element_text(size = 8.16, colour = TINTA),
        strip.text = element_text(size = 8.16, colour = TINTA, hjust = 0), legend.position = "bottom",
        legend.title = element_blank(), legend.text = element_text(size = 7.56, colour = TINTA2),
        plot.tag = element_text(size = 9.6, face = "bold", colour = TINTA), plot.margin = margin(3, 4, 3, 3),
        plot.title = element_text(size = 8.64, face = "bold", colour = TINTA, margin = margin(b = 1)),
        plot.subtitle = element_text(size = 7.68, colour = TINTA2, margin = margin(b = 3)), plot.title.position = "plot")
linea_objetivo <- function(y = 0.90) geom_hline(yintercept = y, linetype = "22", colour = TINTA2, linewidth = 0.3)
Y_COBERTURA <- "Share of districts whose census value\nfell inside the 90% set"

fig2 <- function() {
  # Simulación: cobertura frente a la distancia, modelo correcto (tres phi) y dos modelos mal especificados.
  lee <- function(a) read.csv(file.path(RESULTADOS, "simulacion", a))
  sim <- rbind(lee("cobertura_por_distancia.csv"), lee("cobertura_por_distancia_mal.csv"))
  sim <- sim[sim$diseno == "bloque145", ]
  sim <- rbind(data.frame(phi = sim$phi, tramo = sim$tramo_delta, regla = "nominal", cobertura = sim$cub_nominal),
               data.frame(phi = sim$phi, tramo = sim$tramo_delta, regla = "bloque", cobertura = sim$cub_bloque))
  etiquetas <- c("0.2" = "phi == 0.2", "0.5" = "phi == 0.5", "0.9" = "phi == 0.9",
                 colas = '"Heavy-tailed effect"', no_estacionario = '"Non-stationary, highlands"',
                 no_estacionario_este = '"Non-stationary, east-west"')
  sim$panel <- factor(etiquetas[sim$phi], levels = etiquetas)
  sim$grupo <- factor(ifelse(sim$phi %in% c("colas", "no_estacionario", "no_estacionario_este"), "mal", "bien"), levels = c("bien", "mal"))
  sim$tramo <- factor(TRAMOS[sim$tramo], levels = TRAMOS)
  sim$regla <- factor(NOMBRE_REGLA[sim$regla], levels = names(COL_REGLA))
  usadas <- c("Bayesian 90% set", "Block conformal (PIT score)")
  esc <- scale_colour_manual(values = COL_REGLA[usadas], limits = usadas)
  dodge <- position_dodge(width = 0.55)
  panel <- function(datos, tag, titulo, subtitulo, eje_y = TRUE) {
    ggplot(datos, aes(tramo, cobertura, colour = regla)) + linea_objetivo() +
      geom_point(size = 1.4, position = dodge) + facet_wrap(~panel, nrow = 1, labeller = label_parsed) + esc +
      scale_y_continuous(limits = c(0.75, 1.00), breaks = seq(0.75, 1, 0.05), labels = function(v) paste0(100 * v, "%")) +
      labs(x = "Distance to the nearest training district (km)", y = if (eje_y) "Share of held-out districts whose simulated\ncount fell inside the 90% set" else NULL, tag = tag,
           title = titulo, subtitle = subtitulo) +
      tema_ejes + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  }
  a <- panel(sim[sim$grupo == "bien", ], "a", "Simulated data, model correctly specified",
             "Both sets stay on the 90% target (dashed line) at every distance")
  b <- panel(sim[sim$grupo == "mal", ], "b", "Simulated data, model misspecified (spatial share 0.5)",
             "The Bayesian set over-covers with heavy tails and under-covers with a non-stationary effect; calibration corrects both")
  g <- (a / b + plot_layout(guides = "collect")) & theme(legend.position = "bottom")
  guardar(g, "Fig2_cobertura_simulacion", 190, 150)
}

fig3 <- function() {
  # Censo: cobertura por distancia (a) y cobertura por región en simulación y censo (b).
  cen <- read.csv(file.path(RESULTADOS, "aplicacion", "cobertura_censo_distancia.csv"))
  cen <- cen[cen$modelo == "M1" & cen$regla %in% c("nominal", "bloque", "geoconformal", "tramos_di", "mondrian"), ]
  cen$tramo <- factor(TRAMOS[cen$tramo], levels = TRAMOS)
  cen$regla <- factor(NOMBRE_REGLA[cen$regla], levels = names(COL_REGLA))
  usadas <- c("Bayesian 90% set", "Block conformal (PIT score)", "Block conformal by region", "GeoConformal", "Conformal by DI bins")
  esc <- scale_colour_manual(values = COL_REGLA[usadas], limits = usadas)
  dodge <- position_dodge(width = 0.65)
  n_tramo <- unique(cen[, c("tramo", "n")])
  a <- ggplot(cen, aes(tramo, cobertura, colour = regla)) + linea_objetivo() +
    geom_linerange(aes(ymin = ic_lo, ymax = ic_hi), position = dodge, linewidth = 0.4, show.legend = FALSE) +
    geom_point(size = 1.3, position = dodge) + esc +
    geom_text(data = n_tramo, aes(tramo, 0.31, label = paste0("n = ", n)), inherit.aes = FALSE,
              size = 2.28, colour = TINTA2, family = "Latin Modern Roman") +
    scale_y_continuous(limits = c(0.30, 1.00), breaks = seq(0.3, 1, 0.1), labels = function(v) paste0(100 * v, "%")) +
    annotate("text", x = 4.55, y = 0.42, label = "only 14 districts:\nwide intervals", hjust = 1, size = 2.28,
             colour = TINTA2, family = "Latin Modern Roman", lineheight = 0.9) +
    labs(x = "Distance to the nearest training district (km)", y = paste(Y_COBERTURA, "(95% CI)"), tag = "a",
         title = "Census 2012: coverage by distance",
         subtitle = "The Bayesian set falls short of 90%; conformal\nsets fix it except beyond 100 km (14 districts)") +
    tema_ejes + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  # Panel b: por región. Simulación (phi = 0,5 y los dos mal especificados) con IC entre réplicas; censo con IC de Clopper-Pearson.
  lee <- function(a) read.csv(file.path(RESULTADOS, "simulacion", a))
  sr <- rbind(lee("cobertura_por_region.csv"), lee("cobertura_por_region_mal.csv"))
  sr <- sr[sr$diseno == "bloque145" & sr$phi %in% c("0.5", "colas", "no_estacionario"), ]
  larga <- function(r, pre) data.frame(caso = sr$phi, region = sr$region, regla = r, cobertura = sr[[pre]],
                                       ic_lo = sr[[paste0(pre, "_lo")]], ic_hi = sr[[paste0(pre, "_hi")]])
  sr <- rbind(larga("nominal", "cub_nominal"), larga("bloque", "cub_bloque"), larga("mondrian", "cub_mondrian"))
  cr <- read.csv(file.path(RESULTADOS, "aplicacion", "cobertura_censo_region.csv"))
  cr <- cr[cr$modelo == "M1" & cr$regla %in% c("nominal", "bloque", "mondrian"), ]
  cr <- data.frame(caso = "censo", region = cr$region, regla = cr$regla, cobertura = cr$cobertura, ic_lo = cr$ic_lo, ic_hi = cr$ic_hi)
  mi <- read.csv(file.path(RESULTADOS, "simulacion", "cobertura_por_mitad_mal.csv"))
  mi <- mi[mi$diseno == "bloque145" & mi$phi == "no_estacionario_este", ]
  me <- function(r, col) data.frame(caso = "este", region = mi$mitad, regla = r, cobertura = mi[[col]], ic_lo = NA, ic_hi = NA)
  r <- rbind(sr, cr, me("nominal", "cub_nominal"), me("bloque", "cub_bloque"), me("mondrian", "cub_mondrian"))
  casos <- c("0.5" = "Simulation,\ncorrect model", colas = "Simulation,\nheavy tails", no_estacionario = "Simulation,\nnon-stationary",
             este = "Simulation,\neast-west variance", censo = "Census 2012\n(n = 1138 and 645)")
  r$caso <- factor(casos[r$caso], levels = casos)
  r$region <- factor(c(sierra = "Highlands", resto = "Rest", este = "East", oeste = "West")[r$region], levels = c("Highlands", "Rest", "East", "West"))
  r$regla <- factor(NOMBRE_REGLA[r$regla], levels = names(COL_REGLA))
  b <- ggplot(r, aes(region, cobertura, colour = regla)) + linea_objetivo() +
    geom_linerange(aes(ymin = ic_lo, ymax = ic_hi), position = dodge, linewidth = 0.4, show.legend = FALSE) +
    geom_point(size = 1.3, position = dodge) + esc + facet_wrap(~caso, nrow = 1, scales = "free_x") +
    scale_y_continuous(limits = c(0.68, 1.005), breaks = seq(0.70, 1, 0.05), labels = function(v) paste0(100 * v, "%")) +
    labs(x = "Region (highlands: mean agricultural elevation of 2300 m or more; east and west: halves of the districts)", y = NULL, tag = "b",
         title = "Coverage within regions",
         subtitle = "One threshold reaches 90% on average but not within regions of unequal error variance. Calibrating highlands and rest\nseparately fixes it when those are the regions that differ, but not when the variance differs between east and west") +
    tema_ejes + theme(panel.spacing.x = unit(3, "mm")) + guides(colour = "none")
  g <- (a + b + plot_layout(widths = c(1.1, 3.4), guides = "collect")) & theme(legend.position = "bottom")
  guardar(g, "Fig3_cobertura_censo_region", 190, 100)
}

fig4 <- function() {
  # Cobertura frente a ancho, para todos los distritos y para los alejados del entrenamiento (50 km o más).
  gl <- read.csv(file.path(RESULTADOS, "aplicacion", "cobertura_censo.csv")); gl <- gl[gl$modelo == "M1", ]
  cd <- read.csv(file.path(RESULTADOS, "aplicacion", "cobertura_censo_distancia.csv"))
  cd <- cd[cd$modelo == "M1" & cd$tramo_delta %in% c("[50,100)", "[100,Inf)"), ]
  lejos <- do.call(rbind, lapply(split(cd, cd$regla), function(g) {
    k <- round(sum(g$cobertura * g$n)); n <- sum(g$n); ci <- binom.test(k, n)$conf.int
    data.frame(regla = g$regla[1], n = n, cobertura = k / n, ic_lo = ci[1], ic_hi = ci[2], ancho_medio = sum(g$ancho_medio * g$n) / n)
  }))
  todos <- gl[, c("regla", "n", "cobertura", "ic_lo", "ic_hi", "ancho_medio")]
  x <- rbind(cbind(todos, panel = sprintf("All districts (n = %d)", todos$n[1])),
             cbind(lejos, panel = sprintf("Districts 50 km or more from any training district (n = %d)", lejos$n[1])))
  x$panel <- factor(x$panel, levels = unique(x$panel))
  x$regla <- factor(NOMBRE_REGLA[x$regla], levels = names(COL_REGLA))
  x$etiqueta <- sprintf("%s\n%.1f%%", x$regla, 100 * x$cobertura)
  # Posición de cada etiqueta, fijada a mano (panel 1: todos; panel 2: alejados).
  pos <- data.frame(p = rep(1:2, each = 6), regla = rep(names(COL_REGLA), 2),
                    lx = c(0.300, 0.392, 0.262, 0.345, 0.395, 0.300,   0.290, 0.335, 0.292, 0.350, 0.360, 0.262),
                    ly = c(0.826, 0.933, 0.872, 0.936, 0.868, 0.928,   0.818, 0.933, 0.868, 0.866, 0.906, 0.936))
  x$p <- as.integer(x$panel); x <- merge(x, pos, by.x = c("p", "regla"), by.y = c("p", "regla"))
  x$regla <- factor(x$regla, levels = names(COL_REGLA))
  x$xd <- x$ancho_medio + c(-0.002, 0, 0, 0.002, 0.004, 0)[as.integer(x$regla)] * (x$p == 2)
  g <- ggplot(x, aes(ancho_medio, cobertura, colour = regla)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.88, ymax = 0.92, fill = "#eaf4ee") +
    geom_hline(yintercept = 0.90, linetype = "22", colour = TINTA2, linewidth = 0.3) +
    geom_segment(aes(x = lx, y = ly, xend = xd, yend = cobertura), linewidth = 0.2, alpha = 0.7) +
    geom_linerange(aes(x = xd, ymin = ic_lo, ymax = ic_hi), linewidth = 0.45) +
    geom_point(aes(x = xd), size = 2.2) +
    geom_label(aes(lx, ly, label = etiqueta), size = 2.4, family = "Latin Modern Roman", lineheight = 0.9, fill = "white",
               linewidth = 0, label.padding = unit(0.5, "mm"), show.legend = FALSE) +
    geom_text(data = data.frame(panel = factor(levels(x$panel), levels = levels(x$panel)), xt = c(0.243, 0.345), yt = 0.9185),
              aes(xt, yt, label = "On target: 88-92%"), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.28,
              colour = "#2f6b45", family = "Latin Modern Roman") +
    facet_wrap(~panel, nrow = 1) +
    scale_colour_manual(values = COL_REGLA, guide = "none") +
    scale_x_continuous(limits = c(0.24, 0.42), breaks = seq(0.25, 0.40, 0.05)) +
    scale_y_continuous(limits = c(0.79, 0.94), breaks = seq(0.80, 0.94, 0.02), labels = function(v) paste0(100 * v, "%")) +
    labs(x = "Mean width of the set (share of crop records with potato); narrower is more precise",
         y = "Share of districts whose census value\nfell inside the 90% set (95% CI)",
         title = "Which set keeps its promise, and at what price in width?",
         subtitle = "Over all districts every conformal set reaches 90%. Far from the training data the block-conformal sets, GeoConformal and the DI bins\nstay on target; the fraction score and the Bayesian set fall short. Calibrating by region gives the same coverage with narrower sets.") +
    tema_ejes + theme(panel.grid.major.x = element_line(colour = "#ebeae6", linewidth = 0.25),
                      panel.spacing.x = unit(6, "mm"), strip.text = element_text(size = 8.16, face = "bold", colour = TINTA, hjust = 0))
  guardar(g, "Fig4_comparacion_conformal", 190, 95)
}

fig5 <- function() {
  ca <- read.csv(file.path(RESULTADOS, "aplicacion", "riesgo_cobertura_alfa.csv"))
  cd <- read.csv(file.path(RESULTADOS, "aplicacion", "riesgo_cobertura.csv"))
  cd <- cd[cd$metodo == "solo DI" & cd$fraccion_retenida >= 0.2, ]
  ab <- read.csv(file.path(RESULTADOS, "aplicacion", "abstencion.csv"))
  ab$letra <- ifelse(grepl("^[a-g]\\.", ab$regla), substr(ab$regla, 1, 1), "none")
  nombres <- c(a = "Bayesian 90% set", b = "Area-of-applicability mask", c = "GeoConformal", d = "Conformal by DI bins",
               e = "Block conformal 90%", f = "Block conformal + DI (plan)", g = "Risk-calibrated level", none = "No abstention")
  ab$nombre <- nombres[ab$letra]
  ab <- ab[ab$letra != "f", ]   # f coincide con e (lambda = infinito en los 20 bloques); se indica en la leyenda
  ab$nombre[ab$letra == "e"] <- "Block conformal 90% (= plan rule)"
  ab$nombre <- ifelse(ab$letra %in% c("e", "g", "none"),
                      sprintf("%s\n%.0f%% predicted, %.1f%% wrong", ab$nombre, 100 * ab$fraccion_retenida, 100 * ab$riesgo), ab$nombre)
  # Posición de cada etiqueta (fijada a mano para que no se encimen).
  pos <- data.frame(letra = c("d", "e", "c", "a", "b", "g", "none"),
                    lx = c(0.30, 0.29, 0.47, 0.64, 0.80, 0.70, 0.90), ly = c(0.024, 0.058, 0.074, 0.034, 0.064, 0.120, 0.124))
  ab <- merge(ab, pos, by = "letra")
  curvas <- rbind(data.frame(x = ca$fraccion_retenida, y = ca$riesgo, curva = "Abstain where the conformal set is unsure (set level varied from 90% to 10%)"),
                  data.frame(x = cd$fraccion_retenida, y = cd$riesgo, curva = "Abstain where conditions are new to the model (dissimilarity threshold varied)"))
  g <- ggplot() +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 0.10, fill = "#eaf4ee") +
    geom_hline(yintercept = 0.10, linetype = "22", colour = TINTA2, linewidth = 0.3) +
    annotate("text", x = 0.205, y = 0.004, label = "Acceptable error: 10% or less", hjust = 0, vjust = 0, size = 2.4, colour = "#2f6b45", family = "Latin Modern Roman") +
    geom_line(data = curvas, aes(x, y, linetype = curva), colour = TINTA2, linewidth = 0.45) +
    geom_linerange(data = ab, aes(x = fraccion_retenida, ymin = riesgo_ic_lo, ymax = riesgo_ic_hi), colour = TINTA, linewidth = 0.35) +
    geom_point(data = ab, aes(fraccion_retenida, riesgo), shape = 21, size = 2, fill = "white", colour = TINTA, stroke = 0.5) +
    geom_segment(data = ab, aes(x = lx, y = ly, xend = fraccion_retenida, yend = riesgo), colour = TINTA2, linewidth = 0.2) +
    geom_label(data = ab, aes(lx, ly, label = nombre), size = 2.4, colour = TINTA, family = "Latin Modern Roman",
               linewidth = 0, label.padding = unit(0.6, "mm"), fill = "white") +
    scale_linetype_manual(values = setNames(c("solid", "11"), unique(curvas$curva))) +
    scale_x_continuous(limits = c(0.2, 1), breaks = seq(0.2, 1, 0.2), labels = function(v) paste0(100 * v, "%")) +
    scale_y_continuous(limits = c(0, 0.135), breaks = seq(0, 0.12, 0.02), labels = function(v) paste0(100 * v, "%")) +
    labs(x = "Districts that receive a prediction (the rest are left as inconclusive)",
         y = "Share of predicted districts misclassified\n(95% CI)",
         title = "Abstaining buys accuracy at the cost of coverage",
         subtitle = "Best rules sit low (few errors) and to the right (many districts predicted)") +
    tema_ejes + theme(panel.grid.major.x = element_line(colour = "#ebeae6", linewidth = 0.25),
                      legend.key.width = unit(6, "mm"), legend.key.height = unit(2.5, "mm"),
                      legend.spacing.y = unit(0, "mm")) + guides(linetype = guide_legend(nrow = 2))
  guardar(g, "Fig5_riesgo_cobertura", 140, 110)
}

fig6 <- function() {
  d <- st_transform(leer_distritos(), CRS_METRICO)
  dec <- read.csv(file.path(RESULTADOS, "aplicacion", "decision_por_distrito.csv"), colClasses = c(ubigeo = "character"))
  d <- merge(d, dec[, c("ubigeo", "pliegue", "error", "nc_riesgo")], by = "ubigeo")
  d$decision <- factor(ifelse(d$nc_riesgo, "Inconclusive", ifelse(d$error == 1, "Predicted, misclassified", "Predicted, correct")),
                       levels = c("Predicted, correct", "Predicted, misclassified", "Inconclusive"))
  rb <- read.csv(file.path(RESULTADOS, "aplicacion", "riesgo_por_bloque.csv"))
  rb <- rb[startsWith(rb$regla, "g."), ]
  bloques <- aggregate(d["pliegue"], by = list(bloque = d$pliegue), FUN = function(x) x[1])
  bloques <- st_make_valid(bloques)
  bloques$riesgo <- rb$error[match(bloques$bloque, rb$pliegue)]
  bloques$clase <- cut(bloques$riesgo, c(0, 0.05, 0.10, 0.20, 1), right = TRUE, include.lowest = TRUE,
                       labels = c("5% or less", "5% to 10%", "10% to 20%", "More than 20%"))
  pais <- st_union(d)
  base_mapa <- function(g) g + capa_relieve() + geom_sf(data = pais, fill = NA, colour = TINTA2, linewidth = 0.3) +
    ggspatial::annotation_scale(location = "bl", width_hint = 0.25, height = unit(1.2, "mm"), text_cex = 0.55,
                                line_width = 0.3, style = "ticks") +
    ggspatial::annotation_north_arrow(location = "tr", height = unit(7, "mm"), width = unit(5, "mm"),
                                      style = ggspatial::north_arrow_orienteering(text_size = 5)) +
    coord_sf(crs = CRS_METRICO, datum = 4326, expand = FALSE) +
    scale_x_continuous(breaks = c(-80, -76, -72)) + scale_y_continuous(breaks = c(-4, -8, -12, -16)) +
    tema + theme(axis.text = element_text(size = 6.6, colour = TINTA2), axis.ticks = element_line(colour = TINTA2, linewidth = 0.2),
                 axis.ticks.length = unit(0.8, "mm"), panel.border = element_rect(fill = NA, colour = TINTA2, linewidth = 0.3))
  leyenda <- guides(fill = guide_legend(title.position = "top", nrow = 2, byrow = TRUE, keywidth = unit(4, "mm"), keyheight = unit(2.5, "mm")))
  a <- base_mapa(ggplot() + geom_sf(data = d, aes(fill = decision), colour = NA) +
    geom_sf(data = bloques, fill = NA, colour = TINTA, linewidth = 0.25) +
    scale_fill_manual(values = c("#c9dcef", "#D55E00", "#8a8984"),
                      name = "Out-of-fold decision, risk-calibrated rule\n(black lines: the 20 validation folds of 145 km blocks)") + leyenda) +
    labs(tag = "a", title = "Where the rule predicts correctly, errs or abstains") + theme(plot.title = element_text(size = 8.64, face = "bold", colour = TINTA))
  b <- base_mapa(ggplot() + geom_sf(data = bloques, aes(fill = clase), colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = c("#fde3cf", "#f7a96b", "#e1682b", "#a63603"), drop = FALSE,
                      name = "Share of predicted districts misclassified, by validation\nfold (each fold groups several 145 km blocks)") + leyenda) +
    labs(tag = "b", title = sprintf("The error is uneven: %.0f%% overall, up to %.0f%% in some folds", 100 * mean(d$error[!d$nc_riesgo]), 100 * max(bloques$riesgo))) + theme(plot.title = element_text(size = 8.64, face = "bold", colour = TINTA))
  guardar(a + b, "Fig6_mapa_abstencion", 190, 135)
}


# Figura 3 (ggplot2 + patchwork): cobertura y ancho del intervalo al 90 % de M1 en el censo, por estrato.
fig3_gg <- function() {
  OI <- c(nominal = "#E69F00", bloque = "#56B4E9", mondrian = "#0072B2")
  NOM <- c(nominal = "Bayesian interval (uncalibrated)", bloque = "Block-conformal calibration",
           mondrian = "Region-wise block-conformal calibration")
  FORMA <- c(nominal = 21, bloque = 24, mondrian = 22)
  gl <- read.csv(file.path(RESULTADOS, "aplicacion", "cobertura_censo.csv")); gl <- gl[gl$modelo == "M1", ]
  rg <- read.csv(file.path(RESULTADOS, "aplicacion", "cobertura_censo_region.csv")); rg <- rg[rg$modelo == "M1", ]
  cd <- read.csv(file.path(RESULTADOS, "aplicacion", "cobertura_censo_distancia.csv"))
  cd <- cd[cd$modelo == "M1" & cd$tramo_delta %in% c("[50,100)", "[100,Inf)"), ]
  fila <- function(estrato, regla, cob, lo, hi, ancho, n) data.frame(estrato, regla, cob, lo, hi, ancho, n)
  x <- list()
  for (r in names(OI)) {
    g <- gl[gl$regla == r, ]; x[[length(x) + 1]] <- fila("All districts", r, g$cobertura, g$ic_lo, g$ic_hi, g$ancho_medio, 1783)
    for (reg in c("sierra", "resto")) {
      g <- rg[rg$regla == r & rg$region == reg, ]
      x[[length(x) + 1]] <- fila(if (reg == "sierra") "Highlands" else "Rest of Peru", r, g$cobertura, g$ic_lo, g$ic_hi, g$ancho_medio, g$n)
    }
    g <- cd[cd$regla == r, ]; n <- sum(g$n); k <- round(sum(g$cobertura * g$n)); ci <- binom.test(k, n)$conf.int
    x[[length(x) + 1]] <- fila("50 km or more from\ntraining data", r, k / n, ci[1], ci[2], sum(g$ancho_medio * g$n) / n, n)
  }
  x <- do.call(rbind, x)
  orden <- c("All districts", "Highlands", "Rest of Peru", "50 km or more from\ntraining data")
  x$etiqueta <- factor(sprintf("%s\n(n = %d)", x$estrato, x$n), levels = rev(unique(sprintf("%s\n(n = %d)", x$estrato, x$n)[match(orden, x$estrato)])))
  x$regla <- factor(x$regla, levels = names(OI))
  dy <- c(nominal = 0.24, bloque = 0, mondrian = -0.24)
  x$y <- as.numeric(x$etiqueta) + dy[as.character(x$regla)]
  zebra <- data.frame(ymin = seq(1, 4, 2) - 0.5, ymax = seq(1, 4, 2) + 0.5)
  base_t <- theme_minimal(base_family = "Latin Modern Roman", base_size = 9) +
    theme(panel.grid = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom",
          legend.title = element_blank(), legend.text = element_text(size = 8.4, colour = TINTA2),
          axis.text = element_text(colour = TINTA2), axis.title = element_text(size = 8.64, colour = TINTA),
          plot.title = element_text(size = 9.6, face = "bold", colour = TINTA), plot.title.position = "plot",
          plot.margin = margin(4, 6, 4, 4))
  esc_y <- scale_y_continuous(breaks = 1:4, labels = levels(x$etiqueta), limits = c(0.45, 4.55), expand = c(0, 0))
  capas <- function(g) g + geom_rect(data = zebra, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
                                     inherit.aes = FALSE, fill = "#f5f5f3")
  a <- capas(ggplot(x, aes(cob, y, colour = regla, fill = regla, shape = regla))) +
    geom_vline(xintercept = 0.90, colour = "#8c8c8c", linetype = "22", linewidth = 0.4) +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.7, orientation = "y") +
    geom_point(size = 2.4, colour = "white", stroke = 0.5) +
    scale_colour_manual(values = OI, labels = NOM) + scale_fill_manual(values = OI, labels = NOM) +
    scale_shape_manual(values = FORMA, labels = NOM) + esc_y +
    scale_x_continuous(limits = c(0.68, 1), breaks = seq(0.70, 1, 0.05), labels = function(v) paste0(round(100 * v), "%")) +
    labs(x = "Districts whose census count fell inside\nthe 90% interval (95% CI)", y = NULL, title = "Coverage", tag = "a)") + base_t +
    theme(axis.text.y = element_text(size = 8.64, colour = TINTA, hjust = 1, lineheight = 0.95),
          panel.grid.major.x = element_line(colour = "#ececea", linewidth = 0.3))
  v <- capas(ggplot(x, aes(1, y, colour = regla))) +
    geom_text(aes(label = sprintf("%.1f%%  [%.1f, %.1f]", 100 * cob, 100 * lo, 100 * hi)), hjust = 1, size = 2.82,
              family = "Latin Modern Roman", show.legend = FALSE) +
    scale_colour_manual(values = c(nominal = "#B07800", bloque = "#2f8fc4", mondrian = "#0072B2")) + esc_y +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = "Coverage [95% CI]") + base_t +
    theme(axis.text = element_blank(), plot.title = element_text(size = 8.64, face = "bold", hjust = 1))
  b <- capas(ggplot(x, aes(ancho, y, colour = regla, fill = regla, shape = regla))) +
    geom_segment(aes(x = 0, xend = ancho, yend = y), linewidth = 0.45, alpha = 0.45, show.legend = FALSE) +
    geom_point(size = 2.4, colour = "white", stroke = 0.5, show.legend = FALSE) +
    scale_colour_manual(values = OI, labels = NOM) + scale_fill_manual(values = OI, labels = NOM) +
    scale_shape_manual(values = FORMA, labels = NOM) + esc_y +
    scale_x_continuous(limits = c(0, 0.62), breaks = seq(0, 0.6, 0.2), expand = c(0, 0)) +
    labs(x = "Mean interval width\n(share of crop records)", y = NULL, title = "Width", tag = "b)") + base_t +
    theme(axis.text.y = element_blank(), panel.grid.major.x = element_line(colour = "#ececea", linewidth = 0.3))
  g <- (a + v + b + plot_layout(widths = c(1.9, 0.75, 1), guides = "collect")) &
    theme(legend.position = "bottom", legend.margin = margin(6, 0, 0, 0),
          plot.tag = element_text(size = 10.8, face = "bold", family = "Latin Modern Roman", colour = "black"))
  guardar(g, "Fig3_cobertura_censo_region", 190, 100)
}


# Figura 4 (ggplot2): cobertura frente a ancho de siete reglas en el censo; cada regla va de todos los
# distritos (punto lleno) a los alejados 50 km o más (punto vacío).
fig4_gg <- function() {
  COL7 <- c(nominal = "#E69F00", bloque = "#56B4E9", mondrian = "#0072B2", geoconformal = "#009E73",
            tramos_di = "#CC79A7", fraccion = "#D55E00", bloque_pesos = "#000000")
  NOM7 <- c(nominal = "Bayesian interval (uncalibrated)", bloque = "Block-conformal calibration",
            mondrian = "Region-wise block-conformal", geoconformal = "GeoConformal", tramos_di = "Conformal by DI bins",
            fraccion = "Conformal on the fraction", bloque_pesos = "Block-weighted (Proposition 1)")
  gl <- read.csv(file.path(RESULTADOS, "aplicacion", "cobertura_censo.csv")); gl <- gl[gl$modelo == "M1", ]
  cd <- read.csv(file.path(RESULTADOS, "aplicacion", "cobertura_censo_distancia.csv"))
  cd <- cd[cd$modelo == "M1" & cd$tramo_delta %in% c("[50,100)", "[100,Inf)"), ]
  x <- do.call(rbind, lapply(names(COL7), function(r) {
    g <- gl[gl$regla == r, ]; h <- cd[cd$regla == r, ]; n <- sum(h$n)
    data.frame(regla = r, c0 = g$cobertura, w0 = g$ancho_medio, c1 = sum(h$cobertura * h$n) / n, w1 = sum(h$ancho_medio * h$n) / n)
  }))
  x$regla <- factor(x$regla, levels = names(COL7))
  x$etiqueta <- sprintf("%s\n%.1f%% all,  %.1f%% far", NOM7[as.character(x$regla)], 100 * x$c0, 100 * x$c1)
  # desplazamiento de cada etiqueta (ancho, cobertura), fijado a mano para no tapar los puntos vacíos
  nud <- data.frame(regla = names(COL7), nx = c(0.045, 0.040, -0.030, 0.000, 0.035, -0.010, -0.020),
                    ny = c(-0.012, 0.024, 0.040, 0.034, -0.020, 0.025, -0.016))
  x <- merge(x, nud, by = "regla")
  g <- ggplot(x) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.88, ymax = 0.92, fill = "#eaf3ec") +
    annotate("text", x = 0.236, y = 0.9185, label = "88-92%", hjust = 0, vjust = 1, size = 2.64, colour = "#3f7a52", family = "Latin Modern Roman") +
    geom_hline(yintercept = 0.90, colour = "#8c8c8c", linetype = "22", linewidth = 0.4) +
    geom_segment(aes(x = w0, y = c0, xend = w1, yend = c1, colour = regla), linewidth = 0.7, alpha = 0.6) +
    geom_point(aes(w1, c1, colour = regla), shape = 21, fill = "white", size = 2.4, stroke = 0.9) +
    geom_point(aes(w0, c0, fill = regla), shape = 21, colour = "white", size = 3, stroke = 0.5) +
    ggrepel::geom_label_repel(aes(w0, c0, label = etiqueta, colour = regla), nudge_x = x$nx, nudge_y = x$ny, size = 2.46, family = "Latin Modern Roman",
                              lineheight = 0.95, label.size = 0, fill = alpha("white", 0.85), label.padding = unit(0.6, "mm"),
                              box.padding = 0.3, point.padding = 0.3, force = 0.4, min.segment.length = 0, segment.size = 0.25,
                              segment.colour = "#a8a8a8", seed = SEMILLA, max.overlaps = Inf) +
    scale_colour_manual(values = COL7, guide = "none") + scale_fill_manual(values = COL7, guide = "none") +
    scale_x_continuous(limits = c(0.23, 0.52), breaks = seq(0.25, 0.50, 0.05)) +
    scale_y_continuous(limits = c(0.825, 0.975), breaks = seq(0.84, 0.96, 0.02), labels = function(v) paste0(round(100 * v), "%")) +
    labs(x = "Mean interval width (share of crop records); narrower is more precise",
         y = "Districts whose census count fell\ninside the 90% interval",
         caption = "Filled circle: all districts (n = 1783).   Open circle: districts 50 km or more from training data (n = 388).") +
    theme_minimal(base_family = "Latin Modern Roman", base_size = 9) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#ececea", linewidth = 0.3),
          axis.text = element_text(colour = TINTA2), axis.title = element_text(size = 8.64, colour = TINTA),
          axis.line = element_line(colour = TINTA2, linewidth = 0.3),
          plot.caption = element_text(size = 7.68, colour = TINTA2, hjust = 0), plot.caption.position = "plot",
          plot.margin = margin(4, 8, 4, 4))
  guardar(g, "Fig4_comparacion_conformal", 160, 108)
}


# Figura 5 (ggplot2 + patchwork): curvas riesgo-cobertura con las reglas comparadas (a) y error de los
# 20 bloques de validación (b).
fig5_gg <- function() {
  ab <- read.csv(file.path(RESULTADOS, "aplicacion", "abstencion.csv"))
  ab$clave <- ifelse(ab$regla == "sin abstención", "z", substr(ab$regla, 1, 1))
  nombres <- c(a = "Bayesian interval", b = "Area-of-applicability mask", c = "GeoConformal",
               d = "Conformal by DI bins", e = "Planned rule (block conformal)", g = "Risk-chosen level (exploratory)",
               z = "No abstention")
  ab <- ab[ab$clave %in% names(nombres), ]
  ab$etiqueta <- sprintf("%s\n%.0f%% classified, %.1f%% wrong", nombres[ab$clave], 100 * ab$fraccion_retenida, 100 * ab$riesgo)
  nud <- data.frame(clave = c("a", "b", "c", "d", "e", "g", "z"),
                    nx = c(0.06, -0.20, 0.03, -0.12, -0.10, -0.14, -0.06), ny = c(-0.030, 0.028, 0.036, -0.022, 0.030, -0.036, 0.024))
  ab <- merge(ab, nud, by = "clave")
  ca <- read.csv(file.path(RESULTADOS, "aplicacion", "riesgo_cobertura_alfa.csv")); ca <- ca[order(ca$fraccion_retenida), ]
  cdi <- read.csv(file.path(RESULTADOS, "aplicacion", "riesgo_cobertura.csv")); cdi <- cdi[cdi$metodo == "solo DI", ]
  curvas <- rbind(data.frame(x = ca$fraccion_retenida, y = ca$riesgo, curva = "Abstain where the calibrated interval contains the threshold"),
                  data.frame(x = cdi$fraccion_retenida, y = cdi$riesgo, curva = "Abstain where conditions are new to the model"))
  curvas$curva <- factor(curvas$curva, levels = c("Abstain where the calibrated interval contains the threshold",
                                                  "Abstain where conditions are new to the model"))
  tema_f5 <- theme_minimal(base_family = "Latin Modern Roman", base_size = 9) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#ececea", linewidth = 0.3),
          axis.text = element_text(colour = TINTA2), axis.title = element_text(size = 8.64, colour = TINTA),
          axis.line = element_line(colour = TINTA2, linewidth = 0.3), legend.position = "bottom",
          legend.title = element_blank(), legend.text = element_text(size = 8.16, colour = TINTA2),
          plot.title = element_text(size = 9.6, face = "bold", colour = TINTA), plot.title.position = "plot",
          plot.margin = margin(4, 8, 4, 4))
  pct <- function(v) paste0(round(100 * v), "%")
  a <- ggplot() +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 0.10, fill = "#eaf3ec") +
    annotate("text", x = 1.0, y = 0.003, label = "Target: 10% or less", hjust = 1, vjust = 0, size = 2.64,
             colour = "#3f7a52", family = "Latin Modern Roman") +
    geom_hline(yintercept = 0.10, colour = "#8c8c8c", linetype = "22", linewidth = 0.4) +
    geom_line(data = curvas, aes(x, y, linetype = curva, colour = curva), linewidth = 0.8) +
    geom_linerange(data = ab, aes(x = fraccion_retenida, ymin = riesgo_ic_lo, ymax = riesgo_ic_hi), colour = TINTA, linewidth = 0.5) +
    geom_point(data = ab, aes(fraccion_retenida, riesgo), shape = 21, fill = "white", colour = TINTA, size = 2.4, stroke = 0.8) +
    ggrepel::geom_label_repel(data = ab, aes(fraccion_retenida, riesgo, label = etiqueta), nudge_x = ab$nx, nudge_y = ab$ny,
                              size = 2.4, family = "Latin Modern Roman", lineheight = 0.95, label.size = NA, fill = alpha("white", 0.92),
                              label.padding = unit(0.5, "mm"), min.segment.length = 0, segment.size = 0.25, segment.colour = "#a8a8a8",
                              box.padding = 0.25, force = 0.3, seed = SEMILLA, max.overlaps = Inf) +
    scale_colour_manual(values = c("#56B4E9", "#5b5a57")) + scale_linetype_manual(values = c("solid", "13")) +
    scale_x_continuous(limits = c(0.18, 1.02), breaks = seq(0.2, 1, 0.2), labels = pct) +
    scale_y_continuous(limits = c(0, 0.165), breaks = seq(0, 0.14, 0.02), labels = pct, expand = c(0, 0)) +
    labs(x = "Districts that receive a classification", y = "Misclassified among classified\ndistricts (95% CI)",
         title = "Risk against coverage", tag = "a)") + tema_f5 +
    theme(legend.position = "inside", legend.position.inside = c(0.01, 0.99), legend.justification = c(0, 1),
          legend.key.height = unit(3, "mm"), legend.spacing.y = unit(0, "mm"), legend.background = element_rect(fill = alpha("white", 0.9), colour = NA))
  rb <- read.csv(file.path(RESULTADOS, "aplicacion", "riesgo_por_bloque.csv"))
  orden <- c("sin abstención" = "No\nabstention", "f. propuesta: conformal por bloques + DI con control de riesgo" = "Planned\nrule",
             "g. exploratoria: nivel del conjunto conformal calibrado al riesgo" = "Risk-chosen\nlevel")
  rb <- rb[rb$regla %in% names(orden), ]; rb$grupo <- factor(orden[rb$regla], levels = orden)
  cuenta <- aggregate(error ~ grupo, rb, function(e) sum(e > 0.10))
  b <- ggplot(rb, aes(grupo, error, colour = grupo)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.10, fill = "#eaf3ec") +
    geom_hline(yintercept = 0.10, colour = "#8c8c8c", linetype = "22", linewidth = 0.4) +
    geom_boxplot(width = 0.45, outlier.shape = NA, fill = NA, colour = "#bdbdbd", linewidth = 0.4) +
    ggbeeswarm::geom_beeswarm(size = 1.6, cex = 2.6, alpha = 0.9) +
    geom_text(data = cuenta, aes(grupo, 0.385, label = sprintf("%d of 20\nabove 10%%", error)), colour = TINTA,
              size = 2.52, family = "Latin Modern Roman", lineheight = 0.95, vjust = 1) +
    scale_colour_manual(values = c("#5b5a57", "#56B4E9", "#0072B2"), guide = "none") +
    scale_y_continuous(limits = c(-0.005, 0.39), breaks = seq(0, 0.3, 0.05), labels = pct) +
    labs(x = NULL, y = "Misclassified in each validation block", title = "Error by block", tag = "b)") + tema_f5
  g <- (a + b + plot_layout(widths = c(1.75, 1))) &
    theme(plot.tag = element_text(size = 10.8, face = "bold", family = "Latin Modern Roman", colour = "black"))
  guardar(g, "Fig5_riesgo_cobertura", 190, 92)
}

# Globo de ubicación como imagen aparte, para el recuadro del PNG de la figura 1 (05_figuras.py).
globo_png <- function() {
  ggsave(file.path(FIGURAS, "globo_ubicacion.png"), globo(), width = 40, height = 40, units = "mm", dpi = 600, bg = "transparent")
  cat("guardado: globo_ubicacion.png\n")
}

pedidas <- commandArgs(TRUE)
# Las figuras 1 a 6 se hacen ahora con Python (03-codigo/05_figuras.py); aquí queda solo el globo.
if (!length(pedidas)) pedidas <- c("globo_png", "fig3_gg", "fig4_gg", "fig5_gg")
for (f in pedidas) get(f)()
