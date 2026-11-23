# Tabla distrital 2012 (etapa 2)

Archivo: `tabla_distrital_2012.csv` (y `.parquet`, ignorado por git). Una fila por distrito del IV CENAGRO 2012: **1 783 distritos, 63 columnas**. Generada el 5 de noviembre de 2026.

## Cómo se reproduce

```bash
data-raw/construir_tabla.sh
```

Corre en orden, en unos 2 minutos: `data-raw/etiqueta_cenagro.py` → `armonizar_ubigeo.py` → `puntos_agricolas.py` → `clima_terraclimate.py` → `tierra_controles.py` → `fenologia_ndvi.py` → `tabla_distrital.py`. El último script valida la tabla antes de guardarla y se detiene si algo falla.

## Decisiones de construcción

| Tema | Decisión |
| --- | --- |
| Unidad | Distrito 2012. Los 107 distritos creados después se disuelven en su distrito de origen (`ubigeo_2023_a_2012.csv`). El README decía 109: el GPKG ya trae los dos de Putumayo con su código nuevo |
| Asignación de distritos nuevos | Vecino censal de la misma provincia con el que comparte más límite. Dos quedaron fuera de su provincia y hay que verificarlos contra su ley de creación: **La Morada (100704) → José Crespo y Castillo (100604)** y **Capazo (210502) → Pisacoma (210405)** |
| Base de agregación | 5 936 756 píxeles de la máscara agrícola MIDAGRI 2024 (0,00125°). 17 distritos sin superficie agrícola (Lima Metropolitana, Callao, cercados urbanos) usan todo su territorio: `base_agregacion = "distrito"` |
| Clima | Métricas calculadas en cada celda de TerraClimate y promediadas por píxeles agrícolas. Campaña c = septiembre de c-1 a mayo de c |
| SPEI-3 | Logística generalizada por L-momentos (como el paquete SPEI de R), referencia 1981–2025. Un primer intento con la fórmula log-logística de momentos daba forma negativa en el 18–53 % de las celdas y se descartó |
| Fenología | Año propio de cada distrito desde su mes de NDVI mínimo. En la sierra empieza en septiembre, el pico cae en marzo y la estación verde dura unos 8 meses |

## Columnas

| Grupo | Columnas | Unidad / nota |
| --- | --- | --- |
| Identificación | `ubigeo`, `departamento`, `provincia` | |
| Respuesta (CENAGRO) | `ha_sembrada`, `ha_papa`, `ha_papa_secano`, `registros_cultivo` (n_d), `registros_papa` (k_d), `fraccion_papa`, `fraccion_registros`, `y_tau0.5`, `y_tau1`, `y_tau5` | ha; fracciones 0–1; y = 1 si la fracción de superficie ≥ τ % |
| Suelo (0–30 cm) | `ph`, `carbono_org_gkg`, `arcilla_pct`, `arena_pct` | SoilGrids 2.0 |
| Relieve | `altitud_media_m`, `altitud_p10_m`, `altitud_p90_m`, `pendiente_grados` | Copernicus GLO-90 |
| Controles | `fraccion_riego`, `viaje_min` | CENAGRO; mediana de minutos a la ciudad (MAP 2015) |
| Máscara | `ha_mascara_agricola`, `base_agregacion` | |
| Clima × 3 periodos | `ppt_mm`, `tmax_c`, `tmin_c`, `meses_helada`, `gdd5`, `et0_mm`, `deficit_mm`, `humedad_suelo_mm`, `spei3_medio`, `spei3_min` con sufijo `_ajuste` (campañas 1991–2012), `_p1` (1982–2000) y `_p2` (2001–2025) | Sumas de la campaña (ppt, ET₀, déficit, GDD sobre 5 °C); medias (temperaturas, humedad); meses con Tmín media ≤ 0 °C |
| Fenología 2001–2012 | `feno_ndvi_max`, `feno_ndvi_min`, `feno_ndvi_amplitud`, `feno_duracion_verde_meses`, `feno_mes_pico`, `feno_inicio_verde`, `feno_mes_inicio_anio` | Meses calendario (1 = enero) |
| Región | `sierra` | 1 si la altitud media agrícola ≥ 2 300 m |

Paneles completos para P3 (no van en git): `clima_campana_distrito.parquet` (1 783 × 44 campañas) y `fenologia_anual_distrito.parquet` (1 764 × 24 años).

## Faltantes

- Suelo: 2 distritos urbanos de Lima (150110 Comas, 150117 Los Olivos).
- Fenología: 19 distritos, los 17 urbanos más Caylloma (040505) y Tisco (040517), en la provincia de Caylloma (Arequipa), donde ningún compuesto MODIS pasó el filtro de peso agrícola ≥ 0,5. Antes de H2 conviene revisar por qué, porque son distritos de altura.
- Clima, relieve y respuesta: completos.

## Validaciones que pasan

- La papa suma 367 692 ha (oficial), con los 34,7 ha excluidos de Satipo.
- Máscara: 11,32 millones de ha frente a 11,25 oficiales (0,6 % de diferencia por la aproximación del área del píxel).
- corr(altitud, Tmín de la campaña) = −0,96.
- SPEI-3: sin celdas faltantes; media 0,00.
- Tmín de la campaña: +0,15 °C de media entre 1982–2000 y 2001–2025 (rango −0,58 a +0,81). Meses con helada en la sierra: −0,07.

## Para conversar con el director: H2 no se puede probar como está escrita

| | Distritos | y = 0 con τ = 1 % | y = 0 con τ = 5 % | Mediana de la fracción de registros con papa |
| --- | --- | --- | --- | --- |
| Sierra | 1 138 | **37** | 187 | 0,21 |
| Resto | 645 | 522 | 592 | 0,00 |

En la sierra casi todos los distritos tienen papa. Con τ = 1 % solo hay 37 negativos, y repartidos en k ≥ 20 bloques espaciales la mayoría de los bloques no tendrá ninguno, así que el AUC por bloque de H2 no se puede calcular. En la sierra la pregunta no es **si** hay papa sino **cuánta**.

Propuesta (no cambia el objetivo, solo la métrica):
- Evaluar H2 con la respuesta binomial que ya es la del modelo principal: **log score (devianza predictiva) por bloque, o RMSE de la fracción**, con la misma Wilcoxon pareada.
- Mantener el AUC solo para H1 a escala nacional, donde las clases están equilibradas (69 % con papa).

Esto conviene decidirlo **antes** de correr los modelos y dejarlo escrito con fecha, para que no parezca que se cambió la métrica después de ver los resultados.

## Clima PISCO (sensibilidad)

Archivo: `clima_pisco_distrito_2012.csv`. Las cuatro variables climáticas de M1 (`tmin_c_pisco`, `ppt_mm_pisco`, `meses_helada_pisco`, `spei3_min_pisco`) calculadas con PISCOp v2.1 y PISCOt v1.1 (SENAMHI), campañas 1991–2012, con el mismo método que TerraClimate (métricas por celda y promedio por píxeles agrícolas). PISCO no trae evapotranspiración: ET0 por Hargreaves-Samani (FAO-56). Se reproduce con `bash data-raw/descargar.sh pisco` y `python3 data-raw/clima_pisco.py`, después de `puntos_agricolas.py`.
