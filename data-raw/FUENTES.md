# Fuentes de datos abiertos

Revisadas y descargadas el 7 de octubre de 2026. `descargar.sh todo` baja todas las fuentes; ninguna pide cuenta. El NDVI se extrae por distrito con `../modelo/ndvi_modis.py`.

| Uso | Fuente | Dónde se descarga | Periodo | Login | Licencia |
| --- | --- | --- | --- | --- | --- |
| **Variable respuesta** | IV CENAGRO 2012, INEI, registro de cultivos REC02 en **SPSS** (módulos 337 a 361 = Amazonas a Ucayali) | https://proyectos.inei.gob.pe/iinei/srienaho/descarga/SPSS/<código>-Modulo231.zip (`descargar_cenagro_spss.py`) | Campaña 2011–2012 | No | Datos abiertos del Estado |
| No usar | IV CENAGRO 2012 en CSV (descarga masiva, 931 MB) | https://www.inei.gob.pe/media/DATOS_ABIERTOS/CENAGRO/DATA/2012.zip | — | No | **Los números vienen corruptos** (valores como 6.01E-154) |
| Tabulados del censo | CENAGRO, tabulados y REDATAM | https://censos.inei.gob.pe/cenagro/tabulados/ | 2012 | No | — |
| Validación de H3 | ENAHO módulo 22, INEI (códigos 498 a 966 = 2015 a 2024) | https://proyectos.inei.gob.pe/iinei/srienaho/descarga/STATA/<código>-Modulo22.zip | 2015–2024. La anual 2025 aún no está en el catálogo (consultado el 7-oct-2026) | No | Microdatos INEI |
| Clima principal | TerraClimate (tmax, tmin, ppt, pet, def, soil) | Recorte del Perú por NCSS con `descargar_terraclimate.py` (la agregación OPeNDAP no responde) | 1981–2025 | No | CC0 |
| Clima, sensibilidad | PISCOp v2.1, SENAMHI (vía IRI) | https://iridl.ldeo.columbia.edu/SOURCES/.SENAMHI/.HSR/.PISCO/ | 1981–2016 estable, hasta 2026 inestable | Navegar pide cuenta gratuita; el enlace `data.nc` directo funcionó sin ella | — |
| Clima, sensibilidad | PISCOt v1.1 (Tmáx, Tmín mensuales, estable), SENAMHI vía IRI | https://iridl.ldeo.columbia.edu/SOURCES/.SENAMHI/.HSR/.PISCO/.Temp/.v1p1/ | 1981–2016 | No (enlace `data.nc` directo) | — |
| Clima, opcional | PISCOt v1.2 diaria (hasta 2020) | figshare: 10.6084/m9.figshare.22712365. **Desde esta red figshare corta la conexión**; probar desde otra red | 1981–2020 | No | CC BY 4.0 |
| Vegetación | MODIS MOD13Q1 (Terra) y MYD13Q1 (Aqua) v061, NDVI de 16 días desfasados 8 días, leído a 500 m; la fenología se deriva de la serie | Microsoft Planetary Computer, colección `modis-13Q1-061`. Resultado: `modis/ndvi_16dias_distrito.parquet` (1 961 485 filas, 1 824 distritos, 1 095 fechas) | ago-2000 a jun-2025, que cubre la campaña 2024–2025 | No | Sin restricciones |
| Suelo | SoilGrids 2.0, ISRIC | https://files.isric.org/soilgrids/latest/data/ | Estático | No | CC BY 4.0 |
| Relieve | Copernicus DEM GLO-90 | AWS público `copernicus-dem-90m`, teselas de 1° (`descargar_suelo_relieve.py`) | 2011–2015 | No | Licencia Copernicus |
| Máscara agrícola | Mapa Nacional de Superficie Agrícola 2024, MIDAGRI (`MNSA_2024.rar`, 402 MB; 1 827 distritos; 11 248 646 ha, igual a la RM 0026-2025-MIDAGRI) | https://siea.midagri.gob.pe/herramientas/superficie-agricola (Google Drive) | 2024 | No | Uso con cita |
| Límites | Distritos INEI 2023 | https://ide.inei.gob.pe/files/Distrito.rar | 2023 | No | — |
| Futuro, sensibilidad | WorldClim 2.1 CMIP6 | https://www.worldclim.org/data/cmip6/cmip6climate.html | 2041–2060 | No | Uso académico |

## Qué se descartó y por qué

- **Estaciones del SENAMHI.** El portal pide sesión o una solicitud. Las grillas PISCO ya incluyen la red de estaciones, así que no hacen falta.
- **ENA (INEI).** Solo es representativa por departamento, excluye las comunidades campesinas y cambió su marco muestral en 2023. No sirve como etiqueta distrital.
- **Serie distrital del SIEA (MIDAGRI).** No hay descarga abierta nacional por distrito. Se puede pedir por la Ley 27806 de acceso a la información pública. Si llega, entra como validación adicional.
- **WorldClim 1970–2000 y ERA5-Land como línea base.** Se reemplazan por TerraClimate, que es una sola fuente para los dos periodos. ERA5-Land además pide cuenta.
- **CHIRTS.** Termina en 2016.

## CENAGRO: código de la papa (verificado)

La papa son los códigos **2610 a 2615** de `P024_03`. Su suma nacional, 367 691,9 ha, coincide con las 367 692 ha de los Resultados Definitivos del INEI (2013, cuadro de cultivos), y es la única combinación de códigos 26xx que coincide. Como control, el código 1203 da 425 416 ha, igual que el café, y el 2107 da 261 577 ha, igual que el maíz amarillo duro. El catálogo de cultivos no viene ni en el diccionario ni en la cédula.

La etiqueta distrital está en `etiqueta_papa_distrito_2012.csv` (la genera `../modelo/etiqueta_cenagro.py`): 1 783 distritos con cultivos, 4 814 476 registros de cultivo y 878 028 de papa; con τ = 1 %, 1 224 son positivos.

Cruce con los límites del INEI 2023 (`limites/DISTRITO.gpkg`, 1 890 distritos, campo `ubigeo`): 1 781 de los 1 784 ubigeos del censo coinciden directo. De los tres que no coinciden, 160109 y 160114 pasaron en 2014 a la provincia de Putumayo y se recodifican a 160801 y 160803; el 120699 (Satipo, sin distrito asignado) se excluye, con 34,7 ha de papa. Los 109 distritos creados después de 2012 heredan la etiqueta del distrito del que se separaron (pendiente: tabla de creaciones de distritos del INEI).

## MODIS: controles de calidad

- Se usan solo píxeles con confiabilidad 0 o 1, NDVI crudo entre −2000 y 10000 y fracción agrícola mayor que 1 %.
- Cada distrito y fecha necesita un peso agrícola de al menos medio píxel de 500 m; el NDVI final está entre −0,04 y 0,95.
- Quedan fuera 66 distritos casi sin tierra agrícola (31 de Lima, 7 del Callao, 8 de Arequipa). Tienen 53,7 ha de papa de 367 657 (0,015 %).
- El catálogo trae compuestos hasta junio de 2025, lo que cubre toda la campaña 2024–2025 (septiembre a mayo).
