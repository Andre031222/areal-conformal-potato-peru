# Data sources and how to cite them

Every variable in `tabla_distrital_2012.csv` comes from one of the sources below. If you use the table, cite this repository **and** the original sources. DOIs were checked at doi.org on 16 November 2026; web pages were accessed on 7 October 2026 unless stated otherwise.

| Variable(s) in the table | Source | Licence / terms |
| --- | --- | --- |
| `registros_cultivo`, `registros_papa`, `ha_*`, `fraccion_*`, `y_tau*`, `fraccion_riego` | IV National Agricultural Census 2012 (INEI), crop record REC02 | INEI microdata terms of use; only district aggregates are shared |
| District polygons (`limites/distritos_2012.gpkg`), `ubigeo` | INEI district boundaries 2023, harmonised to the 2012 census | INEI open data |
| Agricultural land mask (aggregation base) | MIDAGRI National Agricultural Land Map 2024 | Peruvian government open data |
| `*_ajuste`, `*_p1`, `*_p2` (climate) | TerraClimate 1981–2025 | CC0 1.0 |
| `spei3_*` | Computed from TerraClimate with the SPEI method | Method citation |
| Climate sensitivity analysis | PISCOp v2.1 and PISCOt v1.1 (SENAMHI) | SENAMHI open data |
| `ph`, `carbono_org_gkg`, `arcilla_pct`, `arena_pct` | SoilGrids 2.0 (ISRIC) | CC BY 4.0 |
| `altitud_*`, `pendiente_grados` | Copernicus DEM GLO-90 | Copernicus DEM licence (free use with attribution) |
| `feno_*` | MODIS MOD13Q1 and MYD13Q1 v061 (NASA LP DAAC), via Microsoft Planetary Computer | NASA open data |
| `viaje_min` | Travel time to cities 2015 (Malaria Atlas Project) | CC BY 4.0 |
| Expert zoning parameters (model M5) | FAO ECOCROP, via the R package Recocrop | GPL-3 (package) |

## References

**Agricultural census.** Instituto Nacional de Estadística e Informática (INEI), 2013. Resultados definitivos: IV Censo Nacional Agropecuario 2012. INEI, Lima. Microdata: https://proyectos.inei.gob.pe/microdatos/ (IV Censo Nacional Agropecuario 2012, module 231, crop record REC02, SPSS format).

**District boundaries.** Instituto Nacional de Estadística e Informática (INEI), 2023. Límites distritales del Perú [dataset]. https://ide.inei.gob.pe/

**Agricultural land.** Ministerio de Desarrollo Agrario y Riego (MIDAGRI), 2025. Mapa Nacional de Superficie Agrícola del Perú 2024, approved by Resolución Ministerial 0026-2025-MIDAGRI [dataset]. https://www.gob.pe/institucion/midagri/normas-legales/642691

**TerraClimate.** Abatzoglou, J.T., Dobrowski, S.Z., Parks, S.A., Hegewisch, K.C., 2018. TerraClimate, a high-resolution global dataset of monthly climate and climatic water balance from 1958–2015. Scientific Data 5, 170191. https://doi.org/10.1038/sdata.2017.191

**SPEI.** Vicente-Serrano, S.M., Beguería, S., López-Moreno, J.I., 2010. A Multiscalar Drought Index Sensitive to Global Warming: The Standardized Precipitation Evapotranspiration Index. Journal of Climate 23, 1696–1718. https://doi.org/10.1175/2009JCLI2909.1

**PISCO precipitation.** Aybar, C., Fernández, C., Huerta, A., Lavado, W., Vega, F., Felipe-Obando, O., 2020. Construction of a high-resolution gridded rainfall dataset for Peru from 1981 to the present day. Hydrological Sciences Journal 65, 770–785. https://doi.org/10.1080/02626667.2019.1649411

**PISCO temperature.** Huerta, A., Aybar, C., Imfeld, N., Correa, K., Felipe-Obando, O., Rau, P., Drenkhan, F., Lavado-Casimiro, W., 2023. High-resolution grids of daily air temperature for Peru - the new PISCOt v1.2 dataset. Scientific Data 10, 847. https://doi.org/10.1038/s41597-023-02777-w

**SoilGrids.** Poggio, L., de Sousa, L.M., Batjes, N.H., Heuvelink, G.B.M., Kempen, B., Ribeiro, E., Rossiter, D., 2021. SoilGrids 2.0: producing soil information for the globe with quantified spatial uncertainty. SOIL 7, 217–240. https://doi.org/10.5194/soil-7-217-2021

**Copernicus DEM.** European Space Agency, Airbus, 2022. Copernicus DEM [dataset]. Copernicus DEM. https://doi.org/10.5270/ESA-c5d3d65

**MODIS Terra.** Didan, K., 2021. MODIS/Terra Vegetation Indices 16-Day L3 Global 250m SIN Grid V061 [dataset]. NASA Land Processes Distributed Active Archive Center. https://doi.org/10.5067/MODIS/MOD13Q1.061

**MODIS Aqua.** Didan, K., 2021. MODIS/Aqua Vegetation Indices 16-Day L3 Global 250m SIN Grid V061 [dataset]. NASA Land Processes Distributed Active Archive Center. https://doi.org/10.5067/MODIS/MYD13Q1.061

**Travel time.** Weiss, D.J., Nelson, A., Gibson, H.S., Temperley, W., Peedell, S., Lieber, A., Hancher, M., Poyart, E., Belchior, S., Fullman, N., Mappin, B., Dalrymple, U., Rozier, J., Lucas, T.C.D., Howes, R.E., Tusting, L.S., Kang, S.Y., Cameron, E., Bisanzio, D., Battle, K.E., Bhatt, S., Gething, P.W., 2018. A global map of travel time to cities to assess inequalities in accessibility in 2015. Nature 553, 333–336. https://doi.org/10.1038/nature25181

**ECOCROP parameters.** Hijmans, R.J., 2021. Recocrop: Estimating Environmental Suitability for Plants [software], version 0.4-2. CRAN: Contributed Packages. https://doi.org/10.32614/CRAN.package.Recocrop

## This repository

Vilca-Solorzano, R.A., Ibañez-Quispe, V., Torres-Cruz, F., 2026. Data and code for: Agricultural suitability of potato across Peru from a calibrated Bayesian areal model. Mendeley Data, V1. https://doi.org/10.17632/mcp876jbyt.1
