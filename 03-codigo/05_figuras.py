"""P1 · Paso 5 (versión Python): figuras del manuscrito con seaborn y matplotlib.

Lee las tablas de 04-resultados y guarda cada figura en figuras/ como PDF (vectorial), TIFF y PNG
a 600 dpi, a 190 mm de ancho (página completa de Elsevier). Colores de Okabe e Ito.
Uso: .venv/bin/python 03-codigo/05_figuras.py [fig2 ...]
"""
import pathlib
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

RAIZ = pathlib.Path(__file__).resolve().parents[1]
RES = RAIZ / "04-resultados"
FIG = RAIZ / "figuras"
MM = 1 / 25.4

TINTA, TINTA2, GRIS = "#1a1a1a", "#5b5a57", "#d9d8d3"
BANDA = "#e3f1e7"
# Okabe-Ito, igual en todas las figuras: naranja = sin calibrar, celeste = por bloques, azul = por región.
COL = {"nominal": "#E69F00", "bloque": "#56B4E9", "mondrian": "#0072B2",
       "geoconformal": "#009E73", "tramos_di": "#CC79A7", "fraccion": "#D55E00", "bloque_pesos": "#000000"}
NOMBRE = {"nominal": "Bayesian 90% set", "bloque": "Block conformal (PIT score)",
          "mondrian": "Block conformal by region", "geoconformal": "GeoConformal",
          "tramos_di": "Conformal by DI bins", "fraccion": "Conformal on the fraction",
          "bloque_pesos": "Block-weighted (Proposition 1)"}
TRAMOS = ["[0,10)", "[10,25)", "[25,50)", "[50,100)", "[100,Inf)"]
ETIQ_TRAMO = ["0–10", "10–25", "25–50", "50–100", "100 or more"]

sns.set_theme(style="ticks", context="paper")
import glob as _glob, os
from matplotlib import font_manager as _fm
for _f in _glob.glob(str(pathlib.Path(__file__).resolve().parents[1] / "data-raw" / "fuentes_lm" / "lmroman10-*.ttf")):
    _fm.fontManager.addfont(_f)
mpl.rcParams.update({
    "font.family": "Latin Modern Roman", "font.size": 8.4, "axes.titlesize": 9, "axes.labelsize": 8.4,
    "xtick.labelsize": 7.8, "ytick.labelsize": 7.8, "legend.fontsize": 8.4,
    "axes.edgecolor": TINTA2, "axes.labelcolor": TINTA, "xtick.color": TINTA2, "ytick.color": TINTA2,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "pdf.fonttype": 42, "ps.fonttype": 42,
    "savefig.dpi": 600, "figure.dpi": 150,
})


def guardar(fig, nombre):
    for ext in ("pdf", "png", "tiff"):
        kw = {"pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tiff" else {}
        if ext == "pdf": kw["dpi"] = 300
        fig.savefig(FIG / f"{nombre}.{ext}", bbox_inches="tight", pad_inches=0.02, facecolor="white", **kw)
    plt.close(fig)
    print("guardada:", FIG / nombre)


def pct(ax, lo, hi, paso, desde=None):
    ax.set_ylim(lo, hi)
    ax.set_yticks(np.arange(desde if desde is not None else lo, hi + 1e-9, paso))
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{100 * v:.0f}%"))


def banda_objetivo(ax):
    ax.axhspan(0.88, 0.92, color=BANDA, zorder=0, lw=0)
    ax.axhline(0.90, color=TINTA2, lw=0.7, ls=(0, (3, 2)), zorder=1)


def fig2():
    """Simulación: cobertura de cada réplica (a) y ancho medio del intervalo (b), sin calibrar y calibrado,
    con el modelo bien especificado (arriba) y mal especificado (abajo); bloques de 145 km."""
    from matplotlib.gridspec import GridSpec
    rep = pd.concat([pd.read_csv(RES / "simulacion" / "cobertura_por_replica.csv"),
                     pd.read_csv(RES / "simulacion" / "cobertura_por_replica_mal.csv")])
    rep = rep[rep["diseno"] == "bloque145"].copy(); rep["phi"] = rep["phi"].astype(str)
    grupos = [("Correctly specified model", [("0.2", "Spatial share 0.2"), ("0.5", "Spatial share 0.5"), ("0.9", "Spatial share 0.9")]),
              ("Misspecified model (spatial share 0.5)", [("colas", "Heavy-tailed effect"),
                                                          ("no_estacionario", "More variance in highlands"),
                                                          ("no_estacionario_este", "More variance in the east")])]
    mc = 1.96 * np.sqrt(0.9 * 0.1 / 1783)       # error Monte Carlo de una réplica (1783 distritos)
    fig = plt.figure(figsize=(190 * MM, 112 * MM))
    gs = GridSpec(2, 2, width_ratios=[1.65, 1], hspace=0.32, wspace=0.08, figure=fig)
    rng = np.random.default_rng(1)
    off = {"nominal": 0.19, "bloque": -0.19}
    for fila, (titulo, esc) in enumerate(grupos):
        axc = fig.add_subplot(gs[fila, 0]); axw = fig.add_subplot(gs[fila, 1], sharey=axc)
        axc.axvspan(0.90 - mc, 0.90 + mc, color=BANDA, zorder=0, lw=0)
        axc.axvline(0.90, color=TINTA2, lw=0.8, ls=(0, (3, 2)), zorder=1)
        for k, (clave, nombre) in enumerate(esc):
            y = -k
            d = rep[rep["phi"] == clave]
            for r in ("nominal", "bloque"):
                v = d[f"cub_{r}"].to_numpy(); yy = y + off[r]
                axc.boxplot(v, positions=[yy], vert=False, widths=0.26, showfliers=False, patch_artist=True,
                            boxprops=dict(fc=COL[r], alpha=0.18, ec=COL[r], lw=0.8), medianprops=dict(color=COL[r], lw=1.4),
                            whiskerprops=dict(color=COL[r], lw=0.8), capprops=dict(color=COL[r], lw=0.8), zorder=2)
                axc.scatter(v, yy + rng.uniform(-0.08, 0.08, len(v)), s=7, color=COL[r], edgecolor="white", linewidth=0.3, zorder=3)
                axc.plot(v.mean(), yy, "D", ms=4.2, color="white", mec=TINTA, mew=0.8, zorder=4)
                w = d[f"ancho_{r}"].to_numpy()
                axw.boxplot(w, positions=[yy], vert=False, widths=0.26, showfliers=False, patch_artist=True,
                            boxprops=dict(fc=COL[r], alpha=0.18, ec=COL[r], lw=0.8), medianprops=dict(color=COL[r], lw=1.4),
                            whiskerprops=dict(color=COL[r], lw=0.8), capprops=dict(color=COL[r], lw=0.8), zorder=2)
                axw.scatter(w, yy + rng.uniform(-0.08, 0.08, len(w)), s=7, color=COL[r], edgecolor="white", linewidth=0.3, zorder=3)
                axw.plot(w.mean(), yy, "D", ms=4.2, color="white", mec=TINTA, mew=0.8, zorder=4)
        axc.set_yticks([-k for k in range(len(esc))], [n for _, n in esc])
        axc.set_ylim(-len(esc) + 0.45, 0.55)
        axc.set_xlim(0.855, 0.945)
        axc.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{100 * v:.0f}%"))
        axc.tick_params(axis="y", length=0)
        axc.grid(axis="x", color="#efeeea", lw=0.5, zorder=0)
        axc.set_title(titulo, loc="left", fontsize=8.88, fontweight="bold", pad=4)
        axw.set_xlim(0.28, 0.66)
        axw.tick_params(axis="y", length=0, labelleft=False)
        axw.grid(axis="x", color="#efeeea", lw=0.5, zorder=0)
        sns.despine(ax=axc, left=True); sns.despine(ax=axw, left=True)
        if fila == 0:
            axc.annotate("a)", xy=(0, 1), xycoords="axes fraction", xytext=(-8, 14), textcoords="offset points", ha="right",
                         fontsize=10.8, fontweight="bold", va="bottom")
            axc.text(0.0, 1.17, "Coverage of each replicate", transform=axc.transAxes, fontsize=8.88, color=TINTA)
            axw.annotate("b)", xy=(0, 1), xycoords="axes fraction", xytext=(-8, 14), textcoords="offset points",
                         fontsize=10.8, fontweight="bold", va="bottom", ha="right")
            axw.text(0.0, 1.17, "Interval width of each replicate", transform=axw.transAxes, fontsize=8.88, color=TINTA)
        else:
            axc.set_xlabel("Held-out districts inside the 90% interval")
            axw.set_xlabel("Width (share of crop records)")
    manejos = [Patch(fc=COL["nominal"], alpha=0.6, label="Bayesian interval (uncalibrated)"),
               Patch(fc=COL["bloque"], alpha=0.6, label="After block-conformal calibration"),
               Line2D([], [], marker="D", lw=0, ms=4, mfc="white", mec=TINTA, label="Mean of the replicates"),
               Patch(fc=BANDA, label="90% ± Monte Carlo error of one replicate")]
    fig.legend(handles=manejos, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.075),
               fontsize=7.68, columnspacing=2.0, handlelength=1.3)
    guardar(fig, "Fig2_cobertura_simulacion")


def fig3():
    """Censo, gráfico de bosque: cobertura del intervalo al 90 % de M1 por estrato, sin calibrar,
    calibrado por bloques y calibrado por región, con su IC de Clopper-Pearson y una columna de valores."""
    from scipy.stats import beta as _b
    gl = pd.read_csv(RES / "aplicacion" / "cobertura_censo.csv"); gl = gl[gl["modelo"] == "M1"].set_index("regla")
    rg = pd.read_csv(RES / "aplicacion" / "cobertura_censo_region.csv"); rg = rg[rg["modelo"] == "M1"]
    cd = pd.read_csv(RES / "aplicacion" / "cobertura_censo_distancia.csv")
    cd = cd[(cd["modelo"] == "M1") & cd["tramo_delta"].isin(["[50,100)", "[100,Inf)"])]

    def lejos(r):
        g = cd[cd["regla"] == r]; n = int(g["n"].sum()); k = int(round((g["cobertura"] * g["n"]).sum()))
        return k / n, _b.ppf(0.025, k, n - k + 1), _b.ppf(0.975, k + 1, n - k), n

    reglas = ["nominal", "bloque", "mondrian"]
    marca = {"nominal": "o", "bloque": "^", "mondrian": "s"}
    etiq = {"nominal": "Bayesian interval (uncalibrated)", "bloque": "Block-conformal calibration",
            "mondrian": "Block-conformal calibration by region"}
    def ancho_lejos(r):
        g = cd[cd["regla"] == r]; return (g["ancho_medio"] * g["n"]).sum() / g["n"].sum()
    anchos = [{r: gl.loc[r, "ancho_medio"] for r in reglas}]
    for reg in ("sierra", "resto"):
        g = rg[rg["region"] == reg].set_index("regla"); anchos.append({r: g.loc[r, "ancho_medio"] for r in reglas})
    anchos.append({r: ancho_lejos(r) for r in reglas})
    estratos = [("All districts", {r: tuple(gl.loc[r, ["cobertura", "ic_lo", "ic_hi"]]) for r in reglas}, 1783)]
    for reg, nombre in (("sierra", "Highlands"), ("resto", "Rest of Peru")):
        g = rg[rg["region"] == reg].set_index("regla")
        estratos.append((nombre, {r: tuple(g.loc[r, ["cobertura", "ic_lo", "ic_hi"]]) for r in reglas}, int(g.loc["nominal", "n"])))
    estratos.append(("50 km or more from\ntraining data", {r: lejos(r)[:3] for r in reglas}, lejos("nominal")[3]))

    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=(190 * MM, 96 * MM))
    gs = GridSpec(1, 2, width_ratios=[1.9, 1], wspace=0.62, figure=fig)
    ax = fig.add_subplot(gs[0]); axw = fig.add_subplot(gs[1], sharey=ax)
    off = {"nominal": 0.26, "bloque": 0.0, "mondrian": -0.26}
    ax.axvline(0.90, color=TINTA, lw=0.9, ls=(0, (3, 2)), zorder=1)
    yt = []
    for k, (nombre, vals, n) in enumerate(estratos):
        y = -k * 1.25; yt.append(y)
        if k % 2 == 0:
            ax.axhspan(y - 0.55, y + 0.55, color="#f6f5f2", zorder=0, lw=0)
            axw.axhspan(y - 0.55, y + 0.55, color="#f6f5f2", zorder=0, lw=0)
        for r in reglas:
            axw.plot(anchos[k][r], y + off[r], marca[r], ms=5.2, color=COL[r], mec="white", mew=0.6, zorder=3)
            axw.plot([0, anchos[k][r]], [y + off[r]] * 2, color=COL[r], lw=0.6, alpha=0.35, zorder=2)
        for r in reglas:
            m, lo, hi = vals[r]; yy = y + off[r]
            ax.errorbar(m, yy, xerr=[[m - lo], [hi - m]], fmt=marca[r], ms=5.2, color=COL[r], mec="white", mew=0.6,
                        elinewidth=1.1, capsize=0, zorder=3)
            ax.text(1.035, yy, f"{100 * m:.1f}%  [{100 * lo:.1f}, {100 * hi:.1f}]", transform=mpl.transforms.blended_transform_factory(ax.transAxes, ax.transData),
                    ha="left", va="center", fontsize=7.56, color="#2a86bd" if r == "mondrian" else COL[r])
    ax.set_yticks(yt, [f"{e[0]}\n(n = {e[2]})" for e in estratos])
    ax.set_ylim(yt[-1] - 0.7, 0.7)
    ax.set_xlim(0.66, 1.0)
    ticks = [0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00]
    ax.set_xticks(ticks, [f"{100 * t:.0f}%" for t in ticks])
    for t in ax.get_xticklabels():
        if t.get_text() == "90%": t.set_fontweight("bold"); t.set_color(TINTA)
    ax.text(1.035, 0.62, "Coverage [95% CI]", transform=mpl.transforms.blended_transform_factory(ax.transAxes, ax.transData),
            ha="left", va="bottom", fontsize=7.92, fontweight="bold", color=TINTA)
    ax.set_xlabel("Census districts whose count fell inside the 90% interval (Clopper–Pearson 95% CI)")
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color="#efeeea", lw=0.5, zorder=0)
    axw.set_xlim(0, 0.62); axw.set_xticks([0, 0.2, 0.4, 0.6])
    axw.set_xlabel("Mean interval width\n(share of crop records)")
    axw.tick_params(axis="y", length=0, labelleft=False)
    axw.grid(axis="x", color="#efeeea", lw=0.5, zorder=0)
    sns.despine(ax=axw, left=True)
    for a_, l_, tt in ((ax, "a", "Coverage"), (axw, "b", "Width")):
        a_.annotate(l_, xy=(0, 1), xycoords="axes fraction", xytext=(-6, 6), textcoords="offset points",
                    fontsize=10.8, fontweight="bold", ha="right", va="bottom")
        a_.set_title(tt, loc="left", fontsize=8.88, pad=6)
    ax.legend(handles=[Line2D([], [], marker=marca[r], lw=0, ms=5.5, color=COL[r], mec="white", label=etiq[r]) for r in reglas],
              loc="upper center", bbox_to_anchor=(0.85, -0.18), ncol=3, frameon=False, handletextpad=0.3, columnspacing=1.4)
    sns.despine(ax=ax, left=True)
    guardar(fig, "Fig3_cobertura_censo_region")


def fig4():
    """Cobertura frente a ancho en el censo: cada regla es una trayectoria de todos los distritos
    (punto lleno) a los alejados 50 km o más (punto vacío)."""
    gl = pd.read_csv(RES / "aplicacion" / "cobertura_censo.csv"); gl = gl[gl["modelo"] == "M1"].set_index("regla")
    cd = pd.read_csv(RES / "aplicacion" / "cobertura_censo_distancia.csv")
    cd = cd[(cd["modelo"] == "M1") & cd["tramo_delta"].isin(["[50,100)", "[100,Inf)"])]
    lejos = {}
    for r, g in cd.groupby("regla"):
        n = g["n"].sum(); lejos[r] = ((g["cobertura"] * g["n"]).sum() / n, (g["ancho_medio"] * g["n"]).sum() / n)
    n_lejos = int(cd[cd["regla"] == "nominal"]["n"].sum())
    reglas = ["nominal", "fraccion", "mondrian", "geoconformal", "bloque", "tramos_di", "bloque_pesos"]
    # posición de cada etiqueta (ancho, cobertura), junto al punto de todos los distritos
    pos = {"nominal": (0.335, 0.840), "fraccion": (0.240, 0.912), "mondrian": (0.300, 0.936),
           "geoconformal": (0.352, 0.927), "bloque": (0.400, 0.876), "tramos_di": (0.405, 0.915),
           "bloque_pesos": (0.425, 0.968)}
    fig, ax = plt.subplots(figsize=(150 * MM, 100 * MM))
    ax.axhspan(0.88, 0.92, color=BANDA, zorder=0, lw=0)
    ax.axhline(0.90, color=TINTA2, lw=0.7, ls=(0, (3, 2)), zorder=1)
    for r in reglas:
        c0, w0 = gl.loc[r, "cobertura"], gl.loc[r, "ancho_medio"]
        c1, w1 = lejos[r]
        ax.plot([w0, w1], [c0, c1], color=COL[r], lw=1.0, alpha=0.75, zorder=2)
        ax.plot(w0, c0, "o", color=COL[r], ms=6.5, mec="white", mew=0.8, zorder=4)
        ax.plot(w1, c1, "o", mfc="white", mec=COL[r], mew=1.4, ms=5.5, zorder=4)
        lx, ly = pos[r]
        ax.annotate(f"{NOMBRE[r]}\n{100 * c0:.1f}% all, {100 * c1:.1f}% far", xy=(w0, c0), xytext=(lx, ly),
                    fontsize=7.44, color=COL[r] if r != "bloque_pesos" else TINTA, ha="left", va="center",
                    linespacing=1.15, arrowprops=dict(arrowstyle="-", color=GRIS, lw=0.5, shrinkA=1, shrinkB=4))
    ax.set_xlim(0.235, 0.515)
    ax.set_ylim(0.825, 0.98)
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{100 * v:.0f}%"))
    ax.set_xlabel("Mean width of the set (share of crop records with potato); left is narrower, more precise")
    ax.set_ylabel("Districts inside the 90% set")
    ax.grid(color="#efeeea", lw=0.5, zorder=0)
    manejos = [Line2D([], [], color=TINTA2, marker="o", ms=6, lw=0, mec="white", label="All districts (n = 1783)"),
               Line2D([], [], color=TINTA2, marker="o", ms=5, lw=0, mfc="white", mew=1.3,
                      label=f"Districts 50 km or more from training data (n = {n_lejos})"),
               Patch(fc=BANDA, ec="none", label="88–92%")]
    ax.legend(handles=manejos, loc="upper left", frameon=True, framealpha=1, edgecolor=GRIS, fontsize=7.56)
    sns.despine(ax=ax)
    guardar(fig, "Fig4_comparacion_conformal")


def fig5():
    """Abstención: curvas riesgo-cobertura con las reglas comparadas (a) y error por bloque (b)."""
    from matplotlib.gridspec import GridSpec
    ab = pd.read_csv(RES / "aplicacion" / "abstencion.csv")
    ab["clave"] = ab["regla"].str[0].where(ab["regla"] != "sin abstención", "z")
    ca = pd.read_csv(RES / "aplicacion" / "riesgo_cobertura_alfa.csv").sort_values("fraccion_retenida")
    cdi = pd.read_csv(RES / "aplicacion" / "riesgo_cobertura.csv")
    cdi = cdi[cdi["metodo"] == "solo DI"].sort_values("fraccion_retenida")
    rb = pd.read_csv(RES / "aplicacion" / "riesgo_por_bloque.csv")

    fig = plt.figure(figsize=(190 * MM, 100 * MM))
    gs = GridSpec(1, 2, width_ratios=[1.55, 1], wspace=0.28, figure=fig)
    ax = fig.add_subplot(gs[0])
    ax.axhspan(0, 0.10, color=BANDA, zorder=0, lw=0)
    ax.axhline(0.10, color=TINTA2, lw=0.7, ls=(0, (3, 2)), zorder=1)
    ax.plot(ca["fraccion_retenida"], ca["riesgo"], color=COL["bloque"], lw=1.6, zorder=2,
            label="Abstain where the block-conformal set contains the threshold\n(set level varied from 90% to 10%)")
    ax.plot(cdi["fraccion_retenida"], cdi["riesgo"], color=TINTA2, lw=1.0, ls=(0, (1, 1.5)), zorder=2,
            label="Abstain where conditions are new to the model\n(dissimilarity threshold varied)")
    texto = {"a": ("Bayesian 90% set", (0.50, 0.016)), "c": ("GeoConformal", (0.47, 0.084)),
             "d": ("Conformal by DI bins", (0.21, 0.020)), "e": ("Planned rule\n(block conformal 90%)", (0.19, 0.074)),
             "b": ("Area-of-applicability mask", (0.56, 0.124)), "g": ("Risk-chosen set level\n(exploratory)", (0.66, 0.040)),
             "z": ("No abstention", (0.84, 0.131))}
    for clave, (nombre, (lx, ly)) in texto.items():
        f = ab[ab["clave"] == clave].iloc[0]
        x, y = f["fraccion_retenida"], f["riesgo"]
        ax.errorbar(x, y, yerr=[[y - f["riesgo_ic_lo"]], [f["riesgo_ic_hi"] - y]], fmt="o", ms=5, mfc="white",
                    mec=TINTA, mew=1.0, color=TINTA, elinewidth=0.8, capsize=0, zorder=4)
        ax.annotate(f"{nombre}\n{100 * x:.0f}% classified, {100 * y:.1f}% wrong", xy=(x, y), xytext=(lx, ly),
                    fontsize=7.08, color=TINTA, linespacing=1.15, va="center",
                    arrowprops=dict(arrowstyle="-", color=GRIS, lw=0.5, shrinkA=1, shrinkB=4))
    ax.set_xlim(0.18, 1.02); ax.set_ylim(0, 0.142)
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{100 * v:.0f}%"))
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{100 * v:.0f}%"))
    ax.set_xlabel("Districts that receive a classification")
    ax.set_ylabel("Misclassified among classified districts (95% CI)")
    ax.text(0.86, 0.004, "target: 10% or less", fontsize=7.2, color="#2f6b45")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, frameon=False, fontsize=7.2, handlelength=2.2)
    ax.set_title("Risk against coverage", loc="left", pad=4)
    ax.annotate("a)", xy=(0, 1), xycoords="axes fraction", xytext=(-34, 8), textcoords="offset points",
                fontsize=10.8, fontweight="bold")
    ax.grid(color="#efeeea", lw=0.5, zorder=0)
    sns.despine(ax=ax)

    # b: error de cada uno de los 20 bloques (enjambre) con su resumen (caja)
    ax = fig.add_subplot(gs[1])
    orden = {"sin abstención": "No\nabstention", "f. propuesta: conformal por bloques + DI con control de riesgo": "Planned\nrule",
             "g. exploratoria: nivel del conjunto conformal calibrado al riesgo": "Risk-chosen\nlevel"}
    d = rb[rb["regla"].isin(orden)].copy(); d["regla"] = d["regla"].map(orden)
    colores = {"No\nabstention": TINTA2, "Planned\nrule": COL["bloque"], "Risk-chosen\nlevel": COL["mondrian"]}
    ax.axhspan(0, 0.10, color=BANDA, zorder=0, lw=0)
    ax.axhline(0.10, color=TINTA2, lw=0.7, ls=(0, (3, 2)), zorder=1)
    sns.boxplot(data=d, x="regla", y="error", order=list(orden.values()), ax=ax, width=0.5, showfliers=False,
                boxprops=dict(facecolor="none", edgecolor=GRIS, lw=0.8), whiskerprops=dict(color=GRIS, lw=0.8),
                capprops=dict(color=GRIS, lw=0.8), medianprops=dict(color=TINTA, lw=1.2), zorder=2)
    sns.swarmplot(data=d, x="regla", y="error", order=list(orden.values()), hue="regla", palette=colores,
                  size=3.6, edgecolor="white", linewidth=0.4, legend=False, ax=ax, zorder=3)
    for i, r in enumerate(orden.values()):
        k = int((d[d["regla"] == r]["error"] > 0.10).sum())
        ax.text(i, 0.335, f"{k} of 20\nabove 10%", ha="center", va="top", fontsize=7.2, color=TINTA, linespacing=1.15)
    ax.set_ylim(-0.01, 0.34)
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{100 * v:.0f}%"))
    ax.set_xlabel(""); ax.set_ylabel("Misclassified in each validation block")
    ax.set_title("Error varies between blocks", loc="left", pad=4)
    ax.annotate("b)", xy=(0, 1), xycoords="axes fraction", xytext=(-34, 8), textcoords="offset points",
                fontsize=10.8, fontweight="bold")
    ax.grid(axis="y", color="#efeeea", lw=0.5, zorder=0)
    sns.despine(ax=ax)
    guardar(fig, "Fig5_riesgo_cobertura")


def _mapa_base(ax, relieve, pais, bloques=None):
    """Relieve sombreado, borde del país, flecha del norte y escala, en EPSG:32718."""
    import rasterio
    with rasterio.open(relieve) as r:
        sh = r.read(1, masked=True); b = r.bounds
    alfa = np.clip((0.80 - sh.filled(np.nan)) * 0.75, 0, 0.42)
    capa = np.zeros(sh.shape + (4,)); capa[..., 3] = np.nan_to_num(alfa)
    ax.imshow(capa, extent=(b.left, b.right, b.bottom, b.top), zorder=3, interpolation="bilinear")
    if bloques is not None:
        bloques.boundary.plot(ax=ax, color=TINTA, lw=0.35, zorder=4)
    pais.boundary.plot(ax=ax, color=TINTA2, lw=0.6, zorder=5)
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    xs, ys = x0 + 0.06 * (x1 - x0), y0 + 0.22 * (y1 - y0)
    ax.plot([xs, xs + 300e3], [ys, ys], color=TINTA, lw=1.2, solid_capstyle="butt", zorder=6)
    for k in (0, 100e3, 200e3, 300e3):
        ax.plot([xs + k, xs + k], [ys, ys + 12e3], color=TINTA, lw=0.7, zorder=6)
    ax.text(xs + 150e3, ys + 22e3, "300 km", ha="center", va="bottom", fontsize=7.2, color=TINTA)
    ax.annotate("N", xy=(x1 - 0.07 * (x1 - x0), y1 - 0.04 * (y1 - y0)), xytext=(x1 - 0.07 * (x1 - x0), y1 - 0.12 * (y1 - y0)),
                ha="center", va="center", fontsize=8.4, fontweight="bold", color=TINTA,
                arrowprops=dict(arrowstyle="-|>,head_length=0.6,head_width=0.3", color=TINTA, lw=1.0))
    ax.set_axis_off()


def fig6():
    """Mapas de la papa predicha por M1 fuera de bloque: fracción predicha (a), ancho del intervalo
    calibrado (b) y decisión del mapa con abstención (c), sobre el mismo relieve de la figura 1."""
    import geopandas as gpd
    import matplotlib.patheffects as pe
    from pyproj import Transformer
    crs = 32718
    d = gpd.read_file(RAIZ / "data" / "limites" / "distritos_2012.gpkg").to_crs(crs)
    dec = pd.read_csv(RES / "aplicacion" / "decision_por_distrito.csv", dtype={"ubigeo": str})
    an = pd.read_csv(RES / "aplicacion" / "ancho_por_distrito.csv", dtype={"ubigeo": str})
    d = d.merge(dec[["ubigeo", "p_media_m1", "yhat", "error", "nc_riesgo"]], on="ubigeo").merge(an, on="ubigeo")
    d["dep"] = d["ubigeo"].str[:2]
    pais = d.dissolve(); deps = d.dissolve(by="dep")
    pais["geometry"] = pais.simplify(400); deps["geometry"] = deps.simplify(400)
    x0, y0, x1, y1 = pais.total_bounds
    tr = Transformer.from_crs(4326, crs, always_xy=True)
    halo = [pe.withStroke(linewidth=1.5, foreground="white")]
    from matplotlib.gridspec import GridSpec
    NOMBRES = {"01": "Amazonas", "02": "Áncash", "03": "Apurímac", "04": "Arequipa", "05": "Ayacucho",
               "06": "Cajamarca", "08": "Cusco", "09": "Huancavelica", "10": "Huánuco", "11": "Ica",
               "12": "Junín", "13": "La Libertad", "14": "Lambayeque", "15": "Lima", "16": "Loreto",
               "17": "Madre de Dios", "18": "Moquegua", "19": "Pasco", "20": "Piura", "21": "Puno",
               "22": "San Martín", "23": "Tacna", "24": "Tumbes", "25": "Ucayali"}
    fig = plt.figure(figsize=(190 * MM, 200 * MM))
    gs = GridSpec(2, 2, width_ratios=[1.75, 1], height_ratios=[1, 1], wspace=0.06, hspace=0.42, figure=fig)
    ax_a = fig.add_subplot(gs[:, 0]); ax_b = fig.add_subplot(gs[0, 1]); ax_c = fig.add_subplot(gs[1, 1])

    def base(ax, capas, titulo, letra, leyenda, titulo_leyenda, grande=False, ncol=1):
        ax.set_xlim(x0 - 25e3, x1 + 25e3); ax.set_ylim(y0 - 25e3, y1 + 25e3); ax.set_aspect("equal")
        _tinte_relieve(ax, alfa_tinte=0.5)
        for valores, color in capas:
            d[valores].plot(ax=ax, color=color, lw=0, alpha=0.9, zorder=3, rasterized=True)
        deps.boundary.plot(ax=ax, color="white", lw=0.7 if grande else 0.35, zorder=4)
        pais.boundary.plot(ax=ax, color=TINTA2, lw=0.9 if grande else 0.7, zorder=5)
        _marco_geografico(ax, crs)
        if grande:
            for cod, g in deps.iterrows():
                if cod not in NOMBRES: continue
                c = g.geometry.representative_point()
                ax.text(c.x, c.y, NOMBRES[cod], ha="center", va="center", fontsize=6.4, fontstyle="italic",
                        color="#2f2e2b", path_effects=[pe.withStroke(linewidth=1.8, foreground="white")], zorder=6)
            lon, lat, _ = CIUDADES["Lima"]; cx, cy = tr.transform(lon, lat)
            ax.plot(cx, cy, "*", ms=7, mfc="#f0e442", mec=TINTA, mew=0.6, zorder=7)
            _escala_norte(ax, fx=0.40, fy=0.02)
            ax.legend(handles=leyenda, title=titulo_leyenda, loc="lower left", bbox_to_anchor=(0.0, 0.0), ncol=ncol,
                      frameon=False, fontsize=7.2, title_fontsize=7.6, alignment="left", handletextpad=0.4)
            ax.set_title(titulo, loc="left", fontsize=9.6, pad=5)
        else:
            ax.set_xticklabels([]); ax.yaxis.tick_right()
            ax.legend(handles=leyenda, title=titulo_leyenda, loc="upper left", bbox_to_anchor=(0.0, -0.02), ncol=ncol,
                      frameon=False, fontsize=6.6, title_fontsize=7.0, alignment="left", columnspacing=0.8, handletextpad=0.4)
            ax.set_title(titulo, loc="left", fontsize=8.4, pad=4)
        ax.annotate(letra, xy=(0, 1), xycoords="axes fraction", xytext=(-6, 4), textcoords="offset points",
                    fontsize=10.8, fontweight="bold", ha="right", va="bottom")

    # a: fracción predicha, mismas clases y colores que las burbujas de la figura 1 (panel grande)
    cortes = [0, 0.05, 0.15, 0.30, 0.45, 1.01]
    pal = sns.color_palette("rocket_r", len(cortes) - 1)
    k = np.digitize(d["p_media_m1"], cortes) - 1
    capas = [(k == i, pal[i]) for i in range(len(pal))]
    ley = [Patch(fc=pal[i], label=f"{100 * a:.0f}–{100 * b:.0f}%" if b < 1 else f"{100 * a:.0f}% or more")
           for i, (a, b) in enumerate(zip(cortes[:-1], cortes[1:]))]
    base(ax_a, capas, "Predicted potato share", "a)", ley, "Predicted share of crop\nrecords with potato (M1)", grande=True)

    # b: ancho del intervalo calibrado (incertidumbre)
    cortes = [0, 0.05, 0.15, 0.30, 0.50, 1.0001]
    pal = sns.color_palette("mako_r", len(cortes) - 1)
    k = np.digitize(d["ancho.bloque"], cortes) - 1
    capas = [(k == i, pal[i]) for i in range(len(pal))]
    ley = [Patch(fc=pal[i], label=f"{a:.2f}–{b:.2f}" if b < 1 else f"{a:.2f} or more")
           for i, (a, b) in enumerate(zip(cortes[:-1], cortes[1:]))]
    base(ax_b, capas, "Uncertainty of the prediction", "b)", ley, "Width of the calibrated 90% interval", ncol=2)

    # c: decisión
    nc = d["nc_riesgo"].astype(bool); err = d["error"] == 1; pos = d["yhat"] == 1
    clases = [(~nc & ~err & pos, "#2b8a3e", "Potato, correct"), (~nc & ~err & ~pos, "#d9d4c7", "No potato, correct"),
              (~nc & err, COL["nominal"], "Wrong"), (nc, "#5f5e5a", "Inconclusive")]
    ley = [Patch(fc=c, label=f"{t} ({int(m.sum())})") for m, c, t in clases]
    base(ax_c, [(m, c) for m, c, _ in clases], "Decision of the map", "c)", ley, "Decision with abstention (districts)", ncol=2)
    guardar(fig, "Fig6_mapa_abstencion")


CIUDADES = {"Lima": (-77.043, -12.046, "left"), "Cajamarca": (-78.517, -7.162, "right"), "Huancayo": (-75.205, -12.065, "right"),
            "Cusco": (-71.968, -13.532, "right"), "Arequipa": (-71.537, -16.409, "left"), "Puno": (-70.020, -15.840, "right"),
            "Iquitos": (-73.254, -3.744, "right")}


def _tinte_relieve(ax, alfa_tinte=0.85):
    """Tinte hipsométrico (Copernicus GLO-90) con sombreado encima, recortado al Perú."""
    import rasterio
    from matplotlib.colors import LinearSegmentedColormap
    with rasterio.open(RAIZ / "data" / "elevacion_peru.tif") as r:
        z = r.read(1, masked=True); b = r.bounds
    with rasterio.open(RAIZ / "data" / "relieve_peru.tif") as r:
        sh = r.read(1, masked=True)
    hip = LinearSegmentedColormap.from_list("hip", ["#eef3e2", "#dfe8c8", "#ece2bf", "#e2cfa6", "#d4bfa0", "#ece8e3", "#ffffff"])
    norma = mpl.colors.Normalize(0, 5500)
    rgba = hip(norma(z.filled(np.nan))); rgba[..., 3] = np.where(z.mask, 0, alfa_tinte)
    ext = (b.left, b.right, b.bottom, b.top)
    ax.imshow(rgba, extent=ext, zorder=1, interpolation="bilinear", rasterized=True)
    capa = np.zeros(sh.shape + (4,)); capa[..., 3] = np.nan_to_num(np.clip((0.80 - sh.filled(np.nan)) * 0.9, 0, 0.5))
    ax.imshow(capa, extent=ext, zorder=2, interpolation="bilinear", rasterized=True)
    return hip, norma


def _marco_geografico(ax, crs):
    """Marcas de longitud y latitud en el borde del mapa (proyección UTM 18S)."""
    from pyproj import Transformer
    t = Transformer.from_crs(4326, crs, always_xy=True)
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    xs = [t.transform(lon, -10)[0] for lon in (-80, -76, -72)]
    ys = [t.transform(-75, lat)[1] for lat in (-4, -8, -12, -16)]
    ax.set_xticks(xs, ["80°W", "76°W", "72°W"]); ax.set_yticks(ys, ["4°S", "8°S", "12°S", "16°S"])
    ax.tick_params(labelsize=5.8, length=2, width=0.5, colors=TINTA2)
    for lado in ax.spines.values():
        lado.set_visible(True); lado.set_color(TINTA2); lado.set_linewidth(0.6)
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    ax.set_xlabel(""); ax.set_ylabel("")


def _escala_norte(ax, fx=0.06, fy=0.04):
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    xs, ys = x0 + fx * (x1 - x0), y0 + fy * (y1 - y0)
    for k, c in zip(range(3), (TINTA, "white", TINTA)):
        ax.add_patch(mpl.patches.Rectangle((xs + k * 100e3, ys), 100e3, 14e3, fc=c, ec=TINTA, lw=0.5, zorder=8))
    for k, t in ((0, "0"), (300, "300 km")):
        ax.text(xs + k * 1e3, ys + 26e3, t, ha="center", va="bottom", fontsize=6.6, color=TINTA, zorder=8)
    ax.annotate("N", xy=(x1 - 0.07 * (x1 - x0), y1 - 0.035 * (y1 - y0)), xytext=(x1 - 0.07 * (x1 - x0), y1 - 0.11 * (y1 - y0)),
                ha="center", va="center", fontsize=8.4, fontweight="bold", color=TINTA, zorder=8,
                arrowprops=dict(arrowstyle="-|>,head_length=0.6,head_width=0.3", color=TINTA, lw=1.0))


def fig1():
    """Área de estudio: papa en el censo sobre el relieve (a), sierra y bloques de validación (b)
    y fracción de papa frente a la altitud (c)."""
    import geopandas as gpd
    import matplotlib.patheffects as pe
    from matplotlib.gridspec import GridSpec
    from pyproj import Transformer
    crs = 32718
    d = gpd.read_file(RAIZ / "data" / "limites" / "distritos_2012.gpkg").to_crs(crs)
    tab = pd.read_csv(RAIZ / "data" / "tabla_distrital_2012.csv", dtype={"ubigeo": str})
    pl = pd.read_csv(RES / "pliegues.csv", dtype={"ubigeo": str})
    d = d.merge(tab[["ubigeo", "registros_papa", "registros_cultivo", "ha_papa", "altitud_media_m", "sierra"]], on="ubigeo")
    d = d.merge(pl, on="ubigeo")
    d["frac"] = d["registros_papa"] / d["registros_cultivo"]
    d["dep"] = d["ubigeo"].str[:2]
    pais = d.dissolve(); deps = d.dissolve(by="dep"); bloques = d.dissolve(by="pliegue")
    pais["geometry"] = pais.simplify(400); deps["geometry"] = deps.simplify(400); bloques["geometry"] = bloques.simplify(400)
    x0, y0, x1, y1 = pais.total_bounds
    ciudades = {k: Transformer.from_crs(4326, crs, always_xy=True).transform(lon, lat) + (lado,)
                for k, (lon, lat, lado) in CIUDADES.items()}
    halo = [pe.withStroke(linewidth=1.6, foreground="white")]

    NOMBRES = {"01": "Amazonas", "02": "Áncash", "03": "Apurímac", "04": "Arequipa", "05": "Ayacucho",
               "06": "Cajamarca", "08": "Cusco", "09": "Huancavelica", "10": "Huánuco", "11": "Ica",
               "12": "Junín", "13": "La Libertad", "14": "Lambayeque", "15": "Lima", "16": "Loreto",
               "17": "Madre de Dios", "18": "Moquegua", "19": "Pasco", "20": "Piura", "21": "Puno",
               "22": "San Martín", "23": "Tacna", "24": "Tumbes", "25": "Ucayali"}
    halo_dep = [pe.withStroke(linewidth=1.8, foreground="white")]

    fig = plt.figure(figsize=(190 * MM, 200 * MM))
    gs = GridSpec(2, 2, width_ratios=[1.75, 1], height_ratios=[1, 1], wspace=0.06, hspace=0.42, figure=fig)

    # a: mapa principal grande; burbujas de superficie de papa coloreadas por la fracción de registros con papa
    ax = fig.add_subplot(gs[:, 0])
    ax.set_xlim(x0 - 25e3, x1 + 25e3); ax.set_ylim(y0 - 25e3, y1 + 25e3); ax.set_aspect("equal")
    hip, norma = _tinte_relieve(ax)
    deps.boundary.plot(ax=ax, color="white", lw=0.7, zorder=3)
    pais.boundary.plot(ax=ax, color=TINTA2, lw=0.9, zorder=4)
    p = d[d["ha_papa"] > 0].copy(); pt = p.representative_point()
    p = p.assign(px=pt.x, py=pt.y).sort_values("ha_papa", ascending=False)
    cortes = [0, 0.05, 0.15, 0.30, 0.45, 1.01]
    paleta = sns.color_palette("rocket_r", len(cortes) - 1)
    cls = np.digitize(p["frac"], cortes) - 1
    tam = lambda ha: 2.0 + 0.050 * np.power(ha, 0.88)
    ax.scatter(p["px"], p["py"], s=tam(p["ha_papa"].to_numpy()), c=[paleta[i] for i in cls], edgecolor="white",
               linewidth=0.3, alpha=0.9, zorder=5, rasterized=True)
    for cod, g in deps.iterrows():
        if cod not in NOMBRES: continue
        c = g.geometry.representative_point()
        ax.text(c.x, c.y, NOMBRES[cod], ha="center", va="center", fontsize=6.4, fontstyle="italic",
                color="#3d3c39", path_effects=halo_dep, zorder=6)
    lx, ly, _ = ciudades["Lima"]
    ax.plot(lx, ly, "*", ms=7, mfc="#f0e442", mec=TINTA, mew=0.6, zorder=7)
    _marco_geografico(ax, crs); _escala_norte(ax, fx=0.40, fy=0.02)
    leg1 = ax.legend(handles=[Line2D([], [], marker="o", lw=0, ms=5.5, mfc=c, mec="white", mew=0.3,
                                     label=f"{100 * a:.0f}–{100 * b:.0f}%" if b < 1 else f"{100 * a:.0f}% or more")
                              for a, b, c in zip(cortes[:-1], cortes[1:], paleta)],
                     title="Crop records with potato", loc="lower left", bbox_to_anchor=(0.0, 0.0), ncol=1,
                     frameon=False, fontsize=7.2, title_fontsize=7.6, alignment="left", handletextpad=0.3)
    ax.add_artist(leg1)
    ax.legend(handles=[plt.scatter([], [], s=tam(h), color="#bdbdbd", edgecolor=TINTA2, linewidth=0.4, label=f"{h:,} ha")
                       for h in (500, 2000, 6000)], title="Potato area, 2012", loc="lower left", bbox_to_anchor=(0.0, 0.205), ncol=1,
              frameon=False, fontsize=7.2, title_fontsize=7.6, labelspacing=1.3, alignment="left", handletextpad=0.8, borderpad=0.6)
    ax.set_title("Potato in the 2012 census", loc="left", fontsize=9.6, pad=5)
    ax.annotate("a)", xy=(0, 1), xycoords="axes fraction", xytext=(-6, 4), textcoords="offset points",
                fontsize=10.8, fontweight="bold", ha="right", va="bottom")

    # b: sierra y bloques de validación sobre el mismo relieve
    ax = fig.add_subplot(gs[0, 1])
    ax.set_xlim(x0 - 25e3, x1 + 25e3); ax.set_ylim(y0 - 25e3, y1 + 25e3); ax.set_aspect("equal")
    _tinte_relieve(ax, alfa_tinte=0.55)
    d[d["sierra"] == 1].plot(ax=ax, color=COL["bloque"], alpha=0.38, lw=0, zorder=3, rasterized=True)
    bloques.boundary.plot(ax=ax, color=TINTA, lw=0.45, zorder=4)
    pais.boundary.plot(ax=ax, color=TINTA2, lw=0.7, zorder=5)
    _marco_geografico(ax, crs); ax.set_xticklabels([]); ax.yaxis.tick_right()
    sm = mpl.cm.ScalarMappable(norm=norma, cmap=hip)
    cax = ax.inset_axes([0.05, -0.10, 0.45, 0.03])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal"); cb.outline.set_linewidth(0.4)
    cb.set_ticks([0, 2300, 5500]); cb.ax.tick_params(labelsize=5.5, length=2)
    cb.ax.axvline(2300, color=COL["bloque"], lw=1.2)
    cb.set_label("Elevation (m)", fontsize=6.6, labelpad=1)
    ax.legend(handles=[Patch(fc=COL["bloque"], alpha=0.38, ec="none", label="Highlands (2300 m or more)"),
                       Line2D([], [], color=TINTA, lw=0.8, label="145 km blocks (20 folds)")],
              loc="upper left", bbox_to_anchor=(0.0, -0.17), frameon=False, fontsize=6.6)
    ax.set_title("Highlands and validation blocks", loc="left", fontsize=8.4, pad=4)
    ax.annotate("b)", xy=(0, 1), xycoords="axes fraction", xytext=(-6, 4), textcoords="offset points",
                fontsize=10.8, fontweight="bold", ha="right", va="bottom")

    # c: perfil de altitud (altitud en el eje vertical, como el terreno)
    ax = fig.add_subplot(gs[1, 1])
    hb = ax.hexbin(d["frac"], d["altitud_media_m"], gridsize=(14, 22), cmap="mako_r", mincnt=1, linewidths=0.2,
                   edgecolors="white", extent=(0, 0.8, 0, 5000))
    ax.axhline(2300, color=COL["bloque"], lw=0.9, ls=(0, (3, 2)))
    ax.text(0.78, 2380, "Highlands", ha="right", va="bottom", fontsize=6.96, color=COL["bloque"])
    ax.set_xlim(0, 0.8); ax.set_ylim(0, 5000)
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{100 * v:.0f}%"))
    ax.set_xlabel("Crop records with potato"); ax.set_ylabel("Mean elevation of agricultural land (m)")
    ax.yaxis.set_label_position("right"); ax.yaxis.tick_right()
    cax = ax.inset_axes([0.42, 0.10, 0.52, 0.03])
    cb = fig.colorbar(hb, cax=cax, orientation="horizontal"); cb.outline.set_linewidth(0.4)
    cb.ax.tick_params(labelsize=5.5, length=2); cb.set_label("Districts per cell", fontsize=6.6, labelpad=1)
    ax.set_title("Potato rises with elevation", loc="left", fontsize=8.4, pad=4)
    ax.annotate("c)", xy=(0, 1), xycoords="axes fraction", xytext=(-6, 4), textcoords="offset points",
                fontsize=10.8, fontweight="bold", ha="right", va="bottom")
    sns.despine(ax=ax, left=True, right=False)

    for ext in ("pdf", "tiff", "png"):
        kw = {"pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tiff" else ({"dpi": 300} if ext == "pdf" else {})
        fig.savefig(FIG / f"Fig1_area_estudio.{ext}", bbox_inches="tight", pad_inches=0.02, facecolor="white", **kw)
    plt.close(fig); print("guardada:", FIG / "Fig1_area_estudio")


if __name__ == "__main__":
    pedidas = sys.argv[1:] or ["fig1", "fig2", "fig6"]   # las figuras 3, 4 y 5 se hacen con ggplot2 (05_figuras.R)
    for f in pedidas:
        globals()[f]()
