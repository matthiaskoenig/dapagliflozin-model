"""Parameter scan over oral dose."""

from typing import Dict, List
import numpy as np
import matplotlib
import pandas as pd
from matplotlib.ticker import ScalarFormatter
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap, BoundaryNorm, LinearSegmentedColormap, LogNorm
from sbmlsim.simulation import Timecourse, TimecourseSim, ScanSim, Dimension
from sbmlsim.plot.serialization_matplotlib import FigureMPL
from sbmlsim.plot.serialization_matplotlib import plt
from sbmlutils.console import console
from pathlib import Path
from pkdb_models.models.dapagliflozin.experiments.base_experiment import DapagliflozinSimulationExperiment
from pkdb_models.models.dapagliflozin.helpers import run_experiments


class DapagliflozinDoseScan(DapagliflozinSimulationExperiment):
    """Scan the effect of oral dose on pharmacokinetics and pharmacodynamics."""

    study_params_tsv = Path(__file__).resolve().parent.parent.parent / "data/parameters/parameters_dose.tsv"

    tend = 2 * 24 * 60    # [min]
    steps = 5000
    dose_values_mg: List[float] = [1.0, 2.5, 5, 10, 20, 25, 50, 100, 250]
    sids = [
        "[Cve_dap]",
        "[Cve_d3g]",
        "Aurine_dap",
        "Aurine_d3g",
        "Afeces_dap",
    ]

    def _get_dose_color(self, dose: float):
        if hasattr(self, "dose_colors") and dose in self.dose_colors:
            return self.dose_colors[dose]
        cmap = matplotlib.colormaps.get("Oranges")
        x = 0.0 if dose <= 0 else np.log10(dose) / np.log10(max(self.dose_values_mg))
        return cmap(np.clip(x, 0.0, 1.0))

    def simulations(self) -> Dict[str, ScanSim]:
        Q_ = self.Q_
        sim = ScanSim(
            simulation=TimecourseSim(
                Timecourse(
                    start=0,
                    end=self.tend,
                    steps=self.steps,
                    changes={**self.default_changes()},
                )
            ),
            dimensions=[
                Dimension(
                    "dim_dose",
                    changes={"PODOSE_dap": Q_(self.dose_values_mg, "mg")},
                )
            ],
        )
        return {"scan_po_dose": sim}

    def figures_mpl(self) -> Dict[str, FigureMPL]:
        self.pk_dfs = self.calculate_dapagliflozin_pk()
        return {
            # **self.figures_mpl_timecourses(),
            **self.figures_mpl_pharmacokinetics(),
        }

    def figures_mpl_timecourses(self) -> Dict[str, FigureMPL]:
        figures: Dict[str, FigureMPL] = {}
        scan_key = "dose"
        f, axes = plt.subplots(
            nrows=1,
            ncols=len(self.sids),
            figsize=(6 * len(self.sids), 6),
            dpi=600,
            layout="constrained",
        )

        Q_ = self.Q_
        sim_key = "scan_po_dose"
        xres = self.results[f"task_{sim_key}"]
        t_vec = xres.dim_mean("time").to(self.units["time"])
        scandim = xres._redop_dims()[0]
        dose_vec = Q_(xres["PODOSE_dap"].values[0], xres.uinfo["PODOSE_dap"]).to("mg")
        dose_vals = [float(d.magnitude) for d in dose_vec]
        discrete_colors = [self._get_dose_color(d) for d in dose_vals]
        cmap = LinearSegmentedColormap.from_list("dose_cmap", discrete_colors, N=256)
        norm = LogNorm(vmin=min(dose_vals), vmax=max(dose_vals))

        for kcol, sid in enumerate(self.sids):
            ax = axes[kcol]
            for k_dose, dose in enumerate(dose_vec):
                color = cmap(norm(float(dose.magnitude)))
                c_vec = Q_(
                    xres[sid].sel({scandim: k_dose}).values,
                    xres.uinfo[sid],
                ).to(self.units[sid])
                ax.plot(
                    t_vec.magnitude,
                    c_vec.magnitude,
                    color=color,
                    linewidth=2.0,
                )
            ax.set_xlabel(f"{self.label_time} [{self.units['time']}]", fontdict=self.font)
            ax.set_ylabel(f"{self.labels[sid]} [{self.units[sid]}]", fontdict=self.font)
            ax.tick_params(axis="x", labelsize=self.tick_font_size)
            ax.tick_params(axis="y", labelsize=self.tick_font_size)

            if sid in ["[Cve_dap]", "[Cve_d3g]"]:
                ax.set_xlim(0, 21)

        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        cax = inset_axes(
            axes[0],
            width="75%",
            height="10%",
            loc="upper center",
            bbox_to_anchor=(0, 0, 1, 1),
            bbox_transform=axes[0].transAxes,
            borderpad=0.5,
        )

        mappable = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
        cbar = f.colorbar(mappable, cax=cax, orientation="horizontal")

        ticks = [1, 10, 50, 500]
        cbar.set_ticks(ticks)
        cbar.set_ticklabels([str(t) for t in ticks])
        cbar.ax.tick_params(labelsize=self.tick_font_size)
        cbar.ax.tick_params(labelsize=self.tick_font_size)

        figures[f"timecourse__{scan_key}"] = f
        return figures

    def figures_mpl_pharmacokinetics(self) -> Dict[str, FigureMPL]:
        figures: Dict[str, FigureMPL] = {}
        Q_ = self.Q_
        sim_key = "scan_po_dose"
        xres = self.results[f"task_{sim_key}"]
        df_all = self.pk_dfs[sim_key]

        # doses from sim
        dose_vec = Q_(xres["PODOSE_dap"].values[0], xres.uinfo["PODOSE_dap"]).to("mg")
        dose_vals = dose_vec.magnitude.tolist()

        # time vector and index at 24h
        t_vec = xres.dim_mean("time").to(self.units["time"]).magnitude  # minutes
        t24_min = 24 * 60
        idx_24h = int(np.argmin(np.abs(t_vec - t24_min)))

        # figure + axes
        f, axes = plt.subplots(nrows=1, ncols=5, figsize=(30, 6), dpi=600, layout="constrained")
        axes = axes.flatten()

        # styling
        series_colors = {"dap": "black", "d3g": "darkgrey"}
        marker_size_pts = 10.0
        scatter_s = marker_size_pts ** 2
        edge_width = 2.0
        line_width = 2.5

        # collect legend handles
        series_handles = [
            Line2D([0], [0], marker="o", linestyle="",
                   markerfacecolor="white", markeredgecolor=series_colors["dap"],
                   markeredgewidth=edge_width, markersize=marker_size_pts, label="dap"),
            Line2D([0], [0], marker="o", linestyle="",
                   markerfacecolor="white", markeredgecolor=series_colors["d3g"],
                   markeredgewidth=edge_width, markersize=marker_size_pts, label="d3g"),
        ]

        for k, pk_key in enumerate(["aucinf", "cmax", "kel", "thalf"]):
            ax = axes[k]
            for substance in ["dap", "d3g"]:
                df = df_all[df_all.substance == substance].copy()
                yq = Q_(df[pk_key].to_numpy(), df[f"{pk_key}_unit"].values[0]).to(self.pk_units[pk_key])
                y = yq.magnitude
                ax.plot(dose_vals, y, linestyle="-", linewidth=line_width,
                        color=series_colors[substance], alpha=1.0, zorder=2)
                facecols = [self._get_dose_color(d) for d in dose_vals]
                ax.scatter(dose_vals, y, s=scatter_s, facecolors=facecols,
                           edgecolors=series_colors[substance], linewidths=edge_width, zorder=3)

            ax.set_xlabel("Oral dose [mg]", fontdict=self.scan_font)
            ax.set_ylabel(f"{self.pk_labels[pk_key]} [{self.pk_units[pk_key]}]", fontdict=self.scan_font)
            ax.xaxis.set_major_formatter(ScalarFormatter())
            ax.ticklabel_format(style="plain", axis="x")
            ax.tick_params(axis="x", labelsize=self.tick_font_size)
            ax.tick_params(axis="y", labelsize=self.tick_font_size)

        # UGE panel (last)
        uge_sid = "KI__UGE"
        uge_unit = xres.uinfo[uge_sid]
        uge_24h_g = [
            Q_(xres[uge_sid].isel({xres._redop_dims()[0]: i}).values[idx_24h], uge_unit).to("g").magnitude
            for i, _ in enumerate(dose_vec)
        ]
        ax_uge = axes[4]
        ax_uge.plot(dose_vals, uge_24h_g, linestyle="-", linewidth=line_width, color="black", zorder=2)
        facecols = [self._get_dose_color(d) for d in dose_vals]
        ax_uge.scatter(dose_vals, uge_24h_g, s=scatter_s, facecolors=facecols,
                       edgecolors="black", linewidths=edge_width, zorder=3)
        ax_uge.set_xlabel("Oral dose [mg]", fontdict=self.scan_font)
        ax_uge.set_ylabel("UGE (24 hr) [g]", fontdict=self.scan_font)
        ax_uge.xaxis.set_major_formatter(ScalarFormatter())
        ax_uge.ticklabel_format(style="plain", axis="x")
        ax_uge.tick_params(axis="x", labelsize=self.tick_font_size)
        ax_uge.tick_params(axis="y", labelsize=self.tick_font_size)

        # overlay clinical study points
        df_st = pd.read_csv(self.study_params_tsv, sep="\t")
        df_st = df_st[["study", "parameter", "dose", "substance", "mean_unit_conv", "sd_unit_conv"]].copy()

        base_markers = ["s", "^", "v", "<", ">", "D", "d", "p", "h"]
        fillstyles = ["full", "left", "right", "top", "bottom", "none"]
        M, F = len(base_markers), len(fillstyles)

        def style_for_index(i: int):
            m = base_markers[i % M]
            fs = fillstyles[(i // M) % F]
            return m, fs

        studies = pd.unique(df_st["study"])
        study_style = {s: style_for_index(i) for i, s in enumerate(studies)}

        # study legend handles
        study_handles: Dict[str, Line2D] = {}
        legend_mfc = "#bbbbbb"  # neutral facecolor in legend

        param_ax_idx = {"aucinf": 0, "cmax": 1, "thalf": 3, "UGE": 4}

        for pk_param in ["aucinf", "cmax", "thalf"]:
            ax = axes[param_ax_idx[pk_param]]
            for substance in ["dap", "d3g"]:
                dfx = df_st[(df_st["parameter"] == pk_param) & (df_st["substance"] == substance)]
                if dfx.empty:
                    continue
                for _, row in dfx.iterrows():
                    d = float(row["dose"])
                    y = float(row["mean_unit_conv"])
                    yerr = row["sd_unit_conv"]
                    m, fs = study_style[row["study"]]
                    face = self._get_dose_color(d)
                    edge = series_colors[substance]
                    mfc = face if fs != "none" else "white"

                    if pd.notna(yerr) and float(yerr) != 0.0:
                        ax.errorbar([d], [y], yerr=[[yerr], [yerr]],
                                    fmt=m, markersize=marker_size_pts, fillstyle=fs, mfc=mfc,
                                    mec=edge, mew=edge_width, lw=1.5, ecolor=edge, capsize=3, zorder=4)
                    else:
                        ax.plot([d], [y], linestyle="", marker=m, markersize=marker_size_pts,
                                fillstyle=fs, mfc=mfc, mec=edge, mew=edge_width, zorder=4)

                    if row["study"] not in study_handles:
                        lm, lfs = study_style[row["study"]]
                        study_handles[row["study"]] = Line2D(
                            [0], [0], marker=lm, linestyle="",
                            fillstyle=lfs,
                            markerfacecolor=(legend_mfc if lfs != "none" else "white"),
                            markeredgecolor="black",
                            markeredgewidth=edge_width,
                            markersize=marker_size_pts,
                            label=row["study"],
                        )

        # UGE study overlay
        dfx = df_st[df_st["parameter"] == "UGE"]
        if not dfx.empty:
            for _, row in dfx.iterrows():
                d = float(row["dose"])
                y = float(row["mean_unit_conv"])
                yerr = row["sd_unit_conv"]
                m, fs = study_style[row["study"]]
                face = self._get_dose_color(d)
                mfc = face if fs != "none" else "white"

                if pd.notna(yerr) and float(yerr) != 0.0:
                    ax_uge.errorbar([d], [y], yerr=[[yerr], [yerr]],
                                    fmt=m, markersize=marker_size_pts, fillstyle=fs, mfc=mfc,
                                    mec="black", mew=edge_width, lw=1.5, ecolor="black", capsize=3, zorder=4)
                else:
                    ax_uge.plot([d], [y], linestyle="",
                                marker=m, markersize=marker_size_pts,
                                fillstyle=fs, mfc=mfc, mec="black", mew=edge_width, zorder=4)

                if row["study"] not in study_handles:
                    lm, lfs = study_style[row["study"]]
                    study_handles[row["study"]] = Line2D(
                        [0], [0], marker=lm, linestyle="",
                        fillstyle=lfs,
                        markerfacecolor=(legend_mfc if lfs != "none" else "white"),
                        markeredgecolor="black",
                        markeredgewidth=edge_width,
                        markersize=marker_size_pts,
                        label=row["study"],
                    )

        # LEGENDS
        if study_handles:
            combined_handles = series_handles + list(study_handles.values())
        else:
            combined_handles = series_handles

        axes[0].legend(
            handles=combined_handles,
            loc="best",
            frameon=True,
            fontsize=10,
            ncol=2
        )

        figures["pk_dose"] = f
        return figures


if __name__ == "__main__":
    run_experiments(DapagliflozinDoseScan, output_dir=DapagliflozinDoseScan.__name__)