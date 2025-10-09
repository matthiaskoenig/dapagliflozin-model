"""Parameter scans dapagliflozin."""
from typing import Dict, List
import numpy as np
import matplotlib
import matplotlib.axes
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap, Normalize, LogNorm
from matplotlib.ticker import ScalarFormatter
from matplotlib.lines import Line2D
import pandas as pd
from pathlib import Path

from sbmlsim.simulation import Timecourse, TimecourseSim, ScanSim, Dimension
from sbmlsim.plot.serialization_matplotlib import FigureMPL, MatplotlibFigureSerializer
from sbmlsim.plot.serialization_matplotlib import plt
from sbmlutils.console import console

from pkdb_models.models.dapagliflozin.experiments.base_experiment import DapagliflozinSimulationExperiment
from pkdb_models.models.dapagliflozin.helpers import run_experiments


class DapagliflozinParameterScan(DapagliflozinSimulationExperiment):
    """Scan the effect of parameters on pharmacokinetics and pharmacodynamics."""

    parameters_tsv = {
        "food_scan":    Path(__file__).resolve().parent.parent.parent / "data/parameters/parameters_food.tsv",
        "hepatic_scan": Path(__file__).resolve().parent.parent.parent / "data/parameters/parameters_hepatic.tsv",
        "renal_scan":   Path(__file__).resolve().parent.parent.parent / "data/parameters/parameters_renal.tsv",
        "renal_scan_20":Path(__file__).resolve().parent.parent.parent / "data/parameters/parameters_renal.tsv",
    }

    SCAN_CLASS = {
        "food_scan": "food",
        "renal_scan": "renal function",
        "renal_scan_20": "renal function",
        "hepatic_scan": "hepatic function",
    }

    SCAN_MAPPING = {
        "food": {
            "fasted": "fasted",
            "fed": "fed",
        },
        "renal function": {
            "healthy": "Normal renal function",
            "normal": "Normal renal function",
            "mild": "Mild renal impairment",
            "moderate": "Moderate renal impairment",
            "severe": "Severe renal impairment",
        },
        "hepatic function": {
            "healthy": "Control",
            "mild": "Mild cirrhosis",
            "moderate": "Moderate cirrhosis",
            "severe": "Severe cirrhosis",
        },
    }

    tend = 2.0 * 24 * 60
    steps = 4000
    dose_dap = 10
    num_points = 15

    scan_map = {
        "hepatic_scan": {
            "parameter": "f_cirrhosis",
            "default": 0.0,
            "range": np.linspace(0, 0.9, num=num_points),
            "scale": "linear",
            "units": "dimensionless",
            "label": "cirrhosis degree [-]",
            "dose_mg": 10,
        },
        "renal_scan": {
            "parameter": "KI__f_renal_function",
            "default": 1.0,
            "range": np.sort(np.append(np.logspace(-1, 1, num=num_points), [1.0])),  # 0.1..10, include 1
            "scale": "log",
            "units": "dimensionless",
            "label": "renal function [-]",
            "dose_mg": 50,  # PK overlays at 50 mg
        },
        "renal_scan_20": {
            "parameter": "KI__f_renal_function",
            "default": 1.0,
            "range": np.sort(np.append(np.logspace(-1, 1, num=num_points), [1.0])),
            "scale": "log",
            "units": "dimensionless",
            "label": "renal function [-]",
            "dose_mg": 20,  # UGE overlays at 20 mg
        },
        "food_scan": {
            "parameter": "GU__f_absorption",
            "default": 1.0,
            "range": np.sort(np.append(np.logspace(-1, 1, num=num_points), [1.0])),
            "scale": "log",
            "colormap": "seismic_r",
            "units": "dimensionless",
            "label": "absorption activity [-]",
            "dose_mg": 10,
        },
    }

    def simulations(self) -> Dict[str, ScanSim]:
        Q_ = self.Q_
        tcscans = {}
        for scan_key, scan_data in self.scan_map.items():
            base_dose = scan_data.get("dose_mg", self.dose_dap)
            tcscans[f"scan_po_{scan_key}"] = ScanSim(
                simulation=TimecourseSim(
                    Timecourse(
                        start=0,
                        end=self.tend,
                        steps=self.steps,
                        changes={
                            **self.default_changes(),
                            "PODOSE_dap": Q_(base_dose, "mg"),
                        },
                    )
                ),
                dimensions=[
                    Dimension(
                        "dim_scan",
                        changes={
                            scan_data["parameter"]: Q_(
                                scan_data["range"], scan_data["units"]
                            )
                        },
                    ),
                ],
            )
        return tcscans

    def figures_mpl(self) -> Dict[str, FigureMPL]:
        """Matplotlib figures."""
        self.pk_dfs = self.calculate_dapagliflozin_pk()
        return {
            # **self.figures_mpl_timecourses(),
            **self.figures_mpl_pharmacokinetics(),
        }

    def figures_mpl_timecourses(self) -> Dict[str, FigureMPL]:
        sids = [
            "[Cve_dap]",
            "[Cve_d3g]",
            "Aurine_dap",
            "Aurine_d3g",
            "Afeces_dap",
        ]

        figures = {}
        for scan_key, scan_data in self.scan_map.items():
            range_vals = scan_data["range"]
            rmin, rmax = range_vals[0], range_vals[-1]

            if scan_key == "hepatic_scan":
                cmap = LinearSegmentedColormap.from_list(
                    "cirrhosis_blues",
                    [
                        self.cirrhosis_colors["Mild cirrhosis"],
                        self.cirrhosis_colors["Moderate cirrhosis"],
                        self.cirrhosis_colors["Severe cirrhosis"],
                    ],
                )
                norm = Normalize(vmin=0.0, vmax=0.9, clip=False)
            elif scan_key in ["renal_scan", "renal_scan_20"]:
                vmin, vmax = 0.1, 10.0

                def _pos(v):
                    return (np.log10(v) - np.log10(vmin)) / (np.log10(vmax) - np.log10(vmin))

                severe = self.renal_colors["Severe renal impairment"]       # ~0.1
                moderate = self.renal_colors["Moderate renal impairment"]   # ~0.32
                mild = self.renal_colors["Mild renal impairment"]           # ~0.69
                high_ext = "#e5f5e0"  # very light green for >1.0 extension

                cmap = LinearSegmentedColormap.from_list(
                    "renal_greens",
                    [
                        (_pos(0.1), severe),
                        (_pos(0.32), moderate),
                        (_pos(0.69), mild),
                        (_pos(10.0), high_ext),
                    ],
                )
                norm = LogNorm(vmin=vmin, vmax=vmax, clip=False)
            else:
                cmap_str = scan_data["colormap"]
                cmap = matplotlib.colormaps.get_cmap(cmap_str)
                if scan_data["scale"] == "linear":
                    norm = Normalize(vmin=rmin, vmax=rmax, clip=False)
                else:
                    norm = LogNorm(vmin=rmin, vmax=rmax, clip=False)

            f, axes = plt.subplots(
                nrows=1,
                ncols=len(sids),
                figsize=(6 * len(sids), 6),
                dpi=600,
                layout="constrained"
            )

            for kcol, sid in enumerate(sids):
                ax = axes[kcol]
                Q_ = self.Q_
                xres = self.results[f"task_scan_po_{scan_key}"]
                scandim = xres._redop_dims()[0]
                parameter_id = scan_data["parameter"]
                par_vec = Q_(xres[parameter_id].values[0], xres.uinfo[parameter_id])
                t_vec = xres.dim_mean("time").to(self.units["time"])
                t_vec_default = None
                c_vec_default = None
                for k_par, par in enumerate(par_vec):
                    c_vec = Q_(
                        xres[sid].sel({scandim: k_par}).values,
                        xres.uinfo[sid],
                    ).to(self.units[sid])
                    linewidth = 2.0
                    if np.isclose(scan_data["default"], par):
                        color = "black"
                        t_vec_default = t_vec
                        c_vec_default = c_vec
                    else:
                        color = cmap(norm(par.magnitude))
                    ax.plot(
                        t_vec.magnitude,
                        c_vec.magnitude,
                        color=color,
                        linewidth=linewidth,
                    )

                if t_vec_default is not None and c_vec_default is not None:
                    ax.plot(
                        t_vec_default.magnitude,
                        c_vec_default.magnitude,
                        color="black",
                        linewidth=2.0,
                    )

                ax.set_xlabel(f"{self.label_time} [{self.units['time']}]", fontdict=self.font)
                ax.tick_params(axis="x", labelsize=self.tick_font_size)
                ax.tick_params(axis="y", labelsize=self.tick_font_size)
                ax.set_ylabel(f"{self.labels[sid]} [{self.units[sid]}]", fontdict=self.font)

                if sid in ["[Cve_dap]", "[Cve_d3g]"]:
                    if scan_key in ["renal_scan", "renal_scan_20"]:
                        ax.set_xlim(0, 35)
                    else:
                        ax.set_xlim(0, 21)

            cb_ax = f.add_axes(rect=[0.08, 0.85, 0.1, 0.08])
            cb_ax.set_in_layout(True)

            cbar = f.colorbar(
                cm.ScalarMappable(norm=norm, cmap=cmap),
                cax=cb_ax,
                orientation="horizontal",
            )

            if scan_key == "hepatic_scan":
                ticks = [0.0, 0.9]
                cbar.set_ticks(ticks)
                cbar.set_ticklabels(ticks, **{"size": 15, "weight": "medium"})
                cbar.ax.axvline(x=scan_data["default"], color="black", linewidth=2)
            elif scan_key in ["renal_scan", "renal_scan_20"]:
                ticks = [0.1, 10.0]
                cbar.set_ticks(ticks)
                cbar.set_ticklabels(ticks, **{"size": 15, "weight": "medium"})
                cbar.ax.axvline(x=scan_data["default"], color="black", linewidth=2)
            else:
                rmin, rmax = scan_data["range"][0], scan_data["range"][-1]
                ticks = [rmin, rmax]
                if scan_data["default"] not in ticks:
                    ticks.append(scan_data["default"])
                    ticks = sorted(ticks)
                cbar.set_ticks(ticks)
                cbar.set_ticklabels(ticks, **{"size": 15, "weight": "medium"})
                cbar.ax.axvline(x=scan_data["default"], color="black", linewidth=2)

            cbar.ax.set_xlabel(scan_data["label"], **{"size": 15, "weight": "bold"})
            console.print(f"{scan_key} colorbar ticks: {ticks}")

            figures[f"timecourse__{scan_key}"] = f

        return figures

    def figures_mpl_pharmacokinetics(self) -> Dict[str, FigureMPL]:
        Q_ = self.Q_
        figures: Dict[str, FigureMPL] = {}

        parameters_info = {
            "dap": ["aucinf", "cmax", "kel", "thalf"],
            "d3g": ["aucinf", "cmax", "kel", "thalf"],
        }
        series_colors = {"dap": "black", "d3g": "darkgrey"}

        marker_size_pts = 10.0
        edge_width = 2.0

        series_handles = [
            Line2D([0], [0], marker="o", linestyle="", markerfacecolor="white",
                   markeredgecolor=series_colors["dap"], markeredgewidth=edge_width,
                   markersize=marker_size_pts, label="dap"),
            Line2D([0], [0], marker="o", linestyle="", markerfacecolor="white",
                   markeredgecolor=series_colors["d3g"], markeredgewidth=edge_width,
                   markersize=marker_size_pts, label="d3g"),
        ]

        marker_size_pts = 11.0
        scatter_s = marker_size_pts ** 2
        edge_width = 1.8
        line_width = 2.0

        for scan_key, scan_data in self.scan_map.items():
            range_vals = scan_data["range"]
            rmin, rmax = float(range_vals[0]), float(range_vals[-1])

            if scan_key == "hepatic_scan":
                cmap = LinearSegmentedColormap.from_list(
                    "cirrhosis_blues",
                    [
                        self.cirrhosis_colors["Mild cirrhosis"],
                        self.cirrhosis_colors["Moderate cirrhosis"],
                        self.cirrhosis_colors["Severe cirrhosis"],
                    ],
                )
                norm = Normalize(vmin=0.0, vmax=0.9, clip=False)
            elif scan_key in ["renal_scan", "renal_scan_20"]:
                vmin, vmax = 0.1, 10.0

                def _pos(v):
                    return (np.log10(v) - np.log10(vmin)) / (np.log10(vmax) - np.log10(vmin))

                severe = self.renal_colors["Severe renal impairment"]
                moderate = self.renal_colors["Moderate renal impairment"]
                mild = self.renal_colors["Mild renal impairment"]
                high_ext = "#e5f5e0"
                cmap = LinearSegmentedColormap.from_list(
                    "renal_greens",
                    [
                        (_pos(0.1), severe),
                        (_pos(0.32), moderate),
                        (_pos(0.69), mild),
                        (_pos(10.0), high_ext),
                    ],
                )
                norm = LogNorm(vmin=vmin, vmax=vmax, clip=False)
            else:
                cmap_str = scan_data["colormap"]
                cmap = matplotlib.colormaps.get_cmap(cmap_str)
                if scan_data["scale"] == "linear":
                    norm = Normalize(vmin=rmin, vmax=rmax, clip=False)
                else:
                    norm = LogNorm(vmin=rmin, vmax=rmax, clip=False)

            metrics = parameters_info["dap"]
            f, axes = plt.subplots(
                nrows=1, ncols=len(metrics) + 1, figsize=(6 * (len(metrics) + 1), 6),
                dpi=600, layout="constrained"
            )
            axes = axes.flatten()

            sim_key = f"scan_po_{scan_key}"
            xres = self.results[f"task_{sim_key}"]
            df_all = self.pk_dfs[sim_key]

            parameter_id = scan_data["parameter"]
            x_q = Q_(xres[parameter_id].values[0], xres.uinfo[parameter_id])
            x_vals = np.asarray([float(v) for v in x_q.magnitude])

            t_vec = xres.dim_mean("time").to(self.units["time"]).magnitude  # minutes
            t24_min = 24 * 60
            idx_24h = int(np.argmin(np.abs(t_vec - t24_min)))

            # model lines + scatters
            for k, pk_key in enumerate(metrics):
                ax = axes[k]
                ax.axvline(x=scan_data["default"], color="grey", linestyle="--", linewidth=1.2, zorder=1)

                for substance in ["dap", "d3g"]:
                    df = df_all[df_all.substance == substance].copy()
                    yq = Q_(df[pk_key].to_numpy(), df[f"{pk_key}_unit"].values[0]).to(self.pk_units[pk_key])
                    y = yq.magnitude

                    ax.plot(
                        x_vals, y,
                        linestyle="-", linewidth=line_width,
                        color=series_colors[substance], alpha=0.9, zorder=2,
                    )

                    facecols = [cmap(norm(x)) for x in x_vals]
                    ax.scatter(
                        x_vals, y,
                        s=scatter_s, facecolors=facecols,
                        edgecolors=series_colors[substance],
                        linewidths=edge_width, zorder=3,
                    )

                # axes formatting
                ax.tick_params(axis="x", labelsize=self.tick_font_size)
                ax.tick_params(axis="y", labelsize=self.tick_font_size)
                ax.set_xlabel(scan_data["label"], fontdict=self.scan_font)
                ax.set_ylabel(f"{self.pk_labels[pk_key]} [{self.pk_units[pk_key]}]", fontdict=self.scan_font)

                if scan_data["scale"] == "log":
                    ax.set_xscale("log")
                    ax.xaxis.set_major_formatter(ScalarFormatter())
                    ax.ticklabel_format(style="plain", axis="x")

            # UGE panel
            uge_sid = "KI__UGE"
            uge_unit = xres.uinfo[uge_sid]
            dose_dim = xres._redop_dims()[0]
            uge_24h_g = [
                Q_(xres[uge_sid].isel({dose_dim: i}).values[idx_24h], uge_unit).to("g").magnitude
                for i in range(len(x_vals))
            ]
            ax_uge = axes[-1]
            ax_uge.axvline(x=scan_data["default"], color="grey", linestyle="--", linewidth=1.2, zorder=1)
            ax_uge.plot(
                x_vals, uge_24h_g,
                linestyle="-", linewidth=line_width, color="black", zorder=2,
            )
            facecols = [cmap(norm(x)) for x in x_vals]
            ax_uge.scatter(
                x_vals, uge_24h_g,
                s=scatter_s, facecolors=facecols,
                edgecolors="black", linewidths=edge_width, zorder=3,
            )
            ax_uge.tick_params(axis="x", labelsize=self.tick_font_size)
            ax_uge.tick_params(axis="y", labelsize=self.tick_font_size)
            ax_uge.set_xlabel(scan_data["label"], fontdict=self.scan_font)
            ax_uge.set_ylabel("UGE (24 hr) [g]", fontdict=self.scan_font)
            if scan_data["scale"] == "log":
                ax_uge.set_xscale("log")
                ax_uge.xaxis.set_major_formatter(ScalarFormatter())
                ax_uge.ticklabel_format(style="plain", axis="x")

            # overlay clinical study points
            tsv_path = self.parameters_tsv[scan_key]
            df_st = pd.read_csv(tsv_path, sep="\t", comment="#")[[
                "study", "class", "condition", "parameter", "dose",
                "substance", "mean_unit_conv", "sd_unit_conv"
            ]].copy()

            df_st["parameter"] = df_st["parameter"].str.strip()
            df_st["substance"] = df_st["substance"].str.strip().str.lower()
            df_st["class"] = df_st["class"].str.strip().str.lower()
            df_st["condition"] = df_st["condition"].str.strip().str.lower()
            df_st["dose"] = pd.to_numeric(df_st["dose"], errors="coerce")

            if scan_key == "renal_scan":
                params_to_overlay = {"aucinf", "cmax", "thalf"}  # PK only
                df_st = df_st[(df_st["class"] == "renal function") & (df_st["dose"] == 50)]
            elif scan_key == "renal_scan_20":
                params_to_overlay = {"UGE"}  # UGE only
                df_st = df_st[(df_st["class"] == "renal function") & (df_st["dose"] == 20)]
            else:
                params_to_overlay = {"aucinf", "cmax", "thalf", "UGE"}

            cls = self.SCAN_CLASS[scan_key]
            num_map = {
                "food": self.fasting_map,
                "renal function": self.renal_map,
                "hepatic function": self.cirrhosis_map,
            }[cls]

            # marker styles: cycle shapes first, then fillstyles
            base_markers = ["s", "^", "v", "<", ">", "D", "d", "p", "h"]
            fillstyles = ["full", "left", "right", "top", "bottom", "none"]
            M, F = len(base_markers), len(fillstyles)

            def style_for_index(i: int):
                return base_markers[i % M], fillstyles[(i // M) % F]

            studies = pd.unique(df_st["study"])
            study_style = {s: style_for_index(i) for i, s in enumerate(studies)}
            study_handles: Dict[str, Line2D] = {}
            legend_mfc = "#bbbbbb"

            param_ax_idx = {"aucinf": 0, "cmax": 1, "thalf": 3, "UGE": len(axes) - 1}

            for pk_param in ["aucinf", "cmax", "thalf"]:
                if pk_param not in params_to_overlay:
                    continue
                ax = axes[param_ax_idx[pk_param]]
                for substance in ["dap", "d3g"]:
                    dfx = df_st[(df_st["parameter"] == pk_param) & (df_st["substance"] == substance)]
                    if dfx.empty:
                        continue
                    for _, row in dfx.iterrows():
                        cond_raw = row["condition"]
                        cond_canon = self.SCAN_MAPPING[cls][cond_raw]
                        x = float(num_map[cond_canon])
                        y = float(row["mean_unit_conv"])
                        yerr = row["sd_unit_conv"]

                        m, fs = study_style[row["study"]]
                        face = cmap(norm(x))
                        edge = series_colors[substance]
                        mfc = face if fs != "none" else "white"

                        if pd.notna(yerr) and float(yerr) != 0.0:
                            ax.errorbar(
                                [x], [y], yerr=[[yerr], [yerr]],
                                fmt=m, markersize=marker_size_pts,
                                fillstyle=fs, mfc=mfc, mec=edge, mew=edge_width,
                                lw=1.2, ecolor=edge, capsize=3, zorder=4
                            )
                        else:
                            ax.plot(
                                [x], [y],
                                linestyle="", marker=m, markersize=marker_size_pts,
                                fillstyle=fs, mfc=mfc, mec=edge, mew=edge_width, zorder=4
                            )

                        if row["study"] not in study_handles:
                            lm, lfs = study_style[row["study"]]
                            study_handles[row["study"]] = Line2D(
                                [0], [0], marker=lm, linestyle="",
                                fillstyle=lfs,
                                markerfacecolor=(legend_mfc if lfs != "none" else "white"),
                                markeredgecolor="black", markeredgewidth=edge_width,
                                markersize=marker_size_pts, label=row["study"]
                            )

            if "UGE" in params_to_overlay:
                dfx = df_st[df_st["parameter"] == "UGE"]
                if not dfx.empty:
                    for _, row in dfx.iterrows():
                        cond_raw = row["condition"]
                        cond_canon = self.SCAN_MAPPING[cls][cond_raw]
                        x = float(num_map[cond_canon])
                        y = float(row["mean_unit_conv"])
                        yerr = row["sd_unit_conv"]

                        m, fs = study_style[row["study"]]
                        face = cmap(norm(x))
                        mfc = face if fs != "none" else "white"

                        if pd.notna(yerr) and float(yerr) != 0.0:
                            ax_uge.errorbar(
                                [x], [y], yerr=[[yerr], [yerr]],
                                fmt=m, markersize=marker_size_pts,
                                fillstyle=fs, mfc=mfc, mec="black", mew=edge_width,
                                lw=1.2, ecolor="black", capsize=3, zorder=4
                            )
                        else:
                            ax_uge.plot(
                                [x], [y],
                                linestyle="", marker=m, markersize=marker_size_pts,
                                fillstyle=fs, mfc=mfc, mec="black", mew=edge_width, zorder=4
                            )

                        if row["study"] not in study_handles:
                            lm, lfs = study_style[row["study"]]
                            study_handles[row["study"]] = Line2D(
                                [0], [0], marker=lm, linestyle="",
                                fillstyle=lfs,
                                markerfacecolor=(legend_mfc if lfs != "none" else "white"),
                                markeredgecolor="black", markeredgewidth=edge_width,
                                markersize=marker_size_pts, label=row["study"]
                            )

            # LEGEND
            combined_handles = series_handles + list(study_handles.values())
            axes[0].legend(
                handles=combined_handles,
                loc="best",
                frameon=True,
                fontsize=10,
                ncol=2
            )

            figures[f"pk_{scan_key}"] = f

        return figures


if __name__ == "__main__":
    run_experiments(DapagliflozinParameterScan, output_dir=DapagliflozinParameterScan.__name__)