from copy import deepcopy
from typing import Dict
import numpy as np
import matplotlib
from sbmlsim.plot import Axis, Figure, Plot
from sbmlsim.simulation import Timecourse, TimecourseSim
from pkdb_models.models.dapagliflozin.experiments.base_experiment import DapagliflozinSimulationExperiment
from pkdb_models.models.dapagliflozin.helpers import run_experiments


class HepaticRenalImpairment(DapagliflozinSimulationExperiment):
    """Tests hepatic and renal impairment."""

    maps = {
        "hepatic": DapagliflozinSimulationExperiment.cirrhosis_map,
        "renal": DapagliflozinSimulationExperiment.renal_map,
    }
    parameters = {
        "hepatic": "f_cirrhosis",
        "renal": "KI__f_renal_function",
    }
    colors = {
        "hepatic": DapagliflozinSimulationExperiment.cirrhosis_colors,
        "renal": DapagliflozinSimulationExperiment.renal_colors,
    }
    impairments = list(maps.keys())
    dpi = 600
    legend_font_size = 10

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        for impairment in self.impairments:
            map = self.maps[impairment]
            parameter = self.parameters[impairment]
            for group, value in map.items():
                tcsims[f"dap_{impairment}_{group}"] = TimecourseSim(
                    Timecourse(
                        start=0,
                        end=36 * 60,  # [min]
                        steps=5000,
                        changes={
                            **self.default_changes(),
                            f"PODOSE_dap": Q_(10, "mg"),
                            parameter: Q_(value, "dimensionless")
                        },
                    )
                )
        return tcsims

    def figures(self) -> Dict[str, Figure]:
        return {
            **self.figure_pk(),
            **self.figure_pd(),
        }

    def figure_pk(self) -> Dict[str, Figure]:
        figures = {}
        Figure.legend_fontsize = self.legend_font_size
        Figure.fig_dpi = self.dpi
        for impairment in self.impairments:
            map = self.maps[impairment]
            parameter = self.parameters[impairment]
            colors = self.colors[impairment]
            fig = Figure(
                experiment=self,
                sid=f"Fig_{parameter}_pk",
                num_rows=1,
                num_cols=5,
                name=f"Pharmacokinetics: {impairment.title()} impairment",
            )
            plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
            # plots[7].xaxis = None
            sids = [
                # plasma
                "[Cve_dap]",
                "[Cve_d3g]",
                #"[Cve_daptot]",
                # urine
                "Aurine_dap",
                "Aurine_d3g",
                #"Aurine_daptot",
                # feces
                "Afeces_dap",
                # None,
                # "Afeces_daptot",
            ]

            for ksid, sid in enumerate(sids):
                if not sid:
                    continue
                plots[ksid].set_yaxis(label=self.labels[sid], unit=self.units[sid])
                for group, value in map.items():
                    plots[ksid].add_data(
                        task=f"task_dap_{impairment}_{group}",
                        xid="time",
                        yid=sid,
                        label=f"{group}", # if ksid == 0 else None,
                        color=colors[group],
                    )
            figures[fig.sid] = fig
        return figures

    def figure_pd(self) -> Dict[str, Figure]:
        Figure.legend_fontsize = self.legend_font_size
        Figure.fig_dpi = self.dpi
        figures = {}
        for impairment in self.impairments:
            map = self.maps[impairment]
            parameter = self.parameters[impairment]
            colors = self.colors[impairment]
            fig = Figure(
                experiment=self,
                sid=f"Fig_{parameter}_pd",
                num_rows=1,
                num_cols=5,
                name=f"Pharmacodynamics: {impairment.title()} impairment",
            )
            plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
            sids = [
                #"[KI__glc_ext]",
                "KI__RTG",
                "KI__UGE",
            ]
            for ksid, sid in enumerate(sids):
                plots[ksid].set_yaxis(label=self.labels[sid], unit=self.units[sid])
            for ksid, sid in enumerate(sids):
                plots[ksid].set_yaxis(label=self.labels[sid], unit=self.units[sid])
                for group, value in map.items():
                    # simulations
                    plots[ksid].add_data(
                        task=f"task_dap_{impairment}_{group}",
                        xid="time",
                        yid=sid,
                        label=f"{group}",
                        color=colors[group],
                    )
            for ksid in range(2, 5):
                plots[ksid].set_xaxis(label=self.label_dap_plasma, unit=self.unit_dap)
            plots[2].set_yaxis(label=self.label_rtg, unit=self.unit_rtg)
            plots[3].set_yaxis(label=self.labels["KI__GLCEX"], unit=self.units["KI__GLCEX"])
            plots[4].set_yaxis(label=self.label_uge, unit=self.unit_uge)
            for group, value in map.items():
                # simulations
                plots[2].add_data(
                    task=f"task_dap_{impairment}_{group}",
                    xid="[Cve_dap]",
                    yid="KI__RTG",
                    label=f"{group}",
                    color=colors[group],
                )
                plots[3].add_data(
                    task=f"task_dap_{impairment}_{group}",
                    xid="[Cve_dap]",
                    yid="KI__GLCEX",
                    label=f"{group}",
                    color=colors[group],
                )
                plots[4].add_data(
                    task=f"task_dap_{impairment}_{group}",
                    xid="[Cve_dap]",
                    yid="KI__UGE",
                    label=f"{group}",
                    color=colors[group],
                )
            figures[fig.sid] = fig
        return figures


if __name__ == "__main__":
    run_experiments(HepaticRenalImpairment, output_dir=HepaticRenalImpairment.__name__)