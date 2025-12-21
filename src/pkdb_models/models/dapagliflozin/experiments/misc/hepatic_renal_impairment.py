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
    dpi = 300
    legend_font_size = 10

    def simulations(self) -> Dict[str, TimecourseSim]:
        """
        Hepatic: single 10 mg PO
        Renal:
            - 50 mg PO (used for PK plots)
            - 20 mg PO (used for PD plots)
        """
        Q_ = self.Q_
        tcsims: Dict[str, TimecourseSim] = {}
        for impairment in self.impairments:
            map_ = self.maps[impairment]
            parameter = self.parameters[impairment]

            for group, value in map_.items():
                if impairment == "hepatic":
                    # 10 mg PO
                    key = f"dap_hepatic_{group}"
                    tcsims[key] = TimecourseSim(
                        Timecourse(
                            start=0,
                            end=36 * 60,  # [min]
                            steps=5000,
                            changes={
                                **self.default_changes(),
                                "PODOSE_dap": Q_(10, "mg"),
                                parameter: Q_(value, "dimensionless"),
                            },
                        )
                    )
                else:
                    # 50 mg PO (for PK plots)
                    key50 = f"dap_renal_{group}_po50"
                    tcsims[key50] = TimecourseSim(
                        Timecourse(
                            start=0,
                            end=36 * 60,  # [min]
                            steps=5000,
                            changes={
                                **self.default_changes(),
                                "PODOSE_dap": Q_(50, "mg"),
                                parameter: Q_(value, "dimensionless"),
                            },
                        )
                    )
                    # 20 mg PO (for PD plots)
                    key20 = f"dap_renal_{group}_po20"
                    tcsims[key20] = TimecourseSim(
                        Timecourse(
                            start=0,
                            end=36 * 60,  # [min]
                            steps=5000,
                            changes={
                                **self.default_changes(),
                                "PODOSE_dap": Q_(20, "mg"),
                                parameter: Q_(value, "dimensionless"),
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
        figures: Dict[str, Figure] = {}
        Figure.legend_fontsize = self.legend_font_size
        Figure.fig_dpi = self.dpi

        for impairment in self.impairments:
            map_ = self.maps[impairment]
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

            sids = [
                # plasma
                "[Cve_dap]",
                "[Cve_d3g]",
                # urine
                "Aurine_dap",
                "Aurine_d3g",
                # feces
                "Afeces_dap",
            ]

            for ksid, sid in enumerate(sids):
                if not sid:
                    continue
                plots[ksid].set_yaxis(label=self.labels[sid], unit=self.units[sid])

                for group, _ in map_.items():
                    if impairment == "hepatic":
                        task_id = f"task_dap_hepatic_{group}"  # 10 mg
                    else:
                        task_id = f"task_dap_renal_{group}_po50"  # 50 mg for PK

                    plots[ksid].add_data(
                        task=task_id,
                        xid="time",
                        yid=sid,
                        label=f"{group}",
                        color=colors[group],
                    )

            figures[fig.sid] = fig
        return figures

    def figure_pd(self) -> Dict[str, Figure]:
        Figure.legend_fontsize = self.legend_font_size
        Figure.fig_dpi = self.dpi
        figures: Dict[str, Figure] = {}

        for impairment in self.impairments:
            map_ = self.maps[impairment]
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

            sids_time = [
                "KI__RTG",
                "KI__UGE",
            ]
            for ksid, sid in enumerate(sids_time):
                plots[ksid].set_yaxis(label=self.labels[sid], unit=self.units[sid])

            # Exposure-response panels
            for ksid in range(2, 5):
                plots[ksid].set_xaxis(label=self.label_dap_plasma, unit=self.unit_dap)
            plots[2].set_yaxis(label=self.label_rtg, unit=self.unit_rtg)
            plots[3].set_yaxis(label=self.labels["KI__GLCEX"], unit=self.units["KI__GLCEX"])
            plots[4].set_yaxis(label=self.label_uge, unit=self.unit_uge)

            for group, _ in map_.items():
                if impairment == "hepatic":
                    task_time = f"task_dap_hepatic_{group}"    # 10 mg
                    task_er   = f"task_dap_hepatic_{group}"    # 10 mg
                else:
                    task_time = f"task_dap_renal_{group}_po20"  # 20 mg
                    task_er   = f"task_dap_renal_{group}_po20"  # 20 mg

                color = colors[group]

                plots[0].add_data(
                    task=task_time,
                    xid="time",
                    yid="KI__RTG",
                    label=f"{group}",
                    color=color,
                )
                plots[1].add_data(
                    task=task_time,
                    xid="time",
                    yid="KI__UGE",
                    label=f"{group}",
                    color=color,
                )

                plots[2].add_data(
                    task=task_er,
                    xid="[Cve_dap]",
                    yid="KI__RTG",
                    label=f"{group}",
                    color=color,
                )
                plots[3].add_data(
                    task=task_er,
                    xid="[Cve_dap]",
                    yid="KI__GLCEX",
                    label=f"{group}",
                    color=color,
                )
                plots[4].add_data(
                    task=task_er,
                    xid="[Cve_dap]",
                    yid="KI__UGE",
                    label=f"{group}",
                    color=color,
                )

            figures[fig.sid] = fig
        return figures


if __name__ == "__main__":
    run_experiments(HepaticRenalImpairment, output_dir=HepaticRenalImpairment.__name__)