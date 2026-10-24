"""Hepatic and renal impairment, the groups of the studies of dapagliflozin.

The degrees of cirrhosis (`f_cirrhosis`) and of renal impairment
(`KI__f_renal_function`) of the groups of the studies are simulated after a
single oral dose of dapagliflozin.
"""

from collections.abc import Mapping
from typing import ClassVar, override

from sbmlsim import Q
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Simulation

from pkdb_models.models.dapagliflozin import MODEL
from pkdb_models.models.dapagliflozin.experiments.base_experiment import (
    DapagliflozinSimulationExperiment,
)


class HepaticRenalImpairment(DapagliflozinSimulationExperiment):
    """Tests hepatic and renal impairment.

    Hepatic impairment after a single oral dose of 10 mg; renal impairment
    after a single oral dose of 50 mg (used for the PK plots) and of 20 mg (used
    for the PD plots).
    """

    maps: ClassVar[dict[str, Mapping[str, float]]] = {
        "hepatic": DapagliflozinSimulationExperiment.cirrhosis_map,
        "renal": DapagliflozinSimulationExperiment.renal_map,
    }
    parameters: ClassVar[dict[str, str]] = {
        "hepatic": "f_cirrhosis",
        "renal": "KI__f_renal_function",
    }
    colors: ClassVar[dict[str, Mapping[str, str]]] = {
        "hepatic": DapagliflozinSimulationExperiment.cirrhosis_colors,
        "renal": DapagliflozinSimulationExperiment.renal_colors,
    }
    impairments: ClassVar[list[str]] = list(maps)
    #: oral doses of the figures of the pharmacokinetics and pharmacodynamics [mg]
    doses_pk: ClassVar[dict[str, float]] = {"hepatic": 10, "renal": 50}
    doses_pd: ClassVar[dict[str, float]] = {"hepatic": 10, "renal": 20}

    def _simulation(self, dose: float, parameter: str, value: float) -> Simulation:
        """Single oral dose of dapagliflozin with an impairment.

        Args:
            dose: oral dose of dapagliflozin [mg].
            parameter: parameter of the impairment.
            value: value of the parameter.
        """
        return Simulation(
            time_unit="hr",
            end=96,
            preinit_changes={
                **self.default_changes(),
                "PODOSE_dap": Q(dose, "mg"),
                parameter: Q(value, "dimensionless"),
            },
        )

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations: dict[str, Simulation] = {}
        for impairment in self.impairments:
            parameter = self.parameters[impairment]

            for group, value in self.maps[impairment].items():
                if impairment == "hepatic":
                    # 10 mg PO
                    simulations[self.to_sid(f"dap_hepatic_{group}")] = self._simulation(
                        dose=10, parameter=parameter, value=value
                    )
                else:
                    # 50 mg PO (for PK plots)
                    simulations[self.to_sid(f"dap_renal_{group}_po50")] = (
                        self._simulation(dose=50, parameter=parameter, value=value)
                    )
                    # 20 mg PO (for PD plots)
                    simulations[self.to_sid(f"dap_renal_{group}_po20")] = (
                        self._simulation(dose=20, parameter=parameter, value=value)
                    )
        return simulations

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self.figure_pk(),
            **self.figure_pd(),
        }

    def figure_pk(self) -> dict[str, Figure]:
        """Timecourses of dapagliflozin and its glucuronide in plasma, urine and feces."""
        figures: list[Figure] = []

        for impairment in self.impairments:
            parameter = self.parameters[impairment]
            colors = self.colors[impairment]

            fig = Figure(
                experiment=self,
                sid=self.to_sid(f"Fig_{parameter}_pk"),
                num_rows=1,
                num_cols=5,
                name=(
                    f"Pharmacokinetics: {impairment.title()} impairment, "
                    f"{self.doses_pk[impairment]} mg PO"
                ),
            )
            plots = fig.create_plots(
                xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
            )

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
                plots[ksid].set_yaxis(label=self.labels[sid], unit=self.units[sid])

                for group in self.maps[impairment]:
                    if impairment == "hepatic":
                        task_id = f"task_dap_hepatic_{group}"  # 10 mg
                    else:
                        task_id = f"task_dap_renal_{group}_po50"  # 50 mg for PK

                    plots[ksid].add_data(
                        task=self.to_sid(task_id),
                        xid="time",
                        yid=sid,
                        label=group,
                        color=colors[group],
                    )

            figures.append(fig)
        return self.figures_by_sid(*figures)

    def figure_pd(self) -> dict[str, Figure]:
        """Timecourses of RTG and UGE and their dependency on dapagliflozin in plasma."""
        figures: list[Figure] = []

        for impairment in self.impairments:
            parameter = self.parameters[impairment]
            colors = self.colors[impairment]

            fig = Figure(
                experiment=self,
                sid=self.to_sid(f"Fig_{parameter}_pd"),
                num_rows=1,
                num_cols=5,
                name=(
                    f"Pharmacodynamics: {impairment.title()} impairment, "
                    f"{self.doses_pd[impairment]} mg PO"
                ),
            )
            plots = fig.create_plots(
                xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
            )

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
            plots[3].set_yaxis(
                label=self.labels["KI__GLCEX"], unit=self.units["KI__GLCEX"]
            )
            plots[4].set_yaxis(label=self.label_uge, unit=self.unit_uge)

            for group in self.maps[impairment]:
                if impairment == "hepatic":
                    task = f"task_dap_hepatic_{group}"  # 10 mg
                else:
                    task = f"task_dap_renal_{group}_po20"  # 20 mg

                for ksid, xid, yid in [
                    (0, "time", "KI__RTG"),
                    (1, "time", "KI__UGE"),
                    (2, "[Cve_dap]", "KI__RTG"),
                    (3, "[Cve_dap]", "KI__GLCEX"),
                    (4, "[Cve_dap]", "KI__UGE"),
                ]:
                    plots[ksid].add_data(
                        task=self.to_sid(task),
                        xid=xid,
                        yid=yid,
                        label=group,
                        color=colors[group],
                    )

            figures.append(fig)
        return self.figures_by_sid(*figures)


if __name__ == "__main__":
    MODEL.run_experiments(
        HepaticRenalImpairment, output_dir=HepaticRenalImpairment.__name__
    )
