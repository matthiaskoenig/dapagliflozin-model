"""Dependency of the pharmacokinetics and pharmacodynamics on the dose and glucose.

Single oral doses of dapagliflozin (0 - 500 mg) at a plasma glucose of 5 mM, and
plasma glucose levels from 5 to 15 mM without dapagliflozin and with a single
oral dose of 10 mg.
"""

from typing import ClassVar, override

from sbmlsim import Q
from sbmlsim.plot import Axis, Figure, Plot
from sbmlsim.simulation import Simulation

from pkdb_models.models.dapagliflozin import MODEL
from pkdb_models.models.dapagliflozin.experiments.base_experiment import (
    DapagliflozinSimulationExperiment,
)


class DoseDependencyExperiment(DapagliflozinSimulationExperiment):
    """Test different doses of dapagliflozin and different glucose concentrations."""

    dapagliflozin_doses: ClassVar[list[float]] = [
        0,
        2.5,
        5,
        10,
        20,
        50,
        100,
        250,
        500,
    ]  # [mg]
    glucose_values: ClassVar[list[float]] = [5, 7, 9, 11, 13, 15]  # [mM]
    #: oral doses of dapagliflozin of the glucose simulations [mg]
    glucose_doses: ClassVar[list[float]] = [0, 10]
    colors_glc: ClassVar[list[str]] = [
        "black",
        "#d4b9da",
        "#c994c7",
        "#df65b0",
        "#dd1c77",
        "#980043",
    ]

    def _simulation(self, dose: float, glc: float) -> Simulation:
        """Single oral dose of dapagliflozin at a plasma glucose.

        Args:
            dose: oral dose of dapagliflozin [mg].
            glc: plasma glucose [mM].
        """
        return Simulation(
            time_unit="hr",
            end=96,
            preinit_changes={
                **self.default_changes(),
                "PODOSE_dap": Q(dose, "mg"),
                "[KI__fpg]": Q(glc, "mM"),  # Glucose
            },
        )

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations: dict[str, Simulation] = {}
        # dapagliflozin
        for dap in self.dapagliflozin_doses:
            simulations[self.to_sid(f"dapagliflozin_glc5_dap{dap}")] = self._simulation(
                dose=dap, glc=5
            )
        # glucose
        for glc in self.glucose_values:
            for dose in self.glucose_doses:
                simulations[self.to_sid(f"glucose_glc{glc}_dap{dose}")] = (
                    self._simulation(dose=dose, glc=glc)
                )

        return simulations

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self.figure_pk(),
            **self.figure_pd(),
        }

    def _add_curves(
        self, plot: Plot, key: str, xid: str, yid: str, by_dose: bool = True
    ) -> None:
        """Add the curves of the doses of dapagliflozin or of the glucose levels.

        Args:
            plot: plot to add the curves to.
            key: `dapagliflozin` (the doses of dapagliflozin at 5 mM glucose) or
                `glucose` (the glucose levels without and with 10 mg
                dapagliflozin, dashed without dapagliflozin).
            xid: id of the x data.
            yid: id of the y data.
            by_dose: order of the glucose curves, all glucose levels of a dose
                after another, otherwise the doses of a glucose level after
                another.
        """
        if key == "dapagliflozin":
            for dose in self.dapagliflozin_doses:
                plot.add_data(
                    task=self.to_sid(f"task_dapagliflozin_glc5_dap{dose}"),
                    xid=xid,
                    yid=yid,
                    label=f"{dose} mg",
                    color=self.dose_colors[dose],
                )
            return

        curves = [
            (k, dose)
            for dose in self.glucose_doses
            for k in range(len(self.glucose_values))
        ]
        if not by_dose:
            curves = [
                (k, dose)
                for k in range(len(self.glucose_values))
                for dose in self.glucose_doses
            ]
        for k, dose in curves:
            glc = self.glucose_values[k]
            plot.add_data(
                task=self.to_sid(f"task_glucose_glc{glc}_dap{dose}"),
                xid=xid,
                yid=yid,
                label=f"{glc} mM at {dose} mg",
                color=self.colors_glc[k],
                linestyle="-" if dose == 10 else "--",
            )

    def figure_pk(self) -> dict[str, Figure]:
        """Timecourses of dapagliflozin and its glucuronide in plasma, urine and feces."""
        figures: list[Figure] = []
        for key, name in [
            ("dapagliflozin", "Pharmacokinetics: Dapagliflozin dose dependency"),
            ("glucose", "Pharmacokinetics: Glucose dependency"),
        ]:
            fig = Figure(
                experiment=self,
                sid=self.to_sid(f"Fig_{key}_dependency_pk"),
                num_rows=1,
                num_cols=5,
                name=name,
                settings={"legend_scope": "row"},
            )
            plots = fig.create_plots(xaxis=Axis("time", unit="hr"), legend=True)
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
                self._add_curves(plots[ksid], key, xid="time", yid=sid)
            figures.append(fig)
        return self.figures_by_sid(*figures)

    def figure_pd(self) -> dict[str, Figure]:
        """Timecourses of RTG and UGE and their dependency on dapagliflozin in plasma."""
        figures: list[Figure] = []
        for key, name in [
            ("dapagliflozin", "Pharmacodynamics: Dapagliflozin dose dependency"),
            ("glucose", "Pharmacodynamics: Glucose dependency"),
        ]:
            fig = Figure(
                experiment=self,
                sid=self.to_sid(f"Fig_{key}_dependency_pd"),
                num_rows=1,
                num_cols=5,
                name=name,
                settings={"legend_scope": "row"},
            )
            plots = fig.create_plots(xaxis=Axis("time", unit="hr"), legend=True)
            sids = [
                "KI__RTG",
                "KI__UGE",
            ]
            for ksid, sid in enumerate(sids):
                plots[ksid].set_yaxis(label=self.labels[sid], unit=self.units[sid])

            for ksid, sid in enumerate(sids):
                self._add_curves(plots[ksid], key, xid="time", yid=sid)

            # exposure-response panels
            for ksid in range(2, 5):
                plots[ksid].set_xaxis(label=self.label_dap_plasma, unit=self.unit_dap)
            plots[2].set_yaxis(label=self.label_rtg, unit=self.unit_rtg)
            plots[3].set_yaxis(
                label=self.labels["KI__GLCEX"], unit=self.units["KI__GLCEX"]
            )
            plots[4].set_yaxis(label=self.label_uge, unit=self.unit_uge)

            for ksid, yid in [(2, "KI__RTG"), (3, "KI__GLCEX"), (4, "KI__UGE")]:
                self._add_curves(
                    plots[ksid], key, xid="[Cve_dap]", yid=yid, by_dose=False
                )
            figures.append(fig)
        return self.figures_by_sid(*figures)


if __name__ == "__main__":
    MODEL.run_experiments(
        DoseDependencyExperiment, output_dir=DoseDependencyExperiment.__name__
    )
