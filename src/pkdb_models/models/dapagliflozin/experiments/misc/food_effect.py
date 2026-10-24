"""Effect of food on dapagliflozin, as a scan of the absorption.

The effect of food is modelled as a change of the absorption rate of
dapagliflozin in the intestine in the fed state (`GU__fed` 1, the factor
`GU__FED_f_absorption` is scanned) after a single oral dose of 10 mg
dapagliflozin.
"""

from typing import ClassVar, Literal, override

import matplotlib
import numpy as np
from sbmlsim import Q
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Simulation

from pkdb_models.models.dapagliflozin import MODEL
from pkdb_models.models.dapagliflozin.experiments.base_experiment import (
    DapagliflozinSimulationExperiment,
)


class FoodEffect(DapagliflozinSimulationExperiment):
    """Tests food intake effects."""

    ddi_parameters: ClassVar[list[str]] = ["GU__FED_f_absorption"]
    parameter_values: ClassVar[list[float]] = sorted(
        [*np.logspace(start=-1, stop=1, num=10).tolist(), 1.0]
    )  # [0.1, 10]

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations: dict[str, Simulation] = {}
        for parameter in self.ddi_parameters:
            for value in self.parameter_values:
                simulations[self.to_sid(f"dap_{parameter}_{value}")] = Simulation(
                    time_unit="hr",
                    end=96,
                    preinit_changes={
                        **self.default_changes(),
                        "PODOSE_dap": Q(10, "mg"),
                        "GU__fed": Q(1, "dimensionless"),
                        parameter: Q(value, "dimensionless"),
                    },
                )
        return simulations

    @staticmethod
    def get_color_for_absorption(
        v: float,
        scale: Literal["linear", "log"] = "log",
        vmin: float = 0.1,
        vmax: float = 10,
    ) -> str:
        """Get the color of an absorption activity, black for the reference 1.0.

        The color is a hex string: the plots over the draws of the fitted
        parameters (`--uncertainty`) take no tuple of numpy floats.

        Args:
            v: absorption activity.
            scale: scale of the colors between `vmin` and `vmax`.
            vmin: absorption activity of the first color.
            vmax: absorption activity of the last color.

        Raises:
            ValueError: for an unsupported scale.
        """
        cmap = matplotlib.colormaps["seismic_r"]
        if np.isclose(v, 1.0):
            return "black"
        if scale == "linear":
            cvalue = (v - vmin) / np.abs(vmax - vmin)
        elif scale == "log":
            cvalue = (np.log10(v) - np.log10(vmin)) / np.abs(
                np.log10(vmax) - np.log10(vmin)
            )
        else:
            raise ValueError(f"Unsupported scale: '{scale}'")
        return matplotlib.colors.to_hex(cmap(cvalue))

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self.figures_pk(),
            **self.figures_pd(),
        }

    def figures_pk(self) -> dict[str, Figure]:
        """Timecourses of dapagliflozin and its glucuronide in plasma, urine and feces."""
        figures: list[Figure] = []
        for parameter in self.ddi_parameters:
            fig = Figure(
                experiment=self,
                sid=self.to_sid(f"Fig_{parameter}_pk"),
                num_rows=1,
                num_cols=5,
                name="Pharmacokinetics: Food effect",
            )
            plots = fig.create_plots(
                xaxis=Axis(self.label_time, self.unit_time), legend=True
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
                for value in self.parameter_values:
                    plots[ksid].add_data(
                        task=self.to_sid(f"task_dap_{parameter}_{value}"),
                        xid="time",
                        yid=sid,
                        label=f"{value:.1f}",
                        color=self.get_color_for_absorption(v=value),
                    )
            figures.append(fig)
        return self.figures_by_sid(*figures)

    def figures_pd(self) -> dict[str, Figure]:
        """Timecourses of RTG and UGE and their dependency on dapagliflozin in plasma."""
        figures: list[Figure] = []
        for parameter in self.ddi_parameters:
            fig = Figure(
                experiment=self,
                sid=self.to_sid(f"Fig_{parameter}_pd"),
                num_rows=1,
                num_cols=5,
                name="Pharmacodynamics: Food effect",
            )
            plots = fig.create_plots(
                xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
            )
            sids = [
                "KI__RTG",
                "KI__UGE",
            ]
            for ksid, sid in enumerate(sids):
                plots[ksid].set_yaxis(label=self.labels[sid], unit=self.units[sid])
            for ksid, sid in enumerate(sids):
                for value in self.parameter_values:
                    plots[ksid].add_data(
                        task=self.to_sid(f"task_dap_{parameter}_{value}"),
                        xid="time",
                        yid=sid,
                        label=f"{value:.1f}",
                        color=self.get_color_for_absorption(v=value),
                    )

            # exposure-response panels
            for ksid in range(2, 5):
                plots[ksid].set_xaxis(label=self.label_dap_plasma, unit=self.unit_dap)
            plots[2].set_yaxis(label=self.label_rtg, unit=self.unit_rtg)
            plots[3].set_yaxis(
                label=self.labels["KI__GLCEX"], unit=self.units["KI__GLCEX"]
            )
            plots[4].set_yaxis(label=self.label_uge, unit=self.unit_uge)
            for value in self.parameter_values:
                color = self.get_color_for_absorption(v=value)
                for ksid, yid in [(2, "KI__RTG"), (3, "KI__GLCEX"), (4, "KI__UGE")]:
                    plots[ksid].add_data(
                        task=self.to_sid(f"task_dap_{parameter}_{value}"),
                        xid="[Cve_dap]",
                        yid=yid,
                        label=f"{value:.1f}",
                        color=color,
                    )
            figures.append(fig)
        return self.figures_by_sid(*figures)


if __name__ == "__main__":
    MODEL.run_experiments(FoodEffect, output_dir=FoodEffect.__name__)
