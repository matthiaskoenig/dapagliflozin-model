"""Simulation experiment of the study Kasichayanula2011b."""

from typing import ClassVar, override

from pkpdsim.metadata import (
    ApplicationForm,
    Dosing,
    Fasting,
    Health,
    PKPDMappingMetaData,
    Route,
    Tissue,
)
from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.fit import FitData, FitMapping
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Simulation

from pkdb_models.models.dapagliflozin import MODEL
from pkdb_models.models.dapagliflozin.experiments.base_experiment import (
    DapagliflozinSimulationExperiment,
)


class Kasichayanula2011b(DapagliflozinSimulationExperiment):
    """Simulation experiment of Kasichayanula2011b.

    Food effect: 14 healthy subjects (body weight not reported), a single 10 mg
    tablet after a fast of 10 hr or 5 min after a high-fat meal (crossover);
    plasma dapagliflozin and D3G 0-72 hr (mean and SD). The meal lowers Cmax to
    0.69 and delays tmax by 1 hr at an unchanged AUC (1.02); the D3G AUC falls to
    0.86, which a slower absorption alone does not give. The crossover is the
    data of the food effect of the model (`GU__fed`, `GU__FED_f_absorption`).
    The fasted dapagliflozin at 72 hr (below 1 ng/ml in the figure, hidden by
    the fed marker) is not in the data.
    """

    conditions: ClassVar[list[str]] = ["fasted", "fed"]

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                # unit conversion to mole/l
                if label.startswith("dapagliflozin_"):
                    dset.unit_conversion("mean", 1 / self.Mr.dap)
                elif label.startswith("dapagliflozin 3-o-glucuronide"):
                    dset.unit_conversion("mean", 1 / self.Mr.d3g)
                dsets[self.to_sid(label)] = dset
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for condition in self.conditions:
            simulations[self.to_sid(f"po_dap10_{condition}")] = Simulation(
                time_unit="hr",
                end=80,
                preinit_changes={
                    **self.default_changes(),
                    # physiological changes
                    "BW": Q(self.bodyweight_default, "kg"),
                    "[KI__fpg]": Q(self.fpg_healthy, "mM"),
                    "GU__fed": Q(self.fasting_map[condition], "dimensionless"),
                    "f_cirrhosis": Q(self.cirrhosis_map["Control"], "dimensionless"),
                    "KI__f_renal_function": Q(
                        self.renal_map["Normal renal function"], "dimensionless"
                    ),
                    # dose (IVDOSE, PODOSE)
                    "PODOSE_dap": Q(10, "mg"),
                },
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for condition in self.conditions:
            for sid, name in [
                ("dap", "dapagliflozin"),
                ("d3g", "dapagliflozin 3-o-glucuronide"),
            ]:
                mappings[self.to_sid(f"fm_po_dap10_{condition}_{sid}")] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=self.to_sid(f"{name}_DAP10_{condition}"),
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=self.to_sid(f"task_po_dap10_{condition}"),
                        xid="time",
                        yid=f"[Cve_{sid}]",
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED
                        if condition == "fasted"
                        else Fasting.FED,
                    ),
                )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        plots[1].set_yaxis(self.label_d3g_plasma, unit=self.unit_d3g)

        for condition in self.conditions:
            dose_label = f"{condition} - 10 mg PO"
            # simulation
            plots[0].add_data(
                task=self.to_sid(f"task_po_dap10_{condition}"),
                xid="time",
                yid="[Cve_dap]",
                label=dose_label,
                color=self.fasting_colors[condition],
            )
            plots[1].add_data(
                task=self.to_sid(f"task_po_dap10_{condition}"),
                xid="time",
                yid="[Cve_d3g]",
                label=dose_label,
                color=self.fasting_colors[condition],
            )
            # data
            plots[0].add_data(
                dataset=self.to_sid(f"dapagliflozin_DAP10_{condition}"),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=dose_label,
                color=self.fasting_colors[condition],
            )
            plots[1].add_data(
                dataset=self.to_sid(f"dapagliflozin 3-o-glucuronide_DAP10_{condition}"),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=dose_label,
                color=self.fasting_colors[condition],
            )
        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Kasichayanula2011b, output_dir=Kasichayanula2011b.__name__)
