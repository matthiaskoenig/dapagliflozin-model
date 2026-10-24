"""Simulation experiment of the study Shah2019a."""

from typing import ClassVar, override

from pkpdsim.metadata import (
    ApplicationForm,
    Coadministration,
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


class Shah2019a(DapagliflozinSimulationExperiment):
    """Simulation experiment of Shah2019a.

    Single dose of the fixed-dose combination of 5 mg dapagliflozin and 500 mg
    metformin XR in the fasted and fed state (high-fat breakfast 30 min before
    the dose, crossover) in 24 healthy Indian men. Metformin does not change the
    exposure of dapagliflozin (Kasichayanula2011c), no interaction factor. The
    error bars of Fig 3 are about constant over the profile, also at the predose
    0, and are not SDs of the means (Cmax SD 13.4 ng/ml in Tab 7, bar 5.7 at the
    peak): the data have no SD.
    """

    conditions: ClassVar[list[str]] = ["fasted", "fed"]
    #: body weight [kg], the midpoint of the range 54.6-82.2 kg, no mean
    #: reported
    bodyweight = (54.6 + 82.2) / 2

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig3"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                # unit conversion to mole/l
                if label.startswith("dapagliflozin_"):
                    dset.unit_conversion("mean", 1 / self.Mr.dap)
                dsets[self.to_sid(label)] = dset
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for condition in self.conditions:
            simulations[self.to_sid(f"po_dap5_{condition}")] = Simulation(
                time_unit="hr",
                end=100,
                preinit_changes={
                    **self.default_changes(),
                    # physiological changes
                    "BW": Q(self.bodyweight, "kg"),
                    "[KI__fpg]": Q(self.fpg_healthy, "mM"),
                    "GU__fed": Q(self.fasting_map[condition], "dimensionless"),
                    "f_cirrhosis": Q(self.cirrhosis_map["Control"], "dimensionless"),
                    "KI__f_renal_function": Q(
                        self.renal_map["Normal renal function"], "dimensionless"
                    ),
                    # dose (IVDOSE, PODOSE)
                    "PODOSE_dap": Q(5, "mg"),
                },
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for condition in self.conditions:
            mappings[self.to_sid(f"fm_po_dap5_{condition}")] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=self.to_sid(f"dapagliflozin_DAP5, MET500_{condition}"),
                    xid="time",
                    yid="mean",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task=self.to_sid(f"task_po_dap5_{condition}"),
                    xid="time",
                    yid="[Cve_dap]",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    coadministration=Coadministration.METFORMIN,
                    fasting=Fasting.FASTED if condition == "fasted" else Fasting.FED,
                ),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig3",
            num_rows=1,
            num_cols=1,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        for condition in self.conditions:
            # simulation
            plots[0].add_data(
                task=self.to_sid(f"task_po_dap5_{condition}"),
                xid="time",
                yid="[Cve_dap]",
                label=f"{condition} - 5 mg PO",
                color=self.fasting_colors[condition],
            )
            # data
            plots[0].add_data(
                dataset=self.to_sid(f"dapagliflozin_DAP5, MET500_{condition}"),
                xid="time",
                yid="mean",
                count="count",
                label=f"{condition} - 5 mg PO + MET500",
                color=self.fasting_colors[condition],
            )
        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Shah2019a, output_dir=Shah2019a.__name__)
