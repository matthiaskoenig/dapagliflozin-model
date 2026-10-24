"""Simulation experiment of the study Cho2021."""

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


class Cho2021(DapagliflozinSimulationExperiment):
    """Simulation experiment of Cho2021.

    Bioequivalence of a cocrystal tablet (RC-T4, test, `RC10`) and Forxiga
    10 mg (reference, `DAP10`) in a randomized two-way crossover: 37 healthy
    Korean men (72.8 kg) took the tablet after a fast of at least 10 hr. The
    caption gives mean +- SE, the bars are the SD (52 ng/ml at 1 hr, SD of Cmax
    54.8 ng/ml in Table 5) and are curated as SD. Both curves are the same
    subjects: the test formulation is bioequivalent (AUC ratio 1.04, Cmax
    1.01), its mapping counts the subjects a second time (aggregate).
    """

    interventions: ClassVar[list[str]] = ["DAP10", "RC10"]  # reference and test
    bodyweight = 72.8  # [kg]

    colors: ClassVar[list[str]] = ["black", "tab:blue"]

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig6"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                # unit conversion
                if label.startswith("dapagliflozin_"):
                    dset.unit_conversion("mean", 1 / self.Mr.dap)
                dsets[label] = dset
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        simulations["po_dap10"] = Simulation(
            time_unit="hr",
            end=25,
            preinit_changes={
                **self.default_changes(),
                # physiological changes
                "BW": Q(self.bodyweight, "kg"),
                "[KI__fpg]": Q(self.fpg_healthy, "mM"),
                "GU__fed": Q(self.fasting_map["fasted"], "dimensionless"),
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
        for intervention in self.interventions:
            mappings[f"fm_dap10_{intervention}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"dapagliflozin_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task="task_po_dap10",
                    xid="time",
                    yid="[Cve_dap]",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
                ),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig6",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        # linear and log scale as the source
        plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        plots[1].set_yaxis(self.label_dap_plasma, unit=self.unit_dap, scale="log")
        for plot in plots:
            # simulation
            plot.add_data(
                task="task_po_dap10",
                xid="time",
                yid="[Cve_dap]",
                label="10 mg PO",
                color="black",
            )
            # data
            for ks, intervention in enumerate(self.interventions):
                label = "Ref 10 mg PO" if intervention == "DAP10" else "Test 10 mg PO"
                plot.add_data(
                    dataset=f"dapagliflozin_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=label,
                    color=self.colors[ks],
                )
        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Cho2021, output_dir=Cho2021.__name__)
