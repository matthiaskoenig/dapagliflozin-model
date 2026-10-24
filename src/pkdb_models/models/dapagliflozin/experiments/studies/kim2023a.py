"""Simulation experiment of the study Kim2023a."""

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


class Kim2023a(DapagliflozinSimulationExperiment):
    """Simulation experiment of Kim2023a.

    Single dose of 10 mg dapagliflozin propanediol monohydrate (Forxiga, DAP10)
    and of 10 mg dapagliflozin formate (DAP-FOR10, an ester prodrug converted
    during the absorption) in healthy Korean men (crossover, 71.49 kg). The
    fasting state at the dose is not reported (the model takes it as fasted).
    Whether the 10 mg of the formate is the ester or its dapagliflozin
    equivalent is not stated; the bioequivalence (GMR AUCinf 0.98) suggests
    the equivalent, both are simulated as 10 mg dapagliflozin. The sampling
    times are the nominal times; 10-48 hr are from the log panel (Fig2B).
    """

    interventions: ClassVar[list[str]] = ["DAP10", "DAP-FOR10"]
    bodyweight = 71.49  # [kg]
    fpg = 5.28  # (95.48/18 + 94.63/18)/2  # fasting plasma glucose [mM]

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                # unit conversion
                if label.startswith("dapagliflozin_"):
                    dset.unit_conversion("mean", 1 / self.Mr.dap)
                dsets[self.to_sid(label)] = dset
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        simulations["po_dap10"] = Simulation(
            time_unit="hr",
            end=50,
            preinit_changes={
                **self.default_changes(),
                # physiological changes
                "BW": Q(self.bodyweight, "kg"),
                "[KI__fpg]": Q(self.fpg, "mM"),
                "GU__fed": Q(self.fasting_map["not reported"], "dimensionless"),
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
            mappings[self.to_sid(f"fm_dap10_{intervention}")] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=self.to_sid(f"dapagliflozin_{intervention}"),
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
                    fasting=Fasting.NR,
                ),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2",
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        # simulation
        plots[0].add_data(
            task="task_po_dap10",
            xid="time",
            yid="[Cve_dap]",
            label="10 mg PO",
            color="black",
        )
        # data
        for intervention in self.interventions:
            label = "10 mg PO" if intervention == "DAP10" else "DAP-FOR 10 mg PO"
            plots[0].add_data(
                dataset=self.to_sid(f"dapagliflozin_{intervention}"),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=label,
                color="black" if intervention == "DAP10" else "tab:blue",
            )
        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Kim2023a, output_dir=Kim2023a.__name__)
