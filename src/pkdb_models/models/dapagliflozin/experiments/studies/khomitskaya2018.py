"""Simulation experiment of the study Khomitskaya2018."""

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


class Khomitskaya2018(DapagliflozinSimulationExperiment):
    """Simulation experiment of Khomitskaya2018.

    Bioequivalence of a single dose of dapagliflozin 10 mg with metformin XR
    1000 mg in healthy subjects (crossover, 71.9 kg), 30 min after a standard
    breakfast (about 600 kcal): the fixed combination (FCDP, Xigduo XR
    10/1000 mg, `DAP10, MET1000`) and the individual tablets (Forxiga 10 mg
    with 2 x 500 mg metformin XR, `DAP10`). Both arms are coadministered with
    metformin, which does not change the exposure of dapagliflozin
    (Kasichayanula2011c), so both are simulated without an interaction
    factor; the late points (24-72 hr) are from the log inset of Fig2A.
    """

    interventions: ClassVar[list[str]] = ["DAP10", "DAP10, MET1000"]
    labels_intervention: ClassVar[dict[str, str]] = {
        "DAP10": "10 mg + MET XR 1000 mg",
        "DAP10, MET1000": "FCDP 10/1000 mg",
    }
    bodyweight = 71.9  # [kg]

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
            end=80,
            preinit_changes={
                **self.default_changes(),
                # physiological changes
                "BW": Q(self.bodyweight, "kg"),
                "[KI__fpg]": Q(self.fpg_healthy, "mM"),
                "GU__fed": Q(self.fasting_map["fed"], "dimensionless"),
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
                    fasting=Fasting.FED,
                    # both arms with metformin XR 1000 mg
                    coadministration=Coadministration.METFORMIN,
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
            label="10 mg PO (fed)",
            color="black",
        )
        # data
        for intervention in self.interventions:
            label = self.labels_intervention[intervention]
            plots[0].add_data(
                dataset=self.to_sid(f"dapagliflozin_{intervention}"),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=label,
                color="tab:blue" if "MET" in intervention else "black",
            )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Khomitskaya2018, output_dir=Khomitskaya2018.__name__)
