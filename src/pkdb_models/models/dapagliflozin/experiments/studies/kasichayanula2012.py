"""Simulation experiment of the study Kasichayanula2012."""

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


class Kasichayanula2012(DapagliflozinSimulationExperiment):
    """Simulation experiment of Kasichayanula2012.

    Single dose of 20 mg dapagliflozin (tablet, fasted) alone and with
    simvastatin (40 mg) or valsartan (320 mg) in healthy subjects (study 1, a
    crossover of the same 24 subjects; body weight not reported). Neither drug
    changes the exposure of dapagliflozin (AUC ratios 0.98 with simvastatin,
    1.02 with valsartan, Cmax ratios 0.98 and 0.88), so all arms are simulated
    without an interaction factor; the coadministration arms are validation
    data of the interactions. Study 2 (warfarin, digoxin) has no data of
    dapagliflozin and is not simulated.
    """

    interventions: ClassVar[list[str]] = ["DAP20", "DAP20, SIM40", "DAP20, VAL320"]
    colors: ClassVar[dict[str, str]] = {
        "DAP20": "black",
        "DAP20, SIM40": "tab:blue",
        "DAP20, VAL320": "tab:orange",
    }

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
        simulations["po_dap20"] = Simulation(
            time_unit="hr",
            end=25,
            preinit_changes={
                **self.default_changes(),
                # physiological changes
                "BW": Q(self.bodyweight_default, "kg"),
                "[KI__fpg]": Q(self.fpg_healthy, "mM"),
                "GU__fed": Q(self.fasting_map["fasted"], "dimensionless"),
                "f_cirrhosis": Q(self.cirrhosis_map["Control"], "dimensionless"),
                "KI__f_renal_function": Q(
                    self.renal_map["Normal renal function"], "dimensionless"
                ),
                # dose (IVDOSE, PODOSE)
                "PODOSE_dap": Q(20, "mg"),
            },
        )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for intervention in self.interventions:
            mappings[self.to_sid(f"fm_dap20_{intervention}")] = FitMapping(
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
                    task="task_po_dap20",
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
                    coadministration=Coadministration.VALSARTAN
                    if "VAL" in intervention
                    else Coadministration.SIMVASTATIN
                    if "SIM" in intervention
                    else Coadministration.NONE,
                ),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2",
            num_rows=1,
            num_cols=1,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        # simulation, the same for all arms (no interaction)
        plots[0].add_data(
            task="task_po_dap20",
            xid="time",
            yid="[Cve_dap]",
            label="20 mg PO",
            color="black",
        )
        # data
        plots[0].add_data(
            dataset="dapagliflozin_DAP20",
            xid="time",
            yid="mean",
            yid_sd="mean_sd",
            count="count",
            label="20 mg PO",
            color=self.colors["DAP20"],
        )
        plots[0].add_data(
            dataset="dapagliflozin_DAP20__SIM40",
            xid="time",
            yid="mean",
            yid_sd="mean_sd",
            count="count",
            label="20 mg PO + SIM40",
            color=self.colors["DAP20, SIM40"],
        )
        plots[0].add_data(
            dataset="dapagliflozin_DAP20__VAL320",
            xid="time",
            yid="mean",
            yid_sd="mean_sd",
            count="count",
            label="20 mg PO + VAL320",
            color=self.colors["DAP20, VAL320"],
        )
        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Kasichayanula2012, output_dir=Kasichayanula2012.__name__)
