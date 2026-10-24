"""Simulation experiment of the study Imamura2013."""

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


class Imamura2013(DapagliflozinSimulationExperiment):
    """Simulation experiment of Imamura2013.

    22 Japanese patients with type 2 diabetes (median 67.3 kg) on voglibose
    0.2 mg t.i.d., a single 10 mg tablet of dapagliflozin with voglibose, then
    after 7 days without both drugs alone; fasting not reported; plasma 0-72 hr
    (Fig2 0-24 hr, mean and SD). Voglibose is not absorbed and does not change
    the exposure (AUC ratio 1.009, Cmax 1.040), both arms are one simulation.
    No baseline glucose is reported (inclusion FPG <= 180 mg/dl, HbA1c <=
    8.5 %), the glucose is the default of type 2 diabetes; there are no PD data.
    """

    interventions: ClassVar[list[str]] = ["DAP10", "DAP10, VOG"]
    bodyweight = 67.3  # [kg]

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
            end=25,
            preinit_changes={
                **self.default_changes(),
                # physiological changes
                "BW": Q(self.bodyweight, "kg"),
                "[KI__fpg]": Q(self.fpg_t2dm, "mM"),
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
                    health=Health.T2DM,
                    fasting=Fasting.NR,
                    coadministration=Coadministration.VOGLIBOSE
                    if "VOG" in intervention
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
            name=f"{self.__class__.__name__} (T2DM)",
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
        plots[0].add_data(
            dataset="dapagliflozin_DAP10",
            xid="time",
            yid="mean",
            yid_sd="mean_sd",
            count="count",
            label="10 mg DAP",
            color="black",
        )
        plots[0].add_data(
            dataset="dapagliflozin_DAP10__VOG",
            xid="time",
            yid="mean",
            yid_sd="mean_sd",
            count="count",
            label="10 mg DAP + VOG",
            color="tab:blue",
        )
        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Imamura2013, output_dir=Imamura2013.__name__)
