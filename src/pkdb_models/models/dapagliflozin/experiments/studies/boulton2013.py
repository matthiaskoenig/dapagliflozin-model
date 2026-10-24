"""Simulation experiment of the study Boulton2013."""

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
from sbmlsim.simulation import Change, Simulation

from pkdb_models.models.dapagliflozin import MODEL
from pkdb_models.models.dapagliflozin.experiments.base_experiment import (
    DapagliflozinSimulationExperiment,
)


class Boulton2013(DapagliflozinSimulationExperiment):
    """Simulation experiment of Boulton2013.

    Absolute bioavailability: 7 healthy men (body weight not reported) took
    10 mg dapagliflozin orally after an overnight fast of at least 10 hr and
    received about 80 µg [14C]-dapagliflozin as an intravenous infusion over
    1 min at 1 hr (tmax of the oral dose). The unlabeled drug of the oral dose
    is measured by LC-MS/MS (LLOQ 1 ng/ml), the intravenous dose by AMS (LLOQ
    9 pg/ml); the AMS concentrations are expressed for the whole intravenous
    dose (80 µg / AUCinf 6.8 ng*h/ml = 196 ml/min, reported CL 207 ml/min), so
    the full 80 µg are simulated. The intravenous dose adds 0.8 % to the oral
    curve. The intravenous tail of 36-49 hr is flat at about twice the LLOQ of
    AMS (assay limited). Reported: F 78 %, CL 207 ml/min, Vss 118 l, t1/2
    12.2 hr (i.v.) and 13.7 hr (oral).
    """

    simulation_keys: ClassVar[list[str]] = ["po_iv_dap10", "iv_dap008"]
    interventions: ClassVar[list[str]] = ["DAP10", "C14DAP008"]
    colors: ClassVar[list[str]] = ["black", "tab:blue"]

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1"]:
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
        # po followed by iv
        simulations["po_iv_dap10"] = Simulation(
            time_unit="hr",
            end=50,
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
                "PODOSE_dap": Q(10, "mg"),
            },
            changes=[
                Change(
                    1,
                    {
                        "IVDOSE_dap": Q(80, "µg"),  # 1.6 % labeled; 98.4% unlabeled
                    },
                ),
            ],
        )
        # only iv at 1 hr
        simulations["iv_dap008"] = Simulation(
            time_unit="hr",
            end=50,
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
            },
            changes=[
                Change(
                    1,
                    {
                        "IVDOSE_dap": Q(80, "µg"),  # 1.6 % labeled; 98.4% unlabeled
                    },
                ),
            ],
        )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for ks, sim_key in enumerate(self.simulation_keys):
            intervention = self.interventions[ks]
            mappings[f"fm_{intervention}_dapagliflozin"] = FitMapping(
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
                    task=f"task_{sim_key}",
                    xid="time",
                    yid="[Cve_dap]",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO if intervention == "DAP10" else Route.IV,
                    # the oral curve is the unlabeled tablet dose (LC-MS/MS)
                    application_form=ApplicationForm.TABLET
                    if intervention == "DAP10"
                    else ApplicationForm.SOLUTION,
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
            sid="Fig1",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        # log scale as the source, the data spans 0.02-143 ng/ml
        plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap, scale="log")
        plots[1].set_yaxis(self.label_dap_plasma, unit=self.unit_dap, scale="log")

        for ks, sim_key in enumerate(self.simulation_keys):
            dose_label = "10 mg PO" if sim_key == "po_iv_dap10" else "80 µg IV"
            # simulation
            plots[ks].add_data(
                task=f"task_{sim_key}",
                xid="time",
                yid="[Cve_dap]",
                label=dose_label,
                color=self.colors[ks],
            )
            # data
            plots[ks].add_data(
                dataset=f"dapagliflozin_{self.interventions[ks]}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=dose_label,
                color=self.colors[ks],
            )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Boulton2013, output_dir=Boulton2013.__name__)
