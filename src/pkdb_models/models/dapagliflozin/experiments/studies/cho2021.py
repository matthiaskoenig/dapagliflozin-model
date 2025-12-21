from typing import Dict
from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console
from pkdb_models.models.dapagliflozin.experiments.base_experiment import DapagliflozinSimulationExperiment
from pkdb_models.models.dapagliflozin.experiments.metadata import Tissue, Route, Dosing, ApplicationForm, Health, \
    Fasting, DapagliflozinMappingMetaData
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim
import pkdb_models.models.dapagliflozin as dapagliflozin
from pathlib import Path
from pkdb_models.models.dapagliflozin.helpers import run_experiments


class Cho2021(DapagliflozinSimulationExperiment):
    """Simulation experiment of Cho2021."""

    interventions = ["DAP10", "RC10"] # reference and test formulation
    bodyweight = 72.8  # [kg]

    colors = ["black", "tab:blue"]

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig6"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                # unit conversion
                if label.startswith("dapagliflozin_"):
                    dset.unit_conversion("mean", 1 / self.Mr.dap)
                dsets[f"{label}"] = dset
        # console.print(dsets)
        # console.print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        tcsims[f"po_dap10"] = TimecourseSim(
            [Timecourse(
                start=0,
                end=25 * 60,  # [min]
                steps=500,
                changes={
                    **self.default_changes(),
                    # physiological changes
                    "BW": Q_(self.bodyweight, "kg"),
                    "[KI__glc_ext]": Q_(self.fpg_healthy, "mM"),
                    "GU__f_absorption": Q_(self.fasting_map["fasted"], "dimensionless"),
                    "f_cirrhosis": Q_(self.cirrhosis_map["Control"], "dimensionless"),
                    "KI__f_renal_function": Q_(self.renal_map["Normal renal function"], "dimensionless"),
                    # dose (IVDOSE, PODOSE)
                    "PODOSE_dap": Q_(10, "mg"),
                },
            )]
        )

        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
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
                    self, task=f"task_po_dap10", xid="time", yid=f"[Cve_dap]",
                ),
                metadata=DapagliflozinMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
                ),
            )
        # console.print(mappings)
        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig6",
            num_rows=1,
            num_cols=1,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        # simulation
        plots[0].add_data(
            task=f"task_po_dap10",
            xid="time",
            yid=f"[Cve_dap]",
            label="10 mg PO",
            color="black",
        )
        # data
        for ks, intervention in enumerate(self.interventions):
            if intervention == "DAP10":
                label = "Ref 10 mg PO"
            else:
                label = "Test 10 mg PO"
            plots[0].add_data(
                dataset=f"dapagliflozin_{intervention}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=label,
                color=self.colors[ks],
            )
        return {
            fig.sid: fig,
        }


if __name__ == "__main__":
    out = dapagliflozin.RESULTS_PATH_SIMULATION / Cho2021.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Cho2021, output_dir=out)