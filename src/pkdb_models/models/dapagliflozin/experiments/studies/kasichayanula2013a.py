from typing import Dict
from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console
from pkdb_models.models.dapagliflozin.experiments.base_experiment import DapagliflozinSimulationExperiment
from pkdb_models.models.dapagliflozin.experiments.metadata import Tissue, Route, Dosing, ApplicationForm, Health, \
    Fasting, DapagliflozinMappingMetaData
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim
from pkdb_models.models.dapagliflozin.helpers import run_experiments
from pathlib import Path
import pkdb_models.models.dapagliflozin as dapagliflozin


class Kasichayanula2013a(DapagliflozinSimulationExperiment):
    """Simulation experiment of Kasichayanula2013a."""

    groups = ["RIF", "MFA"]
    colors = {
        "RIF": "tab:blue",
        "MFA": "tab:orange",
    }
    info = [
        ("[Cve_dap]", "dapagliflozin"),
        ("Aurine_dap", "dapagliflozin_cumulative amount"),
        ("[Cve_d3g]", "dapagliflozin 3-o-glucuronide"),
        ("Aurine_d3g", "dapagliflozin 3-o-glucuronide_cumulative amount"),
    ]

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Tab1A"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                # unit conversion to mole/l
                if label.startswith("dapagliflozin_"):
                    dset.unit_conversion("mean", 1 / self.Mr.dap)
                elif label.startswith("dapagliflozin 3-o-glucuronide"):
                    dset.unit_conversion("mean", 1 / self.Mr.d3g)
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
                    "BW": Q_(self.bodyweight_default, "kg"),
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
        for kp, data in enumerate(self.info):
            sid, prefix = data[0], data[1]
            for group in self.groups:
                if "cumulative" in prefix:
                    tissue = Tissue.URINE
                else:
                    tissue = Tissue.PLASMA
                mappings[f"fm_po_dap10_{prefix}_{group}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{prefix}_DAP10_{group}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self, task=f"task_po_dap10", xid="time", yid=sid,
                    ),
                    metadata=DapagliflozinMappingMetaData(
                        tissue=tissue,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                    ),
                )
        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_rows=2,
            num_cols=2,
            name=f"{self.__class__.__name__}",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        plots[1].set_yaxis(self.label_dap_urine, unit=self.unit_dap_urine)
        plots[2].set_yaxis(self.label_d3g_plasma, unit=self.unit_d3g)
        plots[3].set_yaxis(self.label_d3g_urine, unit=self.unit_d3g_urine)
        for kp, data in enumerate(self.info):
            sid, prefix = data[0], data[1]
            # simulation
            plots[kp].add_data(
                task=f"task_po_dap10",
                xid="time",
                yid=sid,
                label="10 mg PO",
                color="black",
            )
            # data
            for group in self.groups:
                plots[kp].add_data(
                    dataset=f"{prefix}_DAP10_{group}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"{group} - 10 mg PO",
                    color=self.colors[group],
                    linestyle="" if sid in ["Aurine_dap", "Aurine_d3g"] else "--",
                )
        return {
            fig.sid: fig,
        }


if __name__ == "__main__":
    # run_experiments(Kasichayanula2013a, output_dir=Kasichayanula2013a.__name__)
    out = dapagliflozin.RESULTS_PATH_SIMULATION / Kasichayanula2013a.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Kasichayanula2013a, output_dir=out)