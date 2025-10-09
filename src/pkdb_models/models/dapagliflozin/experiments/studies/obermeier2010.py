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


class Obermeier2010(DapagliflozinSimulationExperiment):
    """Simulation experiment of Obermeier2010."""
    # FIXME find supplemental for bodyweight

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig6"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                # unit conversion
                if label.startswith("dapagliflozin_"):
                    dset.unit_conversion("mean", 1 / self.Mr.dap)
                elif label.startswith("total dapagliflozin_"):
                    dset.unit_conversion("mean", 1 / self.Mr.daptot)
                dsets[f"{label}"] = dset
        # console.print(dsets)
        # console.print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        tcsims[f"po_dap50"] = TimecourseSim(
            [Timecourse(
                start=0,
                end=30 * 60,  # [min]
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
                    "PODOSE_dap": Q_(50, "mg"),
                },
            )]
        )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        infos = [
            ("[Cve_dap]", "dapagliflozin"),
            ("[Cve_daptot]", "total dapagliflozin"),
        ]
        for sid, name in infos:
            mappings[f"fm_dap50_DAP50"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"{name}_DAP50",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self, task=f"task_po_dap50", xid="time", yid=sid,
                ),
                metadata=DapagliflozinMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.SOLUTION,
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
            name=f"{self.__class__.__name__}",
        )
        Figure.legend_fontsize = 11
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis("Plasma concentration", unit=self.unit_dap)
        substance = [
            ("[Cve_dap]", "dap", "black"),
            ("[Cve_daptot]", "total dap", "tab:blue"),
        ]
        # simulation
        for sid, label_analyte, color in substance:
            plots[0].add_data(
                task="task_po_dap50",
                xid="time",
                yid=sid,
                label=f"{label_analyte}, 50 mg PO",
                color=color,
            )
        # data
        for _, label_analyte, color in substance:
            dataset_name = "dapagliflozin" if label_analyte == "DAP" \
                else "total dapagliflozin"
            plots[0].add_data(
                dataset=f"{dataset_name}_DAP50",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{label_analyte}, 50 mg PO",
                color=color,
            )
        return {fig.sid: fig}


if __name__ == "__main__":
    # run_experiments(Obermeier2010, output_dir=Obermeier2010.__name__)
    out = dapagliflozin.RESULTS_PATH_SIMULATION / Obermeier2010.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Obermeier2010, output_dir=out)