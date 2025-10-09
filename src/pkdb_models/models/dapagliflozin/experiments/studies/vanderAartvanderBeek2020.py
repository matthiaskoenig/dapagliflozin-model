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


class vanderAartvanderBeek2020(DapagliflozinSimulationExperiment):
    """
    Simulation experiment of vanderAartvanderBeek2020.
    Renal disease & macroalbunimurea.
    """

    subjects_markers = {
        "S1": "o",  # circle
        "S2": "s",  # square
        "S3": "D",  # diamond
        "S4": "^",  # triangle up
        "S5": "v",  # triangle down
        "S6": "<",  # triangle left
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig3"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                # unit conversion
                if label.startswith("dapagliflozin_"):
                    dset.unit_conversion("value", 1 / self.Mr.dap)
                dsets[f"{label}"] = dset
        # console.print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        tc0 = Timecourse(
            start=0,
            end=24 * 60,  # [min]
            steps=500,
            changes={
                **self.default_changes(),
                # physiological changes
                "BW": Q_(self.bodyweight_default, "kg"),
                "[KI__glc_ext]": Q_(self.fpg_healthy, "mM"),
                "GU__f_absorption": Q_(self.fasting_map["not reported"], "dimensionless"),
                "f_cirrhosis": Q_(self.cirrhosis_map["Control"], "dimensionless"),
                "KI__f_renal_function": Q_(self.renal_map["Mild renal impairment"], "dimensionless"),
                # dose (IVDOSE, PODOSE)
                "PODOSE_dap": Q_(10, "mg"),
            },
        )
        tc1 = Timecourse(
            start=0,
            end=24 * 60,  # [min]
            steps=500,
            changes={
                "PODOSE_dap": Q_(10, "mg"),
            },
        )
        tcsims[f"po_dap10"] = TimecourseSim(
            [tc0] + [tc1 for _ in range(9)],
            time_offset=-9 * 24 * 60,
        )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for subject in self.subjects_markers.keys():
            mappings[f"fm_dap10_{subject}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"dapagliflozin_DAP_{subject}",
                    xid="time",
                    yid="value",
                    yid_sd=None,
                    count="count",
                ),
                observable=FitData(
                    self, task=f"task_po_dap10", xid="time", yid=f"[Cve_dap]",
                ),
                metadata=DapagliflozinMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.MULTIPLE,
                    health=Health.RENAL_IMPAIRMENT,
                    fasting=Fasting.NR,
                ),
            )
        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig3",
            name=f"{self.__class__.__name__}",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        plots[0].xaxis.min = -3
        plots[0].xaxis.max = 24
        # simulation
        plots[0].add_data(
            task="task_po_dap10",
            xid="time",
            yid="[Cve_dap]",
            label="10 mg PO",
            color="black",
        )
        # data: individual subjects
        for subject, marker in self.subjects_markers.items():
            plots[0].add_data(
                dataset=f"dapagliflozin_DAP_{subject}",
                xid="time",
                yid="value",
                yid_sd=None,
                count="count",
                color="#2ca25f",
                label=f"10 mg PO ({subject})",
                marker=marker,
            )
        return {fig.sid: fig}


if __name__ == "__main__":
    # run_experiments(vanderAartvanderBeek2020, output_dir=vanderAartvanderBeek2020.__name__)
    out = dapagliflozin.RESULTS_PATH_SIMULATION / vanderAartvanderBeek2020.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(vanderAartvanderBeek2020, output_dir=out)