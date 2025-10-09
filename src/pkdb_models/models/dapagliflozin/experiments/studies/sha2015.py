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


class Sha2015(DapagliflozinSimulationExperiment):
    """Simulation experiment of Sha2015."""

    colors = {10: "black"}
    bodyweight = 78.9  # [kg]
    gfr = 96.8  # [ml/min/1.73*m^2]

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2", "Fig3"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                dsets[f"{label}"] = dset
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        # single dose
        tcsims["po_dap10"] = TimecourseSim(
            Timecourse(
                start=0,
                end=25 * 60,  # [min]
                steps=500,
                changes={
                    **self.default_changes(),
                    # physiological changes
                    "BW": Q_(self.bodyweight, "kg"),
                    "[KI__glc_ext]": Q_(self.fpg_healthy, "mM"),
                    "KI__f_renal_function": Q_(self.gfr / 100, "dimensionless"),  # [0,1] ↔ [0,100] gfr
                    "GU__f_absorption": Q_(self.fasting_map["fasted"], "dimensionless"),
                    "f_cirrhosis": Q_(self.cirrhosis_map["Control"], "dimensionless"),
                    # dose (IVDOSE, PODOSE)
                    "PODOSE_dap": Q_(10, "mg"),
                },
            )
        )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for name, sid in [
            ("dapagliflozin", "[Cve_dap]"),
            ("glucose_cumulative_amount", "KI__UGE"),
            ("renal_threshold", "KI__RTG"),
        ]:
            mappings[f"task_po_dap10_{name}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"{name}_DAP10",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self, task="task_po_dap10", xid="time", yid=sid,
                ),
                metadata=DapagliflozinMappingMetaData(
                    tissue=Tissue.URINE if "cumulative" in name else Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
                ),
            )
        return mappings

    def figures(self) -> Dict[str, Figure]:
        figs: Dict[str, Figure] = {}

        # Figure plasma
        fig_plasma = Figure(
            experiment=self,
            sid="Sha2015_plasma",
            num_rows=1,
            num_cols=1,
            name=f"{self.__class__.__name__}",
        )
        p_plasma = fig_plasma.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )[0]
        p_plasma.set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        # simulation
        p_plasma.add_data(
            task="task_po_dap10",
            xid="time",
            yid="[Cve_dap]",
            label="10 mg PO",
            color=self.colors[10],
        )
        # data
        p_plasma.add_data(
            dataset="dapagliflozin_DAP10",
            xid="time",
            yid="mean",
            yid_sd="mean_sd",
            count="count",
            label="10 mg PO",
            color=self.colors[10],
        )
        figs[fig_plasma.sid] = fig_plasma

        # Figure UGE
        fig_uge = Figure(
            experiment=self,
            sid="Sha2015_uge",
            num_rows=1,
            num_cols=1,
            name=f"{self.__class__.__name__}",
        )
        p_uge = fig_uge.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )[0]
        p_uge.set_yaxis(self.label_uge, unit=self.unit_uge)
        # simulation
        p_uge.add_data(
            task="task_po_dap10",
            xid="time",
            yid="KI__UGE",
            label="10 mg PO",
            color=self.colors[10],
        )
        # data
        p_uge.add_data(
            dataset="glucose_cumulative_amount_DAP10",
            xid="time",
            yid="mean",
            yid_sd="mean_sd",
            count="count",
            label="10 mg PO",
            color=self.colors[10],
        )
        figs[fig_uge.sid] = fig_uge

        # Figure RTG
        fig_rtg = Figure(
            experiment=self,
            sid="Sha2015_rtg",
            num_rows=1,
            num_cols=1,
            name=f"{self.__class__.__name__}",
        )
        p_rtg = fig_rtg.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )[0]
        p_rtg.set_yaxis(self.label_rtg, unit=self.unit_rtg)
        # simulation
        p_rtg.add_data(
            task="task_po_dap10",
            xid="time",
            yid="KI__RTG",
            label="10 mg PO",
            color=self.colors[10],
        )
        # data
        p_rtg.add_data(
            dataset="renal_threshold_DAP10",
            xid="time",
            yid="mean",
            yid_sd="mean_sd",
            count="count",
            label="10 mg PO",
            color=self.colors[10],
        )
        figs[fig_rtg.sid] = fig_rtg

        return figs


if __name__ == "__main__":
    # run_experiments(Sha2015, output_dir=Sha2015.__name__)
    out = dapagliflozin.RESULTS_PATH_SIMULATION / Sha2015.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Sha2015, output_dir=out)