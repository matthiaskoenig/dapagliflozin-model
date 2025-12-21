from typing import Dict
from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console

import pkdb_models
from pkdb_models.models import dapagliflozin
from pkdb_models.models.dapagliflozin.experiments.base_experiment import DapagliflozinSimulationExperiment
from pkdb_models.models.dapagliflozin.experiments.metadata import Tissue, Route, Dosing, ApplicationForm, Health, \
    Fasting, DapagliflozinMappingMetaData
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim
from pkdb_models.models.dapagliflozin.helpers import run_experiments
from pathlib import Path
import pkdb_models.models.dapagliflozin as dapagliflozin


class LaCreta2016(DapagliflozinSimulationExperiment):
    """Simulation experiment of LaCreta2016."""

    conditions = ["fasted", "fed"]
    groups = ["NH", "HS"]
    doses = [2.5, 10]
    bodyweights = {
        2.5: 73.5,
        10: 70.1,
    }  # [kg]

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                # unit conversion to mole/l
                if label.startswith("dapagliflozin_"):
                    dset.unit_conversion("mean", 1 / self.Mr.dap)
                dsets[f"{label}"] = dset
        # console.print(dsets)
        # console.print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        for condition in self.conditions:
            for dose in [2.5, 10]:
                tcsims[f"po_dap{dose}_{condition}"] = TimecourseSim(
                    [Timecourse(
                        start=0,
                        end=50 * 60,  # [min]
                        steps=500,
                        changes={
                            **self.default_changes(),
                            # physiological changes
                            "BW": Q_(self.bodyweights[dose], "kg"),
                            "[KI__glc_ext]": Q_(self.fpg_healthy, "mM"),
                            "GU__f_absorption": Q_(self.fasting_map[condition], "dimensionless"),
                            "f_cirrhosis": Q_(self.cirrhosis_map["Control"], "dimensionless"),
                            "KI__f_renal_function": Q_(self.renal_map["Normal renal function"], "dimensionless"),
                            # dose (IVDOSE, PODOSE)
                            "PODOSE_dap": Q_(dose, "mg"),
                        },
                    )]
                )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for condition in self.conditions:
            for group in self.groups:
                for dose in self.doses:
                    if dose == 10 and condition == "fed" and group == "NH":
                        continue
                    mappings[f"fm_po_dap{dose}_{group}_{condition}"] = FitMapping(
                        self,
                        reference=FitData(
                            self,
                            dataset=f"dapagliflozin_DAP{dose}_{group}_{condition}",
                            xid="time",
                            yid="mean",
                            yid_sd="mean_sd",
                            count="count",
                        ),
                        observable=FitData(
                            self, task=f"task_po_dap{dose}_{condition}", xid="time", yid=f"[Cve_dap]",
                        ),
                        metadata=DapagliflozinMappingMetaData(
                            tissue=Tissue.PLASMA,
                            route=Route.PO,
                            application_form=ApplicationForm.TABLET,
                            dosing=Dosing.SINGLE,
                            health=Health.HEALTHY,
                            fasting=Fasting.FASTED if condition == "fasted" else Fasting.FED,
                        ),
                    )
        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        Figure.legend_fontsize = 11
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        plots[1].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        for k, dose in enumerate(self.doses):
            for condition in self.conditions:
                for group in self.groups:
                    if dose == 10 and condition == "fed" and group == "NH":
                        continue
                    # simulation
                    plots[k].add_data(
                        task=f"task_po_dap{dose}_{condition}",
                        xid="time",
                        yid=f"[Cve_dap]",
                        label=f"{condition} - {dose} mg PO",
                        color=self.fasting_colors[condition],
                    )
                    # data
                    plots[k].add_data(
                        dataset=f"dapagliflozin_DAP{dose}_{group}_{condition}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                        label=f"{condition}, {group} - {dose} mg PO",
                        color=self.fasting_colors[condition],
                        marker="s" if group == "NH" else "D",
                    )
        return {fig.sid: fig}


if __name__ == "__main__":
    out = dapagliflozin.RESULTS_PATH_SIMULATION / LaCreta2016.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(LaCreta2016, output_dir=out)