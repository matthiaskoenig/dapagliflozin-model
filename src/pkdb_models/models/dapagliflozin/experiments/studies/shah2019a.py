from typing import Dict
from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console
from pkdb_models.models.dapagliflozin.experiments.base_experiment import DapagliflozinSimulationExperiment
from pkdb_models.models.dapagliflozin.experiments.metadata import Tissue, Route, Dosing, ApplicationForm, Health, \
    Fasting, DapagliflozinMappingMetaData, Coadministration
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim
from pkdb_models.models.dapagliflozin.helpers import run_experiments
from pathlib import Path
import pkdb_models.models.dapagliflozin as dapagliflozin

class Shah2019a(DapagliflozinSimulationExperiment):
    """Simulation experiment of Shah2019a."""

    conditions = ["fasted", "fed"]
    bodyweight = (54.6 + 82.2)/2  # [kg] mean out of (min 54.6; max 82.2)

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig3"]:
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
            tcsims[f"po_dap5_{condition}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=100 * 60,  # [min]
                    steps=500,
                    changes={
                        **self.default_changes(),

                        # physiological changes
                        "BW": Q_(self.bodyweight, "kg"),
                        "[KI__glc_ext]": Q_(self.fpg_healthy, "mM"),
                        "GU__f_absorption": Q_(self.fasting_map[condition], "dimensionless"),
                        "f_cirrhosis": Q_(self.cirrhosis_map["Control"], "dimensionless"),
                        "KI__f_renal_function": Q_(self.renal_map["Normal renal function"], "dimensionless"),

                        # dose (IVDOSE, PODOSE)
                        "PODOSE_dap": Q_(5, "mg"),

                    },
                )]
            )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for condition in self.conditions:
                mappings[f"fm_po_dap5_{condition}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"dapagliflozin_DAP5, MET500_{condition}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self, task=f"task_po_dap5_{condition}", xid="time", yid=f"[Cve_dap]",
                    ),
                    metadata=DapagliflozinMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        coadministration=Coadministration.METFORMIN,
                        fasting=Fasting.FASTED if condition == "fasted" else Fasting.FED,
                    ),
                )
        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_rows=1,
            num_cols=1,
            name=f"{self.__class__.__name__}",
        )
        Figure.legend_fontsize = 11
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)

        for condition in self.conditions:
            # simulation
            plots[0].add_data(
                task=f"task_po_dap5_{condition}",
                xid="time",
                yid="[Cve_dap]",
                label=f"{condition} - 5 mg PO",
                color=self.fasting_colors[condition],
            )
            # data
            plots[0].add_data(
                dataset=f"dapagliflozin_DAP5, MET500_{condition}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{condition} - 5 mg PO + MET500",
                color=self.fasting_colors[condition],
            )

        return {fig.sid: fig}


if __name__ == "__main__":
    # run_experiments(Shah2019a, output_dir=Shah2019a.__name__)
    out = dapagliflozin.RESULTS_PATH_SIMULATION / Shah2019a.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Shah2019a, output_dir=out)