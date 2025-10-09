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

class Kasichayanula2011(DapagliflozinSimulationExperiment):
    """
    Simulation experiment of Kasichayanula2011.
    Hepatic impairment.
    """

    conditions = {
        "healthy": "Control",
        "mild": "Mild cirrhosis",
        "moderate": "Moderate cirrhosis",
        "severe": "Severe cirrhosis",
    }
    bodyweights = {
        "healthy": 81.8,
        "mild": 80.9,
        "moderate": 87.5,
        "severe": 89.2,
    }  # [kg]
    renal_functions = {  # creatine clearance  # FIXME: Check why values so high
        "healthy":  138.3 / 100,
        "mild": 127.8 / 100,
        "moderate": 158.3 / 100,
        "severe": 142.8/ 100,
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig2"]:
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
        for condition, cirrhosis_class in self.conditions.items():
            tcsims[f"po_dap10_{condition}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=50 * 60,  # [min]
                    steps=500,
                    changes={
                        **self.default_changes(),
                        # physiological changes
                        "BW": Q_(self.bodyweights[condition], "kg"),
                        "[KI__glc_ext]": Q_(self.fpg_healthy, "mM"),
                        "GU__f_absorption": Q_(self.fasting_map["fasted"], "dimensionless"),
                        "f_cirrhosis": Q_(self.cirrhosis_map[cirrhosis_class], "dimensionless"),
                        "KI__f_renal_function": Q_(self.renal_functions[condition], "dimensionless"),
                        # dose (IVDOSE, PODOSE)
                        "PODOSE_dap": Q_(10, "mg"),
                    },
                )]
            )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for condition in self.conditions:
            health = Health.HEALTHY if condition == "healthy" else Health.CIRRHOSIS
            for sid, name in [("dap", "dapagliflozin"), ("d3g", "dapagliflozin 3-o-glucuronide")]:
                mappings[f"fm_dap10_{condition}_{sid}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_DAP10_{condition}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self, task=f"task_po_dap10_{condition}", xid="time", yid=f"[Cve_{sid}]",
                    ),
                    metadata=DapagliflozinMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=health,
                        fasting=Fasting.FASTED,
                    ),
                )
        # console.print(mappings)
        return mappings

    def figures(self) -> Dict[str, Figure]:
        colors = {k: self.cirrhosis_colors[v] for k, v in self.conditions.items()}
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__}",
        )
        Figure.legend_fontsize=10
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        plots[1].set_yaxis(self.label_d3g_plasma, unit=self.unit_d3g)

        for condition in self.conditions:
            dose_label = f"{condition} - 10 mg PO"
            # simulation: dapagliflozin plasma
            plots[0].add_data(
                task=f"task_po_dap10_{condition}",
                xid="time",
                yid="[Cve_dap]",
                label=f"{dose_label}",
                color=colors[condition],
            )
            # simulation: d3g plasma
            plots[1].add_data(
                task=f"task_po_dap10_{condition}",
                xid="time",
                yid="[Cve_d3g]",
                label=f"{dose_label}",
                color=colors[condition],
            )
            # data: dapagliflozin plasma
            plots[0].add_data(
                dataset=f"dapagliflozin_DAP10_{condition}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{dose_label}",
                color=colors[condition],
            )
            # data: d3g plasma
            plots[1].add_data(
                dataset=f"dapagliflozin 3-o-glucuronide_DAP10_{condition}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{dose_label}",
                color=colors[condition],
            )
        return {fig.sid: fig}


if __name__ == "__main__":
    # run_experiments(Kasichayanula2011, output_dir=Kasichayanula2011.__name__)
    out = dapagliflozin.RESULTS_PATH_SIMULATION / Kasichayanula2011.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Kasichayanula2011, output_dir=out)