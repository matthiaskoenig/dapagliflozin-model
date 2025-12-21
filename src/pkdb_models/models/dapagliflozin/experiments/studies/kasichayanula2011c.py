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


class Kasichayanula2011c(DapagliflozinSimulationExperiment):
    """Simulation experiment of Kasichayanula2011c."""

    group_doses = {
        "S1": 50,
        "S2": 20,
        "S3": 20,
    }
    group_interventions = {
        "S1": ["DAP50", "DAP50, PIO45"],
        "S2": ["DAP20", "DAP20, MET1000"],
        "S3": ["DAP20", "DAP20, GLI4", "DAP20, SIT100"],
    }
    intervention_colors = {
        "DAP50": "black",
        "DAP50, PIO45": "tab:blue",
        "DAP20": "black",
        "DAP20, MET1000": "tab:blue",
        "DAP20, GLI4": "tab:blue",
        "DAP20, SIT100": "tab:orange",
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig3"]:
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
        for dose in set(self.group_doses.values()):
            tcsims[f"po_dap{dose}"] = TimecourseSim(
                    [Timecourse(
                        start=0,
                        end=75 * 60,  # [min]
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
                            "PODOSE_dap": Q_(dose, "mg")
                        },
                    )]
                )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for group, interventions in self.group_interventions.items():
            dose = self.group_doses[group]
            for intervention in interventions:
                mappings[f"fm_po_{intervention}_{group}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"dapagliflozin_{intervention}_{group}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self, task=f"task_po_dap{dose}", xid="time", yid=f"[Cve_dap]",
                    ),
                    metadata=DapagliflozinMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        coadministration=(
                            Coadministration.PIOGLITAZONE if "PIO" in intervention
                            else Coadministration.METFORMIN if "MET" in intervention
                            else Coadministration.GLIMEPIRIDE if "GLI" in intervention
                            else Coadministration.SITAGLIPTIN if "SIT" in intervention
                            else Coadministration.NONE
                        )
                    ),
                )
        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig3",
            num_rows=1,
            num_cols=3,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        Figure.legend_fontsize = 11
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        for idx, (group, interventions) in enumerate(self.group_interventions.items()):
            plots[idx].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
            dose = self.group_doses[group]
            dose_label = f"{dose} mg PO"
            # simulation
            plots[idx].add_data(
                task=f"task_po_dap{dose}",
                xid="time",
                yid="[Cve_dap]",
                label=f"{dose_label}",
                color="black",
            )
            # data
            for intervention in interventions:
                if intervention == f"DAP{dose}":
                    label = f"{dose_label}"
                else:
                    label = f"{dose} mg PO + {intervention.split(', ')[1]}"
                plots[idx].add_data(
                    dataset=f"dapagliflozin_{intervention}_{group}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=label,
                    color=self.intervention_colors[intervention],
                )
        return {fig.sid: fig}


if __name__ == "__main__":
    out = dapagliflozin.RESULTS_PATH_SIMULATION / Kasichayanula2011c.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Kasichayanula2011c, output_dir=out)