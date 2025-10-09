from typing import Dict
from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console

import pkdb_models
from pkdb_models.models import dapagliflozin
from pkdb_models.models.dapagliflozin.experiments.base_experiment import DapagliflozinSimulationExperiment
from pkdb_models.models.dapagliflozin.experiments.metadata import Tissue, Route, Dosing, ApplicationForm, Health, \
    Fasting, DapagliflozinMappingMetaData, Coadministration
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim
from pkdb_models.models.dapagliflozin.helpers import run_experiments
from pathlib import Path
import pkdb_models.models.dapagliflozin as dapagliflozin

class Watada2019(DapagliflozinSimulationExperiment):
    """Simulation experiment of Watada2019."""

    info = [
        ("[Cve_dap]", "dapagliflozin"),
        ("[Cve_d3g]", "dapagliflozin 3-o-glucuronide"),
        ("KI__UGE", "uge"),
    ]
    interventions = ["DAP0", "DAP5", "DAP10"]
    doses = {
        "DAP0": 0,
        "DAP5": 5,
        "DAP10": 10,
    }
    bodyweights = {
        "DAP0": 57.2,
        "DAP5": 61.6,
        "DAP10": 59.8,
    }  # [kg]
    fpgs = {  # fasting plasma glucose
        "DAP0": 134/18,  # [mM]
        "DAP5": 142.9/18,  # [mM]
        "DAP10": 133.4/18,  # [mM]
    }
    gfrs = {  # fasting plasma glucose
        "DAP0": 95.4/100,  # [ml/min/(1.73*m^2)]
        "DAP5": 91.6/100,  # [ml/min/(1.73*m^2)]
        "DAP10": 94.6/100,  # [ml/min/(1.73*m^2)]
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Tab3"]:
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
        for intervention, dose in self.doses.items():
            # baseline simulation for
            tc0 = Timecourse(
                start=0,
                end=25 * 60,  # [min]
                steps=500,
                changes={
                    **self.default_changes(),
                    # physiological changes
                    "BW": Q_(self.bodyweights[intervention], "kg"),
                    "[KI__glc_ext]": Q_(self.fpgs[intervention], "mM"),
                    "KI__f_renal_function": Q_(self.gfrs[intervention], "dimensionless"),
                    "GU__f_absorption": Q_(self.fasting_map["fasted"], "dimensionless"),
                    "f_cirrhosis": Q_(self.cirrhosis_map["Control"], "dimensionless"),
                    # dose (IVDOSE, PODOSE)
                    "PODOSE_dap": Q_(0, "mg"),
                },
            )
            tc1 = Timecourse(
                start=0,
                end=24 * 60,  # [min]
                steps=500,
                changes={
                    "KI__glc_urine": Q_(0, "mmole"),  # reset UGE
                    "PODOSE_dap": Q_(dose, "mg"),
                },
            )
            tc2 = Timecourse(
                start=0,
                end=48 * 60,  # [min]
                steps=500,
                changes={
                    "KI__glc_urine": Q_(0, "mmole"),  # reset UGE
                    "PODOSE_dap": Q_(dose, "mg"),
                },
            )
            tcsims[f"po_{intervention}"] = TimecourseSim(
                [tc0] + [tc1 for _ in range(6)] + [tc2],
                time_offset=-7 * 24 * 60,
            )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for intervention, dose in self.doses.items():
            for k, (sid, name) in enumerate(self.info):
                if dose == 0 and name != "uge":
                    continue
                tissue = Tissue.URINE if name == "uge" else Tissue.PLASMA
                mappings[f"fm_po_{intervention}_{name}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_{intervention}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self, task=f"task_po_{intervention}", xid="time", yid=sid,
                    ),
                    metadata=DapagliflozinMappingMetaData(
                        tissue=tissue,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.MULTIPLE,
                        health=Health.T1DM,
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.INSULIN
                    ),
                )
        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_rows=1,
            num_cols=3,
            name=f"{self.__class__.__name__}",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        plots[1].set_yaxis(self.label_d3g_plasma, unit=self.unit_d3g)
        plots[2].set_yaxis(self.label_uge, unit=self.unit_uge)
        plots[0].xaxis.min = -24
        plots[0].xaxis.max = 30
        plots[1].xaxis.min = -24
        plots[1].xaxis.max = 30

        for intervention, dose in self.doses.items():
            for k, (sid, name) in enumerate(self.info):
                # simulation
                if dose == 0 and k in [0, 1]:
                    pass
                else:
                    plots[k].add_data(
                        task=f"task_po_{intervention}",
                        xid="time",
                        yid=sid,
                        label=f"{dose} mg",
                        color=self.dose_colors[dose],
                    )
                # data
                if dose == 0 and name != "uge":
                    continue
                plots[k].add_data(
                    dataset=f"{name}_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"{dose} mg",
                    linestyle="" if name == "uge" else "--",
                    color=self.dose_colors[dose],
                )
        return {fig.sid: fig,}


if __name__ == "__main__":
    # run_experiments(Watada2019, output_dir=Watada2019.__name__)
    out = dapagliflozin.RESULTS_PATH_SIMULATION / Watada2019.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Watada2019, output_dir=out)