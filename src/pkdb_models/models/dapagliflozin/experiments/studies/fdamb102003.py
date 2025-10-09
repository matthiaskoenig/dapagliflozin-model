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

class FDAMB102003(DapagliflozinSimulationExperiment):
    """Simulation experiment of FDA_MB102003."""

    doses = [5, 25, 100]  # placebo uge, i.e. 0 dose removed

    info = {
        "[Cve_dap]": "dapagliflozin",
        "Aurine_dap": "dapagliflozin_urine",
        "KI__UGE": "uge",
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        # Fig13: dap
        # Tab8A: dap_urine
        # Fig18: UGE
        for fig_id in ["Fig13", "Fig18", "Tab8A"]:
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
        for dose in self.doses:
            tc0 = Timecourse(
                start=0,
                end=24 * 60,  # [min]
                steps=500,
                changes={
                    **self.default_changes(),
                    # physiological changes
                    "BW": Q_(self.bodyweight_default, "kg"),
                    "[KI__glc_ext]": Q_(self.fpg_t2dm, "mM"),
                    "GU__f_absorption": Q_(self.fasting_map["fasted"], "dimensionless"),
                    "f_cirrhosis": Q_(self.cirrhosis_map["Control"], "dimensionless"),
                    "KI__f_renal_function": Q_(self.renal_map["Normal renal function"], "dimensionless"),
                    # dose (IVDOSE, PODOSE)
                    "PODOSE_dap": Q_(dose, "mg"),
                },
            )
            tc1 = Timecourse(
                start=0,
                end=24 * 60,  # [min]
                steps=500,
                changes={
                    "KI__glc_urine": Q_(0, "mmole"),  # reset UGE
                    "Aurine_dap": Q_(0, "mmole"),  # reset urinary dap
                    "PODOSE_dap": Q_(dose, "mg"),
                },
            )
            tc2 = Timecourse(
                start=0,
                end=25 * 60,  # [min]
                steps=500,
                changes={
                    "KI__glc_urine": Q_(0, "mmole"),  # reset UGE
                    "Aurine_dap": Q_(0, "mmole"),  # reset urinary dap
                    "PODOSE_dap": Q_(dose, "mg"),
                },
            )
            tcsims[f"po_DAP{dose}"] = TimecourseSim(
                [tc0] + [tc1 for _ in range(12)] + [tc2],
                # time_offset=-13*24*60,
            )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for dose in self.doses:
            for kp, sid in enumerate(self.info):
                name = self.info[sid]
                if name.startswith("dapagliflozin") and dose == 0:
                    continue
                mappings[f"fm_dap{dose}_{name}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_MAD{dose}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self, task=f"task_po_DAP{dose}", xid="time", yid=sid,
                    ),
                    metadata=DapagliflozinMappingMetaData(
                        tissue=Tissue.PLASMA if sid == "[Cve_dap]" else Tissue.URINE,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.MULTIPLE,
                        health=Health.T2DM,
                        fasting=Fasting.FASTED,
                        outlier=True if sid == "KI__UGE" else None  # FIXME: uncontrolled T2DM (probably much higher glc)
                    ),
                )
        return mappings

    def figures(self) -> Dict[str, Figure]:
        figures = {}
        subplots = ["all", "start", "end"]
        for subplot in subplots:
            fig = Figure(
                experiment=self,
                sid=f"Fig1_{subplot}",
                num_rows=1,
                num_cols=3,
                name=f"{self.__class__.__name__} ({subplot})"
            )
            plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
            plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
            plots[1].set_yaxis(self.label_dap_urine, unit=self.unit_dap_urine)
            plots[2].set_yaxis(self.label_uge, unit=self.unit_uge)

            for kp in [0, 1, 2]:
                if subplot == "start":
                    plots[kp].xaxis.min = -3
                    plots[kp].xaxis.max = 30
                elif subplot == "end":
                    plots[kp].xaxis.min = 310
                    plots[kp].xaxis.max = 340

            for dose in self.doses:
                for kp, sid in enumerate(self.info):
                    name = self.info[sid]
                    if name.startswith("dapagliflozin") and dose == 0:
                        continue
                    dose_label = f"{dose} mg PO"
                    # simulation
                    plots[kp].add_data(
                        task=f"task_po_DAP{dose}",
                        xid="time",
                        yid=sid,
                        label=f"{dose_label}",
                        color=self.dose_colors[dose],
                    )
                    # data
                    plots[kp].add_data(
                        dataset=f"{name}_MAD{dose}",
                        xid="time",
                        yid="mean",
                        yid_sd=None if sid == "KI__UGE" else "mean_sd",
                        label=f"{dose_label}",
                        color=self.dose_colors[dose],
                        linestyle="--" if name == "dapagliflozin" else ""
                    )
                    # print(sid)
            figures[fig.sid] = fig

        return figures


if __name__ == "__main__":
    # run_experiments(FDAMB102003, output_dir=FDAMB102003.__name__)
    out = dapagliflozin.RESULTS_PATH_SIMULATION / FDAMB102003.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(FDAMB102003, output_dir=out)