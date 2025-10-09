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


class Kasichayanula2011a(DapagliflozinSimulationExperiment):
    """Simulation experiment of Kasichayanula2011a."""

    doses_mad = [0, 2.5, 10, 20]
    bodyweights = {
        # "SAD_0": 64.3,
        # "SAD_2.5": 64.7,
        # "SAD_10": 60,
        # "SAD_20": 64.8,
        # "SAD_50": 67.4,
        "MAD0": 61.0,
        "MAD2.5": 60.8,
        "MAD10": 60.2,
        "MAD20": 67.4,
    }  # [kg]
    groups = list(bodyweights.keys())

    #  creatinine clearance overestimates GFR by approximately 10-20% due to tubular secretion of creatinine.
    # A simple approximation often used in clinical practice is:
    # GFR ≈ CrCl×0.85
    f_gfr = 0.85
    renal_functions = {  # eGFR from creatinine clearance
        # "SAD_0": 135.6 / 100 * f_gfr,
        # "SAD_2.5": 151.6 / 100 * f_gfr,
        # "SAD_10": 150.2 / 100 * f_gfr,
        # "SAD_20": 147.7 / 100 * f_gfr,
        # "SAD_50": 154.3 / 100 * f_gfr,
        "MAD0": 103 / 100 * f_gfr,
        "MAD2.5": 112 / 100 * f_gfr,
        "MAD10": 103 / 100 * f_gfr,
        "MAD20": 103 / 100 * f_gfr,
    }  # [ml/min]

    # assuming normal glucose values for healthy subjects

    # fpgs = {
    #     "SAD_0": fpg_healthy,
    #     "SAD_2.5": fpg_healthy,
    #     "SAD_10": fpg_healthy,
    #     "SAD_20": fpg_healthy,
    #     "SAD_50": fpg_healthy,
    #     "MAD_0": fpg_t2dm,
    #     "MAD_2.5": fpg_t2dm,
    #     "MAD_10": fpg_t2dm,
    #     "MAD_20": fpg_t2dm,
    # }  # [mM]

    info = {
        "[Cve_dap]": "dapagliflozin",
        "[Cve_d3g]": "dapagliflozin 3-o-glucuronide",
        "KI__UGE": "uge",
    }


    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        # Fig2: dap & d3g
        # Fig3: UGE
        # Tab2A: dap recovery urine
        # Fig4: UGE individual
        for fig_id in ["Fig2", "Fig3", "Tab2A"]:  #"Tab2A", "Fig4"
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
        # multiple dose
        for dose in self.doses_mad:
            tc0 = Timecourse(
                start=0,
                end=24 * 60,  # [min]
                steps=500,
                changes={
                    **self.default_changes(),
                    # physiological changes
                    "BW": Q_(self.bodyweights[f"MAD{dose}"], "kg"),
                    "[KI__glc_ext]": Q_(self.fpg_t2dm, "mM"),
                    "GU__f_absorption": Q_(self.fasting_map["fasted"], "dimensionless"),
                    "f_cirrhosis": Q_(self.cirrhosis_map["Control"], "dimensionless"),
                    "KI__f_renal_function": Q_(self.renal_functions[f"MAD{dose}"], "dimensionless"),
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
            tcsims[f"po_MAD{dose}"] = TimecourseSim(
                [tc0] + [tc1 for _ in range(12)] + [tc2],
                # time_offset=-13*24*60,
            )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for dose in self.doses_mad:
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
                        self, task=f"task_po_MAD{dose}", xid="time", yid=sid,
                    ),
                    metadata=DapagliflozinMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.T2DM,
                        fasting=Fasting.FASTED,
                    ),
                )
        return mappings

    def figures(self) -> Dict[str, Figure]:
        return {
            **self.figure_mad(),
        }

    def figure_mad(self) -> Dict[str, Figure]:
        figures = {}
        subplots = ["all", "start", "end"]
        for subplot in subplots:
            fig = Figure(
                experiment=self,
                sid=f"Fig1_{subplot}",
                num_rows=1,
                num_cols=3,
                name=f"{self.__class__.__name__} ({subplot})",
            )
            plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
            plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
            plots[1].set_yaxis(self.label_d3g_plasma, unit=self.unit_d3g)
            plots[2].set_yaxis(self.label_uge, unit=self.unit_uge)

            for kp in [0, 1, 2]:
                if subplot == "start":
                    plots[kp].xaxis.min = -3
                    plots[kp].xaxis.max = 30
                elif subplot == "end":
                    plots[kp].xaxis.min = 310
                    plots[kp].xaxis.max = 340

            for dose in self.doses_mad:
                for kp, sid in enumerate(self.info):
                    name = self.info[sid]
                    if name.startswith("dapagliflozin") and dose == 0:
                        continue
                    dose_label = f"{dose} mg PO"
                    # simulation
                    plots[kp].add_data(
                        task=f"task_po_MAD{dose}",
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
                        yid_sd="mean_sd",
                        label=f"{dose_label}",
                        color=self.dose_colors[dose],
                        linestyle="" if name == "uge" else "--"
                    )
            figures[fig.sid] = fig

        return figures


if __name__ == "__main__":
    # run_experiments(Kasichayanula2011a, output_dir=Kasichayanula2011a.__name__)
    out = dapagliflozin.RESULTS_PATH_SIMULATION / Kasichayanula2011a.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Kasichayanula2011a, output_dir=out)