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


class Yang2013(DapagliflozinSimulationExperiment):
    """Simulation experiment of Yang2013."""

    doses_single = [5, 10]
    doses_multi = [5, 10]
    info = [
        ("[Cve_dap]", "dapagliflozin"),
        ("[Cve_d3g]", "dapagliflozin 3-o-glucuronide"),
        ("KI__UGE", "glucose_cumulative amount")
    ]
    bodyweights = {
        5: 62.8,
        10: 62.2,
    }  # [kg]

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
        # single dose simulation
        for dose in self.doses_single:
            tcsims[f"po_dap{dose}_single"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=100 * 60,  # [min]
                    steps=1000,
                    changes={
                        **self.default_changes(),
                        # physiological changes
                        "BW": Q_(self.bodyweights[dose], "kg"),
                        "[KI__glc_ext]": Q_(self.fpg_healthy, "mM"),
                        "GU__f_absorption": Q_(self.fasting_map["fasted"], "dimensionless"),
                        "f_cirrhosis": Q_(self.cirrhosis_map["Control"], "dimensionless"),
                        "KI__f_renal_function": Q_(self.renal_map["Normal renal function"], "dimensionless"),
                        # dose
                        "PODOSE_dap": Q_(dose, "mg"),
                    },
                )]
            )
        # multi dose simulation
        for dose in self.doses_multi:
            tc0 = Timecourse(
                start=0,
                end=24 * 60,  # [min]
                steps=50,
                changes={
                    **self.default_changes(),
                    ## physiological changes
                    "[KI__glc_ext]": Q_(self.fpg_healthy, "mM"),
                    "GU__f_absorption": Q_(self.fasting_map["fasted"], "dimensionless"),
                    "f_cirrhosis": Q_(self.cirrhosis_map["Control"], "dimensionless"),
                    "KI__f_renal_function": Q_(self.renal_map["Normal renal function"], "dimensionless"),
                    # dose
                    "PODOSE_dap": Q_(dose, "mg"),
                },
            )
            tc1 = Timecourse(
                start=0,
                end=24 * 60,  # [min]
                steps=50,
                changes={
                    "KI__glc_urine": Q_(0, "mmole"),  # reset UGE
                    "PODOSE_dap": Q_(dose, "mg"),
                },
            )
            tc2 = Timecourse(
                start=0,
                end=100 * 60,  # [min]
                steps=2000,
                changes={
                    "KI__glc_urine": Q_(0, "mmole"),  # reset UGE
                    "PODOSE_dap": Q_(dose, "mg"),
                },
            )
            tcsims[f"po_dap{dose}_multi"] = TimecourseSim(
                [tc0] + [tc1 for _ in range(8)] + [tc2],
                time_offset=-9 * 24 * 60,
            )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        # Fig1
        for kp, data in enumerate(self.info):
            sid, prefix = data[0], data[1]
            for dose in self.doses_single:
                tissue = Tissue.URINE if "cumulative" in prefix else Tissue.PLASMA
                mappings[f"task_po_dap{dose}_{sid}_single"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{prefix}_SD_DAP{dose}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self, task=f"task_po_dap{dose}_single", xid="time", yid=sid,
                    ),
                    metadata=DapagliflozinMappingMetaData(
                        tissue=tissue,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                    ),
                )
        # Fig2
        for kp, data in enumerate(self.info):
            sid, prefix = data[0], data[1]
            for dose in self.doses_multi:
                tissue = Tissue.URINE if "cumulative" in prefix else Tissue.PLASMA
                mappings[f"task_po_dap{dose}_{sid}_multi"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{prefix}_MD_DAP{dose}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self, task=f"task_po_dap{dose}_multi", xid="time", yid=sid,
                    ),
                    metadata=DapagliflozinMappingMetaData(
                        tissue=tissue,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.MULTIPLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                    ),
                )
        # console.print(mappings)
        return mappings

    def figures(self) -> Dict[str, Figure]:
        return {
            **self.figure_single(),
            **self.figure_multi(),
        }

    def figure_single(self) -> Dict[str, Figure]:
        figs: Dict[str, Figure] = {}

        # Fig1 — DAP plasma
        fig1_plasma = Figure(
            experiment=self,
            sid="Fig1_plasma",
            num_rows=1,
            num_cols=1,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        p = fig1_plasma.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)[0]
        p.set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        p.xaxis.min = -2
        p.xaxis.max = 30
        for dose in self.doses_single:
            p.add_data(
                task=f"task_po_dap{dose}_single",
                xid="time",
                yid="[Cve_dap]",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
            )
            p.add_data(
                dataset=f"dapagliflozin_SD_DAP{dose}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
            )
        figs[fig1_plasma.sid] = fig1_plasma

        # Fig1 — D3G plasma
        fig1_d3g = Figure(
            experiment=self,
            sid="Fig1_d3g",
            num_rows=1,
            num_cols=1,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        p = fig1_d3g.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)[0]
        p.set_yaxis(self.label_d3g_plasma, unit=self.unit_d3g)
        p.xaxis.min = -2
        p.xaxis.max = 30
        for dose in self.doses_single:
            p.add_data(
                task=f"task_po_dap{dose}_single",
                xid="time",
                yid="[Cve_d3g]",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
            )
            p.add_data(
                dataset=f"dapagliflozin 3-o-glucuronide_SD_DAP{dose}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
            )
        figs[fig1_d3g.sid] = fig1_d3g

        # Fig1 — UGE
        fig1_uge = Figure(
            experiment=self,
            sid="Fig1_uge",
            num_rows=1,
            num_cols=1,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        p = fig1_uge.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)[0]
        p.set_yaxis(self.label_uge, unit=self.unit_uge)
        for dose in self.doses_single:
            p.add_data(
                task=f"task_po_dap{dose}_single",
                xid="time",
                yid="KI__UGE",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
            )
            p.add_data(
                dataset=f"glucose_cumulative amount_SD_DAP{dose}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
            )
        figs[fig1_uge.sid] = fig1_uge

        return figs

    def figure_multi(self) -> Dict[str, Figure]:
        figs: Dict[str, Figure] = {}

        # Fig2 — DAP plasma
        fig2_plasma = Figure(
            experiment=self,
            sid="Fig2_plasma",
            num_rows=1,
            num_cols=1,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        p = fig2_plasma.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)[0]
        p.set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        p.xaxis.min = -24
        p.xaxis.max = 30
        for dose in self.doses_multi:
            p.add_data(
                task=f"task_po_dap{dose}_multi",
                xid="time",
                yid="[Cve_dap]",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
            )
            p.add_data(
                dataset=f"dapagliflozin_MD_DAP{dose}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
            )
        figs[fig2_plasma.sid] = fig2_plasma

        # Fig2 — D3G plasma
        fig2_d3g = Figure(
            experiment=self,
            sid="Fig2_d3g",
            num_rows=1,
            num_cols=1,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        p = fig2_d3g.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)[0]
        p.set_yaxis(self.label_d3g_plasma, unit=self.unit_d3g)
        p.xaxis.min = -24
        p.xaxis.max = 30
        for dose in self.doses_multi:
            p.add_data(
                task=f"task_po_dap{dose}_multi",
                xid="time",
                yid="[Cve_d3g]",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
            )
            p.add_data(
                dataset=f"dapagliflozin 3-o-glucuronide_MD_DAP{dose}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
            )
        figs[fig2_d3g.sid] = fig2_d3g

        # Fig2 — UGE
        fig2_uge = Figure(
            experiment=self,
            sid="Fig2_uge",
            num_rows=1,
            num_cols=1,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        p = fig2_uge.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)[0]
        p.set_yaxis(self.label_uge, unit=self.unit_uge)
        for dose in self.doses_multi:
            p.add_data(
                task=f"task_po_dap{dose}_multi",
                xid="time",
                yid="KI__UGE",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
            )
            p.add_data(
                dataset=f"glucose_cumulative amount_MD_DAP{dose}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
            )
        figs[fig2_uge.sid] = fig2_uge

        return figs


if __name__ == "__main__":
    out = dapagliflozin.RESULTS_PATH_SIMULATION / Yang2013.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Yang2013, output_dir=out)