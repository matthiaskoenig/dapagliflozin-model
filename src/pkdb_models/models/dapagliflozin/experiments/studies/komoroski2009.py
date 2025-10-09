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
from pathlib import  Path
import pkdb_models.models.dapagliflozin as dapagliflozin


class Komoroski2009(DapagliflozinSimulationExperiment):
    """Simulation experiment of Komoroski2009."""
    #FIXME: Fig6

    conditions = ["fasted", "fed"]
    doses_multi = [2.5, 10, 20, 50, 100]
    doses_uge_single = [0, 2.5, 5, 10, 20, 50, 100, 250, 500]
    doses_uge_multi = [2.5, 10, 20, 50, 100]
    colors = {
        0: "black",
        2.5: "tab:blue",
        5: "tab:orange",
        10: "tab:green",
        20: "tab:red",
        50: "tab:purple",
        100: "tab:brown",
        250: "tab:pink",
        500: "tab:grey",
    }
    bodyweight = (74 + 84)/2  # [kg] range [74, 84]

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig3", "Tab1A", "Fig4", "Fig5", "Fig6"]:
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
        # fasting simulation
        for condition in self.conditions:
            tcsims[f"po_dap_250_{condition}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=25 * 60,  # [min]
                    steps=500,
                    changes={
                        **self.default_changes(),
                        # physiological changes
                        "BW": Q_(self.bodyweight, "kg"),
                        "GU__f_absorption": Q_(self.fasting_map[condition], "dimensionless"),
                        "[KI__glc_ext]": Q_(self.fpg_healthy, "mM"),
                        "f_cirrhosis": Q_(self.cirrhosis_map["Control"], "dimensionless"),
                        "KI__f_renal_function": Q_(self.renal_map["Normal renal function"], "dimensionless"),
                        # dose (IVDOSE, PODOSE)
                        "PODOSE_dap": Q_(250, "mg"),
                    },
                )]
            )
        # single dose simulation
        for dose in self.doses_uge_single:
            tcsims[f"po_dap_{dose}_single"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=125 * 60,  # [min]
                    steps=1000,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(self.bodyweight, "kg"),
                        "GU__f_absorption": Q_(self.fasting_map["fasted"], "dimensionless"),
                        "[KI__glc_ext]": Q_(self.fpg_healthy, "mM"),
                        "f_cirrhosis": Q_(self.cirrhosis_map["Control"], "dimensionless"),
                        "KI__f_renal_function": Q_(self.renal_map["Normal renal function"], "dimensionless"),
                        "PODOSE_dap": Q_(dose, "mg"),
                    },
                )]
            )
        # multi dose simulation
        for dose in self.doses_multi:
            tc0 = Timecourse(
                start=0,
                end=24 * 60,  # [min]
                steps=500,
                changes={
                    **self.default_changes(),
                    "BW": Q_(self.bodyweight, "kg"),
                    "GU__f_absorption": Q_(self.fasting_map["fasted"], "dimensionless"),
                    "[KI__glc_ext]": Q_(self.fpg_healthy, "mM"),
                    "f_cirrhosis": Q_(self.cirrhosis_map["Control"], "dimensionless"),
                    "KI__f_renal_function": Q_(self.renal_map["Normal renal function"], "dimensionless"),
                    "PODOSE_dap": Q_(dose, "mg"),
                },
            )
            tc1 = Timecourse(
                start=0,
                end=24 * 60,  # [min]
                steps=500,
                changes={
                    "PODOSE_dap": Q_(dose, "mg"),
                },
            )
            tc2 = Timecourse(
                start=0,
                end=30 * 60,  # [min]
                steps=500,
                changes={
                    "PODOSE_dap": Q_(dose, "mg"),
                },
            )
            tcsims[f"po_dap_{dose}_multi"] = TimecourseSim(
                [tc0] + [tc1 for _ in range(12)] + [tc2],
                time_offset=-24*60*13
            )
            # multi dose simulation
            for dose in self.doses_uge_multi:
                tc0 = Timecourse(
                    start=0,
                    end=24 * 60,  # [min]
                    steps=500,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(self.bodyweight, "kg"),
                        "GU__f_absorption": Q_(self.fasting_map["fasted"], "dimensionless"),
                        "[KI__glc_ext]": Q_(self.fpg_healthy, "mM"),
                        "f_cirrhosis": Q_(self.cirrhosis_map["Control"], "dimensionless"),
                        "KI__f_renal_function": Q_(self.renal_map["Normal renal function"], "dimensionless"),
                        "PODOSE_dap": Q_(dose, "mg"),
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
                    end=30 * 60,  # [min]
                    steps=500,
                    changes={
                        "KI__glc_urine": Q_(0, "mmole"),  # reset UGE
                        "PODOSE_dap": Q_(dose, "mg"),
                    },
                )
                tcsims[f"po_dap_{dose}_uge_multi"] = TimecourseSim(
                    [tc0] + [tc1 for _ in range(12)] + [tc2],
                    time_offset=0
                )

        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        # fasting
        for condition in self.conditions:
            mappings[f"fm_po_dap250_{condition}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"dapagliflozin_DAP250_{condition}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self, task=f"task_po_dap_250_{condition}", xid="time", yid=f"[Cve_dap]",
                ),
                metadata=DapagliflozinMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.CAPSULE,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED if condition == "fasted" else Fasting.FED,
                ),
            )
        for condition in self.conditions:
                mappings[f"fm_po_dap250_{condition}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"dapagliflozin_cumulative amount_{condition}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self, task=f"task_po_dap_250_{condition}", xid="time", yid=f"Aurine_dap",
                    ),
                    metadata=DapagliflozinMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.CAPSULE,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED if condition == "fasted" else Fasting.FED,
                    ),
                )
        # FIXME: plot0 uge
        for dose in self.doses_uge_single:
            mappings[f"fm_po_dap_{dose}_single"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"glucose_cumulative amount_DAP{dose}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self, task=f"task_po_dap_{dose}_single", xid="time", yid=f"KI__UGE",
                ),
                metadata=DapagliflozinMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.CAPSULE,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
                ),
            )
        # multi
        for dose in self.doses_multi:
                mappings[f"fm_po_dap{dose}_MD_single"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"dapagliflozin_DAP{dose}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self, task=f"task_po_dap_{dose}_single", xid="time", yid=f"[Cve_dap]",
                    ),
                    metadata=DapagliflozinMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.CAPSULE,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                    ),
                )
        for dose in self.doses_multi:
            mappings[f"fm_po_dap{dose}_MD_multi"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"dapagliflozin_DAPM{dose}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self, task=f"task_po_dap_{dose}_multi", xid="time", yid=f"[Cve_dap]",
                ),
                metadata=DapagliflozinMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.CAPSULE,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
                ),
            )
        for dose in self.doses_uge_multi:
            mappings[f"fm_po_dap_{dose}_uge_multi"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"glucose_cumulative amount_DAPM{dose}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self, task=f"task_po_dap_{dose}_uge_multi", xid="time", yid=f"KI__UGE",
                ),
                metadata=DapagliflozinMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.CAPSULE,
                    dosing=Dosing.MULTIPLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
                ),
            )
        # console.print(mappings)
        return mappings

    def figures(self) -> Dict[str, Figure]:
        return {
            **self.figure_fasting(),
            **self.figure_uge_single(),
            **self.figure_multi(),
            **self.figure_uge_multi(),
        }

    def figure_fasting(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig_fasting",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__}",
        )
        Figure.legend_fontsize = 11
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        plots[1].set_yaxis(self.label_dap_urine, unit=self.unit_dap_urine)
        for condition in self.conditions:
            # simulation
            plots[0].add_data(
                task=f"task_po_dap_250_{condition}",
                xid="time",
                yid=f"[Cve_dap]",
                label=f"{condition} - 250 mg PO",
                color=self.fasting_colors[condition],
            )
            plots[1].add_data(
                task=f"task_po_dap_250_{condition}",
                xid="time",
                yid=f"Aurine_dap",
                label=f"{condition} - 250 mg PO",
                color=self.fasting_colors[condition],
            )
            # data
            plots[0].add_data(
                dataset=f"dapagliflozin_DAP250_{condition}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{condition} - 250 mg PO",
                color=self.fasting_colors[condition],
            )
            plots[1].add_data(
                dataset=f"dapagliflozin_cumulative amount_{condition}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{condition} - 250 mg PO",
                color=self.fasting_colors[condition],
                linestyle="",
            )
        return {fig.sid: fig}

    def figure_uge_single(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig_uge_single",
            name=f"{self.__class__.__name__}",
        )
        Figure.legend_fontsize = 7
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_uge, unit=self.unit_uge)
        for dose in self.doses_uge_single:
            plots[0].add_data(
                task=f"task_po_dap_{dose}_single",
                xid="time",
                yid=f"KI__UGE",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
            )
        for dose in self.doses_uge_single:
            plots[0].add_data(
                dataset=f"glucose_cumulative amount_DAP{dose}",
                xid="time",
                yid="mean",
                yid_sd=None,
                count="count",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
            )
        return {fig.sid: fig}

    def figure_multi(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig_multi",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__}",
        )
        Figure.legend_fontsize = 9
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        plots[1].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        plots[0].xaxis.min = -3;
        plots[0].xaxis.max = 30
        plots[1].xaxis.min = -24;
        plots[1].xaxis.max = 30

        for dose in self.doses_multi:
            # simulation
            plots[0].add_data(
                task=f"task_po_dap_{dose}_single",
                xid="time",
                yid=f"[Cve_dap]",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
            )
            plots[1].add_data(
                task=f"task_po_dap_{dose}_multi",
                xid="time",
                yid=f"[Cve_dap]",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
            )
            # data
            plots[0].add_data(
                dataset=f"dapagliflozin_DAP{dose}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
            )
            plots[1].add_data(
                dataset=f"dapagliflozin_DAPM{dose}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
            )
        return {fig.sid: fig}

    def figure_uge_multi(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig_uge_multi",
            name=f"{self.__class__.__name__}",
        )
        Figure.legend_fontsize = 9
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_uge, unit=self.unit_uge)
        for dose in self.doses_uge_multi:
            plots[0].add_data(
                task=f"task_po_dap_{dose}_uge_multi",
                xid="time",
                yid=f"KI__UGE",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
            )
        for dose in self.doses_uge_multi:
            plots[0].add_data(
                dataset=f"glucose_cumulative amount_DAPM{dose}",
                xid="time",
                yid="mean",
                yid_sd=None,
                count="count",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
                linestyle="",
            )
        return {fig.sid: fig}


if __name__ == "__main__":
    # run_experiments(Komoroski2009, output_dir=Komoroski2009.__name__)
    out = dapagliflozin.RESULTS_PATH_SIMULATION / Komoroski2009.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Komoroski2009, output_dir=out)