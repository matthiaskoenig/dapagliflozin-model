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

class FDAMB102007(DapagliflozinSimulationExperiment):
    """Simulation experiment of FDAMB102007. Renal impairment and T2DMM."""
    groups = [
        "healthy",
        "normal",
        "mild",
        "moderate",
        "severe",
    ]
    colors = {
        "healthy": DapagliflozinSimulationExperiment.renal_colors["Normal renal function"],
        "normal": DapagliflozinSimulationExperiment.renal_colors["Normal renal function"],
        "mild": DapagliflozinSimulationExperiment.renal_colors["Mild renal impairment"],
        "moderate": DapagliflozinSimulationExperiment.renal_colors["Moderate renal impairment"],
        "severe": DapagliflozinSimulationExperiment.renal_colors["Severe renal impairment"],
    }
    renal_functions = {
        "healthy": DapagliflozinSimulationExperiment.renal_map["Normal renal function"],
        "normal": DapagliflozinSimulationExperiment.renal_map["Normal renal function"],
        "mild": DapagliflozinSimulationExperiment.renal_map["Mild renal impairment"],
        "moderate": DapagliflozinSimulationExperiment.renal_map["Moderate renal impairment"],
        "severe": DapagliflozinSimulationExperiment.renal_map["Severe renal impairment"],
    }
    info = [
        ("Aurine_dap", "dapagliflozin_cumulative amount"),
        ("Aurine_d3g", "dapagliflozin 3-o-glucuronide_cumulative amount"),
        ("KI__UGE", "uge"),
    ]

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig3", "Tab13A", "Tab14A"]:
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
        for group in self.groups:
            changes = {
                # physiological changes
                "BW": Q_(self.bodyweight_default, "kg"),
                "[KI__glc_ext]": Q_(self.fpg_t2dm, "mM"),  # healthy not reported ???
                "GU__f_absorption": Q_(self.fasting_map["not reported"], "dimensionless"),
                "f_cirrhosis": Q_(self.cirrhosis_map["Control"], "dimensionless"),
                "KI__f_renal_function": Q_(self.renal_functions[group], "dimensionless"),
            }
            # single dose simulation DAP50
            tcsims[f"po_dap50_{group}"] = TimecourseSim(
                Timecourse(
                    start=0,
                    end=65 * 60,  # [min]
                    steps=500,
                    changes={
                        **self.default_changes(),
                        **changes,
                        # dose (IVDOSE, PODOSE)
                        "PODOSE_dap": Q_(50, "mg"),
                    },
                )
            )
            # multi dose simulation DAP20
            tc0 = Timecourse(
                start=0,
                end=24 * 60,  # [min]
                steps=500,
                changes={
                    **self.default_changes(),
                    **changes,
                    # dose
                    "PODOSE_dap": Q_(20, "mg"),
                },
            )
            tc1 = Timecourse(
                start=0,
                end=24 * 60,  # [min]
                steps=500,
                changes={
                    "KI__glc_urine": Q_(0, "mmole"),  # reset urinary data
                    "PODOSE_dap": Q_(20, "mg"),
                },
            )
            tc2 = Timecourse(
                start=0,
                end=48 * 60,  # [min]
                steps=500,
                changes={
                    "KI__glc_urine": Q_(0, "mmole"),  # reset urinary data
                    "PODOSE_dap": Q_(20, "mg"),
                },
            )
            tcsims[f"po_dap20_{group}"] = TimecourseSim(
                [tc0] + [tc1 for _ in range(5)] + [tc2],
                time_offset=-6*24*60,
            )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        # Fig3
        for group in self.groups:
            if group == "healthy":
                continue
            health = Health.T2DM if group == "normal" else Health.RENAL_IMPAIRMENT_T2DM
            mappings[f"fm_dap20_{group}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"uge_DAP20_{group}",
                    xid="time",
                    yid="mean",
                    yid_sd=None,
                    count="count",
                ),
                observable=FitData(
                    self, task=f"task_po_dap20_{group}", xid="time", yid="KI__UGE",
                ),
                metadata=DapagliflozinMappingMetaData(
                    tissue=Tissue.URINE,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.MULTIPLE,
                    health=health,
                    fasting=Fasting.NR,
                ),
            )
        # Tab13A, Tab14A
        for group in self.groups:
            health = Health.HEALTHY if group == "healthy" else Health.T2DM if group == "normal" else Health.RENAL_IMPAIRMENT_T2DM
            # simulation
            for sid, substance in [
                ("Aurine_dap", "dapagliflozin_cumulative amount"),
                ("Aurine_d3g", "dapagliflozin 3-o-glucuronide_cumulative amount"),
            ]:
                mappings[f"fm_dap50_{group}_{substance}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{substance}_DAP50_{group}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self, task=f"task_po_dap50_{group}", xid="time", yid=sid,
                    ),
                    metadata=DapagliflozinMappingMetaData(
                        tissue=Tissue.URINE,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=health,
                        fasting=Fasting.NR,
                    ),
                )
        return mappings

    def figures(self) -> Dict[str, Figure]:
        return {
            **self.figure_fig3(),
            **self.figure_tab13_14(),
        }

    def figure_fig3(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig3",
            name=f"{self.__class__.__name__}"
        )
        Figure.legend_fontsize = 10
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_uge, unit=self.unit_uge)
        for group in self.groups:
            if group == "healthy":
                continue
            dose_label = f"{group} - 20 mg PO"
            # simulation
            plots[0].add_data(
                task=f"task_po_dap20_{group}",
                xid="time",
                yid="KI__UGE",
                label=f"{dose_label}",
                color=self.colors[group],
            )
            # data
            plots[0].add_data(
                dataset=f"uge_DAP20_{group}",
                xid="time",
                yid="mean",
                yid_sd=None,
                count="count",
                label=f"{dose_label}",
                color=self.colors[group],
            )

        return {
            fig.sid: fig,
        }

    def figure_tab13_14(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Tab13_14",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__}"
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_dap_urine, unit=self.unit_dap_urine)
        plots[1].set_yaxis(self.label_d3g_urine, unit=self.unit_d3g_urine)

        for group in self.groups:
            dose_label_dap = f"{group} - 50 mg PO"
            dose_label_d3g = f"{group} - 50 mg PO"
            # simulation
            plots[0].add_data(
                task=f"task_po_dap50_{group}",
                xid="time",
                yid="Aurine_dap",
                label=f"{dose_label_dap}",
                color=self.colors[group],
            )
            plots[1].add_data(
                task=f"task_po_dap50_{group}",
                xid="time",
                yid="Aurine_d3g",
                label=f"{dose_label_d3g}",
                color=self.colors[group],
            )
            # data
            plots[0].add_data(
                dataset=f"dapagliflozin_cumulative amount_DAP50_{group}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{dose_label_dap}",
                linestyle="",
                color=self.colors[group],
            )
            plots[1].add_data(
                dataset=f"dapagliflozin 3-o-glucuronide_cumulative amount_DAP50_{group}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{dose_label_d3g}",
                linestyle="",
                color=self.colors[group],
            )

        return {
            fig.sid: fig,
        }

if __name__ == "__main__":
    # run_experiments(FDAMB102007, output_dir=FDAMB102007.__name__)
    out = dapagliflozin.RESULTS_PATH_SIMULATION / FDAMB102007.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(FDAMB102007, output_dir=out)