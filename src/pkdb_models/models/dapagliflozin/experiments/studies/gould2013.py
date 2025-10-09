from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console

from pkdb_models.models.dapagliflozin.experiments.base_experiment import (
    DapagliflozinSimulationExperiment,
)
from pkdb_models.models.dapagliflozin.experiments.metadata import Tissue, Route, Dosing, ApplicationForm, Health, \
    Fasting, DapagliflozinMappingMetaData

from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.dapagliflozin.helpers import run_experiments


class Gould2013(DapagliflozinSimulationExperiment):
    """Simulation experiment of Gould2013."""

    doses = [0, 0.001, 0.01, 0.1, 0.3, 1, 2.5, 5, 10, 20, 50, 100, 250, 500]
    colors = {
        0: "black",
        0.001: "tab:blue",
        0.01: "tab:orange",
        0.1: "tab:green",
        0.3: "tab:red",
        1: "tab:purple",
        2.5: "tab:brown",
        5: "tab:pink",
        10: "tab:grey",
        20: "tab:olive",
        50: "tab:cyan",
        100: "lightcoral",
        250: "gold",
        500: "lightgreen"
    }
    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1"]:
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
        for dose in self.doses:
            tcsims[f"po_dap{dose}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=25 * 60,  # [min]
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
                        "PODOSE_dap": Q_(dose, "mg"),
                    },
                )]
            )

            return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:

        mappings = {}
        for dose in self.doses:
            applicationform = ApplicationForm.TABLET if dose == [0, 1, 2.5, 5, 10, 20, 50, 100, 250, 500] else ApplicationForm.SOLUTION
            mappings[f"fm_dap{dose}"] = FitMapping(
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
                    self, task=f"task_po_dap{dose}", xid="time", yid=f"KI__UGE",
                ),
                metadata=DapagliflozinMappingMetaData(
                    tissue=Tissue.URINE,
                    route=Route.PO,
                    application_form=applicationform,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
                ),
            )
        # console.print(mappings)
        return mappings

    def figures(self) -> Dict[str, Figure]:

        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_rows=1,
            num_cols=1,
            name=f"{self.__class__.__name__}",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_uge, unit=self.unit_uge)

        # simulation
        for dose in self.doses:
            plots[0].add_data(
                task=f"task_po_dap{dose}",
                xid="time",
                yid=f"KI__UGE",
                label=f"Sim {dose} mg",
                color=self.colors[dose],
            )
            # data
            plots[0].add_data(
                dataset=f"glucose_cumulative amount_DAP{dose}",
                xid="time",
                yid="mean",
                yid_sd=None,
                count="count",
                label=f"UGE {dose} mg",
                color=self.colors[dose],
            )

        return {
            fig.sid: fig,
        }


if __name__ == "__main__":
    run_experiments(Gould2013, output_dir=Gould2013.__name__)
