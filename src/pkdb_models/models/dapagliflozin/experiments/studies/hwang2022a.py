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

class Hwang2022a(DapagliflozinSimulationExperiment):
    """
    Simulation experiment of Hwang2022a.
    Data from Fig1 is also reported in Fig2 for C3435T_TT. C3435T_CC+CT from Fig1 not used (included as CC and CT data
    from Fig2). Fig1 data not used.
    """

    bodyweights = {
        "C1236T_CC": 71.03,
        "C1236T_CT": 71.79,
        "C1236T_TT": 72.44,
        "G2677TA_GG": 71.18,
        "G2677TA_GT": 71,
        "G2677TA_GA": 73.61,
        "G2677TA_TT": 71.43,
        "G2677TA_TA": 72.78,
        "G2677TA_AA": 67.05,
        "C3435T_CC": 71.13,
        "C3435T_CT": 73.44,
        "C3435T_TT": 68.09,
    }  # [kg]

    markers = {
        "C1236T_CC": "o",  # circle
        "C1236T_CT": "s",  # square
        "C1236T_TT": "D",  # diamond
        "G2677TA_GG": "^",  # triangle up
        "G2677TA_GT": "v",  # triangle down
        "G2677TA_GA": "<",  # triangle left
        "G2677TA_TT": ">",  # triangle right
        "G2677TA_TA": "+",  # plus
        "G2677TA_AA": "x",  # cross
        "C3435T_CC": "*",  # star
        "C3435T_CT": "|",  # vertical line
        "C3435T_TT": "_",  # horizontal line
    }

    groups = list(bodyweights.keys())

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2"]:
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
        for group in self.groups:
            tcsims[f"po_dap10_{group}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=50 * 60,  # [min]
                    steps=1000,
                    changes={
                        # assuming no changes due to ABCB1 genotype
                        **self.default_changes(),
                        # physiological changes
                        "BW": Q_(self.bodyweights[group], "kg"),
                        "[KI__glc_ext]": Q_(self.fpg_healthy, "mM"),
                        "GU__f_absorption": Q_(self.fasting_map["fasted"], "dimensionless"),
                        "f_cirrhosis": Q_(self.cirrhosis_map["Control"], "dimensionless"),
                        "KI__f_renal_function": Q_(self.renal_map["Normal renal function"], "dimensionless"),
                        # dose (IVDOSE, PODOSE)
                        "PODOSE_dap": Q_(10, "mg"),
                    },
                )]
            )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for group in self.groups:
            mappings[f"fm_dap10_{group}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"dapagliflozin_DAP10_{group}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self, task=f"task_po_dap10_{group}", xid="time", yid=f"[Cve_dap]",
                ),
                metadata=DapagliflozinMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
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
            sid="Fig2",
            num_rows=1,
            num_cols=1,
            name=f"{self.__class__.__name__}",
        )
        Figure.legend_fontsize = 10
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)

        for group in self.groups:
            # simulation
            plots[0].add_data(
                task=f"task_po_dap10_{group}",
                xid="time",
                yid=f"[Cve_dap]",
                label=None,
                color="black",
                linewidths=0.5,
            )
            # data
            plots[0].add_data(
                dataset=f"dapagliflozin_DAP10_{group}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=group,
                color="tab:blue",
                linewidths=0.5,
                markersize=9,
                markeredgewidth=0.5,
                marker=self.markers[group],
            )

        return {
            fig.sid: fig,
        }


if __name__ == "__main__":
    # run_experiments(Hwang2022a, output_dir=Hwang2022a.__name__)
    out = dapagliflozin.RESULTS_PATH_SIMULATION / Hwang2022a.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Hwang2022a, output_dir=out)