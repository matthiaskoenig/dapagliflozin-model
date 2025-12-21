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


class Boulton2013(DapagliflozinSimulationExperiment):
    """Simulation experiment of Boulton2013."""

    simulation_keys = ["po_iv_dap10", "iv_dap008"]
    interventions = ["DAP10", "C14DAP008"]
    colors = ["black", "tab:blue"]

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
        # po followed by iv
        tc0 = Timecourse(
            start=0,
            end=60,  # [min]
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
                "PODOSE_dap": Q_(10, "mg"),
            },
        )
        tc1 = Timecourse(
            start=0,
            end=49 * 60,  # [min]
            steps=500,
            changes={
                "IVDOSE_dap": Q_(80, "µg"),  # 1.6 % labeled; 98.4% unlabeled
            },
        )
        tcsims[f"po_iv_dap10"] = TimecourseSim(
            [tc0, tc1],
        )
        # only iv at 1 hr
        tc0_iv = Timecourse(
            start=0,
            end=60,  # [min]
            steps=500,
            changes={
                **self.default_changes(),
                # physiological changes
                "BW": Q_(self.bodyweight_default, "kg"),
                "[KI__glc_ext]": Q_(self.fpg_healthy, "mM"),
                "GU__f_absorption": Q_(self.fasting_map["fasted"], "dimensionless"),
                "f_cirrhosis": Q_(self.cirrhosis_map["Control"], "dimensionless"),
                "KI__f_renal_function": Q_(self.renal_map["Normal renal function"], "dimensionless"),
            },
        )
        tc1_iv = Timecourse(
            start=0,
            end=49 * 60,  # [min]
            steps=500,
            changes={
                "IVDOSE_dap": Q_(80, "µg") #Q_(80*0.016, "µg"),  # 1.6 % labeled; 98.4% unlabeled
            },
        )
        tcsims[f"iv_dap008"] = TimecourseSim(
            [tc0_iv, tc1_iv]
        )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for ks, sim_key in enumerate(self.simulation_keys):
            intervention = self.interventions[ks]
            mappings[f"fm_{intervention}_dapagliflozin"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"dapagliflozin_{self.interventions[ks]}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self, task=f"task_{sim_key}", xid="time", yid=f"[Cve_dap]",
                ),
                metadata=DapagliflozinMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO if intervention == "DAP10" else Route.IV,
                    application_form=ApplicationForm.MIXED if intervention == "DAP10" else ApplicationForm.SOLUTION,
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
            num_cols=2,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        plots[1].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)

        for ks, sim_key in enumerate(self.simulation_keys):
            dose_label = "10 mg PO" if sim_key == "po_iv_dap10" else "80 µg IV"
            # simulation
            plots[ks].add_data(
                task=f"task_{sim_key}",
                xid="time",
                yid=f"[Cve_dap]",
                label=f"{dose_label}",
                color=self.colors[ks],
            )
            # data
            plots[ks].add_data(
                dataset=f"dapagliflozin_{self.interventions[ks]}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{dose_label}",
                color=self.colors[ks],
            )

        return {fig.sid: fig,}


if __name__ == "__main__":
    out = dapagliflozin.RESULTS_PATH_SIMULATION / Boulton2013.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Boulton2013, output_dir=out)