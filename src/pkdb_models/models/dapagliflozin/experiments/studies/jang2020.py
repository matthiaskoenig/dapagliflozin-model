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

class Jang2020(DapagliflozinSimulationExperiment):
    """Simulation experiment of Jang2020."""

    interventions = ["DAP10", "DAP10, LOB05"]
    bodyweight = 71.45,  # [kg]

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
        tc0 = Timecourse(
            start=0,
            end=24 * 60,  # [min]
            steps=500,
            changes={
                **self.default_changes(),
                # physiological changes
                "BW": Q_(self.bodyweight, "kg"),
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
            end=24 * 60,  # [min]
            steps=500,
            changes={
                "PODOSE_dap": Q_(10, "mg"),
            },
        )
        tc2 = Timecourse(
            start=0,
            end=48 * 60,  # [min]
            steps=500,
            changes={
                "PODOSE_dap": Q_(10, "mg"),
            },
        )
        tcsims[f"po_dap10"] = TimecourseSim(
            [tc0] + [tc1 for _ in range(3)] + [tc2],
            time_offset=-4 * 24 * 60,
        )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for intervention in self.interventions:
            mappings[f"fm_dap10_{intervention}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"dapagliflozin_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self, task=f"task_po_dap10", xid="time", yid=f"[Cve_dap]",
                ),
                metadata=DapagliflozinMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.MULTIPLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
                    coadministration=Coadministration.LOBEGLITAZONE if "LOB" in intervention else Coadministration.NONE
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
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        plots[0].xaxis.min = -24
        plots[0].xaxis.max = 50

        label = "10 mg PO"
        # simulation
        plots[0].add_data(
            task=f"task_po_dap10",
            xid="time",
            yid=f"[Cve_dap]",
            label=label,
            color="black",
        )

        # data
        plots[0].add_data(
            dataset="dapagliflozin_DAP10",
            xid="time",
            yid="mean",
            yid_sd="mean_sd",
            count="count",
            label=label,
            color="black",
        )
        plots[0].add_data(
            dataset="dapagliflozin_DAP10, LOB05",
            xid="time",
            yid="mean",
            yid_sd="mean_sd",
            count="count",
            label="LOB05 + 10 mg PO",
            color="tab:blue",
        )

        return {fig.sid: fig}


if __name__ == "__main__":
    # run_experiments(Jang2020, output_dir=Jang2020.__name__)
    out = dapagliflozin.RESULTS_PATH_SIMULATION / Jang2020.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Jang2020, output_dir=out)