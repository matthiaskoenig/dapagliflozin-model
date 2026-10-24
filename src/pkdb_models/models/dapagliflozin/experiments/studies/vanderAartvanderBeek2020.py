"""Simulation experiment of the study vanderAartvanderBeek2020."""

from typing import ClassVar, override

from pkpdsim.metadata import (
    ApplicationForm,
    Dosing,
    Fasting,
    Health,
    PKPDMappingMetaData,
    Route,
    Tissue,
)
from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.fit import FitData, FitMapping
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Simulation

from pkdb_models.models.dapagliflozin import MODEL
from pkdb_models.models.dapagliflozin.experiments.base_experiment import (
    DapagliflozinSimulationExperiment,
)


class vanderAartvanderBeek2020(DapagliflozinSimulationExperiment):
    """Simulation experiment of vanderAartvanderBeek2020.

    10 mg dapagliflozin once daily for 6 weeks in patients with non-diabetic
    chronic kidney disease and macroalbuminuria: the individual plasma profiles
    of six patients at steady state, predose and 30 to 240 min after the dose
    (Fig 3, LC-MS/MS method paper). Body weight, eGFR, food and comedication
    are not reported: the default body weight and the renal function of the
    cohort of the DIAMOND trial (same group, design and dose; mean eGFR 58
    ml/min/1.73 m2, Cherney2020), which is the likely source of the samples.
    """

    #: daily doses until the profile, 6 weeks
    n_doses = 42
    #: eGFR of the DIAMOND cohort [ml/min/1.73 m2]
    egfr = 58.0

    subjects_markers: ClassVar[dict[str, str]] = {
        "S1": "o",  # circle
        "S2": "s",  # square
        "S3": "D",  # diamond
        "S4": "^",  # triangle up
        "S5": "v",  # triangle down
        "S6": "<",  # triangle left
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig3"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                # unit conversion
                if label.startswith("dapagliflozin_"):
                    dset.unit_conversion("value", 1 / self.Mr.dap)
                dsets[self.to_sid(label)] = dset
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        simulations["po_dap10"] = Simulation(
            time_unit="hr",
            end=(self.n_doses - 1) * 24 + 6,
            preinit_changes={
                **self.default_changes(),
                # physiological changes
                "BW": Q(self.bodyweight_default, "kg"),
                "[KI__fpg]": Q(self.fpg_healthy, "mM"),
                "GU__fed": Q(self.fasting_map["not reported"], "dimensionless"),
                "f_cirrhosis": Q(self.cirrhosis_map["Control"], "dimensionless"),
                "KI__f_renal_function": Q(
                    self.egfr / self.gfr_healthy, "dimensionless"
                ),
                # dose (IVDOSE, PODOSE)
                "PODOSE_dap": Q(10, "mg"),
            },
            changes=[
                self.add_doses(
                    [24 * k for k in range(1, self.n_doses)],
                    {"PODOSE_dap": Q(10, "mg")},
                ),
            ],
            time_shift=-(self.n_doses - 1) * 24,
        )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for subject in self.subjects_markers:
            mappings[self.to_sid(f"fm_dap10_{subject}")] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=self.to_sid(f"dapagliflozin_DAP_{subject}"),
                    xid="time",
                    yid="value",
                    yid_sd=None,
                    count="count",
                ),
                observable=FitData(
                    self,
                    task="task_po_dap10",
                    xid="time",
                    yid="[Cve_dap]",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.MULTIPLE,
                    health=Health.RENAL_IMPAIRMENT,
                    fasting=Fasting.NR,
                ),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig3",
            name=f"{self.__class__.__name__} (CKD)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, min=-0.5, max=6),
            legend=True,
        )
        plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        # simulation
        plots[0].add_data(
            task="task_po_dap10",
            xid="time",
            yid="[Cve_dap]",
            label="10 mg PO, day 42",
            color="black",
        )
        # data
        for subject, marker in self.subjects_markers.items():
            plots[0].add_data(
                dataset=self.to_sid(f"dapagliflozin_DAP_{subject}"),
                xid="time",
                yid="value",
                yid_sd=None,
                count="count",
                color="black",
                label=f"10 mg PO ({subject})",
                marker=marker,
            )
        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(
        vanderAartvanderBeek2020, output_dir=vanderAartvanderBeek2020.__name__
    )
