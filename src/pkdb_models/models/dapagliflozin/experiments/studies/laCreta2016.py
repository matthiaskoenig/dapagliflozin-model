"""Simulation experiment of the study LaCreta2016."""

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


class LaCreta2016(DapagliflozinSimulationExperiment):
    """Simulation experiment of LaCreta2016.

    Two randomized crossover bioequivalence studies in healthy subjects: single
    doses of 10 mg (study 1: NH fasted, HS fasted, HS fed) and 2.5 mg (study 2:
    NH and HS, fasted and fed) of non-heat-stressed (NH) and heat-stressed (HS)
    dapagliflozin tablets; fed is a high-fat high-calorie breakfast 30 min
    before the dose. The plasma of Fig 2 (mean + SD); the times of the 10 mg
    panel before 4 hr are the sampling times (the curation read them one
    sampling time late, the 0 of the x axis is offset in panel A).
    """

    conditions: ClassVar[list[str]] = ["fasted", "fed"]
    groups: ClassVar[list[str]] = ["NH", "HS"]
    doses: ClassVar[list[float]] = [2.5, 10]
    #: the body weights of the studies: 10 mg study 1, 2.5 mg study 2 [kg]
    bodyweights: ClassVar[dict[float, float]] = {
        2.5: 70.1,
        10: 73.5,
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                # unit conversion to mole/l
                if label.startswith("dapagliflozin_"):
                    dset.unit_conversion("mean", 1 / self.Mr.dap)
                dsets[self.to_sid(label)] = dset
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for condition in self.conditions:
            for dose in self.doses:
                simulations[self.to_sid(f"po_dap{dose}_{condition}")] = Simulation(
                    time_unit="hr",
                    end=50,
                    preinit_changes={
                        **self.default_changes(),
                        # physiological changes
                        "BW": Q(self.bodyweights[dose], "kg"),
                        "[KI__fpg]": Q(self.fpg_healthy, "mM"),
                        "GU__fed": Q(self.fasting_map[condition], "dimensionless"),
                        "f_cirrhosis": Q(
                            self.cirrhosis_map["Control"], "dimensionless"
                        ),
                        "KI__f_renal_function": Q(
                            self.renal_map["Normal renal function"], "dimensionless"
                        ),
                        # dose (IVDOSE, PODOSE)
                        "PODOSE_dap": Q(dose, "mg"),
                    },
                )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for condition in self.conditions:
            for group in self.groups:
                for dose in self.doses:
                    if dose == 10 and condition == "fed" and group == "NH":
                        continue
                    mappings[self.to_sid(f"fm_po_dap{dose}_{group}_{condition}")] = (
                        FitMapping(
                            self,
                            reference=FitData(
                                self,
                                dataset=self.to_sid(
                                    f"dapagliflozin_DAP{dose}_{group}_{condition}"
                                ),
                                xid="time",
                                yid="mean",
                                yid_sd="mean_sd",
                                count="count",
                            ),
                            observable=FitData(
                                self,
                                task=self.to_sid(f"task_po_dap{dose}_{condition}"),
                                xid="time",
                                yid="[Cve_dap]",
                            ),
                            metadata=PKPDMappingMetaData(
                                tissue=Tissue.PLASMA,
                                route=Route.PO,
                                application_form=ApplicationForm.TABLET,
                                dosing=Dosing.SINGLE,
                                health=Health.HEALTHY,
                                fasting=Fasting.FASTED
                                if condition == "fasted"
                                else Fasting.FED,
                            ),
                        )
                    )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        plots[1].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        for k, dose in enumerate(self.doses):
            for condition in self.conditions:
                # simulation
                plots[k].add_data(
                    task=self.to_sid(f"task_po_dap{dose}_{condition}"),
                    xid="time",
                    yid="[Cve_dap]",
                    label=f"{condition} - {dose} mg PO",
                    color=self.fasting_colors[condition],
                )
                for group in self.groups:
                    if dose == 10 and condition == "fed" and group == "NH":
                        continue
                    # data
                    plots[k].add_data(
                        dataset=self.to_sid(
                            f"dapagliflozin_DAP{dose}_{group}_{condition}"
                        ),
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                        label=f"{condition}, {group} - {dose} mg PO",
                        color=self.fasting_colors[condition],
                        marker="s" if group == "NH" else "D",
                    )
        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(LaCreta2016, output_dir=LaCreta2016.__name__)
