"""Simulation experiment of the study Yang2013."""

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
from sbmlsim.simulation import Change, Simulation

from pkdb_models.models.dapagliflozin import MODEL
from pkdb_models.models.dapagliflozin.experiments.base_experiment import (
    DapagliflozinSimulationExperiment,
)


class Yang2013(DapagliflozinSimulationExperiment):
    """Simulation experiment of Yang2013.

    Single and multiple doses of 5 mg and 10 mg dapagliflozin in healthy Chinese
    subjects (two studies, 14 subjects each): a single dose on day 1, then once
    daily on days 5 to 10, fasted. Plasma of dapagliflozin and D3G after day 1
    and day 10 (Figure A-D, mean (SD)), urinary glucose 0-24 and 0-96 hr after
    day 1 and day 10 (Tab III). The day 1 series have no predose point (the
    predose values of day 10 were copied into them in the curation).
    """

    doses_single: ClassVar[list[float]] = [5, 10]
    doses_multi: ClassVar[list[float]] = [5, 10]
    info: ClassVar[list[tuple[str, str]]] = [
        ("[Cve_dap]", "dapagliflozin"),
        ("[Cve_d3g]", "dapagliflozin 3-o-glucuronide"),
        ("KI__UGE", "glucose_cumulative amount"),
    ]
    bodyweights: ClassVar[dict[float, float]] = {
        5: 62.8,
        10: 62.2,
    }  # [kg]

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Tab3"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                if label.startswith("glucose_cumulative"):
                    # the 0 hr points are not data of the paper, and the 0 hr
                    # of day 10 is the reset of the urine at the dose
                    df_label = df_label[df_label["time"] > 0]
                dset = DataSet.from_df(df_label, self.ureg)
                # unit conversion to mole/l
                if label.startswith("dapagliflozin_"):
                    dset.unit_conversion("mean", 1 / self.Mr.dap)
                elif label.startswith("dapagliflozin 3-o-glucuronide"):
                    dset.unit_conversion("mean", 1 / self.Mr.d3g)
                dsets[self.to_sid(label)] = dset
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        # single dose simulation
        for dose in self.doses_single:
            simulations[self.to_sid(f"po_dap{dose}_single")] = Simulation(
                time_unit="hr",
                end=100,
                preinit_changes={
                    **self.default_changes(),
                    # physiological changes
                    "BW": Q(self.bodyweights[dose], "kg"),
                    "[KI__fpg]": Q(self.fpg_healthy, "mM"),
                    "GU__fed": Q(self.fasting_map["fasted"], "dimensionless"),
                    "f_cirrhosis": Q(self.cirrhosis_map["Control"], "dimensionless"),
                    "KI__f_renal_function": Q(
                        self.renal_map["Normal renal function"], "dimensionless"
                    ),
                    # dose
                    "PODOSE_dap": Q(dose, "mg"),
                },
            )
        # multi dose simulation: the single dose of day 1, then daily doses on
        # days 5 to 10 (TabInterventions); the data is of day 10 at time 0
        for dose in self.doses_multi:
            simulations[self.to_sid(f"po_dap{dose}_multi")] = Simulation(
                time_unit="hr",
                end=9 * 24 + 100,
                preinit_changes={
                    **self.default_changes(),
                    # physiological changes
                    "BW": Q(self.bodyweights[dose], "kg"),
                    "[KI__fpg]": Q(self.fpg_healthy, "mM"),
                    "GU__fed": Q(self.fasting_map["fasted"], "dimensionless"),
                    "f_cirrhosis": Q(self.cirrhosis_map["Control"], "dimensionless"),
                    "KI__f_renal_function": Q(
                        self.renal_map["Normal renal function"], "dimensionless"
                    ),
                    # dose
                    "PODOSE_dap": Q(dose, "mg"),
                },
                changes=[
                    Change(
                        [24 * k for k in range(4, 10)],
                        {"KI__glc_urine": Q(0, "mmole")},  # reset UGE
                    ),
                    self.add_doses(
                        [24 * k for k in range(4, 10)], {"PODOSE_dap": Q(dose, "mg")}
                    ),
                ],
                time_shift=-9 * 24,
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        # Fig1
        for sid, prefix in self.info:
            for dose in self.doses_single:
                tissue = Tissue.URINE if "cumulative" in prefix else Tissue.PLASMA
                mappings[self.to_sid(f"fm_po_dap{dose}_{sid.strip('[]')}_single")] = (
                    FitMapping(
                        self,
                        reference=FitData(
                            self,
                            dataset=self.to_sid(f"{prefix}_SD_DAP{dose}"),
                            xid="time",
                            yid="mean",
                            yid_sd="mean_sd",
                            count="count",
                        ),
                        observable=FitData(
                            self,
                            task=self.to_sid(f"task_po_dap{dose}_single"),
                            xid="time",
                            yid=sid,
                        ),
                        metadata=PKPDMappingMetaData(
                            tissue=tissue,
                            route=Route.PO,
                            application_form=ApplicationForm.TABLET,
                            dosing=Dosing.SINGLE,
                            health=Health.HEALTHY,
                            fasting=Fasting.FASTED,
                        ),
                    )
                )
        # Fig2
        for sid, prefix in self.info:
            for dose in self.doses_multi:
                tissue = Tissue.URINE if "cumulative" in prefix else Tissue.PLASMA
                mappings[self.to_sid(f"fm_po_dap{dose}_{sid.strip('[]')}_multi")] = (
                    FitMapping(
                        self,
                        reference=FitData(
                            self,
                            dataset=self.to_sid(f"{prefix}_MD_DAP{dose}"),
                            xid="time",
                            yid="mean",
                            yid_sd="mean_sd",
                            count="count",
                        ),
                        observable=FitData(
                            self,
                            task=self.to_sid(f"task_po_dap{dose}_multi"),
                            xid="time",
                            yid=sid,
                        ),
                        metadata=PKPDMappingMetaData(
                            tissue=tissue,
                            route=Route.PO,
                            application_form=ApplicationForm.TABLET,
                            dosing=Dosing.MULTIPLE,
                            health=Health.HEALTHY,
                            fasting=Fasting.FASTED,
                        ),
                    )
                )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        return {**self.figure_plasma(), **self.figure_uge()}

    def figure_plasma(self) -> dict[str, Figure]:
        """Dapagliflozin and its glucuronide of day 1 and day 10 (Figure A-D)."""
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_rows=2,
            num_cols=2,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, min=-1, max=25),
            legend=True,
        )
        # A: 5 mg dapagliflozin, B: 5 mg D3G, C: 10 mg dapagliflozin, D: 10 mg D3G
        for k, (dose, (sid, prefix)) in enumerate(
            (dose, info) for dose in self.doses_single for info in self.info[:2]
        ):
            plot = plots[k]
            label = (
                self.label_dap_plasma if sid == "[Cve_dap]" else self.label_d3g_plasma
            )
            # down to about the LLOQ of 0.30 ng/ml
            plot.set_yaxis(label, unit=self.unit_dap, scale="log", min=5e-4)
            for dosing, color, marker in [
                ("single", "darkgrey", "o"),
                ("multi", self.dose_colors[dose], "s"),
            ]:
                day = "day 1" if dosing == "single" else "day 10"
                plot.add_data(
                    task=self.to_sid(f"task_po_dap{dose}_{dosing}"),
                    xid="time",
                    yid=sid,
                    label=f"{dose} mg PO, {day}",
                    color=color,
                )
                plot.add_data(
                    dataset=self.to_sid(
                        f"{prefix}_{'SD' if dosing == 'single' else 'MD'}_DAP{dose}"
                    ),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"{dose} mg PO, {day}",
                    color=color,
                    marker=marker,
                )
        return self.figures_by_sid(fig)

    def figure_uge(self) -> dict[str, Figure]:
        """Urinary glucose excretion after day 1 and day 10 (Tab III)."""
        fig = Figure(
            experiment=self,
            sid="Tab3",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, min=-2, max=100),
            legend=True,
        )
        for k, dosing in enumerate(["single", "multi"]):
            day = "day 1" if dosing == "single" else "day 10"
            plots[k].set_title(f"{day}")
            plots[k].set_yaxis(self.label_uge, unit=self.unit_uge)
            for dose in self.doses_single:
                plots[k].add_data(
                    task=self.to_sid(f"task_po_dap{dose}_{dosing}"),
                    xid="time",
                    yid="KI__UGE",
                    label=f"{dose} mg PO, {day}",
                    color=self.dose_colors[dose],
                )
                plots[k].add_data(
                    dataset=self.to_sid(
                        "glucose_cumulative amount_"
                        f"{'SD' if dosing == 'single' else 'MD'}_DAP{dose}"
                    ),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"{dose} mg PO, {day}",
                    color=self.dose_colors[dose],
                    linestyle="",
                )
        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Yang2013, output_dir=Yang2013.__name__)
