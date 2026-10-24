"""Simulation experiment of the study FDAMB102002."""

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


class FDAMB102002(DapagliflozinSimulationExperiment):
    """Simulation experiment of FDAMB102002.

    Study MB102002 of the clinical pharmacology review of NDA 202293, the
    multiple ascending dose study of Komoroski2009: healthy subjects (n=6 per
    dose, body weight not reported) took 2.5, 10, 20, 50 or 100 mg
    dapagliflozin as a capsule (Komoroski2009) once daily for 14 days, fasted
    on the days of the pharmacokinetics. Plasma on days 1, 7 and 14 (Fig12),
    the cumulative urinary glucose excretion over 0-24 hr on days 1 and 14
    (Fig17, without SD). The urinary glucose is collected per day and reset
    with every dose; a data point at the time of a reset is compared with the
    state before it. The fasting plasma glucose is the healthy value: the
    subjects had meals during the urine collections (a limitation of all
    healthy UGE studies). The predose sample of day 1 (10 mg, below the LLOQ,
    drawn at the bottom of the log axis) is removed from the data; the points
    at 0.2 hr of day 1 at the floor of the axis are kept (lag of the capsule).
    """

    days: ClassVar[dict[int, float]] = {1: 0.0, 7: 144.0, 14: 312.0}

    doses: ClassVar[list[float]] = [2.5, 10, 20, 50, 100]

    info: ClassVar[dict[str, str]] = {
        "[Cve_dap]": "dapagliflozin",
        "KI__UGE": "uge",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        # Fig12: dap
        # Fig17: UGE
        for fig_id in ["Fig12", "Fig17"]:
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
        for dose in self.doses:
            simulations[self.to_sid(f"po_MAD{dose}")] = Simulation(
                time_unit="hr",
                end=13 * 24 + 25,
                preinit_changes={
                    **self.default_changes(),
                    # physiological changes
                    "BW": Q(self.bodyweight_default, "kg"),
                    "[KI__fpg]": Q(self.fpg_healthy, "mM"),
                    "GU__fed": Q(self.fasting_map["fasted"], "dimensionless"),
                    "f_cirrhosis": Q(self.cirrhosis_map["Control"], "dimensionless"),
                    "KI__f_renal_function": Q(
                        self.renal_map["Normal renal function"], "dimensionless"
                    ),
                    # dose (IVDOSE, PODOSE)
                    "PODOSE_dap": Q(dose, "mg"),
                },
                changes=[
                    Change(
                        [24 * k for k in range(1, 14)],
                        {"KI__glc_urine": Q(0, "mmole")},  # reset UGE
                    ),
                    self.add_doses(
                        [24 * k for k in range(1, 14)], {"PODOSE_dap": Q(dose, "mg")}
                    ),
                ],
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for dose in self.doses:
            for sid, name in self.info.items():
                mappings[self.to_sid(f"fm_dap{dose}_{name}")] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=self.to_sid(f"{name}_MAD{dose}"),
                        xid="time",
                        yid="mean",
                        # the urinary glucose has no SD
                        yid_sd=None if name == "uge" else "mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=self.to_sid(f"task_po_MAD{dose}"),
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.URINE if "UGE" in sid else Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.CAPSULE,
                        dosing=Dosing.MULTIPLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                    ),
                )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self._figure_days(
                sid="Fig12", key="[Cve_dap]", days=[1, 7, 14], log=True
            ),
            **self._figure_days(sid="Fig17", key="KI__UGE", days=[1, 14], log=False),
        }

    def _figure_days(
        self, sid: str, key: str, days: list[int], log: bool
    ) -> dict[str, Figure]:
        """Plasma concentration (Fig12) or urinary glucose (Fig17) per study day."""
        name = self.info[key]
        fig = Figure(
            experiment=self,
            sid=sid,
            num_rows=1,
            num_cols=len(days),
            name=f"{self.__class__.__name__} (Healthy)",
            # one legend of the row, right of the grid
            settings={"legend_scope": "row"},
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        for plot, day in zip(plots, days, strict=True):
            start = self.days[day]
            plot.set_title(f"day {day}")
            plot.set_xaxis(
                self.label_time, unit=self.unit_time, min=start - 1, max=start + 25
            )
            if key == "KI__UGE":
                plot.set_yaxis(self.label_uge, unit=self.unit_uge)
            else:
                plot.set_yaxis(
                    self.label_dap_plasma,
                    unit=self.unit_dap,
                    scale="log" if log else "linear",
                )
            for dose in self.doses:
                dose_label = f"{dose} mg PO"
                # simulation
                plot.add_data(
                    task=self.to_sid(f"task_po_MAD{dose}"),
                    xid="time",
                    yid=key,
                    label=dose_label,
                    color=self.dose_colors[dose],
                )
                # data
                plot.add_data(
                    dataset=self.to_sid(f"{name}_MAD{dose}"),
                    xid="time",
                    yid="mean",
                    yid_sd=None if name == "uge" else "mean_sd",
                    label=f"{dose_label} data",
                    color=self.dose_colors[dose],
                    linestyle="",
                )
        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(FDAMB102002, output_dir=FDAMB102002.__name__)
