"""Simulation experiment of the study FDAMB102003."""

import itertools
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


class FDAMB102003(DapagliflozinSimulationExperiment):
    """Simulation experiment of FDAMB102003.

    Study MB102003 of the clinical pharmacology review of NDA 202293 (T2DM
    study of Komoroski 2009, Clin Pharmacol Ther 85:513): subjects with type 2
    diabetes took 5 (n=11), 25 (n=12) or 100 mg (n=16) dapagliflozin once daily
    for 14 days (the "2.5, 10 mg" of the review text is a typo). The review
    reports neither the formulation (curated as tablet), the fasting state,
    the body weight nor the glucose: the T2DM defaults of the model
    (`fpg_t2dm`, `bodyweight_default`) are used until the publication is
    curated. Plasma on days 1 and 14 (Fig13, day 1 n=11/12/15, day 14
    n=10/12/14), the urinary dapagliflozin of 0-24 hr (Table 8, %UR as amount)
    and the cumulative urinary glucose excretion over 0-24 hr (Fig18, without
    SD) on days 1 and 14. The urine is collected per day and reset with every
    dose; a data point at the time of a reset is compared with the state
    before it. The placebo group is not simulated.
    """

    doses: ClassVar[list[float]] = [5, 25, 100]
    days: ClassVar[dict[int, float]] = {1: 0.0, 14: 312.0}

    info: ClassVar[dict[str, str]] = {
        "[Cve_dap]": "dapagliflozin",
        "Aurine_dap": "dapagliflozin_urine",
        "KI__UGE": "uge",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        # Fig13: dap
        # Tab8A: dap_urine
        # Fig18: UGE
        for fig_id in ["Fig13", "Fig18", "Tab8A"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                if "_recovery_" in label:
                    # the recovery in % of the dose repeats the amounts
                    continue
                # one dataset per day, the days differ in the number of
                # subjects
                for day, start in self.days.items():
                    df_day = df_label[
                        (df_label.time >= start) & (df_label.time <= start + 24)
                    ]
                    if df_day.empty:
                        continue
                    dset = DataSet.from_df(df_day, self.ureg)
                    if label.startswith("dapagliflozin_"):
                        dset.unit_conversion("mean", 1 / self.Mr.dap)
                    dsets[f"{label}_day{day}"] = dset
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for dose in self.doses:
            simulations[f"po_DAP{dose}"] = Simulation(
                time_unit="hr",
                end=13 * 24 + 25,
                preinit_changes={
                    **self.default_changes(),
                    # physiological changes
                    "BW": Q(self.bodyweight_default, "kg"),
                    "[KI__fpg]": Q(self.fpg_t2dm, "mM"),
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
                        {
                            "KI__glc_urine": Q(0, "mmole"),  # reset UGE
                            "Aurine_dap": Q(0, "mmole"),  # reset urinary dap
                        },
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
            for (sid, name), day in itertools.product(self.info.items(), self.days):
                mappings[f"fm_dap{dose}_{name}_day{day}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_MAD{dose}_day{day}",
                        xid="time",
                        yid="mean",
                        # the urinary glucose has no SD
                        yid_sd=None if name == "uge" else "mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=f"task_po_DAP{dose}",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.PLASMA if sid == "[Cve_dap]" else Tissue.URINE,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.MULTIPLE,
                        health=Health.T2DM,
                        fasting=Fasting.NR,
                    ),
                )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self._figure_days(sid="Fig13", key="[Cve_dap]"),
            **self._figure_days(sid="Tab8", key="Aurine_dap"),
            **self._figure_days(sid="Fig18", key="KI__UGE"),
        }

    def _figure_days(self, sid: str, key: str) -> dict[str, Figure]:
        """Plasma (Fig13), urinary dapagliflozin (Tab8) or glucose (Fig18) per day."""
        name = self.info[key]
        fig = Figure(
            experiment=self,
            sid=sid,
            num_rows=1,
            num_cols=len(self.days),
            name=f"{self.__class__.__name__} (T2DM)",
            # one legend of the row, right of the grid
            settings={"legend_scope": "row"},
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        for plot, (day, start) in zip(plots, self.days.items(), strict=True):
            plot.set_title(f"day {day}")
            plot.set_xaxis(
                self.label_time, unit=self.unit_time, min=start - 1, max=start + 25
            )
            if key == "[Cve_dap]":
                # log scale as the source
                plot.set_yaxis(self.label_dap_plasma, unit=self.unit_dap, scale="log")
            elif key == "Aurine_dap":
                plot.set_yaxis(self.label_dap_urine, unit=self.unit_dap_urine)
            else:
                plot.set_yaxis(self.label_uge, unit=self.unit_uge)
            for dose in self.doses:
                dose_label = f"{dose} mg PO"
                # simulation
                plot.add_data(
                    task=f"task_po_DAP{dose}",
                    xid="time",
                    yid=key,
                    label=dose_label,
                    color=self.dose_colors[dose],
                )
                # data
                plot.add_data(
                    dataset=f"{name}_MAD{dose}_day{day}",
                    xid="time",
                    yid="mean",
                    yid_sd=None if name == "uge" else "mean_sd",
                    # no count in the label: the days differ in n, the legend
                    # of the row shows a dose once
                    label=f"{dose_label} data",
                    color=self.dose_colors[dose],
                    linestyle="",
                )
        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(FDAMB102003, output_dir=FDAMB102003.__name__)
