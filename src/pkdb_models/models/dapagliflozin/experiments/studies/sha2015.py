"""Simulation experiment of the study Sha2015."""

from typing import ClassVar, override

import numpy as np
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


class Sha2015(DapagliflozinSimulationExperiment):
    """Simulation experiment of Sha2015.

    10 mg dapagliflozin once daily for 4 days in healthy subjects, the
    canagliflozin arm of the crossover is not simulated. The dose of day 4 is
    given 20 min before a mixed meal tolerance test (MMTT, 600 kcal, 75 g
    glucose) after an overnight fast. Time 0 is the dose of day 4 (Fig 2, plasma
    at the nominal times); the urine (UGE, RTG of Fig 3A/B) and the plasma
    glucose (Fig 3C) are timed from the start of the MMTT and shifted by the
    20 min.
    """

    colors: ClassVar[dict[float, str]] = {10: "black"}
    bodyweight = 78.9  # [kg]
    #: start of the MMTT after the dose of day 4 [hr]
    t_mmtt = 1 / 3
    #: fasting plasma glucose, the mean of the samples before the MMTT [mM]
    fpg = 4.86
    #: the mean plasma glucose of the day above the fasting glucose [mM]
    pg_delta_default = 0.85

    #: the intervals of the urine of the RTG, by the end (curation) [hr]
    rtg_midpoints: ClassVar[dict[float, float]] = {
        4.0: 2.0,
        10.0: 7.0,
        14.0: 12.0,
        24.0: 19.0,
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2", "Fig3"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                if "DAP10" not in label:
                    continue
                df_label = df_label.copy()
                if fig_id == "Fig3":
                    if label.startswith("renal_threshold"):
                        # the RTG is the mean of a urine interval: at its
                        # midpoint, not at its end
                        df_label["time"] = df_label["time"].map(self.rtg_midpoints)
                    elif label.startswith("glucose_cumulative_amount"):
                        # no information: the collection starts at the MMTT
                        # with the reset of the urine
                        df_label = df_label[df_label["time"] > 0]
                    # times from the start of the MMTT to times from the dose
                    df_label["time"] = df_label["time"] + self.t_mmtt
                dsets[self.to_sid(label)] = DataSet.from_df(df_label, self.ureg)
        return dsets

    def pg_changes(self) -> list[Change]:
        """The plasma glucose of the MMTT as the mean glucose of the model.

        The fasting glucose in the overnight fast (8 hr) before the MMTT, the
        measured glucose (Fig 3C, 0 to 2 hr of the MMTT) linearly interpolated
        in steps of 3 min, as the difference to the fasting glucose
        (`KI__PG_delta`), and from 2 hr after the start of the MMTT the mean
        glucose of the day of the model (`PG_delta` 0.85 mM, the other meals
        were not measured).

        Returns:
            The changes in the time of the simulation (dose of day 4 at 72 hr).
        """
        df = self.load_dataframe("Fig3")
        df = df[(df["label"] == "glucose_DAP10") & (df["time"] >= 0)]
        times = df["time"].to_numpy(dtype=float)
        values = df["mean"].to_numpy(dtype=float)
        steps = np.arange(0.0, times[-1], 0.05)
        changes = [
            Change([72 + self.t_mmtt - 8], {"KI__PG_delta": Q(0, "mM")}),
        ]
        changes += [
            Change(
                [72 + self.t_mmtt + t],
                {"KI__PG_delta": Q(float(pg) - self.fpg, "mM")},
            )
            for t, pg in zip(steps, np.interp(steps, times, values), strict=True)
        ]
        changes.append(
            Change(
                [72 + self.t_mmtt + times[-1]],
                {"KI__PG_delta": Q(self.pg_delta_default, "mM")},
            )
        )
        return changes

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        # 4 daily doses, the dose of day 4 at time 0; the urine is collected from
        # the start of the MMTT
        simulations["po_dap10"] = Simulation(
            time_unit="hr",
            end=3 * 24 + 25,
            preinit_changes={
                **self.default_changes(),
                # physiological changes
                "BW": Q(self.bodyweight, "kg"),
                "[KI__fpg]": Q(self.fpg, "mM"),
                "KI__f_renal_function": Q(
                    self.renal_map["Normal renal function"], "dimensionless"
                ),
                "GU__fed": Q(self.fasting_map["fasted"], "dimensionless"),
                "f_cirrhosis": Q(self.cirrhosis_map["Control"], "dimensionless"),
                # dose (IVDOSE, PODOSE)
                "PODOSE_dap": Q(10, "mg"),
            },
            changes=[
                # the urine per day (not measured before day 4)
                Change([24, 48, 72], {"KI__glc_urine": Q(0, "mmole")}),
                self.add_doses([24, 48, 72], {"PODOSE_dap": Q(10, "mg")}),
                Change([72 + self.t_mmtt], {"KI__glc_urine": Q(0, "mmole")}),
                *self.pg_changes(),
            ],
            time_shift=-3 * 24,
        )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for name, sid, tissue in [
            ("dapagliflozin", "[Cve_dap]", Tissue.PLASMA),
            ("glucose_cumulative_amount", "KI__UGE", Tissue.URINE),
            ("renal_threshold", "KI__RTG", Tissue.NONE),
        ]:
            mappings[self.to_sid(f"fm_po_dap10_{name}")] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=self.to_sid(f"{name}_DAP10"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task="task_po_dap10",
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
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        figures = []
        infos = [
            (
                "Fig2",
                "dapagliflozin",
                "[Cve_dap]",
                self.label_dap_plasma,
                self.unit_dap,
            ),
            (
                "Fig3A",
                "glucose_cumulative_amount",
                "KI__UGE",
                self.label_uge,
                self.unit_uge,
            ),
            ("Fig3B", "renal_threshold", "KI__RTG", self.label_rtg, self.unit_rtg),
            ("Fig3C", "glucose", "KI__pg", self.label_pg, self.unit_glc),
        ]
        for fig_id, name, sid, label, unit in infos:
            fig = Figure(
                experiment=self,
                sid=fig_id,
                num_rows=1,
                num_cols=1,
                name=f"{self.__class__.__name__} (Healthy)",
            )
            plot = fig.create_plots(
                xaxis=Axis(self.label_time, unit=self.unit_time, min=-1, max=25),
                legend=True,
            )[0]
            plot.set_yaxis(label, unit=unit)
            plot.add_data(
                task="task_po_dap10",
                xid="time",
                yid=sid,
                label="10 mg PO, day 4",
                color=self.colors[10],
            )
            plot.add_data(
                dataset=self.to_sid(f"{name}_DAP10"),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label="10 mg PO, day 4",
                color=self.colors[10],
            )
            figures.append(fig)
        return self.figures_by_sid(*figures)


if __name__ == "__main__":
    MODEL.run_experiments(Sha2015, output_dir=Sha2015.__name__)
