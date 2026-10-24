"""Simulation experiment of the study Gould2013."""

from typing import ClassVar, override

import matplotlib.pyplot as plt
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
from pkpdsim.scan import point_band, point_timecourse
from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.fit import FitData, FitMapping
from sbmlsim.plot import Axis, Figure
from sbmlsim.plot.serialization_matplotlib import FigureMPL
from sbmlsim.simulation import Simulation

from pkdb_models.models.dapagliflozin import MODEL
from pkdb_models.models.dapagliflozin.experiments.base_experiment import (
    DapagliflozinSimulationExperiment,
)


class Gould2013(DapagliflozinSimulationExperiment):
    """Simulation experiment of Gould2013.

    Dose response of the urinary glucose excretion (UGE) over 0-24 hr after a
    single dose in healthy subjects (Fig1). The doses of 0.001-2.5 mg are the
    low dose study of the publication: 36 healthy men and women (body weight
    not reported), an oral solution of 0.001, 0.01, 0.1 or 0.3 mg (n=5) or a
    tablet of 1 or 2.5 mg (n=6), the fasting state not reported; the values
    are the UGE of day 1 of Table 3 (not corrected for the baseline of 51-119
    mg/day). Placebo and 5-500 mg are the pooled standard clinical trials (153
    subjects on dapagliflozin, 5-40 per dose, 22 on placebo, fasted, several
    formulations, partly the subjects of Komoroski2009): the baseline corrected
    means of Fig1 plus the placebo 0.085 g, the SEM converted to an SD with n=6
    (n unknown). The model has no glucose excretion at the baseline glucose:
    placebo and 0.001-0.1 mg (0.04-0.09 g, the baseline glucosuria) are below
    its resolution.
    """

    #: the low dose study of the publication, the other doses are pooled trials
    doses_low: ClassVar[list[float]] = [0.001, 0.01, 0.1, 0.3, 1, 2.5]

    doses: ClassVar[list[float]] = [
        0,
        0.001,
        0.01,
        0.1,
        0.3,
        1,
        2.5,
        5,
        10,
        20,
        50,
        100,
        250,
        500,
    ]

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                # unit conversion
                if label.startswith("dapagliflozin_"):
                    dset.unit_conversion("mean", 1 / self.Mr.dap)
                dsets[self.to_sid(label)] = dset
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for dose in self.doses:
            simulations[self.to_sid(f"po_dap{dose}")] = Simulation(
                time_unit="hr",
                end=25,
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
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for dose in self.doses:
            # solution 0.001-0.3 mg, tablet 1 and 2.5 mg, the pooled trials
            # used capsules, solutions and tablets
            if dose in [0.001, 0.01, 0.1, 0.3]:
                applicationform = ApplicationForm.SOLUTION
            elif dose in self.doses_low:
                applicationform = ApplicationForm.TABLET
            else:
                applicationform = ApplicationForm.MIXED
            mappings[self.to_sid(f"fm_dap{dose}")] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=self.to_sid(f"glucose_cumulative amount_DAP{dose}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task=self.to_sid(f"task_po_dap{dose}"),
                    xid="time",
                    yid="KI__UGE",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.URINE,
                    route=Route.PO,
                    application_form=applicationform,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.NR if dose in self.doses_low else Fasting.FASTED,
                ),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1_timecourses",
            num_rows=1,
            num_cols=1,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        # the legend is outside of the panels, inside it covers the data
        fig.settings = {**fig.settings, "legend_position": "outside"}
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_uge, unit=self.unit_uge)
        # simulation
        for dose in self.doses:
            plots[0].add_data(
                task=self.to_sid(f"task_po_dap{dose}"),
                xid="time",
                yid="KI__UGE",
                label=f"Sim {dose} mg",
                color=self.dose_colors[dose],
            )
            # data
            plots[0].add_data(
                dataset=self.to_sid(f"glucose_cumulative amount_DAP{dose}"),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"UGE {dose} mg",
                color=self.dose_colors[dose],
                linestyle="",
            )
        return self.figures_by_sid(fig)

    @override
    def figures_mpl(self) -> dict[str, FigureMPL]:
        return {"Fig1": self._figure_dose_response()}

    def _uge_24hr(self, dose: float) -> tuple[float, float, float]:
        """Get the simulated UGE at 24 hr of a dose [g].

        Returns:
            The median and the band of the draws of the fitted parameters in a
            run over them (`uncertainty`), the value three times without draws.
        """
        res = self.results[self.to_sid(f"task_po_dap{dose}")]
        factor = Q(1, res.units["KI__UGE"]).to(self.unit_uge).magnitude
        # the time of the result is in the time unit of the model
        t24 = Q(24.0, "hr").to(res.units["time"]).magnitude
        run = self.uncertainty
        if run is None:
            time, values = point_timecourse(res, "KI__UGE")
            value = float(np.interp(t24, time, values)) * factor
            return value, value, value
        time, lower, median, upper = point_band(
            res, "KI__UGE", across=run.id, quantiles=run.quantiles
        )
        low, mid, high = (
            float(np.interp(t24, time, values)) * factor
            for values in (lower, median, upper)
        )
        return mid, low, high

    def _figure_dose_response(self) -> FigureMPL:
        """UGE over 0-24 hr against the dose, log scale of the dose as the source."""
        df = self.load_dataframe("Fig1")
        df = df[df.time == 24.0]
        doses = [dose for dose in self.doses if dose > 0]
        simulated = np.array([self._uge_24hr(dose) for dose in doses])

        f, ax = plt.subplots(figsize=(8, 6), dpi=150, layout="constrained")
        ax.fill_between(
            doses, simulated[:, 1], simulated[:, 2], color="black", alpha=0.2
        )
        ax.plot(doses, simulated[:, 0], color="black", label="simulation")
        placebo_sim = self._uge_24hr(0)[0]
        ax.axhline(
            placebo_sim, color="black", linestyle=":", label="simulation placebo"
        )
        for dose in self.doses:
            row = df[df.intervention == f"DAP{dose:g}"]
            mean = float(row["mean"].iloc[0])
            sd = float(row["sd"].iloc[0]) if not row["sd"].isna().iloc[0] else 0.0
            low_dose = dose in self.doses_low
            if dose == 0:
                ax.axhline(mean, color="tab:gray", linestyle="--", label="data placebo")
                continue
            ax.errorbar(
                x=dose,
                y=mean,
                yerr=sd,
                color=self.dose_colors[dose],
                marker="o" if low_dose else "s",
                markersize=9,
                markeredgecolor="black",
                capsize=4,
                linestyle="none",
                label=(
                    "data low dose study (Tab3)"
                    if dose == self.doses_low[0]
                    else "data pooled trials (Fig1)"
                    if dose == 5
                    else None
                ),
            )
        ax.set_xscale("log")
        ax.set_xlabel("dose [mg]", fontdict=self.font)
        ax.set_ylabel(f"{self.label_uge} 0-24 hr [{self.unit_uge}]", fontdict=self.font)
        ax.tick_params(labelsize=self.tick_font_size)
        ax.set_ylim(bottom=0)
        ax.legend(fontsize=12)
        f.suptitle(f"{self.__class__.__name__} (Healthy)", fontsize=20)
        return f


if __name__ == "__main__":
    MODEL.run_experiments(Gould2013, output_dir=Gould2013.__name__)
