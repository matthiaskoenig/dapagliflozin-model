"""Simulation experiment of the study FDAMB102007."""

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


class FDAMB102007(DapagliflozinSimulationExperiment):
    """Simulation experiment of FDAMB102007.

    Study MB102007 of the clinical pharmacology review of NDA 202293, the renal
    impairment study of Kasichayanula2013 (same trial): healthy subjects and
    subjects with type 2 diabetes and normal renal function or mild, moderate
    or severe renal impairment took a single dose of 50 mg dapagliflozin on day
    1, and 20 mg once daily on days 4-10, fasted for at least 10 hr (food
    6 hr after the dose). The review groups the subjects by the MDRD eGFR
    (healthy n=5, normal 9, mild 9, moderate 8, severe 6), the publication by
    the Cockcroft-Gault creatinine clearance (8/12/8/8/4); the renal function
    is the class of `renal_map`, the body weight and the fasting plasma glucose
    are the means of the groups of the publication (Kasichayanula2013 Table 1,
    an approximation for the regrouped subjects). Urinary dapagliflozin (Table
    13: CLR x AUCinf, i.e. the amount to infinity, not reported) and D3G
    (Table 14: %UR) of the collection over 0-72 hr after the 50 mg dose, the
    urinary glucose excretion of 0-24 hr on day 10 (Fig3, medians of the T2DM
    groups, box plots, n=3/3/6/5).
    """

    groups: ClassVar[list[str]] = [
        "healthy",
        "normal",
        "mild",
        "moderate",
        "severe",
    ]
    colors: ClassVar[dict[str, str]] = {
        # the healthy subjects, the normal renal function of the T2DM groups is black
        "healthy": "tab:gray",
        "normal": DapagliflozinSimulationExperiment.renal_colors[
            "Normal renal function"
        ],
        "mild": DapagliflozinSimulationExperiment.renal_colors["Mild renal impairment"],
        "moderate": DapagliflozinSimulationExperiment.renal_colors[
            "Moderate renal impairment"
        ],
        "severe": DapagliflozinSimulationExperiment.renal_colors[
            "Severe renal impairment"
        ],
    }
    renal_functions: ClassVar[dict[str, float]] = {
        "healthy": DapagliflozinSimulationExperiment.renal_map["Normal renal function"],
        "normal": DapagliflozinSimulationExperiment.renal_map["Normal renal function"],
        "mild": DapagliflozinSimulationExperiment.renal_map["Mild renal impairment"],
        "moderate": DapagliflozinSimulationExperiment.renal_map[
            "Moderate renal impairment"
        ],
        "severe": DapagliflozinSimulationExperiment.renal_map[
            "Severe renal impairment"
        ],
    }
    #: body weight [kg] and fasting plasma glucose [mM] of the groups
    #: [Kasichayanula2013, Table 1]
    bodyweights: ClassVar[dict[str, float]] = {
        "healthy": 77.9,
        "normal": 97.2,
        "mild": 73.9,
        "moderate": 77.4,
        "severe": 78.4,
    }
    fpgs: ClassVar[dict[str, float]] = {
        "healthy": 4.94,
        "normal": 9.32,
        "mild": 7.75,
        "moderate": 6.79,
        "severe": 7.23,
    }
    info: ClassVar[list[tuple[str, str]]] = [
        ("Aurine_dap", "dapagliflozin_cumulative amount"),
        ("Aurine_d3g", "dapagliflozin 3-o-glucuronide_cumulative amount"),
        ("KI__UGE", "uge"),
    ]

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig3", "Tab13A", "Tab14A"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                if "_recovery_" in label:
                    # the recovery in % of the dose repeats the amounts
                    continue
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
        for group in self.groups:
            changes = {
                # physiological changes
                "BW": Q(self.bodyweights[group], "kg"),
                "[KI__fpg]": Q(self.fpgs[group], "mM"),
                "GU__fed": Q(self.fasting_map["fasted"], "dimensionless"),
                "f_cirrhosis": Q(self.cirrhosis_map["Control"], "dimensionless"),
                "KI__f_renal_function": Q(self.renal_functions[group], "dimensionless"),
            }
            # single dose simulation DAP50, urine collected over 0-72 hr
            simulations[self.to_sid(f"po_dap50_{group}")] = Simulation(
                time_unit="hr",
                end=73,
                preinit_changes={
                    **self.default_changes(),
                    **changes,
                    # dose (IVDOSE, PODOSE)
                    "PODOSE_dap": Q(50, "mg"),
                },
            )
            # 50 mg on day 1, 20 mg on days 4-10 (72-216 hr), the urinary
            # glucose of 0-24 hr on day 10 (time 0 is the last dose)
            simulations[self.to_sid(f"po_dap20_{group}")] = Simulation(
                time_unit="hr",
                end=9 * 24 + 48,
                preinit_changes={
                    **self.default_changes(),
                    **changes,
                    # dose
                    "PODOSE_dap": Q(50, "mg"),
                },
                changes=[
                    # the urine of every day is collected on its own
                    Change(
                        [24, 48],
                        {"KI__glc_urine": Q(0, "mmole")},
                    ),
                    Change(
                        [24 * k for k in range(3, 10)],
                        {"KI__glc_urine": Q(0, "mmole")},  # reset urinary data
                    ),
                    self.add_doses(
                        [24 * k for k in range(3, 10)], {"PODOSE_dap": Q(20, "mg")}
                    ),
                ],
                time_shift=-9 * 24,
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        # Fig3
        for group in self.groups:
            if group == "healthy":
                continue
            health = Health.T2DM if group == "normal" else Health.T2DM_RENAL_IMPAIRMENT
            mappings[self.to_sid(f"fm_dap20_{group}")] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=self.to_sid(f"uge_DAP20_{group}"),
                    xid="time",
                    yid="mean",
                    yid_sd=None,
                    count="count",
                ),
                observable=FitData(
                    self,
                    task=self.to_sid(f"task_po_dap20_{group}"),
                    xid="time",
                    yid="KI__UGE",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.URINE,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.MULTIPLE,
                    health=health,
                    fasting=Fasting.FASTED,
                ),
            )
        # Tab13A, Tab14A
        for group in self.groups:
            health = (
                Health.HEALTHY
                if group == "healthy"
                else Health.T2DM
                if group == "normal"
                else Health.T2DM_RENAL_IMPAIRMENT
            )
            # simulation
            for sid, substance in [
                ("Aurine_dap", "dapagliflozin_cumulative amount"),
                ("Aurine_d3g", "dapagliflozin 3-o-glucuronide_cumulative amount"),
            ]:
                mappings[self.to_sid(f"fm_dap50_{group}_{substance}")] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=self.to_sid(f"{substance}_DAP50_{group}"),
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=self.to_sid(f"task_po_dap50_{group}"),
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.URINE,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=health,
                        fasting=Fasting.FASTED,
                    ),
                )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self.figure_fig3(),
            **self.figure_tab13_14(),
        }

    def figure_fig3(self) -> dict[str, Figure]:
        """Urinary glucose excretion after multiple doses of 20 mg (Fig3)."""
        fig = Figure(
            experiment=self,
            sid="Fig3",
            name=f"{self.__class__.__name__} (T2DM)",
            # the legend inside the panel covers the simulations
            settings={"legend_position": "outside"},
        )
        # the last two doses, time 0 is the dose of day 10
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, min=-26, max=26),
            legend=True,
        )
        plots[0].set_yaxis(self.label_uge, unit=self.unit_uge)
        for group in self.groups:
            if group == "healthy":
                continue
            dose_label = f"{group} - 20 mg PO"
            # simulation
            plots[0].add_data(
                task=self.to_sid(f"task_po_dap20_{group}"),
                xid="time",
                yid="KI__UGE",
                label=dose_label,
                color=self.colors[group],
            )
            # data
            plots[0].add_data(
                dataset=self.to_sid(f"uge_DAP20_{group}"),
                xid="time",
                yid="mean",
                yid_sd=None,
                count="count",
                label=f"{dose_label} median",
                color=self.colors[group],
                linestyle="",
            )

        return self.figures_by_sid(fig)

    def figure_tab13_14(self) -> dict[str, Figure]:
        """Urinary dapagliflozin and D3G after a single dose of 50 mg (Tab13, Tab14)."""
        fig = Figure(
            experiment=self,
            sid="Tab13_14",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__} (T2DM)",
            # one legend of the row, right of the grid
            settings={"legend_scope": "row"},
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_dap_urine, unit=self.unit_dap_urine)
        plots[1].set_yaxis(self.label_d3g_urine, unit=self.unit_d3g_urine)

        for group in self.groups:
            dose_label_dap = f"{group} - 50 mg PO"
            dose_label_d3g = f"{group} - 50 mg PO"
            # simulation
            plots[0].add_data(
                task=self.to_sid(f"task_po_dap50_{group}"),
                xid="time",
                yid="Aurine_dap",
                label=dose_label_dap,
                color=self.colors[group],
            )
            plots[1].add_data(
                task=self.to_sid(f"task_po_dap50_{group}"),
                xid="time",
                yid="Aurine_d3g",
                label=dose_label_d3g,
                color=self.colors[group],
            )
            # data
            plots[0].add_data(
                dataset=self.to_sid(f"dapagliflozin_cumulative amount_DAP50_{group}"),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=dose_label_dap,
                linestyle="",
                color=self.colors[group],
            )
            plots[1].add_data(
                dataset=self.to_sid(
                    f"dapagliflozin 3-o-glucuronide_cumulative amount_DAP50_{group}"
                ),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=dose_label_d3g,
                linestyle="",
                color=self.colors[group],
            )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(FDAMB102007, output_dir=FDAMB102007.__name__)
