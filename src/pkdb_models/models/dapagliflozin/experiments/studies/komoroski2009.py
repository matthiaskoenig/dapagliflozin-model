"""Simulation experiment of the study Komoroski2009."""

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


class Komoroski2009(DapagliflozinSimulationExperiment):
    """Simulation experiment of Komoroski2009.

    Single (0 - 500 mg, 250 mg fasted and fed) and multiple (2.5 - 100 mg once
    daily for 14 days) doses of dapagliflozin in healthy men (group means
    74-84 kg). SAD (MB102001): 2.5-5 mg oral solution, 20-500 mg capsule (10 mg
    not stated), after a 10 hr fast, urine for at least 120 hr (48 hr for 5
    and 20 mg); the 250 mg panel took a second dose after a high-fat
    breakfast (period 2). MAD (MB102002): capsules once daily for 14 days,
    plasma 0-24 hr on day 1 and 14, urinary glucose 0-24 hr on day 1 and 14
    (the amounts of the text of the results, Fig6 are bars).

    Curation of the copies (digitization): the times of the day 14 profiles
    are recalibrated (offset of the x axis, t = 1.019 t_dig - 0.46 hr) and set
    to the sampling times, their points at the time axis before 0.5 hr
    (1.1-3.1 ng/ml, far below the trough of 1.3-63 ng/ml at 24 hr) are
    removed; the urine of 250 mg (Table 1) is of 0-120 hr; the 5 mg UGE after
    48 hr is removed (urine collected 0-48 hr).
    """

    conditions: ClassVar[list[str]] = ["fasted", "fed"]
    doses_multi: ClassVar[list[float]] = [2.5, 10, 20, 50, 100]
    doses_uge_single: ClassVar[list[float]] = [0, 2.5, 5, 10, 20, 50, 100, 250, 500]
    doses_uge_multi: ClassVar[list[float]] = [2.5, 10, 20, 50, 100]
    colors: ClassVar[dict[float, str]] = {
        0: "black",
        2.5: "tab:blue",
        5: "tab:orange",
        10: "tab:green",
        20: "tab:red",
        50: "tab:purple",
        100: "tab:brown",
        250: "tab:pink",
        500: "tab:grey",
    }
    bodyweight = (74 + 84) / 2  # [kg] range [74, 84]

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig3", "Tab1A", "Fig4", "Fig5", "Fig6"]:
            df = self.load_dataframe(fig_id)
            # the recovery in % of the dose duplicates the cumulative amounts
            df = df[df["measurement_type"] != "recovery"]
            if fig_id == "Fig6":
                # the urinary glucose of day 14 starts with the reset at the 14th
                # dose: a data point at the time of a change is compared with
                # the state before it, i.e., the 0 with the glucose of day 13
                df = df[df["time"] != 13 * 24]
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
        # fasted and fed, urine 0-120 hr
        for condition in self.conditions:
            simulations[self.to_sid(f"po_dap_250_{condition}")] = Simulation(
                time_unit="hr",
                end=121,
                preinit_changes={
                    **self.default_changes(),
                    # physiological changes
                    "BW": Q(self.bodyweight, "kg"),
                    "GU__fed": Q(self.fasting_map[condition], "dimensionless"),
                    "[KI__fpg]": Q(self.fpg_healthy, "mM"),
                    "f_cirrhosis": Q(self.cirrhosis_map["Control"], "dimensionless"),
                    "KI__f_renal_function": Q(
                        self.renal_map["Normal renal function"], "dimensionless"
                    ),
                    # dose (IVDOSE, PODOSE)
                    "PODOSE_dap": Q(250, "mg"),
                },
            )
        # single dose simulation
        for dose in self.doses_uge_single:
            simulations[self.to_sid(f"po_dap_{dose}_single")] = Simulation(
                time_unit="hr",
                end=125,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweight, "kg"),
                    "GU__fed": Q(self.fasting_map["fasted"], "dimensionless"),
                    "[KI__fpg]": Q(self.fpg_healthy, "mM"),
                    "f_cirrhosis": Q(self.cirrhosis_map["Control"], "dimensionless"),
                    "KI__f_renal_function": Q(
                        self.renal_map["Normal renal function"], "dimensionless"
                    ),
                    "PODOSE_dap": Q(dose, "mg"),
                },
            )
        # multi dose simulation, the profile after the 14th dose
        for dose in self.doses_multi:
            simulations[self.to_sid(f"po_dap_{dose}_multi")] = Simulation(
                time_unit="hr",
                end=13 * 24 + 30,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweight, "kg"),
                    "GU__fed": Q(self.fasting_map["fasted"], "dimensionless"),
                    "[KI__fpg]": Q(self.fpg_healthy, "mM"),
                    "f_cirrhosis": Q(self.cirrhosis_map["Control"], "dimensionless"),
                    "KI__f_renal_function": Q(
                        self.renal_map["Normal renal function"], "dimensionless"
                    ),
                    "PODOSE_dap": Q(dose, "mg"),
                },
                changes=[
                    self.add_doses(
                        [24 * k for k in range(1, 14)], {"PODOSE_dap": Q(dose, "mg")}
                    ),
                ],
                time_shift=-13 * 24,
            )
        # multi dose simulation, the urinary glucose of every day
        for dose in self.doses_uge_multi:
            simulations[self.to_sid(f"po_dap_{dose}_uge_multi")] = Simulation(
                time_unit="hr",
                end=13 * 24 + 30,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweight, "kg"),
                    "GU__fed": Q(self.fasting_map["fasted"], "dimensionless"),
                    "[KI__fpg]": Q(self.fpg_healthy, "mM"),
                    "f_cirrhosis": Q(self.cirrhosis_map["Control"], "dimensionless"),
                    "KI__f_renal_function": Q(
                        self.renal_map["Normal renal function"], "dimensionless"
                    ),
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
        # fasting
        for condition in self.conditions:
            mappings[self.to_sid(f"fm_po_dap250_{condition}")] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=self.to_sid(f"dapagliflozin_DAP250_{condition}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task=self.to_sid(f"task_po_dap_250_{condition}"),
                    xid="time",
                    yid="[Cve_dap]",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.CAPSULE,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED if condition == "fasted" else Fasting.FED,
                ),
            )
        for condition in self.conditions:
            mappings[self.to_sid(f"fm_po_dap250_{condition}_urine")] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=self.to_sid(f"dapagliflozin_cumulative amount_{condition}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task=self.to_sid(f"task_po_dap_250_{condition}"),
                    xid="time",
                    yid="Aurine_dap",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.URINE,
                    route=Route.PO,
                    application_form=ApplicationForm.CAPSULE,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED if condition == "fasted" else Fasting.FED,
                ),
            )
        # single doses, urinary glucose
        for dose in self.doses_uge_single:
            mappings[self.to_sid(f"fm_po_dap_{dose}_single")] = FitMapping(
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
                    task=self.to_sid(f"task_po_dap_{dose}_single"),
                    xid="time",
                    yid="KI__UGE",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.URINE,
                    route=Route.PO,
                    # the doses of 0 to 5 mg were a solution
                    application_form=(
                        ApplicationForm.SOLUTION
                        if dose in [0, 2.5, 5]
                        else ApplicationForm.CAPSULE
                    ),
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
                ),
            )
        # multi
        for dose in self.doses_multi:
            mappings[self.to_sid(f"fm_po_dap{dose}_MD_single")] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=self.to_sid(f"dapagliflozin_DAP{dose}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task=self.to_sid(f"task_po_dap_{dose}_single"),
                    xid="time",
                    yid="[Cve_dap]",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.CAPSULE,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
                ),
            )
        for dose in self.doses_multi:
            mappings[self.to_sid(f"fm_po_dap{dose}_MD_multi")] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=self.to_sid(f"dapagliflozin_DAPM{dose}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task=self.to_sid(f"task_po_dap_{dose}_multi"),
                    xid="time",
                    yid="[Cve_dap]",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.CAPSULE,
                    dosing=Dosing.MULTIPLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
                ),
            )
        for dose in self.doses_uge_multi:
            mappings[self.to_sid(f"fm_po_dap_{dose}_uge_multi")] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=self.to_sid(f"glucose_cumulative amount_DAPM{dose}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task=self.to_sid(f"task_po_dap_{dose}_uge_multi"),
                    xid="time",
                    yid="KI__UGE",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.URINE,
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
            **self.figure_fig3(),
            **self.figure_fig5(),
            **self.figure_fig4(),
            **self.figure_fig6(),
        }

    def figure_fig3(self) -> dict[str, Figure]:
        """Single dose of 250 mg in the fasted and fed state (Fig3, Table 1)."""
        fig = Figure(
            experiment=self,
            sid="Fig3",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        plots[0].set_xaxis(self.label_time, unit=self.unit_time, min=-1, max=25)
        plots[1].set_yaxis(self.label_dap_urine, unit=self.unit_dap_urine)
        for condition in self.conditions:
            # simulation
            plots[0].add_data(
                task=self.to_sid(f"task_po_dap_250_{condition}"),
                xid="time",
                yid="[Cve_dap]",
                label=f"{condition} - 250 mg PO",
                color=self.fasting_colors[condition],
            )
            plots[1].add_data(
                task=self.to_sid(f"task_po_dap_250_{condition}"),
                xid="time",
                yid="Aurine_dap",
                label=f"{condition} - 250 mg PO",
                color=self.fasting_colors[condition],
            )
            # data
            plots[0].add_data(
                dataset=self.to_sid(f"dapagliflozin_DAP250_{condition}"),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{condition} - 250 mg PO",
                color=self.fasting_colors[condition],
            )
            plots[1].add_data(
                dataset=self.to_sid(f"dapagliflozin_cumulative amount_{condition}"),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{condition} - 250 mg PO",
                color=self.fasting_colors[condition],
                linestyle="",
            )
        return self.figures_by_sid(fig)

    def figure_fig5(self) -> dict[str, Figure]:
        """Urinary glucose excretion after the single doses (Fig5)."""
        fig = Figure(
            experiment=self,
            sid="Fig5",
            name=f"{self.__class__.__name__} (Healthy)",
        )
        # the legend is outside of the panels, inside it covers the data
        fig.settings = {**fig.settings, "legend_position": "outside"}
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_uge, unit=self.unit_uge)
        for dose in self.doses_uge_single:
            plots[0].add_data(
                task=self.to_sid(f"task_po_dap_{dose}_single"),
                xid="time",
                yid="KI__UGE",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
            )
        for dose in self.doses_uge_single:
            plots[0].add_data(
                dataset=self.to_sid(f"glucose_cumulative amount_DAP{dose}"),
                xid="time",
                yid="mean",
                yid_sd=None,
                count="count",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
            )
        return self.figures_by_sid(fig)

    def figure_fig4(self) -> dict[str, Figure]:
        """Plasma profiles after the first and the 14th dose (Fig4c)."""
        fig = Figure(
            experiment=self,
            sid="Fig4",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        plots[1].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        plots[0].set_xaxis(self.label_time, unit=self.unit_time, min=-3, max=30)
        plots[1].set_xaxis(self.label_time, unit=self.unit_time, min=-24, max=30)

        for dose in self.doses_multi:
            # simulation
            plots[0].add_data(
                task=self.to_sid(f"task_po_dap_{dose}_single"),
                xid="time",
                yid="[Cve_dap]",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
            )
            plots[1].add_data(
                task=self.to_sid(f"task_po_dap_{dose}_multi"),
                xid="time",
                yid="[Cve_dap]",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
            )
            # data
            plots[0].add_data(
                dataset=self.to_sid(f"dapagliflozin_DAP{dose}"),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
            )
            plots[1].add_data(
                dataset=self.to_sid(f"dapagliflozin_DAPM{dose}"),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
            )
        return self.figures_by_sid(fig)

    def figure_fig6(self) -> dict[str, Figure]:
        """Urinary glucose of day 1 and day 14 of the multiple dosing (Fig6)."""
        fig = Figure(
            experiment=self,
            sid="Fig6",
            name=f"{self.__class__.__name__} (Healthy)",
        )
        # the legend is outside of the panels, inside it covers the data
        fig.settings = {**fig.settings, "legend_position": "outside"}
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_uge, unit=self.unit_uge)
        for dose in self.doses_uge_multi:
            plots[0].add_data(
                task=self.to_sid(f"task_po_dap_{dose}_uge_multi"),
                xid="time",
                yid="KI__UGE",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
            )
        for dose in self.doses_uge_multi:
            plots[0].add_data(
                dataset=self.to_sid(f"glucose_cumulative amount_DAPM{dose}"),
                xid="time",
                yid="mean",
                yid_sd=None,
                count="count",
                label=f"{dose} mg PO",
                color=self.dose_colors[dose],
                linestyle="",
            )
        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Komoroski2009, output_dir=Komoroski2009.__name__)
