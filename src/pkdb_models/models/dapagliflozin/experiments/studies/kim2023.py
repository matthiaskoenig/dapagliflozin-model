"""Simulation experiment of the study Kim2023."""

from typing import ClassVar, override

from pkpdsim.metadata import (
    ApplicationForm,
    Coadministration,
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


class Kim2023(DapagliflozinSimulationExperiment):
    """Simulation experiment of Kim2023.

    Multiple doses of 10 mg dapagliflozin once daily for 5 days alone and with
    5 mg evogliptin (EV5) in healthy Korean men (arm 2, 69.1 kg), fasted;
    plasma 0-24 hr and urinary glucose 0-24 hr after the 5th dose (TabS2: LS
    means, the SD from the 90 % CI). Evogliptin does not change the exposure
    of dapagliflozin (GMR AUC 1.02, Cmax 1.09) nor the urinary glucose: no
    interaction factor, the evogliptin arm is validation data.
    """

    interventions: ClassVar[list[str]] = ["DAP10", "EV5, DAP10"]
    colors: ClassVar[dict[str, str]] = {
        "DAP10": "black",
        "EV5, DAP10": "tab:blue",
    }
    info: ClassVar[dict[str, str]] = {
        "[Cve_dap]": "dapagliflozin",
        "KI__UGE": "uge",
    }
    bodyweight = 69.1  # [kg] Group2, Arm2

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2", "TabS2"]:
            df = self.load_dataframe(fig_id)
            if fig_id == "TabS2":
                # the urinary glucose of day 5 starts with the reset at the 5th
                # dose: a data point at the time of a change is compared with
                # the state before it, i.e., the 0 with the glucose of day 4
                df = df[df["time"] > 0]
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
        simulations["po_dap10"] = Simulation(
            time_unit="hr",
            end=4 * 24 + 48,
            preinit_changes={
                **self.default_changes(),
                # physiological changes
                "BW": Q(self.bodyweight, "kg"),
                "[KI__fpg]": Q(self.fpg_healthy, "mM"),
                "GU__fed": Q(self.fasting_map["fasted"], "dimensionless"),
                "f_cirrhosis": Q(self.cirrhosis_map["Control"], "dimensionless"),
                "KI__f_renal_function": Q(
                    self.renal_map["Normal renal function"], "dimensionless"
                ),
                # dose (IVDOSE, PODOSE)
                "PODOSE_dap": Q(10, "mg"),
            },
            changes=[
                Change(
                    [24 * k for k in range(1, 5)],
                    {"KI__glc_urine": Q(0, "mmole")},  # reset UGE
                ),
                self.add_doses(
                    [24 * k for k in range(1, 5)], {"PODOSE_dap": Q(10, "mg")}
                ),
            ],
            time_shift=-4 * 24,
        )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for intervention in self.interventions:
            for sid, name in self.info.items():
                mappings[self.to_sid(f"fm_dap10_{intervention}_{name}")] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=self.to_sid(f"{name}_{intervention}"),
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
                        tissue=Tissue.URINE if name == "uge" else Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.MULTIPLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.EVOGLIPTIN
                        if "EV5" in intervention
                        else Coadministration.NONE,
                    ),
                )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2_S2",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, min=-5, max=25),
            legend=True,
        )
        plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        plots[1].set_yaxis(self.label_uge, unit=self.unit_uge)
        # simulation
        for k, (sid, name) in enumerate(self.info.items()):
            linestyle = "--" if sid == "[Cve_dap]" else ""
            plots[k].add_data(
                task="task_po_dap10",
                xid="time",
                yid=sid,
                label="10 mg PO",
                color="black",
            )
            # data
            for intervention in self.interventions:
                label = "10 mg PO" if intervention == "DAP10" else "EV5 + 10 mg PO"
                plots[k].add_data(
                    dataset=self.to_sid(f"{name}_{intervention}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=label,
                    color=self.colors[intervention],
                    linestyle=linestyle,
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Kim2023, output_dir=Kim2023.__name__)
