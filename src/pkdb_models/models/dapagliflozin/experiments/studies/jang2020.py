"""Simulation experiment of the study Jang2020."""

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
from sbmlsim.simulation import Simulation

from pkdb_models.models.dapagliflozin import MODEL
from pkdb_models.models.dapagliflozin.experiments.base_experiment import (
    DapagliflozinSimulationExperiment,
)


class Jang2020(DapagliflozinSimulationExperiment):
    """Simulation experiment of Jang2020.

    Part I: 16 healthy Korean men (71.45 kg), 10 mg dapagliflozin tablets once
    daily on days 1-5, alone or with lobeglitazone 0.5 mg (crossover), after an
    overnight fast; plasma after the dose of day 5 (0-48 hr, mean and SD, time 0
    is the dose of day 5). Lobeglitazone does not change the exposure (AUCtau
    ratio 0.99, Cmax 0.92) and is not simulated, both arms are one simulation.
    """

    interventions: ClassVar[list[str]] = ["DAP10", "DAP10, LOB05"]
    bodyweight = 71.45  # [kg]

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2"]:
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
            mappings[self.to_sid(f"fm_dap10_{intervention}")] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=self.to_sid(f"dapagliflozin_{intervention}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
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
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
                    coadministration=Coadministration.LOBEGLITAZONE
                    if "LOB" in intervention
                    else Coadministration.NONE,
                ),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2",
            num_rows=1,
            num_cols=1,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        plots[0].set_xaxis(self.label_time, unit=self.unit_time, min=-24, max=50)

        label = "10 mg PO"
        # simulation
        plots[0].add_data(
            task="task_po_dap10",
            xid="time",
            yid="[Cve_dap]",
            label=label,
            color="black",
        )

        # data
        plots[0].add_data(
            dataset="dapagliflozin_DAP10",
            xid="time",
            yid="mean",
            yid_sd="mean_sd",
            count="count",
            label=label,
            color="black",
        )
        plots[0].add_data(
            dataset="dapagliflozin_DAP10__LOB05",
            xid="time",
            yid="mean",
            yid_sd="mean_sd",
            count="count",
            label="10 mg PO + LOB 0.5 mg",
            color="tab:blue",
        )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Jang2020, output_dir=Jang2020.__name__)
