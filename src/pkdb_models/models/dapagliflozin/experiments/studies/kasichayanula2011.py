"""Simulation experiment of the study Kasichayanula2011."""

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


class Kasichayanula2011(DapagliflozinSimulationExperiment):
    """Simulation experiment of Kasichayanula2011.

    Hepatic impairment: 6 healthy subjects and 6 each with mild, moderate and
    severe hepatic impairment (Child-Pugh A, B, C), a single 10 mg tablet after
    a fast of 6 hr; plasma dapagliflozin and D3G 0-96 hr (Fig1, Fig2 0-48 hr,
    mean with error bars of an unstated type, curated as SE).

    The estimated (Cockcroft-Gault) creatinine clearance of the groups is 138,
    128, 158 and 143 ml/min. Cockcroft-Gault overestimates the GFR in cirrhosis
    (body weight with ascites and edema, low creatinine production), and a renal
    function above 1 would raise the renal glucuronidation (60 % of the
    metabolism) and mask the hepatic effect: all groups have the normal renal
    function.
    """

    conditions: ClassVar[dict[str, str]] = {
        "healthy": "Control",
        "mild": "Mild cirrhosis",
        "moderate": "Moderate cirrhosis",
        "severe": "Severe cirrhosis",
    }
    bodyweights: ClassVar[dict[str, float]] = {
        "healthy": 81.8,
        "mild": 80.9,
        "moderate": 87.5,
        "severe": 89.2,
    }  # [kg]

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig2"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
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
        for condition, cirrhosis_class in self.conditions.items():
            simulations[self.to_sid(f"po_dap10_{condition}")] = Simulation(
                time_unit="hr",
                end=50,
                preinit_changes={
                    **self.default_changes(),
                    # physiological changes
                    "BW": Q(self.bodyweights[condition], "kg"),
                    "[KI__fpg]": Q(self.fpg_healthy, "mM"),
                    "GU__fed": Q(self.fasting_map["fasted"], "dimensionless"),
                    "f_cirrhosis": Q(
                        self.cirrhosis_map[cirrhosis_class], "dimensionless"
                    ),
                    "KI__f_renal_function": Q(
                        self.renal_map["Normal renal function"], "dimensionless"
                    ),
                    # dose (IVDOSE, PODOSE)
                    "PODOSE_dap": Q(10, "mg"),
                },
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for condition in self.conditions:
            health = Health.HEALTHY if condition == "healthy" else Health.CIRRHOSIS
            for sid, name in [
                ("dap", "dapagliflozin"),
                ("d3g", "dapagliflozin 3-o-glucuronide"),
            ]:
                mappings[self.to_sid(f"fm_dap10_{condition}_{sid}")] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=self.to_sid(f"{name}_DAP10_{condition}"),
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=self.to_sid(f"task_po_dap10_{condition}"),
                        xid="time",
                        yid=f"[Cve_{sid}]",
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.PLASMA,
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
        colors = {k: self.cirrhosis_colors[v] for k, v in self.conditions.items()}
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__} (Hepatic)",
        )
        # the legends inside the panels cover data
        fig.settings = {**fig.settings, "legend_position": "outside"}
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        plots[1].set_yaxis(self.label_d3g_plasma, unit=self.unit_d3g)

        for condition in self.conditions:
            dose_label = f"{condition} - 10 mg PO"
            # simulation: dapagliflozin plasma
            plots[0].add_data(
                task=self.to_sid(f"task_po_dap10_{condition}"),
                xid="time",
                yid="[Cve_dap]",
                label=dose_label,
                color=colors[condition],
            )
            # simulation: d3g plasma
            plots[1].add_data(
                task=self.to_sid(f"task_po_dap10_{condition}"),
                xid="time",
                yid="[Cve_d3g]",
                label=dose_label,
                color=colors[condition],
            )
            # data: dapagliflozin plasma
            plots[0].add_data(
                dataset=self.to_sid(f"dapagliflozin_DAP10_{condition}"),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=dose_label,
                color=colors[condition],
            )
            # data: d3g plasma
            plots[1].add_data(
                dataset=self.to_sid(f"dapagliflozin 3-o-glucuronide_DAP10_{condition}"),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=dose_label,
                color=colors[condition],
            )
        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Kasichayanula2011, output_dir=Kasichayanula2011.__name__)
