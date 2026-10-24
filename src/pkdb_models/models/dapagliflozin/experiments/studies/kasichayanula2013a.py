"""Simulation experiment of the study Kasichayanula2013a."""

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


class Kasichayanula2013a(DapagliflozinSimulationExperiment):
    """Simulation experiment of Kasichayanula2013a.

    Two cohorts of healthy subjects (fasted, body weight not reported). RIF
    (n=14): 10 mg dapagliflozin on day 1, rifampicin 600 mg once daily on days
    4-11 and 10 mg dapagliflozin with rifampicin on day 9. MFA (n=16): 10 mg
    dapagliflozin on day 1, mefenamic acid 500 mg and then 250 mg every 6 hr
    from day 4 and 10 mg dapagliflozin with mefenamic acid on day 5. Plasma
    0-12 hr (Fig1), urine 0-72 hr (Table 1).

    The interactions are the UGT1A9 activities of the base class
    (`interactions`), set against the observed AUC ratios of dapagliflozin of
    this study (rifampicin 0.78, mefenamic acid 1.51), never fitted: the
    interaction arms are validation data. The dose of day 1 is simulated as a
    separate single dose (the residual of day 1 on day 5 is 0.4 ng/ml).
    Mefenamic acid has no member of the coadministration vocabulary of
    pkpdsim, its arm is `Coadministration.COMEDICATION`.
    """

    groups: ClassVar[list[str]] = ["RIF", "MFA"]
    #: the coadministered drug of the cohort, a key of `interactions`
    drugs: ClassVar[dict[str, str]] = {
        "RIF": "rifampicin",
        "MFA": "mefenamic acid",
    }
    coadministrations: ClassVar[dict[str, Coadministration]] = {
        "RIF": Coadministration.RIFAMPICIN,
        "MFA": Coadministration.COMEDICATION,
    }
    colors: ClassVar[dict[str, str]] = {
        "RIF": "tab:blue",
        "MFA": "tab:orange",
    }
    info: ClassVar[list[tuple[str, str]]] = [
        ("[Cve_dap]", "dapagliflozin"),
        ("Aurine_dap", "dapagliflozin_cumulative amount"),
        ("[Cve_d3g]", "dapagliflozin 3-o-glucuronide"),
        ("Aurine_d3g", "dapagliflozin 3-o-glucuronide_cumulative amount"),
    ]

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Tab1A"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                # unit conversion to mole/l and mole; the urinary glucuronide is
                # the recovery in % of the dose times the dose, i.e., mg of
                # dapagliflozin equivalents
                if label.startswith(
                    ("dapagliflozin_", "dapagliflozin 3-o-glucuronide_cumulative")
                ):
                    dset.unit_conversion("mean", 1 / self.Mr.dap)
                elif label.startswith("dapagliflozin 3-o-glucuronide"):
                    dset.unit_conversion("mean", 1 / self.Mr.d3g)
                dsets[self.to_sid(label)] = dset
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for sim_id, changes in [
            ("po_dap10", {}),
            *[
                (f"po_dap10_{group}", self.interaction_changes(drug))
                for group, drug in self.drugs.items()
            ],
        ]:
            simulations[sim_id] = Simulation(
                time_unit="hr",
                end=73,
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
                    # interaction
                    **changes,
                    # dose (IVDOSE, PODOSE)
                    "PODOSE_dap": Q(10, "mg"),
                },
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for sid, prefix in self.info:
            tissue = Tissue.URINE if "cumulative" in prefix else Tissue.PLASMA
            for group in self.groups:
                for interaction in [False, True]:
                    if interaction:
                        key = f"fm_po_dap10_{prefix}_{group}_{self.drugs[group]}"
                        dset_id = f"{prefix}_DAP10, {group}_{group}"
                        task = f"task_po_dap10_{group}"
                        coadministration = self.coadministrations[group]
                    else:
                        key = f"fm_po_dap10_{prefix}_{group}"
                        dset_id = f"{prefix}_DAP10_{group}"
                        task = "task_po_dap10"
                        coadministration = Coadministration.NONE
                    mappings[self.to_sid(key)] = FitMapping(
                        self,
                        reference=FitData(
                            self,
                            dataset=self.to_sid(dset_id),
                            xid="time",
                            yid="mean",
                            yid_sd="mean_sd",
                            count="count",
                        ),
                        observable=FitData(
                            self,
                            task=self.to_sid(task),
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
                            coadministration=coadministration,
                        ),
                    )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        figures = {}
        for group in self.groups:
            drug = self.drugs[group]
            fig = Figure(
                experiment=self,
                sid=f"Fig1_{group}",
                num_rows=2,
                num_cols=2,
                name=f"{self.__class__.__name__} (Healthy, {drug})",
            )
            plots = fig.create_plots(
                xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
            )
            plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
            plots[1].set_yaxis(self.label_dap_urine, unit=self.unit_dap_urine)
            plots[2].set_yaxis(self.label_d3g_plasma, unit=self.unit_d3g)
            plots[3].set_yaxis(self.label_d3g_urine, unit=self.unit_d3g_urine)
            for kp, (sid, prefix) in enumerate(self.info):
                for task, dset_id, label, color in [
                    ("task_po_dap10", f"{prefix}_DAP10_{group}", "10 mg PO", "black"),
                    (
                        f"task_po_dap10_{group}",
                        f"{prefix}_DAP10, {group}_{group}",
                        f"10 mg PO + {drug}",
                        self.colors[group],
                    ),
                ]:
                    # simulation
                    plots[kp].add_data(
                        task=self.to_sid(task),
                        xid="time",
                        yid=sid,
                        label=label,
                        color=color,
                    )
                    # data
                    plots[kp].add_data(
                        dataset=self.to_sid(dset_id),
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                        label=label,
                        color=color,
                        linestyle="" if sid in ["Aurine_dap", "Aurine_d3g"] else "--",
                    )
            figures.update(self.figures_by_sid(fig))
        return figures


if __name__ == "__main__":
    MODEL.run_experiments(Kasichayanula2013a, output_dir=Kasichayanula2013a.__name__)
