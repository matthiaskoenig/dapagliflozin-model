"""Simulation experiment of the study Kasichayanula2011c."""

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


class Kasichayanula2011c(DapagliflozinSimulationExperiment):
    """Simulation experiment of Kasichayanula2011c.

    Three crossover studies in healthy subjects (body weight not reported),
    single doses after a fast of 10 hr: S1 50 mg with or without pioglitazone
    45 mg, S2 20 mg with or without metformin 1000 mg, S3 20 mg alone, with
    glimepiride 4 mg or with sitagliptin 100 mg; plasma 0-72 hr (Fig3, mean and
    error bars of an unstated type, curated as SD). In S3 every dose was taken
    with 240 ml of a 20 % glucose solution against hypoglycaemia, which may
    delay the absorption (tmax 1.5 hr); it is simulated fasted, the glucose
    acts only on the glucose excretion, which was not measured. None of the
    drugs changes the exposure (AUC ratios 0.99-1.08, Cmax 0.93-1.09): the
    combinations are the simulation of the dose alone.
    """

    group_doses: ClassVar[dict[str, int]] = {
        "S1": 50,
        "S2": 20,
        "S3": 20,
    }
    group_interventions: ClassVar[dict[str, list[str]]] = {
        "S1": ["DAP50", "DAP50, PIO45"],
        "S2": ["DAP20", "DAP20, MET1000"],
        "S3": ["DAP20", "DAP20, GLI4", "DAP20, SIT100"],
    }
    intervention_colors: ClassVar[dict[str, str]] = {
        "DAP50": "black",
        "DAP50, PIO45": "tab:blue",
        "DAP20": "black",
        "DAP20, MET1000": "tab:blue",
        "DAP20, GLI4": "tab:blue",
        "DAP20, SIT100": "tab:orange",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig3"]:
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
        for dose in set(self.group_doses.values()):
            simulations[self.to_sid(f"po_dap{dose}")] = Simulation(
                time_unit="hr",
                end=75,
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
        for group, interventions in self.group_interventions.items():
            dose = self.group_doses[group]
            for intervention in interventions:
                mappings[self.to_sid(f"fm_po_{intervention}_{group}")] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=self.to_sid(f"dapagliflozin_{intervention}_{group}"),
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=self.to_sid(f"task_po_dap{dose}"),
                        xid="time",
                        yid="[Cve_dap]",
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        coadministration=(
                            Coadministration.PIOGLITAZONE
                            if "PIO" in intervention
                            else Coadministration.METFORMIN
                            if "MET" in intervention
                            else Coadministration.GLIMEPIRIDE
                            if "GLI" in intervention
                            else Coadministration.SITAGLIPTIN
                            if "SIT" in intervention
                            else Coadministration.NONE
                        ),
                    ),
                )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig3",
            num_rows=1,
            num_cols=3,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        for idx, (group, interventions) in enumerate(self.group_interventions.items()):
            plots[idx].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
            dose = self.group_doses[group]
            dose_label = f"{dose} mg PO"
            # simulation
            plots[idx].add_data(
                task=self.to_sid(f"task_po_dap{dose}"),
                xid="time",
                yid="[Cve_dap]",
                label=dose_label,
                color="black",
            )
            # data
            for intervention in interventions:
                if intervention == f"DAP{dose}":
                    label = dose_label
                else:
                    label = f"{dose} mg PO + {intervention.split(', ')[1]}"
                plots[idx].add_data(
                    dataset=self.to_sid(f"dapagliflozin_{intervention}_{group}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=label,
                    color=self.intervention_colors[intervention],
                )
        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Kasichayanula2011c, output_dir=Kasichayanula2011c.__name__)
