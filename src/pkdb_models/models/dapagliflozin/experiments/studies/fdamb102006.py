"""Simulation experiment of the study FDAMB102006."""

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


class FDAMB102006(DapagliflozinSimulationExperiment):
    """Simulation experiment of FDAMB102006.

    Study MB102006 of the clinical pharmacology review of NDA 202293, the mass
    balance study of Obermeier2010 (the plasma curves are the experiment
    Obermeier2010): 6 healthy men took 50 mg [14C]-dapagliflozin as an oral
    solution (PEG/water) in the fasted state; urine and feces over 13 days
    (Table 11). The total radioactivity (TRA) in urine (75.2 % of the dose) is
    dapagliflozin (1.2 %), D3G (61 %) and other metabolites, which the model has
    (`Aurine_daptot` = dap + d3g + dapmet). The TRA in feces (21 %) contains
    metabolites besides dapagliflozin (15.4 %) which the model does not have
    (`Afeces_daptot` = `Afeces_dap`): a limitation. The amounts of the copy
    are in µmol of the dose (D3G in molar equivalents).
    """

    colors: ClassVar[dict[str, str]] = {
        "Aurine_daptot": "black",
        "Aurine_dap": "tab:blue",
        "Aurine_d3g": "tab:orange",
        "Afeces_daptot": "black",
        "Afeces_dap": "tab:blue",
    }

    info: ClassVar[dict[str, str]] = {
        "Aurine_daptot": "total dapagliflozin_urine",
        "Aurine_dap": "dapagliflozin_urine",
        "Aurine_d3g": "dapagliflozin 3-o-glucuronide_urine",
        "Afeces_daptot": "total dapagliflozin_feces",
        "Afeces_dap": "dapagliflozin_feces",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Tab11A"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                dsets[self.to_sid(label)] = dset
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        simulations["po_DAP50"] = Simulation(
            time_unit="hr",
            end=14 * 24,
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
                "PODOSE_dap": Q(50, "mg"),
            },
        )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.info.items():
            mappings[self.to_sid(f"fm_DAP50_{name}")] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=self.to_sid(name),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task="task_po_DAP50",
                    xid="time",
                    yid=sid,
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.URINE if "urine" in sid else Tissue.FECES,
                    route=Route.PO,
                    application_form=ApplicationForm.SOLUTION,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
                ),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Tab11",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit="day"), legend=True)
        plots[0].set_yaxis("Excretion urine", unit=self.unit_daptot_urine)
        plots[1].set_yaxis("Excretion feces", unit=self.unit_daptot_feces)

        for sid, name in self.info.items():
            if "urine" in sid:
                kp = 0
            elif "feces" in sid:
                kp = 1
            # labels
            if sid == "Aurine_daptot":
                label = "total, 50 mg PO "
            elif sid == "Aurine_dap":
                label = "dap, 50 mg PO "
            elif sid == "Aurine_d3g":
                label = "d3g, 50 mg PO "
            elif sid == "Afeces_daptot":
                label = "total, 50 mg PO "
            elif sid == "Afeces_dap":
                label = "dap, 50 mg PO "
            else:
                label = name

            # simulation
            plots[kp].add_data(
                task="task_po_DAP50",
                xid="time",
                yid=sid,
                label=label,
                color=self.colors[sid],
            )
            # data
            plots[kp].add_data(
                dataset=self.to_sid(name),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                label=label,
                color=self.colors[sid],
                linestyle="",
            )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(FDAMB102006, output_dir=FDAMB102006.__name__)
