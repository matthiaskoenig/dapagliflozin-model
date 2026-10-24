"""Simulation experiment of the study Obermeier2010."""

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


class Obermeier2010(DapagliflozinSimulationExperiment):
    """Simulation experiment of Obermeier2010.

    Human ADME study: single oral dose of 50 mg [14C]dapagliflozin as a solution
    in 6 fasted healthy men. Fig 6: plasma dapagliflozin and total radioactivity
    (TRA, ng equivalents of dapagliflozin); Fig 7: dapagliflozin, its
    3-O-glucuronide (M15, quantified by radioactivity, ng equivalents of
    dapagliflozin), their sum and TRA in the samples of the metabolite profiling
    at 1, 4 and 12 hr.
    """

    #: body weight [kg], not reported: CL/F 4.9 ml/min/kg and AUC 2.43 µg*hr/ml
    #: (50 mg / 2.43 µg*hr/ml = 343 ml/min) give 70 kg
    bodyweight = 70.0
    #: the analytes of the plasma: observable, label of the data, key, label
    #: of the figures, color
    analytes: ClassVar[list[tuple[str, str, str, str, str]]] = [
        ("[Cve_dap]", "dapagliflozin", "dap", "dapagliflozin", "black"),
        (
            "[Cve_d3g]",
            "dapagliflozin 3-o-glucuronide",
            "d3g",
            "D3G (M15)",
            "tab:orange",
        ),
        (
            "[Cve_daptot]",
            "total dapagliflozin",
            "daptot",
            "DAP + metabolites (TRA)",
            "tab:blue",
        ),
    ]

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig6", "Fig7"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                # unit conversion to mole/l: the glucuronide (M15) and the
                # total radioactivity are in ng equivalents of dapagliflozin
                dset.unit_conversion("mean", 1 / self.Mr.dap)
                dsets[self.to_sid(f"{fig_id}_{label}")] = dset
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        simulations["po_dap50"] = Simulation(
            time_unit="hr",
            end=30,
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
                "PODOSE_dap": Q(50, "mg"),
            },
        )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        # Fig 6: dapagliflozin and TRA; Fig 7: the analytes of the profiling
        # samples (dapagliflozin and TRA are the subjects of Fig 6 again)
        infos = [
            ("Fig6", "[Cve_dap]", "dapagliflozin", "dap"),
            ("Fig6", "[Cve_daptot]", "total dapagliflozin", "daptot"),
            ("Fig7", "[Cve_dap]", "dapagliflozin", "dap"),
            ("Fig7", "[Cve_d3g]", "dapagliflozin 3-o-glucuronide", "d3g"),
            ("Fig7", "[Cve_daptot]", "total dapagliflozin", "daptot"),
        ]
        for fig_id, sid, name, key in infos:
            suffix = "" if fig_id == "Fig6" else "_profiling"
            mappings[self.to_sid(f"fm_dap50_{key}{suffix}")] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=self.to_sid(f"{fig_id}_{name}_DAP50"),
                    xid="time",
                    yid="mean",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task="task_po_dap50",
                    xid="time",
                    yid=sid,
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.PLASMA,
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
        figures = {}
        for fig_id in ["Fig6", "Fig7"]:
            fig = Figure(
                experiment=self,
                sid=fig_id,
                num_rows=1,
                num_cols=2,
                name=f"{self.__class__.__name__} (Healthy)",
            )
            plots = fig.create_plots(
                xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
            )
            plots[0].set_yaxis("plasma concentration", unit=self.unit_dap)
            plots[1].set_yaxis("plasma concentration", unit=self.unit_dap, scale="log")
            for plot in plots:
                for sid, name, key, label, color in self.analytes:
                    if fig_id == "Fig6" and key == "d3g":
                        continue
                    plot.add_data(
                        task="task_po_dap50",
                        xid="time",
                        yid=sid,
                        label=f"{label}, 50 mg PO",
                        color=color,
                    )
                    plot.add_data(
                        dataset=self.to_sid(f"{fig_id}_{name}_DAP50"),
                        xid="time",
                        yid="mean",
                        count="count",
                        label=label,
                        color=color,
                    )
                if fig_id == "Fig7":
                    # the measured sum of dapagliflozin and M15
                    plot.add_data(
                        dataset=self.to_sid(
                            "Fig7_dapagliflozin + dapagliflozin-3-o-glucuronide_DAP50"
                        ),
                        xid="time",
                        yid="mean",
                        count="count",
                        label="dapagliflozin + D3G",
                        color="tab:green",
                    )
            figures.update(self.figures_by_sid(fig))
        return figures


if __name__ == "__main__":
    MODEL.run_experiments(Obermeier2010, output_dir=Obermeier2010.__name__)
