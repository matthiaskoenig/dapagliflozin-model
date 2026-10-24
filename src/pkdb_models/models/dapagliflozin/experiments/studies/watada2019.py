"""Simulation experiment of the study Watada2019."""

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


class Watada2019(DapagliflozinSimulationExperiment):
    """Simulation experiment of Watada2019.

    Placebo, 5 and 10 mg dapagliflozin once daily for 7 days as add-on to
    insulin in Japanese patients with type 1 diabetes. Plasma of day 7 (Fig 1,
    60 min before and 0.5 to 24 hr after the dose) and the 24 hr urinary glucose
    of day -1 and day 7 (Tab 3). The urinary glucose is driven by the glucose of
    the day, the mean of the 7-point SMBG (9.5-10.3 mM), not by the FPG
    (7.4-7.9 mM). Fasting before the dose is not reported.
    """

    info: ClassVar[list[tuple[str, str]]] = [
        ("[Cve_dap]", "dapagliflozin"),
        ("[Cve_d3g]", "dapagliflozin 3-o-glucuronide"),
        ("KI__UGE", "uge"),
    ]
    interventions: ClassVar[list[str]] = ["DAP0", "DAP5", "DAP10"]
    doses: ClassVar[dict[str, float]] = {
        "DAP0": 0,
        "DAP5": 5,
        "DAP10": 10,
    }
    bodyweights: ClassVar[dict[str, float]] = {
        "DAP0": 57.2,
        "DAP5": 61.6,
        "DAP10": 59.8,
    }  # [kg]

    markers: ClassVar[dict[str, str]] = {"DAP0": "o", "DAP5": "s", "DAP10": "D"}
    #: mg/dl of glucose in mM
    mgdl = 10 / 180.16
    #: fasting plasma glucose at baseline (Tab 1) [mg/dl]
    fpgs: ClassVar[dict[str, float]] = {
        "DAP0": 134.0,
        "DAP5": 142.9,
        "DAP10": 133.4,
    }
    #: the mean of the 7-point SMBG at baseline (Tab 1) and on day 7 (baseline
    #: and the change of the text) [mg/dl]: the glucose of the day which is
    #: filtered, as `KI__PG_delta` above the FPG
    smbgs_baseline: ClassVar[dict[str, float]] = {
        "DAP0": 178.70,
        "DAP5": 184.32,
        "DAP10": 177.01,
    }
    smbgs_day7: ClassVar[dict[str, float]] = {
        "DAP0": 178.70 + 6.54,
        "DAP5": 184.32 - 1.99,
        "DAP10": 177.01 - 5.68,
    }
    #: eGFR (Tab 1) [ml/min/1.73 m2]
    egfrs: ClassVar[dict[str, float]] = {
        "DAP0": 95.4,
        "DAP5": 91.6,
        "DAP10": 94.6,
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Tab3"]:
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
        for intervention, dose in self.doses.items():
            fpg = self.fpgs[intervention] * self.mgdl
            # the 24 hr urine of day -1 (baseline) without dapagliflozin, then
            # doses on days 1 to 7, the dose of day 7 at time 0; the urine is
            # collected per day from the doses
            simulations[self.to_sid(f"po_{intervention}")] = Simulation(
                time_unit="hr",
                end=24 + 6 * 24 + 48,
                preinit_changes={
                    **self.default_changes(),
                    # physiological changes
                    "BW": Q(self.bodyweights[intervention], "kg"),
                    "[KI__fpg]": Q(fpg, "mM"),
                    "KI__PG_delta": Q(
                        self.smbgs_baseline[intervention] * self.mgdl - fpg, "mM"
                    ),
                    "KI__f_renal_function": Q(
                        self.egfrs[intervention] / self.gfr_healthy, "dimensionless"
                    ),
                    "GU__fed": Q(self.fasting_map["not reported"], "dimensionless"),
                    "f_cirrhosis": Q(self.cirrhosis_map["Control"], "dimensionless"),
                    # dose (IVDOSE, PODOSE)
                    "PODOSE_dap": Q(0, "mg"),
                },
                changes=[
                    Change(
                        [24 + 24 * k for k in range(7)],
                        {"KI__glc_urine": Q(0, "mmole")},  # reset UGE
                    ),
                    self.add_doses(
                        [24 + 24 * k for k in range(7)], {"PODOSE_dap": Q(dose, "mg")}
                    ),
                    # the glucose of the days on treatment
                    Change(
                        [24],
                        {
                            "KI__PG_delta": Q(
                                self.smbgs_day7[intervention] * self.mgdl - fpg,
                                "mM",
                            )
                        },
                    ),
                ],
                time_shift=-(24 + 6 * 24),
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for intervention, dose in self.doses.items():
            for sid, name in self.info:
                if dose == 0 and name != "uge":
                    continue
                tissue = Tissue.URINE if name == "uge" else Tissue.PLASMA
                mappings[self.to_sid(f"fm_po_{intervention}_{name}")] = FitMapping(
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
                        task=self.to_sid(f"task_po_{intervention}"),
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=tissue,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.MULTIPLE,
                        health=Health.T1DM,
                        fasting=Fasting.NR,
                        coadministration=Coadministration.INSULIN,
                    ),
                )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_rows=1,
            num_cols=3,
            name=f"{self.__class__.__name__} (T1DM)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
        plots[1].set_yaxis(self.label_d3g_plasma, unit=self.unit_d3g)
        plots[2].set_yaxis(self.label_uge, unit=self.unit_uge)
        plots[0].set_xaxis(self.label_time, unit=self.unit_time, min=-2, max=26)
        plots[1].set_xaxis(self.label_time, unit=self.unit_time, min=-2, max=26)

        for intervention, dose in self.doses.items():
            for k, (sid, name) in enumerate(self.info):
                # simulation
                if not (dose == 0 and k in [0, 1]):
                    plots[k].add_data(
                        task=self.to_sid(f"task_po_{intervention}"),
                        xid="time",
                        yid=sid,
                        label=f"{dose} mg",
                        color=self.dose_colors[dose],
                    )
                # data
                if dose == 0 and name != "uge":
                    continue
                plots[k].add_data(
                    dataset=self.to_sid(f"{name}_{intervention}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"{dose} mg",
                    linestyle="" if name == "uge" else "--",
                    color=self.dose_colors[dose],
                    marker=self.markers[intervention],
                )
        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Watada2019, output_dir=Watada2019.__name__)
