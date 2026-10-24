"""Simulation experiment of the study Hwang2022a."""

from typing import ClassVar, override

from pkpdsim.metadata import (
    ApplicationForm,
    Coadministration,
    Dosing,
    Fasting,
    Genotype,
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


class Hwang2022a(DapagliflozinSimulationExperiment):
    """Simulation experiment of Hwang2022a.

    100 healthy Korean subjects (94 men) of four crossover trials, a single
    10 mg tablet of dapagliflozin together with 100 mg sitagliptin, fasting not
    reported; plasma 0-48 hr. Fig2 gives the mean plasma concentrations by the
    genotypes of the ABCB1 SNPs C1236T, G2677T/A and C3435T (no error bars),
    each SNP a partition of the same 100 subjects. Sitagliptin does not change
    the exposure of dapagliflozin (AUC ratio 1.08, Cmax 0.96,
    Kasichayanula2011c) and is not simulated, neither is the genotype (only
    3435TT had a higher AUC, 1.17, p=0.049). The C3435T partition is the data
    of the study; the C1236T and G2677T/A curves count the subjects again
    (aggregates). Fig1 (C3435T CC+CT vs TT) shows the same subjects again and
    is not used.
    """

    bodyweights: ClassVar[dict[str, float]] = {
        "C1236T_CC": 71.03,
        "C1236T_CT": 71.79,
        "C1236T_TT": 72.44,
        "G2677TA_GG": 71.18,
        "G2677TA_GT": 71,
        "G2677TA_GA": 73.61,
        "G2677TA_TT": 71.43,
        "G2677TA_TA": 72.78,
        "G2677TA_AA": 67.05,
        "C3435T_CC": 71.13,
        "C3435T_CT": 73.44,
        "C3435T_TT": 68.09,
    }  # [kg]

    groups: ClassVar[list[str]] = list(bodyweights.keys())

    #: the SNPs of the panels of Fig2
    snps: ClassVar[tuple[str, ...]] = ("C1236T", "G2677TA", "C3435T")

    genotypes: ClassVar[dict[str, Genotype]] = {
        "C3435T_CC": Genotype.ABCB1_3435CC,
        "C3435T_CT": Genotype.ABCB1_3435CT,
        "C3435T_TT": Genotype.ABCB1_3435TT,
    }

    colors: ClassVar[dict[str, str]] = {
        "CC": "tab:blue",
        "CT": "tab:orange",
        "TT": "tab:red",
        "GG": "tab:blue",
        "GT": "tab:orange",
        "GA": "tab:green",
        "TA": "tab:purple",
        "AA": "tab:brown",
    }

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
                dsets[label] = dset
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for group in self.groups:
            simulations[f"po_dap10_{group}"] = Simulation(
                time_unit="hr",
                end=50,
                preinit_changes={
                    # assuming no changes due to ABCB1 genotype
                    **self.default_changes(),
                    # physiological changes
                    "BW": Q(self.bodyweights[group], "kg"),
                    "[KI__fpg]": Q(self.fpg_healthy, "mM"),
                    "GU__fed": Q(self.fasting_map["not reported"], "dimensionless"),
                    "f_cirrhosis": Q(self.cirrhosis_map["Control"], "dimensionless"),
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
        for group in self.groups:
            mappings[f"fm_dap10_{group}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"dapagliflozin_DAP10_{group}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task=f"task_po_dap10_{group}",
                    xid="time",
                    yid="[Cve_dap]",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.NR,
                    coadministration=Coadministration.SITAGLIPTIN,
                    genotype=self.genotypes.get(group, Genotype.NR),
                ),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2",
            num_rows=1,
            num_cols=3,
            name=f"{self.__class__.__name__} (Healthy, with sitagliptin)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        for k, snp in enumerate(self.snps):
            plots[k].set_title(snp.replace("TA", "T/A"))
            plots[k].set_yaxis(self.label_dap_plasma, unit=self.unit_dap)
            for group in self.groups:
                if not group.startswith(f"{snp}_"):
                    continue
                genotype = group.split("_")[1]
                plots[k].add_data(
                    task=f"task_po_dap10_{group}",
                    xid="time",
                    yid="[Cve_dap]",
                    label=None,
                    color=self.colors[genotype],
                )
                plots[k].add_data(
                    dataset=f"dapagliflozin_DAP10_{group}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=genotype,
                    color=self.colors[genotype],
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Hwang2022a, output_dir=Hwang2022a.__name__)
