"""Sensitivity analysis."""
from pathlib import Path

import numpy as np
import xarray as xr
from roadrunner._roadrunner import NamedArray


from sbmlsim.sensitivity.analysis import (
    SensitivitySimulation,
    LocalSensitivityAnalysis, SensitivityOutput,
)

from sbmlsim.sensitivity.parameters import (
    SensitivityParameter,
    parameters_for_sensitivity_analysis,
)
from sbmlutils.console import console
from pkdb_analysis.pk.pharmacokinetics import TimecoursePK
from pint import UnitRegistry


class DapagliflozinSensitivitySimulation(SensitivitySimulation):
    """Simulation for sensitivity calculation."""
    tend = 2 * 24 * 60  # [min]


    def simulate(self, changes: dict[str, float]) -> dict[str, float]:
        # apply changes and simulate
        all_changes = {
            **self.changes_simulation,
            **changes
        }
        self.apply_changes(all_changes, reset_all=True)
        s: NamedArray = self.rr.simulate(start=0, end=self.tend, steps=1000)
        # self._plot(s)

        # pharmacokinetic parameters
        d: dict[str, float] = {}

        # pharmacokinetics
        ureg = UnitRegistry()
        Q_ = ureg.Quantity

        # dapagliflozin
        time = Q_(s["time"], "min")
        for sid in ["dap", "d3g"]:
            tcpk = TimecoursePK(
                time=time,
                concentration=Q_(s[f"[Cve_{sid}]"], "mM"),
                substance="dapagliflozin",
                ureg=ureg,
                dose=None,
                min_treshold=1E4,
            )
            pk_dict = tcpk.pk.to_dict()
            for pk_key in [
                "aucinf",
                "cmax",
                "thalf",
                # "kel",
            ]:
                d[f"{sid}_{pk_key}"] = pk_dict[pk_key]

        # pharmacodynamics
        t_idx = np.argmin(np.abs(time - Q_(24*60, "min")))
        d["uge24"] = s["KI__UGE"][t_idx]

        # console.print(d)

        return d

    def _plot(self, s: NamedArray) -> None:

        # plotting
        from matplotlib import pyplot as plt
        plt.plot(
            s["time"], s["[Cve_dap]"],
            # df["time"], df["[Cve_los]"],
            marker="o",
            markeredgecolor="black",
            color="tab:blue",
        )
        plt.show()


if __name__ == "__main__":
    from pkdb_models.models.dapagliflozin import MODEL_PATH, MODEL_BASE_PATH

    sensitivity_simulation = DapagliflozinSensitivitySimulation(
        model_path=MODEL_PATH,
        selections=[
            "time",
            "[Cve_dap]",
            "[Cve_d3g]",
            "KI__UGE"
        ],
        changes_simulation = {
            "PODOSE_dap": 10.0,  # [mg]
        },
        outputs = [
            SensitivityOutput(uid='dap_aucinf', name='DAP AUC∞'),
            SensitivityOutput(uid='dap_cmax', name='DAP Cmax'),
            SensitivityOutput(uid='dap_thalf', name='DAP Half-life'),
            SensitivityOutput(uid='d3g_aucinf', name='D3G AUC∞'),
            SensitivityOutput(uid='d3g_cmax', name='D3G Cmax'),
            SensitivityOutput(uid='d3g_thalf', name='D3G Half-life'),
            SensitivityOutput(uid='uge24', name='UGE (24h)'),
        ]
    )
    console.print(sensitivity_simulation.outputs)

    # parameters for sensitivity analysis
    parameters: list[SensitivityParameter] = parameters_for_sensitivity_analysis(
        sbml_path=MODEL_PATH,
        exclude_ids={
            # conversion factors
            "conversion_min_per_day",
            "KI__cf_mg_per_g",
            "KI__cf_ml_per_l",

            # molecular weights
            "Mr_dap",  # molecular weight
            "Mr_d3g",  # molecular weight
            "KI__Mr_glc",
            
            "PODOSE_dap", # dose
            "Vurine", # unused volume
            "Vfeces",  # unused volume
            "GU__Vstomach",  # unused volume
        },
        exclude_na=True,
        exclude_zero=True,
    )
    # parameters = [parameters[k] for k in range(5)]

    """ Local Sensitivity Analysis """
    difference = 0.01
    sa = LocalSensitivityAnalysis(
        sensitivity_simulation=sensitivity_simulation,
        parameters=parameters,
        difference=difference,
    )
    sa.create_samples()
    console.print(sa.samples)
    sa.simulate_samples()
    console.print(sa.results)
    sa.calculate_sensitivity()
    console.print(sa.sensitivity)
    sa.plot_sensitivity()
    df = sa.sensitivity_df

    from pkdb_models.models.dapagliflozin import RESULTS_PATH
    sensitivity_dir = RESULTS_PATH / "sensitivity"
    Path(sensitivity_dir).mkdir(parents=True, exist_ok=True)
    sensitivity_fname = f"parameter_local_sensitivity_{difference}"

    df.to_csv(sensitivity_dir / f"{sensitivity_fname}.tsv", sep="\t")

    from matplotlib import pyplot as plt
    plt.savefig(
        sensitivity_dir / f"{sensitivity_fname}.png", dpi=300, bbox_inches="tight"
    )
    plt.show()


    #
    # console.print(parameters)

    #
    # SamplingSensitivityAnalysis
    # sa = GlobalSobolSensitivityAnalysis(
    #     sensitivity_simulation=sensitivity_simulation,
    #     parameters=list(parameters.keys()),
    # )



