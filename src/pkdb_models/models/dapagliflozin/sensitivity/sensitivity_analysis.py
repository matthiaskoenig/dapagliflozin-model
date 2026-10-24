"""Sensitivity analysis of the dapagliflozin model.

The outputs are the pharmacokinetic parameters of dapagliflozin and its
glucuronide (D3G) in plasma after a single oral dose of 10 mg (the
non-compartmental analysis of `pkpdutils`) and the urinary glucose excretion
(UGE) in the 24 hr after the dose. The parameters are the parameters of the
model with bounds of +-15 % around their value. The analyses are run as a
script, the local analysis by default:

    python -m pkdb_models.models.dapagliflozin.sensitivity.sensitivity_analysis
    python -m pkdb_models.models.dapagliflozin.sensitivity.sensitivity_analysis local sampling

The analyses run on `pkpdsim.sensitivity.SensitivityRunner`, which writes the
table of the parameters into `results/sensitivity` of the model and the tables
and figures of an analysis into `results/sensitivity/<analysis>`.
An analysis does not vary a parameter which one of its groups sets
(`pkpdsim.sensitivity.free_parameters`): `sampling`, over the seven groups,
leaves out `KI__f_renal_function` (`f_cirrhosis`, 0 in the model, is no
parameter of the analyses). The environment variable `PKPDSIM_SENSITIVITY_N`
replaces the numbers of samples of the global analyses for a quick run, see
`pkpdsim.sensitivity.sample_size`.

The model without dapagliflozin is stationary, so the simulation of a group
starts with the dose. The dose is followed for `END_HR` hours, long enough for
the terminal phase of dapagliflozin (a half-life of about 12 hr, and the slow
pool of the rest tissue).
"""

from collections.abc import Mapping

import numpy as np
from pkpdsim.sensitivity import (
    FastAnalysis,
    LocalAnalysis,
    ParameterType,
    SamplingAnalysis,
    SensitivityGroup,
    SensitivityOutput,
    SensitivityParameter,
    SensitivityRunner,
    SobolAnalysis,
    parameters_from_sbml,
    set_bounds,
    set_relative_bounds,
)
from pkpdutils import NCAOptions
from sbmlsim import Q
from sbmlsim.simulation import PK, Formula, Observable, Simulation

from pkdb_models.models.dapagliflozin import MODEL, MODEL_PATH
from pkdb_models.models.dapagliflozin.experiments.base_experiment import (
    DapagliflozinSimulationExperiment,
)
from pkdb_models.models.dapagliflozin.fitting.parameters import (
    parameters_pk_pd as fit_parameters,
)

#: oral dose of the analysis [mg]
DOSE_DAP = 10.0
#: the hours after the dose
END_HR = 168
#: the limit of quantification of the non-compartmental analysis [mM], about
#: 1e-4 of the maximal plasma concentration after 10 mg (3e-4 mM), so that the
#: noise of the integrator is not the terminal phase
LLOQ = 3e-8
#: fraction of the bounds of a parameter relative to its value
BOUNDS_FRACTION = 0.15

#: outputs of the non-compartmental analysis -> name and unit, by id of the
#: concentration in plasma (`[Cve_<id>]`); the glucuronide has no dose, so no
#: clearance and volume
PK_OUTPUTS: dict[str, dict[str, tuple[str, str]]] = {
    "dap": {
        "auc_inf_obs": ("AUC∞", "mM*min"),
        "cmax": ("Cmax", "mM"),
        "tmax": ("tmax", "min"),
        "thalf": ("Half-life", "min"),
        "vz_f": ("Vz/F", "l"),
        "cl_f": ("CL/F", "ml/min"),
    },
    "d3g": {
        "auc_inf_obs": ("AUC∞", "mM*min"),
        "cmax": ("Cmax", "mM"),
        "tmax": ("tmax", "min"),
        "thalf": ("Half-life", "min"),
    },
}
#: the names of the substances
SUBSTANCES: dict[str, str] = {"dap": "DAP", "d3g": "D3G"}
#: the substance with a dose, the glucuronide has none
DOSED = "dap"
#: the outputs whose unit has the amount of the dose: the dose is a mass
#: (`PODOSE_dap` in mg), so they are divided by the molecular weight (its value
#: at the start, a value per simulation like the parameter of the analysis)
PK_PER_DOSE_MASS: frozenset[str] = frozenset({"vz_f", "cl_f"})
#: the time of the urinary glucose excretion [min], in the time unit of the
#: model (`at` of sbmlsim)
UGE_TIME = 24 * 60

#: the outputs with their ids and names
OUTPUTS: list[SensitivityOutput] = [
    *(
        SensitivityOutput(f"{sid}_{key}", f"{SUBSTANCES[sid]} {name}")
        for sid, outputs in PK_OUTPUTS.items()
        for key, (name, _) in outputs.items()
    ),
    SensitivityOutput("uge24", "UGE (24h)"),
]

#: the observables of a simulation, the outputs with their units; the id of the
#: analysis of a substance (`pk_<id>`) is no species of the model
OBSERVABLES: list[Observable] = [
    *(
        PK(
            f"pk_{sid}",
            f"[Cve_{sid}]",
            dose="PODOSE_dap" if sid == DOSED else None,
            route="oral" if sid == DOSED else None,
            options=NCAOptions(lloq=LLOQ),
            parameters=list(outputs),
        )
        for sid, outputs in PK_OUTPUTS.items()
    ),
    *(
        Formula(
            f"{sid}_{key}",
            f"pk_{sid}.{key} / at(Mr_dap, 0)"
            if key in PK_PER_DOSE_MASS
            else f"pk_{sid}.{key}",
            unit=unit,
        )
        for sid, outputs in PK_OUTPUTS.items()
        for key, (_, unit) in outputs.items()
    ),
    # the pharmacodynamics: the urinary glucose excretion in 24 hr
    Formula("uge24", f"at(KI__UGE, {UGE_TIME})", unit="g"),
]

#: the output times [hr]: every 0.02 hr in the first day after the dose
#: (absorption, maximum and the UGE at 24 hr), every 0.5 hr after it
TIMES_HR: tuple[float, ...] = tuple(
    np.unique(
        np.concatenate(
            [
                np.linspace(0, 24, 1201),
                np.linspace(24, END_HR, 2 * (END_HR - 24) + 1),
            ]
        )
    ).tolist()
)


def sensitivity_simulation(changes: Mapping[str, float]) -> Simulation:
    """Get the simulation of a group: a single oral dose of dapagliflozin at time 0."""
    return Simulation(
        start=0,
        end=END_HR,
        times=TIMES_HR,
        time_unit="hr",
        preinit_changes={**changes, "PODOSE_dap": Q(DOSE_DAP, "mg")},
    )


#: the groups the analyses are run for: the control and the stages of renal
#: and hepatic impairment of the simulation experiments
SENSITIVITY_GROUPS: list[SensitivityGroup] = [
    SensitivityGroup("control", "Control", sensitivity_simulation({}), "dimgrey"),
    *(
        SensitivityGroup(
            uid,
            name,
            sensitivity_simulation(
                {
                    "KI__f_renal_function": (
                        DapagliflozinSimulationExperiment.renal_map[name]
                    )
                }
            ),
            DapagliflozinSimulationExperiment.renal_colors[name],
        )
        for uid, name in [
            ("mildRI", "Mild renal impairment"),
            ("modRI", "Moderate renal impairment"),
            ("sevRI", "Severe renal impairment"),
        ]
    ),
    *(
        SensitivityGroup(
            uid,
            name,
            sensitivity_simulation(
                {"f_cirrhosis": DapagliflozinSimulationExperiment.cirrhosis_map[name]}
            ),
            DapagliflozinSimulationExperiment.cirrhosis_colors[name],
        )
        for uid, name in [
            ("mildHI", "Mild cirrhosis"),
            ("modHI", "Moderate cirrhosis"),
            ("sevHI", "Severe cirrhosis"),
        ]
    ),
]


def sensitivity_parameters() -> list[SensitivityParameter]:
    """Get the parameters of the analyses with their bounds."""
    parameters = parameters_from_sbml(
        MODEL_PATH,
        exclude_ids={
            # the reference body weight of the tissue exchange
            "BW_ref",
            # unchangeable values: the whole cardiac output passes the lungs,
            # the volumes of the portal and hepatic veins are no parameters of
            # the distribution
            "FQlu",
            "FVhv",
            "FVpo",
            # conversion factors
            "KI__cf_mg_per_g",
            # molecular weights
            "Mr_dap",
            "KI__Mr_glc",
            # dosing parameters
            "PODOSE_dap",
            "IVDOSE_dap",
            "ti_dap",
            # unused volumes
            "Vurine",
            "Vfeces",
            "GU__Vstomach",
        },
        exclude_na=True,
        exclude_zero=True,
    )

    # the fitted parameters, their bounds are set relative to the value below
    set_bounds(
        parameters,
        [(fp.pid, np.nan, np.nan, ParameterType.FIT) for fp in fit_parameters],
    )
    # the scaling factors of the organ functions
    set_bounds(
        parameters,
        [
            (uid, 1 - BOUNDS_FRACTION, 1 + BOUNDS_FRACTION, ParameterType.SCALING)
            for uid in [
                "f_cardiac_function",
                "KI__f_renal_function",
                "f_ugt1a9",
            ]
        ],
    )

    # every parameter without bounds varies around its value
    set_relative_bounds(parameters, BOUNDS_FRACTION)
    return parameters


#: the analyses of the script, with their groups and numbers of samples
SENSITIVITY = SensitivityRunner(
    model=MODEL,
    parameters=sensitivity_parameters,
    outputs=OUTPUTS,
    observables=OBSERVABLES,
    analyses={
        "local": LocalAnalysis(groups=SENSITIVITY_GROUPS[:1], delta=0.01),
        "sampling": SamplingAnalysis(groups=SENSITIVITY_GROUPS, n=1000),
        "fast": FastAnalysis(groups=SENSITIVITY_GROUPS[:1], n=1000),
        "sobol": SobolAnalysis(groups=SENSITIVITY_GROUPS[:1], n=4096),
    },
    default=("local",),
)


if __name__ == "__main__":
    SENSITIVITY.main()
