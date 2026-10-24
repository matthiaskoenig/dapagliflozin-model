"""Pharmacokinetic parameters of the dapagliflozin scans.

The parameters are the non-compartmental analysis of `pkpdutils`, named like
the variables of its `NCAResult` (e.g. `auc_last`, `auc_inf_obs`, `cmax`,
`tmax`, `thalf`, `lambda_z`, `cl_f`, `vz_f`), with the unit of every parameter
in `<name>_unit`, for dapagliflozin, its glucuronide (D3G) and total
dapagliflozin (dapagliflozin and its metabolites) in plasma. The dose of
dapagliflozin and of total dapagliflozin is the oral dose of dapagliflozin;
D3G has no dose, so it has no clearance and volume. Dapagliflozin has the
amount excreted in the urine (`Aurine_dap`) and the renal clearance
(`cl_renal`). The concentrations below `LLOQ_FRACTION` of the maximum are below
the limit of quantification, so that the noise of the integrator is not the
terminal phase.

A scan has one dimension, the scanned parameter; a run over the draws of the
fitted parameters (`run_simulation pkdb_models.models.dapagliflozin
--uncertainty`) adds the dimension of the draws. The parameters are calculated
for every point of every dimension, a row per point with its position along
each dimension (the columns `point` and, with draws, `draw`), see `summarize`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from pkpdsim.pk import nca_frame, scan_nca
from pkpdsim.scan import end_values, scan_dimensions, start_values
from pkpdutils import NCAResult, Route, ureg
from sbmlsim import Q
from sbmlsim.result import ScanResult

if TYPE_CHECKING:
    from pkdb_models.models.dapagliflozin.experiments.base_experiment import (
        DapagliflozinSimulationExperiment,
    )

#: oral dose of the scans
DOSE = "PODOSE_dap"
#: the substances by id of the concentration in plasma (`[Cve_<id>]`)
SUBSTANCES: dict[str, str] = {
    "dap": "dapagliflozin",
    "d3g": "dapagliflozin-3-O-glucuronide",
    "daptot": "total dapagliflozin (dapagliflozin and metabolites)",
}
#: the substances whose dose is the oral dose of dapagliflozin
DOSED_SUBSTANCES: set[str] = {"dap", "daptot"}
#: cumulative amount of dapagliflozin excreted in the urine
URINE = "Aurine_dap"
#: the limit of quantification relative to the maximal concentration
LLOQ_FRACTION = 1e-4


def _dimensions(res: ScanResult, across: str | None) -> list[str]:
    """Get the dimension of a scan and the draws `across`, in the order of the rows.

    Raises:
        ValueError: if the scan has not one dimension besides the draws.
    """
    dims = scan_dimensions(res, across)
    if len(dims) != 1:
        raise ValueError(
            f"A scan of one dimension is required, not of {dims} (draws: {across})."
        )
    return [dims[0], *([across] if across is not None else [])]


def _analysis(
    experiment: DapagliflozinSimulationExperiment, res: ScanResult, sid: str
) -> NCAResult:
    """Analyse a substance at every point of a scan at once.

    Args:
        experiment: experiment of the scan, for the molecular weight.
        res: result of the scan.
        sid: id of the substance, a key of `SUBSTANCES`.

    Returns:
        The analysis over the dimensions of the scan.
    """
    concentration = f"[Cve_{sid}]"
    if sid in DOSED_SUBSTANCES:
        # the dose in the substance unit of the concentration, which `cl_f`
        # and `vz_f` require
        dose = Q(start_values(res, DOSE).values, res.units[DOSE]) / experiment.Mr.dap
        return scan_nca(
            res,
            concentration,
            dose=dose.to("mmole"),
            route=Route.ORAL,
            substance=SUBSTANCES[sid],
            lloq_fraction=LLOQ_FRACTION,
        )
    # the parameters of the curve alone
    return scan_nca(
        res,
        concentration,
        dose=None,
        route=None,
        substance=SUBSTANCES[sid],
        lloq_fraction=LLOQ_FRACTION,
    )


def calculate_dapagliflozin_pk(
    experiment: DapagliflozinSimulationExperiment,
    res: ScanResult,
    across: str | None = None,
) -> pd.DataFrame:
    """Calculate the pharmacokinetic parameters of every simulation of a scan.

    Every substance is analysed at all points of the scan at once
    (`scan_nca`). The renal clearance of dapagliflozin is the amount excreted
    in the urine until the end of the simulation over `auc_last`.

    Args:
        experiment: experiment of the scan.
        res: result of the scan, a scan of the oral dose.
        across: the dimension of the draws of the fitted parameters, `None`
            for a scan without draws.

    Returns:
        One row per simulation of the scan and substance, with its position
        (`point`, and `draw` with draws), the id of the substance in
        `substance`, the parameters and their units (`<name>_unit`), the
        substances of a simulation one after the other.
    """
    names = dict(zip(_dimensions(res, across), ["point", "draw"], strict=False))
    positions = list(names.values())
    dfs = []
    for sid in SUBSTANCES:
        df = nca_frame(_analysis(experiment, res, sid)).drop(columns="nca_flags")
        if sid == "dap":
            position = tuple(df[dim].to_numpy() for dim in res.dims)
            aurine = ureg.Quantity(
                end_values(res, URINE).values[position], str(res.units[URINE])
            )
            auc_last = ureg.Quantity(
                df["auc_last"].to_numpy(), df["auc_last_unit"].iloc[0]
            )
            for key, value in {
                URINE: aurine,
                "cl_renal": (aurine / auc_last).to("ml/min"),
            }.items():
                df[key] = value.magnitude
                df[f"{key}_unit"] = str(value.units)
        df = df.rename(columns=names)
        df.insert(len(positions), "substance", sid)
        dfs.append(df[[*positions, *df.columns.drop(positions)]])
    df = pd.concat(dfs, ignore_index=True)
    return df.sort_values(positions, kind="stable").reset_index(drop=True)


def summarize(
    df: pd.DataFrame, key: str, quantiles: tuple[float, float]
) -> pd.DataFrame:
    """Get the median and the quantiles of a parameter over the draws of every point.

    Args:
        df: the parameters with the columns `point` and `draw`.
        key: the parameter.
        quantiles: the lower and the upper quantile.

    Returns:
        One row per point (the index), the columns `lower`, `median`, `upper`.
    """
    low, high = quantiles
    grouped = df.groupby("point")[key]
    return pd.DataFrame(
        {
            "lower": grouped.quantile(low),
            "median": grouped.median(),
            "upper": grouped.quantile(high),
        }
    ).astype(np.float64)
