"""Dapagliflozin pharmacokinetics."""
import pandas as pd

from pkdb_analysis.pk.pharmacokinetics import TimecoursePK
from sbmlsim.result import XResult
from sbmlutils.log import get_logger

logger = get_logger(__name__)


def calculate_dapagliflozin_pk(
    experiment: "DapagliflozinSimulationExperiment",
    xres: XResult,
) -> pd.DataFrame:
    """Calculate PK parameters.

    Only works for 1D-scans.
    Currently only supporting po scans.
    """
    Q_ = experiment.Q_

    # scanned dimension
    scandim = xres._redop_dims()[0]
    dose_vec = Q_(xres["PODOSE_dap"].values[0], xres.uinfo["PODOSE_dap"])

    pk_dicts = []

    for k_dose, dose in enumerate(dose_vec):
        # common time vector
        t_vec = Q_(xres.dim_mean("time").magnitude, xres.uinfo["time"])

        # parent
        c_dap = Q_(
            xres["[Cve_dap]"].sel({scandim: k_dose}).values,
            xres.uinfo["[Cve_dap]"],
        )
        dose_mmole = dose / experiment.Mr.dap
        pk_dap = TimecoursePK(
            time=t_vec,
            concentration=c_dap,
            substance="dapagliflozin",
            dose=dose_mmole,
            ureg=experiment.ureg,
        ).pk.to_dict()
        pk_dap["substance"] = "dap"
        pk_dicts.append(pk_dap)

        # metabolite
        c_d3g = Q_(
            xres["[Cve_d3g]"].sel({scandim: k_dose}).values,
            xres.uinfo["[Cve_d3g]"],
        )
        pk_d3g = TimecoursePK(
            time=t_vec,
            concentration=c_d3g,
            substance="dapagliflozin-3-o-glucuronide",
            dose=None,
            ureg=experiment.ureg,
        ).pk.to_dict()
        pk_d3g["substance"] = "d3g"
        pk_dicts.append(pk_d3g)

        # total = parent + metabolite
        conc_unit = c_dap.units
        c_d3g_same = c_d3g.to(conc_unit)
        c_tot = Q_(c_dap.magnitude + c_d3g_same.magnitude, conc_unit)

        pk_tot = TimecoursePK(
            time=t_vec,
            concentration=c_tot,
            substance="total dapagliflozin (dap + d3g)",
            dose=dose_mmole,
            ureg=experiment.ureg,
        ).pk.to_dict()
        pk_tot["substance"] = "daptot"
        pk_dicts.append(pk_tot)

    return pd.DataFrame(pk_dicts)
