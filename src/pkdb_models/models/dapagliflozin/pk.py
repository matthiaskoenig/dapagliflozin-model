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

    # calculate dap pharmacokinetic parameters
    pk_dicts = list()
    for k_dose, dose in enumerate(dose_vec):
        t_vec: Q_ = xres.dim_mean("time")
        t_vec = Q_(t_vec.magnitude, xres.uinfo["time"])

        c_vec = Q_(
            xres["[Cve_dap]"].sel({scandim: k_dose}).values,
            xres.uinfo["[Cve_dap]"]
        )

        dose_mmole = dose / experiment.Mr.dap
        tcpk = TimecoursePK(
            time=t_vec,
            concentration=c_vec,
            substance="dapagliflozin",
            dose=dose_mmole,
            ureg=experiment.ureg,
        )

        pk_dict = tcpk.pk.to_dict()
        # print(pk_dict)

        # pk = tcpk.pk
        # # Calculate renal clearance as amount in urine/AUC plasma
        # aurine_vec = Q_(
        #     xres["Aurine_ena"].sel({scandim: k_dose}).values, xres.uinfo["Aurine_ena"]
        # )
        # pk_dict["Aurine_ena"] = aurine_vec.magnitude[-1]
        # pk_dict["Aurine_ena_unit"] = aurine_vec.units
        # cl_renal = aurine_vec[-1] / pk.auc  # [mmole] / [mmole/l*min] = [l/min]
        # pk_dict["cl_renal"] = cl_renal.magnitude
        # pk_dict["cl_renal_unit"] = cl_renal.units
        #
        # # # Calculate hepatic clearance as amount in cl_total - cl_renal
        # cl_hepatic = pk.cl - cl_renal  # [l/min]
        # pk_dict["cl_hepatic"] = cl_hepatic.magnitude
        # pk_dict["cl_hepatic_unit"] = cl_hepatic.units

        pk_dict["substance"] = "dap"
        pk_dicts.append(pk_dict)

    # calculate d3g pharmacokinetics
    for k_dose, dose in enumerate(dose_vec):
        t_vec: Q_ = xres.dim_mean("time")
        t_vec = Q_(t_vec.magnitude, xres.uinfo["time"])
        c_vec = Q_(
            xres["[Cve_d3g]"].sel({scandim: k_dose}).values,
            xres.uinfo["[Cve_d3g]"]
        )

        tcpk = TimecoursePK(
            time=t_vec,
            concentration=c_vec,
            substance="dapagliflozin-3-o-glucuronide",
            dose=None,
            ureg=experiment.ureg,
        )

        pk_dict = tcpk.pk.to_dict()
        # print(pk_dict)

        # # Calculate renal clearance as amount in urine/AUC plasma
        # aurine_vec = Q_(
        #     xres["Aurine_eat"].sel({scandim: k_dose}).values, xres.uinfo["Aurine_eat"]
        # )
        # pk_dict["Aurine_eat"] = aurine_vec.magnitude[-1]
        # pk_dict["Aurine_eat_unit"] = aurine_vec.units
        # cl_renal = aurine_vec[-1] / pk.auc  # [mmole] / [mmole/l*min] = [l/min]
        # pk_dict["cl_renal"] = cl_renal.magnitude
        # pk_dict["cl_renal_unit"] = cl_renal.units

        pk_dict["substance"] = "d3g"
        pk_dicts.append(pk_dict)

    return pd.DataFrame(pk_dicts)
