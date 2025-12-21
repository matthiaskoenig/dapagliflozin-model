"""Parameter fit problems for dapagliflozin."""
from typing import Dict, List
from sbmlsim.fit.helpers import f_fitexp, filter_empty
from sbmlutils.console import console
from sbmlutils.log import get_logger

from sbmlsim.fit import FitExperiment, FitMapping

from pkdb_models.models.dapagliflozin import DAPAGLIFLOZIN_PATH, DATA_PATHS
from pkdb_models.models.dapagliflozin.experiments.metadata import (
    Tissue, Route, Dosing, ApplicationForm, Health, Coadministration,
    Fasting, DapagliflozinMappingMetaData
)
from pkdb_models.models.dapagliflozin.experiments.studies import *


logger = get_logger(__name__)


# --- Experiment classes ---
experiment_classes = [
    Boulton2013,
    Cho2021,
    FDAMB102002,
    FDAMB102003,
    FDAMB102006,
    FDAMB102007,
    Gould2013,
    Hwang2022a,
    Imamura2013,
    Jang2020,
    Kasichayanula2011,
    Kasichayanula2011a,
    Kasichayanula2011b,
    Kasichayanula2011c,
    Kasichayanula2012,
    Kasichayanula2013,
    Kasichayanula2013a,
    Khomitskaya2018,
    Kim2023,
    Kim2023a,
    Komoroski2009,
    LaCreta2016,
    Obermeier2010,
    Sha2015,
    Shah2019a,
    vanderAartvanderBeek2020,
    Watada2019,
    Yang2013
]

# --- Filters ---
def filter_control(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """Return control experiments/mappings."""

    metadata: DapagliflozinMappingMetaData = fit_mapping.metadata

    # only PO and IV (no SL, MU, RE)
    if metadata.route not in {Route.PO, Route.IV}:
        return False

    # filter multiple dosing (only single dosing)
    # if metadata.dosing == Dosing.MULTIPLE:
    #     return False

    # no coadministration
    if metadata.coadministration != Coadministration.NONE:
        return False

    # filter health (no renal, cardiac impairment, hepatic impairment, ...)
    if metadata.health not in {Health.HEALTHY, Health.T1DM, Health.T2DM, Health.HYPERTENSION}:
        return False

    # only fasted subjects
    if metadata.fasting not in {Fasting.FASTED, Fasting.NR}:
        return False

    # remove outliers
    if metadata.outlier is True:
        return False

    return True


def filter_iv(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """Only iv data."""
    return fit_mapping.metadata.route == Route.IV


def filter_po(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """Only po data."""
    return fit_mapping.metadata.route == Route.PO


def filter_dap(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """Only dap data."""
    yid = "__".join(fit_mapping.observable.y.sid.split("__")[1:])
    if yid not in {
        "Cve_dap",
        "Aurine_dap",
        "Afeces_dap",
        }:
        return False
    return True


def filter_d3g(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """Only d3g data."""
    yid = "__".join(fit_mapping.observable.y.sid.split("__")[1:])
    if yid not in {
        "Cve_d3g",
        "Aurine_d3g",
        }:
        return False
    return True

def filter_pharmacokinetics(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """Only d3g data."""
    yid = "__".join(fit_mapping.observable.y.sid.split("__")[1:])
    if yid not in {
        "Cve_dap",
        "Cve_d3g",
        "Cve_daptot",

        "Aurine_dap",
        "Aurine_d3g",
        "Aurine_daptot",

        "Afeces_dap",
        "Afeces_daptot",
    }:
        return False
    return True

def filter_pharmacodynamics(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """Only pharmacodynamics data."""
    yid = "__".join(fit_mapping.observable.y.sid.split("__")[1:])
    if yid not in {
        "KI__RTG",
        "KI__UGE"
    }:
        return False

    return True

kwargs_f_fitexp = dict(
    experiment_classes=experiment_classes,
    base_path=DAPAGLIFLOZIN_PATH,
    data_path=DATA_PATHS,
)

# --- Fit experiments ---
def f_fitexp_all():
    """All data."""
    return f_fitexp(metadata_filters=filter_empty, **kwargs_f_fitexp)


def f_fitexp_control() -> Dict[str, List[FitExperiment]]:
    """Control data."""
    return f_fitexp(metadata_filters=filter_control, **kwargs_f_fitexp)


def f_fitexp_pharmacokinetics() -> Dict[str, List[FitExperiment]]:
    """Pharmacodynamics data."""
    return f_fitexp(metadata_filters=[filter_control, filter_pharmacokinetics], **kwargs_f_fitexp)


def f_fitexp_pharmacodynamics() -> Dict[str, List[FitExperiment]]:
    """Pharmacodynamics data."""
    return f_fitexp(metadata_filters=[filter_control, filter_pharmacodynamics], **kwargs_f_fitexp)


if __name__ == "__main__":
    """Test construction of FitExperiments."""

    for f in [
        # f_fitexp_all,
        # f_fitexp_control,
        # f_fitexp_pharmacokinetics,
        f_fitexp_pharmacodynamics,
    ]:
        console.rule(style="white")
        console.print(f"{f.__name__}")
        fitexp = f()
