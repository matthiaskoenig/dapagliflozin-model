from dataclasses import dataclass
from enum import Enum

from sbmlsim.fit.objects import MappingMetaData


class Tissue(str, Enum):
    PLASMA = "plasma"
    SERUM = "serum"
    URINE = "urine"
    FECES = "feces"


class Route(str, Enum):
    PO = "po"
    IV = "iv"


class Dosing(str, Enum):
    SINGLE = "single"
    MULTIPLE = "multiple"
    CONSTANT_INFUSION = "infusion"


class ApplicationForm(str, Enum):
    TABLET = "tablet"
    SOLUTION = "solution"
    CAPSULE = "capsule"
    MIXED = "mixed"  # mix of forms, e.g. po and iv
    NR = "not reported"


class Health(str, Enum):
    HEALTHY = "healthy"
    T1DM = "t1dm"
    T2DM = "t2dm"
    HYPERTENSION = "hypertension"
    CIRRHOSIS = "cirrhosis"  # affects PK
    RENAL_IMPAIRMENT = "renal impairment"  # affects PK
    RENAL_IMPAIRMENT_T2DM = "renal impairment and t2dm"  # affects PK
    CHF = "congestive heart failure"  # affects PK


class Fasting(str, Enum):
    NR = "not reported"
    FASTED = "fasted"
    FED = "fed"


class Coadministration(str, Enum):
    NONE = "none"
    CYCLOSPORINE = "cyclosporine"
    EVOGLIPTIN = "evogliptin"
    GLIMEPIRIDE = "glimepiride"
    HCTZ = "hydrochlorothiazide"
    INSULIN = "insulin"
    LOBEGLITAZONE = "lobeglitazone"
    METFORMIN = "metformin"
    PIOGLITAZONE = "pioglitazone"
    PROBENECID = "probenecid"
    RIFAMPICIN = "rifampicin"
    SIMVASTATIN = "simvastatin"
    SITAGLIPTIN = "siztagliptin"
    VOGLIBOSE = "voglibose"
    VALSARTAN = "valsartan"


@dataclass
class DapagliflozinMappingMetaData(MappingMetaData):
    """Metadata for fitting experiment."""
    tissue: Tissue
    route: Route
    application_form: ApplicationForm
    dosing: Dosing
    health: Health
    fasting: Fasting
    coadministration: Coadministration = Coadministration.NONE
    outlier: bool = False

    def to_dict(self):
        return {
            "tissue": self.tissue.name,
            "route": self.route.name,
            "application_form": self.application_form.name,
            "dosing": self.dosing.name,
            "health": self.health.name,
            "fasting": self.fasting.name,
            "coadministration": self.coadministration.name,
            "outlier": self.outlier,
        }
