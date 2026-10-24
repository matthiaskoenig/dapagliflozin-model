"""Base class of the dapagliflozin simulation experiments."""

from collections.abc import Iterable, Mapping
from typing import ClassVar, NamedTuple

import pandas as pd
from pkpdsim.experiment import PKPDSimulationExperiment
from pkpdsim.impairment import RENAL_COLORS_END_STAGE, RENAL_FUNCTION_GFR
from sbmlsim import Q
from sbmlsim.units import Quantity

from pkdb_models.models.dapagliflozin import MODEL
from pkdb_models.models.dapagliflozin.dapagliflozin_pk import (
    calculate_dapagliflozin_pk,
)


class MolecularWeights(NamedTuple):
    """Molecular weights for the conversion of the data."""

    dap: Quantity
    d3g: Quantity
    glc: Quantity


class DapagliflozinSimulationExperiment(PKPDSimulationExperiment):
    """Base class of the dapagliflozin simulation experiments."""

    pkpd_model = MODEL
    selections = (
        "time",
        # dosing
        "IVDOSE_dap",
        "PODOSE_dap",
        # venous
        "[Cve_daptot]",
        "[Cve_dap]",
        "[Cve_d3g]",
        "[Cve_dapmet]",
        # urine
        "Aurine_daptot",
        "Aurine_dap",
        "Aurine_d3g",
        "Aurine_dapmet",
        # feces
        "Afeces_dap",
        "Afeces_daptot",
        # glucose
        "[KI__fpg]",
        "KI__pg",
        "KI__glc_urine",
        "KI__UGE",
        "KI__RTG",
        "KI__GLCEX",
        # renal excretion rate
        "KI__D3GEX",
        "KI__DAPEX",
        "KI__METEX",
        # cases
        "KI__f_renal_function",
        "KI__GFR",
        "f_cirrhosis",
        "GU__fed",
        "GU__FED_f_absorption",
        "f_ugt1a9",
    )

    font: ClassVar[dict[str, str | int]] = {"weight": "bold", "size": 22}
    scan_font: ClassVar[dict[str, str | int]] = {"weight": "bold", "size": 20}
    tick_font_size = 15

    # healthy reference values
    #: fasting plasma glucose of healthy subjects [mM], `KI__fpg_healthy` of
    #: the model; the filtered glucose is the mean plasma glucose of the day,
    #: `KI__fpg` + `KI__PG_delta` (0.85 mM)
    fpg_healthy = 5.0
    #: fasting plasma glucose of type 2 diabetes if a study reports none [mM]
    fpg_t2dm = 7.5
    #: GFR of healthy kidneys at 75 kg [ml/min], `GFR_base` of the model; a
    #: measured creatinine clearance or GFR of X ml/min is the renal function
    #: X/101
    gfr_healthy = 101.0
    #: body weight if a study reports none [kg]
    bodyweight_default = 75

    # labels
    label_time = "time"

    label_dap = "dapagliflozin"
    label_d3g = "dapagliflozin-3-O-\nglucuronide"
    label_dapmet = "other metabolites"
    label_daptot = "DAP + metabolites"

    label_dap_plasma = f"{label_dap}\nplasma"
    label_d3g_plasma = f"{label_d3g} plasma"
    label_daptot_plasma = f"{label_daptot}\nplasma"

    label_glc_plasma = "glucose plasma"
    label_pg = "mean glucose plasma"
    label_uge = "UGE"
    label_rtg = "RTG"
    label_dap_urine = label_dap + "\nurine"
    label_d3g_urine = label_d3g + " urine"
    label_daptot_urine = label_daptot + "\nurine"

    label_glc_urine = "glucose urine"
    label_dap_feces = label_dap + "\nfeces"
    label_daptot_feces = label_daptot + "\nfeces"

    labels: ClassVar[Mapping[str, str]] = {
        "time": "time",
        "[Cve_dap]": label_dap_plasma,
        "[Cve_d3g]": label_d3g_plasma,
        "[Cve_daptot]": label_daptot_plasma,
        "[Cve_dapmet]": f"{label_dapmet}\nplasma",
        "[KI__fpg]": label_glc_plasma,
        "KI__pg": label_pg,
        "Aurine_dap": label_dap_urine,
        "Aurine_d3g": label_d3g_urine,
        "Aurine_dapmet": f"{label_dapmet}\nurine",
        "Aurine_daptot": label_daptot_urine,
        "KI__glc_urine": label_glc_urine,
        "KI__UGE": label_uge,
        "KI__RTG": label_rtg,
        "Afeces_dap": label_dap_feces,
        "Afeces_daptot": label_daptot_feces,
        "KI__D3GEX": f"{label_d3g} excretion urine",
        "KI__DAPEX": f"{label_dap} excretion\nurine",
        "KI__METEX": f"{label_dapmet} excretion\nurine",
        "KI__GLCEX": "glucose excretion\nurine",
    }

    # units
    unit_time = "hr"
    unit_metabolite = "µM"
    unit_metabolite_urine = "µmole"
    unit_metabolite_feces = "µmole"

    unit_dap = unit_metabolite
    unit_d3g = unit_metabolite
    unit_glc = "mM"
    unit_dap_urine = unit_metabolite_urine
    unit_d3g_urine = unit_metabolite_urine
    unit_daptot_urine = unit_metabolite_urine
    unit_glc_urine = "mole"
    unit_uge = "g"
    unit_rtg = "mM"
    unit_dap_feces = unit_metabolite_feces
    unit_daptot_feces = unit_metabolite_feces

    units: ClassVar[Mapping[str, str]] = {
        "time": unit_time,
        "[Cve_dap]": unit_metabolite,
        "[Cve_d3g]": unit_metabolite,
        "[Cve_daptot]": unit_metabolite,
        "[Cve_dapmet]": unit_metabolite,
        "[KI__fpg]": unit_glc,
        "KI__pg": unit_glc,
        "Aurine_dap": unit_metabolite_urine,
        "Aurine_d3g": unit_metabolite_urine,
        "Aurine_dapmet": unit_metabolite_urine,
        "Aurine_daptot": unit_metabolite_urine,
        "KI__glc_urine": unit_glc_urine,
        "KI__UGE": unit_uge,
        "KI__RTG": unit_rtg,
        "Afeces_dap": unit_dap_feces,
        "Afeces_daptot": unit_daptot_feces,
        "KI__D3GEX": "µmole/min",
        "KI__DAPEX": "µmole/min",
        "KI__METEX": "µmole/min",
        "KI__GLCEX": "µmole/min",
    }

    # ----------- Fasting/food -----
    #: the fed state of the intestine (`GU__fed`): food slows the absorption
    #: (`GU__FED_f_absorption`), Cmax 0.69, tmax +1 hr, AUC unchanged
    #: [Kasichayanula2011b]; a dose without a reported state is fasted
    fasting_map: ClassVar[dict[str, float]] = {
        "not reported": 0.0,
        "fasted": 0.0,
        "fed": 1.0,
    }
    fasting_colors: ClassVar[dict[str, str]] = {
        "fasted": "black",
        "fed": "tab:red",
    }

    # ----------- Renal function --------------
    #: the classes of the GFR with the end stage renal disease; a measured
    #: creatinine clearance or GFR of X ml/min is the renal function X/101
    renal_map = RENAL_FUNCTION_GFR
    renal_colors = RENAL_COLORS_END_STAGE

    # ----------- Interactions --------------
    #: the activity of UGT1A9 (`f_ugt1a9`, the glucuronidation to D3G of liver
    #: and kidney) with a coadministered drug, set from the mechanism and the
    #: observed change of the exposure, never fitted; the arms are validation
    #: data [Kasichayanula2013a]:
    #: - rifampicin 600 mg QD (inducer of UGT1A9, CYP3A4, UGT2B7): the AUC of
    #:   dapagliflozin falls to 0.78, the D3G AUC to 0.86 (model with the start
    #:   values: AUC 0.78, D3G 1.07; the induction of the other pathways is not
    #:   described)
    #: - mefenamic acid 500 mg, then 250 mg q6h (inhibitor of UGT1A9): the AUC
    #:   of dapagliflozin rises to 1.51, the D3G AUC falls to 0.78, D3G Cmax to
    #:   0.58 (model: AUC 1.51, D3G 0.83, D3G Cmax 0.68)
    #: The factors are set by hand against the observed AUC ratios of
    #: dapagliflozin. Metformin, pioglitazone, glimepiride, sitagliptin,
    #: voglibose and lobeglitazone do not change the exposure of dapagliflozin
    #: (AUC ratios 0.9-1.1, Kasichayanula2011c, Imamura2013, Jang2020): no
    #: factors.
    interactions: ClassVar[dict[str, dict[str, float]]] = {
        "rifampicin": {"f_ugt1a9": 1.37},
        "mefenamic acid": {"f_ugt1a9": 0.55},
    }

    def interaction_changes(self, drug: str) -> dict[str, Quantity]:
        """Get the changes of the model of a coadministered drug.

        Args:
            drug: the drug, a key of `interactions`.

        Returns:
            The activities of the enzymes by parameter id.
        """
        return {
            pid: Q(value, "dimensionless")
            for pid, value in self.interactions[drug].items()
        }

    # ----------- Cirrhosis map --------------
    cirrhosis_colors: ClassVar[Mapping[str, str]] = {
        "Control": "black",
        "Mild cirrhosis": "#C4DFF2",  # CPT A
        "Moderate cirrhosis": "#2b8cbe",  # CPT B
        "Severe cirrhosis": "#00446E",  # CPT C
    }

    # ----------- Doses map --------------
    # from light to dark with the hue turning from yellow over red to violet:
    # the lightness (OKLab) falls by 0.053 from dose to dose, adjacent doses
    # differ by about 7 (OKLab distance x 100)
    dose_colors: ClassVar[dict[float, str]] = {
        0: "black",
        0.001: "#fae776",
        0.01: "#ffcb57",
        0.1: "#feaf4c",
        0.3: "#fd9243",
        1: "#f4774f",
        2.5: "#e55f60",
        5: "#d04b6f",
        10: "#b63c7c",
        20: "#993185",
        25: "#792989",
        50: "#592387",
        100: "#3a1d7e",
        250: "#1a176e",
        500: "#001a45",
    }

    # molecular weights
    Mr: ClassVar[MolecularWeights] = MolecularWeights(
        dap=Q(408.873, "g/mole"),
        d3g=Q(585, "g/mole"),
        glc=Q(180.16, "g/mole"),
    )

    # --- Pharmacokinetic parameters ---
    def calculate_dapagliflozin_pk(
        self, scans: Iterable[str] | None = None
    ) -> dict[str, pd.DataFrame]:
        """Calculate the pharmacokinetic parameters of the scans.

        Args:
            scans: keys of the simulations, all simulations by default.

        Returns:
            The parameters of every scan by simulation key, per draw of the
            fitted parameters in a run over them (`uncertainty`).
        """
        keys = self._simulations if scans is None else scans
        across = None if self.uncertainty is None else self.uncertainty.id
        return {
            key: calculate_dapagliflozin_pk(
                experiment=self, res=self.results[f"task_{key}"], across=across
            )
            for key in keys
        }
