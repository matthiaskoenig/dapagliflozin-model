"""
Reusable functionality for multiple simulation experiments.
"""
from collections import namedtuple
from typing import Dict
import pandas as pd

from pkdb_models.models.dapagliflozin.pk import calculate_dapagliflozin_pk
from pkdb_models.models.dapagliflozin import MODEL_PATH
from sbmlsim.experiment import SimulationExperiment
from sbmlsim.model import AbstractModel
from sbmlsim.task import Task


# Constants for conversion
MolecularWeights = namedtuple("MolecularWeights", "dap d3g daptot glc")


class DapagliflozinSimulationExperiment(SimulationExperiment):
    """Base class for all SimulationExperiments."""

    font = {"weight": "bold", "size": 22}
    scan_font = {"weight": "bold", "size": 20}
    tick_font_size = 15
    legend_font_size = 13
    suptitle_font_size = 25

    # labels
    label_time = "time"

    label_dap = "dapagliflozin"
    label_d3g = "dapagliflozin-3-O-\nglucuronide"
    label_daptot = "DAP + D3G"

    label_dap_plasma = f"{label_dap}\nplasma"
    label_d3g_plasma = f"{label_d3g} plasma"
    label_daptot_plasma = f"{label_daptot}\nplasma"

    label_glc_plasma = "glucose plasma"
    label_uge = "UGE"
    label_rtg = "RTG"
    label_dap_urine = label_dap + "\nurine"
    label_d3g_urine = label_d3g + " urine"
    label_daptot_urine = label_daptot + "\nurine"

    label_glc_urine = "glucose urine"
    label_dap_feces = label_dap + "\nfeces"
    label_daptot_feces = label_daptot + "\nfeces"

    labels: Dict[str, str] = {
        "time": "time",
        "[Cve_dap]": label_dap_plasma,
        "[Cve_d3g]": label_d3g_plasma,
        "[Cve_daptot]": label_daptot_plasma,
        "[KI__glc_ext]": label_glc_plasma,
        "Aurine_dap": label_dap_urine,
        "Aurine_d3g": label_d3g_urine,
        "Aurine_daptot": label_daptot_urine,
        "KI__glc_urine": label_glc_urine,
        "KI__UGE": label_uge,
        "KI__RTG": label_rtg,
        "Afeces_dap": label_dap_feces,
        "Afeces_daptot": label_daptot_feces,
        "KI__D3GEX": f"{label_d3g}excretion urine",
        "KI__DAPEX": f"{label_dap} excretion\nurine",
        "KI__GLCEX": f"glucose excretion\nurine",
    }

    # units
    unit_time = "hr"
    unit_metabolite = "µM"
    unit_metabolite_urine = "µmole"
    unit_metabolite_feces = "µmole"

    unit_dap = unit_metabolite
    unit_d3g = unit_metabolite
    unit_daptot = unit_metabolite
    unit_glc = "mM"
    unit_dap_urine = unit_metabolite_urine
    unit_d3g_urine = unit_metabolite_urine
    unit_daptot_urine = unit_metabolite_urine
    unit_glc_urine = "mole"
    unit_uge = "g"
    unit_rtg = "mM"
    unit_dap_feces = unit_metabolite_feces
    unit_daptot_feces = unit_metabolite_feces

    units: Dict[str, str] = {
        "time": unit_time,
        "[Cve_dap]": unit_metabolite,
        "[Cve_d3g]": unit_metabolite,
        "[Cve_daptot]": unit_metabolite,
        "[KI__glc_ext]": unit_glc,
        "Aurine_dap": unit_metabolite_urine,
        "Aurine_d3g": unit_metabolite_urine,
        "Aurine_daptot": unit_metabolite_urine,
        "KI__glc_urine": unit_glc_urine,
        "KI__UGE": unit_uge,
        "KI__RTG": unit_rtg,
        "Afeces_dap": unit_dap_feces,
        "Afeces_daptot": unit_daptot_feces,
        "KI__D3GEX": "µmole/min",
        "KI__DAPEX": "µmole/min",
        "KI__GLCEX": "µmole/min",
    }

    # ----------- Body weight -----
    # fall backs if no information is available
    bodyweight_default = 75  # [kg]
    bodyweight_asian = 65  # [kg]

    # ----------- Fasting plasma glucose -----
    # fallbacks if no information is available
    fpg_healthy = 5.0  # [mM]
    fpg_t2dm = 7.5  # [mM]
    fpg_t1dm = 7.5  # [mM]

    # ----------- Fasting/food -----
    fasting_map = {
        "not reported": 1.0,  # assuming fasted state if nothing is reported
        "fasted": 1.0,
        "fed": 0.3,
    }
    fasting_colors = {
        "fasted": "black",
        "fed": "tab:red",
    }

    # ----------- Renal map --------------
    renal_map = {
        "Normal renal function": 101.0 / 101.0,  # 1.0,
        "Mild renal impairment": 69.5 / 101.0,  # 0.69
        "Moderate renal impairment": 32.5 / 101.0,  # 0.32
        "Severe renal impairment": 19.5 / 101.0,  # 0.19
    }
    renal_colors = {
        "Normal renal function": "black",
        "Mild renal impairment": "#66c2a4",
        "Moderate renal impairment": "#2ca25f",
        "Severe renal impairment": "#006d2c",
    }

    # ----------- Cirrhosis map --------------
    cirrhosis_map = {
        "Control": 0,
        "Mild cirrhosis": 0.3994897959183674,  # CPT A
        "Moderate cirrhosis": 0.6979591836734694,  # CPT B
        "Severe cirrhosis": 0.8127551020408164,  # CPT C
    }
    cirrhosis_colors = {
        "Control": "black",
        "Mild cirrhosis": "#C4DFF2",  # CPT A
        "Moderate cirrhosis": "#2b8cbe",  # CPT B
        "Severe cirrhosis": "#00446E",  # CPT C
    }

    # ----------- Doses map --------------

    dose_colors = {
        0: "black",
        0.001: "#ffeb3b",
        0.01: "#ffd54f",
        0.1: "#ffb74d",
        0.3: "#ffa726",
        1: "#ff9800",
        2.5: "#fb8c00",
        5: "#f57c00",
        10: "#ef6c00",
        20: "#e65100",
        25: "#d84315",
        50: "#bf360c",
        100: "#a52a00",
        250: "#872200",
        500: "#5e1800",
    }

    def models(self) -> Dict[str, AbstractModel]:
        Q_ = self.Q_
        return {
            "model": AbstractModel(
                source=MODEL_PATH,
                language_type=AbstractModel.LanguageType.SBML,
                changes={},
            )
        }

    @staticmethod
    def _default_changes(Q_):
        """Default changes to simulations."""

        changes = {
        }

        return changes

    def default_changes(self: SimulationExperiment) -> Dict:
        """Default changes to simulations."""
        return DapagliflozinSimulationExperiment._default_changes(Q_=self.Q_)

    def tasks(self) -> Dict[str, Task]:
        if self.simulations():
            return {
                f"task_{key}": Task(model="model", simulation=key)
                for key in self.simulations()
            }
        return {}

    def data(self) -> Dict:
        self.add_selections_data(
            selections=[
                "time",
                # dosing
                "IVDOSE_dap",
                "PODOSE_dap",

                # venous
                "[Cve_daptot]",
                "[Cve_dap]",
                "[Cve_d3g]",
                "[KI__glc_ext]",

                # liver
                "[LI__dap]",
                "[LI__d3g]",
                
                # urine
                "Aurine_daptot",
                "Aurine_dap",
                "Aurine_d3g",
                "KI__glc_urine",
                "KI__UGE",
                "KI__RTG",

                # feces
                "Afeces_dap",
                "Afeces_daptot",

                # renal excretion rate
                "KI__D3GEX",
                "KI__DAPEX",
                "KI__GLCEX",

                # cases
                'KI__f_renal_function',
                'f_cirrhosis',
                'GU__f_absorption',
                'f_ugt1a9',
            ]
        )
        return {}

    @property
    def Mr(self):
        return MolecularWeights(
            dap=self.Q_(408.873, "g/mole"),
            d3g=self.Q_(585, "g/mole"),
            daptot=self.Q_((408.873 + 585)/2.0, "g/mole"),  # FIXME: this is only approximation
            glc=self.Q_(180, "g/mole"),
        )

    # --- Pharmacokinetic parameters ---
    pk_labels = {
        "auc": "AUCend",
        "aucinf": "AUC",
        "cl": "Total clearance",
        "cl_renal": "Renal clearance",
        "cl_hepatic": "Hepatic clearance",
        "cmax": "Cmax",
        "thalf": "Half-life",
        "kel": "kel",
        "vd": "vd",
        "Aurine_eat": "Enalaprilat urine",
    }

    pk_units = {
        "auc": "µmole/l*hr",
        "aucinf": "µmole/l*hr",
        "cl": "ml/min",
        "cl_renal": "ml/min",
        "cl_hepatic": "ml/min",
        "cmax": "µmole/l",
        "thalf": "hr",
        "kel": "1/hr",
        "vd": "l",
        "Aurine_eat": "µmole",
    }

    def calculate_dapagliflozin_pk(self, scans: list = []) -> Dict[str, pd.DataFrame]:
       """Calculate pk parameters for simulations (scans)"""
       pk_dfs = {}
       if scans:
           for sim_key in scans:
               xres = self.results[f"task_{sim_key}"]
               df = calculate_dapagliflozin_pk(experiment=self, xres=xres)
               pk_dfs[sim_key] = df
       else:
           for sim_key in self._simulations.keys():
               xres = self.results[f"task_{sim_key}"]
               df = calculate_dapagliflozin_pk(experiment=self, xres=xres)
               pk_dfs[sim_key] = df
       return pk_dfs

