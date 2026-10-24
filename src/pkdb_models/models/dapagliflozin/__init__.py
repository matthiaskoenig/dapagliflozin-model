"""Dapagliflozin: PBPK/PD model of dapagliflozin and the urinary glucose excretion.

The paths of the model and its model object `MODEL`, which the commands
`run_simulation`, `run_fit`, `run_report` and `run_identifiability` of pkpdsim
run by its package `pkdb_models.models.dapagliflozin`.
"""

from pathlib import Path

from pkpdsim.model import FigureSettings, PKPDModel, SimulatorSettings
from pkpdsim.release.settings import Author, ReleaseSettings

from pkdb_models.models.dapagliflozin.tolerances import ABSOLUTE_TOLERANCE

DAPAGLIFLOZIN_PATH = Path(__file__).parent

MODEL_BASE_PATH = DAPAGLIFLOZIN_PATH / "models" / "output" / "models"
MODEL_PATH = MODEL_BASE_PATH / "dapagliflozin_body_flat.xml"

# the subset of the PK-DB studies the model uses, see `data/collect_data.py`
DATA_PATH_BASE = DAPAGLIFLOZIN_PATH / "data"
DATA_PATH_DAPAGLIFLOZIN = DATA_PATH_BASE / "dapagliflozin"
DATA_PATH_CANAGLIFLOZIN = DATA_PATH_BASE / "canagliflozin"
DATA_PATHS = [
    DATA_PATH_DAPAGLIFLOZIN,
    DATA_PATH_CANAGLIFLOZIN,
]

MODEL = PKPDModel(
    name="dapagliflozin",
    package=__name__,
    title="Dapagliflozin PBPK/PD model",
    path=DAPAGLIFLOZIN_PATH,
    model_path=MODEL_PATH,
    data_paths=tuple(DATA_PATHS),
    archive_path=MODEL_BASE_PATH.parent / "dapagliflozin_model.omex",
    factory=".models.factory:create_dapagliflozin_models",
    experiments=".simulations:EXPERIMENTS",
    fit_definitions=".fitting.fitting:FIT_DEFINITIONS",
    fits=("PK_20261020_002954__89bab", "PD_20261020_014612__f2826"),
    fit_identifiability=True,
    simulator=SimulatorSettings(
        absolute_tolerance=ABSOLUTE_TOLERANCE, relative_tolerance=1e-10
    ),
    figures=FigureSettings(
        formats=("svg",), fig_dpi=300, fig_titlesize=20, legend_fontsize=10
    ),
    # the release of the model repository
    release=ReleaseSettings(
        title=(
            "Physiologically based pharmacokinetic/pharmacodynamic (PBPK/PD) model"
            " of dapagliflozin"
        ),
        repository="matthiaskoenig/dapagliflozin-model",
        doi="10.5281/zenodo.13987865",
        authors=(
            Author(
                "Nemitz, Nike",
                "0009-0006-4892-3279",
                "Humboldt-University Berlin, Institute for Theoretical Biology, Berlin",
            ),
            Author(
                "Elias, Michelle",
                "0009-0005-7959-3164",
                "Humboldt-University Berlin, Institute for Theoretical Biology, Berlin",
            ),
            Author(
                "König, Matthias",
                "0000-0003-1725-179X",
                (
                    "Humboldt-Universität zu Berlin, Faculty of Life Sciences, Department of Biology, Institute of Theoretical Biology, Systems Medicine of the Liver; University of Stuttgart, Institute of Structural Mechanics and Dynamics in Aerospace Engineering; University of Southern Denmark, Clinical Pharmacology, Pharmacy and Environmental Medicine, Department of Public Health; University Hospital Schleswig-Holstein, Campus Lübeck, First Department of Medicine, Lübeck, Germany"
                ),
            ),
        ),
        keywords=(
            "model",
            "dapagliflozin",
            "SGLT2 inhibitor",
            "urinary glucose excretion",
            "type 2 diabetes",
            "pharmacokinetics",
            "pharmacodynamics",
            "PBPK/PD",
            "physiologically based pharmacokinetic/ pharmacodynamic model",
            "SBML",
        ),
        docker_image="matthiaskoenig/dapagliflozin",
        funding=(
            "Matthias König was supported by the Federal Ministry of Research, Technology and Space (BMFTR, Germany) within ATLAS by grant number 031L0304B and by the German Research Foundation (DFG) within the Research Unit Program FOR 5151 QuaLiPerF (Quantifying Liver Perfusion-Function Relationship in Complex Resection - A Systems Medicine Approach) by grant number 436883643 and by grant number 465194077 (Priority Programme SPP 2311, Subproject SimLivA). This work was supported by the BMBF-funded de.NBI Cloud within the German Network for Bioinformatics Infrastructure (de.NBI) (031A537B, 031A533A, 031A538A, 031A533B, 031A535A, 031A537C, 031A534A, 031A532B)."
        ),
        description=(
            "This repository provides the dapagliflozin physiologically based"
            " pharmacokinetic/pharmacodynamic (PBPK/PD) model: the absorption,"
            " the tissue distribution, the hepatic and renal metabolism to the"
            " 3-O-glucuronide and the minor metabolites and the renal excretion"
            " of dapagliflozin, and the inhibition of SGLT2, which lowers the"
            " renal threshold of glucose and increases the urinary glucose"
            " excretion."
        ),
        first_year=2024,
    ),
)
"""The dapagliflozin model: its commands, see `pkpdsim.model`."""
