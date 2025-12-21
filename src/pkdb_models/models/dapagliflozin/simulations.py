"""Run all simulation experiments."""
import shutil
from typing import List
from pathlib import Path
from sbmlutils.console import console
from pkdb_models.models.dapagliflozin.helpers import run_experiments
from pkdb_models.models.dapagliflozin.experiments.studies import *
from pkdb_models.models.dapagliflozin.experiments.misc import *
from pkdb_models.models.dapagliflozin.experiments.scans import *
import pkdb_models.models.dapagliflozin as dapagliflozin
from sbmlutils import log
from sbmlsim.plot import Figure

logger = log.get_logger(__name__)

EXPERIMENTS = {
    "studies": [
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
        Yang2013,
    ],
    "pharmacodynamics": [
        FDAMB102002,
        FDAMB102003,
        FDAMB102007,
        Gould2013,
        Kasichayanula2011a,
        Kim2023,
        Komoroski2009,
        Sha2015,
        Watada2019,
        Yang2013,
    ],
    "dose_dependency": [
        FDAMB102002,
        FDAMB102003,
        Gould2013,
        Kasichayanula2011a,
        Komoroski2009,
        Watada2019,
        Yang2013,
    ],
    "food": [
        Kasichayanula2011b,
        Komoroski2009,
        LaCreta2016,
        Shah2019a,
    ],
    "hepatic_impairment": [
        Kasichayanula2011,
    ],
    "renal_impairment": [
        Kasichayanula2013,
        FDAMB102007,
    ],
    "misc": [
        DoseDependencyExperiment,
        FoodEffect,
        HepaticRenalImpairment,
    ],
    "scan": [
        DapagliflozinParameterScan,
    ]
}

EXPERIMENTS["all"] = EXPERIMENTS["studies"] + EXPERIMENTS["misc"] + EXPERIMENTS["scan"]


def run_simulation_experiments(
    selected: str = None,
    experiment_classes: List = None,
    output_dir: Path = None
) -> None:
    """Run dapagliflozin simulation experiments."""

    Figure.fig_dpi = 600
    Figure.legend_fontsize = 10

    # Determine which experiments to run
    if experiment_classes is not None:
        experiments_to_run = experiment_classes
        if output_dir is None:
            output_dir = dapagliflozin.RESULTS_PATH_SIMULATION / "custom_selection"
    elif selected:
        # Using the 'selected' parameter
        if selected not in EXPERIMENTS:
            console.rule(style="red bold")
            console.print(
                f"[red]Error: Unknown group '{selected}'. Valid groups: {', '.join(EXPERIMENTS.keys())}[/red]"
            )
            console.rule(style="red bold")
            return
        experiments_to_run = EXPERIMENTS[selected]
        if output_dir is None:
            output_dir = dapagliflozin.RESULTS_PATH_SIMULATION / selected
    else:
        console.print("\n[red bold]Error: No experiments specified![/red bold]")
        console.print("[yellow]Use selected='all' or selected='studies' or provide experiment_classes=[...][/yellow]\n")
        return

    # Run the experiments
    run_experiments(experiment_classes=experiments_to_run, output_dir=output_dir)

    # Collect figures into one folder
    figures_dir = output_dir / "_figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    for f in output_dir.glob("**/*.png"):
        if f.parent == figures_dir:
            continue
        try:
            shutil.copy2(f, figures_dir / f.name)
        except Exception as err:
            print(f"file {f.name} in {f.parent} fails, skipping. Error: {err}")
    console.print(f"Figures copied to: file://{figures_dir}", style="info")


if __name__ == "__main__":
    """
    # Run experiments
    
    # selected = "all"
    # selected = "misc"
    # selected = "studies"
    # selected = "pharmacodynamics"
    # selected = "dose_dependency"
    # selected = "food"
    # selected = "hepatic_impairment"
    # selected = "renal_impairment"
    # selected = "scan"
    """

    run_simulation_experiments(selected="all")