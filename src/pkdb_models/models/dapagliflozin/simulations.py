"""Run all simulation experiments."""
import shutil
from typing import List
from sbmlutils.console import console
from pkdb_models.models.dapagliflozin.experiments.scans.scan_dose import DapagliflozinDoseScan
from pkdb_models.models.dapagliflozin.helpers import run_experiments
from pkdb_models.models.dapagliflozin.experiments.studies import *
from pkdb_models.models.dapagliflozin.experiments.misc import *
from pkdb_models.models.dapagliflozin.experiments.scans import *
import pkdb_models.models.dapagliflozin as dapagliflozin
from sbmlutils import log
from sbmlsim.plot import Figure

Figure.legend_fontsize = 10
# Figure.fig_dpi = 600

logger = log.get_logger(__name__)

def run_simulation_experiments(selected: str = None, specific_experiments: List[str] = None, list_only: bool = False) -> None:
    """Run simulation experiments."""

    experiments = {
        "studies": [
            Boulton2013,
            Cho2021,
            FDAMB102002,
            FDAMB102003,
            FDAMB102006,
            FDAMB102007,
            # Gould2013,
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
            DapagliflozinDoseScan,
        ]

    }
    experiments["all"] = experiments["studies"] + experiments["misc"] + experiments["scan"]

    # experiment name collector for run_dapagliflozin.py
    if list_only:
        console.rule("[bold cyan]Available Simulation Experiments[/bold cyan]", style="cyan")
        console.print("\n[bold]You can use these group names:[/bold]")
        console.print(f"  {', '.join([g for g in experiments.keys()])}")
        console.print("\n[bold]Or these individual experiment names:[/bold]")
        for group_name in ["studies", "misc", "scan"]:
            if group_name in experiments and experiments[group_name]:
                console.print(f"\n[yellow]{group_name}:[/yellow]")
                for exp in experiments[group_name]:
                    console.print(f"{exp.__name__}")
        console.print("\n[dim]Use '--experiments' with comma-separated names to run specific experiments.[/dim]")
        console.print('[dim]Example: run_dapagliflozin --action simulate --experiments "misc,LaCreta2016"[/dim]')
        console.print('[dim]Or use "all" to run all experiments: run_dapagliflozin --action simulate --experiments all[/dim]\n')
        return

    # Determine which experiments to run
    if specific_experiments:
        # dictionary of all available experiments
        all_available_exp = {}
        for category, exp_list in experiments.items():
            if category != "all":
                for exp in exp_list:
                    all_available_exp[exp.__name__] = exp

        # Validate and collect requested experiments
        experiment_classes = []
        not_found = []
        for exp_name in specific_experiments:
            # group name
            if exp_name in experiments:
                experiment_classes.extend(experiments[exp_name])
            elif exp_name in all_available_exp:
                # individual experiment
                experiment_classes.append(all_available_exp[exp_name])
            else:
                not_found.append(exp_name)

        # Report any experiments that weren't found
        if not_found:
            console.rule(style="red bold")
            console.print(f"[red]Warning: The following experiments were not found: {', '.join(not_found)}[/red]")
            console.rule(style="red bold")
        if not experiment_classes:
            console.rule(style="red bold")
            console.print("[red]Error: No valid experiments to run![/red]")
            console.rule(style="red bold")
            return

        # output directory for custom selection
        output_dir = dapagliflozin.RESULTS_PATH_SIMULATION / "custom_selection"

    else:
        if not selected:
            console.print("\n[red bold]Error: No experiments specified![/red bold]")
            console.print("[yellow]For example, use selected='all' or selected='studies' or specific_experiments=[...][/yellow]\n")
            return
        if selected not in experiments:
            console.rule(style="red bold")
            console.print(
                f"[red]Error: Unknown group '{selected}'. Valid groups: {', '.join(experiments.keys())}[/red]")
            console.rule(style="red bold")
            return

        experiment_classes = experiments[selected]
        output_dir = dapagliflozin.RESULTS_PATH_SIMULATION / selected

    run_experiments(
        experiment_classes=experiment_classes,
        output_dir=output_dir,
    )
    # collect figures
    figures_dir = output_dir / "_figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    for f in output_dir.glob("**/*.png"):
        try:
            shutil.copy2(f, figures_dir / f.name)
        except shutil.SameFileError as err:
            pass
    console.print(f"Figures copied to: file://{figures_dir}", style="info")


if __name__ == "__main__":
    """
    # Run experiments
    
    run_simulation_experiments(specific_experiments=["Jang2020", "Kim2023"])
    run_simulation_experiments(selected="dose_dependency")
    
    # selected = "all"
    # selected = "misc"
    selected = "studies"
    # selected = "pharmacodynamics"
    # selected = "food"
    # selected = "hepatic_impairment"
    # selected = "renal_impairment"
    # selected = "scan"
    # selected = "dose_dependency"
    """

    run_simulation_experiments(selected="dose_dependency")
    # run_simulation_experiments(specific_experiments=["hepatic_impairment"])