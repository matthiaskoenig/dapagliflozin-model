"""Tool to run the dapagliflozin model factory and simulation experiments."""

import sys
import os
import subprocess
from enum import Enum
import optparse
from pathlib import Path
from pkdb_models.models.dapagliflozin import DAPAGLIFLOZIN_PATH
from pkdb_models.models.dapagliflozin.simulations import run_simulation_experiments
from sbmlutils.console import console

FACTORY_SCRIPT_PATH = DAPAGLIFLOZIN_PATH / "models" / "factory.py"

class Action(str, Enum):
    # simulations
    SIMULATE = "simulate"
    LIST_EXPERIMENTS = "list_experiments"
    # factory
    FACTORY = "factory"
    # run all
    ALL = "all"

def _setup_custom_results_paths(results_dir: str):
    """Override default results paths with custom directory for custom figure output directory."""
    import pkdb_models.models.dapagliflozin as dapagliflozin
    custom_path = Path(results_dir).resolve()
    custom_path.mkdir(parents=True, exist_ok=True)
    # Override the module paths
    dapagliflozin.RESULTS_PATH = custom_path
    dapagliflozin.RESULTS_PATH_SIMULATION = custom_path / "simulation"
    console.print(f"Figure output directory set to: [cyan]{custom_path}[/cyan]")
    return custom_path

def _get_current_results_path():
    """Get the current results path (default or custom)."""
    from pkdb_models.models.dapagliflozin import RESULTS_PATH
    return RESULTS_PATH

def _run_factory():
    """Executes the model factory script."""
    console.rule("[bold cyan]Running Model Factory[/bold cyan]", style="cyan")
    subprocess.run(
        [sys.executable, str(FACTORY_SCRIPT_PATH)],
        cwd=FACTORY_SCRIPT_PATH.parent,
        check=True,
    )
    console.print("[bold green]Factory finished.[/bold green]")
    console.rule(style="green")

def main() -> None:
    parser = optparse.OptionParser()
    parser.add_option(
        "-a", "--action",
        dest="action",
        help=f"The main action to perform. Choices: {[a.value for a in Action]} (required)",
    )
    parser.add_option(
        "-r", "--results-dir",
        dest="results_dir",
        help="Optional: Custom directory to save all results/figures (default: ./results)",
    )
    parser.add_option(
        "-e", "--experiments",
        dest="experiments",
        help="Comma-separated list of simulation experiments and/or groups (for '--action simulate'). "
             "Use '--action list_experiments' to see all available options.",
    )

    console.rule("[bold cyan]DAPAGLIFLOZIN PBPK/PD MODEL[/bold cyan]", style="cyan")

    options, args = parser.parse_args()

    def _parser_message(text: str) -> None:
        console.print(f"[bold red]Error: {text}[/bold red]")
        parser.print_help()
        console.rule(style="red")
        sys.exit(1)

    if not options.action:
        _parser_message("Required argument '--action' is missing.")

    try:
        action = Action(options.action.lower())
    except ValueError:
        _parser_message(f"Invalid action '{options.action}'. Please choose from {[a.value for a in Action]}.")

    # Setup custom results directory if provided
    if options.results_dir:
        _setup_custom_results_paths(options.results_dir)
    else:
        import pkdb_models.models.dapagliflozin as dapagliflozin
        default_path = _get_current_results_path()
        if options.action and options.action.lower() == "factory":
            console.print(f"[cyan]Factory will use model base path: {dapagliflozin.MODEL_BASE_PATH}[/cyan]")
        else:
            console.print(f"[cyan]Using figure results directory: {default_path}[/cyan]")

    if action == Action.FACTORY:
        _run_factory()

    elif action == Action.LIST_EXPERIMENTS:
        run_simulation_experiments(list_only=True)

    elif action == Action.SIMULATE:

        if not options.experiments:
            _parser_message("For '--action simulate', the '--experiments' argument is required.")
        exp_list = [e.strip() for e in options.experiments.split(",")]
        results_path = _get_current_results_path()
        console.rule("[bold cyan]Running Simulations[/bold cyan]", style="cyan")
        run_simulation_experiments(specific_experiments=exp_list)
        console.print("[bold green]Simulations finished.[/bold green]")
        console.print(f"[bold green]Results saved to: {results_path / 'simulation'}[/bold green]")

    elif action == Action.ALL:
        console.rule("[bold cyan]Running: Factory and all simulations.[/bold cyan]", style="cyan")
        _run_factory()                              # Run factory
        run_simulation_experiments(selected="all")  # Run all simulations
        console.print("\n[bold green]All scripts completed successfully![/bold green]")

    console.rule(style="white")


if __name__ == "__main__":
    """
    This script is intended to be run from the command line as 'run_dapagliflozin'.

    -------------------------------------------------
    Usage Examples (to be run in your terminal):
    -------------------------------------------------

    1. Help:
       Shows all available options and actions.
       $ run_dapagliflozin --help

    2. Run the Model Factory:
       Generates all SBML model files.
       $ run_dapagliflozin --action factory

    3. Run Simulations:
       List available experiments:
       $ run_dapagliflozin --action list_experiments

       Run all experiments:
       $ run_dapagliflozin --action simulate --experiments all

    5. Run Everything:
       Runs factory and all simulations.
       $ run_dapagliflozin --action all
       
       With custom results directory for figures:
       $ run_dapagliflozin --action all --results-dir '/path/to/my/results'
    """
    main()