# Dapagliflozin

Physiologically based pharmacokinetic/pharmacodynamic (PBPK/PD) model of the SGLT2 inhibitor dapagliflozin: its oral absorption (dissolution and absorption, a slower absorption with food), its distribution into the tissues (a fast pool and a slow pool of the rest tissue), its glucuronidation by UGT1A9 in the liver and the kidney to dapagliflozin-3-O-glucuronide (D3G) and the minor metabolites, and the renal excretion of dapagliflozin and its metabolites relative to the glomerular filtration rate. The plasma dapagliflozin inhibits SGLT2 and lowers the renal threshold of glucose, which increases the urinary glucose excretion (UGE) depending on the plasma glucose. The model describes the effects of the dose, food, renal and hepatic impairment, type 2 diabetes and the UGT1A9 interactions of rifampicin and mefenamic acid.

The published version of the model is [dapagliflozin-model](https://github.com/matthiaskoenig/dapagliflozin-model) ([10.5281/zenodo.13987865](https://doi.org/10.5281/zenodo.13987865)).

## Structure

- `__init__.py`: the paths of the model and its model object `MODEL` with its release metadata and its adopted fits, which the commands below run (pkpdsim).
- `models/`: the tissue models (`model_intestine.py`, `model_liver.py`, `model_kidney.py`) and the whole-body model (`model_body.py`). `models/factory.py` creates the SBML models, the flattened whole-body model which is simulated, the ODEs and the COMBINE archive in `models/output/`; the fitted values come from the adopted fits (`models/fits.py`, `models/output/fits/`).
- `experiments/`: the simulation experiments of the clinical studies (`studies/`), the dose dependency, the food effect and the hepatic and renal impairment (`misc/`) and the parameter scans of the dose, food, cirrhosis and renal function at plasma glucose levels of 5-15 mM with the pharmacokinetic parameters and the UGE of the studies (`scans/`, data in `data/parameters/`). The fit mappings are described with the vocabulary of `pkpdsim.metadata`.
- `fitting/`: the fit problems of the pharmacokinetics (`PK`) and of the pharmacodynamics on the adopted pharmacokinetics (`PD`), their parameters and the selection of the training, validation, aggregate and outlier data.
- `sensitivity/`: the sensitivity analysis of the pharmacokinetics of dapagliflozin and D3G and of the UGE in 24 hr.
- `dapagliflozin_pk.py`: the pharmacokinetic parameters of the scans (the non-compartmental analysis of `pkpdutils`), per draw of the fitted parameters in a run over them.
- `tolerances.py`: the absolute tolerances of the integrator per state.
- `simulations.py`: the simulation experiments and their groups, `MODEL` runs them.
- `data/`: the data of the studies, a subset of the PK-DB studies (`data/collect_data.py`); the copies are corrected against the publications (`TODO.md`).
- `food/`, `latex/`, `literature/`: the food of the studies, the tables of the studies and models of the publication and the literature search.
- `TODO.md`: the open issues of the model.

## Commands

```bash
# models, simulations
run_simulation pkdb_models.models.dapagliflozin --action factory
run_simulation pkdb_models.models.dapagliflozin --action list_experiments
run_simulation pkdb_models.models.dapagliflozin --action simulate --experiments all
run_simulation pkdb_models.models.dapagliflozin --action simulate --experiments all --uncertainty

# fits (the PD fit on the adopted PK), adoption of a fit, reports, identifiability
run_fit pkdb_models.models.dapagliflozin --subset=PK --runs=3 --cores=3 --method=DE
run_fit pkdb_models.models.dapagliflozin --subset=PD --runs=3 --cores=3 --method=DE
run_fit pkdb_models.models.dapagliflozin --adopt PK_<fit id>
run_report pkdb_models.models.dapagliflozin --subset=PK results/fit/PK_<fit id>/parameters.json
run_identifiability pkdb_models.models.dapagliflozin --subset=PK results/fit/PK_<fit id>/parameters.json

# sensitivity analysis (local by default)
python -m pkdb_models.models.dapagliflozin.sensitivity.sensitivity_analysis
```

The results are written into `results/` of the model.

## Adopted fits and uncertainty

The values of the fitted parameters come from the adopted fits of the model,
not from numbers in the code: `PK_20261020_002954__89bab` (8 parameters of the
dissolution, the absorption, the food effect, the metabolism, the renal
excretion and the distribution, all identifiable) and
`PD_20261020_014612__f2826` (the half maximal concentration, the Hill
coefficient and the maximal reduction of the renal threshold of glucose and
its increase with the fasting plasma glucose; the Hill coefficient is not
bounded upwards). `run_fit pkdb_models.models.dapagliflozin --adopt <fit id>`
copies the record of a finished fit (its parameters, problem, profile
likelihood and Fisher information) into `models/output/fits/<fit id>/`, which
is tracked; the ids of the adopted fits are listed in `fits` of the model
object (`__init__.py`). The factory takes the fitted values and their
uncertainty from the records (`models/fits.py`) and writes the uncertainty into
the SBML (`distrib:uncertainty`). `run_simulation ... --uncertainty` runs the
experiments over draws of these distributions and draws the simulations as the
median with a band of quantiles into `results/simulation_uncertainty`; the
parameter scans run on a coarser grid of output times then. A refit is: fit,
adopt, replace the id in `fits`, rebuild the model. See the documentation of
the repository on fitting (`docs/fitting.md`).
