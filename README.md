[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13987865.svg)](https://doi.org/10.5281/zenodo.13987865)
[![Release](https://img.shields.io/github/v/release/matthiaskoenig/dapagliflozin-model)](https://github.com/matthiaskoenig/dapagliflozin-model/releases)
[![Documentation](https://github.com/matthiaskoenig/dapagliflozin-model/actions/workflows/docs.yml/badge.svg)](https://matthiaskoenig.github.io/dapagliflozin-model/)
[![CI](https://github.com/matthiaskoenig/dapagliflozin-model/actions/workflows/ci.yml/badge.svg)](https://github.com/matthiaskoenig/dapagliflozin-model/actions/workflows/ci.yml)
[![Docker](https://img.shields.io/docker/v/matthiaskoenig/dapagliflozin?sort=semver&label=docker)](https://hub.docker.com/r/matthiaskoenig/dapagliflozin)
[![Python](https://img.shields.io/badge/python-3.14-blue)](https://www.python.org)
[![Code license: MIT](https://img.shields.io/badge/code%20license-MIT-blue)](https://opensource.org/license/MIT)
[![Model license: CC BY-SA 4.0](https://img.shields.io/badge/model%20license-CC%20BY--SA%204.0-lightgrey)](https://creativecommons.org/licenses/by-sa/4.0/)

# Physiologically based pharmacokinetic/pharmacodynamic (PBPK/PD) model of dapagliflozin

This repository provides the dapagliflozin physiologically based pharmacokinetic/pharmacodynamic (PBPK/PD) model: the absorption, the tissue distribution, the hepatic and renal metabolism to the 3-O-glucuronide and the minor metabolites and the renal excretion of dapagliflozin, and the inhibition of SGLT2, which lowers the renal threshold of glucose and increases the urinary glucose excretion.

Documentation: <https://matthiaskoenig.github.io/dapagliflozin-model/>

## Quick start

```bash
git clone https://github.com/matthiaskoenig/dapagliflozin-model.git && cd dapagliflozin-model
uv sync
uv run run_simulation pkdb_models.models.dapagliflozin --action all -r results
```

## How to cite

> Nemitz, N., Elias, M. & König, M. (2026). *Physiologically based pharmacokinetic/pharmacodynamic (PBPK/PD) model of dapagliflozin*. Zenodo. https://doi.org/10.5281/zenodo.13987865

## License

* Source Code: [MIT](https://opensource.org/license/MIT)
* Documentation: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
* Models: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)

&copy; 2024-2026 Nike Nemitz, Michelle Elias and Matthias König
