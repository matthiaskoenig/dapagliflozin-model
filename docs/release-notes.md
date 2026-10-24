# Release notes

## Release notes for dapagliflozin-model 0.10.0

Made from pkdb_models 3a0b49144b with pkpdsim 0.11.0.

The changes since 0.9.9.

### Core libraries

pkpdsim 0.11.0: not used by the previous release.

sbmlsim 0.3.0 → 0.13.0 (latest 5 of 428):

- **The experiments of a run run in a pool of worker processes** ([#387](https://github.com/matthiaskoenig/sbmlsim/pull/387)).
- **A figure of the serializer is saved without a second drawing** ([#388](https://github.com/matthiaskoenig/sbmlsim/pull/388)).
- **Formulas are parsed without importing the package petab** ([#389](https://github.com/matthiaskoenig/sbmlsim/pull/389)).
- **The points of a PK observable are analysed in padded groups** ([#390](https://github.com/matthiaskoenig/sbmlsim/pull/390)).
- **The derived units of a model are matched without `areIdentical` of libsbml** ([#383](https://github.com/matthiaskoenig/sbmlsim/pull/383)).

sbmlutils 0.9.4 → 0.17.0 (latest 5 of 269):

- `certifi>=2026.10.19` is a dependency: sbmlutils reaches the network through requests, which verifies with the CA bundle of certifi, and certifi 2026.10.19 removes the distrusted emSign root (GHSA-6q48-6c7w-9gw2) (#510)
- continuous integration, the documentation and the release are built with uv 0.13.4, whose bundled CA store removes the emSign root (GHSA-hq85-5cpj-mr62) (#510)
- model creation returns immutable, grouped construction diagnostics separately from the final validation; the opt-in `strict=True` rejects reported input losses, failed checked libSBML operations and undeclared targets of rules and initial assignments before writing, preserving an existing output file; the default continues to create missing target parameters; `ConstructionError` exposes the findings, annotation resources retained verbatim are informational (#509)
- the SBML parser reports grouped known round-trip losses through `Model.preservation_diagnostics`; factory results carry them independently from the construction and validation findings, also after merging models; the opt-in `strict_preservation=True` of parsing or model creation raises `PreservationError`, and a strict write rejects known source losses before touching its destination; default parsing and writing remain permissive (#509)
- `FactoryResult.validation` exposes the validation result of the final file, after annotations; `create_model(raise_on_error=True)` raises on validation errors while leaving the file available for diagnosis (#509)

### Libraries

| Package | Previous | New |
| --- | --- | --- |
| `libroadrunner` | 2.8.0 | 2.11.0 |
| `matplotlib` | 3.10.8 | 3.11.2 |
| `numpy` | 2.4.0 | 2.5.3 |
| `pandas` | 2.3.3 | 3.0.6 |
| `petab` | 0.8.0 | 0.9.0 |
| `pillow` | 12.0.0 | 12.3.0 |
| `pint` | 0.25.2 | 0.26.1 |
| `pkpdsim` | - | 0.11.0 |
| `pkpdutils` | - | 1.3.0 |
| `pymetadata` | 0.5.9 | 0.8.1 |
| `python-libsbml` | 5.20.5 | 5.21.2 |
| `rich` | 14.2.0 | 15.0.0 |
| `sbmlode` | - | 0.3.0 |
| `sbmlsim` | 0.3.0 | 0.13.0 |
| `sbmlutils` | 0.9.4 | 0.17.0 |
| `scikit-learn` | 1.8.0 | 1.9.1 |
| `scipy` | 1.16.3 | 1.18.1 |
| `typst` | - | 0.15.0 |
| `xarray` | 2025.12.0 | 2026.9.0 |

107 other packages changed.

<details>
<summary>All changed packages</summary>

| Package | Previous | New |
| --- | --- | --- |
| `altair` | 6.0.0 | - |
| `annotated-doc` | 0.0.4 | - |
| `annotated-types` | 0.7.0 | 0.8.0 |
| `antimony` | 3.1.1 | 3.1.3 |
| `anyio` | 4.12.0 | - |
| `asttokens` | 3.0.1 | - |
| `attrs` | 25.4.0 | 26.1.0 |
| `blosc2` | 3.12.2 | - |
| `bottleneck` | 1.6.0 | - |
| `cachetools` | 6.2.4 | - |
| `certifi` | 2025.11.12 | 2026.10.19 |
| `chardet` | 5.2.0 | 7.6.0 |
| `charset-normalizer` | 3.4.4 | 3.5.2 |
| `click` | 8.3.1 | 8.5.0 |
| `cloudpickle` | - | 3.1.2 |
| `coloredlogs` | 15.0.1 | - |
| `contourpy` | 1.3.3 | 1.4.0 |
| `decorator` | 5.2.1 | 5.3.1 |
| `deepmerge` | - | 3.0.1 |
| `depinfo` | 2.2.0 | - |
| `dill` | 0.4.0 | 0.4.1 |
| `elementpath` | 5.0.4 | - |
| `executing` | 2.2.1 | - |
| `fastapi` | 0.127.0 | - |
| `fastobo` | 0.14.1 | - |
| `fonttools` | 4.61.1 | 4.66.1 |
| `google-auth` | 2.45.0 | - |
| `google-auth-oauthlib` | 1.2.2 | - |
| `gspread` | 5.12.4 | - |
| `gspread-pandas` | 3.3.0 | - |
| `h11` | 0.16.0 | - |
| `h5netcdf` | - | 1.8.1 |
| `h5py` | - | 3.16.0 |
| `humanfriendly` | 10.0 | - |
| `idna` | 3.11 | 3.20 |
| `ipython` | 9.8.0 | - |
| `ipython-pygments-lexers` | 1.1.1 | - |
| `jedi` | 0.19.2 | - |
| `joblib` | 1.5.3 | 1.6.0 |
| `jsonschema` | 4.25.1 | 4.26.0 |
| `kiwisolver` | 1.4.9 | 1.5.1 |
| `latex2mathml` | - | 3.81.1 |
| `libroadrunner` | 2.8.0 | 2.11.0 |
| `lxml` | 6.0.2 | 6.1.3 |
| `markdown` | - | 3.11 |
| `markdown-it-py` | 4.0.0 | 4.2.0 |
| `markupsafe` | 3.0.3 | 3.0.4 |
| `matplotlib` | 3.10.8 | 3.11.2 |
| `matplotlib-inline` | 0.2.1 | - |
| `msgpack` | 1.1.2 | - |
| `multiprocess` | 0.70.18 | 0.70.19 |
| `narwhals` | 2.14.0 | 2.26.0 |
| `ndindex` | 1.10.1 | - |
| `networkx` | 3.6.1 | 3.7 |
| `numexpr` | 2.14.1 | - |
| `numpy` | 2.4.0 | 2.5.3 |
| `oauthlib` | 3.3.1 | - |
| `packaging` | 25.0 | 26.3 |
| `pandas` | 2.3.3 | 3.0.6 |
| `parso` | 0.8.5 | - |
| `pathspec` | - | 1.1.1 |
| `patsy` | 1.0.2 | - |
| `petab` | 0.8.0 | 0.9.0 |
| `pexpect` | 4.9.0 | - |
| `pillow` | 12.0.0 | 12.3.0 |
| `pint` | 0.25.2 | 0.26.1 |
| `pkdb-analysis` | 0.3.1 | - |
| `pkpdsim` | - | 0.11.0 |
| `pkpdutils` | - | 1.3.0 |
| `platformdirs` | 4.5.1 | 4.12.4 |
| `plotly` | 6.5.0 | - |
| `prompt-toolkit` | 3.0.52 | - |
| `pronto` | 2.7.2 | - |
| `psutil` | 7.1.3 | - |
| `ptyprocess` | 0.7.0 | - |
| `pure-eval` | 0.2.3 | - |
| `py-cpuinfo` | 9.0.0 | - |
| `py4cytoscape` | 1.12.0 | 1.13.0 |
| `pyasn1` | 0.6.1 | - |
| `pyasn1-modules` | 0.4.2 | - |
| `pydantic` | 2.12.5 | 2.14.1 |
| `pydantic-core` | 2.41.5 | 2.50.1 |
| `pydoe` | 0.3.8 | - |
| `pygments` | 2.19.2 | 2.21.0 |
| `pymdown-extensions` | - | 12.1 |
| `pymetadata` | 0.5.9 | 0.8.1 |
| `pyparsing` | 3.2.5 | 3.3.3 |
| `pyreadline3` | 3.5.4 | - |
| `python-libnuml` | 1.1.7 | - |
| `python-libsbml` | 5.20.5 | 5.21.2 |
| `python-libsedml` | 2.0.33 | - |
| `python-multipart` | 0.0.21 | - |
| `pytz` | 2025.2 | - |
| `requests` | 2.32.5 | 2.34.2 |
| `requests-oauthlib` | 2.0.0 | - |
| `rich` | 14.2.0 | 15.0.0 |
| `rpds-py` | 0.30.0 | 2026.9.1 |
| `rsa` | 4.9.1 | - |
| `salib` | 1.5.2 | 1.6.0 |
| `sbmlmath` | - | 0.4.1 |
| `sbmlode` | - | 0.3.0 |
| `sbmlsim` | 0.3.0 | 0.13.0 |
| `sbmlutils` | 0.9.4 | 0.17.0 |
| `scikit-learn` | 1.8.0 | 1.9.1 |
| `scipy` | 1.16.3 | 1.18.1 |
| `setproctitle` | 1.3.7 | - |
| `setuptools` | 80.9.0 | - |
| `stack-data` | 0.6.3 | - |
| `starlette` | 0.50.0 | - |
| `statsmodels` | 0.14.6 | - |
| `tables` | 3.10.2 | - |
| `threadpoolctl` | 3.6.0 | 3.7.0 |
| `tomli` | - | 2.5.0 |
| `traitlets` | 5.14.3 | - |
| `typing-extensions` | 4.15.0 | 4.16.0 |
| `typing-inspection` | 0.4.2 | 0.4.4 |
| `typst` | - | 0.15.0 |
| `tzdata` | 2025.3 | 2026.5 |
| `urllib3` | 2.6.2 | 2.8.0 |
| `uvicorn` | 0.40.0 | - |
| `wcwidth` | 0.2.14 | - |
| `xarray` | 2025.12.0 | 2026.9.0 |
| `xlsxwriter` | 3.2.9 | - |
| `xmlschema` | 4.2.0 | - |
| `xmltodict` | 1.0.2 | 1.0.4 |
| `zensical` | - | 0.0.70 |

</details>

### Data

- Changed: `canagliflozin/Sha2015`, `dapagliflozin/Boulton2013`, `dapagliflozin/FDAMB102002`, `dapagliflozin/FDAMB102003`, `dapagliflozin/FDAMB102006`, `dapagliflozin/FDAMB102007`, `dapagliflozin/Gould2013`, `dapagliflozin/Hwang2022a`, `dapagliflozin/Imamura2013`, `dapagliflozin/Jang2020`, `dapagliflozin/Kasichayanula2011`, `dapagliflozin/Kasichayanula2011a`, `dapagliflozin/Kasichayanula2011b`, `dapagliflozin/Kasichayanula2011c`, `dapagliflozin/Kasichayanula2012`, `dapagliflozin/Kasichayanula2013`, `dapagliflozin/Kasichayanula2013a`, `dapagliflozin/Khomitskaya2018`, `dapagliflozin/Kim2023`, `dapagliflozin/Kim2023a`, `dapagliflozin/Komoroski2009`, `dapagliflozin/LaCreta2016`, `dapagliflozin/Obermeier2010`, `dapagliflozin/Shah2019a`, `dapagliflozin/Watada2019`, `dapagliflozin/Yang2013`, `dapagliflozin/vanderAartvanderBeek2020`

### Analysis

#### Models

- `dapagliflozin_body.xml`: added 11 species, 15 reactions, 9 parameters; removed 2 parameters; changed values `Kp_dap` 25.5174 → 1.20018, `ftissue_dap` 0.0100001 → 0.411373
- `dapagliflozin_body_flat.xml`: added 12 species, 18 reactions, 17 parameters; removed 5 species, 4 reactions, 15 parameters; changed values `GU__DAPABS_k` 0.0594648 → 0.0444919, `GU__Ka_dis_dap` 0.84842 → 1.09176, `KI__Mr_glc` 180 → 180.16, `KI__RTG_E50` 9.38382e-06 → 5.7563e-06, `KI__RTG_base` 7.00191 → 11, `KI__RTG_gamma` 1.03599 → 2.80108, `KI__RTG_m_fpg` 1.65321 → 1.99099, `KI__RTG_max_inhibition` 0.712547 → 0.640014, `Kp_dap` 25.5174 → 1.20018, `ftissue_dap` 0.0100001 → 0.411373
- `dapagliflozin_intestine.xml`: added 2 parameters; changed values `DAPABS_k` 0.0594648 → 0.0444919, `Ka_dis_dap` 0.84842 → 1.09176
- `dapagliflozin_kidney.xml`: added 3 species, 2 reactions, 9 parameters; removed 3 species, 2 reactions, 11 parameters; changed values `Mr_glc` 180 → 180.16, `RTG_E50` 9.38382e-06 → 5.7563e-06, `RTG_base` 7.00191 → 11, `RTG_gamma` 1.03599 → 2.80108, `RTG_m_fpg` 1.65321 → 1.99099, `RTG_max_inhibition` 0.712547 → 0.640014
- `dapagliflozin_liver.xml`: added 1 species, 1 reaction, 2 parameters; removed 2 species, 2 reactions, 6 parameters

#### Experiments

No changes.

#### Fits

- `PD`: `PD_20261020_014612__f2826`, no fit record in the previous release
- `PK`: `PK_20261020_002954__89bab`, no fit record in the previous release

## Release notes for dapagliflozin-model v0.9.9

We are pleased to release a new version of the dapagliflozin-model including the 
following changes:

- model publication and reproducibility publication added
- FAIR assessment added
- uv update
- update to python 14
- updated github actions
- bugfix uv environment fixed

Your dapagliflozin-model team

## Release notes for dapagliflozin-model v0.9.8

We are pleased to release a new version of the dapagliflozin-model including the 
following changes:

- sensitivity analysis
- injected all model parameters
- removed unused parameters
- improved names for parameters

Your dapagliflozin-model team

## Release notes for dapagliflozin-model v0.9.7

We are pleased to release a new version of the dapagliflozin-model including the 
following changes:

- updated model annotations

Your dapagliflozin-model team

## Release notes for dapagliflozin-model v0.9.6

We are pleased to release a new version of the dapagliflozin-model including the 
following changes:

- python package
- docker container

Your dapagliflozin-model team

## Release notes for dapagliflozin-model v0.9.5

We are pleased to release a new version of the dapagliflozin-model including the 
following changes:

- updated pharmacodynamics model (RTG dependency on dapagliflozin)
- parameter optimization for dapagliflozin

Your dapagliflozin-model team

## Release notes for dapagliflozin-model v0.9.4

We are pleased to release a new version of the dapagliflozin-model including the 
following changes:

- updated pharmacodynamics model
- parameter optimization for dapagliflozin

Your dapagliflozin-model team

## Release notes for dapagliflozin-model v0.9.3

We are pleased to release a new version of the dapagliflozin-model including the 
following changes:

- added the omex file
- updated metadata
- licensing under MIT

Your dapagliflozin-model team

## Release notes for dapagliflozin-model v0.9.2

We are pleased to release a new version of the dapagliflozin-model including the 
following changes:

- added the omex file
- parameter optimization

Your dapagliflozin-model team

## Release notes for dapagliflozin-model v0.9.1

We are pleased to release a new version of the dapagliflozin-model including the 
following changes:

- initial model release
- zenodo integration

Your dapagliflozin-model team

## Release notes for dapagliflozin-model v0.9.0

We are pleased to release a new version of the dapagliflozin-model including the 
following changes:

- initial model release

Your dapagliflozin-model team

