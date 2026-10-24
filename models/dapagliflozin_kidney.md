# Model for the renal excretion of dapagliflozin and its metabolites and the urinary glucose excretion.

Model `dapagliflozin_kidney`, SBML Level 3 Version 2, read from dapagliflozin\_kidney.xml, written by sbmlode 0.3.0.

Model for the renal excretion of dapagliflozin and its metabolites, the renal metabolism of dapagliflozin and the urinary glucose excretion.

1\.2-1.6 % of an oral dose of dapagliflozin is excreted unchanged in the urine with a renal clearance of 3.5-5.6 ml/min, below the filtration of the unbound fraction (fu 0.09, 9 ml/min): net reabsorption \[Obermeier2010, Kasichayanula2013\].

61 % of the dose is excreted as dapagliflozin-3-O-glucuronide (D3G) in the urine with a renal clearance of 151 ml/min in healthy subjects, above the GFR (tubular secretion) \[Kasichayanula2013\]; 75.2 % of the radioactivity is in the urine \[Obermeier2010\].

The renal excretion is relative to the glomerular filtration rate (GFR).

The kidney glucuronidates dapagliflozin (UGT1A9, 8-fold mRNA and 3-fold D3G formation per mg microsomal protein of the liver \[Kasichayanula2013\]); the renal metabolism scales with the renal function (the functional renal mass), which makes the exposure of dapagliflozin rise in renal impairment (AUC ratio 1.60, 2.07, 1.95 mild, moderate, severe \[Kasichayanula2013\]).

Dapagliflozin lowers the renal threshold of glucose (RTG) by the inhibition of SGLT2; glucose above the RTG is excreted with the GFR.

The content of this model has been carefully created in a manual research effort.

Terms of use

Copyright © 2026 Matthias König.

This work is licensed under a Creative Commons Attribution-ShareAlike 4.0 International License.

Redistribution and use of any part of this model, with or without modification, are permitted provided that the following conditions are met:

Redistributions of this SBML file must retain the above copyright notice, this list of conditions and the following disclaimer.

Redistributions in a different form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

This model is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

## Units

| Time | Substance | Extent | Volume | Area | Length |
| --- | --- | --- | --- | --- | --- |
| min | mmol | mmol | l | m<sup>2</sup> | m |

## Compartments

| Symbol | Id | Name | Size | Unit | Constant |
| --- | --- | --- | --- | --- | :---: |
| $\mathrm{Vext}$ | `Vext` | plasma | $1.5$ | l | ✓ |
| $\mathrm{Vki}$ | `Vki` | kidney | $0.3$ | l | ✓ |
| $\mathrm{Vurine}$ | `Vurine` | urine | $1$ | l | ✓ |

## Species

| Symbol | Id | Name | Compartment | Value | Unit | Properties |
| --- | --- | --- | --- | --- | --- | --- |
| $\mathrm{dap}_{\mathrm{ext}}$ | `dap_ext` | dapagliflozin (plasma) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |
| $\mathrm{dap}_{\mathrm{urine}}$ | `dap_urine` | dapagliflozin (urine) | $\mathrm{Vurine}$ |  | mmol | amount |
| $\mathrm{d3g}_{\mathrm{ext}}$ | `d3g_ext` | dapagliflozin-3-O-glucuronide (plasma) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |
| $\mathrm{d3g}_{\mathrm{urine}}$ | `d3g_urine` | dapagliflozin-3-O-glucuronide (urine) | $\mathrm{Vurine}$ |  | mmol | amount |
| $\mathrm{dapmet}_{\mathrm{ext}}$ | `dapmet_ext` | other dapagliflozin metabolites (plasma) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |
| $\mathrm{dapmet}_{\mathrm{urine}}$ | `dapmet_urine` | other dapagliflozin metabolites (urine) | $\mathrm{Vurine}$ |  | mmol | amount |
| $\mathrm{fpg}$ | `fpg` | fasting plasma glucose (FPG) | $\mathrm{Vext}$ | $5$ | mmol/l | concentration, boundary |
| $\mathrm{glc}_{\mathrm{urine}}$ | `glc_urine` | glucose (urine) | $\mathrm{Vurine}$ |  | mmol | amount |

## Parameters

| Symbol | Id | Name | Value | Unit | Constant |
| --- | --- | --- | --- | --- | :---: |
| $f_{\mathrm{renal\_function}}$ | `f_renal_function` | renal function | $1$ | \- | ✓ |
| $\mathrm{GFR}_{\mathrm{base}}$ | `GFR_base` | glomerular filtration rate (healthy) | $101$ | ml/min | ✓ |
| $\mathrm{fu}_{\mathrm{dap}}$ | `fu_dap` | fraction unbound in plasma dapagliflozin | $0.09$ | \- | ✓ |
| $\mathrm{DAPEX}_{\mathrm{f}}$ | `DAPEX_f` | renal clearance of dapagliflozin relative to its filtration | $0.460206039705467$ | \- | ✓ |
| $\mathrm{D3GEX}_{\mathrm{f}}$ | `D3GEX_f` | renal clearance of the metabolites relative to the GFR | $1.80545550669025$ | \- | ✓ |
| $\mathrm{DAPMET}_{\mathrm{k}}$ | `DAPMET_k` | rate of the renal metabolism of dapagliflozin | $0.1$ | 1/min | ✓ |
| $f_{\mathrm{d3g}}$ | `f_d3g` | fraction of the metabolism to D3G | $0.77$ | \- | ✓ |
| $f_{\mathrm{ugt1a9}}$ | `f_ugt1a9` | scaling factor UGT1A9 activity | $1$ | \- | ✓ |
| $\mathrm{Mr}_{\mathrm{glc}}$ | `Mr_glc` | molecular weight glucose | $180.16$ | g/mol | ✓ |
| $\mathrm{cf}_{\mathrm{mg\_per\_g}}$ | `cf_mg_per_g` | conversion factor mg per g | $1000$ | mg/g | ✓ |
| $\mathrm{fpg}_{\mathrm{healthy}}$ | `fpg_healthy` | fasting plasma glucose (healthy) | $5$ | mmol/l | ✓ |
| $\mathrm{PG}_{\mathrm{delta}}$ | `PG_delta` | mean plasma glucose of the day minus the fasting plasma glucose | $0.85$ | mmol/l | ✓ |
| $\mathrm{RTG}_{\mathrm{E50}}$ | `RTG_E50` | EC50 reduction in RTG | $5.75629727963379 \times 10^{-6}$ | mmol/l | ✓ |
| $\mathrm{RTG}_{\mathrm{gamma}}$ | `RTG_gamma` | hill coefficient reduction in RTG | $2.80107895262435$ | \- | ✓ |
| $\mathrm{RTG}_{\mathrm{base}}$ | `RTG_base` | baseline RTG | $11$ | mmol/l | ✓ |
| $\mathrm{RTG}_{\mathrm{m\_fpg}}$ | `RTG_m_fpg` | FPG effect on RTG | $1.99098619575577$ | \- | ✓ |
| $\mathrm{RTG}_{\mathrm{max\_inhibition}}$ | `RTG_max_inhibition` | maximal relative reduction of the RTG | $0.640013975351734$ | \- | ✓ |
| $\mathrm{RTG}_{\mathrm{splay}}$ | `RTG_splay` | splay of the glucose excretion around the RTG | $1$ | mmol/l | ✓ |
| $\mathrm{GFR}$ | `GFR` | glomerular filtration rate |  | ml/min |  |
| $\mathrm{pg}$ | `pg` | mean plasma glucose of the day |  | mmol/l |  |
| $\mathrm{RTG}_{\mathrm{fpg}}$ | `RTG_fpg` | RTG without SGLT2 inhibition |  | mmol/l |  |
| $\mathrm{RTG}_{\mathrm{delta}}$ | `RTG_delta` | maximal reduction of the RTG |  | mmol/l |  |
| $\mathrm{RTG}$ | `RTG` | renal threshold glucose (RTG) |  | mmol/l |  |
| $\mathrm{UGE}$ | `UGE` | urinary glucose excretion (UGE) |  | gram |  |

## Initial assignments and assignment rules

The initial assignments set the values at $t = 0$:

$$
\begin{aligned}
\mathrm{dap}_{\mathrm{urine}} &= 0 \cdot \mathrm{Vurine} \\[1ex]
\mathrm{d3g}_{\mathrm{urine}} &= 0 \cdot \mathrm{Vurine} \\[1ex]
\mathrm{dapmet}_{\mathrm{urine}} &= 0 \cdot \mathrm{Vurine} \\[1ex]
\mathrm{glc}_{\mathrm{urine}} &= 0 \cdot \mathrm{Vurine}
\end{aligned}
$$

The assignment rules hold at every time $t$:

$$
\begin{aligned}
\mathrm{GFR} &= f_{\mathrm{renal\_function}} \cdot \mathrm{GFR}_{\mathrm{base}} \\[1ex]
\mathrm{pg} &= \mathrm{fpg} + \mathrm{PG}_{\mathrm{delta}} \\[1ex]
\mathrm{RTG}_{\mathrm{fpg}} &= \mathrm{RTG}_{\mathrm{base}} + \mathrm{RTG}_{\mathrm{m\_fpg}} \cdot \mathopen{}\left(\mathrm{fpg} - \mathrm{fpg}_{\mathrm{healthy}}\right) \\[1ex]
\mathrm{RTG}_{\mathrm{delta}} &= \mathrm{RTG}_{\mathrm{fpg}} \cdot \mathrm{RTG}_{\mathrm{max\_inhibition}} \\[1ex]
\mathrm{RTG} &= \mathrm{RTG}_{\mathrm{fpg}} - \frac{\mathrm{RTG}_{\mathrm{delta}} \cdot \mathrm{dap}_{\mathrm{ext}}^{\mathrm{RTG}_{\mathrm{gamma}}}}{\mathrm{RTG}_{\mathrm{E50}}^{\mathrm{RTG}_{\mathrm{gamma}}} + \mathrm{dap}_{\mathrm{ext}}^{\mathrm{RTG}_{\mathrm{gamma}}}} \\[1ex]
\mathrm{UGE} &= \frac{\mathrm{glc}_{\mathrm{urine}} \cdot \mathrm{Mr}_{\mathrm{glc}}}{\mathrm{cf}_{\mathrm{mg\_per\_g}}}
\end{aligned}
$$

## Reactions

| Rate | Id | Name | Equation | Modifiers |
| --- | --- | --- | --- | --- |
| $v_{\mathrm{DAP2D3G}}$ | `DAP2D3G` | dapagliflozin glucuronidation (DAP2D3G, UGT1A9) | $\mathrm{dap}_{\mathrm{ext}} \longrightarrow \mathrm{d3g}_{\mathrm{ext}}$ |  |
| $v_{\mathrm{DAP2MET}}$ | `DAP2MET` | dapagliflozin other metabolism (DAP2MET) | $\mathrm{dap}_{\mathrm{ext}} \longrightarrow \mathrm{dapmet}_{\mathrm{ext}}$ |  |
| $v_{\mathrm{DAPEX}}$ | `DAPEX` | dapagliflozin excretion (DAPEX) | $\mathrm{dap}_{\mathrm{ext}} \longrightarrow \mathrm{dap}_{\mathrm{urine}}$ |  |
| $v_{\mathrm{D3GEX}}$ | `D3GEX` | D3G excretion (D3GEX) | $\mathrm{d3g}_{\mathrm{ext}} \longrightarrow \mathrm{d3g}_{\mathrm{urine}}$ |  |
| $v_{\mathrm{METEX}}$ | `METEX` | other metabolites excretion (METEX) | $\mathrm{dapmet}_{\mathrm{ext}} \longrightarrow \mathrm{dapmet}_{\mathrm{urine}}$ |  |
| $v_{\mathrm{GLCEX}}$ | `GLCEX` | glucose excretion (GLCEX) | $\mathrm{fpg} \longrightarrow \mathrm{glc}_{\mathrm{urine}}$ | $\mathrm{dap}_{\mathrm{ext}}$ |

The rates of the reactions are:

$$
\begin{aligned}
v_{\mathrm{DAP2D3G}} &= f_{\mathrm{renal\_function}} \cdot f_{\mathrm{ugt1a9}} \cdot f_{\mathrm{d3g}} \cdot \mathrm{DAPMET}_{\mathrm{k}} \cdot \mathrm{Vki} \cdot \mathrm{dap}_{\mathrm{ext}} \\[1ex]
v_{\mathrm{DAP2MET}} &= f_{\mathrm{renal\_function}} \cdot \mathopen{}\left(1 - f_{\mathrm{d3g}}\right) \cdot \mathrm{DAPMET}_{\mathrm{k}} \cdot \mathrm{Vki} \cdot \mathrm{dap}_{\mathrm{ext}} \\[1ex]
v_{\mathrm{DAPEX}} &= \mathrm{GFR} \cdot 0.001 \cdot \mathrm{fu}_{\mathrm{dap}} \cdot \mathrm{DAPEX}_{\mathrm{f}} \cdot \mathrm{dap}_{\mathrm{ext}} \\[1ex]
v_{\mathrm{D3GEX}} &= \mathrm{GFR} \cdot 0.001 \cdot \mathrm{D3GEX}_{\mathrm{f}} \cdot \mathrm{d3g}_{\mathrm{ext}} \\[1ex]
v_{\mathrm{METEX}} &= \mathrm{GFR} \cdot 0.001 \cdot \mathrm{D3GEX}_{\mathrm{f}} \cdot \mathrm{dapmet}_{\mathrm{ext}} \\[1ex]
v_{\mathrm{GLCEX}} &= \mathrm{GFR} \cdot 0.001 \cdot \mathrm{RTG}_{\mathrm{splay}} \cdot \ln\mathopen{}\left(1 + \exp\mathopen{}\left(\frac{\mathrm{pg} - \mathrm{RTG}}{\mathrm{RTG}_{\mathrm{splay}}}\right)\right)
\end{aligned}
$$

## ODE system

The states change in time with the rates of the reactions:

$$
\begin{aligned}
\frac{\mathrm{d} \,\mathrm{dap}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{DAP2D3G}} - v_{\mathrm{DAP2MET}} - v_{\mathrm{DAPEX}}}{\mathrm{Vext}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{dap}_{\mathrm{urine}}}{\mathrm{d} t} &= v_{\mathrm{DAPEX}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{d3g}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{v_{\mathrm{DAP2D3G}} - v_{\mathrm{D3GEX}}}{\mathrm{Vext}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{d3g}_{\mathrm{urine}}}{\mathrm{d} t} &= v_{\mathrm{D3GEX}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{dapmet}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{v_{\mathrm{DAP2MET}} - v_{\mathrm{METEX}}}{\mathrm{Vext}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{dapmet}_{\mathrm{urine}}}{\mathrm{d} t} &= v_{\mathrm{METEX}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{glc}_{\mathrm{urine}}}{\mathrm{d} t} &= v_{\mathrm{GLCEX}}
\end{aligned}
$$
