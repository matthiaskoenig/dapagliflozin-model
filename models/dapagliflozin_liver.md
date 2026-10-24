# Model for the hepatic metabolism of dapagliflozin.

Model `dapagliflozin_liver`, SBML Level 3 Version 2, read from dapagliflozin\_liver.xml, written by sbmlode 0.3.0.

Model for the hepatic metabolism of dapagliflozin.

Dapagliflozin is metabolized mainly by UGT1A9 to the inactive dapagliflozin-3-O-glucuronide (D3G, 61 % of an oral dose in the urine), and to a minor extent to the 2-O-glucuronide (UGT2B4/2B7) and oxidative metabolites (CYP3A4), about 13 % of the dose in the urine and 6 % in the feces \[Obermeier2010\]. Shares of the metabolism: UGT1A9 79.9 %, CYP3A4 9.8 %, UGT2B4/2B7 8.9 % \[Jo2021\].

The metabolism is a first-order clearance of the plasma of the liver: the pharmacokinetics are dose proportional from 0.1 to 500 mg \[Komoroski2009\], far below the Km of UGT1A9 (479 µM in kidney microsomes \[Kasichayanula2013\]). It scales with the volume of the liver tissue, so with the tissue loss in cirrhosis.

The fraction f\_d3g of the metabolism forms D3G, the rest the other metabolites (dapmet).

The exposure rises with the Child-Pugh class (AUC ratio 1.03, 1.36, 1.67 \[Kasichayanula2011\]).

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
| $\mathrm{Vli}$ | `Vli` | liver | $1.5$ | l | ✓ |

## Species

| Symbol | Id | Name | Compartment | Value | Unit | Properties |
| --- | --- | --- | --- | --- | --- | --- |
| $\mathrm{dap}_{\mathrm{ext}}$ | `dap_ext` | dapagliflozin (plasma) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |
| $\mathrm{d3g}_{\mathrm{ext}}$ | `d3g_ext` | dapagliflozin-3-O-glucuronide (plasma) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |
| $\mathrm{dapmet}_{\mathrm{ext}}$ | `dapmet_ext` | other dapagliflozin metabolites (plasma) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |

## Parameters

| Symbol | Id | Name | Value | Unit | Constant |
| --- | --- | --- | --- | --- | :---: |
| $\mathrm{DAPMET}_{\mathrm{k}}$ | `DAPMET_k` | rate of the hepatic metabolism of dapagliflozin | $0.0658712610363756$ | 1/min | ✓ |
| $f_{\mathrm{d3g}}$ | `f_d3g` | fraction of the metabolism to D3G | $0.77$ | \- | ✓ |
| $f_{\mathrm{ugt1a9}}$ | `f_ugt1a9` | scaling factor UGT1A9 activity | $1$ | \- | ✓ |

## Reactions

| Rate | Id | Name | Equation |
| --- | --- | --- | --- |
| $v_{\mathrm{DAP2D3G}}$ | `DAP2D3G` | dapagliflozin glucuronidation (DAP2D3G, UGT1A9) | $\mathrm{dap}_{\mathrm{ext}} \longrightarrow \mathrm{d3g}_{\mathrm{ext}}$ |
| $v_{\mathrm{DAP2MET}}$ | `DAP2MET` | dapagliflozin other metabolism (DAP2MET) | $\mathrm{dap}_{\mathrm{ext}} \longrightarrow \mathrm{dapmet}_{\mathrm{ext}}$ |

The rates of the reactions are:

$$
\begin{aligned}
v_{\mathrm{DAP2D3G}} &= f_{\mathrm{ugt1a9}} \cdot f_{\mathrm{d3g}} \cdot \mathrm{DAPMET}_{\mathrm{k}} \cdot \mathrm{Vli} \cdot \mathrm{dap}_{\mathrm{ext}} \\[1ex]
v_{\mathrm{DAP2MET}} &= \mathopen{}\left(1 - f_{\mathrm{d3g}}\right) \cdot \mathrm{DAPMET}_{\mathrm{k}} \cdot \mathrm{Vli} \cdot \mathrm{dap}_{\mathrm{ext}}
\end{aligned}
$$

## ODE system

The states change in time with the rates of the reactions:

$$
\begin{aligned}
\frac{\mathrm{d} \,\mathrm{dap}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{DAP2D3G}} - v_{\mathrm{DAP2MET}}}{\mathrm{Vext}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{d3g}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{v_{\mathrm{DAP2D3G}}}{\mathrm{Vext}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{dapmet}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{v_{\mathrm{DAP2MET}}}{\mathrm{Vext}}
\end{aligned}
$$
