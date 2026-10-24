# Model for dapagliflozin absorption in the small intestine

Model `dapagliflozin_intestine`, SBML Level 3 Version 2, read from dapagliflozin\_intestine.xml, written by sbmlode 0.3.0.

Model for dapagliflozin absorption

Dissolution of the oral dose and absorption of dapagliflozin; the fraction absorbed F\_dap\_abs of the dapagliflozin of the lumen is absorbed, the rest is excreted in the feces. After an oral solution of 50 mg of 14C-dapagliflozin 21 % of the radioactivity is recovered in the feces, 15.4 % of the dose as unchanged dapagliflozin \[Obermeier2010\]; the absolute bioavailability is 78 % \[Boulton2013\].

The enterohepatic circulation is minimal and not modeled.

Food (fed) slows the absorption (FED\_f\_absorption): a high fat meal lowers Cmax by 31 % and delays tmax by about 1 hr, the AUC is unchanged \[Kasichayanula2011b\].

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
| $\mathrm{Vext}$ | `Vext` | plasma | $1$ | l | ✓ |
| $\mathrm{Vgu}$ | `Vgu` | intestine | $1.2825$ | l | ✓ |
| $\mathrm{Vlumen}$ | `Vlumen` | intestinal lumen (inner part of intestine) | $1.15425$ | l |  |
| $\mathrm{Vfeces}$ | `Vfeces` | feces | $1$ | l | ✓ |
| $\mathrm{Vapical}$ | `Vapical` | apical membrane (intestinal membrane enterocytes) | $\mathrm{NaN}$ | m<sup>2</sup> | ✓ |
| $\mathrm{Vstomach}$ | `Vstomach` | stomach | $1$ | l | ✓ |

## Species

| Symbol | Id | Name | Compartment | Value | Unit | Properties |
| --- | --- | --- | --- | --- | --- | --- |
| $\mathrm{dap}_{\mathrm{stomach}}$ | `dap_stomach` | dapagliflozin (stomach) | $\mathrm{Vstomach}$ |  | mmol | amount, boundary |
| $\mathrm{dap}_{\mathrm{lumen}}$ | `dap_lumen` | dapagliflozin (intestinal lumen) | $\mathrm{Vlumen}$ | $0$ | mmol/l | concentration |
| $\mathrm{dap}_{\mathrm{ext}}$ | `dap_ext` | dapagliflozin (plasma) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |
| $\mathrm{dap}_{\mathrm{feces}}$ | `dap_feces` | dapagliflozin (feces) | $\mathrm{Vfeces}$ |  | mmol | amount |

## Parameters

| Symbol | Id | Name | Value | Unit | Constant |
| --- | --- | --- | --- | --- | :---: |
| $F_{\mathrm{dap\_abs}}$ | `F_dap_abs` | fraction absorbed dapagliflozin | $0.84$ | \- | ✓ |
| $\mathrm{DAPABS}_{\mathrm{k}}$ | `DAPABS_k` | rate of dapagliflozin absorption | $0.0444918821768688$ | 1/min | ✓ |
| $\mathrm{fed}$ | `fed` | fed state (0 fasted, 1 fed) | $0$ | \- | ✓ |
| $\mathrm{FED}_{\mathrm{f\_absorption}}$ | `FED_f_absorption` | fed effect on the absorption rate | $0.189289061510997$ | \- | ✓ |
| $f_{\mathrm{absorption}}$ | `f_absorption` | factor of the absorption rate |  | \- |  |
| $\mathrm{PODOSE}_{\mathrm{dap}}$ | `PODOSE_dap` | oral dose dapagliflozin \[mg\] | $0$ | mg |  |
| $\mathrm{Ka}_{\mathrm{dis\_dap}}$ | `Ka_dis_dap` | dissolution rate dapagliflozin | $1.09176322622891$ | 1/hr | ✓ |
| $\mathrm{Mr}_{\mathrm{dap}}$ | `Mr_dap` | Molecular weight dapagliflozin \[g/mole\] | $408.873$ | g/mol | ✓ |
| $\mathrm{absorption}$ | `absorption` | dapagliflozin leaving the lumen |  | mmol/min |  |

## Initial assignments and assignment rules

The initial assignments set the values at $t = 0$:

$$
\begin{aligned}
\mathrm{dap}_{\mathrm{stomach}} &= 0 \cdot \mathrm{Vstomach} \\[1ex]
\mathrm{dap}_{\mathrm{feces}} &= 0 \cdot \mathrm{Vfeces}
\end{aligned}
$$

The assignment rules hold at every time $t$:

$$
\begin{aligned}
f_{\mathrm{absorption}} &= 1 + \mathrm{fed} \cdot \mathopen{}\left(\mathrm{FED}_{\mathrm{f\_absorption}} - 1\right) \\[1ex]
\mathrm{absorption} &= f_{\mathrm{absorption}} \cdot \mathrm{DAPABS}_{\mathrm{k}} \cdot \mathrm{Vlumen} \cdot \mathrm{dap}_{\mathrm{lumen}}
\end{aligned}
$$

## Reactions

| Rate | Id | Name | Equation |
| --- | --- | --- | --- |
| $v_{\mathrm{DAPABS}}$ | `DAPABS` | absorption dapagliflozin (DAPABS) | $\mathrm{dap}_{\mathrm{lumen}} \longrightarrow \mathrm{dap}_{\mathrm{ext}}$ |
| $v_{\mathrm{DAPEXC}}$ | `DAPEXC` | excretion dapagliflozin (feces) | $\mathrm{dap}_{\mathrm{lumen}} \longrightarrow \mathrm{dap}_{\mathrm{feces}}$ |
| $v_{\mathrm{dissolution\_dap}}$ | `dissolution_dap` | dissolution dapagliflozin | $\mathrm{dap}_{\mathrm{stomach}} \longrightarrow \mathrm{dap}_{\mathrm{lumen}}$ |

The rates of the reactions are:

$$
\begin{aligned}
v_{\mathrm{DAPABS}} &= F_{\mathrm{dap\_abs}} \cdot \mathrm{absorption} \\[1ex]
v_{\mathrm{DAPEXC}} &= \mathopen{}\left(1 - F_{\mathrm{dap\_abs}}\right) \cdot \mathrm{absorption} \\[1ex]
v_{\mathrm{dissolution\_dap}} &= \frac{\frac{\mathrm{Ka}_{\mathrm{dis\_dap}}}{60} \cdot \mathrm{PODOSE}_{\mathrm{dap}}}{\mathrm{Mr}_{\mathrm{dap}}}
\end{aligned}
$$

## ODE system

The states change in time with the rates of the reactions and the rate rules:

$$
\begin{aligned}
\frac{\mathrm{d} \,\mathrm{dap}_{\mathrm{lumen}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{DAPABS}} - v_{\mathrm{DAPEXC}} + v_{\mathrm{dissolution\_dap}}}{\mathrm{Vlumen}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{dap}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{v_{\mathrm{DAPABS}}}{\mathrm{Vext}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{dap}_{\mathrm{feces}}}{\mathrm{d} t} &= v_{\mathrm{DAPEXC}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{PODOSE}_{\mathrm{dap}}}{\mathrm{d} t} &= -v_{\mathrm{dissolution\_dap}} \cdot \mathrm{Mr}_{\mathrm{dap}} \qquad \text{(rate rule)}
\end{aligned}
$$
