# Dapagliflozin body model

Model `dapagliflozin_body`, SBML Level 3 Version 2, read from dapagliflozin\_body\_flat.xml, written by sbmlode 0.3.0.

Physiologically based pharmacokinetic/ pharmacodynamic (PBPK/PD) model of dapagliflozin

Dapagliflozin is an SGLT2 inhibitor prescribed for the management of type 2 diabetes mellitus. The drug lowers blood glucose levels by increasing urinary glucose excretion. Despite established efficacy, dapagliflozin demonstrates significant inter-individual variability in pharmacokinetics and pharmacodynamics, with potential impact on treatment outcomes. To evaluate the sources of variability, we developed a digital twin of dapagliflozin pharmacokinetics and pharmacodynamics built as a PBPK/PD model using curated data from 28 clinical studies. This constitutes a full mechanistic model that integrates absorption, metabolism, excretion, and pharmacodynamics. It accounts for key determinants of variability including renal and hepatic function, and food effects.

Substances:

dap: dapagliflozin

d3g: dapagliflozin-3-O-glucuronide (D3G, about 61 % of an oral dose in the urine \[Obermeier2010\]), formed by UGT1A9 in liver and kidney

dapmet: the other metabolites, the 2-O-glucuronide (about 5 %) and the oxidative metabolites (about 10 %) \[Obermeier2010, Callegari2021\], lumped into one metabolite

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
| $\mathrm{Vre}$ | `Vre` | rest of body |  | l |  |
| $\mathrm{Vgu}$ | `Vgu` | gut |  | l |  |
| $\mathrm{Vlumen}$ | `Vlumen` | intestinal lumen (inner part of intestine) |  | l |  |
| $\mathrm{Vki}$ | `Vki` | kidney |  | l |  |
| $\mathrm{Vli}$ | `Vli` | liver |  | l |  |
| $\mathrm{Vlu}$ | `Vlu` | lung |  | l |  |
| $\mathrm{Vve}$ | `Vve` | venous blood |  | l |  |
| $\mathrm{Var}$ | `Var` | arterial blood |  | l |  |
| $\mathrm{Vurine}$ | `Vurine` | urine | $1$ | l | ✓ |
| $\mathrm{Vfeces}$ | `Vfeces` | feces | $1$ | l | ✓ |
| $\mathrm{Vpo}$ | `Vpo` | portal plasma |  | l |  |
| $\mathrm{Vhv}$ | `Vhv` | hepatic venous plasma |  | l |  |
| $\mathrm{Vki}_{\mathrm{tissue}}$ | `Vki_tissue` | ki tissue |  | l |  |
| $\mathrm{Vki}_{\mathrm{plasma}}$ | `Vki_plasma` | ki plasma |  | l |  |
| $\mathrm{Vli}_{\mathrm{tissue}}$ | `Vli_tissue` | li tissue |  | l |  |
| $\mathrm{Vli}_{\mathrm{plasma}}$ | `Vli_plasma` | li plasma |  | l |  |
| $\mathrm{Vlu}_{\mathrm{tissue}}$ | `Vlu_tissue` | lu tissue |  | l |  |
| $\mathrm{Vlu}_{\mathrm{plasma}}$ | `Vlu_plasma` | lu plasma |  | l |  |
| $\mathrm{Vgu}_{\mathrm{tissue}}$ | `Vgu_tissue` | gu tissue |  | l |  |
| $\mathrm{Vgu}_{\mathrm{plasma}}$ | `Vgu_plasma` | gu plasma |  | l |  |
| $\mathrm{Vre}_{\mathrm{tissue}}$ | `Vre_tissue` | re tissue |  | l |  |
| $\mathrm{Vre}_{\mathrm{plasma}}$ | `Vre_plasma` | re plasma |  | l |  |
| $\mathrm{GU}_{\mathrm{\_Vapical}}$ | `GU__Vapical` | apical membrane (intestinal membrane enterocytes) | $\mathrm{NaN}$ | m<sup>2</sup> | ✓ |
| $\mathrm{GU}_{\mathrm{\_Vstomach}}$ | `GU__Vstomach` | stomach | $1$ | l | ✓ |

## Species

| Symbol | Id | Name | Compartment | Value | Unit | Properties |
| --- | --- | --- | --- | --- | --- | --- |
| $\mathrm{Cki}_{\mathrm{plasma\_dap}}$ | `Cki_plasma_dap` | dapagliflozin (kidney plasma) | $\mathrm{Vki}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cli}_{\mathrm{plasma\_dap}}$ | `Cli_plasma_dap` | dapagliflozin (liver plasma) | $\mathrm{Vli}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Clu}_{\mathrm{plasma\_dap}}$ | `Clu_plasma_dap` | dapagliflozin (lung plasma) | $\mathrm{Vlu}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cgu}_{\mathrm{plasma\_dap}}$ | `Cgu_plasma_dap` | dapagliflozin (gut plasma) | $\mathrm{Vgu}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cre}_{\mathrm{plasma\_dap}}$ | `Cre_plasma_dap` | dapagliflozin (rest plasma) | $\mathrm{Vre}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Car}_{\mathrm{dap}}$ | `Car_dap` | dapagliflozin (arterial blood plasma) | $\mathrm{Var}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cve}_{\mathrm{dap}}$ | `Cve_dap` | dapagliflozin (venous blood plasma) | $\mathrm{Vve}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cpo}_{\mathrm{dap}}$ | `Cpo_dap` | dapagliflozin (portal vein plasma) | $\mathrm{Vpo}$ | $0$ | mmol/l | concentration |
| $\mathrm{Chv}_{\mathrm{dap}}$ | `Chv_dap` | dapagliflozin (hepatic vein plasma) | $\mathrm{Vhv}$ | $0$ | mmol/l | concentration |
| $\mathrm{Clu}_{\mathrm{dap}}$ | `Clu_dap` | dapagliflozin (lung) | $\mathrm{Vlu}_{\mathrm{tissue}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cre}_{\mathrm{dap}}$ | `Cre_dap` | dapagliflozin (rest) | $\mathrm{Vre}_{\mathrm{tissue}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cre}_{\mathrm{deep\_dap}}$ | `Cre_deep_dap` | dapagliflozin slow pool (rest) | $\mathrm{Vre}_{\mathrm{tissue}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Aurine}_{\mathrm{dap}}$ | `Aurine_dap` | dapagliflozin (urine) | $\mathrm{Vurine}$ |  | mmol | amount |
| $\mathrm{Afeces}_{\mathrm{dap}}$ | `Afeces_dap` | dapagliflozin (feces) | $\mathrm{Vfeces}$ |  | mmol | amount |
| $\mathrm{Cgu}_{\mathrm{dap}}$ | `Cgu_dap` | dapagliflozin (intestinal lumen) | $\mathrm{Vlumen}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cki}_{\mathrm{plasma\_d3g}}$ | `Cki_plasma_d3g` | dapagliflozin-3-O-glucuronide (kidney plasma) | $\mathrm{Vki}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cli}_{\mathrm{plasma\_d3g}}$ | `Cli_plasma_d3g` | dapagliflozin-3-O-glucuronide (liver plasma) | $\mathrm{Vli}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Clu}_{\mathrm{plasma\_d3g}}$ | `Clu_plasma_d3g` | dapagliflozin-3-O-glucuronide (lung plasma) | $\mathrm{Vlu}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cgu}_{\mathrm{plasma\_d3g}}$ | `Cgu_plasma_d3g` | dapagliflozin-3-O-glucuronide (gut plasma) | $\mathrm{Vgu}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cre}_{\mathrm{plasma\_d3g}}$ | `Cre_plasma_d3g` | dapagliflozin-3-O-glucuronide (rest plasma) | $\mathrm{Vre}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Car}_{\mathrm{d3g}}$ | `Car_d3g` | dapagliflozin-3-O-glucuronide (arterial blood plasma) | $\mathrm{Var}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cve}_{\mathrm{d3g}}$ | `Cve_d3g` | dapagliflozin-3-O-glucuronide (venous blood plasma) | $\mathrm{Vve}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cpo}_{\mathrm{d3g}}$ | `Cpo_d3g` | dapagliflozin-3-O-glucuronide (portal vein plasma) | $\mathrm{Vpo}$ | $0$ | mmol/l | concentration |
| $\mathrm{Chv}_{\mathrm{d3g}}$ | `Chv_d3g` | dapagliflozin-3-O-glucuronide (hepatic vein plasma) | $\mathrm{Vhv}$ | $0$ | mmol/l | concentration |
| $\mathrm{Aurine}_{\mathrm{d3g}}$ | `Aurine_d3g` | dapagliflozin-3-O-glucuronide (urine) | $\mathrm{Vurine}$ |  | mmol | amount |
| $\mathrm{Cki}_{\mathrm{plasma\_dapmet}}$ | `Cki_plasma_dapmet` | other dapagliflozin metabolites (kidney plasma) | $\mathrm{Vki}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cli}_{\mathrm{plasma\_dapmet}}$ | `Cli_plasma_dapmet` | other dapagliflozin metabolites (liver plasma) | $\mathrm{Vli}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Clu}_{\mathrm{plasma\_dapmet}}$ | `Clu_plasma_dapmet` | other dapagliflozin metabolites (lung plasma) | $\mathrm{Vlu}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cgu}_{\mathrm{plasma\_dapmet}}$ | `Cgu_plasma_dapmet` | other dapagliflozin metabolites (gut plasma) | $\mathrm{Vgu}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cre}_{\mathrm{plasma\_dapmet}}$ | `Cre_plasma_dapmet` | other dapagliflozin metabolites (rest plasma) | $\mathrm{Vre}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Car}_{\mathrm{dapmet}}$ | `Car_dapmet` | other dapagliflozin metabolites (arterial blood plasma) | $\mathrm{Var}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cve}_{\mathrm{dapmet}}$ | `Cve_dapmet` | other dapagliflozin metabolites (venous blood plasma) | $\mathrm{Vve}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cpo}_{\mathrm{dapmet}}$ | `Cpo_dapmet` | other dapagliflozin metabolites (portal vein plasma) | $\mathrm{Vpo}$ | $0$ | mmol/l | concentration |
| $\mathrm{Chv}_{\mathrm{dapmet}}$ | `Chv_dapmet` | other dapagliflozin metabolites (hepatic vein plasma) | $\mathrm{Vhv}$ | $0$ | mmol/l | concentration |
| $\mathrm{Aurine}_{\mathrm{dapmet}}$ | `Aurine_dapmet` | other dapagliflozin metabolites (urine) | $\mathrm{Vurine}$ |  | mmol | amount |
| $\mathrm{Cve}_{\mathrm{daptot}}$ | `Cve_daptot` | total dapagliflozin (plasma) | $\mathrm{Vve}$ |  | mmol/l | concentration |
| $\mathrm{Aurine}_{\mathrm{daptot}}$ | `Aurine_daptot` | total dapagliflozin (urine) | $\mathrm{Vurine}$ |  | mmol | amount |
| $\mathrm{Afeces}_{\mathrm{daptot}}$ | `Afeces_daptot` | total dapagliflozin (feces) | $\mathrm{Vfeces}$ |  | mmol | amount |
| $\mathrm{KI}_{\mathrm{\_fpg}}$ | `KI__fpg` | fasting plasma glucose (FPG) | $\mathrm{Vki}_{\mathrm{plasma}}$ | $5$ | mmol/l | concentration, boundary |
| $\mathrm{KI}_{\mathrm{\_glc\_urine}}$ | `KI__glc_urine` | glucose (urine) | $\mathrm{Vurine}$ |  | mmol | amount |
| $\mathrm{GU}_{\mathrm{\_dap\_stomach}}$ | `GU__dap_stomach` | dapagliflozin (stomach) | $\mathrm{GU}_{\mathrm{\_Vstomach}}$ |  | mmol | amount, boundary |

## Parameters

| Symbol | Id | Name | Value | Unit | Constant |
| --- | --- | --- | --- | --- | :---: |
| $\mathrm{BW}$ | `BW` | body weight | $75$ | kg | ✓ |
| $\mathrm{BW}_{\mathrm{ref}}$ | `BW_ref` | reference body weight \[kg\] | $75$ | kg | ✓ |
| $\mathrm{HEIGHT}$ | `HEIGHT` | height | $170$ | cm | ✓ |
| $\mathrm{BSA}$ | `BSA` | body surface area |  | m<sup>2</sup> |  |
| $\mathrm{COBW}$ | `COBW` | cardiac output per bodyweight | $1.548$ | ml/s/kg | ✓ |
| $f_{\mathrm{cardiac\_function}}$ | `f_cardiac_function` | heart function | $1$ | \- | ✓ |
| $\mathrm{CO}$ | `CO` | cardiac output \[ml/s\] |  | ml/s |  |
| $\mathrm{QC}$ | `QC` | cardiac output \[l/min\] |  | l/min |  |
| $\mathrm{Fblood}$ | `Fblood` | blood fraction of organ volume | $0.02$ | \- | ✓ |
| $\mathrm{HCT}$ | `HCT` | hematocrit | $0.51$ | \- | ✓ |
| $\mathrm{FVgu}$ | `FVgu` | gut fractional tissue volume | $0.0171$ | l/kg | ✓ |
| $\mathrm{FVki}$ | `FVki` | kidney fractional tissue volume | $0.0044$ | l/kg | ✓ |
| $\mathrm{FVli}$ | `FVli` | liver fractional tissue volume | $0.021$ | l/kg | ✓ |
| $\mathrm{FVlu}$ | `FVlu` | lung fractional tissue volume | $0.0076$ | l/kg | ✓ |
| $\mathrm{FVre}$ | `FVre` | rest of body fractional tissue volume |  | l/kg |  |
| $\mathrm{FVve}$ | `FVve` | venous fractional tissue volume | $0.0514$ | l/kg | ✓ |
| $\mathrm{FVar}$ | `FVar` | arterial fractional tissue volume | $0.0257$ | l/kg | ✓ |
| $\mathrm{FVpo}$ | `FVpo` | portal fractional tissue volume | $0.001$ | l/kg | ✓ |
| $\mathrm{FVhv}$ | `FVhv` | hepatic venous fractional tissue volume | $0.001$ | l/kg | ✓ |
| $\mathrm{FQgu}$ | `FQgu` | gut fractional tissue blood flow | $0.18$ | \- | ✓ |
| $\mathrm{FQki}$ | `FQki` | kidney fractional tissue blood flow | $0.19$ | \- | ✓ |
| $\mathrm{FQh}$ | `FQh` | hepatic (venous side) fractional tissue blood flow | $0.215$ | \- | ✓ |
| $\mathrm{FQlu}$ | `FQlu` | lung fractional tissue blood flow | $1$ | \- | ✓ |
| $\mathrm{FQre}$ | `FQre` | rest of body fractional tissue blood flow |  | \- |  |
| $f_{\mathrm{cirrhosis}}$ | `f_cirrhosis` | severity of cirrhosis \[0, 0.95\] | $0$ | \- | ✓ |
| $f_{\mathrm{shunts}}$ | `f_shunts` | fraction of portal venous blood shunted by the liver |  | \- |  |
| $f_{\mathrm{tissue\_loss}}$ | `f_tissue_loss` | fraction of lost parenchymal liver volume |  | \- |  |
| $f_{\mathrm{GFR\_Qki}}$ | `f_GFR_Qki` | relationship GFR and Qki | $76.3105006271061$ | ml/l | ✓ |
| $\mathrm{GFR}_{\mathrm{base}}$ | `GFR_base` | glomerular filtration rate (healthy) |  | ml/min |  |
| $f_{\mathrm{ugt1a9}}$ | `f_ugt1a9` | scaling factor UGT1A9 activity | $1$ | \- | ✓ |
| $\mathrm{DAPMET}_{\mathrm{k}}$ | `DAPMET_k` | rate of the hepatic metabolism of dapagliflozin | $0.0658712610363756$ | 1/min | ✓ |
| $f_{\mathrm{met\_renal}}$ | `f_met_renal` | renal fraction of the metabolism | $0.6$ | \- | ✓ |
| $\mathrm{DAPMET}_{\mathrm{k\_renal}}$ | `DAPMET_k_renal` | rate of the renal metabolism of dapagliflozin |  | 1/min |  |
| $f_{\mathrm{d3g}}$ | `f_d3g` | fraction of the metabolism to D3G | $0.77$ | \- | ✓ |
| $\mathrm{PODOSE}_{\mathrm{dap}}$ | `PODOSE_dap` | oral dose dap \[mg\] | $0$ | mg |  |
| $\mathrm{Mr}_{\mathrm{dap}}$ | `Mr_dap` | Molecular weight dap \[g/mole\] | $408.873$ | g/mol | ✓ |
| $\mathrm{ftissue}_{\mathrm{dap}}$ | `ftissue_dap` | tissue distribution dap | $0.411373243178683$ | l/min | ✓ |
| $\mathrm{Kp}_{\mathrm{dap}}$ | `Kp_dap` | tissue/plasma partition coefficient dap | $1.20017936506384$ | \- | ✓ |
| $\mathrm{Kpdeep}_{\mathrm{dap}}$ | `Kpdeep_dap` | slow pool/tissue partition coefficient dap (rest) | $0.5$ | \- | ✓ |
| $\mathrm{kdeep}_{\mathrm{dap}}$ | `kdeep_dap` | rate of the slow pool dap (rest) | $0.002$ | 1/min | ✓ |
| $\mathrm{IVDOSE}_{\mathrm{dap}}$ | `IVDOSE_dap` | IV bolus dose dap \[mg\] | $0$ | mg |  |
| $\mathrm{ti}_{\mathrm{dap}}$ | `ti_dap` | injection time dap \[s\] | $10$ | s | ✓ |
| $\mathrm{Ki}_{\mathrm{dap}}$ | `Ki_dap` | Ki \[1/min\] injection dap |  | 1/min |  |
| $\mathrm{Ri}_{\mathrm{dap}}$ | `Ri_dap` | Ri \[mg/min\] rate of infusion dap | $0$ | mg/min | ✓ |
| $\mathrm{Qgu}$ | `Qgu` | gut blood flow |  | l/min |  |
| $\mathrm{Qki}$ | `Qki` | kidney blood flow |  | l/min |  |
| $\mathrm{Qh}$ | `Qh` | hepatic (venous side) blood flow |  | l/min |  |
| $\mathrm{Qha}$ | `Qha` | hepatic artery blood flow |  | l/min |  |
| $\mathrm{Qlu}$ | `Qlu` | lung blood flow |  | l/min |  |
| $\mathrm{Qre}$ | `Qre` | rest of body blood flow |  | l/min |  |
| $\mathrm{Qpo}$ | `Qpo` | portal blood flow |  | l/min |  |
| $\mathrm{KI}_{\mathrm{\_f\_renal\_function}}$ | `KI__f_renal_function` | renal function | $1$ | \- | ✓ |
| $\mathrm{KI}_{\mathrm{\_fu\_dap}}$ | `KI__fu_dap` | fraction unbound in plasma dapagliflozin | $0.09$ | \- | ✓ |
| $\mathrm{KI}_{\mathrm{\_DAPEX\_f}}$ | `KI__DAPEX_f` | renal clearance of dapagliflozin relative to its filtration | $0.460206039705467$ | \- | ✓ |
| $\mathrm{KI}_{\mathrm{\_D3GEX\_f}}$ | `KI__D3GEX_f` | renal clearance of the metabolites relative to the GFR | $1.80545550669025$ | \- | ✓ |
| $\mathrm{KI}_{\mathrm{\_Mr\_glc}}$ | `KI__Mr_glc` | molecular weight glucose | $180.16$ | g/mol | ✓ |
| $\mathrm{KI}_{\mathrm{\_cf\_mg\_per\_g}}$ | `KI__cf_mg_per_g` | conversion factor mg per g | $1000$ | mg/g | ✓ |
| $\mathrm{KI}_{\mathrm{\_fpg\_healthy}}$ | `KI__fpg_healthy` | fasting plasma glucose (healthy) | $5$ | mmol/l | ✓ |
| $\mathrm{KI}_{\mathrm{\_PG\_delta}}$ | `KI__PG_delta` | mean plasma glucose of the day minus the fasting plasma glucose | $0.85$ | mmol/l | ✓ |
| $\mathrm{KI}_{\mathrm{\_RTG\_E50}}$ | `KI__RTG_E50` | EC50 reduction in RTG | $5.75629727963379 \times 10^{-6}$ | mmol/l | ✓ |
| $\mathrm{KI}_{\mathrm{\_RTG\_gamma}}$ | `KI__RTG_gamma` | hill coefficient reduction in RTG | $2.80107895262435$ | \- | ✓ |
| $\mathrm{KI}_{\mathrm{\_RTG\_base}}$ | `KI__RTG_base` | baseline RTG | $11$ | mmol/l | ✓ |
| $\mathrm{KI}_{\mathrm{\_RTG\_m\_fpg}}$ | `KI__RTG_m_fpg` | FPG effect on RTG | $1.99098619575577$ | \- | ✓ |
| $\mathrm{KI}_{\mathrm{\_RTG\_max\_inhibition}}$ | `KI__RTG_max_inhibition` | maximal relative reduction of the RTG | $0.640013975351734$ | \- | ✓ |
| $\mathrm{KI}_{\mathrm{\_RTG\_splay}}$ | `KI__RTG_splay` | splay of the glucose excretion around the RTG | $1$ | mmol/l | ✓ |
| $\mathrm{KI}_{\mathrm{\_GFR}}$ | `KI__GFR` | glomerular filtration rate |  | ml/min |  |
| $\mathrm{KI}_{\mathrm{\_pg}}$ | `KI__pg` | mean plasma glucose of the day |  | mmol/l |  |
| $\mathrm{KI}_{\mathrm{\_RTG\_fpg}}$ | `KI__RTG_fpg` | RTG without SGLT2 inhibition |  | mmol/l |  |
| $\mathrm{KI}_{\mathrm{\_RTG\_delta}}$ | `KI__RTG_delta` | maximal reduction of the RTG |  | mmol/l |  |
| $\mathrm{KI}_{\mathrm{\_RTG}}$ | `KI__RTG` | renal threshold glucose (RTG) |  | mmol/l |  |
| $\mathrm{KI}_{\mathrm{\_UGE}}$ | `KI__UGE` | urinary glucose excretion (UGE) |  | gram |  |
| $\mathrm{GU}_{\mathrm{\_F\_dap\_abs}}$ | `GU__F_dap_abs` | fraction absorbed dapagliflozin | $0.84$ | \- | ✓ |
| $\mathrm{GU}_{\mathrm{\_DAPABS\_k}}$ | `GU__DAPABS_k` | rate of dapagliflozin absorption | $0.0444918821768688$ | 1/min | ✓ |
| $\mathrm{GU}_{\mathrm{\_fed}}$ | `GU__fed` | fed state (0 fasted, 1 fed) | $0$ | \- | ✓ |
| $\mathrm{GU}_{\mathrm{\_FED\_f\_absorption}}$ | `GU__FED_f_absorption` | fed effect on the absorption rate | $0.189289061510997$ | \- | ✓ |
| $\mathrm{GU}_{\mathrm{\_f\_absorption}}$ | `GU__f_absorption` | factor of the absorption rate |  | \- |  |
| $\mathrm{GU}_{\mathrm{\_Ka\_dis\_dap}}$ | `GU__Ka_dis_dap` | dissolution rate dapagliflozin | $1.09176322622891$ | 1/hr | ✓ |
| $\mathrm{GU}_{\mathrm{\_absorption}}$ | `GU__absorption` | dapagliflozin leaving the lumen |  | mmol/min |  |

## Initial assignments and assignment rules

The initial assignments set the values at $t = 0$:

$$
\begin{aligned}
\mathrm{Aurine}_{\mathrm{dap}} &= 0 \cdot \mathrm{Vurine} \\[1ex]
\mathrm{Afeces}_{\mathrm{dap}} &= 0 \cdot \mathrm{Vfeces} \\[1ex]
\mathrm{Aurine}_{\mathrm{d3g}} &= 0 \cdot \mathrm{Vurine} \\[1ex]
\mathrm{Aurine}_{\mathrm{dapmet}} &= 0 \cdot \mathrm{Vurine} \\[1ex]
\mathrm{KI}_{\mathrm{\_glc\_urine}} &= 0 \cdot \mathrm{Vurine} \\[1ex]
\mathrm{GU}_{\mathrm{\_dap\_stomach}} &= 0 \cdot \mathrm{GU}_{\mathrm{\_Vstomach}}
\end{aligned}
$$

The assignment rules hold at every time $t$:

$$
\begin{aligned}
f_{\mathrm{shunts}} &= f_{\mathrm{cirrhosis}} \\[1ex]
f_{\mathrm{tissue\_loss}} &= f_{\mathrm{cirrhosis}} \\[1ex]
\mathrm{Cve}_{\mathrm{daptot}} &= \mathrm{Cve}_{\mathrm{dap}} + \mathrm{Cve}_{\mathrm{d3g}} + \mathrm{Cve}_{\mathrm{dapmet}} \\[1ex]
\mathrm{Aurine}_{\mathrm{daptot}} &= \mathrm{Aurine}_{\mathrm{dap}} + \mathrm{Aurine}_{\mathrm{d3g}} + \mathrm{Aurine}_{\mathrm{dapmet}} \\[1ex]
\mathrm{Afeces}_{\mathrm{daptot}} &= \mathrm{Afeces}_{\mathrm{dap}} \\[1ex]
\mathrm{DAPMET}_{\mathrm{k\_renal}} &= \frac{\frac{\mathrm{DAPMET}_{\mathrm{k}} \cdot \mathrm{FVli}}{\mathrm{FVki}} \cdot f_{\mathrm{met\_renal}}}{1 - f_{\mathrm{met\_renal}}} \\[1ex]
\mathrm{FVre} &= 1 - \mathrm{FVgu} - \mathrm{FVki} - \mathrm{FVli} \\
&\quad {} - \mathrm{FVlu} - \mathrm{FVve} - \mathrm{FVar} - \mathrm{FVpo} \\
&\quad {} - \mathrm{FVhv} \\[1ex]
\mathrm{FQre} &= 1 - \mathrm{FQki} - \mathrm{FQh} \\[1ex]
\mathrm{BSA} &= 0.024265 \cdot \mathopen{}\left(\frac{\mathrm{BW}}{1}\right)^{0.5378} \cdot \mathopen{}\left(\frac{\mathrm{HEIGHT}}{1}\right)^{0.3964} \\[1ex]
\mathrm{CO} &= f_{\mathrm{cardiac\_function}} \cdot \mathrm{BW} \cdot \mathrm{COBW} \\[1ex]
\mathrm{QC} &= \frac{\mathrm{CO}}{1000} \cdot 60 \\[1ex]
\mathrm{Vgu} &= \mathrm{BW} \cdot \mathrm{FVgu} \\[1ex]
\mathrm{Vki} &= \mathrm{BW} \cdot \mathrm{FVki} \\[1ex]
\mathrm{Vli} &= \mathrm{BW} \cdot \mathrm{FVli} \\[1ex]
\mathrm{Vlu} &= \mathrm{BW} \cdot \mathrm{FVlu} \\[1ex]
\mathrm{Vre} &= \mathrm{BW} \cdot \mathrm{FVre} \\[1ex]
\mathrm{Vve} &= \mathrm{BW} \cdot \mathrm{FVve} - \frac{\mathrm{FVve}}{\mathrm{FVar} + \mathrm{FVve}} \cdot \mathrm{BW} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{FVve} - \mathrm{FVar}\right) \\[1ex]
\mathrm{Var} &= \mathrm{BW} \cdot \mathrm{FVar} - \frac{\mathrm{FVar}}{\mathrm{FVar} + \mathrm{FVve}} \cdot \mathrm{BW} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{FVve} - \mathrm{FVar}\right) \\[1ex]
\mathrm{Vpo} &= \mathopen{}\left(1 - \mathrm{HCT}\right) \cdot \mathopen{}\left(\mathrm{BW} \cdot \mathrm{FVpo} - \frac{\mathrm{FVpo}}{\mathrm{FVar} + \mathrm{FVve} + \mathrm{FVpo} + \mathrm{FVhv}} \cdot \mathrm{BW} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{FVar} - \mathrm{FVve} - \mathrm{FVpo} - \mathrm{FVhv}\right)\right) \\[1ex]
\mathrm{Vhv} &= \mathopen{}\left(1 - \mathrm{HCT}\right) \cdot \mathopen{}\left(\mathrm{BW} \cdot \mathrm{FVhv} - \frac{\mathrm{FVhv}}{\mathrm{FVar} + \mathrm{FVve} + \mathrm{FVpo} + \mathrm{FVhv}} \cdot \mathrm{BW} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{FVar} - \mathrm{FVve} - \mathrm{FVpo} - \mathrm{FVhv}\right)\right) \\[1ex]
\mathrm{Qgu} &= \mathrm{QC} \cdot \mathrm{FQgu} \\[1ex]
\mathrm{Qki} &= \mathrm{QC} \cdot \mathrm{FQki} \\[1ex]
\mathrm{Qh} &= \mathrm{QC} \cdot \mathrm{FQh} \\[1ex]
\mathrm{Qha} &= \mathrm{Qh} - \mathrm{Qgu} \\[1ex]
\mathrm{Qlu} &= \mathrm{QC} \cdot \mathrm{FQlu} \\[1ex]
\mathrm{Qre} &= \mathrm{QC} \cdot \mathrm{FQre} \\[1ex]
\mathrm{Qpo} &= \mathrm{Qgu} \\[1ex]
\mathrm{GFR}_{\mathrm{base}} &= f_{\mathrm{GFR\_Qki}} \cdot \mathrm{Qki} \\[1ex]
\mathrm{Vlumen} &= 0.9 \cdot \mathrm{Vgu} \\[1ex]
\mathrm{Vki}_{\mathrm{plasma}} &= \mathrm{Vki} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{HCT}\right) \\[1ex]
\mathrm{Vki}_{\mathrm{tissue}} &= \mathrm{Vki} \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Vli}_{\mathrm{plasma}} &= \mathrm{Vli} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{HCT}\right) \\[1ex]
\mathrm{Vli}_{\mathrm{tissue}} &= \mathrm{Vli} \cdot \mathopen{}\left(1 - f_{\mathrm{tissue\_loss}}\right) \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Vlu}_{\mathrm{plasma}} &= \mathrm{Vlu} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{HCT}\right) \\[1ex]
\mathrm{Vlu}_{\mathrm{tissue}} &= \mathrm{Vlu} \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Vgu}_{\mathrm{plasma}} &= \mathrm{Vgu} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{HCT}\right) \\[1ex]
\mathrm{Vgu}_{\mathrm{tissue}} &= \mathrm{Vgu} \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Vre}_{\mathrm{plasma}} &= \mathrm{Vre} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{HCT}\right) \\[1ex]
\mathrm{Vre}_{\mathrm{tissue}} &= \mathrm{Vre} \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Ki}_{\mathrm{dap}} &= \frac{0.693}{\mathrm{ti}_{\mathrm{dap}}} \cdot 60 \\[1ex]
\mathrm{KI}_{\mathrm{\_GFR}} &= \mathrm{KI}_{\mathrm{\_f\_renal\_function}} \cdot \mathrm{GFR}_{\mathrm{base}} \\[1ex]
\mathrm{KI}_{\mathrm{\_pg}} &= \mathrm{KI}_{\mathrm{\_fpg}} + \mathrm{KI}_{\mathrm{\_PG\_delta}} \\[1ex]
\mathrm{KI}_{\mathrm{\_RTG\_fpg}} &= \mathrm{KI}_{\mathrm{\_RTG\_base}} + \mathrm{KI}_{\mathrm{\_RTG\_m\_fpg}} \cdot \mathopen{}\left(\mathrm{KI}_{\mathrm{\_fpg}} - \mathrm{KI}_{\mathrm{\_fpg\_healthy}}\right) \\[1ex]
\mathrm{KI}_{\mathrm{\_RTG\_delta}} &= \mathrm{KI}_{\mathrm{\_RTG\_fpg}} \cdot \mathrm{KI}_{\mathrm{\_RTG\_max\_inhibition}} \\[1ex]
\mathrm{KI}_{\mathrm{\_RTG}} &= \mathrm{KI}_{\mathrm{\_RTG\_fpg}} - \frac{\mathrm{KI}_{\mathrm{\_RTG\_delta}} \cdot \mathrm{Cki}_{\mathrm{plasma\_dap}}^{\mathrm{KI}_{\mathrm{\_RTG\_gamma}}}}{\mathrm{KI}_{\mathrm{\_RTG\_E50}}^{\mathrm{KI}_{\mathrm{\_RTG\_gamma}}} + \mathrm{Cki}_{\mathrm{plasma\_dap}}^{\mathrm{KI}_{\mathrm{\_RTG\_gamma}}}} \\[1ex]
\mathrm{KI}_{\mathrm{\_UGE}} &= \frac{\mathrm{KI}_{\mathrm{\_glc\_urine}} \cdot \mathrm{KI}_{\mathrm{\_Mr\_glc}}}{\mathrm{KI}_{\mathrm{\_cf\_mg\_per\_g}}} \\[1ex]
\mathrm{GU}_{\mathrm{\_f\_absorption}} &= 1 + \mathrm{GU}_{\mathrm{\_fed}} \cdot \mathopen{}\left(\mathrm{GU}_{\mathrm{\_FED\_f\_absorption}} - 1\right) \\[1ex]
\mathrm{GU}_{\mathrm{\_absorption}} &= \mathrm{GU}_{\mathrm{\_f\_absorption}} \cdot \mathrm{GU}_{\mathrm{\_DAPABS\_k}} \cdot \mathrm{Vlumen} \cdot \mathrm{Cgu}_{\mathrm{dap}}
\end{aligned}
$$

## Reactions

| Rate | Id | Name | Equation | Modifiers |
| --- | --- | --- | --- | --- |
| $v_{\mathrm{transport\_lu\_dap}}$ | `transport_lu_dap` | transport dapagliflozin | $\mathrm{Clu}_{\mathrm{plasma\_dap}} \rightleftharpoons \mathrm{Clu}_{\mathrm{dap}}$ |  |
| $v_{\mathrm{transport\_re\_dap}}$ | `transport_re_dap` | transport dapagliflozin | $\mathrm{Cre}_{\mathrm{plasma\_dap}} \rightleftharpoons \mathrm{Cre}_{\mathrm{dap}}$ |  |
| $v_{\mathrm{deep\_re\_dap}}$ | `deep_re_dap` | slow pool dapagliflozin (tissue) | $\mathrm{Cre}_{\mathrm{dap}} \rightleftharpoons \mathrm{Cre}_{\mathrm{deep\_dap}}$ |  |
| $v_{\mathrm{iv\_dap}}$ | `iv_dap` | iv dapagliflozin | $\varnothing \longrightarrow \mathrm{Cve}_{\mathrm{dap}}$ |  |
| $v_{\mathrm{Flow\_ar\_ki\_dap}}$ | `Flow_ar_ki_dap` | inflow kidney dapagliflozin | $\mathrm{Car}_{\mathrm{dap}} \longrightarrow \mathrm{Cki}_{\mathrm{plasma\_dap}}$ |  |
| $v_{\mathrm{Flow\_ki\_ve\_dap}}$ | `Flow_ki_ve_dap` | outflow kidney dapagliflozin | $\mathrm{Cki}_{\mathrm{plasma\_dap}} \longrightarrow \mathrm{Cve}_{\mathrm{dap}}$ |  |
| $v_{\mathrm{Flow\_arli\_li\_dap}}$ | `Flow_arli_li_dap` | arterial inflow liver dapagliflozin | $\mathrm{Car}_{\mathrm{dap}} \longrightarrow \mathrm{Cli}_{\mathrm{plasma\_dap}}$ |  |
| $v_{\mathrm{Flow\_arli\_hv\_dap}}$ | `Flow_arli_hv_dap` | flow arterial shunts dapagliflozin | $\mathrm{Car}_{\mathrm{dap}} \longrightarrow \mathrm{Chv}_{\mathrm{dap}}$ |  |
| $v_{\mathrm{Flow\_po\_li\_dap}}$ | `Flow_po_li_dap` | portal inflow liver dapagliflozin | $\mathrm{Cpo}_{\mathrm{dap}} \longrightarrow \mathrm{Cli}_{\mathrm{plasma\_dap}}$ |  |
| $v_{\mathrm{Flow\_po\_hv\_dap}}$ | `Flow_po_hv_dap` | flow portal shunts dapagliflozin | $\mathrm{Cpo}_{\mathrm{dap}} \longrightarrow \mathrm{Chv}_{\mathrm{dap}}$ |  |
| $v_{\mathrm{Flow\_li\_hv\_dap}}$ | `Flow_li_hv_dap` | outflow liver dapagliflozin | $\mathrm{Cli}_{\mathrm{plasma\_dap}} \longrightarrow \mathrm{Chv}_{\mathrm{dap}}$ |  |
| $v_{\mathrm{Flow\_hv\_ve\_dap}}$ | `Flow_hv_ve_dap` | outflow hepatic vein dapagliflozin | $\mathrm{Chv}_{\mathrm{dap}} \longrightarrow \mathrm{Cve}_{\mathrm{dap}}$ |  |
| $v_{\mathrm{Flow\_ve\_lu\_dap}}$ | `Flow_ve_lu_dap` | inflow lung dapagliflozin | $\mathrm{Cve}_{\mathrm{dap}} \longrightarrow \mathrm{Clu}_{\mathrm{plasma\_dap}}$ |  |
| $v_{\mathrm{Flow\_lu\_ar\_dap}}$ | `Flow_lu_ar_dap` | outflow lung dapagliflozin | $\mathrm{Clu}_{\mathrm{plasma\_dap}} \longrightarrow \mathrm{Car}_{\mathrm{dap}}$ |  |
| $v_{\mathrm{Flow\_ar\_gu\_dap}}$ | `Flow_ar_gu_dap` | inflow gut dapagliflozin | $\mathrm{Car}_{\mathrm{dap}} \longrightarrow \mathrm{Cgu}_{\mathrm{plasma\_dap}}$ |  |
| $v_{\mathrm{Flow\_gu\_po\_dap}}$ | `Flow_gu_po_dap` | outflow gut dapagliflozin | $\mathrm{Cgu}_{\mathrm{plasma\_dap}} \longrightarrow \mathrm{Cpo}_{\mathrm{dap}}$ |  |
| $v_{\mathrm{Flow\_ar\_re\_dap}}$ | `Flow_ar_re_dap` | inflow rest dapagliflozin | $\mathrm{Car}_{\mathrm{dap}} \longrightarrow \mathrm{Cre}_{\mathrm{plasma\_dap}}$ |  |
| $v_{\mathrm{Flow\_re\_ve\_dap}}$ | `Flow_re_ve_dap` | outflow rest dapagliflozin | $\mathrm{Cre}_{\mathrm{plasma\_dap}} \longrightarrow \mathrm{Cve}_{\mathrm{dap}}$ |  |
| $v_{\mathrm{Flow\_ar\_ki\_d3g}}$ | `Flow_ar_ki_d3g` | inflow kidney dapagliflozin-3-O-glucuronide | $\mathrm{Car}_{\mathrm{d3g}} \longrightarrow \mathrm{Cki}_{\mathrm{plasma\_d3g}}$ |  |
| $v_{\mathrm{Flow\_ki\_ve\_d3g}}$ | `Flow_ki_ve_d3g` | outflow kidney dapagliflozin-3-O-glucuronide | $\mathrm{Cki}_{\mathrm{plasma\_d3g}} \longrightarrow \mathrm{Cve}_{\mathrm{d3g}}$ |  |
| $v_{\mathrm{Flow\_arli\_li\_d3g}}$ | `Flow_arli_li_d3g` | arterial inflow liver dapagliflozin-3-O-glucuronide | $\mathrm{Car}_{\mathrm{d3g}} \longrightarrow \mathrm{Cli}_{\mathrm{plasma\_d3g}}$ |  |
| $v_{\mathrm{Flow\_arli\_hv\_d3g}}$ | `Flow_arli_hv_d3g` | flow arterial shunts dapagliflozin-3-O-glucuronide | $\mathrm{Car}_{\mathrm{d3g}} \longrightarrow \mathrm{Chv}_{\mathrm{d3g}}$ |  |
| $v_{\mathrm{Flow\_po\_li\_d3g}}$ | `Flow_po_li_d3g` | portal inflow liver dapagliflozin-3-O-glucuronide | $\mathrm{Cpo}_{\mathrm{d3g}} \longrightarrow \mathrm{Cli}_{\mathrm{plasma\_d3g}}$ |  |
| $v_{\mathrm{Flow\_po\_hv\_d3g}}$ | `Flow_po_hv_d3g` | flow portal shunts dapagliflozin-3-O-glucuronide | $\mathrm{Cpo}_{\mathrm{d3g}} \longrightarrow \mathrm{Chv}_{\mathrm{d3g}}$ |  |
| $v_{\mathrm{Flow\_li\_hv\_d3g}}$ | `Flow_li_hv_d3g` | outflow liver dapagliflozin-3-O-glucuronide | $\mathrm{Cli}_{\mathrm{plasma\_d3g}} \longrightarrow \mathrm{Chv}_{\mathrm{d3g}}$ |  |
| $v_{\mathrm{Flow\_hv\_ve\_d3g}}$ | `Flow_hv_ve_d3g` | outflow hepatic vein dapagliflozin-3-O-glucuronide | $\mathrm{Chv}_{\mathrm{d3g}} \longrightarrow \mathrm{Cve}_{\mathrm{d3g}}$ |  |
| $v_{\mathrm{Flow\_ve\_lu\_d3g}}$ | `Flow_ve_lu_d3g` | inflow lung dapagliflozin-3-O-glucuronide | $\mathrm{Cve}_{\mathrm{d3g}} \longrightarrow \mathrm{Clu}_{\mathrm{plasma\_d3g}}$ |  |
| $v_{\mathrm{Flow\_lu\_ar\_d3g}}$ | `Flow_lu_ar_d3g` | outflow lung dapagliflozin-3-O-glucuronide | $\mathrm{Clu}_{\mathrm{plasma\_d3g}} \longrightarrow \mathrm{Car}_{\mathrm{d3g}}$ |  |
| $v_{\mathrm{Flow\_ar\_gu\_d3g}}$ | `Flow_ar_gu_d3g` | inflow gut dapagliflozin-3-O-glucuronide | $\mathrm{Car}_{\mathrm{d3g}} \longrightarrow \mathrm{Cgu}_{\mathrm{plasma\_d3g}}$ |  |
| $v_{\mathrm{Flow\_gu\_po\_d3g}}$ | `Flow_gu_po_d3g` | outflow gut dapagliflozin-3-O-glucuronide | $\mathrm{Cgu}_{\mathrm{plasma\_d3g}} \longrightarrow \mathrm{Cpo}_{\mathrm{d3g}}$ |  |
| $v_{\mathrm{Flow\_ar\_re\_d3g}}$ | `Flow_ar_re_d3g` | inflow rest dapagliflozin-3-O-glucuronide | $\mathrm{Car}_{\mathrm{d3g}} \longrightarrow \mathrm{Cre}_{\mathrm{plasma\_d3g}}$ |  |
| $v_{\mathrm{Flow\_re\_ve\_d3g}}$ | `Flow_re_ve_d3g` | outflow rest dapagliflozin-3-O-glucuronide | $\mathrm{Cre}_{\mathrm{plasma\_d3g}} \longrightarrow \mathrm{Cve}_{\mathrm{d3g}}$ |  |
| $v_{\mathrm{Flow\_ar\_ki\_dapmet}}$ | `Flow_ar_ki_dapmet` | inflow kidney other dapagliflozin metabolites | $\mathrm{Car}_{\mathrm{dapmet}} \longrightarrow \mathrm{Cki}_{\mathrm{plasma\_dapmet}}$ |  |
| $v_{\mathrm{Flow\_ki\_ve\_dapmet}}$ | `Flow_ki_ve_dapmet` | outflow kidney other dapagliflozin metabolites | $\mathrm{Cki}_{\mathrm{plasma\_dapmet}} \longrightarrow \mathrm{Cve}_{\mathrm{dapmet}}$ |  |
| $v_{\mathrm{Flow\_arli\_li\_dapmet}}$ | `Flow_arli_li_dapmet` | arterial inflow liver other dapagliflozin metabolites | $\mathrm{Car}_{\mathrm{dapmet}} \longrightarrow \mathrm{Cli}_{\mathrm{plasma\_dapmet}}$ |  |
| $v_{\mathrm{Flow\_arli\_hv\_dapmet}}$ | `Flow_arli_hv_dapmet` | flow arterial shunts other dapagliflozin metabolites | $\mathrm{Car}_{\mathrm{dapmet}} \longrightarrow \mathrm{Chv}_{\mathrm{dapmet}}$ |  |
| $v_{\mathrm{Flow\_po\_li\_dapmet}}$ | `Flow_po_li_dapmet` | portal inflow liver other dapagliflozin metabolites | $\mathrm{Cpo}_{\mathrm{dapmet}} \longrightarrow \mathrm{Cli}_{\mathrm{plasma\_dapmet}}$ |  |
| $v_{\mathrm{Flow\_po\_hv\_dapmet}}$ | `Flow_po_hv_dapmet` | flow portal shunts other dapagliflozin metabolites | $\mathrm{Cpo}_{\mathrm{dapmet}} \longrightarrow \mathrm{Chv}_{\mathrm{dapmet}}$ |  |
| $v_{\mathrm{Flow\_li\_hv\_dapmet}}$ | `Flow_li_hv_dapmet` | outflow liver other dapagliflozin metabolites | $\mathrm{Cli}_{\mathrm{plasma\_dapmet}} \longrightarrow \mathrm{Chv}_{\mathrm{dapmet}}$ |  |
| $v_{\mathrm{Flow\_hv\_ve\_dapmet}}$ | `Flow_hv_ve_dapmet` | outflow hepatic vein other dapagliflozin metabolites | $\mathrm{Chv}_{\mathrm{dapmet}} \longrightarrow \mathrm{Cve}_{\mathrm{dapmet}}$ |  |
| $v_{\mathrm{Flow\_ve\_lu\_dapmet}}$ | `Flow_ve_lu_dapmet` | inflow lung other dapagliflozin metabolites | $\mathrm{Cve}_{\mathrm{dapmet}} \longrightarrow \mathrm{Clu}_{\mathrm{plasma\_dapmet}}$ |  |
| $v_{\mathrm{Flow\_lu\_ar\_dapmet}}$ | `Flow_lu_ar_dapmet` | outflow lung other dapagliflozin metabolites | $\mathrm{Clu}_{\mathrm{plasma\_dapmet}} \longrightarrow \mathrm{Car}_{\mathrm{dapmet}}$ |  |
| $v_{\mathrm{Flow\_ar\_gu\_dapmet}}$ | `Flow_ar_gu_dapmet` | inflow gut other dapagliflozin metabolites | $\mathrm{Car}_{\mathrm{dapmet}} \longrightarrow \mathrm{Cgu}_{\mathrm{plasma\_dapmet}}$ |  |
| $v_{\mathrm{Flow\_gu\_po\_dapmet}}$ | `Flow_gu_po_dapmet` | outflow gut other dapagliflozin metabolites | $\mathrm{Cgu}_{\mathrm{plasma\_dapmet}} \longrightarrow \mathrm{Cpo}_{\mathrm{dapmet}}$ |  |
| $v_{\mathrm{Flow\_ar\_re\_dapmet}}$ | `Flow_ar_re_dapmet` | inflow rest other dapagliflozin metabolites | $\mathrm{Car}_{\mathrm{dapmet}} \longrightarrow \mathrm{Cre}_{\mathrm{plasma\_dapmet}}$ |  |
| $v_{\mathrm{Flow\_re\_ve\_dapmet}}$ | `Flow_re_ve_dapmet` | outflow rest other dapagliflozin metabolites | $\mathrm{Cre}_{\mathrm{plasma\_dapmet}} \longrightarrow \mathrm{Cve}_{\mathrm{dapmet}}$ |  |
| $v_{\mathrm{KI\_\_DAP2D3G}}$ | `KI__DAP2D3G` | dapagliflozin glucuronidation (DAP2D3G, UGT1A9) | $\mathrm{Cki}_{\mathrm{plasma\_dap}} \longrightarrow \mathrm{Cki}_{\mathrm{plasma\_d3g}}$ |  |
| $v_{\mathrm{KI\_\_DAP2MET}}$ | `KI__DAP2MET` | dapagliflozin other metabolism (DAP2MET) | $\mathrm{Cki}_{\mathrm{plasma\_dap}} \longrightarrow \mathrm{Cki}_{\mathrm{plasma\_dapmet}}$ |  |
| $v_{\mathrm{KI\_\_DAPEX}}$ | `KI__DAPEX` | dapagliflozin excretion (DAPEX) | $\mathrm{Cki}_{\mathrm{plasma\_dap}} \longrightarrow \mathrm{Aurine}_{\mathrm{dap}}$ |  |
| $v_{\mathrm{KI\_\_D3GEX}}$ | `KI__D3GEX` | D3G excretion (D3GEX) | $\mathrm{Cki}_{\mathrm{plasma\_d3g}} \longrightarrow \mathrm{Aurine}_{\mathrm{d3g}}$ |  |
| $v_{\mathrm{KI\_\_METEX}}$ | `KI__METEX` | other metabolites excretion (METEX) | $\mathrm{Cki}_{\mathrm{plasma\_dapmet}} \longrightarrow \mathrm{Aurine}_{\mathrm{dapmet}}$ |  |
| $v_{\mathrm{KI\_\_GLCEX}}$ | `KI__GLCEX` | glucose excretion (GLCEX) | $\mathrm{KI}_{\mathrm{\_fpg}} \longrightarrow \mathrm{KI}_{\mathrm{\_glc\_urine}}$ | $\mathrm{Cki}_{\mathrm{plasma\_dap}}$ |
| $v_{\mathrm{LI\_\_DAP2D3G}}$ | `LI__DAP2D3G` | dapagliflozin glucuronidation (DAP2D3G, UGT1A9) | $\mathrm{Cli}_{\mathrm{plasma\_dap}} \longrightarrow \mathrm{Cli}_{\mathrm{plasma\_d3g}}$ |  |
| $v_{\mathrm{LI\_\_DAP2MET}}$ | `LI__DAP2MET` | dapagliflozin other metabolism (DAP2MET) | $\mathrm{Cli}_{\mathrm{plasma\_dap}} \longrightarrow \mathrm{Cli}_{\mathrm{plasma\_dapmet}}$ |  |
| $v_{\mathrm{GU\_\_DAPABS}}$ | `GU__DAPABS` | absorption dapagliflozin (DAPABS) | $\mathrm{Cgu}_{\mathrm{dap}} \longrightarrow \mathrm{Cgu}_{\mathrm{plasma\_dap}}$ |  |
| $v_{\mathrm{GU\_\_DAPEXC}}$ | `GU__DAPEXC` | excretion dapagliflozin (feces) | $\mathrm{Cgu}_{\mathrm{dap}} \longrightarrow \mathrm{Afeces}_{\mathrm{dap}}$ |  |
| $v_{\mathrm{GU\_\_dissolution\_dap}}$ | `GU__dissolution_dap` | dissolution dapagliflozin | $\mathrm{GU}_{\mathrm{\_dap\_stomach}} \longrightarrow \mathrm{Cgu}_{\mathrm{dap}}$ |  |

The rates of the reactions are:

$$
\begin{aligned}
v_{\mathrm{transport\_lu\_dap}} &= \frac{\mathrm{ftissue}_{\mathrm{dap}} \cdot \mathrm{BW}}{\mathrm{BW}_{\mathrm{ref}}} \cdot \mathopen{}\left(\mathrm{Clu}_{\mathrm{plasma\_dap}} \cdot \mathrm{Kp}_{\mathrm{dap}} - \mathrm{Clu}_{\mathrm{dap}}\right) \\[1ex]
v_{\mathrm{transport\_re\_dap}} &= \frac{\mathrm{ftissue}_{\mathrm{dap}} \cdot \mathrm{BW}}{\mathrm{BW}_{\mathrm{ref}}} \cdot \mathopen{}\left(\mathrm{Cre}_{\mathrm{plasma\_dap}} \cdot \mathrm{Kp}_{\mathrm{dap}} - \mathrm{Cre}_{\mathrm{dap}}\right) \\[1ex]
v_{\mathrm{deep\_re\_dap}} &= \mathrm{kdeep}_{\mathrm{dap}} \cdot \mathrm{Vre}_{\mathrm{tissue}} \cdot \mathopen{}\left(\mathrm{Kpdeep}_{\mathrm{dap}} \cdot \mathrm{Cre}_{\mathrm{dap}} - \mathrm{Cre}_{\mathrm{deep\_dap}}\right) \\[1ex]
v_{\mathrm{iv\_dap}} &= \frac{\mathrm{Ki}_{\mathrm{dap}} \cdot \mathrm{IVDOSE}_{\mathrm{dap}}}{\mathrm{Mr}_{\mathrm{dap}}} \\[1ex]
v_{\mathrm{Flow\_ar\_ki\_dap}} &= \mathrm{Qki} \cdot \mathrm{Car}_{\mathrm{dap}} \\[1ex]
v_{\mathrm{Flow\_ki\_ve\_dap}} &= \mathrm{Qki} \cdot \mathrm{Cki}_{\mathrm{plasma\_dap}} \\[1ex]
v_{\mathrm{Flow\_arli\_li\_dap}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathrm{Qha} \cdot \mathrm{Car}_{\mathrm{dap}} \\[1ex]
v_{\mathrm{Flow\_arli\_hv\_dap}} &= f_{\mathrm{shunts}} \cdot \mathrm{Qha} \cdot \mathrm{Car}_{\mathrm{dap}} \\[1ex]
v_{\mathrm{Flow\_po\_li\_dap}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathrm{Qpo} \cdot \mathrm{Cpo}_{\mathrm{dap}} \\[1ex]
v_{\mathrm{Flow\_po\_hv\_dap}} &= f_{\mathrm{shunts}} \cdot \mathrm{Qpo} \cdot \mathrm{Cpo}_{\mathrm{dap}} \\[1ex]
v_{\mathrm{Flow\_li\_hv\_dap}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathopen{}\left(\mathrm{Qpo} + \mathrm{Qha}\right) \cdot \mathrm{Cli}_{\mathrm{plasma\_dap}} \\[1ex]
v_{\mathrm{Flow\_hv\_ve\_dap}} &= \mathrm{Qh} \cdot \mathrm{Chv}_{\mathrm{dap}} \\[1ex]
v_{\mathrm{Flow\_ve\_lu\_dap}} &= \mathrm{Qlu} \cdot \mathrm{Cve}_{\mathrm{dap}} \\[1ex]
v_{\mathrm{Flow\_lu\_ar\_dap}} &= \mathrm{Qlu} \cdot \mathrm{Clu}_{\mathrm{plasma\_dap}} \\[1ex]
v_{\mathrm{Flow\_ar\_gu\_dap}} &= \mathrm{Qgu} \cdot \mathrm{Car}_{\mathrm{dap}} \\[1ex]
v_{\mathrm{Flow\_gu\_po\_dap}} &= \mathrm{Qgu} \cdot \mathrm{Cgu}_{\mathrm{plasma\_dap}} \\[1ex]
v_{\mathrm{Flow\_ar\_re\_dap}} &= \mathrm{Qre} \cdot \mathrm{Car}_{\mathrm{dap}} \\[1ex]
v_{\mathrm{Flow\_re\_ve\_dap}} &= \mathrm{Qre} \cdot \mathrm{Cre}_{\mathrm{plasma\_dap}} \\[1ex]
v_{\mathrm{Flow\_ar\_ki\_d3g}} &= \mathrm{Qki} \cdot \mathrm{Car}_{\mathrm{d3g}} \\[1ex]
v_{\mathrm{Flow\_ki\_ve\_d3g}} &= \mathrm{Qki} \cdot \mathrm{Cki}_{\mathrm{plasma\_d3g}} \\[1ex]
v_{\mathrm{Flow\_arli\_li\_d3g}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathrm{Qha} \cdot \mathrm{Car}_{\mathrm{d3g}} \\[1ex]
v_{\mathrm{Flow\_arli\_hv\_d3g}} &= f_{\mathrm{shunts}} \cdot \mathrm{Qha} \cdot \mathrm{Car}_{\mathrm{d3g}} \\[1ex]
v_{\mathrm{Flow\_po\_li\_d3g}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathrm{Qpo} \cdot \mathrm{Cpo}_{\mathrm{d3g}} \\[1ex]
v_{\mathrm{Flow\_po\_hv\_d3g}} &= f_{\mathrm{shunts}} \cdot \mathrm{Qpo} \cdot \mathrm{Cpo}_{\mathrm{d3g}} \\[1ex]
v_{\mathrm{Flow\_li\_hv\_d3g}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathopen{}\left(\mathrm{Qpo} + \mathrm{Qha}\right) \cdot \mathrm{Cli}_{\mathrm{plasma\_d3g}} \\[1ex]
v_{\mathrm{Flow\_hv\_ve\_d3g}} &= \mathrm{Qh} \cdot \mathrm{Chv}_{\mathrm{d3g}} \\[1ex]
v_{\mathrm{Flow\_ve\_lu\_d3g}} &= \mathrm{Qlu} \cdot \mathrm{Cve}_{\mathrm{d3g}} \\[1ex]
v_{\mathrm{Flow\_lu\_ar\_d3g}} &= \mathrm{Qlu} \cdot \mathrm{Clu}_{\mathrm{plasma\_d3g}} \\[1ex]
v_{\mathrm{Flow\_ar\_gu\_d3g}} &= \mathrm{Qgu} \cdot \mathrm{Car}_{\mathrm{d3g}} \\[1ex]
v_{\mathrm{Flow\_gu\_po\_d3g}} &= \mathrm{Qgu} \cdot \mathrm{Cgu}_{\mathrm{plasma\_d3g}} \\[1ex]
v_{\mathrm{Flow\_ar\_re\_d3g}} &= \mathrm{Qre} \cdot \mathrm{Car}_{\mathrm{d3g}} \\[1ex]
v_{\mathrm{Flow\_re\_ve\_d3g}} &= \mathrm{Qre} \cdot \mathrm{Cre}_{\mathrm{plasma\_d3g}} \\[1ex]
v_{\mathrm{Flow\_ar\_ki\_dapmet}} &= \mathrm{Qki} \cdot \mathrm{Car}_{\mathrm{dapmet}} \\[1ex]
v_{\mathrm{Flow\_ki\_ve\_dapmet}} &= \mathrm{Qki} \cdot \mathrm{Cki}_{\mathrm{plasma\_dapmet}} \\[1ex]
v_{\mathrm{Flow\_arli\_li\_dapmet}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathrm{Qha} \cdot \mathrm{Car}_{\mathrm{dapmet}} \\[1ex]
v_{\mathrm{Flow\_arli\_hv\_dapmet}} &= f_{\mathrm{shunts}} \cdot \mathrm{Qha} \cdot \mathrm{Car}_{\mathrm{dapmet}} \\[1ex]
v_{\mathrm{Flow\_po\_li\_dapmet}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathrm{Qpo} \cdot \mathrm{Cpo}_{\mathrm{dapmet}} \\[1ex]
v_{\mathrm{Flow\_po\_hv\_dapmet}} &= f_{\mathrm{shunts}} \cdot \mathrm{Qpo} \cdot \mathrm{Cpo}_{\mathrm{dapmet}} \\[1ex]
v_{\mathrm{Flow\_li\_hv\_dapmet}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathopen{}\left(\mathrm{Qpo} + \mathrm{Qha}\right) \cdot \mathrm{Cli}_{\mathrm{plasma\_dapmet}} \\[1ex]
v_{\mathrm{Flow\_hv\_ve\_dapmet}} &= \mathrm{Qh} \cdot \mathrm{Chv}_{\mathrm{dapmet}} \\[1ex]
v_{\mathrm{Flow\_ve\_lu\_dapmet}} &= \mathrm{Qlu} \cdot \mathrm{Cve}_{\mathrm{dapmet}} \\[1ex]
v_{\mathrm{Flow\_lu\_ar\_dapmet}} &= \mathrm{Qlu} \cdot \mathrm{Clu}_{\mathrm{plasma\_dapmet}} \\[1ex]
v_{\mathrm{Flow\_ar\_gu\_dapmet}} &= \mathrm{Qgu} \cdot \mathrm{Car}_{\mathrm{dapmet}} \\[1ex]
v_{\mathrm{Flow\_gu\_po\_dapmet}} &= \mathrm{Qgu} \cdot \mathrm{Cgu}_{\mathrm{plasma\_dapmet}} \\[1ex]
v_{\mathrm{Flow\_ar\_re\_dapmet}} &= \mathrm{Qre} \cdot \mathrm{Car}_{\mathrm{dapmet}} \\[1ex]
v_{\mathrm{Flow\_re\_ve\_dapmet}} &= \mathrm{Qre} \cdot \mathrm{Cre}_{\mathrm{plasma\_dapmet}} \\[1ex]
v_{\mathrm{KI\_\_DAP2D3G}} &= \mathrm{KI}_{\mathrm{\_f\_renal\_function}} \cdot f_{\mathrm{ugt1a9}} \cdot f_{\mathrm{d3g}} \cdot \mathrm{DAPMET}_{\mathrm{k\_renal}} \cdot \mathrm{Vki}_{\mathrm{tissue}} \cdot \mathrm{Cki}_{\mathrm{plasma\_dap}} \\[1ex]
v_{\mathrm{KI\_\_DAP2MET}} &= \mathrm{KI}_{\mathrm{\_f\_renal\_function}} \cdot \mathopen{}\left(1 - f_{\mathrm{d3g}}\right) \cdot \mathrm{DAPMET}_{\mathrm{k\_renal}} \cdot \mathrm{Vki}_{\mathrm{tissue}} \cdot \mathrm{Cki}_{\mathrm{plasma\_dap}} \\[1ex]
v_{\mathrm{KI\_\_DAPEX}} &= \mathrm{KI}_{\mathrm{\_GFR}} \cdot 0.001 \cdot \mathrm{KI}_{\mathrm{\_fu\_dap}} \cdot \mathrm{KI}_{\mathrm{\_DAPEX\_f}} \cdot \mathrm{Cki}_{\mathrm{plasma\_dap}} \\[1ex]
v_{\mathrm{KI\_\_D3GEX}} &= \mathrm{KI}_{\mathrm{\_GFR}} \cdot 0.001 \cdot \mathrm{KI}_{\mathrm{\_D3GEX\_f}} \cdot \mathrm{Cki}_{\mathrm{plasma\_d3g}} \\[1ex]
v_{\mathrm{KI\_\_METEX}} &= \mathrm{KI}_{\mathrm{\_GFR}} \cdot 0.001 \cdot \mathrm{KI}_{\mathrm{\_D3GEX\_f}} \cdot \mathrm{Cki}_{\mathrm{plasma\_dapmet}} \\[1ex]
v_{\mathrm{KI\_\_GLCEX}} &= \mathrm{KI}_{\mathrm{\_GFR}} \cdot 0.001 \cdot \mathrm{KI}_{\mathrm{\_RTG\_splay}} \cdot \ln\mathopen{}\left(1 + \exp\mathopen{}\left(\frac{\mathrm{KI}_{\mathrm{\_pg}} - \mathrm{KI}_{\mathrm{\_RTG}}}{\mathrm{KI}_{\mathrm{\_RTG\_splay}}}\right)\right) \\[1ex]
v_{\mathrm{LI\_\_DAP2D3G}} &= f_{\mathrm{ugt1a9}} \cdot f_{\mathrm{d3g}} \cdot \mathrm{DAPMET}_{\mathrm{k}} \cdot \mathrm{Vli}_{\mathrm{tissue}} \cdot \mathrm{Cli}_{\mathrm{plasma\_dap}} \\[1ex]
v_{\mathrm{LI\_\_DAP2MET}} &= \mathopen{}\left(1 - f_{\mathrm{d3g}}\right) \cdot \mathrm{DAPMET}_{\mathrm{k}} \cdot \mathrm{Vli}_{\mathrm{tissue}} \cdot \mathrm{Cli}_{\mathrm{plasma\_dap}} \\[1ex]
v_{\mathrm{GU\_\_DAPABS}} &= \mathrm{GU}_{\mathrm{\_F\_dap\_abs}} \cdot \mathrm{GU}_{\mathrm{\_absorption}} \\[1ex]
v_{\mathrm{GU\_\_DAPEXC}} &= \mathopen{}\left(1 - \mathrm{GU}_{\mathrm{\_F\_dap\_abs}}\right) \cdot \mathrm{GU}_{\mathrm{\_absorption}} \\[1ex]
v_{\mathrm{GU\_\_dissolution\_dap}} &= \frac{\frac{\mathrm{GU}_{\mathrm{\_Ka\_dis\_dap}}}{60} \cdot \mathrm{PODOSE}_{\mathrm{dap}}}{\mathrm{Mr}_{\mathrm{dap}}}
\end{aligned}
$$

## ODE system

The states change in time with the rates of the reactions and the rate rules:

$$
\begin{aligned}
\frac{\mathrm{d} \,\mathrm{Cki}_{\mathrm{plasma\_dap}}}{\mathrm{d} t} &= \frac{1}{\mathrm{Vki}_{\mathrm{plasma}}} \Bigl( v_{\mathrm{Flow\_ar\_ki\_dap}} - v_{\mathrm{Flow\_ki\_ve\_dap}} - v_{\mathrm{KI\_\_DAP2D3G}} - v_{\mathrm{KI\_\_DAP2MET}} \\
&\quad {} - v_{\mathrm{KI\_\_DAPEX}} \Bigr) \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cli}_{\mathrm{plasma\_dap}}}{\mathrm{d} t} &= \frac{1}{\mathrm{Vli}_{\mathrm{plasma}}} \Bigl( v_{\mathrm{Flow\_arli\_li\_dap}} + v_{\mathrm{Flow\_po\_li\_dap}} - v_{\mathrm{Flow\_li\_hv\_dap}} - v_{\mathrm{LI\_\_DAP2D3G}} \\
&\quad {} - v_{\mathrm{LI\_\_DAP2MET}} \Bigr) \\[1ex]
\frac{\mathrm{d} \,\mathrm{Clu}_{\mathrm{plasma\_dap}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{transport\_lu\_dap}} + v_{\mathrm{Flow\_ve\_lu\_dap}} - v_{\mathrm{Flow\_lu\_ar\_dap}}}{\mathrm{Vlu}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cgu}_{\mathrm{plasma\_dap}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_gu\_dap}} - v_{\mathrm{Flow\_gu\_po\_dap}} + v_{\mathrm{GU\_\_DAPABS}}}{\mathrm{Vgu}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cre}_{\mathrm{plasma\_dap}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{transport\_re\_dap}} + v_{\mathrm{Flow\_ar\_re\_dap}} - v_{\mathrm{Flow\_re\_ve\_dap}}}{\mathrm{Vre}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Car}_{\mathrm{dap}}}{\mathrm{d} t} &= \frac{1}{\mathrm{Var}} \Bigl( -v_{\mathrm{Flow\_ar\_ki\_dap}} - v_{\mathrm{Flow\_arli\_li\_dap}} - v_{\mathrm{Flow\_arli\_hv\_dap}} + v_{\mathrm{Flow\_lu\_ar\_dap}} \\
&\quad {} - v_{\mathrm{Flow\_ar\_gu\_dap}} - v_{\mathrm{Flow\_ar\_re\_dap}} \Bigr) \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cve}_{\mathrm{dap}}}{\mathrm{d} t} &= \frac{1}{\mathrm{Vve}} \Bigl( v_{\mathrm{iv\_dap}} + v_{\mathrm{Flow\_ki\_ve\_dap}} + v_{\mathrm{Flow\_hv\_ve\_dap}} - v_{\mathrm{Flow\_ve\_lu\_dap}} \\
&\quad {} + v_{\mathrm{Flow\_re\_ve\_dap}} \Bigr) \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cpo}_{\mathrm{dap}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{Flow\_po\_li\_dap}} - v_{\mathrm{Flow\_po\_hv\_dap}} + v_{\mathrm{Flow\_gu\_po\_dap}}}{\mathrm{Vpo}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Chv}_{\mathrm{dap}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_arli\_hv\_dap}} + v_{\mathrm{Flow\_po\_hv\_dap}} + v_{\mathrm{Flow\_li\_hv\_dap}} - v_{\mathrm{Flow\_hv\_ve\_dap}}}{\mathrm{Vhv}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Clu}_{\mathrm{dap}}}{\mathrm{d} t} &= \frac{v_{\mathrm{transport\_lu\_dap}}}{\mathrm{Vlu}_{\mathrm{tissue}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cre}_{\mathrm{dap}}}{\mathrm{d} t} &= \frac{v_{\mathrm{transport\_re\_dap}} - v_{\mathrm{deep\_re\_dap}}}{\mathrm{Vre}_{\mathrm{tissue}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cre}_{\mathrm{deep\_dap}}}{\mathrm{d} t} &= \frac{v_{\mathrm{deep\_re\_dap}}}{\mathrm{Vre}_{\mathrm{tissue}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Aurine}_{\mathrm{dap}}}{\mathrm{d} t} &= v_{\mathrm{KI\_\_DAPEX}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Afeces}_{\mathrm{dap}}}{\mathrm{d} t} &= v_{\mathrm{GU\_\_DAPEXC}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cgu}_{\mathrm{dap}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{GU\_\_DAPABS}} - v_{\mathrm{GU\_\_DAPEXC}} + v_{\mathrm{GU\_\_dissolution\_dap}}}{\mathrm{Vlumen}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cki}_{\mathrm{plasma\_d3g}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_ki\_d3g}} - v_{\mathrm{Flow\_ki\_ve\_d3g}} + v_{\mathrm{KI\_\_DAP2D3G}} - v_{\mathrm{KI\_\_D3GEX}}}{\mathrm{Vki}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cli}_{\mathrm{plasma\_d3g}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_arli\_li\_d3g}} + v_{\mathrm{Flow\_po\_li\_d3g}} - v_{\mathrm{Flow\_li\_hv\_d3g}} + v_{\mathrm{LI\_\_DAP2D3G}}}{\mathrm{Vli}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Clu}_{\mathrm{plasma\_d3g}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ve\_lu\_d3g}} - v_{\mathrm{Flow\_lu\_ar\_d3g}}}{\mathrm{Vlu}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cgu}_{\mathrm{plasma\_d3g}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_gu\_d3g}} - v_{\mathrm{Flow\_gu\_po\_d3g}}}{\mathrm{Vgu}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cre}_{\mathrm{plasma\_d3g}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_re\_d3g}} - v_{\mathrm{Flow\_re\_ve\_d3g}}}{\mathrm{Vre}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Car}_{\mathrm{d3g}}}{\mathrm{d} t} &= \frac{1}{\mathrm{Var}} \Bigl( -v_{\mathrm{Flow\_ar\_ki\_d3g}} - v_{\mathrm{Flow\_arli\_li\_d3g}} - v_{\mathrm{Flow\_arli\_hv\_d3g}} + v_{\mathrm{Flow\_lu\_ar\_d3g}} \\
&\quad {} - v_{\mathrm{Flow\_ar\_gu\_d3g}} - v_{\mathrm{Flow\_ar\_re\_d3g}} \Bigr) \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cve}_{\mathrm{d3g}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ki\_ve\_d3g}} + v_{\mathrm{Flow\_hv\_ve\_d3g}} - v_{\mathrm{Flow\_ve\_lu\_d3g}} + v_{\mathrm{Flow\_re\_ve\_d3g}}}{\mathrm{Vve}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cpo}_{\mathrm{d3g}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{Flow\_po\_li\_d3g}} - v_{\mathrm{Flow\_po\_hv\_d3g}} + v_{\mathrm{Flow\_gu\_po\_d3g}}}{\mathrm{Vpo}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Chv}_{\mathrm{d3g}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_arli\_hv\_d3g}} + v_{\mathrm{Flow\_po\_hv\_d3g}} + v_{\mathrm{Flow\_li\_hv\_d3g}} - v_{\mathrm{Flow\_hv\_ve\_d3g}}}{\mathrm{Vhv}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Aurine}_{\mathrm{d3g}}}{\mathrm{d} t} &= v_{\mathrm{KI\_\_D3GEX}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cki}_{\mathrm{plasma\_dapmet}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_ki\_dapmet}} - v_{\mathrm{Flow\_ki\_ve\_dapmet}} + v_{\mathrm{KI\_\_DAP2MET}} - v_{\mathrm{KI\_\_METEX}}}{\mathrm{Vki}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cli}_{\mathrm{plasma\_dapmet}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_arli\_li\_dapmet}} + v_{\mathrm{Flow\_po\_li\_dapmet}} - v_{\mathrm{Flow\_li\_hv\_dapmet}} + v_{\mathrm{LI\_\_DAP2MET}}}{\mathrm{Vli}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Clu}_{\mathrm{plasma\_dapmet}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ve\_lu\_dapmet}} - v_{\mathrm{Flow\_lu\_ar\_dapmet}}}{\mathrm{Vlu}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cgu}_{\mathrm{plasma\_dapmet}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_gu\_dapmet}} - v_{\mathrm{Flow\_gu\_po\_dapmet}}}{\mathrm{Vgu}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cre}_{\mathrm{plasma\_dapmet}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_re\_dapmet}} - v_{\mathrm{Flow\_re\_ve\_dapmet}}}{\mathrm{Vre}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Car}_{\mathrm{dapmet}}}{\mathrm{d} t} &= \frac{1}{\mathrm{Var}} \Bigl( -v_{\mathrm{Flow\_ar\_ki\_dapmet}} - v_{\mathrm{Flow\_arli\_li\_dapmet}} - v_{\mathrm{Flow\_arli\_hv\_dapmet}} + v_{\mathrm{Flow\_lu\_ar\_dapmet}} \\
&\quad {} - v_{\mathrm{Flow\_ar\_gu\_dapmet}} - v_{\mathrm{Flow\_ar\_re\_dapmet}} \Bigr) \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cve}_{\mathrm{dapmet}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ki\_ve\_dapmet}} + v_{\mathrm{Flow\_hv\_ve\_dapmet}} - v_{\mathrm{Flow\_ve\_lu\_dapmet}} + v_{\mathrm{Flow\_re\_ve\_dapmet}}}{\mathrm{Vve}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cpo}_{\mathrm{dapmet}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{Flow\_po\_li\_dapmet}} - v_{\mathrm{Flow\_po\_hv\_dapmet}} + v_{\mathrm{Flow\_gu\_po\_dapmet}}}{\mathrm{Vpo}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Chv}_{\mathrm{dapmet}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_arli\_hv\_dapmet}} + v_{\mathrm{Flow\_po\_hv\_dapmet}} + v_{\mathrm{Flow\_li\_hv\_dapmet}} - v_{\mathrm{Flow\_hv\_ve\_dapmet}}}{\mathrm{Vhv}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Aurine}_{\mathrm{dapmet}}}{\mathrm{d} t} &= v_{\mathrm{KI\_\_METEX}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{KI}_{\mathrm{\_glc\_urine}}}{\mathrm{d} t} &= v_{\mathrm{KI\_\_GLCEX}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{PODOSE}_{\mathrm{dap}}}{\mathrm{d} t} &= -v_{\mathrm{GU\_\_dissolution\_dap}} \cdot \mathrm{Mr}_{\mathrm{dap}} \qquad \text{(rate rule)} \\[1ex]
\frac{\mathrm{d} \,\mathrm{IVDOSE}_{\mathrm{dap}}}{\mathrm{d} t} &= -v_{\mathrm{iv\_dap}} \cdot \mathrm{Mr}_{\mathrm{dap}} + \mathrm{Ri}_{\mathrm{dap}} \qquad \text{(rate rule)}
\end{aligned}
$$
