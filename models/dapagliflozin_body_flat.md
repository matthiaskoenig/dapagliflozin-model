# model: dapagliflozin_body
Autogenerated ODE System from SBML with [sbmlutils](https://github.com/matthiaskoenig/sbmlutils).
```
time: [min]
substance: [mmol]
extent: [mmol]
volume: [l]
area: [m^2]
length: [m]
```

## Parameters `p`
```
BW = 75.0  # [kg] body weight [kg]  
COBW = 1.548  # [ml/s/kg] cardiac output per bodyweight [ml/s/kg]  
COHRI = 150.0  # [ml] increase of cardiac output per heartbeat [ml/min*min]  
FQgu = 0.18  # [-] gut fractional tissue blood flow  
FQh = 0.215  # [-] hepatic (venous side) fractional tissue blood flow  
FQki = 0.19  # [-] kidney fractional tissue blood flow  
FQlu = 1.0  # [-] lung fractional tissue blood flow  
FVar = 0.0257  # [l/kg] arterial fractional tissue volume  
FVgu = 0.0171  # [l/kg] gut fractional tissue volume  
FVhv = 0.001  # [l/kg] hepatic venous fractional tissue volume  
FVki = 0.0044  # [l/kg] kidney fractional tissue volume  
FVli = 0.021  # [l/kg] liver fractional tissue volume  
FVlu = 0.0076  # [l/kg] lung fractional tissue volume  
FVpo = 0.001  # [l/kg] portal fractional tissue volume  
FVve = 0.0514  # [l/kg] venous fractional tissue volume  
Fblood = 0.02  # [-] blood fraction of organ volume  
GU__DAPABS_k = 0.1  # [1/min] rate of dapagliflozin absorption  
GU__F_dap_abs = 0.84  # [-] fraction absorbed dapagliflozin  
GU__Ka_dis_dap = 2.0  # [1/hr] Ka_dis [1/hr] dissolution dapagliflozin  
GU__Mr_dap = 408.873  # [g/mol] Molecular weight dapagliflozin [g/mole]  
GU__Vapical = nan  # [m^2] apical membrane (intestinal membrane enterocytes)  
GU__Vbaso = nan  # [m^2] basolateral membrane (intestinal membrane enterocytes)  
GU__Ventero = 1.0  # [l] intestinal lining (enterocytes)  
GU__Vlumen = 1.15425  # [l] intestinal lumen (inner part of intestine)  
GU__Vstomach = 1.0  # [l] stomach  
GU__f_absorption = 1.0  # [-] scaling factor for absorption rate  
HCT = 0.51  # [-] hematocrit  
HEIGHT = 170.0  # [cm] height [cm]  
HR = 70.0  # [1/min] heart rate [1/min]  
HRrest = 70.0  # [1/min] heart rate [1/min]  
KI__D3GEX_k = 1.0  # [1/min] rate urinary excretion of dapagliflozin-3-o-glucuronide  
KI__DAPEX_k = 1.0  # [1/min] rate urinary excretion of dapagliflozin  
KI__GFR_healthy = 100.0  # [ml/min] Glomerular filtration rate (healthy)  
KI__Mr_glc = 180.0  # [g/mol] Molecular weight glc [g/mole]  
KI__RTG_E50 = 2.5e-06  # [mmol/l] EC50 reduction in RTG  
KI__RTG_base = 12.5  # [mmol/l] Baseline RTG value  
KI__RTG_delta = 9.0  # [mmol/l] RTG value  
KI__RTG_gamma = 1.0  # [-] hill coefficient reduction in RTG  
KI__Vmem = nan  # [m^2] plasma membrane  
KI__cf_mg_per_g = 1000.0  # [mg/g] Conversion factor mg per g  
KI__cf_ml_per_l = 1000.0  # [ml/l] Conversion factor ml per l  
KI__f_renal_function = 1.0  # [-] parameter for renal function  
Kp_dap = 10.0  # [-] tissue/plasma partition coefficient dap  
LI__D3GEX_Vmax = 1000.0  # [1/min] Vmax dapagliflozinat export  
LI__DAP2D3G_Km_dap = 0.479  # [mmol/l] Km dapagliflozin UGT1A9  
LI__DAP2D3G_Vmax = 0.04  # [mmol/min/l] Vmax dapagliflozin conversion  
LI__DAPIM_Vmax = 1000.0  # [1/min] Vmax dapagliflozin import  
LI__Vmem = nan  # [m^2] plasma membrane  
LI__f_ugt1a9 = 1.0  # [-] scaling factor UGT1A9 activity  
MAP = 100.0  # [133.32239 N/m^2] mean arterial pressure [mmHg]  
Mr_d3g = 585.0  # [g/mol] Molecular weight d3g [g/mole]  
Mr_dap = 408.873  # [g/mol] Molecular weight dap [g/mole]  
Ri_dap = 0.0  # [mg/min] Ri [mg/min] rate of infusion dap  
Vfeces = 1.0  # [l] feces  
Vstomach = 1.0  # [l] stomach  
Vurine = 1.0  # [l] urine  
conversion_min_per_day = 1440.0  # [min/day] Conversion factor min to hours  
f_cardiac_function = 1.0  # [-] heart function  
f_cirrhosis = 0.0  # [-] severity of cirrhosis [0, 0.95]  
ftissue_dap = 1.0  # [l/min] tissue distribution dap  
ti_dap = 10.0  # [s] injection time dap [s]  
```

## Initial conditions `x0`
```
Afeces_dap = 0.0  # [mmol] dapagliflozin (feces) in Vfeces  
Aurine_d3g = 0.0  # [mmol] dapagliflozin-3-o-glucuronide (urine) in Vurine  
Aurine_dap = 0.0  # [mmol] dapagliflozin (urine) in Vurine  
Car_d3g = 0.0  # [mmol/l] dapagliflozin-3-o-glucuronide (arterial blood plasma) in Var  
Car_dap = 0.0  # [mmol/l] dapagliflozin (arterial blood plasma) in Var  
Cgu_dap = 0.0  # [mmol/l] dapagliflozin (gut) in Vgu  
Cgu_plasma_d3g = 0.0  # [mmol/l] dapagliflozin-3-o-glucuronide (gut plasma) in Vgu_plasma  
Cgu_plasma_dap = 0.0  # [mmol/l] dapagliflozin (gut plasma) in Vgu_plasma  
Chv_d3g = 0.0  # [mmol/l] dapagliflozin-3-o-glucuronide (hepatic vein plasma) in Vhv  
Chv_dap = 0.0  # [mmol/l] dapagliflozin (hepatic vein plasma) in Vhv  
Cki_plasma_d3g = 0.0  # [mmol/l] dapagliflozin-3-o-glucuronide (kidney plasma) in Vki_plasma  
Cki_plasma_dap = 0.0  # [mmol/l] dapagliflozin (kidney plasma) in Vki_plasma  
Cli_plasma_d3g = 0.0  # [mmol/l] dapagliflozin-3-o-glucuronide (liver plasma) in Vli_plasma  
Cli_plasma_dap = 0.0  # [mmol/l] dapagliflozin (liver plasma) in Vli_plasma  
Clu_dap = 0.0  # [mmol/l] dapagliflozin (lung) in Vlu_tissue  
Clu_plasma_d3g = 0.0  # [mmol/l] dapagliflozin-3-o-glucuronide (lung plasma) in Vlu_plasma  
Clu_plasma_dap = 0.0  # [mmol/l] dapagliflozin (lung plasma) in Vlu_plasma  
Cpo_d3g = 0.0  # [mmol/l] dapagliflozin-3-o-glucuronide (portal vein plasma) in Vpo  
Cpo_dap = 0.0  # [mmol/l] dapagliflozin (portal vein plasma) in Vpo  
Cre_dap = 0.0  # [mmol/l] dapagliflozin (rest) in Vre_tissue  
Cre_plasma_d3g = 0.0  # [mmol/l] dapagliflozin-3-o-glucuronide (rest plasma) in Vre_plasma  
Cre_plasma_dap = 0.0  # [mmol/l] dapagliflozin (rest plasma) in Vre_plasma  
Cve_d3g = 0.0  # [mmol/l] dapagliflozin-3-o-glucuronide (venous blood plasma) in Vve  
Cve_dap = 0.0  # [mmol/l] dapagliflozin (venous blood plasma) in Vve  
GU__dap_stomach = 0.0  # [mmol] dapagliflozin (stomach) in GU__Vstomach  
IVDOSE_dap = 0.0  # [mg] IV bolus dose dap [mg]  
KI__glc_ext = 5.0  # [mmol/l] glucose (plasma) in Vki_plasma  
KI__glc_urine = 0.0  # [mmol] glucose (urine) in Vurine  
LI__d3g = 0.0  # [mmol/l] dapagliflozin-3-o-glucuronide (liver) in Vli_tissue  
LI__dap = 0.0  # [mmol/l] dapagliflozin (liver) in Vli_tissue  
PODOSE_dap = 0.0  # [mg] oral dose dap [mg]  
cum_dose_dap = 0.0  # [mg] Cumulative dose due to infusion dap  
```

## ODE system
```
# y
Aurine_daptot = Aurine_dap + Aurine_d3g  # [mmol] Sum of dapagliflozin and d3g urine  
BSA = 0.024265 * (BW / 1)**0.5378 * (HEIGHT / 1)**0.3964  # [m^2] body surface area [m^2]  
CO = f_cardiac_function * BW * COBW + (HR - HRrest) * COHRI / 60  # [ml/s] cardiac output [ml/s]  
Cve_daptot = Cve_dap + Cve_d3g  # [mmol/l] Sum of dapagliflozin and d3g  
FQre = 1 - (FQki + FQh)  # [-] rest of body fractional tissue blood flow  
FVre = 1 - (FVgu + FVki + FVli + FVlu + FVve + FVar)  # [l/kg] rest of body fractional tissue volume  
GU__dissolution_dap = (GU__Ka_dis_dap / 60) * PODOSE_dap / GU__Mr_dap  # [mmol/min] dissolution dapagliflozin  
KI__GFR = KI__f_renal_function * KI__GFR_healthy  # [ml/min] glomerular filtration rate  
KI__RTG = KI__RTG_base - KI__RTG_delta * Cki_plasma_dap**KI__RTG_gamma / (KI__RTG_E50**KI__RTG_gamma + Cki_plasma_dap**KI__RTG_gamma)  # [mmol/l] renal threshold glucose (RTG)  
KI__UGE = KI__glc_urine * KI__Mr_glc / KI__cf_mg_per_g  # [gram] urinary glucose excretion (UGE)  
Ki_dap = (0.693 / ti_dap) * 60  # [1/min] injection rate IV  
Var = BW * FVar - (FVar / (FVar + FVve)) * BW * Fblood * (1 - FVve - FVar)  # [l] arterial blood  
Vgu = BW * FVgu  # [l] gut  
Vhv = (1 - HCT) * (BW * FVhv - (FVhv / (FVar + FVve + FVpo + FVhv)) * BW * Fblood * (1 - (FVar + FVve + FVpo + FVhv)))  # [l] hepatic venous plasma  
Vki = BW * FVki  # [l] kidney  
Vli = BW * FVli  # [l] liver  
Vlu = BW * FVlu  # [l] lung  
Vpo = (1 - HCT) * (BW * FVpo - (FVpo / (FVar + FVve + FVpo + FVhv)) * BW * Fblood * (1 - (FVar + FVve + FVpo + FVhv)))  # [l] portal plasma  
Vve = BW * FVve - (FVve / (FVar + FVve)) * BW * Fblood * (1 - FVve - FVar)  # [l] venous blood  
Xfeces_dap = Afeces_dap * Mr_dap  # [mg] dapagliflozin amount (feces)  
Xurine_d3g = Aurine_d3g * Mr_d3g  # [mg] dapagliflozin-3-o-glucuronide amount (urine) [mg]  
Xurine_dap = Aurine_dap * Mr_dap  # [mg] dapagliflozin amount (urine) [mg]  
f_shunts = f_cirrhosis  # [-] fraction of portal venous blood shunted by the liver  
f_tissue_loss = f_cirrhosis  # [-] fraction of lost parenchymal liver volume  
transport_lu_dap = ftissue_dap * (Clu_plasma_dap * Kp_dap - Clu_dap)  # [mmol/min] transport dapagliflozin  
transport_re_dap = ftissue_dap * (Cre_plasma_dap * Kp_dap - Cre_dap)  # [mmol/min] transport dapagliflozin  
Aar_d3g = Car_d3g * Var  # [mmol] dapagliflozin-3-o-glucuronide amount (arterial blood) [mmole]  
Aar_dap = Car_dap * Var  # [mmol] dapagliflozin amount (arterial blood) [mmole]  
Ahv_d3g = Chv_d3g * Vhv  # [mmol] dapagliflozin-3-o-glucuronide amount (hepatic vein) [mmole]  
Ahv_dap = Chv_dap * Vhv  # [mmol] dapagliflozin amount (hepatic vein) [mmole]  
Apo_d3g = Cpo_d3g * Vpo  # [mmol] dapagliflozin-3-o-glucuronide amount (portal vein) [mmole]  
Apo_dap = Cpo_dap * Vpo  # [mmol] dapagliflozin amount (portal vein) [mmole]  
Ave_d3g = Cve_d3g * Vve  # [mmol] dapagliflozin-3-o-glucuronide amount (venous blood) [mmole]  
Ave_dap = Cve_dap * Vve  # [mmol] dapagliflozin amount (venous blood) [mmole]  
GU__absorption = GU__f_absorption * GU__DAPABS_k * Vgu * Cgu_dap  # [mmol/min] absorption dapagliflozin  
KI__GLCEX = piecewise((KI__GFR / KI__cf_ml_per_l) * (KI__glc_ext - KI__RTG), KI__glc_ext > KI__RTG, 0)  # [mmol/min] glucose excretion (GLCEX)  
QC = (CO / 1000) * 60  # [l/min] cardiac output [L/hr]  
Vgu_plasma = Vgu * Fblood * (1 - HCT)  # [l] plasma volume of gut  
Vgu_tissue = Vgu * (1 - Fblood)  # [l] tissue volume of gut  
Vki_plasma = Vki * Fblood * (1 - HCT)  # [l] plasma volume of kidney  
Vki_tissue = Vki * (1 - Fblood)  # [l] tissue volume of kidney  
Vli_plasma = Vli * Fblood * (1 - HCT)  # [l] plasma volume of liver  
Vli_tissue = Vli * (1 - f_tissue_loss) * (1 - Fblood)  # [l] tissue volume of liver  
Vlu_plasma = Vlu * Fblood * (1 - HCT)  # [l] plasma volume of lung  
Vlu_tissue = Vlu * (1 - Fblood)  # [l] tissue volume of lung  
Vre = BW * FVre  # [l] rest of body  
iv_dap = Ki_dap * IVDOSE_dap / Mr_dap  # [mmol/min] iv dapagliflozin  
Agu_plasma_d3g = Cgu_plasma_d3g * Vgu_plasma  # [mmol] dapagliflozin-3-o-glucuronide amount (gut) [mmole]  
Agu_plasma_dap = Cgu_plasma_dap * Vgu_plasma  # [mmol] dapagliflozin amount (gut) [mmole]  
Aki_plasma_d3g = Cki_plasma_d3g * Vki_plasma  # [mmol] dapagliflozin-3-o-glucuronide amount (kidney) [mmole]  
Aki_plasma_dap = Cki_plasma_dap * Vki_plasma  # [mmol] dapagliflozin amount (kidney) [mmole]  
Ali_plasma_d3g = Cli_plasma_d3g * Vli_plasma  # [mmol] dapagliflozin-3-o-glucuronide amount (liver) [mmole]  
Ali_plasma_dap = Cli_plasma_dap * Vli_plasma  # [mmol] dapagliflozin amount (liver) [mmole]  
Alu_plasma_d3g = Clu_plasma_d3g * Vlu_plasma  # [mmol] dapagliflozin-3-o-glucuronide amount (lung) [mmole]  
Alu_plasma_dap = Clu_plasma_dap * Vlu_plasma  # [mmol] dapagliflozin amount (lung) [mmole]  
GU__DAPABS = GU__F_dap_abs * GU__absorption  # [mmol/min] absorption dapagliflozin  
GU__DAPEXC = (1 - GU__F_dap_abs) * GU__absorption  # [mmol/min] excretion dapagliflozin (feces)  
KI__D3GEX = KI__f_renal_function * KI__D3GEX_k * Vki_tissue * Cki_plasma_d3g  # [mmol/min] dapagliflozin-3-o-glucuronide excretion (D3GEX)  
KI__DAPEX = KI__f_renal_function * KI__DAPEX_k * Vki_tissue * Cki_plasma_dap  # [mmol/min] dapagliflozin excretion (DAPEX)  
LI__D3GEX = LI__D3GEX_Vmax * Vli_tissue * (LI__d3g - Cli_plasma_d3g)  # [mmol/min] dapagliflozin-3-o-glucuronide export (D3GEX)  
LI__DAP2D3G = LI__f_ugt1a9 * LI__DAP2D3G_Vmax * Vli_tissue * LI__dap / (LI__dap + LI__DAP2D3G_Km_dap)  # [mmol/min] dapagliflozin conversion (DAP2D3G) UGT1A9  
LI__DAPIM = LI__DAPIM_Vmax * Vli_tissue * Cli_plasma_dap  # [mmol/min] dapagliflozin import (DAPIM) OTP1B1  
Mar_d3g = (Aar_d3g / Var) * Mr_d3g  # [mg/l] dapagliflozin-3-o-glucuronide concentration (arterial blood) [mg/l]  
Mar_dap = (Aar_dap / Var) * Mr_dap  # [mg/l] dapagliflozin concentration (arterial blood) [mg/l]  
Mhv_d3g = (Ahv_d3g / Vhv) * Mr_d3g  # [mg/l] dapagliflozin-3-o-glucuronide concentration (hepatic vein) [mg/l]  
Mhv_dap = (Ahv_dap / Vhv) * Mr_dap  # [mg/l] dapagliflozin concentration (hepatic vein) [mg/l]  
Mpo_d3g = (Apo_d3g / Vpo) * Mr_d3g  # [mg/l] dapagliflozin-3-o-glucuronide concentration (portal vein) [mg/l]  
Mpo_dap = (Apo_dap / Vpo) * Mr_dap  # [mg/l] dapagliflozin concentration (portal vein) [mg/l]  
Mve_d3g = (Ave_d3g / Vve) * Mr_d3g  # [mg/l] dapagliflozin-3-o-glucuronide concentration (venous blood) [mg/l]  
Mve_dap = (Ave_dap / Vve) * Mr_dap  # [mg/l] dapagliflozin concentration (venous blood) [mg/l]  
Qgu = QC * FQgu  # [l/min] gut blood flow  
Qh = QC * FQh  # [l/min] hepatic (venous side) blood flow  
Qki = QC * FQki  # [l/min] kidney blood flow  
Qlu = QC * FQlu  # [l/min] lung blood flow  
Qre = QC * FQre  # [l/min] rest of body blood flow  
Vre_plasma = Vre * Fblood * (1 - HCT)  # [l] plasma volume of rest  
Vre_tissue = Vre * (1 - Fblood)  # [l] tissue volume of rest  
Xar_d3g = Aar_d3g * Mr_d3g  # [mg] dapagliflozin-3-o-glucuronide amount (arterial blood) [mg]  
Xar_dap = Aar_dap * Mr_dap  # [mg] dapagliflozin amount (arterial blood) [mg]  
Xhv_d3g = Ahv_d3g * Mr_d3g  # [mg] dapagliflozin-3-o-glucuronide amount (hepatic vein) [mg]  
Xhv_dap = Ahv_dap * Mr_dap  # [mg] dapagliflozin amount (hepatic vein) [mg]  
Xpo_d3g = Apo_d3g * Mr_d3g  # [mg] dapagliflozin-3-o-glucuronide amount (portal vein) [mg]  
Xpo_dap = Apo_dap * Mr_dap  # [mg] dapagliflozin amount (portal vein) [mg]  
Xve_d3g = Ave_d3g * Mr_d3g  # [mg] dapagliflozin-3-o-glucuronide amount (venous blood) [mg]  
Xve_dap = Ave_dap * Mr_dap  # [mg] dapagliflozin amount (venous blood) [mg]  
Are_plasma_d3g = Cre_plasma_d3g * Vre_plasma  # [mmol] dapagliflozin-3-o-glucuronide amount (rest) [mmole]  
Are_plasma_dap = Cre_plasma_dap * Vre_plasma  # [mmol] dapagliflozin amount (rest) [mmole]  
Flow_ar_gu_d3g = Qgu * Car_d3g  # [mmol/min] inflow gut dapagliflozin-3-o-glucuronide  
Flow_ar_gu_dap = Qgu * Car_dap  # [mmol/min] inflow gut dapagliflozin  
Flow_ar_ki_d3g = Qki * Car_d3g  # [mmol/min] inflow kidney dapagliflozin-3-o-glucuronide  
Flow_ar_ki_dap = Qki * Car_dap  # [mmol/min] inflow kidney dapagliflozin  
Flow_ar_re_d3g = Qre * Car_d3g  # [mmol/min] inflow rest dapagliflozin-3-o-glucuronide  
Flow_ar_re_dap = Qre * Car_dap  # [mmol/min] inflow rest dapagliflozin  
Flow_gu_po_d3g = Qgu * Cgu_plasma_d3g  # [mmol/min] outflow gut dapagliflozin-3-o-glucuronide  
Flow_gu_po_dap = Qgu * Cgu_plasma_dap  # [mmol/min] outflow gut dapagliflozin  
Flow_hv_ve_d3g = Qh * Chv_d3g  # [mmol/min] outflow hepatic vein dapagliflozin-3-o-glucuronide  
Flow_hv_ve_dap = Qh * Chv_dap  # [mmol/min] outflow hepatic vein dapagliflozin  
Flow_ki_ve_d3g = Qki * Cki_plasma_d3g  # [mmol/min] outflow kidney dapagliflozin-3-o-glucuronide  
Flow_ki_ve_dap = Qki * Cki_plasma_dap  # [mmol/min] outflow kidney dapagliflozin  
Flow_lu_ar_d3g = Qlu * Clu_plasma_d3g  # [mmol/min] outflow lung dapagliflozin-3-o-glucuronide  
Flow_lu_ar_dap = Qlu * Clu_plasma_dap  # [mmol/min] outflow lung dapagliflozin  
Flow_re_ve_d3g = Qre * Cre_plasma_d3g  # [mmol/min] outflow rest dapagliflozin-3-o-glucuronide  
Flow_re_ve_dap = Qre * Cre_plasma_dap  # [mmol/min] outflow rest dapagliflozin  
Flow_ve_lu_d3g = Qlu * Cve_d3g  # [mmol/min] inflow lung dapagliflozin-3-o-glucuronide  
Flow_ve_lu_dap = Qlu * Cve_dap  # [mmol/min] inflow lung dapagliflozin  
Mgu_plasma_d3g = (Agu_plasma_d3g / Vgu_plasma) * Mr_d3g  # [mg/l] dapagliflozin-3-o-glucuronide concentration (gut) [mg/l]  
Mgu_plasma_dap = (Agu_plasma_dap / Vgu_plasma) * Mr_dap  # [mg/l] dapagliflozin concentration (gut) [mg/l]  
Mki_plasma_d3g = (Aki_plasma_d3g / Vki_plasma) * Mr_d3g  # [mg/l] dapagliflozin-3-o-glucuronide concentration (kidney) [mg/l]  
Mki_plasma_dap = (Aki_plasma_dap / Vki_plasma) * Mr_dap  # [mg/l] dapagliflozin concentration (kidney) [mg/l]  
Mli_plasma_d3g = (Ali_plasma_d3g / Vli_plasma) * Mr_d3g  # [mg/l] dapagliflozin-3-o-glucuronide concentration (liver) [mg/l]  
Mli_plasma_dap = (Ali_plasma_dap / Vli_plasma) * Mr_dap  # [mg/l] dapagliflozin concentration (liver) [mg/l]  
Mlu_plasma_d3g = (Alu_plasma_d3g / Vlu_plasma) * Mr_d3g  # [mg/l] dapagliflozin-3-o-glucuronide concentration (lung) [mg/l]  
Mlu_plasma_dap = (Alu_plasma_dap / Vlu_plasma) * Mr_dap  # [mg/l] dapagliflozin concentration (lung) [mg/l]  
Qha = Qh - Qgu  # [l/min] hepatic artery blood flow  
Qpo = Qgu  # [l/min] portal blood flow  
Xgu_plasma_d3g = Agu_plasma_d3g * Mr_d3g  # [mg] dapagliflozin-3-o-glucuronide amount (gut) [mg]  
Xgu_plasma_dap = Agu_plasma_dap * Mr_dap  # [mg] dapagliflozin amount (gut) [mg]  
Xki_plasma_d3g = Aki_plasma_d3g * Mr_d3g  # [mg] dapagliflozin-3-o-glucuronide amount (kidney) [mg]  
Xki_plasma_dap = Aki_plasma_dap * Mr_dap  # [mg] dapagliflozin amount (kidney) [mg]  
Xli_plasma_d3g = Ali_plasma_d3g * Mr_d3g  # [mg] dapagliflozin-3-o-glucuronide amount (liver) [mg]  
Xli_plasma_dap = Ali_plasma_dap * Mr_dap  # [mg] dapagliflozin amount (liver) [mg]  
Xlu_plasma_d3g = Alu_plasma_d3g * Mr_d3g  # [mg] dapagliflozin-3-o-glucuronide amount (lung) [mg]  
Xlu_plasma_dap = Alu_plasma_dap * Mr_dap  # [mg] dapagliflozin amount (lung) [mg]  
Flow_arli_hv_d3g = f_shunts * Qha * Car_d3g  # [mmol/min] flow arterial shunts  
Flow_arli_hv_dap = f_shunts * Qha * Car_dap  # [mmol/min] flow arterial shunts  
Flow_arli_li_d3g = (1 - f_shunts) * Qha * Car_d3g  # [mmol/min] arterial inflow liver dapagliflozin-3-o-glucuronide  
Flow_arli_li_dap = (1 - f_shunts) * Qha * Car_dap  # [mmol/min] arterial inflow liver dapagliflozin  
Flow_li_hv_d3g = (1 - f_shunts) * (Qpo + Qha) * Cli_plasma_d3g  # [mmol/min] outflow liver dapagliflozin-3-o-glucuronide  
Flow_li_hv_dap = (1 - f_shunts) * (Qpo + Qha) * Cli_plasma_dap  # [mmol/min] outflow liver dapagliflozin  
Flow_po_hv_d3g = f_shunts * Qpo * Cpo_d3g  # [mmol/min] flow portal shunts  
Flow_po_hv_dap = f_shunts * Qpo * Cpo_dap  # [mmol/min] flow portal shunts  
Flow_po_li_d3g = (1 - f_shunts) * Qpo * Cpo_d3g  # [mmol/min] outflow po dapagliflozin-3-o-glucuronide  
Flow_po_li_dap = (1 - f_shunts) * Qpo * Cpo_dap  # [mmol/min] outflow po dapagliflozin  
Mre_plasma_d3g = (Are_plasma_d3g / Vre_plasma) * Mr_d3g  # [mg/l] dapagliflozin-3-o-glucuronide concentration (rest) [mg/l]  
Mre_plasma_dap = (Are_plasma_dap / Vre_plasma) * Mr_dap  # [mg/l] dapagliflozin concentration (rest) [mg/l]  
Xre_plasma_d3g = Are_plasma_d3g * Mr_d3g  # [mg] dapagliflozin-3-o-glucuronide amount (rest) [mg]  
Xre_plasma_dap = Are_plasma_dap * Mr_dap  # [mg] dapagliflozin amount (rest) [mg]  

# odes
d Afeces_dap/dt = GU__DAPEXC  # [mmol/min] dapagliflozin (feces)  
d Aurine_d3g/dt = KI__D3GEX  # [mmol/min] dapagliflozin-3-o-glucuronide (urine)  
d Aurine_dap/dt = KI__DAPEX  # [mmol/min] dapagliflozin (urine)  
d Car_d3g/dt = (-Flow_ar_ki_d3g / Var - Flow_arli_li_d3g / Var - Flow_arli_hv_d3g / Var) + Flow_lu_ar_d3g / Var - Flow_ar_gu_d3g / Var - Flow_ar_re_d3g / Var  # [mmol/l/min] dapagliflozin-3-o-glucuronide (arterial blood plasma)  
d Car_dap/dt = (-Flow_ar_ki_dap / Var - Flow_arli_li_dap / Var - Flow_arli_hv_dap / Var) + Flow_lu_ar_dap / Var - Flow_ar_gu_dap / Var - Flow_ar_re_dap / Var  # [mmol/l/min] dapagliflozin (arterial blood plasma)  
d Cgu_dap/dt = (-GU__DAPABS / Vgu - GU__DAPEXC / Vgu) + GU__dissolution_dap / Vgu  # [mmol/l/min] dapagliflozin (gut)  
d Cgu_plasma_d3g/dt = Flow_ar_gu_d3g / Vgu_plasma - Flow_gu_po_d3g / Vgu_plasma  # [mmol/l/min] dapagliflozin-3-o-glucuronide (gut plasma)  
d Cgu_plasma_dap/dt = (Flow_ar_gu_dap / Vgu_plasma - Flow_gu_po_dap / Vgu_plasma) + GU__DAPABS / Vgu_plasma  # [mmol/l/min] dapagliflozin (gut plasma)  
d Chv_d3g/dt = Flow_arli_hv_d3g / Vhv + Flow_po_hv_d3g / Vhv + Flow_li_hv_d3g / Vhv - Flow_hv_ve_d3g / Vhv  # [mmol/l/min] dapagliflozin-3-o-glucuronide (hepatic vein plasma)  
d Chv_dap/dt = Flow_arli_hv_dap / Vhv + Flow_po_hv_dap / Vhv + Flow_li_hv_dap / Vhv - Flow_hv_ve_dap / Vhv  # [mmol/l/min] dapagliflozin (hepatic vein plasma)  
d Cki_plasma_d3g/dt = Flow_ar_ki_d3g / Vki_plasma - Flow_ki_ve_d3g / Vki_plasma - KI__D3GEX / Vki_plasma  # [mmol/l/min] dapagliflozin-3-o-glucuronide (kidney plasma)  
d Cki_plasma_dap/dt = Flow_ar_ki_dap / Vki_plasma - Flow_ki_ve_dap / Vki_plasma - KI__DAPEX / Vki_plasma  # [mmol/l/min] dapagliflozin (kidney plasma)  
d Cli_plasma_d3g/dt = (Flow_arli_li_d3g / Vli_plasma + Flow_po_li_d3g / Vli_plasma - Flow_li_hv_d3g / Vli_plasma) + LI__D3GEX / Vli_plasma  # [mmol/l/min] dapagliflozin-3-o-glucuronide (liver plasma)  
d Cli_plasma_dap/dt = Flow_arli_li_dap / Vli_plasma + Flow_po_li_dap / Vli_plasma - Flow_li_hv_dap / Vli_plasma - LI__DAPIM / Vli_plasma  # [mmol/l/min] dapagliflozin (liver plasma)  
d Clu_dap/dt = transport_lu_dap / Vlu_tissue  # [mmol/l/min] dapagliflozin (lung)  
d Clu_plasma_d3g/dt = Flow_ve_lu_d3g / Vlu_plasma - Flow_lu_ar_d3g / Vlu_plasma  # [mmol/l/min] dapagliflozin-3-o-glucuronide (lung plasma)  
d Clu_plasma_dap/dt = -transport_lu_dap / Vlu_plasma + Flow_ve_lu_dap / Vlu_plasma - Flow_lu_ar_dap / Vlu_plasma  # [mmol/l/min] dapagliflozin (lung plasma)  
d Cpo_d3g/dt = (-Flow_po_li_d3g / Vpo - Flow_po_hv_d3g / Vpo) + Flow_gu_po_d3g / Vpo  # [mmol/l/min] dapagliflozin-3-o-glucuronide (portal vein plasma)  
d Cpo_dap/dt = (-Flow_po_li_dap / Vpo - Flow_po_hv_dap / Vpo) + Flow_gu_po_dap / Vpo  # [mmol/l/min] dapagliflozin (portal vein plasma)  
d Cre_dap/dt = transport_re_dap / Vre_tissue  # [mmol/l/min] dapagliflozin (rest)  
d Cre_plasma_d3g/dt = Flow_ar_re_d3g / Vre_plasma - Flow_re_ve_d3g / Vre_plasma  # [mmol/l/min] dapagliflozin-3-o-glucuronide (rest plasma)  
d Cre_plasma_dap/dt = -transport_re_dap / Vre_plasma + Flow_ar_re_dap / Vre_plasma - Flow_re_ve_dap / Vre_plasma  # [mmol/l/min] dapagliflozin (rest plasma)  
d Cve_d3g/dt = (Flow_ki_ve_d3g / Vve + Flow_hv_ve_d3g / Vve - Flow_ve_lu_d3g / Vve) + Flow_re_ve_d3g / Vve  # [mmol/l/min] dapagliflozin-3-o-glucuronide (venous blood plasma)  
d Cve_dap/dt = (iv_dap / Vve + Flow_ki_ve_dap / Vve + Flow_hv_ve_dap / Vve - Flow_ve_lu_dap / Vve) + Flow_re_ve_dap / Vve  # [mmol/l/min] dapagliflozin (venous blood plasma)  
d GU__dap_stomach/dt = 0  # [mmol/min] dapagliflozin (stomach)  
d IVDOSE_dap/dt = -iv_dap * Mr_dap + Ri_dap  # [mg/min] IV bolus dose dap [mg]  
d KI__glc_ext/dt = 0  # [mmol/l/min] glucose (plasma)  
d KI__glc_urine/dt = KI__GLCEX  # [mmol/min] glucose (urine)  
d LI__d3g/dt = LI__DAP2D3G / Vli_tissue - LI__D3GEX / Vli_tissue  # [mmol/l/min] dapagliflozin-3-o-glucuronide (liver)  
d LI__dap/dt = LI__DAPIM / Vli_tissue - LI__DAP2D3G / Vli_tissue  # [mmol/l/min] dapagliflozin (liver)  
d PODOSE_dap/dt = -GU__dissolution_dap * GU__Mr_dap  # [mg/min] oral dose dap [mg]  
d cum_dose_dap/dt = Ri_dap  # [mg/min] Cumulative dose due to infusion dap  
```