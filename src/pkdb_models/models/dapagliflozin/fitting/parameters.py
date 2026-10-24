"""Parameters to optimize in the dapagliflozin fits, with their bounds.

The start values are the values of the model, chosen by hand against the
literature (Boulton2013, Obermeier2010, Kasichayanula2013, Kasichayanula2011b,
Komoroski2009, Gould2013, Yang2013): IV clearance 216 ml/min (207), Vss 92 l
(118), half-life 11.2 hr (12.2), 10 mg Cmax 155 ng/ml (136-223), tmax 1.0 hr
(1.0-1.25), AUCinf 647 ng*h/ml (420-940), urine 1.3 % dapagliflozin (1.2-1.6)
and 63 % D3G (61-63), UGE 0-24 hr of healthy subjects 43 g with 10 mg (36-47).
Not fitted: the slow pool of the rest tissue (`Kpdeep_dap` 0.5, `kdeep_dap`
0.002 1/min), which the data of the fits do not identify (the profiles of
PK_20261019_231348__82b65 are flat), set to the terminal half-life of 12-13 hr
and the accumulation of 10 mg once daily [Boulton2013, Komoroski2009]; the
fraction absorbed (`GU__F_dap_abs` 0.84, the feces of
Obermeier2010), the fraction unbound (`KI__fu_dap`), the renal fraction of the
metabolism (`f_met_renal`, which the data of healthy subjects cannot tell from
the hepatic one; it is set from the AUC ratios of renal and hepatic impairment,
Kasichayanula2013 and Kasichayanula2011, which are validation data), the
fraction of the metabolism to D3G (`f_d3g`, the mass balance of
Obermeier2010), the factor of the interactions (`f_ugt1a9`), the RTG without
inhibition of healthy subjects (`KI__RTG_base` 11 mM, 10-11 mM in the
literature; PD_20261020_012618__f9531 does not identify it above 10 mM), the
splay of the glucose excretion
(`KI__RTG_splay`) and the difference of the mean to the fasting plasma glucose
(`KI__PG_delta`, set by the studies).
"""

from sbmlsim.fit import FitParameter

#: the absorption and the distribution of dapagliflozin: tmax 1.0-1.25 hr
#: fasted [Boulton2013, Kasichayanula2013], a high fat meal lowers Cmax by
#: 31 % and delays tmax by about 1 hr, the AUC is unchanged
#: [Kasichayanula2011b]; Vss 118 l and t1/2 12.2 hr after an intravenous
#: microdose [Boulton2013], t1/2 12-14 hr after oral doses: a fast pool of the
#: tissues (`ftissue_dap`, `Kp_dap`) and a slow pool of the rest tissue
#: (`Kpdeep_dap`, `kdeep_dap`, not fitted); a single pool fits the fast decline
#: of the first 12 hr only with a terminal half-life of 7 hr
parameters_pk_absorption_distribution: list[FitParameter] = [
    FitParameter(
        pid="ftissue_dap",
        start_value=0.3,
        lower_bound=1e-3,
        upper_bound=10,
        unit="l/min",
    ),
    FitParameter(
        pid="Kp_dap",
        start_value=0.8,
        lower_bound=0.2,
        upper_bound=10,
        unit="dimensionless",
    ),
    FitParameter(
        pid="GU__Ka_dis_dap",
        start_value=1.0,
        lower_bound=0.05,
        upper_bound=20,
        unit="1/hr",
    ),
    FitParameter(
        pid="GU__DAPABS_k",
        start_value=0.035,
        lower_bound=1e-3,
        upper_bound=1,
        unit="1/min",
    ),
    FitParameter(
        pid="GU__FED_f_absorption",
        start_value=0.4,
        lower_bound=0.05,
        upper_bound=1.0,
        unit="dimensionless",
    ),
]
#: the metabolism of dapagliflozin (hepatic and, `f_met_renal` of it, renal)
#: and the renal excretion relative to the GFR: clearance 207 ml/min
#: [Boulton2013], renal clearance of dapagliflozin 3.5-5.6 ml/min (below the
#: filtration of the unbound fraction, 9 ml/min), urine 1.2-1.6 % of the dose
#: [Obermeier2010, Kasichayanula2013], renal clearance of D3G 151 ml/min,
#: molar AUC ratio D3G/dapagliflozin 1.12 [Kasichayanula2013]
parameters_pk_elimination: list[FitParameter] = [
    FitParameter(
        pid="DAPMET_k",
        start_value=0.0565,
        lower_bound=1e-3,
        upper_bound=1,
        unit="1/min",
    ),
    FitParameter(
        pid="KI__DAPEX_f",
        start_value=0.42,
        lower_bound=0.05,
        upper_bound=5,
        unit="dimensionless",
    ),
    FitParameter(
        pid="KI__D3GEX_f",
        start_value=1.5,
        lower_bound=0.1,
        upper_bound=10,
        unit="dimensionless",
    ),
]
#: the pharmacodynamics, the renal threshold of glucose (RTG): the urinary
#: glucose excretion of healthy subjects rises from 0.6 g/24 hr with 0.3 mg to
#: 36-47 g with 10 mg and 69-75 g with 250-500 mg [Gould2013, Komoroski2009,
#: Yang2013], of type 2 diabetes 38-77 g with 2.5-20 mg [Kasichayanula2011a,
#: FDAMB102003]; the IC50 of SGLT2 is about 1.1 nM (unbound), the bounds of the
#: EC50 are 0.04-40 ng/ml (total)
parameters_pd: list[FitParameter] = [
    FitParameter(
        pid="KI__RTG_E50",
        start_value=4e-6,
        lower_bound=1e-7,
        upper_bound=1e-4,
        unit="mM",
    ),
    FitParameter(
        pid="KI__RTG_gamma",
        start_value=1.5,
        # at least 1: an exponent below 1 makes the integration stiff at zero
        # dapagliflozin
        lower_bound=1.0,
        upper_bound=4.0,
        unit="dimensionless",
    ),
    FitParameter(
        pid="KI__RTG_m_fpg",
        start_value=1.5,
        # positive for the search in log10 space
        lower_bound=0.1,
        upper_bound=5,
        unit="dimensionless",
    ),
    FitParameter(
        pid="KI__RTG_max_inhibition",
        start_value=0.70,
        lower_bound=0.3,
        upper_bound=0.95,
        unit="dimensionless",
    ),
]

#: the pharmacokinetics
parameters_pk: list[FitParameter] = (
    parameters_pk_absorption_distribution + parameters_pk_elimination
)
#: the pharmacokinetics and the pharmacodynamics
parameters_pk_pd: list[FitParameter] = parameters_pk + parameters_pd
