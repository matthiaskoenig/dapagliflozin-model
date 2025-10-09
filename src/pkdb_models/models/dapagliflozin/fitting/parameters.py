"""FitParameters for dapagliflozin fitting."""
import copy

from sbmlsim.fit import FitParameter


parameters_dap_iv = [
    # tissue distribution
    FitParameter(
        pid="ftissue_dap",
        lower_bound=0.01,
        start_value=1,
        upper_bound=100,
        unit="l/min",
    ),
    FitParameter(
        pid="Kp_dap",
        lower_bound=1,
        start_value=10,
        upper_bound=50,
        unit="dimensionless",
    ),

    # hepatic transport
    # FitParameter(
    #     pid="LI__DAPIM_Vmax",
    #     lower_bound=1E-2,
    #     start_value=1.0,
    #     upper_bound=1E4,
    #     unit="mmole/min/l",
    # ),
    # FitParameter(
    #     pid="LI__D3GEX_Vmax",
    #     lower_bound=1E-3,
    #     start_value=1.0,
    #     upper_bound=1E3,
    #     unit="mmole/min/l",
    # ),

    # hepatic metabolism
    FitParameter(
        pid="DAP2D3G_Vmax",
        lower_bound=1E-3,
        start_value=1.0,
        upper_bound=100,
        unit="mmole/min/l",
    ),

    # kidney transport
    # FitParameter(
    #     pid="KI__DAPIM_Vmax",
    #     lower_bound=1E-2,
    #     start_value=1.0,
    #     upper_bound=1E4,
    #     unit="mmole/min/l",
    # ),
    # FitParameter(
    #     pid="KI__D3GIM_Vmax",
    #     lower_bound=1E-3,
    #     start_value=1.0,
    #     upper_bound=1E3,
    #     unit="mmole/min/l",
    # ),

    # kidney metabolism (activity relative to liver)
    FitParameter(
        pid="KI__f_DAP2D3G",
        lower_bound=0.1,
        start_value=1.0,
        upper_bound=10,
        unit="dimensionless",
    ),
    # kidney removal
    FitParameter(
        pid="KI__DAPEX_k",
        lower_bound=1E-4,
        start_value=0.1,
        upper_bound=10,
        unit="1/min",
    ),
    FitParameter(
        pid="KI__D3GEX_k",
        start_value=1.0,
        lower_bound=1E-1,
        upper_bound=10,
        unit="1/min",
    ),
]

parameters_dap_po = [
    # dissolution rate
    FitParameter(
        pid="GU__Ka_dis_dap",
        lower_bound=1E-3,
        start_value=2,
        upper_bound=100,
        unit="1/hr",
    ),
    # absorption rate
    FitParameter(
        pid="GU__DAPABS_k",
        lower_bound=1E-5,
        start_value=0.02,
        upper_bound=10,
        unit="1/min",
    ),
    # fraction absorbed
    # FitParameter(
    #     pid="GU__F_dap_abs",
    #     lower_bound=0.80,
    #     start_value=0.84,
    #     upper_bound=0.85,
    #     unit="dimensionless",
    # ),
]


# --- pharmacodynamics ---
parameters_pd = [
    FitParameter(
        pid="KI__RTG_E50",
        lower_bound=1E-6,
        start_value=1E-3,
        upper_bound=1E-1,
        unit="mM",
    ),
    FitParameter(
        pid="KI__RTG_base",
        lower_bound=8,
        start_value=12.5,
        upper_bound=14,
        unit="mM",
    ),
    # FitParameter(
    #     pid="KI__RTG_delta",
    #     lower_bound=7,
    #     start_value=9,
    #     upper_bound=11,
    #     unit="mM",
    # ),
    FitParameter(
        pid="KI__RTG_gamma",
        lower_bound=1,
        start_value=1.01,
        upper_bound=5,
        unit="dimensionless",
    ),
    FitParameter(
        pid="KI__RTG_max_inhibition",
        lower_bound=0.2,
        start_value=0.75,
        upper_bound=1.0,
        unit="dimensionless",
    ),
    FitParameter(
        pid="KI__RTG_m_fpg",
        lower_bound=0.2,
        start_value=1,
        upper_bound=3,
        unit="dimensionless",
    ),
]

parameters_pk = parameters_dap_iv + parameters_dap_po


