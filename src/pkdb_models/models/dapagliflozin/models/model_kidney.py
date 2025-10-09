"""Kidney model for the SGLT2 inhibitors.

"""

import numpy as np
import pandas as pd
from sbmlutils.factory import *
from sbmlutils.metadata import *

from pkdb_models.models.dapagliflozin.models import annotations
from pkdb_models.models.dapagliflozin.models import templates


class U(templates.U):
    """UnitDefinitions"""

    mg_per_g = UnitDefinition("mg_per_g", "mg/g")
    ml_per_l = UnitDefinition("ml_per_l", "ml/l")
    ml_per_min = UnitDefinition("ml_per_min", "ml/min")


mid = "dapagliflozin_kidney"
version = 3

_m = Model(
    sid=mid,
    name="Model for renal dapagliflozin and dapagliflozin-3-o-glucuronide excretion.",
    notes=f"""
    Model for renal dapagliflozin and dapagliflozin-3-o-glucuronide excretion.
    
    - ~82 % in urine (2% dap, 80% metabolites)[Callegari2021]; (2% dap, 72% metabolites)[Jo2021]
    - ~18 % in feces (16% dap, 2% metabolites)[Callegari2021]; (1.7% metabolites)[Jo2021]
    - kidney impairment has strong effect on dap and d3g plasma concentration (increase)
    
    => mainly renal excretion of d3g, small excretion of dag (~1-2%)
    
    **version** {version}
    
    ## Changelog
    
    **version 1**
    
    - initial model
        
    """ + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
)

_m.compartments = [
    Compartment(
        "Vext",
        value=1.5,
        unit=U.liter,
        name="plasma",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["plasma"],
        port=True
    ),
    Compartment(
        "Vki",
        value=0.3,  # 0.4 % of bodyweight
        unit=U.liter,
        name="kidney",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["ki"],
        port=True
    ),
    Compartment(
        "Vmem",
        value=np.nan,
        unit=U.m2,
        name="plasma membrane",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["plasma membrane"],
        spatialDimensions=2,
    ),
    Compartment(
        "Vurine",
        1.0,
        name="urine",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        port=True,
        annotations=annotations.compartments["urine"],
    ),

]

_m.species = [
    Species(
        "dap_ext",
        name="dapagliflozin (plasma)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["dap"],
        port=True
    ),
    Species(
        "dap",
        name="dapagliflozin (kidney)",
        initialConcentration=0.0,
        compartment="Vki",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["dap"],
    ),
    Species(
        "dap_urine",
        name="dapagliflozin (urine)",
        initialConcentration=0.0,
        compartment="Vurine",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["dap"],
        port=True
    ),
    Species(
        "d3g_ext",
        name="dapagliflozin-3-o-glucuronide (plasma)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["d3g"],
        port=True
    ),
    Species(
        "d3g",
        name="dapagliflozin-3-o-glucuronide (kidney)",
        initialConcentration=0.0,
        compartment="Vki",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["d3g"],
    ),
    Species(
        "d3g_urine",
        name="dapagliflozin-3-o-glucuronide (urine)",
        initialConcentration=0.0,
        compartment="Vurine",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["d3g"],
        port=True
    ),

]

_m.parameters.extend([
    Parameter(
        "f_renal_function",
        name="parameter for renal function",
        value=1.0,
        unit=U.dimensionless,
        sboTerm=SBO.KINETIC_CONSTANT,
        notes="""scaling factor for renal function. 1.0: normal renal function; 
        <1.0: reduced renal function
        """
    ),
    Parameter(
        "DAP2D3G_Vmax",
        0.04,
        U.mmole_per_min_l,
        name="Vmax dapagliflozin conversion",
        sboTerm=SBO.MAXIMAL_VELOCITY,
        port=True,
    ),
    Parameter(
        "DAP2D3G_Km_dap",
        0.479,
        U.mM,
        name="Km dapagliflozin UGT1A9",
        sboTerm=SBO.MICHAELIS_CONSTANT,
        port=True,
        notes="""
        Km apparent in kidney microsomes for dapagliflozin; 479 µM [Kasichaynula2013]
        Ki for UGT1A9 with other substrates 12-15 µM [Obermeier2010a]."""
    ),
    Parameter(
        "f_ugt1a9",
        1,
        U.dimensionless,
        name="scaling factor UGT1A9 activity",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        port=True,
        notes="""Scaling factor to vary UGT1A9 activity.
    1.0: unchanged activity; < 1.0 decreased activity; >1.0 increased activity.
    """
    )
])

_m.reactions = [
    Reaction(
        sid="DAPIM",
        name="DAPIM",  # dapagliflozin import
        equation="dap_ext <-> dap",
        compartment="Vmem",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "DAPIM_Vmax",
                10.0,
                U.mmole_per_min_l,
                name="Vmax dapagliflozin import",
                sboTerm=SBO.MAXIMAL_VELOCITY,
            ),
            Parameter(
                "DAPIM_Km_dap",
                0.033,
                U.mM,
                name="Km dapagliflozin import",
                sboTerm=SBO.MICHAELIS_CONSTANT,
                notes="""Possibly catalyzed by OAT3 with Ki=33µM. [FDA]"""
            ),
        ],
        formula=(
            "f_renal_function * DAPIM_Vmax/DAPIM_Km_dap * Vki * (dap_ext - dap)/(1 dimensionless + dap_ext/DAPIM_Km_dap + dap/DAPIM_Km_dap)"
        )
    ),
    Reaction(
        sid="D3GIM",
        name="D3GIM",  # dapagliflozin-3-O-glucuronide import
        equation="d3g_ext <-> d3g",
        compartment="Vmem",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "D3GIM_Vmax",
                10.0,
                U.mmole_per_min_l,
                name="Vmax dapagliflozin-3-O-glucuronide import",
                sboTerm=SBO.MAXIMAL_VELOCITY,
            ),
            Parameter(
                "D3GIM_Km_d3g",
                0.033,
                U.mM,
                name="Km dapagliflozin-3-O-glucuronide import",
                sboTerm=SBO.MICHAELIS_CONSTANT,
                notes="""Using Km identical to dap"""
            ),
        ],
        formula=(
            "f_renal_function * D3GIM_Vmax/D3GIM_Km_d3g * Vki * (d3g_ext - d3g)/(1 dimensionless + d3g_ext/D3GIM_Km_d3g + d3g/D3GIM_Km_d3g)"
        ),
    ),
    Reaction(
        sid="DAP2D3G",
        name="UGT1A9 (dap -> d3g)",
        equation="dap -> d3g",
        compartment="Vki",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        pars=[
            Parameter(
                "f_DAP2D3G",
                9.999990451401281,
                U.dimensionless,
                name="scaling factor relative to liver activity",
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                notes="""Scaling factor relative to liver activity.
                1.0: equal activity to liver; < 1.0 decreased activity; >1.0 increased activity.
                """
            )
        ],
        formula=(
            "f_renal_function * f_ugt1a9 * f_DAP2D3G * DAP2D3G_Vmax * Vki * dap/(dap + DAP2D3G_Km_dap)"
        ),
    ),
    Reaction(
        sid="D3GEX",
        name="D3GEX",  # dapagliflozin-3-o-glucuronide excretion
        equation="d3g_ext -> d3g_urine",
        compartment="Vki",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "D3GEX_k",
                0.45035618074418376,
                U.per_min,
                name="rate urinary excretion of dapagliflozin-3-o-glucuronide",
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ],
        formula=(
            "f_renal_function * D3GEX_k * Vki * d3g_ext"
        )
    ),
    Reaction(
        sid="DAPEX",
        name="DAPEX",  # dapagliflozin excretion
        equation="dap_ext -> dap_urine",
        compartment="Vki",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "DAPEX_k",
                0.01815179124844871,
                U.per_min,
                name="rate urinary excretion of dapagliflozin",
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ],
        formula=(
            "f_renal_function * DAPEX_k * Vki * dap_ext"
        )
    ),
]

# ---------------------------------------------------------------------------------------------------------------------
# Pharmacodynamics
# ---------------------------------------------------------------------------------------------------------------------
_m.species.extend([
    Species(
        "glc_ext",
        name="glucose (plasma)",
        initialConcentration=5.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        boundaryCondition=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["glc"],
    ),
    Species(
        "glc_urine",
        name="glucose (urine)",
        initialConcentration=0.0,
        compartment="Vurine",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["glc"],
    ),
])

_m.parameters.extend([
    Parameter(
        "Mr_glc",
        180,
        U.g_per_mole,
        name=f"Molecular weight glc [g/mole]",
        sboTerm=SBO.MOLECULAR_MASS,
    ),
    Parameter(
        "cf_mg_per_g",
        1000,
        U.mg_per_g,
        name=f"Conversion factor mg per g",
    ),
    Parameter(
        "cf_ml_per_l",
        1000,
        U.ml_per_l,
        name=f"Conversion factor ml per l",
    ),
    Parameter(
        "RTG_E50",
        6.493779238072141e-06,  # ~1 µmole/l [mmole/l]
        U.mM,
        name="EC50 reduction in RTG",
        sboTerm=SBO.KINETIC_CONSTANT,
        notes="""
        ~10 ng/ml graphical estimate for low doses [Komoroski2009]
        Half-maximal pharmacodynamic effect of dapagliflozin.
        """
    ),
    Parameter(
        "RTG_gamma",
        1,  # [1 - 4] for optimization
        U.dimensionless,
        name="hill coefficient reduction in RTG",
        sboTerm=SBO.KINETIC_CONSTANT,
        notes="""gamma effect of dapagliflozin pharmacodynamics."""
    ),
    Parameter(
        "RTG_base",
        8.000008431987332,  # [10.5 - 14] for optimization
        U.mM,
        name=f"Baseline RTG value",
        notes="typical RTG value without SGLT2 inhibitors in healthy subjects"
    ),
    # Parameter(
    #     "RTG_delta",
    #     7.007231048860945,  # [7 - 10] for optimization
    #     U.mM,
    #     name=f"RTG value",
    #     notes="ΔRTG, maximum reduction in RTG due to SGLT2 inhibition"
    # ),
    Parameter(
        "RTG_m_fpg",
        1.2533715089157764,  # [0.2 - 1] for optimization
        U.dimensionless,
        name=f"FPG effect on RTG",
        notes="""Effect of the FPG on the change in RTG_base."""
    ),
    Parameter(
        "RTG_max_inhibition",
        0.7067308400163536,  # [0 - 1] for optimization
        U.dimensionless,
        name=f"RTG maximum inhibition",
        notes="maximum inhibition of RTG via SGLT2"
    ),
    Parameter(
        "fpg_healthy",
        5,
        U.mM,
        name=f"fasting plasma glucose (healthy)",
    ),

    Parameter(
        "GFR_healthy",
        100,
        U.ml_per_min,
        name=f"Glomerular filtration rate (healthy)",
    ),
])

_m.rules.extend([
    AssignmentRule(
        "RTG_fpg",
        "RTG_base + RTG_m_fpg * (glc_ext - fpg_healthy)",
        U.mM,
        name=f"RTG value (FPG)",
        notes="""
    FPG dependent base RTG value. SGLT2 is induced in T2DM [Rahmoune2005] 
    """
    ),
    AssignmentRule(
        "RTG_delta",
        "RTG_fpg * RTG_max_inhibition",  # [7 - 10]
        U.mM,
        name=f"RTG value",
        notes="""
    ΔRTG, maximum reduction in RTG due to SGLT2 inhibition.
    Normally around 7 - 10 mM. 
    """
    ),
    AssignmentRule(
        variable="RTG",
        value="RTG_fpg - RTG_delta * power(dap_ext, RTG_gamma)/ (power(RTG_E50, RTG_gamma) + power(dap, RTG_gamma))",
        unit=U.mM,
        name="renal threshold glucose (RTG)",
        notes="""
        Renal threshold glucose.

        The renal threshold for glucose (RTG) is the plasma glucose concentration at which tubular reabsorption of 
        glucose begins to saturate; glucose is excreted into the urine in direct proportion 
        to the glucose concentration above this threshold.
        12.3 - 12.7 mM RTG (placebo) [Devineni2012]
        """
    ),
    AssignmentRule(
        variable="GFR",
        value="f_renal_function * GFR_healthy",
        unit=U.ml_per_min,
        name="glomerular filtration rate",
    ),
])


# Glucose excretion (UGE)
_m.reactions.extend([
    Reaction(
        sid="GLCEX",
        name="glucose excretion (GLCEX)",
        equation="glc_ext -> glc_urine [dap]",
        compartment="Vki",
        sboTerm=SBO.TRANSPORT_REACTION,
        #  piecewise      | x1, y1, [x2, y2,] [...] [z] | A piecewise function: if (y1), x1.  Otherwise, if (y2), x2, etc.  Otherwise, z.
        formula=(
            # [ml/min]/[ml/l] * [mmole/l] = [mmole/min]
            "piecewise(GFR/cf_ml_per_l * (glc_ext - RTG), glc_ext > RTG, 0 mmole_per_min)"   # FIXME: no scaling with liver volume Vki? (GFR should scale with kidney volume)
        )
    ),
])

_m.rules.extend([
    AssignmentRule(
        "UGE", "glc_urine * Mr_glc/cf_mg_per_g", unit=U.gram,
        name="urinary glucose excretion (UGE)",
        notes="""
        Urinary glucose excretion is calculated from cumulative amount of glucose in urine.
        """
    )
])

model_kidney = _m


def dapagliflozin_layout(dx=200, dy=200) -> pd.DataFrame:
    """Layout definition."""

    delta_y = 0.5 * dy
    delta_x = 0.7 * dx

    positions = [
        # sid, x, y
        ["glc_ext", 0, 0],
        ["dap_ext", delta_x, 0],
        ["d3g_ext", 3 * delta_x, 0],


        ["DAPIM", delta_x, delta_y],
        ["D3GIM", 3 * delta_x, delta_y],

        ["dap", delta_x, 2*delta_y],
        ["d3g", 3 * delta_x, 2*delta_y],

        ["GLCEX", 0, 3*delta_y],
        ["DAPEX", delta_x, 3* delta_y],
        ["DAP2D3G", 2*delta_x, 2.5* delta_y],
        ["D3GEX", 3*delta_x, 3* delta_y],

        ["glc_urine", 0, 4*delta_y],
        ["dap_urine", delta_x, 4*delta_y],
        ["d3g_urine", 3 * delta_x, 4*delta_y],
    ]

    df = pd.DataFrame(positions, columns=["id", "x", "y"])
    df.set_index("id", inplace=True)
    return df


def dapagliflozin_annotations(dx=200, dy=200) -> list:
    COLOR_BLOOD = "#FF796C"
    COLOR_CELL = "#FFFFFF"
    COLOR_URINE = "#FF7F0E"
    delta_y = 0.5 * dy

    kwargs = {
        "type": cyviz.AnnotationShapeType.ROUND_RECTANGLE,
        "opacity": 20,
        "border_color": "#000000",
        "border_thickness": 2,
    }
    xpos = -0.5 * dx
    width = 3.5 * dx

    annotations = [
        cyviz.AnnotationShape(
            x_pos=xpos, y_pos=-0.5*delta_y, width=width, height=1.5 * delta_y,
            fill_color=COLOR_BLOOD, **kwargs
        ),
        cyviz.AnnotationShape(
            x_pos=xpos, y_pos=delta_y, width=width, height=2* delta_y,
            fill_color=COLOR_CELL, **kwargs
        ),
        cyviz.AnnotationShape(
            x_pos=xpos, y_pos=3 * delta_y, width=width, height=1.5 * delta_y,
            fill_color=COLOR_URINE, **kwargs
        )
    ]
    return annotations


if __name__ == "__main__":
    from pkdb_models.models.dapagliflozin import MODEL_BASE_PATH

    from sbmlutils import cytoscape as cyviz
    results: FactoryResult = create_model(
        model=model_kidney,
        filepath=MODEL_BASE_PATH / f"{model_kidney.sid}.xml",
        sbml_level=3, sbml_version=2,
    )

    cyviz.visualize_sbml(sbml_path=results.sbml_path, delete_session=True)
    cyviz.apply_layout(layout=dapagliflozin_layout())
    cyviz.add_annotations(annotations=dapagliflozin_annotations())

