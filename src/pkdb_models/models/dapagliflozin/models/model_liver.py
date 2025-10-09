"""Liver model for SGLT2 inhibitors."""

import numpy as np
import pandas as pd
from sbmlutils.factory import *
from sbmlutils.metadata import *

from pkdb_models.models.dapagliflozin.models import annotations
from pkdb_models.models.dapagliflozin.models import templates


class U(templates.U):
    """UnitDefinitions"""

    pass


mid = "dapagliflozin_liver"
version = 1

_m = Model(
    sid=mid,
    name="Model for hepatic dapagliflozin metabolism.",
    notes=f"""
    Model for dapagliflozin metabolism.
    
    - UGT1A9: 90% [Callegari2021]; 79.9% [Jo2021]
    - UGT2B4/2B7: 10% [Callegari2021]; 8.9% [Jo2021]
    - CYP3A4: 9.8% [Jo2021]
    - in hepatic impairment: increase in AUC & Cmax of dap and d3g [Kasichayanula2011]
    => only conversion dap -> d3g [UGT1A9], no d2g in first version; only UGT1A9 included
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
        "Vli",
        value=1.5,
        unit=U.liter,
        name="liver",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["li"],
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
    )

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
        "dap",
        name="dapagliflozin (liver)",
        initialConcentration=0.0,
        compartment="Vli",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["dap"],
    ),
    Species(
        "d3g",
        name="dapagliflozin-3-o-glucuronide (liver)",
        initialConcentration=0.0,
        compartment="Vli",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["d3g"],
    ),
]
_m.parameters = [
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
]

_m.reactions = [
    Reaction(
        sid="DAPIM",
        name="dapagliflozin import (DAPIM, OAT3)",
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
                name="Vmax dapagliflozin import",
                sboTerm=SBO.MICHAELIS_CONSTANT,
                notes="""Possibly catalyzed by OAT3 with Ki=33µM. [FDA]"""
            ),
        ],
        formula=(
            "DAPIM_Vmax/DAPIM_Km_dap * Vli * (dap_ext - dap)/(1 dimensionless + dap_ext/DAPIM_Km_dap + dap/DAPIM_Km_dap)"
        ),
        notes="""Possibly catalyzed by OAT3 with Ki=33µM. [FDA]"""
    ),
    Reaction(
        sid="DAP2D3G",
        name="dapagliflozin conversion (DAP2D3G) UGT1A9",
        equation="dap -> d3g",
        compartment="Vli",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        pars=[],
        formula=(
            "f_ugt1a9 * DAP2D3G_Vmax * Vli * dap/(dap + DAP2D3G_Km_dap)"
        ),
    ),
    Reaction(
        sid="D3GEX",
        name="dapagliflozin-3-o-glucuronide export (D3GEX)",
        equation="d3g <-> d3g_ext",
        compartment="Vmem",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "D3GEX_Vmax",
                10.0,
                U.mmole_per_min_l,
                name="Vmax dapagliflozinat export",
                sboTerm=SBO.MAXIMAL_VELOCITY,
            ),
            Parameter(
                "D3GEX_Km_d3g",
                0.115,
                U.mM,
                name="Vmax dapagliflozinat export",
                sboTerm=SBO.MICHAELIS_CONSTANT,
                notes="""Substrate of hOAT3 and its Km value for hOAT3-mediated uptake was 115µM [FDA]."""
            ),
        ],
        formula=(
            "D3GEX_Vmax/D3GEX_Km_d3g * Vli * (d3g - d3g_ext)/(1 dimensionless + d3g/D3GEX_Km_d3g + d3g_ext/D3GEX_Km_d3g)"
        ),
        notes="""Possibly catalyzed by OAT3 with Km=115µM. [FDA]"""
    ),
]


model_liver = _m

def dapagliflozin_layout(dx=200, dy=200) -> pd.DataFrame:
    """Layout definition."""

    delta_y = 0.5 * dy
    delta_x = 0.7 * dx

    positions = [
        # sid, x, y

        ["dap_ext", 0, 0],
        ["d3g_ext", 2 * delta_x, 0],

        ["DAPIM", 0, delta_y],
        ["D3GEX", 2 * delta_x, delta_y],

        ["dap", 0, 2*delta_y],
        ["d3g", 2 * delta_x, 2*delta_y],

        ["DAP2D3G", 1*delta_x, 2.5* delta_y],
    ]

    df = pd.DataFrame(positions, columns=["id", "x", "y"])
    df.set_index("id", inplace=True)
    return df


def dapagliflozin_annotations(dx=200, dy=200) -> list:
    COLOR_BLOOD = "#FF796C"
    COLOR_CELL = "#FFFFFF"
    delta_y = 0.5 * dy

    kwargs = {
        "type": cyviz.AnnotationShapeType.ROUND_RECTANGLE,
        "opacity": 20,
        "border_color": "#000000",
        "border_thickness": 2,
    }
    xpos = -0.5 * dx
    width = 2.5 * dx

    annotations = [
        cyviz.AnnotationShape(
            x_pos=xpos, y_pos=-0.5*delta_y, width=width, height=1.5 * delta_y,
            fill_color=COLOR_BLOOD, **kwargs
        ),
        cyviz.AnnotationShape(
            x_pos=xpos, y_pos=delta_y, width=width, height=2* delta_y,
            fill_color=COLOR_CELL, **kwargs
        ),
    ]
    return annotations


if __name__ == "__main__":
    from sbmlutils import cytoscape as cyviz
    from pkdb_models.models.dapagliflozin import MODEL_BASE_PATH
    results: FactoryResult = create_model(
        model=model_liver,
        filepath=MODEL_BASE_PATH / f"{model_liver.sid}.xml",
        sbml_level=3, sbml_version=2,
    )
    cyviz.visualize_sbml(sbml_path=results.sbml_path, delete_session=True)
    cyviz.apply_layout(layout=dapagliflozin_layout())
    cyviz.add_annotations(annotations=dapagliflozin_annotations())
