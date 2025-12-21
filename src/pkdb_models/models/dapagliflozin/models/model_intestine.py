"""Dapagliflozin intestine model."""
import numpy as np
import pandas as pd

from sbmlutils.factory import *
from sbmlutils.metadata import *

from pkdb_models.models.dapagliflozin.models import annotations
from pkdb_models.models.dapagliflozin.models import templates


class U(templates.U):
    """UnitDefinitions"""

    per_hr = UnitDefinition("per_hr", "1/hr")
    mg_per_min = UnitDefinition("mg_per_min", "mg/min")


_m = Model(
    "dapagliflozin_intestine",
    name="Model for dapagliflozin absorption in the small intestine",
    notes="""
    # Model for dapagliflozin absorption

    - absorption dapagliflozin (dap)
    - fraction absorbed: 0.84 (16% in feces) [Callegari2021]
    - enterohepatische circulation (minimal) -> not modeled
    - slower absorption when fed (compared to fasted) [Kasichayanula2011b] -> via absorption rate depending on food
    => absorption dap, no enterohepatic circulation
    """
    + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=[
        # tissue
        (BQB.OCCURS_IN, "fma/FMA:45615"),  # gut
        (BQB.OCCURS_IN, "bto/BTO:0000545"),  # gut
        (BQB.OCCURS_IN, "ncit/C12736"),  # intestine
        (BQB.OCCURS_IN, "fma/FMA:7199"),  # intestine
        (BQB.OCCURS_IN, "bto/BTO:0000648"),  # intestine

        (BQB.HAS_PROPERTY, "ncit/C79369"),  # Pharmacokinetics: Absorption
        (BQB.HAS_PROPERTY, "ncit/C79372"),  # Pharmacokinetics: Excretion
    ] + templates.model_annotations
)

_m.compartments = [
    Compartment(
        "Vext",
        1.0,
        name="plasma",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        port=True,
        annotations=annotations.compartments["plasma"],
    ),
    Compartment(
        "Vgu",
        1.2825,  # 0.0171 [l/kg] * 75 kg
        name="intestine",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        port=True,
        annotations=annotations.compartments["gu"],
    ),
    Compartment(
        "Vlumen",
        1.2825 * 0.9,  # 0.0171 [l/kg] * 75 kg * 0.9,
        name="intestinal lumen (inner part of intestine)",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        port=True,
        annotations=annotations.compartments["gu_lumen"],
    ),
    Compartment(
        "Vfeces",
        metaId="meta_Vfeces",
        value=1,
        unit=U.liter,
        name="feces",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        port=True,
        annotations=annotations.compartments["feces"],
    ),
    # Compartment(
    #     "Ventero",
    #     1.0,
    #     name="intestinal lining (enterocytes)",
    #     sboTerm=SBO.PHYSICAL_COMPARTMENT,
    #     unit=U.liter,
    # ),
    Compartment(
        "Vapical",
        np.nan,
        name="apical membrane (intestinal membrane enterocytes)",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.m2,
        spatialDimensions=2,
        annotations=annotations.compartments["apical"],
    ),
    # Compartment(
    #     "Vbaso",
    #     np.nan,
    #     name="basolateral membrane (intestinal membrane enterocytes)",
    #     sboTerm=SBO.PHYSICAL_COMPARTMENT,
    #     unit=U.m2,
    #     spatialDimensions=2,
    #     annotations=annotations.compartments["basolateral"],
    # ),
    Compartment(
        "Vstomach",
        metaId="meta_Vstomach",
        value=1,
        unit=U.liter,
        constant=True,
        name="stomach",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        port=True,
        annotations=annotations.compartments["stomach"],
    ),
]


_m.species = [
    Species(
        f"dap_stomach",
        metaId=f"meta_dap_stomach",
        initialConcentration=0.0,
        compartment="Vstomach",
        substanceUnit=U.mmole,
        name=f"dapagliflozin (stomach)",
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["dap"],
        boundaryCondition=True,
    ),
    Species(
        "dap_lumen",
        initialConcentration=0.0,
        name="dapagliflozin (intestinal volume)",
        compartment="Vlumen",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["dap"],
        port=True,
    ),
    Species(
        "dap_ext",
        initialConcentration=0.0,
        name="dapagliflozin (plasma)",
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["dap"],
        port=True,
    ),
    Species(
        "dap_feces",
        initialConcentration=0.0,
        name="dapagliflozin (feces)",
        compartment="Vfeces",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["dap"],
        port=True,
    ),
]

_m.parameters = [
    Parameter(
        f"F_dap_abs",
        0.84,
        U.dimensionless,
        constant=True,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        name=f"fraction absorbed dap",
        notes="""
        Fraction absorbed, i.e., only a fraction of the dapagliflozin in the intestinal lumen
        is absorbed. This parameter determines how much of the dapagliflozin is excreted.
        
        `F_dap_abs` of dose is absorbed. `(1-F_dap_abs)` is excreted in feces.
        
        ~16 % in feces
        """,
    ),
    Parameter(
        "DAPABS_k",
        0.059464824495600456,
        unit=U.per_min,
        name="absorption rate dap",
        sboTerm=SBO.KINETIC_CONSTANT,
    ),
    Parameter(
        "f_absorption",
        1,
        unit=U.dimensionless,
        name="scaling factor absorption rate dap",
        sboTerm=SBO.KINETIC_CONSTANT,
        notes="""1.0: normal absorption corresponding to tablet under fasting conditions.
        
        food decreases the absorption rate, i.e. < 1.0.
        """
    ),
]

_m.rules.append(
    AssignmentRule(
        "absorption",
        value="f_absorption * DAPABS_k * Vgu * dap_lumen",
        unit=U.mmole_per_min,
        name="absorption dapagliflozin",
    ),
)

_m.reactions = [
    Reaction(
        "DAPABS",
        name="absorption dapagliflozin",
        equation="dap_lumen -> dap_ext",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vapical",
        formula=("F_dap_abs * absorption", U.mmole_per_min),
    ),

    Reaction(
        sid="DAPEXC",
        name=f"excretion dapagliflozin (feces)",
        compartment="Vlumen",
        equation=f"dap_lumen -> dap_feces",
        sboTerm=SBO.TRANSPORT_REACTION,
        formula=(
            f"(1 dimensionless - F_dap_abs) * absorption",
            U.mmole_per_min,
        ),
    ),
]


_m.parameters.extend([
    Parameter(
        f"PODOSE_dap",
        0,
        U.mg,
        constant=False,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        name=f"oral dose dapagliflozin [mg]",
        port=True,
    ),
    Parameter(
        f"Ka_dis_dap",
        0.8484201414414877,
        U.per_hr,
        constant=True,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        name=f"dissolution rate dap",
        port=True
    ),
    Parameter(
        f"Mr_dap",
        408.873,
        U.g_per_mole,
        constant=True,
        name=f"Molecular weight dapagliflozin [g/mole]",
        sboTerm=SBO.MOLECULAR_MASS,
        port=True,
    ),
])

# -------------------------------------
# Dissolution of tablet/dose in stomach
# -------------------------------------
_m.reactions.extend(
    [
        # fraction dose available for absorption from stomach
        Reaction(
            sid=f"dissolution_dap",
            name=f"dissolution dapagliflozin",
            formula=(
                f"Ka_dis_dap/60 min_per_hr * PODOSE_dap/Mr_dap",
                U.mmole_per_min,
            ),
            equation=f"dap_stomach -> dap_lumen",
            compartment="Vgu",
            notes="""Swallowing, dissolution of tablet, and transport into intestine.
            Overall process describing the rates of this processes.
            """
        ),
    ]
)
_m.rate_rules.append(
    RateRule(f"PODOSE_dap", f"-dissolution_dap * Mr_dap", U.mg_per_min),
)

model_intestine = _m


def dapagliflozin_layout(dx=200, dy=200) -> pd.DataFrame:
    """Layout definition."""

    delta_y = 0.5 * dy
    delta_x = 1.0 * dx

    positions = [
        # sid, x, y
        ["dap_stomach", 0, 0],
        ["dissolution_dap", 0, delta_y],
        ["dap_lumen", 0, 2 * delta_y],
        ["DAPABS", 0, 3 * delta_y],
        ["dap_ext", 0, 4 * delta_y],

        ["DAPEXC", 1 * delta_x, 3 * delta_y],
        ["dap_feces", 1 * delta_x, 4 * delta_y],

    ]
    df = pd.DataFrame(positions, columns=["id", "x", "y"])
    df.set_index("id", inplace=True)
    return df


def dapagliflozin_annotations(dx=200, dy=200) -> list:
    COLOR_STOMACH = "#1f77b4"   # tab:blue
    COLOR_INTESTINE = "#FFFFFF"
    COLOR_BLOOD = "#FF796C"
    COLOR_FECES = "#8c564b"  # tab:brown

    kwargs = {
        "type": cyviz.AnnotationShapeType.ROUND_RECTANGLE,
        "opacity": 20,
        "border_color": "#000000",
        "border_thickness": 2,
    }

    dy = 0.5 * dy
    dx = 1.0 * dx

    annotations = [
        cyviz.AnnotationShape(
            x_pos=-0.5 * dx, y_pos=-0.5 * dy, width=dx, height=2 * dy,
            fill_color=COLOR_STOMACH, **kwargs
        ),
        cyviz.AnnotationShape(
            x_pos=-0.5 *dx, y_pos=1.5 * dy, width=dx, height=2* dy,
            fill_color=COLOR_INTESTINE, **kwargs
        ),
        cyviz.AnnotationShape(
            x_pos=-0.5 * dx, y_pos=3.5 * dy, width=dx, height=1 * dy,
            fill_color=COLOR_BLOOD, **kwargs
        ),
        cyviz.AnnotationShape(
            x_pos=0.5 * dx, y_pos=2.5 * dy, width=dx, height=2 * dy,
            fill_color=COLOR_FECES, **kwargs
        ),
    ]
    return annotations

if __name__ == "__main__":
    from pkdb_models.models.dapagliflozin import MODEL_BASE_PATH

    from sbmlutils import cytoscape as cyviz
    results = create_model(
        filepath=MODEL_BASE_PATH / f"{model_intestine.sid}.xml",
        model=model_intestine, sbml_level=3, sbml_version=2
    )

    cyviz.visualize_sbml(sbml_path=results.sbml_path, delete_session=True)
    cyviz.apply_layout(layout=dapagliflozin_layout())
    cyviz.add_annotations(annotations=dapagliflozin_annotations())