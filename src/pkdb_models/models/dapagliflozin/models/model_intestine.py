"""Intestine model of dapagliflozin: dissolution, absorption and fecal excretion."""

import numpy as np
from sbmlutils.factory import (
    AssignmentRule,
    Compartment,
    Model,
    Parameter,
    RateRule,
    Reaction,
    Species,
    UnitDefinition,
    create_model,
)
from sbmlutils.metadata import (
    BQB,
    SBO,
)

from pkdb_models.models.dapagliflozin.models import annotations, templates
from pkdb_models.models.dapagliflozin.models.fits import FITS


class U(templates.U):
    """Units of the intestine model."""

    per_hr = UnitDefinition("per_hr", "1/hr")
    mg_per_min = UnitDefinition("mg_per_min", "mg/min")


_m = Model(
    "dapagliflozin_intestine",
    name="Model for dapagliflozin absorption in the small intestine",
    notes="""
    # Model for dapagliflozin absorption

    - Dissolution of the oral dose and absorption of dapagliflozin; the
      fraction absorbed `F_dap_abs` of the dapagliflozin of the lumen is
      absorbed, the rest is excreted in the feces. After an oral solution of
      50 mg of 14C-dapagliflozin 21 % of the radioactivity is recovered in the
      feces, 15.4 % of the dose as unchanged dapagliflozin [Obermeier2010];
      the absolute bioavailability is 78 % [Boulton2013].
    - The enterohepatic circulation is minimal and not modeled.
    - Food (`fed`) slows the absorption (`FED_f_absorption`): a high fat meal
      lowers Cmax by 31 % and delays tmax by about 1 hr, the AUC is unchanged
      [Kasichayanula2011b].
    """
    + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=[
        *annotations.model,
        # tissue
        (BQB.OCCURS_IN, "bto/BTO:0000545"),  # gut
        (BQB.OCCURS_IN, "NCIT:C12736"),  # intestine
        (BQB.OCCURS_IN, "fma/FMA:7199"),  # intestine
        (BQB.OCCURS_IN, "bto/BTO:0000648"),  # intestine
        (BQB.HAS_PROPERTY, "NCIT:C79369"),  # Pharmacokinetics: Absorption
        (BQB.HAS_PROPERTY, "NCIT:C79372"),  # Pharmacokinetics: Excretion
    ],
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
        constant=False,
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
    Compartment(
        "Vapical",
        np.nan,
        name="apical membrane (intestinal membrane enterocytes)",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.m2,
        spatialDimensions=2,
        annotations=annotations.compartments["apical"],
    ),
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
        "dap_stomach",
        metaId="meta_dap_stomach",
        initialConcentration=0.0,
        compartment="Vstomach",
        substanceUnit=U.mmole,
        name="dapagliflozin (stomach)",
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["dap"],
        boundaryCondition=True,
    ),
    Species(
        "dap_lumen",
        initialConcentration=0.0,
        name="dapagliflozin (intestinal lumen)",
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
        "F_dap_abs",
        0.84,
        U.dimensionless,
        constant=True,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        name="fraction absorbed dapagliflozin",
        notes="""
        Fraction absorbed: `F_dap_abs` of the dapagliflozin of the intestinal
        lumen is absorbed, `(1 - F_dap_abs)` is excreted in the feces. 15.4 % of
        an oral dose is excreted unchanged in the feces [Obermeier2010]; with
        the hepatic first pass the bioavailability is 0.78 [Boulton2013]. Not
        fitted: only the feces of Obermeier2010 identify it.
        """,
    ),
    Parameter(
        "DAPABS_k",
        **FITS.parameter("GU__DAPABS_k", unit=U.per_min),
        name="rate of dapagliflozin absorption",
        sboTerm=SBO.KINETIC_CONSTANT,
        notes="""
        First-order rate of the dapagliflozin leaving the lumen, absorbed
        (`F_dap_abs`) or to the feces: tmax 1.0-1.25 hr fasted [Boulton2013,
        Kasichayanula2013].
        """,
    ),
    Parameter(
        "fed",
        0.0,
        unit=U.dimensionless,
        name="fed state (0 fasted, 1 fed)",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        notes="""
        The dose is taken in the fasted state (0) or with a meal (1). The effect
        of food on the absorption rate is `FED_f_absorption`.
        """,
    ),
    Parameter(
        "FED_f_absorption",
        **FITS.parameter("GU__FED_f_absorption", unit=U.dimensionless),
        name="fed effect on the absorption rate",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        notes="""
        Factor of the absorption rate with food: a high fat meal lowers Cmax by
        31 % (10 mg) and delays tmax by about 1 hr, the AUC is unchanged
        [Kasichayanula2011b]; Cmax -45 to -51 % with 2.5 mg [LaCreta2016]. The
        fraction absorbed is not changed.
        """,
    ),
    Parameter(
        "f_absorption",
        np.nan,
        unit=U.dimensionless,
        constant=False,
        name="factor of the absorption rate",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
    ),
]
_m.rules.extend(
    [
        AssignmentRule(
            "f_absorption",
            "1 dimensionless + fed * (FED_f_absorption - 1 dimensionless)",
            unit=U.dimensionless,
        ),
        AssignmentRule(
            "absorption",
            value="f_absorption * DAPABS_k * Vlumen * dap_lumen",
            unit=U.mmole_per_min,
            name="dapagliflozin leaving the lumen",
        ),
    ]
)

_m.reactions = [
    Reaction(
        "DAPABS",
        name="absorption dapagliflozin (DAPABS)",
        equation="dap_lumen -> dap_ext",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vapical",
        formula=("F_dap_abs * absorption", U.mmole_per_min),
    ),
    Reaction(
        sid="DAPEXC",
        name="excretion dapagliflozin (feces)",
        compartment="Vlumen",
        equation="dap_lumen -> dap_feces",
        sboTerm=SBO.TRANSPORT_REACTION,
        formula=(
            "(1 dimensionless - F_dap_abs) * absorption",
            U.mmole_per_min,
        ),
    ),
]

_m.parameters.extend(
    [
        Parameter(
            "PODOSE_dap",
            0,
            U.mg,
            constant=False,
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            name="oral dose dapagliflozin [mg]",
            port=True,
        ),
        Parameter(
            "Ka_dis_dap",
            **FITS.parameter("GU__Ka_dis_dap", unit=U.per_hr),
            constant=True,
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            name="dissolution rate dapagliflozin",
            port=True,
        ),
        Parameter(
            "Mr_dap",
            408.873,
            U.g_per_mole,
            constant=True,
            name="Molecular weight dapagliflozin [g/mole]",
            sboTerm=SBO.MOLECULAR_MASS,
            port=True,
        ),
    ]
)

# -------------------------------------
# Dissolution of tablet/dose in stomach
# -------------------------------------
_m.reactions.extend(
    [
        # fraction dose available for absorption from stomach
        Reaction(
            sid="dissolution_dap",
            name="dissolution dapagliflozin",
            formula=(
                "Ka_dis_dap/60 min_per_hr * PODOSE_dap/Mr_dap",
                U.mmole_per_min,
            ),
            equation="dap_stomach -> dap_lumen",
            compartment="Vgu",
            notes="""Swallowing, dissolution of tablet, and transport into intestine.
            Overall process describing the rates of this processes.
            """,
        ),
    ]
)
_m.rate_rules.append(
    RateRule("PODOSE_dap", "-dissolution_dap * Mr_dap", U.mg_per_min),
)

model_intestine = _m


if __name__ == "__main__":
    from sbmlode import OdeSystem

    from pkdb_models.models.dapagliflozin import MODEL_BASE_PATH

    results = create_model(
        filepath=MODEL_BASE_PATH / f"{model_intestine.sid}.xml",
        model=model_intestine,
        sbml_level=3,
        sbml_version=2,
    )

    # ODE equations
    ode_system = OdeSystem.from_sbml(results.sbml_path)
    ode_system.write(results.sbml_path.parent / f"{results.sbml_path.stem}.md")
