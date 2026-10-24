"""Liver model of dapagliflozin: hepatic metabolism."""

from sbmlutils.factory import (
    Compartment,
    Model,
    Parameter,
    Reaction,
    Species,
    create_model,
)
from sbmlutils.metadata import (
    BQB,
    SBO,
)

from pkdb_models.models.dapagliflozin.models import annotations, templates
from pkdb_models.models.dapagliflozin.models.fits import FITS


class U(templates.U):
    """Units of the liver model."""


mid = "dapagliflozin_liver"

_m = Model(
    sid=mid,
    name="Model for the hepatic metabolism of dapagliflozin.",
    notes="""
    Model for the hepatic metabolism of dapagliflozin.

    - Dapagliflozin is metabolized mainly by UGT1A9 to the inactive
      dapagliflozin-3-O-glucuronide (D3G, 61 % of an oral dose in the urine),
      and to a minor extent to the 2-O-glucuronide (UGT2B4/2B7) and oxidative
      metabolites (CYP3A4), about 13 % of the dose in the urine and 6 % in the
      feces [Obermeier2010]. Shares of the metabolism: UGT1A9 79.9 %, CYP3A4
      9.8 %, UGT2B4/2B7 8.9 % [Jo2021].
    - The metabolism is a first-order clearance of the plasma of the liver:
      the pharmacokinetics are dose proportional from 0.1 to 500 mg
      [Komoroski2009], far below the Km of UGT1A9 (479 µM in kidney microsomes
      [Kasichayanula2013]). It scales with the volume of the liver tissue, so
      with the tissue loss in cirrhosis.
    - The fraction `f_d3g` of the metabolism forms D3G, the rest the other
      metabolites (`dapmet`).
    - The exposure rises with the Child-Pugh class (AUC ratio 1.03, 1.36, 1.67
      [Kasichayanula2011]).
    """
    + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=[
        *annotations.model,
        # tissue
        (BQB.OCCURS_IN, "fma/FMA:7197"),
        (BQB.OCCURS_IN, "bto/BTO:0000759"),
        (BQB.OCCURS_IN, "NCIT:C12392"),
        (BQB.HAS_PROPERTY, "NCIT:C79371"),  # Pharmacokinetics: Metabolism
    ],
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
        port=True,
    ),
    Compartment(
        "Vli",
        value=1.5,
        unit=U.liter,
        name="liver",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["li"],
        port=True,
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
        port=True,
    ),
    Species(
        "d3g_ext",
        name="dapagliflozin-3-O-glucuronide (plasma)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["d3g"],
        port=True,
    ),
    Species(
        "dapmet_ext",
        name="other dapagliflozin metabolites (plasma)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["dapmet"],
        port=True,
    ),
]

_m.parameters = [
    Parameter(
        "DAPMET_k",
        **FITS.parameter("DAPMET_k", unit=U.per_min),
        name="rate of the hepatic metabolism of dapagliflozin",
        sboTerm=SBO.KINETIC_CONSTANT,
        notes="""First-order rate of the hepatic metabolism relative to the
        volume of the liver tissue; the whole-body model sets it.""",
        port=True,
    ),
    Parameter(
        "f_d3g",
        0.77,
        U.dimensionless,
        name="fraction of the metabolism to D3G",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        notes="""Fraction of the metabolism of dapagliflozin which forms D3G;
        the whole-body model sets it.""",
        port=True,
    ),
    Parameter(
        "f_ugt1a9",
        1,
        U.dimensionless,
        name="scaling factor UGT1A9 activity",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        notes="""Scaling factor of the UGT1A9 activity of the liver (and the
        kidney, the parameter of the whole-body model replaces both).
        1.0: normal activity; <1.0 decreased activity; >1.0 increased activity.
        """,
        port=True,
    ),
]

_m.reactions = [
    Reaction(
        sid="DAP2D3G",
        name="dapagliflozin glucuronidation (DAP2D3G, UGT1A9)",
        equation="dap_ext -> d3g_ext",
        compartment="Vli",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        notes="""Uptake of the dapagliflozin of the plasma of the liver and its
        glucuronidation to D3G by UGT1A9; the export of D3G into the plasma is
        fast.""",
        formula=(
            "f_ugt1a9 * f_d3g * DAPMET_k * Vli * dap_ext",
            U.mmole_per_min,
        ),
        annotations=[
            (BQB.IS, "uniprot/O60656"),  # UGT1A9
        ],
    ),
    Reaction(
        sid="DAP2MET",
        name="dapagliflozin other metabolism (DAP2MET)",
        equation="dap_ext -> dapmet_ext",
        compartment="Vli",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        notes="""The minor pathways of the metabolism of dapagliflozin, the
        2-O-glucuronide (UGT2B4/2B7) and the oxidative metabolites (CYP3A4),
        lumped into one metabolite.""",
        formula=(
            "(1 dimensionless - f_d3g) * DAPMET_k * Vli * dap_ext",
            U.mmole_per_min,
        ),
    ),
]

model_liver = _m


if __name__ == "__main__":
    from sbmlode import OdeSystem

    from pkdb_models.models.dapagliflozin import MODEL_BASE_PATH

    results = create_model(
        model=model_liver,
        filepath=MODEL_BASE_PATH / f"{model_liver.sid}.xml",
        sbml_level=3,
        sbml_version=2,
    )

    # ODE equations
    ode_system = OdeSystem.from_sbml(results.sbml_path)
    ode_system.write(results.sbml_path.parent / f"{results.sbml_path.stem}.md")
