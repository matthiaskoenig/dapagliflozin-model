"""Kidney model of dapagliflozin: renal excretion, renal metabolism and glucose excretion."""

from sbmlutils.factory import (
    AssignmentRule,
    Compartment,
    Model,
    Parameter,
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
    """Units of the kidney model."""

    mg_per_g = UnitDefinition("mg_per_g", "mg/g")
    ml_per_min = UnitDefinition("ml_per_min", "ml/min")
    l_per_ml = UnitDefinition("l_per_ml", "l/ml")


mid = "dapagliflozin_kidney"

_m = Model(
    sid=mid,
    name="Model for the renal excretion of dapagliflozin and its metabolites and the urinary glucose excretion.",
    notes="""
    Model for the renal excretion of dapagliflozin and its metabolites, the
    renal metabolism of dapagliflozin and the urinary glucose excretion.

    - 1.2-1.6 % of an oral dose of dapagliflozin is excreted unchanged in the
      urine with a renal clearance of 3.5-5.6 ml/min, below the filtration of
      the unbound fraction (fu 0.09, 9 ml/min): net reabsorption
      [Obermeier2010, Kasichayanula2013].
    - 61 % of the dose is excreted as dapagliflozin-3-O-glucuronide (D3G) in
      the urine with a renal clearance of 151 ml/min in healthy subjects, above
      the GFR (tubular secretion) [Kasichayanula2013]; 75.2 % of the
      radioactivity is in the urine [Obermeier2010].
    - The renal excretion is relative to the glomerular filtration rate (GFR).
    - The kidney glucuronidates dapagliflozin (UGT1A9, 8-fold mRNA and 3-fold
      D3G formation per mg microsomal protein of the liver [Kasichayanula2013]);
      the renal metabolism scales with the renal function (the functional renal
      mass), which makes the exposure of dapagliflozin rise in renal impairment
      (AUC ratio 1.60, 2.07, 1.95 mild, moderate, severe [Kasichayanula2013]).
    - Dapagliflozin lowers the renal threshold of glucose (RTG) by the
      inhibition of SGLT2; glucose above the RTG is excreted with the GFR.
    """
    + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=[
        *annotations.model,
        # tissue
        (BQB.OCCURS_IN, "fma/FMA:7203"),  # kidney
        (BQB.OCCURS_IN, "bto/BTO:0000671"),  # kidney
        (BQB.OCCURS_IN, "NCIT:C12415"),  # kidney
        (BQB.HAS_PROPERTY, "NCIT:C79372"),  # Pharmacokinetics: Excretion
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
        "Vki",
        value=0.3,  # 0.4 % of bodyweight
        unit=U.liter,
        name="kidney",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["ki"],
        port=True,
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

# ---------------------------------------------------------------------------------------------------------------------
# Pharmacokinetics
# ---------------------------------------------------------------------------------------------------------------------
#: the substances of the plasma and the urine: id, name
_SUBSTANCES: list[tuple[str, str]] = [
    ("dap", "dapagliflozin"),
    ("d3g", "dapagliflozin-3-O-glucuronide"),
    ("dapmet", "other dapagliflozin metabolites"),
]
_m.species = []
for _sid, _name in _SUBSTANCES:
    _m.species.extend(
        [
            Species(
                f"{_sid}_ext",
                name=f"{_name} (plasma)",
                initialConcentration=0.0,
                compartment="Vext",
                substanceUnit=U.mmole,
                hasOnlySubstanceUnits=False,  # this is a concentration
                sboTerm=SBO.SIMPLE_CHEMICAL,
                annotations=annotations.species[_sid],
                port=True,
            ),
            Species(
                f"{_sid}_urine",
                name=f"{_name} (urine)",
                initialConcentration=0.0,
                compartment="Vurine",
                substanceUnit=U.mmole,
                hasOnlySubstanceUnits=True,
                sboTerm=SBO.SIMPLE_CHEMICAL,
                annotations=annotations.species[_sid],
                port=True,
            ),
        ]
    )

_m.parameters.extend(
    [
        Parameter(
            "f_renal_function",
            name="renal function",
            value=1.0,
            unit=U.dimensionless,
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            >1.0: increased kidney function
            1.0: normal kidney function (healthy control)
            <1.0: decreased kidney function
            Renal function relative to the healthy glomerular filtration rate
            GFR_base, which is 101 ml/min at 75 kg; a measured creatinine
            clearance or GFR of X ml/min is the renal function X/101. It scales
            the filtration (renal excretion, glucose excretion) and the renal
            metabolism (the functional renal mass).
            """,
        ),
        Parameter(
            "GFR_base",
            101,
            U.ml_per_min,
            name="glomerular filtration rate (healthy)",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            port=True,
        ),
        Parameter(
            "fu_dap",
            0.09,
            U.dimensionless,
            name="fraction unbound in plasma dapagliflozin",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""Dapagliflozin is 91 % bound to plasma proteins
            [Obermeier2010] (92.3-93 % [Kasichayanula2013]), independent of the
            renal and hepatic function.""",
        ),
        Parameter(
            "DAPEX_f",
            **FITS.parameter("KI__DAPEX_f", unit=U.dimensionless),
            name="renal clearance of dapagliflozin relative to its filtration",
            sboTerm=SBO.KINETIC_CONSTANT,
            notes="""
            The renal clearance of dapagliflozin is `DAPEX_f` times the
            filtration of its unbound fraction, fu_dap * GFR (9 ml/min). The
            renal clearance of 3.5-5.6 ml/min [Kasichayanula2013, Obermeier2010]
            is about half of the filtration: net tubular reabsorption.
            """,
        ),
        Parameter(
            "D3GEX_f",
            **FITS.parameter("KI__D3GEX_f", unit=U.dimensionless),
            name="renal clearance of the metabolites relative to the GFR",
            sboTerm=SBO.KINETIC_CONSTANT,
            notes="""
            The renal clearance of D3G is `D3GEX_f` times the GFR: 151 ml/min in
            healthy subjects [Kasichayanula2013], above the GFR (tubular
            secretion, D3G is an OAT3 substrate). The other metabolites (the
            2-O-glucuronide and the oxidative metabolites) are excreted in the
            same way. The plasma kinetics of the metabolites are limited by their
            formation.
            """,
        ),
        Parameter(
            "DAPMET_k",
            0.1,
            U.per_min,
            name="rate of the renal metabolism of dapagliflozin",
            sboTerm=SBO.KINETIC_CONSTANT,
            notes="""
            First-order rate of the renal metabolism relative to the volume of
            the kidney; the whole-body model sets it from the hepatic rate and
            the renal fraction of the metabolism.
            """,
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
            notes="""Scaling factor of the UGT1A9 activity of the kidney (and the
            liver, the parameter of the whole-body model replaces both).
            1.0: normal activity; <1.0 decreased activity; >1.0 increased activity.
            """,
            port=True,
        ),
    ]
)

_m.rules.extend(
    [
        AssignmentRule(
            "GFR",
            "f_renal_function * GFR_base",
            unit=U.ml_per_min,
            name="glomerular filtration rate",
        ),
    ]
)

_m.reactions = [
    Reaction(
        sid="DAP2D3G",
        name="dapagliflozin glucuronidation (DAP2D3G, UGT1A9)",
        equation="dap_ext -> d3g_ext",
        compartment="Vki",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        notes="""The glucuronidation of the dapagliflozin of the renal plasma to
        D3G by UGT1A9; uptake and export are fast. It scales with the renal
        function (the functional renal mass).""",
        formula=(
            "f_renal_function * f_ugt1a9 * f_d3g * DAPMET_k * Vki * dap_ext",
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
        compartment="Vki",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        notes="""The minor pathways of the renal metabolism of dapagliflozin
        (2-O-glucuronide by UGT2B7), lumped into one metabolite.""",
        formula=(
            "f_renal_function * (1 dimensionless - f_d3g) * DAPMET_k * Vki * dap_ext",
            U.mmole_per_min,
        ),
    ),
    Reaction(
        sid="DAPEX",
        name="dapagliflozin excretion (DAPEX)",
        equation="dap_ext -> dap_urine",
        compartment="Vki",
        sboTerm=SBO.TRANSPORT_REACTION,
        # [ml/min] * [l/ml] * [-] * [-] * [mmole/l]
        formula=(
            "GFR * 0.001 l_per_ml * fu_dap * DAPEX_f * dap_ext",
            U.mmole_per_min,
        ),
    ),
    Reaction(
        sid="D3GEX",
        name="D3G excretion (D3GEX)",
        equation="d3g_ext -> d3g_urine",
        compartment="Vki",
        sboTerm=SBO.TRANSPORT_REACTION,
        formula=(
            "GFR * 0.001 l_per_ml * D3GEX_f * d3g_ext",
            U.mmole_per_min,
        ),
    ),
    Reaction(
        sid="METEX",
        name="other metabolites excretion (METEX)",
        equation="dapmet_ext -> dapmet_urine",
        compartment="Vki",
        sboTerm=SBO.TRANSPORT_REACTION,
        formula=(
            "GFR * 0.001 l_per_ml * D3GEX_f * dapmet_ext",
            U.mmole_per_min,
        ),
    ),
]

# ---------------------------------------------------------------------------------------------------------------------
# Pharmacodynamics
# ---------------------------------------------------------------------------------------------------------------------
_m.species.extend(
    [
        Species(
            "fpg",
            name="fasting plasma glucose (FPG)",
            initialConcentration=5.0,
            compartment="Vext",
            substanceUnit=U.mmole,
            hasOnlySubstanceUnits=False,  # this is a concentration
            boundaryCondition=True,
            sboTerm=SBO.SIMPLE_CHEMICAL,
            annotations=annotations.species["glc"],
            notes="""
            Fasting plasma glucose (FPG), set by the studies. It shifts the RTG
            (SGLT2 is induced in type 2 diabetes) and, with the difference of
            the mean plasma glucose of the day to the FPG (`PG_delta`), gives
            the glucose which is filtered.
            """,
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
            port=True,
        ),
    ]
)

_m.parameters.extend(
    [
        Parameter(
            "Mr_glc",
            180.16,
            U.g_per_mole,
            name="molecular weight glucose",
            sboTerm=SBO.MOLECULAR_MASS,
        ),
        Parameter(
            "cf_mg_per_g",
            1000,
            U.mg_per_g,
            name="conversion factor mg per g",
        ),
        Parameter(
            "fpg_healthy",
            5.0,
            U.mM,
            name="fasting plasma glucose (healthy)",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "PG_delta",
            0.85,
            U.mM,
            name="mean plasma glucose of the day minus the fasting plasma glucose",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            The glucose which is filtered over the day is the mean plasma glucose
            of the day, above the fasting plasma glucose by the meals. In healthy
            subjects about 0.85 mM: the urinary glucose excretion of 285 mmol/day
            over a GFR of 139 l/day with a 24-hr mean RTG of 3.79 mM is a mean
            plasma glucose of 5.85 mM [Sha2015]. Studies which report the 24-hr
            mean plasma glucose set the difference to the FPG.
            """,
        ),
        Parameter(
            "RTG_E50",
            **FITS.parameter("KI__RTG_E50", unit=U.mM),
            name="EC50 reduction in RTG",
            sboTerm=SBO.KINETIC_CONSTANT,
            notes="""
            EC50 of the reduction of the RTG by plasma dapagliflozin (total).
            The IC50 of dapagliflozin for human SGLT2 is about 1.1 nM (unbound),
            about 5 ng/ml (1.2e-5 mM) total with fu 0.09; SGLT2 is inhibited
            from the tubular lumen. The urinary glucose excretion rises steeply
            from 0.3 mg (0.6 g/24 hr) to 1 mg (19 g/24 hr) [Gould2013].
            """,
        ),
        Parameter(
            "RTG_gamma",
            **FITS.parameter("KI__RTG_gamma", unit=U.dimensionless),
            name="hill coefficient reduction in RTG",
            sboTerm=SBO.KINETIC_CONSTANT,
            notes="""The steepness of the reduction of the RTG with the
            concentration, at least 1 (an exponent below 1 makes the
            integration stiff at zero dapagliflozin). With an exponent of 1 the
            fits give the steep dose response from 0.3 to 2.5 mg only with an
            RTG at the upper bound of its range (14 mM, PD_20261020_010757__2bd20).""",
        ),
        Parameter(
            "RTG_base",
            11.0,
            U.mM,
            name="baseline RTG",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""RTG without SGLT2 inhibition in healthy subjects (FPG 5 mM),
            which corresponds to the SGLT2 capacity; about 10-11 mM (180-200
            mg/dl). Not fitted: the urinary glucose of the fits bounds it only
            from below (10 mM, PD_20261020_012618__f9531).""",
        ),
        Parameter(
            "RTG_m_fpg",
            **FITS.parameter("KI__RTG_m_fpg", unit=U.dimensionless),
            name="FPG effect on RTG",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""Increase of the RTG per mM of FPG above the healthy FPG:
            SGLT2 is induced in type 2 diabetes [Rahmoune2005]. Phenomenological:
            it gives the urinary glucose excretion of type 2 diabetes (38-77 g
            with 2.5-20 mg [Kasichayanula2011a, FDAMB102003]) with the constant
            glucose of the model, which has no fall of the glucose during the
            treatment.""",
        ),
        Parameter(
            "RTG_max_inhibition",
            **FITS.parameter("KI__RTG_max_inhibition", unit=U.dimensionless),
            name="maximal relative reduction of the RTG",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""Maximal reduction of the RTG by SGLT2 inhibition relative to
            the RTG without inhibition: SGLT1 reabsorbs the glucose which SGLT2
            does not; the RTG under full inhibition is about 2-4 mM.""",
        ),
        Parameter(
            "RTG_splay",
            1.0,
            U.mM,
            name="splay of the glucose excretion around the RTG",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""The thresholds of the nephrons differ, so the glucose
            excretion rises smoothly around the RTG (splay) instead of at a
            sharp threshold; the softplus with a width of 1 mM.""",
        ),
    ]
)

_m.rules.extend(
    [
        AssignmentRule(
            "pg",
            "fpg + PG_delta",
            U.mM,
            name="mean plasma glucose of the day",
            notes="""The glucose which is filtered.""",
        ),
        AssignmentRule(
            "RTG_fpg",
            "RTG_base + RTG_m_fpg * (fpg - fpg_healthy)",
            U.mM,
            name="RTG without SGLT2 inhibition",
            notes="""
            FPG dependent RTG without SGLT2 inhibition. SGLT2 is induced in T2DM
            [Rahmoune2005].
            """,
        ),
        AssignmentRule(
            "RTG_delta",
            "RTG_fpg * RTG_max_inhibition",
            U.mM,
            name="maximal reduction of the RTG",
        ),
        AssignmentRule(
            variable="RTG",
            value="RTG_fpg - RTG_delta * power(dap_ext, RTG_gamma)/ (power(RTG_E50, RTG_gamma) + power(dap_ext, RTG_gamma))",
            unit=U.mM,
            name="renal threshold glucose (RTG)",
            notes="""
            Renal threshold glucose.

            The renal threshold for glucose (RTG) is the plasma glucose
            concentration at which the tubular reabsorption of glucose begins to
            saturate; glucose is excreted into the urine in proportion to the
            glucose concentration above this threshold.
            """,
        ),
    ]
)

# Glucose excretion (UGE)
_m.reactions.extend(
    [
        Reaction(
            sid="GLCEX",
            name="glucose excretion (GLCEX)",
            equation="fpg -> glc_urine [dap_ext]",
            compartment="Vki",
            sboTerm=SBO.TRANSPORT_REACTION,
            notes="""The filtered glucose above the RTG, the softplus of the
            difference with the width `RTG_splay`: the excretion approaches
            GFR * (pg - RTG) above the RTG and zero below it.""",
            # [ml/min] * [l/ml] * [mmole/l] = [mmole/min]
            formula=(
                "GFR * 0.001 l_per_ml * RTG_splay * ln(1 dimensionless + exp((pg - RTG)/RTG_splay))",
                U.mmole_per_min,
            ),
        ),
    ]
)

_m.rules.extend(
    [
        AssignmentRule(
            "UGE",
            "glc_urine * Mr_glc/cf_mg_per_g",
            unit=U.gram,
            name="urinary glucose excretion (UGE)",
            notes="""
            Urinary glucose excretion, the cumulative amount of glucose in the
            urine.
            """,
        )
    ]
)

model_kidney = _m


if __name__ == "__main__":
    from sbmlode import OdeSystem

    from pkdb_models.models.dapagliflozin import MODEL_BASE_PATH

    results = create_model(
        model=model_kidney,
        filepath=MODEL_BASE_PATH / f"{model_kidney.sid}.xml",
        sbml_level=3,
        sbml_version=2,
    )

    # ODE equations
    ode_system = OdeSystem.from_sbml(results.sbml_path)
    ode_system.write(results.sbml_path.parent / f"{results.sbml_path.stem}.md")
