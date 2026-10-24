"""The tolerances of the integrator for the dapagliflozin model."""

from sbmlsim.model.tolerances import AbsoluteTolerance


def absolute_tolerance(substance: float, other: float) -> AbsoluteTolerance:
    """Get the absolute tolerance of the integrator for the dapagliflozin model.

    The species of dapagliflozin and its metabolites are small: the plasma
    concentrations peak at about 1e-6 to 1e-2 mM (0.1 to 500 mg) and fall by
    orders of magnitude in the terminal phase, and the renal threshold for
    glucose depends on the plasma concentration of the kidney with a
    half-maximal effect at about 4e-6 mM (`KI__RTG_E50`); the amounts in the
    urine and the feces are about 1e-6 to 1 mmole, down to 1e-12 mmole at the
    smallest dose of 0.001 mg. The glucose in the urine (about 1e0 to 1e3
    mmole) is controlled by the relative tolerance. The other states are the
    doses in mg (`PODOSE_dap`, `IVDOSE_dap`, the targets of the rate rules).

    Args:
        substance: tolerance of the species, amounts [mmole] and concentrations
            [mM] (multiplied by the volume of the compartment).
        other: tolerance of the targets of the rate rules [mg].

    Returns:
        The absolute tolerance per kind of state.
    """
    return AbsoluteTolerance(amount=substance, concentration=substance, other=other)


#: the absolute tolerance of the simulation experiments
ABSOLUTE_TOLERANCE = absolute_tolerance(substance=1e-14, other=1e-10)

#: the absolute tolerance of the fits, which integrate with a relative
#: tolerance of 1e-6
ABSOLUTE_TOLERANCE_FIT = absolute_tolerance(substance=1e-12, other=1e-6)
