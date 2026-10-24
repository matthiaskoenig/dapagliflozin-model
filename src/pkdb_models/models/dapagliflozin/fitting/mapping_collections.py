"""Parameter fit problems for dapagliflozin.

The data of a fit is selected from the fit mappings of the studies in three
steps, see `sbmlsim.fit.helpers`: the filters on the metadata of the mappings
select the training data of the fit and exclude the rest, the outliers are
named once for all mappings and are the training data which is not usable, and
the validation data is the part of the training data a fit is evaluated on but
not fitted to. The metadata of a mapping describes its curve, what a fit does
with the curve is decided here.

The training data are healthy subjects and patients with type 2 diabetes and a
normal renal function, without a coadministered drug, fasted, fed (the fed
arms identify the effect of food on the absorption) or fasting not reported,
single doses and exactly encoded multiple doses. The validation data has four
parts, each a different question: the interactions (`INTERACTION_VALIDATION`,
the arms with a coadministered drug: rifampicin and mefenamic acid with factors
of UGT1A9 which are not fitted, the other drugs without an interaction), the
doses (`DOSE_VALIDATION`), the populations (`POPULATION_VALIDATION`, the
impaired renal and hepatic function, type 1 diabetes) and the limitations
(`LIMITATION_VALIDATION`, data of a process the model does not describe). The
summaries of data which is fitted one by one and the data of the same subjects
in a second study or arm are excluded (`AGGREGATE_MAPPINGS`).
"""

from pkpdsim.fit import mapping_metadata, observable_filter
from pkpdsim.metadata import Coadministration, Health
from sbmlsim.console import console
from sbmlsim.experiment import SimulationExperiment
from sbmlsim.fit import FitMapping, FitMappingCollection, MappingKind
from sbmlsim.fit.helpers import FitMappings, MappingFilter

from pkdb_models.models.dapagliflozin import DAPAGLIFLOZIN_PATH, DATA_PATHS
from pkdb_models.models.dapagliflozin.experiments.studies import (
    FDAMB102002,
    FDAMB102003,
    FDAMB102006,
    FDAMB102007,
    Boulton2013,
    Cho2021,
    Gould2013,
    Hwang2022a,
    Imamura2013,
    Jang2020,
    Kasichayanula2011,
    Kasichayanula2011a,
    Kasichayanula2011b,
    Kasichayanula2011c,
    Kasichayanula2012,
    Kasichayanula2013,
    Kasichayanula2013a,
    Khomitskaya2018,
    Kim2023,
    Kim2023a,
    Komoroski2009,
    LaCreta2016,
    Obermeier2010,
    Sha2015,
    Shah2019a,
    Watada2019,
    Yang2013,
    vanderAartvanderBeek2020,
)

#: the keys of fit mappings by study: the keys of the mappings are SIds of their
#: experiment, which two studies can share, so a named mapping is named with
#: its study
type NamedMappings = dict[str, set[str]]

EXPERIMENT_CLASSES: list[type[SimulationExperiment]] = [
    Boulton2013,
    Cho2021,
    FDAMB102002,
    FDAMB102003,
    FDAMB102006,
    FDAMB102007,
    Gould2013,
    Hwang2022a,
    Imamura2013,
    Jang2020,
    Kasichayanula2011,
    Kasichayanula2011a,
    Kasichayanula2011b,
    Kasichayanula2011c,
    Kasichayanula2012,
    Kasichayanula2013,
    Kasichayanula2013a,
    Khomitskaya2018,
    Kim2023,
    Kim2023a,
    Komoroski2009,
    LaCreta2016,
    Obermeier2010,
    Sha2015,
    Shah2019a,
    vanderAartvanderBeek2020,
    Watada2019,
    Yang2013,
]

#: health of the training data: healthy subjects and patients with type 2
#: diabetes and a normal renal function (the fasting plasma glucose is set by
#: the experiments)
TRAINING_HEALTH: frozenset[Health] = frozenset({Health.HEALTHY, Health.T2DM})

#: observables of the pharmacokinetics
PK_OBSERVABLES: set[str] = {
    "Cve_dap",
    "Cve_d3g",
    "Cve_daptot",
    "Aurine_dap",
    "Aurine_d3g",
    "Aurine_daptot",
    "Afeces_dap",
    "Afeces_daptot",
}

#: observables of the pharmacodynamics
PD_OBSERVABLES: set[str] = {"KI__RTG", "KI__UGE"}

#: mappings which are not fitted, the data is not usable
OUTLIER_MAPPINGS: NamedMappings = {
    # Fig 1A is labeled dapagliflozin but shows the glucuronide (the values and
    # the time course are the ones of the D3G of Fig 1B and Table 2)
    "Kasichayanula2013": {
        "fm_dap20_d10_normal__Cve_dap_",
        "fm_dap20_d10_mild__Cve_dap_",
        "fm_dap20_d10_moderate__Cve_dap_",
        "fm_dap20_d10_severe__Cve_dap_",
    },
}

#: the validation of the interactions: the arms with a coadministered drug.
#: Rifampicin (UGT1A9 induction, AUC 0.78) and mefenamic acid (UGT1A9
#: inhibition, AUC 1.51) are simulated with factors of UGT1A9 which are not
#: fitted, the other drugs do not change the exposure (AUC ratios 0.98 to
#: 1.08) and are simulated without a factor. `filter_selection` lets the
#: validation data through for this reason, every other coadministration stays
#: out of the selection altogether
INTERACTION_VALIDATION: NamedMappings = {
    # sitagliptin 100 mg, the ABCB1 C3435T groups of the 100 subjects
    "Hwang2022a": {
        "fm_dap10_C3435T_CC",
        "fm_dap10_C3435T_CT",
        "fm_dap10_C3435T_TT",
    },
    # voglibose, not absorbed (AUC ratio 1.01)
    "Imamura2013": {"fm_dap10_DAP10__VOG"},
    # lobeglitazone (AUC ratio 0.99)
    "Jang2020": {"fm_dap10_DAP10__LOB05"},
    # pioglitazone, metformin, glimepiride, sitagliptin (AUC ratios 0.99-1.08)
    "Kasichayanula2011c": {
        "fm_po_DAP50__PIO45_S1",
        "fm_po_DAP20__MET1000_S2",
        "fm_po_DAP20__GLI4_S3",
        "fm_po_DAP20__SIT100_S3",
    },
    # simvastatin and valsartan (AUC ratios 0.98 and 1.02)
    "Kasichayanula2012": {"fm_dap20_DAP20__SIM40", "fm_dap20_DAP20__VAL320"},
    # rifampicin and mefenamic acid, factors of UGT1A9
    "Kasichayanula2013a": {
        "fm_po_dap10_dapagliflozin_RIF_rifampicin",
        "fm_po_dap10_dapagliflozin_cumulative_amount_RIF_rifampicin",
        "fm_po_dap10_dapagliflozin_3_o_glucuronide_RIF_rifampicin",
        "fm_po_dap10_dapagliflozin_3_o_glucuronide_cumulative_amount_RIF_rifampicin",
        "fm_po_dap10_dapagliflozin_MFA_mefenamic_acid",
        "fm_po_dap10_dapagliflozin_cumulative_amount_MFA_mefenamic_acid",
        "fm_po_dap10_dapagliflozin_3_o_glucuronide_MFA_mefenamic_acid",
        "fm_po_dap10_dapagliflozin_3_o_glucuronide_cumulative_amount_MFA_mefenamic_acid",
    },
    # individual tablets of dapagliflozin and metformin, fed
    "Khomitskaya2018": {"fm_dap10_DAP10"},
    # evogliptin
    "Kim2023": {"fm_dap10_EV5__DAP10_dapagliflozin", "fm_dap10_EV5__DAP10_uge"},
    # fixed dose combination with metformin XR, fasted and fed
    "Shah2019a": {"fm_po_dap5_fasted", "fm_po_dap5_fed"},
}

#: the validation of the doses: the urinary dapagliflozin of the single doses
#: of young healthy Japanese men, whose renal function from the
#: Cockcroft-Gault clearance (1.5) lies above the one of the training
DOSE_VALIDATION: NamedMappings = {
    "Kasichayanula2011a": {
        "fm_sad2p5_dapagliflozin_urine",
        "fm_sad10_dapagliflozin_urine",
        "fm_sad20_dapagliflozin_urine",
        "fm_sad50_dapagliflozin_urine",
    },
}

#: the validation of the populations: the impaired renal and hepatic function
#: (the classes of `renal_map` and `cirrhosis_map`) and type 1 diabetes
POPULATION_VALIDATION: NamedMappings = {
    # urinary glucose of patients with type 2 diabetes and renal impairment at
    # steady state, medians (the urine of the same trial is Kasichayanula2013)
    "FDAMB102007": {
        "fm_dap20_normal",
        "fm_dap20_mild",
        "fm_dap20_moderate",
        "fm_dap20_severe",
    },
    # hepatic impairment, Child-Pugh A, B and C
    "Kasichayanula2011": {
        "fm_dap10_mild_dap",
        "fm_dap10_mild_d3g",
        "fm_dap10_moderate_dap",
        "fm_dap10_moderate_d3g",
        "fm_dap10_severe_dap",
        "fm_dap10_severe_d3g",
    },
    # type 2 diabetes with mild, moderate and severe renal impairment
    "Kasichayanula2013": {
        f"fm_{key}"
        for group in ["mild", "moderate", "severe"]
        for key in [
            f"dap50_{group}_dapagliflozin",
            f"dap50_{group}_dapagliflozin_3_o_glucuronide",
            f"dap20_d4_{group}_dapagliflozin",
            f"dap20_d4_{group}_dapagliflozin_3_o_glucuronide",
            f"dap20_d10_{group}_Aurine_dap",
            f"dap20_d10_{group}_Aurine_d3g",
            f"dap20_d10_{group}__Cve_d3g_",
        ]
    },
    # chronic kidney disease at steady state, individual patients
    "vanderAartvanderBeek2020": {f"fm_dap10_S{k}" for k in range(1, 7)},
    # Japanese patients with type 1 diabetes on insulin
    "Watada2019": {
        "fm_po_DAP0_uge",
        "fm_po_DAP5_dapagliflozin",
        "fm_po_DAP10_dapagliflozin",
        "fm_po_DAP5_dapagliflozin_3_o_glucuronide",
        "fm_po_DAP10_dapagliflozin_3_o_glucuronide",
        "fm_po_DAP5_uge",
        "fm_po_DAP10_uge",
    },
}

#: the validation of the limitations of the model: valid data of a process the
#: model does not describe, which would distort the fitted parameters of the
#: processes it does describe
LIMITATION_VALIDATION: NamedMappings = {
    # fecal metabolites (5.6 % of the dose), the model has none
    "FDAMB102006": {"fm_DAP50_total_dapagliflozin_feces"},
    # the baseline glucosuria of 0.04 to 0.09 g a day: the model excretes no
    # glucose at the plasma glucose of healthy subjects
    "Gould2013": {"fm_dap0", "fm_dap0p001", "fm_dap0p01", "fm_dap0p1"},
    # the glucosuria of the placebo of patients with type 2 diabetes, the
    # glucose of the day of the patients is not described
    "Kasichayanula2011a": {"fm_dap0_uge"},
    # the baseline glucosuria of the placebo, 1 to 2 g
    "Komoroski2009": {"fm_po_dap_0_single"},
    # the RTG of the study is the mean of a urine interval, derived from the
    # plasma glucose and the glucose excretion of the interval
    "Sha2015": {"fm_po_dap10_renal_threshold"},
}

#: the summaries of data which the fits use one by one and the data of the same
#: subjects in a second study or arm, excluded from every fit: a mean of
#: subjects would count the subjects twice
AGGREGATE_MAPPINGS: NamedMappings = {
    # the test formulation, the same 37 subjects as the reference
    "Cho2021": {"fm_dap10_RC10"},
    # the urine of trial MB102007, which Kasichayanula2013 maps (Table 2)
    "FDAMB102007": {
        f"fm_dap50_{group}_{substance}_cumulative_amount"
        for group in ["healthy", "normal", "mild", "moderate", "severe"]
        for substance in ["dapagliflozin", "dapagliflozin_3_o_glucuronide"]
    },
    # the urinary glucose of 0-24 hr of the pooled single and multiple dose
    # trials, whose time courses are Komoroski2009 and FDAMB102002
    "Gould2013": {
        "fm_dap5",
        "fm_dap10",
        "fm_dap20",
        "fm_dap50",
        "fm_dap100",
        "fm_dap250",
        "fm_dap500",
    },
    # the ABCB1 C1236T and G2677T/A groups of the 100 subjects of the C3435T
    # groups
    "Hwang2022a": {
        "fm_dap10_C1236T_CC",
        "fm_dap10_C1236T_CT",
        "fm_dap10_C1236T_TT",
        "fm_dap10_G2677TA_GG",
        "fm_dap10_G2677TA_GT",
        "fm_dap10_G2677TA_GA",
        "fm_dap10_G2677TA_TT",
        "fm_dap10_G2677TA_TA",
        "fm_dap10_G2677TA_AA",
    },
    # the fixed dose combination, the same subjects as the individual tablets
    "Khomitskaya2018": {"fm_dap10_DAP10__MET1000"},
    # the test formulation, the same subjects as the reference
    "Kim2023a": {"fm_dap10_DAP_FOR10"},
    # the multiple dose trial of the FDA review FDAMB102002 (days 1, 7 and 14)
    "Komoroski2009": {
        *(
            f"fm_po_dap{dose}_MD_{day}"
            for dose in ["2p5", "10", "20", "50", "100"]
            for day in ["single", "multi"]
        ),
        *(f"fm_po_dap_{dose}_uge_multi" for dose in ["2p5", "10", "20", "50", "100"]),
    },
    # dapagliflozin and the total radioactivity of the profiling samples, the
    # subjects of Fig 6
    "Obermeier2010": {"fm_dap50_dap_profiling", "fm_dap50_daptot_profiling"},
}

#: the validation data by question
VALIDATION_QUESTIONS: dict[str, NamedMappings] = {
    "interactions": INTERACTION_VALIDATION,
    "doses": DOSE_VALIDATION,
    "populations": POPULATION_VALIDATION,
    "limitations": LIMITATION_VALIDATION,
}

#: the validation data
VALIDATION_MAPPINGS: NamedMappings = {
    study: set().union(
        *(mappings.get(study, set()) for mappings in VALIDATION_QUESTIONS.values())
    )
    for mappings in VALIDATION_QUESTIONS.values()
    for study in mappings
}


def study(fit_mapping: FitMapping) -> str:
    """Get the study of a fit mapping, the id of its experiment."""
    return fit_mapping.reference.experiment.sid


def named(
    mappings: NamedMappings, fit_mapping_key: str, fit_mapping: FitMapping
) -> bool:
    """Check whether a fit mapping is one of the named mappings of its study."""
    return fit_mapping_key in mappings.get(study(fit_mapping), set())


def check_named(mappings: NamedMappings, fit_mappings: FitMappings) -> None:
    """Check that the named mappings are fit mappings of their studies.

    Raises:
        ValueError: for a study which is no experiment, or a key which is no fit
            mapping of its study.
    """
    for sid, keys in mappings.items():
        experiment = fit_mappings.experiments.get(sid)
        if experiment is None:
            raise ValueError(f"'{sid}' is no experiment of the fits.")
        if unknown := sorted(keys - set(experiment._fit_mappings)):
            raise ValueError(f"'{unknown}' are no fit mappings of '{sid}'.")


def filter_aggregates(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """No summaries of data which is fitted one by one, see `AGGREGATE_MAPPINGS`."""
    return not named(AGGREGATE_MAPPINGS, fit_mapping_key, fit_mapping)


def filter_selection(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """The training data and the named validation data.

    A mapping is training data without a coadministered drug in healthy
    subjects or patients with type 2 diabetes and a normal renal function;
    every other mapping enters only as named validation data.
    """
    if named(VALIDATION_MAPPINGS, fit_mapping_key, fit_mapping):
        return True
    metadata = mapping_metadata(fit_mapping)
    return (
        metadata.coadministration == Coadministration.NONE
        and metadata.health in TRAINING_HEALTH
    )


#: Only the pharmacokinetic data.
filter_pk = observable_filter(PK_OBSERVABLES)

#: Only the pharmacodynamic data.
filter_pd = observable_filter(PD_OBSERVABLES)


def mapping_collections(
    filters: MappingFilter | list[MappingFilter],
    validation: NamedMappings = VALIDATION_MAPPINGS,
) -> dict[str, list[FitMappingCollection]]:
    """Select the data of a fit from the fit mappings of the studies.

    The filters select the training data, the mappings which fail them are
    excluded, and so are the summaries of data which is fitted one by one
    (`AGGREGATE_MAPPINGS`) and the mappings which are neither training nor
    named validation data (`filter_selection`). The outliers are
    `OUTLIER_MAPPINGS`, the validation data is `VALIDATION_MAPPINGS`. Neither
    kind enters the cost; both are simulated and evaluated with the training
    data, so a report says how the fitted parameters describe data they were
    not fitted on.

    Args:
        filters: filters of the training data.
        validation: the keys of the validation data by study.

    Returns:
        The fit mapping collections of all kinds by experiment id.

    Raises:
        ValueError: if an outlier is a mapping of no study or of several.
    """
    fit_mappings = FitMappings(
        experiment_classes=EXPERIMENT_CLASSES,
        base_path=DAPAGLIFLOZIN_PATH,
        data_path=DATA_PATHS,
    )
    for mappings in [OUTLIER_MAPPINGS, AGGREGATE_MAPPINGS, validation]:
        check_named(mappings, fit_mappings)
    # sbmlsim selects the outliers by their keys alone, which must name one
    # mapping of all studies
    outliers = {key for keys in OUTLIER_MAPPINGS.values() for key in keys}
    for key in outliers:
        studies = [
            sid
            for sid, experiment in fit_mappings.experiments.items()
            if key in experiment._fit_mappings
        ]
        if len(studies) != 1:
            raise ValueError(f"The outlier '{key}' is a mapping of {studies}.")
    selection = fit_mappings.select(
        filters=[
            *(filters if isinstance(filters, list) else [filters]),
            filter_aggregates,
            filter_selection,
        ],
        outliers=outliers,
        validation=lambda key, fit_mapping: named(validation, key, fit_mapping),
        print_info=False,
    )
    # the counts per experiment and kind, not the single mappings
    selection.print(detail=False)
    # the validation mappings by the question they answer
    selected = {
        (sid, key)
        for sid, kinds in selection.kinds.items()
        for key, kind in kinds.items()
        if kind is MappingKind.VALIDATION
    }

    def pairs(mappings: NamedMappings) -> set[tuple[str, str]]:
        return {(sid, key) for sid, keys in mappings.items() for key in keys}

    console.print(
        "validation: "
        + ", ".join(
            f"{len(selected & pairs(mappings))} {question}"
            for question, mappings in VALIDATION_QUESTIONS.items()
        )
    )
    return selection.collections


def f_collections_pk() -> dict[str, list[FitMappingCollection]]:
    """Pharmacokinetic data."""
    return mapping_collections([filter_pk])


def f_collections_pd() -> dict[str, list[FitMappingCollection]]:
    """Pharmacodynamic data."""
    return mapping_collections([filter_pd])


if __name__ == "__main__":
    for f in [f_collections_pk, f_collections_pd]:
        console.rule(title=f.__name__, align="left", style="white")
        f()
