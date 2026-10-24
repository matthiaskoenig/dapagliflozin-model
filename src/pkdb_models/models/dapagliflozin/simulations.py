"""The simulation experiments of the dapagliflozin model and their groups.

They are the registry of `MODEL`, the experiments its commands simulate.
"""

from sbmlsim.experiment import SimulationExperiment

from pkdb_models.models.dapagliflozin import MODEL
from pkdb_models.models.dapagliflozin.experiments.misc import (
    DoseDependencyExperiment,
    FoodEffect,
    HepaticRenalImpairment,
)
from pkdb_models.models.dapagliflozin.experiments.scans import (
    DapagliflozinParameterScan,
)
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

EXPERIMENTS: dict[str, list[type[SimulationExperiment]]] = {
    "studies": [
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
        Sha2015,  # the dapagliflozin groups of the canagliflozin study Sha2015
        Shah2019a,
        vanderAartvanderBeek2020,
        Watada2019,
        Yang2013,
    ],
    "pharmacodynamics": [
        FDAMB102002,
        FDAMB102003,
        FDAMB102007,
        Gould2013,
        Kasichayanula2011a,
        Kim2023,
        Komoroski2009,
        Sha2015,
        Watada2019,
        Yang2013,
    ],
    "dose_dependency": [
        FDAMB102002,
        FDAMB102003,
        Gould2013,
        Kasichayanula2011a,
        Komoroski2009,
        Watada2019,
        Yang2013,
    ],
    "food": [
        Kasichayanula2011b,
        Komoroski2009,
        LaCreta2016,
        Shah2019a,
    ],
    "hepatic_impairment": [
        Kasichayanula2011,
    ],
    "renal_impairment": [
        Kasichayanula2013,
        FDAMB102007,
    ],
    "misc": [
        DoseDependencyExperiment,
        FoodEffect,
        HepaticRenalImpairment,
    ],
    "scans": [
        DapagliflozinParameterScan,
    ],
}
EXPERIMENTS["all"] = EXPERIMENTS["studies"] + EXPERIMENTS["misc"] + EXPERIMENTS["scans"]


if __name__ == "__main__":
    MODEL.run_simulation_experiments(EXPERIMENTS["all"], MODEL.simulation_dir / "all")
