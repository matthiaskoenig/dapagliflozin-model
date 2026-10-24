"""Collect the data of the dapagliflozin studies from the PK-DB study folders.

The studies are copied from a `pkdb_data` checkout (`PKDB_DATA`, else
`../pkdb_data`, see `pkpdsim.data.studies_dir`) into the data folder of the
model, see `pkpdsim.data.collect_tsv_files`.
Besides the dapagliflozin studies, one study of another substance has
dapagliflozin groups: Sha2015 (canagliflozin).
"""

from pathlib import Path

from pkpdsim.data import collect_tsv_files
from pkpdsim.data import studies_dir as pkdb_studies_dir

from pkdb_models.models.dapagliflozin import (
    DATA_PATH_CANAGLIFLOZIN,
    DATA_PATH_DAPAGLIFLOZIN,
)

#: the studies of other substances which the dapagliflozin experiments use
CANAGLIFLOZIN_STUDIES = {"Sha2015"}


def collect_data(studies_dir: Path | None = None) -> None:
    """Copy the data of the dapagliflozin studies into the data folder of the model.

    Args:
        studies_dir: the `studies` folder of pkdb_data, by default the one of
            the checkout `pkpdsim.data.studies_dir` finds.
    """
    if studies_dir is None:
        studies_dir = pkdb_studies_dir()
    collect_tsv_files(
        source_dir=studies_dir / "dapagliflozin", target_dir=DATA_PATH_DAPAGLIFLOZIN
    )
    collect_tsv_files(
        source_dir=studies_dir / "canagliflozin",
        target_dir=DATA_PATH_CANAGLIFLOZIN,
        filter_study=lambda study: study in CANAGLIFLOZIN_STUDIES,
    )


if __name__ == "__main__":
    collect_data()
