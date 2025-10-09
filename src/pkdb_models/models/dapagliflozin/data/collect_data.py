from pathlib import Path

from pkdb_models.models.data import collect_tsv_files

def collect_dapagliflozin_data():
    common_parent: Path = Path(__file__).parents[5]
    source_dir = common_parent / "pkdb_data" / "studies" / "dapagliflozin"
    target_dir = Path(__file__).parent / "dapagliflozin"

    collect_tsv_files(source_dir=source_dir, target_dir=target_dir)

    # collect Sha2015 (canagliflozin)
    def is_Sha2015(study_name) -> bool:
        return study_name == "Sha2015"

    collect_tsv_files(
        source_dir=common_parent / "pkdb_data" / "studies" / "canagliflozin",
        target_dir=Path(__file__).parent / "canagliflozin",
        filter_study=is_Sha2015,
    )

if __name__ == "__main__":
    collect_dapagliflozin_data()

