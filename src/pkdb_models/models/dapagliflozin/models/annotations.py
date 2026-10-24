"""Annotations of the dapagliflozin models: the model and its species.

The compartments are the shared ones of `pkpdsim.annotations`.
"""

from pkpdsim.annotations import compartments
from sbmlutils.metadata import (
    BQB,
)

__all__ = ["compartments", "model", "species"]

model = [
    # taxonomy
    (BQB.HAS_TAXON, "taxonomy/9606"),  # human
    (BQB.HAS_TAXON, "snomedct/337915000"),  # human
    # modelling approach
    (BQB.HAS_PROPERTY, "mamo/MAMO_0000046"),  # ordinary differential equation model
    # biological process explained by model (GO/NCIT)
    (BQB.HAS_PROPERTY, "NCIT:C78126"),  # dapagliflozin
    (BQB.HAS_PROPERTY, "NCIT:C98083"),  # SGLT2 inhibitor
    (BQB.HAS_PROPERTY, "NCIT:C15299"),  # pharmacokinetics
    # model relevance to a particular area
    (BQB.HAS_PROPERTY, "NCIT:C2985"),  # Diabetes Mellitus
    (BQB.HAS_PROPERTY, "NCIT:C99532"),  # Diabetes Therapy
    # model reference
    (BQB.IS_DESCRIBED_BY, "doi/10.3390/pharmaceutics18030287"),
    (BQB.IS_DESCRIBED_BY, "doi/10.36903/physiome.31368556"),
    (BQB.IS_DESCRIBED_BY, "doi/10.5281/zenodo.13987865"),
]

species = {
    "glc": [  # glucose
        (BQB.IS, "CHEBI:4167"),
        (BQB.IS, "NCIT:C2831"),
    ],
    "dap": [  # dapagliflozin
        (BQB.IS, "CHEBI:85078"),
        (BQB.IS, "NCIT:C78126"),
    ],
    "d3g": [  # dapagliflozin-3-o-glucuronide
        (BQB.IS, "pubchem.compound/91617971"),
        (BQB.IS, "inchikey/ZYZULHSUKTZGTR-PTNNFGGUSA-N"),
    ],
    "dapmet": [  # other metabolites: 2-O-glucuronide, oxidative metabolites
    ],
    "daptot": [  # total dapagliflozin
    ],
}
