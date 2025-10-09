"""Module for creating latex table from study table."""

from pathlib import Path
import pandas as pd
from sbmlutils.console import console


def create_latex_table(df: pd.DataFrame) -> str:
    """Transform DataFrame into LaTeX table."""

    #   p{2.6cm}  % Study
    #   p{1.7cm}  % PK-DB ID
    #   p{1.5cm}  % Substance
    #   p{0.9cm}  % Route
    #   p{1.0cm}  % Dosing
    #   p{0.7cm}  % Dose [mg]


    # Preamble:
    latex_header = r"""
\begin{landscape}
\begin{table}[H]
\centering
\tabcolsep=3.5pt\relax
\scriptsize
\begin{threeparttable} 

\caption{\scriptsize{\textbf{Summary of studies for modeling.} 
Overview of study identifiers, PK-DB IDs, administered substance, route, dosing, 
and subject characteristics, including health status (\emph{H}), renal impairment (\emph{RI}), 
hepatic impairment (\emph{HI}), fasting status and urinary glucose excretion (\emph{UGE}) and renal
treshold for glucose (\emph{RTG}). \emph{DAP P} = dapagliflozin plasma, \emph{DAP U} = dapagliflozin urine, 
\emph{DAP F} = dapagliflozin feces, \emph{D3G P} = dapagliflozin-3 glucuronide plasma, 
\emph{D3G U} = dapagliflozin-3 glucuronide urine}},
\label{table:curated_data_overview}

\begin{tabularx}{\textwidth}{
    p{2.6cm}  % Study
    p{1.7cm}  % PK-DB ID
    p{1.5cm}  % Substance
    p{0.9cm}  % Route
    p{1.0cm}  % Dosing
    c % Dose [mg]
    c % dap plasma,
    c % dap urine,
    c % dap feces,
    c % d3g plasma,
    c % d3g urine,
    p{0.4cm}  % H
    p{0.4cm}  % HI
    p{0.4cm}  % RI
    p{0.4cm}  % T1
    p{0.4cm}  % T2
    p{0.4cm}  % Fed
    p{0.4cm}  % Fast
    p{0.6cm}  % UGE
    p{0.6cm}  % RTG
}
\toprule
""".strip('\n')

    # Rename columns
    df = df.rename(columns={
        "study": "Study",
        "pkdb": "PK-DB",
        "substance": "Substance",
        "route": "Route",
        "dosing": "Dosing",
        "dose": "Dose [mg]",

        "healthy": "H",
        "renal impairment": "RI",
        "hepatic impairment": "HI",
        "t1dm": "T1",
        "t2dm": "T2",

        "dap plasma": "DAP P",
        "dap urine": "DAP U",
        "dap feces": "DAP F",
        "d3g plasma": "D3G P",
        "d3g urine": "D3G U",

        "fed": "Fed",
        "fasted": "Fast",

    })
    columns_bool = [
        "H", "RI", "HI", "T1", "T2", "Fed", "Fast", "UGE", "RTG",
        "DAP P",
        "DAP U",
        "DAP F",
        "D3G P",
        "D3G U",
    ]
    # Convert True/False to checkmarks
    for col in columns_bool:
        df[col] = df[col].apply(lambda x: r"\checkmark" if str(x).strip().upper() in "TRUE" else "")


    # console.print(df)

    # Column headers row
    column_names = df.columns.tolist()
    column_headers = " & ".join([f"\\textbf{{{col}}}" for col in column_names]) + r" \\"

    # Genotypes columns "*" -> "\textasteriskcentered"
    # if "Genotype" in df.columns:
    #     df["Genotype"] = df["Genotype"].apply(
    #         lambda x: x.replace("*", r"\textasteriskcentered") if isinstance(x, str) else x
    #     )

    # LaTeX body
    latex_body = f"{column_headers}\n\\midrule\n"

    for _, row in df.iterrows():
        values = list(row.astype(str).values)
        values[0] = f"{values[0]} \\cite{{{values[0]}}}"

        # PK-DB ID to clickable link
        values[1] = f"\\href{{https://identifiers.org/pkdb:{values[1]}}}{{{values[1]}}}"

        # Converting leftover "TRUE" or "-" strings
        row_str = " & ".join(v.replace("TRUE", r"\checkmark").replace("-", "")
                             for v in values) + r" \\"

        # fix for long name
        row_str = row_str.replace("vanderAartvanderBeek2020 ",
                                       r"vanderAart\-vanderBeek2020 ")

        latex_body += row_str + "\n"

    # End of the tabular portion
    latex_footer = r"""
\bottomrule
\end{tabularx}

\end{threeparttable} 
\end{table}
\end{landscape}
""".strip('\n')

    # Combine all parts
    full_latex = latex_header + "\n" + latex_body + latex_footer

    console.rule(style="white")
    console.print(df)
    console.rule(style="white")

    return full_latex


if __name__ == "__main__":
    # Input and output file paths
    tsv_path = Path(__file__).parent / 'dapagliflozin_studies.tsv'
    latex_path = Path(__file__).parent / 'dapagliflozin_studies.tex'

    # Load TSV file
    columns = [
        "study",
        "pkdb",
        # "pmid",
        "substance",
        "route",
        "dosing",
        "dose",

        "healthy",
        "renal impairment",
        "hepatic impairment",
        "t1dm",
        "t2dm",
        "dap plasma",
        "dap urine",
        "dap feces",
        "d3g plasma",
        "d3g urine",
        "fed",
        "fasted",
        "UGE",
        "RTG",
    ]
    df = pd.read_csv(
        tsv_path, sep="\t", skiprows=1, skipfooter=27,
        usecols=columns,
        engine="python",
    )[columns]

    # Replace NaN values with an empty string
    df = df.fillna("")

    # # Create LaTeX string
    latex_str = create_latex_table(df)

    # Save LaTeX table to file
    with open(latex_path, 'w') as f_tex:
        f_tex.write(latex_str)
