"""Module for creating latex table from computational models table."""

from pathlib import Path
import pandas as pd


def create_latex_table(df: pd.DataFrame) -> str:
    """Transform DataFrame into LaTeX table."""

    # Preamble with color definitions
    latex_header = r"""
\begin{landscape}
\begin{table}[H]
\centering
\tabcolsep=3pt
\renewcommand{\arraystretch}{1.3}
\tiny
\begin{threeparttable} 

\caption{\scriptsize{\textbf{Summary of published computational models for dapagliflozin.} 
Overview of published computational models including model type, software/platform, 
reproducibility criteria (open software, open model, open code, open data, 
reproducibility, FAIR, longterm storage), resources, clinical data sources, and model scope.}}
\label{tab:computational_models_overview}

\begin{tabularx}{\linewidth}{
    >{\raggedright\arraybackslash}p{1.6cm}  % Study
    >{\centering\arraybackslash}p{0.8cm}    % PubMed ID
    >{\raggedright\arraybackslash}p{1.0cm}  % Model Type
    >{\raggedright\arraybackslash}p{1.6cm}  % Platform/Software
    >{\centering\arraybackslash}p{0.8cm}    % Open Software
    >{\centering\arraybackslash}p{0.8cm}    % Open Model
    >{\centering\arraybackslash}p{0.8cm}    % Open Code
    >{\centering\arraybackslash}p{0.8cm}    % Open Data
    >{\centering\arraybackslash}p{0.8cm}    % Reproducibility
    >{\centering\arraybackslash}p{0.8cm}    % FAIR
    >{\centering\arraybackslash}p{0.8cm}    % Longterm Storage
    >{\centering\arraybackslash}p{2.3cm}    % Resources
    >{\centering\arraybackslash}p{0.7cm}    % Studies
    >{\raggedright\arraybackslash}p{3.0cm}  % Clinical Data Used
    >{\raggedright\arraybackslash}p{5.0cm}  % Scope
}
\toprule
""".strip('\n')

    # Column headers row
    column_names = df.columns.tolist()
    column_headers = " & ".join([f"\\textbf{{{col}}}" for col in column_names]) + r" \\"

    # LaTeX body
    latex_body = f"{column_headers}\n\\midrule\n"

    # Indices of columns that need conditional coloring
    color_columns = {
        'Open Software': 4,
        'Open Model': 5,
        'Open Code': 6,
        'Open Data': 7,
        'Reproducibility': 8,
        'FAIR': 9,
        'Longterm Storage': 10
    }

    for index, row in df.iterrows():
        values = list(row.astype(str).values)

        # First column: Study with citation
        values[0] = f"{values[0]} \\cite{{{values[0]}}}"

        # Second column: PubMed ID to clickable link
        if values[1] and values[1] != 'nan' and values[1] != '' and values[1] != '-':
            values[1] = f"\\href{{https://pubmed.ncbi.nlm.nih.gov/{values[1]}/}}{{{values[1]}}}"

        # Apply conditional coloring to Yes/No columns
        for col_name, col_idx in color_columns.items():
            if col_idx < len(values):
                val = values[col_idx].strip()
                if val == 'Yes':
                    values[col_idx] = f"\\cellcolor{{green!20}}Yes"
                elif val == 'No':
                    values[col_idx] = f"\\cellcolor{{red!20}}No"

        # Clean up nan values
        values = [v if v != 'nan' else '' for v in values]
        # Join row with proper escaping
        row_str = " & ".join(values) + r" \\"
        # Fix for long name
        row_str = row_str.replace("vanderAartvanderBeek2021 ",
                                  r"vanderAart\-vanderBeek2021 ")

        latex_body += row_str + "\n"
        latex_body += r"\addlinespace[1pt]" + "\n"

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

    return full_latex


if __name__ == "__main__":
    # Input and output file paths
    tsv_path = Path(__file__).parent / 'computational_models.tsv'
    latex_path = Path(__file__).parent / 'computational_models.tex'

    # Load TSV file
    df = pd.read_csv(tsv_path, sep="\t")
    df = df.fillna("")

    # Create LaTeX string
    latex_str = create_latex_table(df)

    # Save LaTeX table to file
    with open(latex_path, 'w') as f_tex:
        f_tex.write(latex_str)