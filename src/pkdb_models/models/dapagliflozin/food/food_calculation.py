from pathlib import Path
from typing import Any

import pandas as pd
import numpy as np

from sbmlutils.console import console
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns

import matplotlib

SMALL_SIZE = 12
MEDIUM_SIZE = 18
BIGGER_SIZE = 25

matplotlib.rc('font', size=SMALL_SIZE)          # controls default text sizes
matplotlib.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
matplotlib.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
matplotlib.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
matplotlib.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
matplotlib.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
matplotlib.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
matplotlib.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


def food_effect(ratio: float) -> str:
    """Calculates the effect size from given ratio of pharmacokinetic parameters.

    E.g. mean(Cmax(fed))/mean(Cmax(fasted))
    """
    if ratio>=5.0:
        return "strong induction"
    elif 5.0 > ratio >= 2.0:
        return "moderate induction"
    elif 2.0 > ratio >= 1.25:
        return "weak induction"
    elif 1.25 > ratio >= 0.8:
        return "no effect"
    elif 0.8 > ratio >= 0.5:
        return "weak inhibition"
    elif 0.5 > ratio >= 0.2:
        return "moderate inhibition"
    elif 0.2 > ratio:
        return "strong inhibition"
    if ratio < 0.0:
        raise ValueError(f"Ratio of pharmacokinetic parameters must be > 0.0, but {ratio=}.")
    else:
        raise ValueError(f"Effect cannot be calculated for {ratio=}.")


def significance_for_pvalue(pvalue: float) -> str:
    """Significance marker for pvalue."""
    if pvalue <= 0.001:
        significance = "***"
    elif pvalue <= 0.01:
        significance = "**"
    elif pvalue <= 0.05:
        significance = "*"
    else:
        significance = " "
    return significance


def food_significance(mean_fasted: float, sd_fasted: float, n_fasted: int, mean_fed: float, sd_fed: float, n_fed: int) -> dict[str, Any]:
    """Calculate significance of effect from mean and standard deviation.

    Works on iterables of values.

    p-values: * = 0.05, ** = 0.01, *** = 0.001
    """
    mean1, sd1, n1 = mean_fasted, sd_fasted, n_fasted
    mean2, sd2, n2 = mean_fed, sd_fed, n_fed

    # Calculate the pooled standard deviation
    pooled_sd = np.sqrt(((n1 - 1) * sd1 ** 2 + (n2 - 1) * sd2 ** 2) / (n1 + n2 - 2))

    # Calculate the t-statistic
    t_statistic = (mean1 - mean2) / (pooled_sd * np.sqrt(1 / n1 + 1 / n2))
    # Calculate degrees of freedom
    df = n1 + n2 - 2

    # Calculate the p-value
    p_value = stats.t.sf(np.abs(t_statistic), df) * 2  # two-tailed test
    significance = [significance_for_pvalue(p) for p in p_value]

    # Calculate effect size
    ratio = mean_fed / mean_fasted
    # Ratio SD (error propagation), worst case approximation
    ratio_sd = ratio * np.sqrt(np.power((sd_fasted / mean_fasted), 2) + np.power((sd_fed / mean_fed), 2))
    effect = [food_effect(ratio=r) for r in ratio]

    return {
        "ratio": ratio,
        "ratio_sd": ratio_sd,
        "effect": effect,
        "p-value": p_value,
        "sig": significance,
        "t-statistic": t_statistic,
        "df": df,
        "method": "t-test two-sided",
    }


def load_food_data(tsv_path: Path, output_path: Path) -> dict[str, pd.DataFrame]:
    """Load the DDI data for Cmax and AUC.

    Data managed in google sheets: https://docs.google.com/spreadsheets/d/1ixcQzusK2JYn9Uz8xerJpZ6qr3HBGOV9nbrMDbUuYQA/edit?gid=790035504#gid=790035504
    """
    df_all = pd.read_csv(tsv_path, sep="\t", comment="#")

    keys = ["cmax", "auc_inf", "auc_end"]
    data: dict[str, pd.DataFrame] = {}
    for key in keys:
        df = df_all[df_all.parameter == key].copy()

        raw_data = dict(
            mean_fasted=df[f"mean_fasted"],
            sd_fasted=df[f"sd_fasted"],
            n_fasted=df[f"n_fasted"],
            mean_fed=df[f"mean_fed"],
            sd_fed=df[f"sd_fed"],
            n_fed=df[f"n_fed"]
        )
        food_ttest = food_significance(**raw_data)
        food_ttest = {
            "study": df["study"],
            "intervention": df["intervention"],
            "dose": df["dose"],
            "dose_unit": df["dose_unit"],
            **raw_data,
            **food_ttest,
        }

        #sorting by values

        df_key = pd.DataFrame(food_ttest)
        #df_key_sorted=df_key.sort_values(by = "ratio")
        data[key] = df_key

        # save
        df_key.to_csv(output_path / f"food_effect_{key}.tsv", sep="\t", index=False)
        df_key.to_latex(output_path / f"food_effect_{key}.tex")

    for key, df in data.items():
        console.rule(f"{key.upper()} Food effect", style="white")
        console.print(df.to_string())

    return data


def plot_food(data: dict[str, pd.DataFrame], output_dir: Path) -> None:
    """Plots the food effect for dapagliflozin."""

    spectral_color = sns.color_palette("Spectral_r", 12, as_cmap=False)

    for key, df_raw in data.items():
        df = df_raw[::-1].reset_index(drop=True)
        #df = df.sort_values(by = "ratio")

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(9, 2), dpi=300, layout="constrained")
        alpha = 0.25
        ax.set_xscale("log")
        ax.axvspan(0.2, 0.5, facecolor=spectral_color[9], alpha=alpha)
        ax.axvspan(0.5, 0.8, facecolor=spectral_color[7], alpha=alpha)
        ax.axvspan(0.8, 1.25, facecolor="w", alpha=alpha)
        ax.axvspan(1.25, 2, facecolor=spectral_color[5], alpha=alpha)
        ax.axvspan(2, 5, facecolor=spectral_color[3], alpha=alpha)

        ax.set_xlim([0.1, 10])
        left, right = ax.get_xlim()
        ax.axvspan(0.2, left, facecolor=spectral_color[11], alpha=alpha)
        ax.axvspan(5, right, facecolor=spectral_color[0], alpha=alpha)

        positions = [0.2, 0.5, 0.8, 1.25, 2, 5]

        # ax2.set_xlim(right=right + 3)
        ax.set_xticks(positions)
        ax.set_xticklabels(positions)

        # set ticks
        ax2 = ax.twiny()
        ax2.set_xscale("log")
        ax2.set_xlim(ax.get_xlim())
        label_positions = [0.1,] + positions + [10.0]
        label_positions = [ (label_positions[k] + label_positions[k+1])/2 for k in range( len(label_positions)-1)]
        ax2.set_xticks(label_positions)
        ax2.set_xticklabels(("strong", "moderate", "weak", "no effect", "weak", "moderate", "strong"))

        # ax.set_title(f"Food effect", fontdict={"weight": "bold"})
        ax.axvline(x=1.0, linestyle="--", color="black")
        ax.errorbar(
            x=df.ratio,
            y=df.study,
            xerr=df.ratio_sd,
            color="black",
            ecolor="darkgrey",
            markeredgecolor="white",
            marker="s",
            # markersize=10,
            linestyle="",
        )
        # significance stars
        k = 0
        for idx, row in df.iterrows():
            significance = row["sig"]
            y = k
            if significance:
                xy = (row["ratio"], y)
                ax.annotate(text=significance, xy=xy)
            k = k + 1

        ax.set_xlabel(f"{key} ratio [-]", fontdict={"weight": "bold"})
        plt.show()

        fig.savefig(output_dir / f"ddi_{key}.png", bbox_inches="tight")


if __name__ == "__main__":
    from pkdb_models.models.dapagliflozin import RESULTS_PATH_SIMULATION

    food_effect_path = RESULTS_PATH_SIMULATION / "food_effect"
    food_effect_path.mkdir(exist_ok=True)

    tsv_path = Path(__file__).parent / "dapagliflozin_studies - food.tsv"
    data: dict[str, pd.DataFrame] = load_food_data(tsv_path, output_path=food_effect_path)
    plot_food(data, output_dir=food_effect_path)

