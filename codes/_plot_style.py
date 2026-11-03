"""Shared plot settings and PNG export utilities."""

# Import libraries -------------------------------------------------------------

from pathlib import Path

import matplotlib.pyplot as plt

# Plot settings ----------------------------------------------------------------

STYLE = {
    "legend.loc": "best",
    "font.size": 7,
    "scatter.edgecolors": "k",
    "scatter.marker": "o",
    "lines.linewidth": 0.5,
    "patch.linewidth": 0.5,
    "savefig.dpi": 600,
    "savefig.transparent": False,
    "savefig.bbox": "tight",
    "mathtext.default": "regular",
    "figure.figsize": (9, 4),
}

C = {
    "figure_bg": "#FFFFFF",
    "axes_bg": "#FFFFFF",

    "blue": "#1455C0",
    "blue1": "#E0EFFB",
    "blue2": "#B4D5F6",
    "blue3": "#73AEF4",
    "blue4": "#347DE0",
    "blue5": "#1455C0",
    "blue6": "#0C3992",
    "blue7": "#0A1E6E",
    "blue8": "#061350",

    "green": "#408335",
    "green1": "#E2F3E5",
    "green2": "#BDDBB9",
    "green3": "#8CBC80",
    "green4": "#66A558",
    "green5": "#408335",
    "green6": "#2A7230",
    "green7": "#165C27",
    "green8": "#154A26",

    "lgreen": "#63A615",
    "lgreen1": "#EBF7DD",
    "lgreen2": "#C9EB9E",
    "lgreen3": "#9FD45F",
    "lgreen4": "#78BE14",
    "lgreen5": "#63A615",
    "lgreen6": "#508B1B",
    "lgreen7": "#44741A",
    "lgreen8": "#375F15",

    "red": "#EC0016",
    "red1": "#FEE6E6",
    "red2": "#FCC8C3",
    "red3": "#FA9090",
    "red4": "#F75056",
    "red5": "#EC0016",
    "red6": "#C50014",
    "red7": "#9B000E",
    "red8": "#740009",

    "yellow": "#FFD800",
    "orange": "#FF7A00",
    "yellow1": "#FFFFDC",
    "yellow2": "#FFFFAF",
    "yellow3": "#FFF876",
    "yellow4": "#FFF000",
    "yellow5": "#FFD800",
    "yellow6": "#FFBB00",
    "yellow7": "#FF9B00",
    "yellow8": "#FF7A00",

    "gray": "#7D7D7D",
    "black": "#000000",
    "white": "#FFFFFF",

    "wgray": "#858379",
    "wgray1": "#F5F4F1",
    "wgray2": "#DDDED6",
    "wgray3": "#BCBBB2",
    "wgray4": "#9C9A8E",
    "wgray5": "#858379",
    "wgray6": "#747067",
    "wgray7": "#4F4B41",
    "wgray8": "#38342F",

    "cgray": "#646973",
    "cgray1": "#F0F3F5",
    "cgray2": "#D7DCE1",
    "cgray3": "#AFB4BB",
    "cgray4": "#878C96",
    "cgray5": "#646973",
    "cgray6": "#3C414B",
    "cgray7": "#282D37",
    "cgray8": "#131821",

    "errorbar": "#cacaca",
    "seawater": "#83CACA",
    "standard_blue": "blue",
    "model_warm": "darkred",
    "model_cold": "darkblue",
    "model_warm_alternative": "hotpink",
    "model_cold_alternative": "cornflowerblue",
}

# Plot utilities ---------------------------------------------------------------

def apply_plot_style() -> None:
    """Apply the original coral figure settings without changing other defaults."""
    plt.rcParams.update(STYLE)


def save_figure(
    fig: plt.Figure, figures_dir: str | Path, filename: str,
    fixed_layout: bool = False, *, tight_layout: bool = False,
) -> None:
    """Optionally adjust spacing, save the PNG, and close only this figure.

    Tight layout adjusts subplot spacing; tight bounding-box cropping is a
    separate export setting. A fixed layout skips adjustment and cropping.
    """
    if tight_layout and not fixed_layout:
        fig.tight_layout()
    figure_path = Path(figures_dir) / f"{filename}.png"
    if fixed_layout:
        with plt.rc_context({"savefig.bbox": None}):
            fig.savefig(figure_path)
    else:
        fig.savefig(figure_path)
    plt.close(fig)
