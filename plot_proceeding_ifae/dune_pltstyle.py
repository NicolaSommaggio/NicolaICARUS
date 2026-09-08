import matplotlib.pyplot as plt
from cycler import cycler

plt.rcParams.update({
    # FONT & TEXT
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "axes.labelsize": 20,
    "axes.titlesize": 24,
    "legend.fontsize": 16,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,

    # COLORS
    "axes.prop_cycle": cycler(
        "color",
        [
            "#cc4d00",
            "#e69900",
            "#009980",
            "#a87fff",
            "#ffb3d9",
            "#1f83e7",
            "#66ffcc",
            "#adb5bd",
            "#1a3a5f",
        ],
    ),
        #"#f2e640",

    # AXES & GRIDS
    "axes.grid": False,
    "axes.facecolor": "white",
    "axes.edgecolor": "black",
    "axes.linewidth": 2.5,

    # TICKS
    "xtick.direction": "in",
    "xtick.major.size": 8,
    "xtick.top": True,
    "ytick.direction": "in",
    "ytick.major.size": 8,
    "ytick.right": True,

    # LINES & PATCHES
    "lines.linewidth": 2.5,
    "patch.linewidth": 2.5,
    "patch.antialiased": True,

    # FIGURE & SAVING
    "figure.facecolor": "white",
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.1,
    #"savefig.dpi": 300,
})

print(plt.rcParams["axes.prop_cycle"].by_key()["color"])