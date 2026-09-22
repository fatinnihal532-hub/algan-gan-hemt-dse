"""One consistent look for every figure, in a light and a dark variant,
matching the style used across the author's other simulation repos."""
import matplotlib as mpl

LIGHT = dict(surface="#fcfcfb", ink="#0b0b0b", secondary="#52514e",
             grid="#e1e0d9", axis="#c3c2b7",
             series=["#2a78d6", "#eb6834", "#2f9e57", "#a24fc9"])
DARK = dict(surface="#1a1a19", ink="#ffffff", secondary="#c3c2b7",
            grid="#2c2c2a", axis="#383835",
            series=["#3987e5", "#d95926", "#3fb86a", "#b968de"])


def apply(theme: dict):
    mpl.rcParams.update({
        "figure.facecolor": theme["surface"],
        "axes.facecolor": theme["surface"],
        "savefig.facecolor": theme["surface"],
        "text.color": theme["ink"],
        "axes.labelcolor": theme["ink"],
        "axes.edgecolor": theme["axis"],
        "xtick.color": theme["secondary"],
        "ytick.color": theme["secondary"],
        "grid.color": theme["grid"],
        "axes.grid": True,
        "grid.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 10,
        "font.family": "DejaVu Sans",
    })


def save_light_dark(fig_fn, path_stem, results_dir="results"):
    """Call fig_fn(theme) to draw onto the current figure/axes, once per
    theme, saving <path_stem>_light.svg and <path_stem>_dark.svg.

    Output is byte-for-byte reproducible: the SVG timestamp is dropped and
    element ids use a fixed hash salt, so CI can regenerate results/ and
    commit it without producing a spurious diff on every run."""
    import matplotlib.pyplot as plt
    mpl.rcParams["svg.hashsalt"] = "reproducible"
    for name, theme in [("light", LIGHT), ("dark", DARK)]:
        apply(theme)
        fig = fig_fn(theme)
        fig.savefig(f"{results_dir}/{path_stem}_{name}.svg", format="svg",
                    metadata={"Date": None})
        plt.close(fig)
