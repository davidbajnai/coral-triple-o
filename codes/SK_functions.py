"""Core utilities for SK triple oxygen isotope calculations.

Notation follows the manuscript convention for d17O, d18O, and ∆'17O.
Code identifiers use ``Dp17O`` for the capital-delta quantity (in ppm).
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(script_dir, "..", "data")
plots_dir = os.path.join(script_dir, "..", "plots")


def apply_mpl_style():
    """Apply the common matplotlib style used across SK figures."""
    plt.rcParams["legend.loc"] = "best"
    plt.rcParams.update({'font.size': 7})
    plt.rcParams['scatter.edgecolors'] = "k"
    plt.rcParams['scatter.marker'] = "o"
    plt.rcParams["lines.linewidth"] = 0.5
    plt.rcParams["patch.linewidth"] = 0.5
    plt.rcParams["figure.figsize"] = (4, 4)
    plt.rcParams["savefig.dpi"] = 600
    plt.rcParams["savefig.bbox"] = "tight"
    plt.rcParams['savefig.transparent'] = False
    plt.rcParams['mathtext.default'] = 'regular'


def prime(x):
    """Return d' notation from d values (per mil), vectorized and stable."""
    x = np.asarray(x, dtype=float)
    return 1000.0 * np.log1p(x / 1000.0)


def unprime(x):
    """Invert d' notation back to d values (per mil), vectorized and stable."""
    x = np.asarray(x, dtype=float)
    return 1000.0 * np.expm1(x / 1000.0)


def Dp17O(d17O, d18O):
    """Calculate ∆'17O (ppm) from d17O and d18O (per mil)."""
    return (prime(d17O) - 0.528 * prime(d18O)) * 1000.0


def d17O(d18O, Dp17O_value):
    """Calculate d17O (per mil) from d18O (per mil) and ∆'17O (ppm)."""
    return unprime(Dp17O_value / 1000.0 + 0.528 * prime(d18O))


def calculate_theta(d18O_A, Dp17O_A, d18O_B, Dp17O_B):
    """Calculate effective theta between two isotope states A and B."""
    a18 = (np.asarray(d18O_B) + 1000.0) / (np.asarray(d18O_A) + 1000.0)
    a17 = (d17O(d18O_B, Dp17O_B) + 1000.0) / (d17O(d18O_A, Dp17O_A) + 1000.0)
    return np.round(np.log(a17) / np.log(a18), 4)


def apply_theta(d18O_A, Dp17O_A, d18O_B=None, shift_d18O=None, theta=None):
    """Apply theta to predict ∆'17O at state B from state A."""
    if d18O_B is None:
        d18O_B = np.asarray(d18O_A) + shift_d18O

    d17O_A = d17O(d18O_A, Dp17O_A)
    a18 = (np.asarray(d18O_B) + 1000.0) / (np.asarray(d18O_A) + 1000.0)
    a17 = a18 ** theta

    d17O_B = a17 * (d17O_A + 1000.0) - 1000.0
    return Dp17O(d17O_B, d18O_B)


def mix_d17O(d18O_A, d17O_A=None, Dp17O_A=None, d18O_B=None, d17O_B=None, Dp17O_B=None, step=100):
    """Mix two reservoirs and return d17O, d18O, and ∆'17O trajectories."""
    ratio_B = np.arange(0.0, 1.0 + 1.0 / step, 1.0 / step)

    if d17O_A is None:
        d17O_A = d17O(d18O_A, Dp17O_A)
    if d17O_B is None:
        d17O_B = d17O(d18O_B, Dp17O_B)

    mix_d18O = ratio_B * float(d18O_B) + (1.0 - ratio_B) * float(d18O_A)
    mix_d17O_values = ratio_B * float(d17O_B) + (1.0 - ratio_B) * float(d17O_A)
    mix_Dp17O = Dp17O(mix_d17O_values, mix_d18O)

    return pd.DataFrame(
        {
            "mix_d17O": mix_d17O_values,
            "mix_d18O": mix_d18O,
            "mix_Dp17O": mix_Dp17O,
            "xB": ratio_B * 100.0,
        }
    )


def B_from_a(a, A):
    """Return reservoir B from fractionation factor ``a`` and reservoir A."""
    return (np.asarray(A) + 1000.0) / np.asarray(a) - 1000.0


def A_from_a(a, B):
    """Return reservoir A from fractionation factor ``a`` and reservoir B."""
    return (np.asarray(B) + 1000.0) * np.asarray(a) - 1000.0


def epsilon(d18O_A, d18O_B):
    """Return epsilon (per mil) between two d18O reservoirs."""
    return ((np.asarray(d18O_A) + 1000.0) / (np.asarray(d18O_B) + 1000.0) - 1.0) * 1000.0
