"""Shared triple oxygen isotope calculations used by the coral scripts."""

# Import libraries -------------------------------------------------------------
import numpy as np
import pandas as pd
from numpy.typing import ArrayLike


# Isotope calculations ---------------------------------------------------------

def prime(x: ArrayLike) -> np.ndarray | np.float64:
    """Convert delta notation to prime notation.

    Parameters
    ----------
    x : array-like or float
        Isotopic composition in per mil.

    Returns
    -------
    numpy.ndarray or float
        Prime-transformed isotope values in per mil.
    """
    x_arr = np.asarray(x, dtype=float)
    return 1000.0 * np.log1p(x_arr / 1000.0)


def unprime(x: ArrayLike) -> np.ndarray | np.float64:
    """Convert prime notation back to delta notation.

    Parameters
    ----------
    x : array-like or float
        Prime-transformed isotope values in per mil.

    Returns
    -------
    numpy.ndarray or float
        Isotopic composition in per mil.
    """
    x_arr = np.asarray(x, dtype=float)
    return np.expm1(x_arr / 1000.0) * 1000.0


def Dp17O(d17O: ArrayLike, d18O: ArrayLike) -> np.ndarray | np.float64:
    """Calculate ∆'17O from d17O and d18O values.

    Parameters
    ----------
    d17O : array-like or float
        d17O values in per mil.
    d18O : array-like or float
        d18O values in per mil.

    Returns
    -------
    numpy.ndarray or float
        ∆'17O values in ppm.
    """
    d17_prime = prime(d17O)
    d18_prime = prime(d18O)
    return (d17_prime - 0.528 * d18_prime) * 1000.0


def d17O(d18O: ArrayLike, Dp17O_value: ArrayLike) -> np.ndarray | np.float64:
    """Calculate d17O (per mil) from d18O (per mil) and ∆'17O (ppm)."""
    return unprime(Dp17O_value / 1000.0 + 0.528 * prime(d18O))


def calculate_theta(
    d18O_A: ArrayLike, Dp17O_A: ArrayLike, d18O_B: ArrayLike, Dp17O_B: ArrayLike,
) -> np.ndarray | np.float64:
    """Calculate effective theta between two isotope states A and B."""
    a18 = (np.asarray(d18O_B) + 1000.0) / (np.asarray(d18O_A) + 1000.0)
    a17 = (d17O(d18O_B, Dp17O_B) + 1000.0) / (d17O(d18O_A, Dp17O_A) + 1000.0)
    return np.round(np.log(a17) / np.log(a18), 4)


def apply_theta(
    d18O_A: ArrayLike, Dp17O_A: ArrayLike, d18O_B: ArrayLike | None = None,
    shift_d18O: ArrayLike | None = None, theta: ArrayLike | None = None,
) -> np.ndarray | np.float64:
    """Apply theta to predict ∆'17O at state B from state A."""
    if d18O_B is None:
        d18O_B = np.asarray(d18O_A) + shift_d18O

    d17O_A = d17O(d18O_A, Dp17O_A)
    a18 = (np.asarray(d18O_B) + 1000.0) / (np.asarray(d18O_A) + 1000.0)
    a17 = a18 ** theta

    d17O_B = a17 * (d17O_A + 1000.0) - 1000.0
    return Dp17O(d17O_B, d18O_B)


def mix_d17O(
    d18O_A: float, d17O_A: float | None = None, Dp17O_A: float | None = None,
    d18O_B: float | None = None, d17O_B: float | None = None,
    Dp17O_B: float | None = None, step: int = 100,
) -> pd.DataFrame:
    """Mix two reservoirs and return d17O, d18O, and ∆'17O trajectories."""
    # Preserve the original sampling grid and endpoint rounding.
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


def B_from_a(a: ArrayLike, A: ArrayLike) -> np.ndarray | np.float64:
    """Return reservoir B from fractionation factor ``a`` and reservoir A."""
    return (np.asarray(A) + 1000.0) / np.asarray(a) - 1000.0


def A_from_a(a: ArrayLike, B: ArrayLike) -> np.ndarray | np.float64:
    """Return reservoir A from fractionation factor ``a`` and reservoir B."""
    return (np.asarray(B) + 1000.0) * np.asarray(a) - 1000.0
