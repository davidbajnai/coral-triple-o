# This code creates a single-panel plot of coral disequilibrium data in
# DeltaDelta'17O vs Deltad18O space (Figure 1b style, points only).
#
# INPUT: SK_Table_S-3_part-2.csv
# OUTPUT: SK_extra.png

# Import libraries
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from SK_functions import *
# all necessary

# Set up directories
script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(script_dir, "..", "data")
plots_dir = os.path.join(script_dir, "..", "plots")

# Apply plot parameters
apply_mpl_style()



def a18_calcite_Coplen(T):
    # Coplen (2007)
    return np.exp((17.4 * 1000 / T - 28.6) / 1000)

def a18_aragonite_GuoZhou(T):
    # Guo and Zhou (2019) - aragonite
    x = 0.0974*(1000/T)**4 - 2.0916*(1000/T)**3 + 15.071*(1000/T)**2 - 25.782*(1000/T) +14.665
    return np.exp(x/1000)

def theta_aragonite_GuoZhou(T):
    # Guo and Zhou (2019) - aragonite
    return 78.1173 / T ** 2 - 1.5152 / T + 0.5299

def theta_aragonite_Wostbrock(T):
    # Wostbrock et al. (2020) – aragonite
    return -1.53 / T + 0.5305


def a17_cc(T_kelvin, a18_func, theta_func):
    return a18_func(T_kelvin) ** theta_func(T_kelvin)


def d18O_cc(T_celsius, d18Ow, a18_func):
    T_kelvin = T_celsius + 273.15
    return a18_func(T_kelvin) * (d18Ow + 1000) - 1000


def d17O_cc(T_celsius, d17Ow, a18_func, theta_func):
    T_kelvin = T_celsius + 273.15
    return a17_cc(T_kelvin, a18_func, theta_func) * (d17Ow + 1000) - 1000


def cc_equilibrium(T, T_err, d18Ow, d18Ow_err, Dp17Ow, Dp17Ow_err, a18_func, theta_func, Sample=None):
    df_equi = pd.DataFrame([])
    if Sample is not None:
        df_equi["SampleName"] = Sample

    df_equi["T"] = T
    df_equi["T_err"] = T_err
    df_equi["d18Ow"] = d18Ow
    df_equi["d18Ow_err"] = d18Ow_err
    df_equi["Dp17Ow"] = Dp17Ow
    df_equi["Dp17Ow_err"] = Dp17Ow_err

    df_equi["d18O_equi"] = d18O_cc(T, d18Ow, a18_func)
    df_equi["d17O_equi"] = d17O_cc(T, unprime((Dp17Ow / 1000) + 0.528 * prime(d18Ow)), a18_func, theta_func)
    df_equi["Dp17O_equi"] = Dp17O(df_equi["d17O_equi"], df_equi["d18O_equi"])

    df_equi["d18O_equi_min"] = d18O_cc(T + T_err, d18Ow - d18Ow_err, a18_func)
    df_equi["d17O_equi_min"] = d17O_cc(T + T_err, unprime((Dp17Ow / 1000) + 0.528 * prime(d18Ow - d18Ow_err)), a18_func, theta_func)
    df_equi["Dp17O_equi_max"] = Dp17O(df_equi["d17O_equi_min"], df_equi["d18O_equi_min"])

    df_equi["d18O_equi_max"] = d18O_cc(T - T_err, d18Ow + d18Ow_err, a18_func)
    df_equi["d17O_equi_max"] = d17O_cc(T - T_err, unprime((Dp17Ow / 1000) + 0.528 * prime(d18Ow + d18Ow_err)), a18_func, theta_func)
    df_equi["Dp17O_equi_min"] = Dp17O(df_equi["d17O_equi_max"], df_equi["d18O_equi_max"])

    df_equi["d18O_equi_err"] = (df_equi["d18O_equi_max"] - df_equi["d18O_equi_min"]) / 2
    df_equi["d17O_equi_err"] = (df_equi["d17O_equi_max"] - df_equi["d17O_equi_min"]) / 2
    df_equi["Dp17O_equi_err"] = np.sqrt(((df_equi["Dp17O_equi_max"] - df_equi["Dp17O_equi_min"]) / 2)**2 + Dp17Ow_err**2)

    return df_equi[["SampleName", "d18O_equi", "d18O_equi_err", "Dp17O_equi", "Dp17O_equi_err"]].rename(
        columns={
            "d18O_equi": "d18O_equilibrium",
            "d18O_equi_err": "d18O_equilibrium_err",
            "Dp17O_equi": "Dp17O_equilibrium",
            "Dp17O_equi_err": "Dp17O_equilibrium_err",
        }
    )

df = pd.read_csv(os.path.join(data_dir, "SK_Table_S-3_part-2.csv"), sep=",")

equation_sets = {
    "guozhou_aragonite_guozhou_aragonite": (a18_aragonite_GuoZhou, theta_aragonite_GuoZhou),
    "coplen_guozhou_aragonite": (a18_calcite_Coplen, theta_aragonite_GuoZhou),
    "daeron_hayles_calcite": (a18_calcite_Coplen, theta_aragonite_Wostbrock),
}

for label, (a18_func, theta_func) in equation_sets.items():
    df_equi = cc_equilibrium(
        T=df["T_database"],
        T_err=df["T_database_err"],
        d18Ow=df["d18Osw_database"],
        d18Ow_err=df["d18Osw_database_err"],
        Dp17Ow=df["Dp17Osw"],
        Dp17Ow_err=df["Dp17Osw_err"],
        a18_func=a18_func,
        theta_func=theta_func,
        Sample=df["SampleName"],
    )
    df[f"d18O_offset_{label}"] = df["d18O_AC"] - df_equi["d18O_equilibrium"].values
    df[f"Dp17O_offset_{label}"] = df["Dp17O_AC"] - df_equi["Dp17O_equilibrium"].values

fig, ax = plt.subplots(1, 1)

plot_colors = {
    "coplen_guozhou_aragonite": "white",
    "guozhou_aragonite_guozhou_aragonite": "black",
    "daeron_hayles_calcite": "blue",
}
plot_labels = {
    "coplen_guozhou_aragonite": r"$^{18}\alpha$: Coplen (2007) | $\theta$: Guo & Zhou (2019) - aragonite",
    "guozhou_aragonite_guozhou_aragonite": r"$^{18}\alpha$: Guo & Zhou (2019) - aragonite | $\theta$: Guo & Zhou (2019) - aragonite",
    "daeron_hayles_calcite": r"$^{18}\alpha$: Coplen (2007) | $\theta$: Wostbrock et al. (2020) - aragonite",
}

for label in equation_sets:
    x = df[f"d18O_offset_{label}"]
    y = df[f"Dp17O_offset_{label}"]
    ax.scatter(
        x,
        y,
        marker="o",
        facecolors=plot_colors[label],
        edgecolors="black",
        label=plot_labels[label],
    )

# Axis properties copied from Figure 1b
ax.set_xlim(-7.5, 0.5)
ax.set_ylim(-35, 25)
ax.set_ylabel("$\\Delta\\prime^{17}$O$_{measured}$ - $\\Delta\\prime^{17}$O$_{expected}$ (ppm)")
ax.set_xlabel("$\\delta^{18}$O$_{measured}$ - $\\delta^{18}$O$_{expected}$ (‰, VSMOW)")
ax.legend(bbox_to_anchor=(-0.1, -0.15), loc='upper left')

plt.savefig(os.path.join(plots_dir, "SK_extra"))
plt.close()