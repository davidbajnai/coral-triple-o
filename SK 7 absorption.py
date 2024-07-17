# This code is used to calcualte and plot the CO2 hydration and hydroxylation slopes

# INPUT: NONE
# OUTPUT: SK Figure S7.png

# >>>>>>>>>

# Import libraries
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Import functions
from functions import *


# Plot parameters
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

# Functions that make life easier

def a18_cc(T):

    # Used for the discussion
    return 0.0201 * (1000 / T) + 0.9642                     # Guo and Zhou (2019) – aragonite

    # Alternative equations: 

    # Hayles et al. (2018) - calcite
    # B_calcite = 7.027321E+14 / T**7 + -1.633009E+13 / T**6 + 1.463936E+11 / T**5 + -5.417531E+08 / T**4 + -4.495755E+05 / T**3  + 1.307870E+04 / T**2 + -5.393675E-01 / T + 1.331245E-04
    # B_water = -6.705843E+15 / T**7 + 1.333519E+14 / T**6 + -1.114055E+12 / T**5 + 5.090782E+09 / T**4 + -1.353889E+07 / T**3 + 2.143196E+04 / T**2 + 5.689300 / T + -7.839005E-03
    # return np.exp(B_calcite) / np.exp(B_water)
    
    # return np.exp((2.84 * 10**6 / T**2 - 2.96) / 1000)    # Wostbrock et al. (2020) – calcite
    # return np.exp((17.88 * 1000 / T - 31.14) / 1000)      # Kim et al. (2007) – aragonite
    # return np.exp((17.57 * 1000 / T - 29.13) / 1000)      # Daeron et al. (2019) – calcite


def theta_cc(T):

    # Used for the discussion
    return 59.1047/T**2 + -1.4089/T + 0.5297                # Guo and Zhou (2019) – aragonite

    # Alternative equations:

    # Hayles et al. (2018) - calcite
    # K_calcite = 1.019124E+09 / T**5 + -2.117501E+07 / T**4 + 1.686453E+05 / T**3 + -5.784679E+02 / T**2 + 1.489666E-01 / T + 0.5304852
    # B_calcite = 7.027321E+14 / T**7 + -1.633009E+13 / T**6 + 1.463936E+11 / T**5 + -5.417531E+08 / T**4 + -4.495755E+05 / T**3  + 1.307870E+04 / T**2 + -5.393675E-01 / T + 1.331245E-04
    # K_water = 7.625734E+06 / T**5 + 1.216102E+06 / T**4 + -2.135774E+04 / T**3 + 1.323782E+02 / T**2 + -4.931630E-01 / T + 0.5306551
    # B_water = -6.705843E+15 / T**7 + 1.333519E+14 / T**6 + -1.114055E+12 / T**5 + 5.090782E+09 / T**4 + -1.353889E+07 / T**3 + 2.143196E+04 / T**2 + 5.689300 / T + -7.839005E-03
    # a18 = np.exp(B_calcite) / np.exp(B_water)
    # return K_calcite + (K_calcite-K_water) * (B_water / np.log(a18))

    # return -1.39 / T + 0.5305                             # Wostbrock et al. (2020) – calcite
    # return -1.53 / T + 0.5305                             # Wostbrock et al. (2020) – aragonite


def a17_cc(T):
    return a18_cc(T)**theta_cc(T)


def d18O_cc(T, d18Ow):
    return a18_cc(T) * (d18Ow+1000) - 1000


def d17O_cc(T, d17Ow):
    return a17_cc(T) * (d17Ow+1000) - 1000


def a18OH(T=273.15+22, eq="Z20-X3LYP"):
    if (eq == "Z20-X3LYP"):
        e18_H2O_OH = (-4.4573 + (10.3255 * 10**3) /
                      (T) + (-0.5976 * 10**6) / (T)**2)
    elif (eq == "Z20-MP2"):
        e18_H2O_OH = (-4.0771 + (9.8350 * 10**3) /
                      (T) + (-0.8729 * 10**6) / (T)**2)
    elif (eq == "BH21_original"):
        e18_H2O_OH = -0.034 * (T-273.15) + 43.4
    elif (eq == "BH21"):
        e18_H2O_OH = -0.035 * (T-273.15) + 40.1

    return e18_H2O_OH / 1000 + 1


def a17OH(T=273.15+22, eq="Z20-X3LYP", theta=0.529):
    return a18OH(T, eq)**theta


def a18CO2(T=273.15+22, eq="GZ19"):
    if (eq == "Beck05"):
        return np.exp(2.52*10**(-3) / (T-273.15)**2 + 0.01212)
    elif (eq == "GZ19"):
        return 0.01795 * 1000/T + 0.98124


def a17CO2(T=273.15+22, eq="GZ19"):
    theta = 209.3282 / (T**2) + -2.3984 / T + 0.5303
    return a18CO2(T, eq)**theta


# CREATE FIGURE S7
fig, ax = plt.subplots()

# Temperature in °C
T = 22

# Seawater
d18O_water = 0
Dp17O_water = -11
d17O_water = d17O(d18O_water, Dp17O_water)
ax.scatter(prime(d18O_water), Dp17O_water,
           marker="*", fc="k", ec="k")
ax.text(prime(d18O_water), Dp17O_water+10,
        f"seawater\n({d18O_water}‰, {Dp17O_water} ppm)",
        ha="center", va="bottom", color="k")

# Carbonate precipitating in equilibrium from seawater
d18Occ = d18O_cc(T + 273.15, d18O_water)
d17Occ = d17O_cc(T + 273.15, d17O(d18O_water, Dp17O_water))
Dp17Occ = Dp17O(d17Occ, d18Occ)
ax.scatter(prime(d18Occ), Dp17Occ,
           marker="o", fc="#858379", ec="k")

# Diffusion
Diff_shift = -3
Diff_theta = (np.log((12+16+16)/(12+17+16)))/(np.log((12+16+16)/(12+18+16)))
d18O_diff = d18Occ + Diff_shift
Dp17O_diff = apply_theta(
    d18Occ, Dp17Occ, shift_d18O=Diff_shift, theta=Diff_theta)

# Dissolved CO2
d18O_CO2 = A_from_a(a18CO2(T + 273.15, eq = "GZ19"), d18O_water)
d17O_CO2 = A_from_a(a17CO2(T + 273.15, eq = "GZ19"), d17O_water)
Dp17O_CO2 = Dp17O(d17O_CO2, d18O_CO2)
ax.scatter(prime(d18O_CO2), Dp17O_CO2,
           marker="D", fc="w", ec="k")
ax.text(prime(d18O_CO2) + 3, Dp17O_CO2,
        f"dissolved CO$_2$",
        ha="left", va="center", color="k")

# CO2 + KIE
CO2_KIE_shift = -3
CO2_KIE_theta = (np.log((12+16+16)/(12+17+16)))/(np.log((12+16+16)/(12+18+16)))
d18O_CO2_KIE = d18O_CO2 + CO2_KIE_shift
Dp17O_CO2_KIE = apply_theta(d18O_CO2, Dp17O_CO2, shift_d18O = CO2_KIE_shift, theta = CO2_KIE_theta)
ax.scatter(prime(d18O_CO2_KIE), Dp17O_CO2_KIE,
           marker="D", fc="k", ec="k")
ax.text(prime(d18O_CO2_KIE) + 3, Dp17O_CO2_KIE,
        r"CO$_2$ $\plus$ KIE",
        ha="left", va="center", color="k")

# Equilibrium OH-
d18O_OH_eq = B_from_a(a18OH(T=273.15 + T), d18O_water)
d17O_OH_eq = B_from_a(a17OH(T=273.15 + T), d17O_water)
Dp17O_OH_eq = Dp17O(d17O_OH_eq, d18O_OH_eq)
ax.scatter(prime(d18O_OH_eq), Dp17O_OH_eq,
           marker="s", fc="w", ec="k")
ax.text(prime(d18O_OH_eq), Dp17O_OH_eq - 10,
        r"OH$_{eq}^{-}$",
        ha="center", va="top", color="k")

# OH- + KIE
d18O_OH_eff = B_from_a(a18OH(T=273.15 + T, eq="BH21"), d18O_water)
d17O_OH_eff = B_from_a(a17OH(T=273.15 + T, eq="BH21", theta=0.523), d17O_water)
Dp17O_OH_eff = Dp17O(d17O_OH_eff, d18O_OH_eff)
ax.scatter(prime(d18O_OH_eff), Dp17O_OH_eff,
           marker="s", fc="k", ec="k")
ax.text(prime(d18O_OH_eff) + 3, Dp17O_OH_eff,
        r"OH$^{-}$ $\plus$ KIE",
        ha="left", va="center", color="k")

# Hydroxylation endmember carbonate
mixdfKIE = mix_d17O(d18O_A=d18O_OH_eff, D17O_A=Dp17O_OH_eff,
                    d18O_B=d18O_CO2_KIE, D17O_B=Dp17O_CO2_KIE)
ax.plot(prime(mixdfKIE["mix_d18O"]), mixdfKIE["mix_Dp17O"],
        color="k", ls=":", zorder=-1)
d18O_hydrox = mixdfKIE["mix_d18O"].iloc[67]
Dp17O_hydrox = mixdfKIE["mix_Dp17O"].iloc[67]
ax.scatter(prime(d18O_hydrox), Dp17O_hydrox,
           marker="o", fc="#EC0016", ec="k")
ax.annotate("", xycoords="data", textcoords="data",
            xy=(prime(d18Occ), Dp17Occ),
            xytext=(prime(d18O_hydrox), Dp17O_hydrox),
            arrowprops=dict(arrowstyle="<|-", color="#EC0016", lw=1.5), zorder=-1)
theta_hydrox = calculate_theta(d18Occ, Dp17Occ, d18O_hydrox, Dp17O_hydrox)
ax.text(prime(d18O_hydrox), Dp17O_hydrox - 10,
        r"$\theta^{}_{hydroxylation}$ = " + f"{theta_hydrox:.3f}",
        ha="center", va="top", color="#EC0016")


# Hydration
mixdf_hydra = mix_d17O(d18O_A=d18O_water, D17O_A=Dp17O_water,
                       d18O_B=d18O_CO2_KIE, D17O_B=Dp17O_CO2_KIE)
ax.plot(prime(mixdf_hydra["mix_d18O"]), mixdf_hydra["mix_Dp17O"],
        color="k", ls=":", zorder=-1)
d18O_hydra = mixdf_hydra["mix_d18O"].iloc[67]
Dp17O_hydra = mixdf_hydra["mix_Dp17O"].iloc[67]
ax.scatter(prime(d18O_hydra), Dp17O_hydra,
           marker="o", fc="#1455C0", ec="k")
ax.annotate("", xycoords="data", textcoords="data",
            xy=(prime(d18Occ), Dp17Occ),
            xytext=(prime(d18O_hydra), Dp17O_hydra),
            arrowprops=dict(arrowstyle="<|-", color="#1455C0", lw=1.5), zorder=-1)
theta_hydra = calculate_theta(d18Occ, Dp17Occ, d18O_hydra, Dp17O_hydra)
ax.text(prime(d18O_hydra) - 5, Dp17O_hydra,
        r"$\theta^{}_{hydration}$ = " + f"{theta_hydra:.3f}",
        bbox=dict(fc='white', ec="None", alpha=0.8, pad=0.2),
        ha="right", va="center", color="#1455C0")

# Arrow to show diffusion KIE
ax.annotate("", xycoords="data", textcoords="data",
            xy=(prime(d18Occ), Dp17Occ),
            xytext=(prime(d18O_diff), Dp17O_diff),
            arrowprops=dict(arrowstyle="<|-", color="#FF9B00", lw=1.5), zorder=-1)
plt.text(prime(d18O_diff), Dp17O_diff,
         "diffusion",
         ha="center", va="center", color="#FF9B00")

ax.set_xlabel("$\delta\prime^{18}$O (‰, VSMOW)")
ax.set_ylabel("$\Delta\prime^{17}$O (ppm)")

ax.set_ylim(-205, 205)
ax.set_xlim(-41, 71)

plt.tight_layout()
plt.savefig(os.path.join(sys.path[0], "SK Figure S7"))
plt.close()