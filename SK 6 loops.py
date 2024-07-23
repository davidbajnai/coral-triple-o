# This code is used to plot the disequilibrium loops for triple oxygen and clumped isotopes

# INPUT: isoDIC_***.csv
# OUTPUT: SK Figure S4.png

# >>>>>>>>>

# Import libraries
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Plot parameters
plt.rcParams["legend.loc"] = "best"
plt.rcParams.update({'font.size': 7})
plt.rcParams['scatter.edgecolors'] = "k"
plt.rcParams['scatter.marker'] = "o"
plt.rcParams["lines.linewidth"] = 0.5
plt.rcParams["patch.linewidth"] = 0.5
plt.rcParams["figure.figsize"] = (9, 4)
plt.rcParams["savefig.dpi"] = 600
plt.rcParams["savefig.bbox"] = "tight"
plt.rcParams['savefig.transparent'] = False
plt.rcParams['mathtext.default'] = 'regular'

# isoDIC models
isoDIC_header = pd.read_csv(os.path.join(sys.path[0], "isoDIC_header.csv"), sep=",")
isoDIC_cwc = pd.read_csv(os.path.join(sys.path[0], "isoDIC_pH8.8_T9.csv"), sep=",")
isoDIC_cwc.columns = isoDIC_header.columns
isoDIC_cwc = isoDIC_cwc - isoDIC_cwc.iloc[0]

isoDIC_wwc = pd.read_csv(os.path.join(sys.path[0], "isoDIC_pH8.5_T27.csv"), sep=",")
isoDIC_wwc.columns = isoDIC_header.columns
isoDIC_wwc = isoDIC_wwc - isoDIC_wwc.iloc[0]

isoDIC_cwc2 = pd.read_csv(os.path.join(sys.path[0], "isoDIC_pH8.5_T9.csv"), sep=",")
isoDIC_cwc2.columns = isoDIC_header.columns
isoDIC_cwc2 = isoDIC_cwc2 - isoDIC_cwc2.iloc[0]

isoDIC_wwc2 = pd.read_csv(os.path.join(sys.path[0], "isoDIC_pH8.8_T27.csv"), sep=",")
isoDIC_wwc2.columns = isoDIC_header.columns
isoDIC_wwc2 = isoDIC_wwc2 - isoDIC_wwc2.iloc[0]


# CREATE FIGURE S4
fig, (ax1, ax2) = plt.subplots(1, 2)

# Subplot a - triple oxygen isotpes

# Get the 15 minute marker position
cwc_target = (isoDIC_cwc["time(s)"] - 15*60).abs().idxmin()
wwc_target = (isoDIC_wwc["time(s)"] - 15*60).abs().idxmin()
cwc2_target = (isoDIC_cwc2["time(s)"] - 15*60).abs().idxmin()
wwc2_target = (isoDIC_wwc2["time(s)"] - 15*60).abs().idxmin()

# Add the isoDIC models
ax1.plot(isoDIC_wwc["d18_CO3"], isoDIC_wwc["D17_CO3"],
         c="darkred", ls="solid", zorder=1, lw=1.5,
         label=r"pH = 8.5, $\it{T}$ = 27 °C")
ax1.plot(isoDIC_cwc["d18_CO3"], isoDIC_cwc["D17_CO3"],
         c="darkblue", ls="solid", zorder=1, lw=1.5,
         label=r"pH = 8.8, $\it{T}$ = 9 °C")
ax1.plot(isoDIC_wwc2["d18_CO3"], isoDIC_wwc2["D17_CO3"],
         c="hotpink", ls="solid", zorder=2, lw=1,
         label=r"pH = 8.8, $\it{T}$ = 27 °C")
ax1.plot(isoDIC_cwc2["d18_CO3"], isoDIC_cwc2["D17_CO3"],
         c="cornflowerblue", ls="solid", zorder=2, lw=1,
         label=r"pH = 8.5, $\it{T}$ = 9 °C")

# Add the arrow markers
ax1.annotate("",
             (isoDIC_wwc["d18_CO3"].iloc[wwc_target],
              isoDIC_wwc["D17_CO3"].iloc[wwc_target]),
             (isoDIC_wwc["d18_CO3"].iloc[wwc_target-1],
              isoDIC_wwc["D17_CO3"].iloc[wwc_target-1]),
             ha="center", va="center", zorder=0,
             arrowprops=dict(arrowstyle="-|>", color="darkred", lw=1.5))
ax1.annotate("",
             (isoDIC_cwc["d18_CO3"].iloc[cwc_target],
              isoDIC_cwc["D17_CO3"].iloc[cwc_target]),
             (isoDIC_cwc["d18_CO3"].iloc[cwc_target-1],
              isoDIC_cwc["D17_CO3"].iloc[cwc_target-1]),
             ha="center", va="center", zorder=0,
             arrowprops=dict(arrowstyle="-|>", color="darkblue", lw=1.5))
ax1.annotate("",
             (isoDIC_wwc2["d18_CO3"].iloc[wwc2_target],
              isoDIC_wwc2["D17_CO3"].iloc[wwc2_target]),
             (isoDIC_wwc2["d18_CO3"].iloc[wwc2_target-1],
              isoDIC_wwc2["D17_CO3"].iloc[wwc2_target-1]),
             ha="center", va="center", zorder=2,
             arrowprops=dict(arrowstyle="-|>", color="hotpink", lw=1))
ax1.annotate("",
             (isoDIC_cwc2["d18_CO3"].iloc[cwc2_target],
              isoDIC_cwc2["D17_CO3"].iloc[cwc2_target]),
             (isoDIC_cwc2["d18_CO3"].iloc[cwc2_target-1],
              isoDIC_cwc2["D17_CO3"].iloc[cwc2_target-1]),
             ha="center", va="center", zorder=2,
             arrowprops=dict(arrowstyle="-|>", color="cornflowerblue", lw=1))

# Add the equilibrium point
ax1.scatter(0, 0, marker="X", color="#747067", label="Equilibrium", zorder=10, lw=1.5)

# Show the commonly achieved measurement precision
ax1.text(0.98, 0.02, "Commonly achieved\nmeasurement precision", size=6, ha="right", va="bottom", transform=ax1.transAxes)
ax1.errorbar(-0.5, -60,
             xerr=0.05, yerr=5,
             fmt="none", color="k", lw =.8)

# Axis properties
ax1.set_xlim(-6.2, 0.2)
ax1.set_ylim(-72, 2)
ax1.set_ylabel("$\Delta\prime ^{17}O$ disequilibrium (ppm)")
ax1.set_xlabel("$\delta^{18}O$ disequilibrium (‰)")
ax1.text(0.02, 0.98, "(a)", size=10, ha="left", va="top",
         transform=ax1.transAxes)

# Subplot b - clumped isotpes

# Add the isoDIC models
ax2.plot(isoDIC_wwc["D64_CO3"], isoDIC_wwc["D63_CO3"],
         c="darkred", ls="solid", zorder=1, lw=1.5, label=r"pH = 8.5, $\it{T}$ = 27 °C")
ax2.plot(isoDIC_cwc["D64_CO3"], isoDIC_cwc["D63_CO3"],
         c="darkblue", ls="solid", zorder=1, lw=1.5, label=r"pH = 8.8, $\it{T}$ = 9 °C")
ax2.plot(isoDIC_wwc2["D64_CO3"], isoDIC_wwc2["D63_CO3"]-isoDIC_wwc2["D63_CO3"].iloc[0],
         c="hotpink", ls="solid", zorder=2, lw=1, label=r"pH = 8.8, $\it{T}$ = 27 °C")
ax2.plot(isoDIC_cwc2["D64_CO3"], isoDIC_cwc2["D63_CO3"],
         c="cornflowerblue", ls="solid", zorder=2, lw=1, label=r"pH = 8.5, $\it{T}$ = 9 °C")

# Add the arrow markers
ax2.annotate("",
             (isoDIC_wwc["D64_CO3"].iloc[wwc_target],
              isoDIC_wwc["D63_CO3"].iloc[wwc_target]),
             (isoDIC_wwc["D64_CO3"].iloc[wwc_target-1],
              isoDIC_wwc["D63_CO3"].iloc[wwc_target-1]),
             ha="center", va="center", zorder=0,
             arrowprops=dict(arrowstyle="-|>", color="darkred", lw=1.5))
ax2.annotate("",
             (isoDIC_cwc["D64_CO3"].iloc[cwc_target],
              isoDIC_cwc["D63_CO3"].iloc[cwc_target]),
             (isoDIC_cwc["D64_CO3"].iloc[cwc_target-1],
              isoDIC_cwc["D63_CO3"].iloc[cwc_target-1]),
             ha="center", va="center", zorder=0,
             arrowprops=dict(arrowstyle="-|>", color="darkblue", lw=1.5))
ax2.annotate("",
             (isoDIC_wwc2["D64_CO3"].iloc[wwc2_target],
              isoDIC_wwc2["D63_CO3"].iloc[wwc2_target]),
             (isoDIC_wwc2["D64_CO3"].iloc[wwc2_target-1],
              isoDIC_wwc2["D63_CO3"].iloc[wwc2_target-1]),
             ha="center", va="center", zorder=2,
             arrowprops=dict(arrowstyle="-|>", color="hotpink", lw=1))
ax2.annotate("",
             (isoDIC_cwc2["D64_CO3"].iloc[cwc2_target],
              isoDIC_cwc2["D63_CO3"].iloc[cwc2_target]),
             (isoDIC_cwc2["D64_CO3"].iloc[cwc2_target-1],
              isoDIC_cwc2["D63_CO3"].iloc[cwc2_target-1]),
             ha="center", va="center", zorder=2,
             arrowprops=dict(arrowstyle="-|>", color="cornflowerblue", lw=1))

# Add the equilibrium point
ax2.scatter(0, 0, marker="X", color="#747067", label="Equilibrium", zorder=10, lw=1.5)

# Show the commonly achieved measurement precision
ax2.text(0.98, 0.02, "Commonly achieved\nmeasurement precision", size=6, ha="right", va="bottom",
         transform=ax2.transAxes)
ax2.errorbar(-0.025, -0.065,
             xerr=0.020, yerr=0.006,
             fmt="none", color="k", lw=.8)

# Legend
ax2.legend(loc='upper right', bbox_to_anchor=(1.45, 1))

# Axis properties
ax2.set_ylabel("$\Delta_{47}$ disequilibrium (‰)")
ax2.set_xlabel("$\Delta_{48}$ disequilibrium (‰)")
ax2.text(0.02, 0.98, "(b)", size=10, ha="left", va="top",
         transform=ax2.transAxes)

plt.savefig(os.path.join(sys.path[0], "SK Figure S4"))
plt.close()