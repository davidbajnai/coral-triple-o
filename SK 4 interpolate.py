# INPUT: SK Table S-3 part-3.csv
# OUTPUT: SK Figure 3.png, SK Figure S6.png, SK Table S-3.csv, SK Table S-3.xlsx

# >>>>>>>>>

# Import libraries
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist

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

# Functions that make life easier
def prime(x):
    return 1000 * np.log(x / 1000 + 1)


def unprime(x):
    return (np.exp(x / 1000) - 1) * 1000


def a18_cc(T):

    # Used for Figure 2:                                   Hayles et al. (2018) - calcite
    B_calcite = 7.027321E+14 / T**7 + -1.633009E+13 / T**6 + 1.463936E+11 / T**5 + -5.417531E+08 / T**4 + -4.495755E+05 / T**3  + 1.307870E+04 / T**2 + -5.393675E-01 / T + 1.331245E-04
    B_water = -6.705843E+15 / T**7 + 1.333519E+14 / T**6 + -1.114055E+12 / T**5 + 5.090782E+09 / T**4 + -1.353889E+07 / T**3 + 2.143196E+04 / T**2 + 5.689300 / T + -7.839005E-03
    return np.exp(B_calcite) / np.exp(B_water)

    # Use this for Figure S5
    # return np.exp((2.84 * 10**6 / T**2 - 2.96) / 1000) # Wostbrock et al. (2020) – calcite
    
    # Alternative equations
    # return np.exp((17.88 * 1000 / T - 31.14) / 1000)   # Kim et al. (2007) – aragonite
    # return np.exp((17.57 * 1000 / T - 29.13) / 1000)   # Daeron et al. (2019) – calcite
    # return 0.0201 * (1000 / T) + 0.9642                # Guo and Zhou (2019) – aragonite


def theta_cc(T):
    # Used for Figure 2:                         Hayles et al. (2018) - calcite
    K_calcite = 1.019124E+09 / T**5 + -2.117501E+07 / T**4 + 1.686453E+05 / T**3 + -5.784679E+02 / T**2 + 1.489666E-01 / T + 0.5304852
    B_calcite = 7.027321E+14 / T**7 + -1.633009E+13 / T**6 + 1.463936E+11 / T**5 + -5.417531E+08 / T**4 + -4.495755E+05 / T**3  + 1.307870E+04 / T**2 + -5.393675E-01 / T + 1.331245E-04
    K_water = 7.625734E+06 / T**5 + 1.216102E+06 / T**4 + -2.135774E+04 / T**3 + 1.323782E+02 / T**2 + -4.931630E-01 / T + 0.5306551
    B_water = -6.705843E+15 / T**7 + 1.333519E+14 / T**6 + -1.114055E+12 / T**5 + 5.090782E+09 / T**4 + -1.353889E+07 / T**3 + 2.143196E+04 / T**2 + 5.689300 / T + -7.839005E-03
    a18 = np.exp(B_calcite) / np.exp(B_water)
    return K_calcite + (K_calcite-K_water) * (B_water / np.log(a18))

    # Use this for Figure S5
    # return -1.39 / T + 0.5305                 # Wostbrock et al. (2020) – calcite

    # Alternative equations
    # return 59.1047/T**2 + -1.4089/T + 0.5297  # Guo and Zhou (2019) – aragonite
    # return -1.53 / T + 0.5305                 # Wostbrock et al. (2020) – aragonite


def a17_cc(T):
    return a18_cc(T)**theta_cc(T)


def d18O_cc(equilibrium_temperatures, d18Ow):
    return a18_cc(equilibrium_temperatures + 273.15) * (d18Ow+1000) - 1000


def get_18O_temp(d18O_coral, d18O_coral_err, d18O_seawater, d18O_seawater_err):

    a18 = (d18O_coral + 1000) / (d18O_seawater + 1000)
    T =  1000 / ((a18 - 0.9642) / 0.0201) - 273.15

    a18_min = (d18O_coral + d18O_coral_err + 1000) / (d18O_seawater - d18O_seawater_err + 1000)
    T_min = 1000 / ((a18_min - 0.9642) / 0.0201) - 273.15

    a18_max = (d18O_coral - d18O_coral_err + 1000) / (d18O_seawater + d18O_seawater_err + 1000)
    T_max = 1000 / ((a18_max - 0.9642) / 0.0201) - 273.15

    return T, (T_max-T_min)/2

def d17O_cc(equilibrium_temperatures, d17Ow):
    return a17_cc(equilibrium_temperatures + 273.15) * (d17Ow+1000) - 1000


def Dp17O(d17O, d18O):
    return (prime(d17O) - 0.528 * prime(d18O)) * 1000


def d17O(d18O, Dp17O):
    return unprime(Dp17O / 1000 + 0.528 * prime(d18O))


def calculate_theta(d18O_A, Dp17O_A, d18O_B, Dp17O_B):
    a18 = (d18O_B + 1000) / (d18O_A + 1000)
    a17 = (d17O(d18O_B, Dp17O_B) + 1000) / (d17O(d18O_A, Dp17O_A) + 1000)
    return np.log(a17) / np.log(a18)


def apply_theta(d18O_A, Dp17O_A, d18O_B=None, shift_d18O=None, theta=None):
    if d18O_B == None:
        d18O_B = d18O_A + shift_d18O

    a18 = (d18O_B + 1000) / (d18O_A + 1000)
    a17 = a18**theta

    d17O_B = a17 * (d17O(d18O_A, Dp17O_A) + 1000) - 1000
    Dp17O_B = Dp17O(d17O_B, d18O_B)

    return Dp17O_B


def cc_equilibrium(T, T_err, d18Ow, d18Ow_err, Dp17Ow, Dp17Ow_err, Sample=None):
    # Calculate aragonite equilibrium with error propagation

    df = pd.DataFrame([])

    if Sample is not None:
        df["SampleName"] = Sample

    df["T"] = T
    df["T_err"] = T_err
    df["d18Ow"] = d18Ow
    df["d18Ow_err"] = d18Ow_err
    df["Dp17Ow"] = Dp17Ow
    df["Dp17Ow_err"] = Dp17Ow_err
    
    df["d18O_equi"] = d18O_cc(T, d18Ow)
    df["d17O_equi"] = d17O_cc(T, unprime((Dp17Ow / 1000) + 0.528 * prime(d18Ow)))
    df["Dp17O_equi"] = Dp17O(df["d17O_equi"], df["d18O_equi"])

    df["d18O_equi_min"] = d18O_cc(T+T_err, d18Ow-d18Ow_err)
    df["d17O_equi_min"] = d17O_cc(T+T_err, unprime((Dp17Ow / 1000) + 0.528 * prime(d18Ow-d18Ow_err)))
    df["Dp17O_equi_max"] = Dp17O(df["d17O_equi_min"], df["d18O_equi_min"])

    df["d18O_equi_max"] = d18O_cc(T-T_err, d18Ow+d18Ow_err)
    df["d17O_equi_max"] = d17O_cc(T-T_err, unprime((Dp17Ow / 1000) + 0.528 * prime(d18Ow+d18Ow_err)))
    df["Dp17O_equi_min"] = Dp17O(df["d17O_equi_max"], df["d18O_equi_max"])

    df["d18O_equi_err"] = (df["d18O_equi_max"] - df["d18O_equi_min"]) / 2
    df["d17O_equi_err"] = (df["d17O_equi_max"] - df["d17O_equi_min"]) / 2
    df["Dp17O_equi_err"] = np.sqrt(((df["Dp17O_equi_max"] - df["Dp17O_equi_min"]) / 2)**2 + Dp17Ow_err**2)

    # Keep only the equilibrium values
    df = df.loc[:, ["SampleName", "d18O_equi",
                    "d18O_equi_err", "Dp17O_equi", "Dp17O_equi_err"]]
    df = df.rename(columns={"d18O_equi": "d18O_equilibrium",
                            "d18O_equi_err": "d18O_equilibrium_err",
                            "Dp17O_equi": "Dp17O_equilibrium",
                            "Dp17O_equi_err": "Dp17O_equilibrium_err"})
    return df


def plot_equilibrium(Dp17Ow, d18Ow, ax, color="k"):

    d17Ow = unprime(0.528 * prime(d18Ow) + Dp17Ow/1000)

    ax.scatter(prime(d18Ow), Dp17O(d17Ow, d18Ow),
            marker="X", fc=color, ec="w", zorder=10)

    # plot equilibrium line, entire T range
    toInf = np.arange(-10, 290, 2)
    d18O_mineral = d18O_cc(toInf, d18Ow)
    d17O_mineral = d17O_cc(toInf, d17Ow)
    mineral_equilibrium = np.array(
        [d18O_mineral, Dp17O(d17O_mineral, d18O_mineral), toInf]).T
    ax.plot(prime(mineral_equilibrium[:, 0]), mineral_equilibrium[:, 1],
            ":", c=color, zorder=3)
    
    # plot equilibrium line, entire T range
    toInf = np.arange(-10, 60, 0.05)
    d18O_mineral = d18O_cc(toInf, d18Ow)
    d17O_mineral = d17O_cc(toInf, d17Ow)
    mineral_equilibrium = np.array(
        [d18O_mineral, Dp17O(d17O_mineral, d18O_mineral), toInf]).T
    ax.plot(prime(mineral_equilibrium[:, 0]), mineral_equilibrium[:, 1],
            ":", c=color, zorder=3)
    
    # Return equilibrium data as a dataframe
    equilibrium_df = pd.DataFrame(mineral_equilibrium)
    equilibrium_df[2] = equilibrium_df[2]
    equilibrium_df = equilibrium_df.rename(
        columns={0: 'd18O', 1: 'Dp17O', 2: 'temperature'})

    return equilibrium_df


def vital_vector(d18O_coral, Dp17O_coral, num_points, shift_d18O, theta_coral):
    new_Dp17O = apply_theta(d18O_A=d18O_coral, Dp17O_A=Dp17O_coral,
                            shift_d18O=shift_d18O, theta=theta_coral)
    x1, y1 = d18O_coral, Dp17O_coral
    x2, y2 = d18O_coral+shift_d18O, new_Dp17O
    slope2 = (y2 - y1) / (x2 - x1)
    intercept2 = y1 - slope2 * x1
    x_values = np.linspace(x1, x2, num_points)
    y_values = slope2 * x_values + intercept2
    return pd.DataFrame({'d18O': x_values, 'Dp17O': y_values})


def get_17O_temp(d18O_coral, d18O_coral_error, Dp17O_coral, Dp17O_coral_error, d18O_seawater, d18O_seawater_err, Dp17O_seawater, Dp17O_seawater_err, theta_coral, ax):

    shift_d18O = 15

    # Get equilibrium values for seawater
    df_eq = plot_equilibrium(Dp17Ow=Dp17O_seawater, d18Ow=d18O_seawater,
                             ax=ax, color="k")
    ax.errorbar(prime(d18O_seawater), Dp17O_seawater,
                    xerr=d18O_seawater_err, yerr=Dp17O_seawater_err,
                    fmt="none", color="k", zorder=-1)
    ax.annotate("carbonate equilibrium", xycoords="data", textcoords="data",
            xy = (prime(df_eq["d18O"]).iloc[-1], df_eq["Dp17O"].iloc[-1]),
            xytext = (prime(df_eq["d18O"]).iloc[-1], df_eq["Dp17O"].iloc[-1]+20),
            ha="left", va="center",
            arrowprops=dict(arrowstyle="->", color="k"))
    ax.text(d18O_seawater-1, Dp17O_seawater-10, "ambient\nseawater", ha="left", va="top")
    ax.text(d18O_coral-2, Dp17O_coral, "coral", ha="right", va="center")
    
    # Get the vital effect line
    df_line = vital_vector(d18O_coral, Dp17O_coral, len(df_eq["d18O"]), shift_d18O, theta_coral)

    # Calculate the distance between points
    df_eq_f = df_eq.loc[:, ["d18O", "Dp17O"]]
    distances = cdist(df_eq_f, df_line)
    min_indices = np.unravel_index(np.argmin(distances), distances.shape)
    closest_point_eq = df_eq_f.iloc[min_indices[0]]

    temp = df_eq.iloc[min_indices[0]]["temperature"]

    # Plot
    ax.plot(prime(df_line["d18O"]), df_line["Dp17O"],"k", ls="-")
    ax.scatter(prime(closest_point_eq.iloc[0]), closest_point_eq.iloc[1], marker=".", fc="k")

    # Maximum temperature
    df_eq = plot_equilibrium(Dp17Ow=Dp17O_seawater - Dp17O_seawater_err, d18Ow=d18O_seawater + d18O_seawater_err,
                             ax=ax, color="red")
    df_line = vital_vector(d18O_coral - d18O_coral_error, Dp17O_coral + Dp17O_coral_error, len(df_eq["d18O"]), shift_d18O, theta_coral)
    df_eq_f = df_eq.loc[:, ["d18O", "Dp17O"]]
    distances = cdist(df_eq_f, df_line)
    min_indices = np.unravel_index(np.argmin(distances), distances.shape)
    closest_point_eq = df_eq_f.iloc[min_indices[0]]
    T_max = df_eq.iloc[min_indices[0]]["temperature"]
    ax.plot(prime(df_line["d18O"]), df_line["Dp17O"],"red", ls="-")
    ax.scatter(prime(closest_point_eq.iloc[0]), closest_point_eq.iloc[1], marker=".", fc="red")


    # Minimum temperature
    df_eq = plot_equilibrium(Dp17Ow=Dp17O_seawater + Dp17O_seawater_err, d18Ow=d18O_seawater - d18O_seawater_err,
                             ax=ax, color="blue")
    df_line = vital_vector(d18O_coral + d18O_coral_error, Dp17O_coral - Dp17O_coral_error, len(df_eq["d18O"]), shift_d18O, theta_coral)
    df_eq_f = df_eq.loc[:, ["d18O", "Dp17O"]]
    distances = cdist(df_eq_f, df_line)
    min_indices = np.unravel_index(np.argmin(distances), distances.shape)
    closest_point_eq = df_eq_f.iloc[min_indices[0]]
    T_min = df_eq.iloc[min_indices[0]]["temperature"]
    ax.plot(prime(df_line["d18O"]), df_line["Dp17O"],"blue", ls="-")
    ax.scatter(prime(closest_point_eq.iloc[0]), closest_point_eq.iloc[1], marker=".", fc="blue")

    return temp, (T_max-T_min)/2


# Read data from CSV files
df = pd.read_csv(sys.path[0] + "/SK Table S-3 part-3.csv", sep=",")


# Assign colors and markers
cat1 = df["Species"].unique()
markers = dict(zip(cat1, ["o", "s", "D", "v", "^", "<", ">", "p", "P", "*"]))
cat2 = df["Type"].unique()
colors = dict(zip(cat2, ["#1455C0", "#EC0016"]))

# Do the sensitivity analysis here (uncomment lines to test different scenarios)
# -> What happens if we change the measurement error
# Dp17O_error = df["Dp17O_error"].mean()
# df["Dp17O_error"] = Dp17O_error - 1
# print("The mean cooral Dp17Oc error is: " + str(round(df["Dp17O_error"].mean(),0)) + " ppm; " + f"∆error = {round(df['Dp17O_error'].mean() - Dp17O_error, 0)}")
# -> What happens if we change the seawater Dp17O error
# Dp17O_error = df["Dp17Osw_err"].mean()
# df["Dp17Osw_err"] = 0
# print("The mean seawater Dp17O error is: " + str(round(df["Dp17Osw_err"].mean(),0)) + " ppm; " + f"∆error = {round(df['Dp17Osw_err'].mean() - Dp17O_error, 0)}")

# Figure S6 - to check if the interpolation works

# Uncomment this line to produce Figure S6
# df = df[df["SampleName"] == "SK-GeoB"]

fig, ax = plt.subplots(1, 1, figsize=(4, 4))

# Get the "vital effect"-corrected temperatures
theta_coral = df["theta_coral"].median()
df['T_17O_tuple'] = df.apply(lambda row: get_17O_temp(d18O_coral=row["d18O_AC"],
                                                      Dp17O_coral=row["Dp17O_AC"],
                                                      d18O_coral_error=row["d18O_error"],
                                                      Dp17O_coral_error=row["Dp17O_error"],
                                                      d18O_seawater=row["d18Osw_modeled"],
                                                      d18O_seawater_err=row["d18Osw_modeled_err"],
                                                      Dp17O_seawater=row["Dp17Osw"],
                                                      Dp17O_seawater_err=row["Dp17Osw_err"],
                                                      theta_coral=theta_coral,
                                                      ax=ax),
                                                    #   theta_coral=row["theta_coral"]), # uncomment this line to use the individual theta values
                             axis=1)

df['T_17O'], df['T_17O_error'] = zip(*df['T_17O_tuple'])
del df['T_17O_tuple']
df["T_18O"], df["T_18O_error"] = get_18O_temp(df["d18O_AC"], df["d18O_error"],
                                              df["d18Osw_modeled"], df["d18Osw_modeled_err"])
print(f'The mean error of the 17O-based temperatures is: {df["T_17O_error"].mean():.1f} °C')
print(f'The mean error of the 18O-based temperatures is: {df["T_18O_error"].mean():.1f} °C')

# Print  the temperature difference - not used in the paper
# print("\nThe average difference between the 17O-based and the MODELED temperatures:")
# print(f'{df[(df["Type"] == "cold-water coral")]["T_17O"].mean() - df[(df["Type"] == "cold-water coral")]["T_modeled"].mean():.1f} °C (COLD-WATER CORALS)')
# print(f'{df[df["Type"] == "warm-water coral"]["T_17O"].mean() - df[df["Type"] == "warm-water coral"]["T_modeled"].mean():.1f} °C (WARM-WATER CORALS)')

# print("\nThe average difference between the 17O-based and the MEASURED temperatures: ")
# print(f'{df[(df["Type"] == "cold-water coral")]["T_17O"].mean() - df[(df["Type"] == "cold-water coral")]["T_measured"].mean():.1f} °C (COLD-WATER CORALS)')
# print(f'{df[df["Type"] == "warm-water coral"]["T_17O"].mean() - df[df["Type"] == "warm-water coral"]["T_measured"].mean():.1f} °C (WARM-WATER CORALS)')

# Create a separate scatter plot for each species
ax.scatter(prime(df["d18O_AC"]), df["Dp17O_AC"],
           marker="o", fc="k")
ax.errorbar(prime(df["d18O_AC"]), df["Dp17O_AC"],
            xerr=df["d18O_error"],
            yerr=df["Dp17O_error"],
            fmt="none", color="k", zorder=-1)

ax.set_ylabel("$\Delta^{\prime 17}$O (ppm)")
ax.set_xlabel("$\delta^{\prime 18}$O (‰, VSMOW)")

plt.savefig(sys.path[0] + "/SK Figure S6.png")
plt.close()


# Create Figure 3
fig, (ax1, ax2) = plt.subplots(1, 2)

# Subplot A

for cat in cat1:
    for dog in cat2:
        data = df[(df["Species"] == cat) & (df["Type"] == dog)]
        if len(data) > 0:
            x = data["T_modeled"]
            y = data["T_18O"]
            xerr = data["T_modeled_err"]
            yerr = data["T_18O_error"]
            ax1.scatter(x, y,
                        marker=markers[cat], fc=colors[dog], label=f"{cat}")
            ax1.errorbar(x, y, xerr=xerr, yerr=yerr,
                            fmt="none", color=colors[dog], zorder=0)

ax1.set_xlim(-1, 31)
ax1.set_ylim(-6, 66)
ylim = ax1.get_ylim()
xlim = ax1.get_xlim()

# 1:1 line
ax1.plot([-10, 100], [-10, 100], c = "k", ls="dashed", zorder = -1)
xmin, xmax = ax1.get_xlim()
ymin, ymax = ax1.get_ylim()
angle = np.arctan((xmax-xmin)/(ymax-ymin)) * 180 / np.pi
ax1.text(20, 18, "1:1", ha="center", va="center", rotation=angle)

ax1.set_ylabel(r"Temperature from $\delta^{18}O$-thermometry (°C)")
ax1.set_xlabel("Growth temperature (°C)")

ax1.text(0.08, 0.98, "a", size=14, ha="right", va="top",
         transform=ax1.transAxes, fontweight="bold")


# Subplot B

for cat in cat1:
    for dog in cat2:
        data = df[(df["Species"] == cat) & (df["Type"] == dog)]
        if len(data) > 0:
            x = data["T_modeled"]
            y = data["T_17O"]
            xerr = data["T_modeled_err"]
            yerr = data["T_17O_error"]
            ax2.scatter(x, y,
                        marker=markers[cat], fc=colors[dog], label=f"{cat}")
            ax2.errorbar(x, y, xerr=xerr, yerr=yerr,
                         fmt="none", color=colors[dog], zorder=0)

# 1:1 line
ax2.plot([-10, 100], [-10, 100], c = "k", ls="dashed", zorder = -1)
ax2.text(20, 18, "1:1", ha="center", va="center", rotation=angle)

# Add legend and format species names to italic
legend = ax2.legend(loc='upper right', bbox_to_anchor=(1.45, 1))
for text in legend.texts:
    text.set_fontsize(5.5)
    text.set_fontstyle('italic')

ax2.set_ylabel(r"Temperature corrected for 'vital effects' ($\it{T}_{\Delta\prime^{17}O}$, °C)")
ax2.set_xlabel("Growth temperature (°C)")

ax2.text(0.08, 0.98, "b", size=14, ha="right", va="top",
         transform=ax2.transAxes, fontweight="bold")

ax2.set_ylim(ylim)
ax2.set_xlim(xlim)

plt.savefig(sys.path[0] + "/SK Figure 3.png")
plt.close("all")

df.to_csv(sys.path[0] + "/SK Table S-3.csv", index=False)
df.to_excel(sys.path[0] + "/SK Table S-3.xlsx", index=False)
