# This code is used to plot the coral data in triple oxygen isotope space,
# and calculate the vital effect theta

# INPUT: SK Table S-3 part-2.csv, seawater.csv, isoDIC_***.csv
# OUTPUT: SK Figure 2.png, (SK Figure S5.png), SK Table S-3 part-3.csv

# >>>>>>>>>

# Import libraries
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import zoom


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

    # Use this for Figure S-4
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

    # Use this for Figure S4
    # return -1.39 / T + 0.5305                 # Wostbrock et al. (2020) – calcite

    # Alternative equations
    # return 59.1047/T**2 + -1.4089/T + 0.5297  # Guo and Zhou (2019) – aragonite
    # return -1.53 / T + 0.5305                 # Wostbrock et al. (2020) – aragonite


def a17_cc(T):
    return a18_cc(T)**theta_cc(T)


def d18O_cc(equilibrium_temperatures, d18Ow):
    return a18_cc(equilibrium_temperatures + 273.15) * (d18Ow+1000) - 1000


def d17O_cc(equilibrium_temperatures, d17Ow):
    return a17_cc(equilibrium_temperatures + 273.15) * (d17Ow+1000) - 1000


def Dp17O(d17O, d18O):
    return (prime(d17O) - 0.528 * prime(d18O)) * 1000


def d17O(d18O, Dp17O):
    return unprime(Dp17O / 1000 + 0.528 * prime(d18O))


def calculate_theta(d18O_A, Dp17O_A, d18O_B, Dp17O_B):
    a18 = (d18O_B + 1000) / (d18O_A + 1000)
    a17 = (d17O(d18O_B, Dp17O_B) + 1000) / (d17O(d18O_A, Dp17O_A) + 1000)

    theta = round(np.log(a17) / np.log(a18), 4)

    return theta


def apply_theta(d18O_A, Dp17O_A, d18O_B=None, shift_d18O=None, theta=None):
    if d18O_B == None:
        d18O_B = d18O_A + shift_d18O

    a18 = (d18O_B + 1000) / (d18O_A + 1000)
    a17 = a18**theta

    d17O_B = a17 * (d17O(d18O_A, Dp17O_A) + 1000) - 1000
    Dp17O_B = Dp17O(d17O_B, d18O_B)

    return Dp17O_B


def cc_equilibrium(T, T_err, d18Ow, d18Ow_err, Dp17Ow, Dp17Ow_err, Sample=None):
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


def plot_equilibrium(Dp17Ow, d18Ow, Tmin, Tmax, ax, fluid_name="precipitating fluid", color="k", highlight=True, mark_water = True):

    d17Ow = unprime(0.528 * prime(d18Ow) + Dp17Ow/1000)

    # mark water
    if mark_water == True:
        ax.scatter(prime(d18Ow), Dp17O(d17Ow, d18Ow),
                marker="X", fc=color, ec="w",
                zorder=10, label=fluid_name)

    # equilibrium line, entire T range
    toInf = np.arange(Tmin, Tmax, 1)
    d18O_mineral = d18O_cc(toInf, d18Ow)
    d17O_mineral = d17O_cc(toInf, d17Ow)
    mineral_equilibrium = np.array(
        [d18O_mineral, Dp17O(d17O_mineral, d18O_mineral), toInf]).T
    ax.plot(prime(mineral_equilibrium[:, 0]), mineral_equilibrium[:, 1],
            ":", c=color, zorder=3)

    # equilibrium points, highlight range
    equilibrium_temperatures = np.arange(Tmin, Tmax, 0.5)
    colors = np.linspace(0, 1, len(equilibrium_temperatures))
    d18O_mineral = d18O_cc(equilibrium_temperatures, d18Ow)
    d17O_mineral = d17O_cc(equilibrium_temperatures, d17Ow)
    mineral_equilibrium = np.array([d18O_mineral, Dp17O(
        d17O_mineral, d18O_mineral), equilibrium_temperatures]).T
    if highlight == True:
        ax.scatter(prime(mineral_equilibrium[:, 0]), mineral_equilibrium[:, 1],
                   marker=".", c=colors, cmap='coolwarm', linewidths=0, zorder=3)
    
    # Return equilibrium data as a dataframe
    equilibrium_df = pd.DataFrame(mineral_equilibrium)
    equilibrium_df[2] = equilibrium_df[2]
    equilibrium_df = equilibrium_df.rename(
        columns={0: 'd18O', 1: 'Dp17O', 2: 'temperature'})

    # equilibrium, highlight range, marker every 10 °C
    equilibrium_temperatures = np.arange(Tmin, Tmax+1, 10)
    d18O_mineral = d18O_cc(equilibrium_temperatures, d18Ow)
    d17O_mineral = d17O_cc(equilibrium_temperatures, d17Ow)
    mineral_equilibrium = np.array([d18O_mineral, Dp17O(
        d17O_mineral, d18O_mineral), equilibrium_temperatures]).T
    if highlight == True:
        ax.scatter(prime(mineral_equilibrium[:, 0]), mineral_equilibrium[:, 1],
                   s=15, marker="o", fc="white", ec=color, zorder=3)

    return equilibrium_df



# Read data from CSV files
swdf = pd.read_csv(sys.path[0] + "/seawater.csv", sep=",")

# isoDIC models
isoDIC_header = pd.read_csv(sys.path[0] + "/isoDIC_header.csv", sep=",")
isoDIC_cwc = pd.read_csv(sys.path[0] + "/isoDIC_pH8.8_T9.csv", sep=",")
isoDIC_cwc.columns = isoDIC_header.columns
isoDIC_wwc = pd.read_csv(sys.path[0] + "/isoDIC_pH8.5_T27.csv", sep=",")
isoDIC_wwc.columns = isoDIC_header.columns

df = pd.read_csv(sys.path[0] + "/SK Table S-3 part-2.csv", sep=",")

# Set Dp17Osw and Dp17Osw_err
df["Dp17Osw"], df["Dp17Osw_err"] = -11, 6

# Calculate equilibrium values using the "measured + modeled" d18Osw and T values
df_equi = cc_equilibrium(T=df["T_modeled"], T_err=df["T_modeled_err"],
                                d18Ow=df["d18Osw_modeled"], d18Ow_err=df["d18Osw_modeled_err"],
                                Dp17Ow=df["Dp17Osw"], Dp17Ow_err=df["Dp17Osw_err"],
                                Sample=df["SampleName"])
df = pd.merge(df, df_equi, on="SampleName", how="left")
df["d18O_offset"] = df["d18O_AC"]-df["d18O_equilibrium"]
df["Dp17O_offset"] = df["Dp17O_AC"]-df["Dp17O_equilibrium"]

# Calculate the effective theta for coral vital effects

# # MONTE CARLO SIMULATION
# monte_carlo_iterations = 10**3

# # SIMULATION FOR COLD-WATER CORALS
# df_cwc = df[df["Type"] == "cold-water coral"].reset_index(drop=True)
# df_cwc["theta"] = np.nan
# df_cwc["theta_err"] = np.nan
# # Plots for testing the monte carlo simulation
# # plt.plot(prime(df_cwc["d18O_equilibrium"]), df_cwc["Dp17O_equilibrium"],
# #          marker="x", color="gray")
# # plt.errorbar(prime(df_cwc["d18O_equilibrium"]), df_cwc["Dp17O_equilibrium"],
# #              xerr=df_cwc["d18O_equilibrium_err"], yerr=df_cwc["Dp17O_equilibrium_err"],
# #              fmt="none", color="gray", elinewidth=0.8)
# # plt.plot(prime(df_cwc["d18O_AC"]), df_cwc["Dp17O_AC"],
# #          marker="x", color="blue")
# # plt.errorbar(prime(df_cwc["d18O_AC"]), df_cwc["Dp17O_AC"],
# #              xerr=df_cwc["d18O_error"]/np.sqrt(df_cwc["Replicates"]), yerr=df_cwc["Dp17O_error"]/np.sqrt(df_cwc["Replicates"]),
# #              fmt="none", color="blue", elinewidth=0.8)
# for i in range(len(df_cwc)):
#     theta_cwc_lst = []
#     for _ in range(monte_carlo_iterations):

#         d18O_equilibrium = np.random.normal(df_cwc["d18O_equilibrium"].iloc[i], df_cwc["d18O_equilibrium_err"].iloc[i])
#         Dp17O_equilibrium = np.random.normal(df_cwc["Dp17O_equilibrium"].iloc[i], df_cwc["Dp17O_equilibrium_err"].iloc[i])
#         d18O_measured = np.random.normal(df_cwc["d18O_AC"].iloc[i], df_cwc["d18O_error"].iloc[i])
#         Dp17O_measured = np.random.normal(df_cwc["Dp17O_AC"].iloc[i], df_cwc["Dp17O_error"].iloc[i])

#         # plt.plot([prime(d18O_equilibrium), prime(d18O_measured)],
#         #             [Dp17O_equilibrium, Dp17O_measured],
#         #             c="lightblue", ls="solid", lw=0.8, zorder=-1, alpha=0.2)

#         result_theta = calculate_theta(d18O_A=d18O_equilibrium, Dp17O_A=Dp17O_equilibrium,
#                                        d18O_B=d18O_measured, Dp17O_B=Dp17O_measured)


#         theta_cwc_lst.append(result_theta)
#     df_cwc.loc[i, "theta"] = np.mean(theta_cwc_lst)
#     df_cwc.loc[i, "theta_err"] = np.std(theta_cwc_lst)

# df_cwc = df_cwc.loc[:, ["SampleName", "theta", "theta_err"]]

# theta_cwc = np.mean(df_cwc["theta"])
# theta_cwc_err = np.std(df_cwc["theta"])
# print(f'The mean effective theta for cold-water corals is {theta_cwc:.3f}(±{theta_cwc_err:.3f})')

# # SIMULATION FOR WARM-WATER CORALS
# df_wwc = df[df["Type"] == "warm-water coral"].reset_index(drop=True)
# df_wwc["theta"] = np.nan
# df_wwc["theta_err"] = np.nan
# # Plots for testing the monte carlo simulation
# # plt.plot(prime(df_wwc["d18O_equilibrium"]), df_wwc["Dp17O_equilibrium"],
# #          marker="x", color="gray")
# # plt.errorbar(prime(df_wwc["d18O_equilibrium"]), df_wwc["Dp17O_equilibrium"],
# #                 xerr=df_wwc["d18O_equilibrium_err"], yerr=df_wwc["Dp17O_equilibrium_err"],
# #                 fmt="none", color="gray", elinewidth=0.8)
# # plt.plot(prime(df_wwc["d18O_AC"]), df_wwc["Dp17O_AC"],
# #         marker="x", color="red", label="Sample")
# # plt.errorbar(prime(df_wwc["d18O_AC"]), df_wwc["Dp17O_AC"],
# #                 xerr=df_wwc["d18O_error"]/np.sqrt(df_wwc["Replicates"]), yerr=df_wwc["Dp17O_error"]/np.sqrt(df_wwc["Replicates"]),
# #                 fmt="none", color="red", elinewidth=0.8)
# for i in range(len(df_wwc)):
#     theta_wwc_lst = []
#     for _ in range(monte_carlo_iterations):

#         d18O_equilibrium = np.random.normal(df_wwc["d18O_equilibrium"].iloc[i], df_wwc["d18O_equilibrium_err"].iloc[i])
#         Dp17O_equilibrium = np.random.normal(df_wwc["Dp17O_equilibrium"].iloc[i], df_wwc["Dp17O_equilibrium_err"].iloc[i])
#         d18O_measured = np.random.normal(df_wwc["d18O_AC"].iloc[i], df_wwc["d18O_error"].iloc[i])
#         Dp17O_measured = np.random.normal(df_wwc["Dp17O_AC"].iloc[i], df_wwc["Dp17O_error"].iloc[i])
        
#         # plt.plot([prime(d18O_equilibrium), prime(d18O_measured)],
#         #          [Dp17O_equilibrium, Dp17O_measured],
#         #          c="pink", ls="solid", lw=0.8, zorder=-1, alpha=0.2)

#         result_theta = calculate_theta(d18O_A=d18O_equilibrium, Dp17O_A=Dp17O_equilibrium,
#                                        d18O_B=d18O_measured, Dp17O_B=Dp17O_measured)
#         theta_wwc_lst.append(result_theta)
#     df_wwc.loc[i, "theta"] = np.mean(theta_wwc_lst)
#     df_wwc.loc[i, "theta_err"] = np.std(theta_wwc_lst)
# df_wwc = df_wwc.loc[:, ["SampleName", "theta", "theta_err"]]

# # plt.savefig(sys.path[0] + "/SK Figure MC test.png", bbox_inches='tight')
# # plt.close()

# theta_wwc = np.mean(df_wwc["theta"])
# theta_wwc_err = np.std(df_wwc["theta"])
# print(f'The mean effective theta for warm-water corals is {theta_wwc:.3f}(±{theta_wwc_err:.3f})')

# # merge the two dataframes with df on SampleName
# df = pd.merge(df, df_cwc, on="SampleName", how="left")
# df = pd.merge(df, df_wwc, on="SampleName", how="left")
# df['theta'] = df['theta_x'].combine_first(df['theta_y'])
# df['theta_err'] = df['theta_err_x'].combine_first(df['theta_err_y'])
# df = df.drop(columns=['theta_x', 'theta_y', 'theta_err_x', 'theta_err_y'])
# theta_coral = np.mean(df["theta"])
# theta_coral_err = np.std(df["theta"])
# print(f'The mean effective theta for coral vital effects is {theta_coral:.3f}(±{theta_coral_err:.3f})')


# Calculate the effective theta for coral vital effects
df["theta_coral"] = calculate_theta(d18O_A=df["d18O_equilibrium"], Dp17O_A=df["Dp17O_equilibrium"],
                                 d18O_B=df["d18O_AC"], Dp17O_B=df["Dp17O_AC"])
theta_coral = np.median(df["theta_coral"])
theta_coral_err = np.std(df["theta_coral"])
print(f'The median effective theta for coral vital effects is {theta_coral:.3f}(±{theta_coral_err:.3f})')
print(f'The effective theta range is {np.min(df["theta_coral"]):.3f} to {np.max(df["theta_coral"]):.3f}')
print("\n")

theta_wwc = df[df["Type"] == "warm-water coral"]["theta_coral"].median()
theta_wwc_err = df[df["Type"] == "warm-water coral"]["theta_coral"].std()
print(f'The median effective theta for warm-water corals is {theta_wwc:.3f}(±{theta_wwc_err:.3f})')

theta_cwc = df[df["Type"] == "cold-water coral"]["theta_coral"].median()
theta_cwc_err = df[df["Type"] == "cold-water coral"]["theta_coral"].std()
print(f'The median effective theta for cold-water corals is {theta_cwc:.3f}(±{theta_cwc_err:.3f})')


# Assign colors and markers
cat1 = df["Species"].unique()
markers = dict(zip(cat1, ["o", "s", "D", "v", "^", "<", ">", "p", "P", "*"]))
cat2 = df["Type"].unique()
colors = dict(zip(cat2, ["#1455C0", "#EC0016"]))


# CREATE FIGURE 2

# Subplot A: Dp17O vs dp18O
fig, (ax1, ax2) = plt.subplots(1, 2)

# Create a separate scatter plot for each species
for cat in cat1:
    for dog in cat2:
        data = df[(df["Species"] == cat) & (df["Type"] == dog)]
        if len(data) > 0:
            x = prime(data["d18O_AC"])
            y = data["Dp17O_AC"]
            xerr = data["d18O_error"] #/np.sqrt(data["Replicates"])
            yerr = data["Dp17O_error"] #/np.sqrt(data["Replicates"])
            ax1.scatter(x, y,
                        marker=markers[cat], fc=colors[dog], label=f"{cat}")
            ax1.errorbar(x, y, xerr=xerr, yerr=yerr,
                            fmt="none", color=colors[dog], zorder=0)
    
# Plot quilibrium points
ax1.scatter(prime(df["d18O_equilibrium"]), df["Dp17O_equilibrium"],
            marker="x", color="#747067", label="Equilibrium", zorder=10, lw=1.5)
ax1.errorbar(prime(df["d18O_equilibrium"]), df["Dp17O_equilibrium"],
             xerr=df["d18O_equilibrium_err"], yerr=df["Dp17O_equilibrium_err"],
             fmt="none", color="#747067", zorder=-1)

# Connect SampleName points with equilibrium points
for i in range(len(df)):
    ax1.plot([prime(df["d18O_AC"][i]), prime(df["d18O_equilibrium"][i])],
             [df["Dp17O_AC"][i], df["Dp17O_equilibrium"][i]],
             c="#cacaca", ls="dashed", lw=0.8, zorder=-1)

ax1.text(0.08, 0.98, "a", size=14, ha="right", va="top",
         transform=ax1.transAxes, fontweight="bold")

ax1.set_ylabel("$\Delta^{\prime 17}$O (ppm)")
ax1.set_xlabel("$\delta^{\prime 18}$O (‰, VSMOW)")

# Subplot B: Disequilibrium Dp17O and d18O values

for cat in cat1:
    for dog in cat2:
        data = df[(df["Species"] == cat) & (df["Type"] == dog)]
        if len(data) > 0:
            x = data["d18O_offset"]
            y = data["Dp17O_offset"]
            xerr = data["d18O_error"] #/np.sqrt(data["Replicates"])
            yerr = data["Dp17O_error"] #/np.sqrt(data["Replicates"])
            ax2.scatter(x, y,
                        marker=markers[cat], fc=colors[dog], label=cat)
            for xi, yi, xerri, yerri in zip(x, y, xerr, yerr):
                if not (np.isnan(xerri) or np.isnan(yerri)):
                    ax2.errorbar(x, y,
                                 xerr=xerr, yerr=yerr,
                                 fmt="none", color=colors[dog], zorder=0)

# Plot segments to equilibrium
slopes = []
for i in range(len(df)):
    x = df["d18O_offset"][i]
    y = df["Dp17O_offset"][i]
    slope = (0 - y) / (0 - x)
    slopes.append(round(slope, 1))
    ax2.plot([x, 0], [y, 0],
             color="#cacaca", ls="dashed", lw=0.8, zorder=-1)

# Plot equilibrium point
mean_xerr = np.mean(df["d18O_equilibrium_err"])
mean_yerr = np.mean(df["Dp17O_equilibrium_err"])
ax2.scatter(0, 0,
            marker="x", color="#747067", lw=1.5, label="Equilibrium", zorder=10)
ax2.errorbar(0, 0,
             xerr=mean_xerr, yerr=mean_yerr,
             fmt="none", color="#747067", zorder=-1)

# Add isoDIC models
cwc_target = (isoDIC_cwc["time(s)"] - 15*60).abs().idxmin()
wwc_target = (isoDIC_wwc["time(s)"] - 15*60).abs().idxmin()

ax2.plot(isoDIC_wwc["d18_CO3"]-isoDIC_wwc["d18_CO3"].iloc[0], isoDIC_wwc["D17_CO3"]-isoDIC_wwc["D17_CO3"].iloc[0],
         c="darkred", ls="solid", zorder=3, lw=1)
ax2.plot(isoDIC_cwc["d18_CO3"]-isoDIC_cwc["d18_CO3"].iloc[0], isoDIC_cwc["D17_CO3"]-isoDIC_cwc["D17_CO3"].iloc[0],
         c="darkblue", ls="solid", zorder=3, lw=1)
# ax2.plot(isoDIC_wwc["d18_CO3"].iloc[wwc_target]-isoDIC_wwc["d18_CO3"].iloc[0], isoDIC_wwc["D17_CO3"].iloc[wwc_target]-isoDIC_wwc["D17_CO3"].iloc[0],
#          c="darkred", marker="|", zorder=3)
# ax2.plot(isoDIC_cwc["d18_CO3"].iloc[cwc_target]-isoDIC_cwc["d18_CO3"].iloc[0], isoDIC_cwc["D17_CO3"].iloc[cwc_target]-isoDIC_cwc["D17_CO3"].iloc[0],
#          c="darkblue", marker="|", zorder=3)
ax2.annotate("",
             (isoDIC_wwc["d18_CO3"].iloc[wwc_target]-isoDIC_wwc["d18_CO3"].iloc[0],
              isoDIC_wwc["D17_CO3"].iloc[wwc_target]-isoDIC_wwc["D17_CO3"].iloc[0]),
             (isoDIC_wwc["d18_CO3"].iloc[wwc_target-1]-isoDIC_wwc["d18_CO3"].iloc[0],
              isoDIC_wwc["D17_CO3"].iloc[wwc_target-1]-isoDIC_wwc["D17_CO3"].iloc[0]),
             ha="center", va="center", zorder=-1,
             arrowprops=dict(arrowstyle="->", color="darkred", lw=1))
ax2.annotate("",
             (isoDIC_cwc["d18_CO3"].iloc[cwc_target]-isoDIC_cwc["d18_CO3"].iloc[0],
              isoDIC_cwc["D17_CO3"].iloc[cwc_target]-isoDIC_cwc["D17_CO3"].iloc[0]),
             (isoDIC_cwc["d18_CO3"].iloc[cwc_target-1]-isoDIC_cwc["d18_CO3"].iloc[0],
              isoDIC_cwc["D17_CO3"].iloc[cwc_target-1]-isoDIC_cwc["D17_CO3"].iloc[0]),
             ha="center", va="center", zorder=-1,
             arrowprops=dict(arrowstyle="->", color="darkblue", lw=1))

ax2.text(-1, -30,
         "CO$_2$ absorbtion",
         ha="center", va="center", color="k")
ax2.text(-1, -35,
         r"$\it{T}$ = 9 °C, pH = 8.8",
         ha="center", va="center", color="darkblue")
ax2.text(-1, -40,
         r"$\it{T}$ = 27 °C, pH = 8.5",
         ha="center", va="center", color="darkred")

# Add diffusion vector
theta_diff = (np.log((12+16+16)/(12+17+16)))/(np.log((12+16+16)/(12+18+16)))
shift_d18O = -1
A = (0, 0)
B = (shift_d18O, apply_theta(0, 0, shift_d18O=shift_d18O, theta=theta_diff))
ax2.annotate("",
             (A[0], A[1]),
             xytext=(B[0], B[1]),
             ha="center", va="center", zorder=-1,
             arrowprops=dict(arrowstyle="<|-", color="#FF7A00", lw=1))
ax2.text(B[0], B[1]-5,
         "diffusion",
         ha="right", va="center", color="#FF7A00")

# Add legend and format species names to italic
legend = ax2.legend(loc='upper right', bbox_to_anchor=(1.45, 1))
for text in legend.texts:
    text.set_fontsize(5.5)
    if 'Equilibrium' not in text.get_text() and 'vent coral' not in text.get_text():
        text.set_fontstyle('italic')

ax2.text(0.08, 0.98, "b", size=14, ha="right", va="top",
         transform=ax2.transAxes, fontweight="bold")

# Axis properties
ax2.set_ylabel("$\Delta^{\prime 17}$O$_{measured}$ - $\Delta^{\prime 17}$O$_{expected}$ (ppm)")
ax2.set_xlabel("$\delta^{18}$O$_{measured}$ - $\delta^{18}$O$_{expected}$ (‰, VSMOW)")

plt.xlim(-7.5, 0.5)
plt.ylim(-70, 25)

plt.savefig(sys.path[0] + "/SK Figure 2.png")
plt.close()

# Save data to CSV file
df.to_csv(sys.path[0] + "/SK Table S-3 part-3.csv", index=False)


# CREATE FIGURE 1
fig, ax = plt.subplots(1, 1, figsize=(4, 4))

# Plot seawater
plt.text(4, 0, "seawater", ha="left", va="center", c="#004B6D")
plt.scatter(prime(swdf["d18O"]), swdf["Dp17O"],
            marker="o", fc="#004B6D", linewidth=0,
            label="seawater")

# Plot equilibrium and annotate the temperature range
df_eq = plot_equilibrium(Dp17Ow=-11, d18Ow=1,
                         Tmin=0, Tmax=100,
                         ax=ax, fluid_name="seawater", color="k", highlight=True)
ax.annotate(str(round(df_eq.iloc[0, 2])) + " °C",
            (prime(df_eq.iloc[0, 0]), df_eq.iloc[0, 1]),
            xytext=(prime(df_eq.iloc[0, 0]+0.3), df_eq.iloc[0, 1]+10),
            ha="center", va="bottom", c="k",
            arrowprops=dict(arrowstyle="-|>", color="k"))
ax.annotate(str(round(df_eq.iloc[-1, 2])) + " °C",
            (prime(df_eq.iloc[-1, 0]), df_eq.iloc[-1, 1]),
            xytext=(prime(df_eq.iloc[-1, 0]+0.3), df_eq.iloc[-1, 1]+10),
            ha="center", va="bottom", c="k",
            arrowprops=dict(arrowstyle="-|>", color="k"))

# Mark the warm-water coral apparent and growth temperatures
mean_d18O = df[df["Type"] == "warm-water coral"]["d18O_AC"].mean()
df_eq_d18O = df_eq.iloc[(df_eq["d18O"]-mean_d18O).abs().argsort()[:1]]
ax.annotate(r"$\it{T}$ from $\delta^{18}O$" + "\n(" + str(round(df_eq_d18O.iloc[0, 2])) + " °C)",
            (prime(df_eq_d18O.iloc[0, 0]), df_eq_d18O.iloc[0, 1]),
            (prime(df_eq_d18O.iloc[0, 0]+0.3), df_eq_d18O.iloc[0, 1]+20),
            ha = "left", va = "bottom",
            arrowprops=dict(arrowstyle="-|>", color="#4F4B41"))
ax.annotate("",
            (prime(mean_d18O), df_eq_d18O.iloc[0, 1]),
            (prime(mean_d18O), -130),
            zorder=-1,
            arrowprops=dict(arrowstyle="-", ls="dashed", color="#4F4B41", lw=1))

mean_T = df[df["Type"] == "warm-water coral"]["T_modeled"].mean()
df_eq_T = df_eq.iloc[(df_eq["temperature"]-mean_T).abs().argsort()[:1]]
ax.annotate("growth $\it{T}$" + "\n(" + str(round(df_eq_T.iloc[0, 2])) + " °C)",
            (prime(df_eq_T.iloc[0, 0]), df_eq_T.iloc[0, 1]),
            xytext=(prime(df_eq_T.iloc[0, 0]+0.3), df_eq_T.iloc[0, 1]+20),
            ha="left", va="top",
            arrowprops=dict(arrowstyle="-|>", color="#4F4B41"))

# Vital effect arrow
shift_d18O = -10
A = (prime(df_eq_T.iloc[0, 0]),
     df_eq_T.iloc[0, 1])
B = (prime(df_eq_T.iloc[0, 0]+shift_d18O),
     apply_theta(df_eq_T.iloc[0, 0], df_eq_T.iloc[0, 1], shift_d18O=shift_d18O, theta=theta_coral))
ax.annotate("",
            (A[0], A[1]),
            (B[0], B[1]),
            ha="center", va="center", zorder=-1,
            arrowprops=dict(arrowstyle="<|-", color="#FF7A00", lw=1.5))
ax.text(B[0], B[1], '"vital effects"\n'+ r'($\theta_{effective}$ = '+ f'{theta_coral:.3f})',
        ha="right", va="center", color="#FF7A00")

# Create a separate scatter plot for each species
for cat in cat1:
    for dog in cat2:
        data = df[(df["Species"] == cat) & (df["Type"] == dog)]
        if len(data) > 0:
            x = prime(data["d18O_AC"])
            y = data["Dp17O_AC"]
            xerr = data["d18O_error"] #/np.sqrt(data["Replicates"])
            yerr = data["Dp17O_error"] #/np.sqrt(data["Replicates"])
            ax.scatter(x, y,
                       marker=markers[cat], fc=colors[dog], label=f"{cat}")
            ax.errorbar(x, y, xerr=xerr, yerr=yerr,
                        fmt="none", color=colors[dog], zorder=0)

# Label the coral groups
ax.text(33, -121, "cold-water\ncorals",
        ha="left", va="center", color="#1455C0")
ax.text(25, -80, "warm-water\ncorals",
        ha="right", va="center", color="#EC0016")

# Axis properties
plt.xlim(-2, 42)
plt.ylim(-130, 10)
plt.ylabel("$\Delta^{\prime 17}$O (ppm)")
plt.xlabel("$\delta^{\prime 18}$O (‰, VSMOW)")

plt.savefig(sys.path[0] + "/SK Figure 1.png")
plt.close()


# Abstract figure

df = df[df["SampleName"] == "SK-SA5"]

plt.rcParams.update({'font.size': 16})
plt.rcParams["figure.figsize"] = (5, 3)
plt.rcParams["lines.linewidth"] = 1  # error bar width
plt.rcParams["patch.linewidth"] = 1  # marker edge width

fig = plt.figure()

ax = fig.add_subplot(111)
ax.set_facecolor('none')

# change axis colors to white
ax.spines['bottom'].set_color('w')
ax.spines['top'].set_color('w')
ax.spines['left'].set_color('w')
ax.spines['right'].set_color('w')
ax.xaxis.label.set_color('w')
ax.yaxis.label.set_color('w')
ax.tick_params(axis='x', colors='w')
ax.tick_params(axis='y', colors='w')

# Plot equilibrium and annotate the temperature range
df_eq = plot_equilibrium(Dp17Ow=-11, d18Ow=1,
                         Tmin=0, Tmax=100,
                         ax=ax, fluid_name="seawater", color="w", highlight=False, mark_water=False)

plt.text(df_eq["d18O"].iloc[-1], df_eq["Dp17O"].iloc[-1], "equilibrium\n(0–100 °C)", ha="right", va="center", c="w")

# Vital effect arrow
shift_d18O = -10
A = (prime(df_eq_T.iloc[0, 0]),
     df_eq_T.iloc[0, 1])
B = (prime(df_eq_T.iloc[0, 0]+shift_d18O),
     apply_theta(df_eq_T.iloc[0, 0], df_eq_T.iloc[0, 1], shift_d18O=shift_d18O, theta=theta_coral))
ax.annotate("",
            (A[0], A[1]),
            (B[0], B[1]),
            ha="center", va="center", zorder=-1,
            arrowprops=dict(arrowstyle="<|-", color="w", lw=1.5))
ax.text(B[0], B[1], r"$\bf{coral}$" +"\n"+ r"$\bf{'vital}$ $\bf{effects'}$" + "\n" + f"($\\theta_{{coral}}$ = {theta_coral:.3f})",
        ha="right", va="center", color="w")

# Mark the warm-water coral apparent and growth temperatures
mean_d18O = df[df["Type"] == "warm-water coral"]["d18O_AC"].mean()
df_eq_d18O = df_eq.iloc[(df_eq["d18O"]-mean_d18O).abs().argsort()[:1]]
ax.annotate(r"$\it{T}$ from $\delta^{18}O$" + "\n(" + str(round(df_eq_d18O.iloc[0, 2])) + " °C)",
            (prime(df_eq_d18O.iloc[0, 0]), df_eq_d18O.iloc[0, 1]),
            (prime(df_eq_d18O.iloc[0, 0])-3, df_eq_d18O.iloc[0, 1]+37),
            ha = "left", va = "bottom", color = "w",
            arrowprops=dict(arrowstyle="-|>", color="w"))
ax.annotate("",
            (prime(mean_d18O), df_eq_d18O.iloc[0, 1]),
            (prime(mean_d18O), -130),
            zorder=-1,
            arrowprops=dict(arrowstyle="-", ls="dashed", color="w", lw=1))
plt.plot(prime(df_eq_d18O.iloc[0, 0]), df_eq_d18O.iloc[0, 1],
         marker=".", c="w", mec="k", mew=0.5, zorder=10)

mean_T = df[df["Type"] == "warm-water coral"]["T_modeled"].mean()
df_eq_T = df_eq.iloc[(df_eq["temperature"]-mean_T).abs().argsort()[:1]]
ax.annotate("growth $\it{T}$" + "\n(" + str(round(df_eq_T.iloc[0, 2])) + " °C)",
            (prime(df_eq_T.iloc[0, 0]), df_eq_T.iloc[0, 1]),
            xytext=(prime(df_eq_T.iloc[0, 0])-2, df_eq_T.iloc[0, 1]+45),
            ha="left", va="top", color = "w",
            arrowprops=dict(arrowstyle="-|>", color="w"))
plt.plot(prime(df_eq_T.iloc[0, 0]), df_eq_T.iloc[0, 1],
            marker=".", c="w", mec="k", mew=0.5, zorder=10)

x = prime(df["d18O_AC"])
y = df["Dp17O_AC"]
xerr = df["d18O_error"]
yerr = df["Dp17O_error"]
ax.scatter(x, y,
           marker="D", fc="#EC0016", ec="w", s=50)
plt.errorbar(x, y,
             xerr=xerr, yerr=yerr,
             fmt="none", color="#EC0016", zorder=2)

# Label the coral
ax.text(np.mean(x)+1, np.mean(y)-5, "coral", ha="left", va="top", color="w")

# Axis properties
plt.xlim(-2, 42)
plt.ylim(-130, 10)
plt.xticks([])
plt.yticks([])
plt.ylabel("$\Delta^{\prime 17}$O")
plt.xlabel("$\delta^{\prime 18}$O")

# Insert background image
im = plt.imread(sys.path[0] + "/background.jpeg")
y_scale = fig.get_size_inches()[1]*600 / im.shape[0]
x_scale = fig.get_size_inches()[0]*600 / im.shape[1]
resized_im = zoom(im, (y_scale, x_scale, 1))
fig.figimage(resized_im, zorder = -12)

plt.savefig(sys.path[0] + "/SK Graphical Abstract.png")
plt.close("all")