# Showing the concept of the paper
# INPUT: SK Table S-3 part-2.csv, seawater.csv, background.jpeg
# OUTPUT: SK Figure 2.png, SK Graphical Abstract.png

# >>>>>>>>>

# Import libraries
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import zoom
from scipy.spatial.distance import cdist

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


def d18O_cc(equilibrium_temperatures, d18Ow):
    return a18_cc(equilibrium_temperatures + 273.15) * (d18Ow+1000) - 1000


def d17O_cc(equilibrium_temperatures, d17Ow):
    return a17_cc(equilibrium_temperatures + 273.15) * (d17Ow+1000) - 1000


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


def plot_equilibrium(Dp17Ow, d18Ow, Tmin, Tmax, ax, color="k", highlight=True, mark_water=True, plot=True):

    d17Ow = unprime(0.528 * prime(d18Ow) + Dp17Ow/1000)

    # mark water
    if mark_water == True and plot == True:
        ax.scatter(prime(d18Ow), Dp17O(d17Ow, d18Ow),
                   marker="X", fc=color, ec=None,
                   zorder=10)

    # equilibrium line, entire T range
    toInf = np.arange(Tmin, Tmax, 1)
    d18O_mineral = d18O_cc(toInf, d18Ow)
    d17O_mineral = d17O_cc(toInf, d17Ow)
    mineral_equilibrium = np.array(
        [d18O_mineral, Dp17O(d17O_mineral, d18O_mineral), toInf]).T
    if plot == True:
        ax.plot(prime(mineral_equilibrium[:, 0]), mineral_equilibrium[:, 1],
                ":", c=color, zorder=3)

    # equilibrium points, highlight range
    equilibrium_temperatures = np.arange(Tmin, Tmax, 0.5)
    colors = np.linspace(0, 1, len(equilibrium_temperatures))
    d18O_mineral = d18O_cc(equilibrium_temperatures, d18Ow)
    d17O_mineral = d17O_cc(equilibrium_temperatures, d17Ow)
    mineral_equilibrium = np.array([d18O_mineral, Dp17O(
        d17O_mineral, d18O_mineral), equilibrium_temperatures]).T
    if highlight == True and plot == True:
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

    # Get equilibrium values for seawater
    df_eq = plot_equilibrium(Dp17Ow=Dp17O_seawater, d18Ow=d18O_seawater,
                             Tmin=0, Tmax=100,
                             ax=ax, color="k", plot=False)
    ax.errorbar(prime(d18O_seawater), Dp17O_seawater,
                xerr=d18O_seawater_err, yerr=Dp17O_seawater_err,
                fmt="none", color="k", zorder=-1)

    # Get the vital effect line
    shift_d18O = 5.9
    df_line = vital_vector(d18O_coral, Dp17O_coral, len(df_eq["d18O"]), shift_d18O, theta_coral)

    # Calculate the distance between points
    df_eq_f = df_eq.loc[:, ["d18O", "Dp17O"]]
    distances = cdist(df_eq_f, df_line)
    min_indices = np.unravel_index(np.argmin(distances), distances.shape)

    temp = df_eq.iloc[min_indices[0]]["temperature"]

    # Plot
    ax.annotate("", xy=(prime(df_line["d18O"].values[0]), df_line["Dp17O"].values[0]),
                xytext=(prime(df_line["d18O"].values[-1]),
                        df_line["Dp17O"].values[-1]),
                arrowprops=dict(arrowstyle="<|-", color="#FF7A00", lw=1.5))

    # Maximum temperature
    shift_d18O = 5
    df_eq = plot_equilibrium(Dp17Ow=Dp17O_seawater, d18Ow=d18O_seawater,
                             Tmin=10, Tmax=100,
                             ax=ax, color="red", plot=False)
    df_line = vital_vector(d18O_coral, Dp17O_coral + Dp17O_coral_error,
                           len(df_eq["d18O"]), shift_d18O, theta_coral)
    df_eq_f = df_eq.loc[:, ["d18O", "Dp17O"]]
    distances = cdist(df_eq_f, df_line)
    min_indices = np.unravel_index(np.argmin(distances), distances.shape)
    T_max = df_eq.iloc[min_indices[0]]["temperature"]
    ax.plot(prime(df_line["d18O"]), df_line["Dp17O"], "#FF7A00", ls=":")

    # Minimum temperature
    shift_d18O = 6.5
    df_eq = plot_equilibrium(Dp17Ow=Dp17O_seawater, d18Ow=d18O_seawater,
                             Tmin=10, Tmax=100,
                             ax=ax, color="blue", plot=False)
    df_line = vital_vector(d18O_coral + d18O_coral_error, Dp17O_coral -
                           Dp17O_coral_error, len(df_eq["d18O"]), shift_d18O, theta_coral)
    df_eq_f = df_eq.loc[:, ["d18O", "Dp17O"]]
    distances = cdist(df_eq_f, df_line)
    min_indices = np.unravel_index(np.argmin(distances), distances.shape)
    T_min = df_eq.iloc[min_indices[0]]["temperature"]
    ax.plot(prime(df_line["d18O"]), df_line["Dp17O"], "#FF7A00", ls=":")

    print(f"Reconstructed temperature: {temp:.1f} ±{(T_max-T_min)/2:.1f} °C")
    return temp, (T_max-T_min)/2


# Read data from CSV files
swdf = pd.read_csv(sys.path[0] + "/seawater.csv", sep=",")
dfAll = pd.read_csv(sys.path[0] + "/SK Table S-3.csv", sep=",")
theta_coral = round(dfAll["theta_coral"].mean(), 3)
print(f"theta_coral: {theta_coral:.3f}")
df = dfAll[dfAll["SampleName"] == "SK-SA5"]

# CREATE FIGURE 2
fig, ax = plt.subplots(1, 1)

# Plot seawater
plt.text(4, 0, "seawater", ha="left", va="center", c="#83CACA")
plt.scatter(prime(swdf["d18O"]), swdf["Dp17O"],
            marker="o", fc="#83CACA", lw=0,
            label="seawater", zorder=-1)

# Plot equilibrium and annotate the temperature range
df_eq = plot_equilibrium(Dp17Ow=df["Dp17Osw"].mean(), d18Ow=df["d18Osw_database"].mean(),
                         Tmin=10, Tmax=100,
                         ax=ax, color="k", highlight=True)

ax.annotate("Carbonate equilibrium\n(10–100 °C)",
            (prime(df_eq.iloc[-1, 0]), df_eq.iloc[-1, 1]),
            xytext=(prime(df_eq.iloc[-1, 0]-3), df_eq.iloc[-1, 1]),
            ha="right", va="center", c="k",
            arrowprops=dict(arrowstyle="->", color="k"))

# Mark the warm-water coral apparent and growth temperatures
mean_d18O = df[df["Type"] == "warm-water coral"]["d18O_AC"].mean()
mean_Dp17O = df[df["Type"] == "warm-water coral"]["Dp17O_AC"].mean()
df_eq_d18O = df_eq.iloc[(df_eq["d18O"]-mean_d18O).abs().argsort()[:1]]
ax.scatter(prime(df_eq_d18O.iloc[0, 0]), df_eq_d18O.iloc[0, 1],
           marker=".", fc="#4F4B41", ec="k", zorder=10)
ax.annotate("apparent $\it{T}$\nfrom $\delta^{18}O$\n(" + f"{df_eq_d18O.iloc[0, 2]:.0f} °C)",
            (prime(df_eq_d18O.iloc[0, 0]), df_eq_d18O.iloc[0, 1]),
            (prime(df_eq_d18O.iloc[0, 0]+0.3), df_eq_d18O.iloc[0, 1]+20),
            ha = "left", va = "bottom",
            arrowprops=dict(arrowstyle="->", color="k"))
ax.annotate("",
            (prime(df_eq_d18O.iloc[0, 0]), df_eq_d18O.iloc[0, 1]),
            (prime(mean_d18O), mean_Dp17O),
            zorder=-1,
            arrowprops=dict(arrowstyle="-|>", color="#4F4B41", lw=1.5))

mean_T = df[df["Type"] == "warm-water coral"]["T_database"].mean()
df_eq_T = df_eq.iloc[(df_eq["temperature"]-mean_T).abs().argsort()[:1]]
ax.scatter(prime(df_eq_T.iloc[0, 0]), df_eq_T.iloc[0, 1],
           marker=".", fc="w", ec="k", zorder=10)
ax.annotate("measured\ngrowth $\it{T}$" + "\n(" + str(round(df_eq_T.iloc[0, 2])) + " °C)",
            (prime(df_eq_T.iloc[0, 0]), df_eq_T.iloc[0, 1]),
            xytext=(prime(df_eq_T.iloc[0, 0]+0.3), df_eq_T.iloc[0, 1]+30),
            ha="left", va="top",
            arrowprops=dict(arrowstyle="->", color="k"))


rT, rTerr = get_17O_temp(d18O_coral=df["d18O_AC"].mean(),
                                Dp17O_coral=df["Dp17O_AC"].mean(),
                                d18O_coral_error=df["d18O_error"].mean(),
                                Dp17O_coral_error=df["Dp17O_error"].mean(),
                                d18O_seawater=df["d18Osw_database"].mean(),
                                d18O_seawater_err=df["d18Osw_database_err"].mean(),
                                Dp17O_seawater=df["Dp17Osw"].mean(),
                                Dp17O_seawater_err=df["Dp17Osw_err"].mean(),
                                theta_coral=theta_coral.mean(),
                                ax=ax)

df_eq_T = df_eq.iloc[(df_eq["temperature"]-rT).abs().argsort()[:1]]
ax.scatter(prime(df_eq_T.iloc[0, 0]), df_eq_T.iloc[0, 1],
           marker=".", fc="#FF7A00", ec="k", zorder=10)
ax.annotate("reconstructed\ngrowth $\it{T}$" + "\n(" + f"{rT:.0f}±{rTerr:.0f} °C)",
            (prime(df_eq_T.iloc[0, 0]), df_eq_T.iloc[0, 1]),
            xytext=(prime(df_eq_T.iloc[0, 0]+4), df_eq_T.iloc[0, 1]+10),
            ha="left", va="top",
            arrowprops=dict(arrowstyle="->", color="k"))

# Assign colors and markers
cat1 = df["Species"].unique()
markers = dict(zip(cat1, ["o", "s", "D", "v", "^", "<", ">", "p", "P", "*"]))
cat2 = df["Type"].unique()
colors = dict(zip(cat2, ["#1455C0", "#EC0016"]))

# Plot the coral
x = prime(df["d18O_AC"])
y = df["Dp17O_AC"]
xerr = df["d18O_error"]
yerr = df["Dp17O_error"]
ax.scatter(x, y,
           marker="P", fc="#EC0016", ec="w", zorder=10)
ax.errorbar(x, y, xerr=xerr, yerr=yerr,
            fmt="none", color="#EC0016", zorder=9)
ax.text(np.mean(x)-1, np.mean(y), "coral",
        ha="right", va="center", c="#EC0016")

# Axis properties
plt.xlim(-2, 47)
plt.ylim(-110, 10)
plt.ylabel("$\Delta^{\prime 17}$O (ppm)")
plt.xlabel("$\delta^{\prime 18}$O (‰, VSMOW)")

plt.savefig(sys.path[0] + "/SK Figure 2.png")
plt.close()


# Abstract figure

df = dfAll[dfAll["SampleName"] == "SK-DS3"]

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
                         Tmin=10, Tmax=100,
                         ax=ax, color="w", highlight=False, mark_water=False)

plt.text(df_eq["d18O"].iloc[-1], df_eq["Dp17O"].iloc[-1], "equilibrium\n(10–100 °C)", ha="right", va="center", c="w")

# Vital effect arrow (plot already calculated values)
df_eq_T = df_eq.iloc[(df_eq["temperature"]-df["T_17O"].mean()).abs().argsort()[:1]]
shift_d18O = -10
A = (prime(df_eq_T.iloc[0, 0]),
     df_eq_T.iloc[0, 1])
B = (prime(df_eq_T.iloc[0, 0]+shift_d18O),
     apply_theta(df_eq_T.iloc[0, 0], df_eq_T.iloc[0, 1], shift_d18O=shift_d18O, theta=theta_coral))
ax.annotate("",
            (A[0], A[1]),
            (B[0], B[1]),
            ha="center", va="center", zorder=-1,
            arrowprops=dict(arrowstyle="-|>", color="w", lw=1.5))
ax.text(B[0], B[1], r"$\bf{coral}$" +"\n"+ r"$\bf{'vital}$ $\bf{effects'}$" + "\n" + f"($\\theta_{{coral}}$ = {theta_coral:.3f})",
        ha="right", va="center", color="w")

# Mark the warm-water coral apparent and growth temperatures
mean_d18O = df["d18O_AC"].mean()
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

mean_T = df[df["Type"] == "warm-water coral"]["T_database"].mean()
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