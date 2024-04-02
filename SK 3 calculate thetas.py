# This code is used to plot the coral data in triple oxygen isotope space,
# and calculate the vital effect theta

# INPUT: SK Table S-3 part-2.csv, isoDIC_***.csv
# OUTPUT: SK Figure 1.png, SK Figure S5.png, (SK Figure S6.png), SK Table S-3 part-3.csv

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
    d17O_A = d17O(d18O_A, Dp17O_A)

    a18 = (d18O_B + 1000) / (d18O_A + 1000)
    a17 = a18**theta

    d17O_B = a17 * (d17O_A + 1000) - 1000
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


# isoDIC models
isoDIC_header = pd.read_csv(sys.path[0] + "/isoDIC_header.csv", sep=",")
isoDIC_cwc = pd.read_csv(sys.path[0] + "/isoDIC_pH8.8_T9.csv", sep=",")
isoDIC_cwc.columns = isoDIC_header.columns
isoDIC_wwc = pd.read_csv(sys.path[0] + "/isoDIC_pH8.5_T27.csv", sep=",")
isoDIC_wwc.columns = isoDIC_header.columns

df = pd.read_csv(sys.path[0] + "/SK Table S-3 part-2.csv", sep=",")

# Calculate equilibrium values using the "measured + database" d18Osw and T values
df_equi = cc_equilibrium(T=df["T_database"], T_err=df["T_database_err"],
                                d18Ow=df["d18Osw_database"], d18Ow_err=df["d18Osw_database_err"],
                                Dp17Ow=df["Dp17Osw"], Dp17Ow_err=df["Dp17Osw_err"],
                                Sample=df["SampleName"])
df = pd.merge(df, df_equi, on="SampleName", how="left")
df["d18O_offset"] = df["d18O_AC"]-df["d18O_equilibrium"]
df["Dp17O_offset"] = df["Dp17O_AC"]-df["Dp17O_equilibrium"]

# Calculate the effective theta for coral vital effects
df["theta_coral"] = calculate_theta(d18O_A=df["d18O_equilibrium"], Dp17O_A=df["Dp17O_equilibrium"],
                                 d18O_B=df["d18O_AC"], Dp17O_B=df["Dp17O_AC"])
theta_coral = round(df["theta_coral"].mean(), 3)
theta_coral_std = np.std(df["theta_coral"])
print(f'The median effective theta for coral vital effects is {theta_coral:.3f}(±{theta_coral_std:.3f})')
print(f'The effective theta range is {np.min(df["theta_coral"]):.3f} to {np.max(df["theta_coral"]):.3f}')

# Calculate the error of the effective theta for coral vital effects
def monte_carlo_simulation_row(row, num_simulations=1000):
    thetas = []
    for _ in range(num_simulations):
        d18O_A_error = np.random.normal(loc=0, scale=row["d18O_error"])
        Dp17O_A_error = np.random.normal(loc=0, scale=row["Dp17O_error"])
        d18O_B_error = np.random.normal(loc=0, scale=row["d18O_equilibrium_err"])
        Dp17O_B_error = np.random.normal(loc=0, scale=row["Dp17O_equilibrium_err"])
        
        d18O_A_with_error = row["d18O_equilibrium"] + d18O_A_error
        Dp17O_A_with_error = row["Dp17O_equilibrium"] + Dp17O_A_error
        d18O_B_with_error = row["d18O_AC"] + d18O_B_error
        Dp17O_B_with_error = row["Dp17O_AC"] + Dp17O_B_error
        
        theta_coral_with_error = calculate_theta(d18O_A_with_error, Dp17O_A_with_error, 
                                                 d18O_B_with_error, Dp17O_B_with_error)
        thetas.append(theta_coral_with_error)
    
    return round(np.std(thetas), 3)

df["theta_coral_error"] = df.apply(monte_carlo_simulation_row, axis=1)
print(f'The error of the individual coral vital effects theta is {df["theta_coral_error"].mean():.3f}')

wwc_df = df[df["Type"] == "warm-water coral"]
theta_wwc = wwc_df["theta_coral"].mean()
theta_wwc_err = wwc_df["theta_coral"].std()
print(f'The mean effective theta for warm-water corals is {theta_wwc:.3f}(±{theta_wwc_err:.3f})')

cwc_df = df[df["Type"] == "cold-water coral"]
theta_cwc = cwc_df["theta_coral"].mean()
theta_cwc_err = cwc_df["theta_coral"].std()
print(f'The mean effective theta for cold-water corals is {theta_cwc:.3f}(±{theta_cwc_err:.3f})')


# Assign colors and markers
cat1 = df["Species"].unique()
markers = dict(zip(cat1, ["o", "s", "D", "v", "^", "<", ">", "p", "P", "*"]))
cat2 = df["Type"].unique()
colors = dict(zip(cat2, ["#1455C0", "#EC0016"]))



# CREATE FIGURE S5

fig, ax = plt.subplots()

# Create a separate scatter plot for each species
for cat in cat1:
    for dog in cat2:
        data = df[(df["Species"] == cat) & (df["Type"] == dog)]
        if len(data) > 0:
            x = data["SampleName"]
            y = data["theta_coral"]
            yerr = data["theta_coral_error"]
            ax.scatter(x, y,
                        marker=markers[cat], fc=colors[dog], label=f"{cat}", zorder=10)
            ax.errorbar(x, y, yerr=yerr,
                            fmt="none", color=colors[dog], zorder=0)
            for xi, yi in zip(x, y):
                ax.text(xi, 0.5265, str(xi), ha="center", va="center", c = colors[dog], fontsize=5)

f = 0.5
ax.set_xlim(0-f, len(df["SampleName"])-1+f)
x_min, x_max = ax.get_xlim()

x_values = np.linspace(x_min, x_max, 2)
ax.axhline(theta_cwc, color="blue", ls="-")
ax.fill_between(x_values, theta_cwc - theta_cwc_err, theta_cwc + theta_cwc_err, color="blue", alpha=0.2, lw = 0)
ax.axhline(theta_wwc, color="red", ls="-")
ax.fill_between(x_values, theta_wwc - theta_wwc_err, theta_wwc + theta_wwc_err, color="red", alpha=0.2, lw = 0)

ax.set_ylabel("$\\theta_{coral}$")
ax.set_xticks([])

plt.savefig(sys.path[0] + "/SK Figure S5.png")
plt.close()



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
         "CO$_2$ absorbtion\n(modelled)",
         ha="center", va="center", color="k")
ax2.text(-1, -37,
         r"$\it{T}$ = 9 °C, pH = 8.8",
         ha="center", va="center", color="darkblue")
ax2.text(-1, -42,
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
             arrowprops=dict(arrowstyle="<|-", color="#63A615", lw=1))
ax2.text(B[0], B[1]-5,
         "diffusion",
         ha="right", va="center", color="#63A615")


# Add revised CO2 absorption vector
theta_diff = 0.532
shift_d18O = -2
A = (0, 0)
B = (shift_d18O, apply_theta(0, 0, shift_d18O=shift_d18O, theta=theta_diff))
ax2.annotate("",
             (A[0], A[1]),
             xytext=(B[0], B[1]),
             ha="center", va="center", zorder=3,
             arrowprops=dict(arrowstyle="<|-", color="#FF7A00", lw=1))
ax2.text(B[0]+0.55, B[1]+6,
         "CO$_2$ absorbtion\n(experimental)",
         bbox=dict(fc='white', ec="None", alpha=0.5, pad=0.1),
         ha="center", va="bottom", color="#FF7A00")


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

plt.savefig(sys.path[0] + "/SK Figure 1.png")
plt.close()

# Save data to CSV file
df.to_csv(sys.path[0] + "/SK Table S-3 part-3.csv", index=False)


# Do some additional calculations for the discussion
coral_slope_range = [df["theta_coral"].min().round(3), df["theta_coral"].max().round(3)]

# Calculate diffusion-induced 'vital effect' percentage
abs_values = [0.538, 0.541] # from Guo and Zhou (2019)
diffusion_slope = (np.log((12+16+16)/(12+17+16)))/(np.log((12+16+16)/(12+18+16)))
for abs in abs_values:
    coral_slope_min = min(coral_slope_range)
    x = (1 - ((coral_slope_min - diffusion_slope) / (abs - diffusion_slope))) * 100
    print(f"If the absorption slope is {abs} and the coral theta is {coral_slope_min},\nthen up to {x:.0f}% of the total ‘vital effect’ is induced by diffusion\n")

# Calculate diffusion-induced 'vital effect' percentage using the revised theta estimates from Bajnai et al. (2023)
abs_values = [0.531, 0.532]
for abs in abs_values:
    coral_slope_min = min(coral_slope_range)
    x = (1 - ((coral_slope_min - diffusion_slope) / (abs - diffusion_slope))) * 100
    print(f"If the absorption slope is {abs} and the coral theta is {coral_slope_min},\nthen up to {x:.0f}% of the total ‘vital effect’ is induced by diffusion\n")