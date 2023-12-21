# The following code is used to scale the TILDAS data
# based on the accepted values for the reference gases

# INPUT: SK Table S-2.csv
# OUTPUT: SK Figure S1.png, SK Figure S2.png, SK Table S-3 part-1.csv

# >>>>>>>>>

# Import libraries
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


# Functions that make life easier
def prime(delta):
    return 1000 * np.log(delta/1000 + 1)


def unprime(dprime):
    return (np.exp(dprime/1000) - 1) * 1000


# This function calculates the ∆'17O value (in ppm) from the d17O and d18O values
def Dp17O(d17O, d18O):
    return (prime(d17O) - 0.528 * prime(d18O)) * 1000


# This function scales the data based on the accepted values of the light and heavy reference gases
# The scaling is done for each measurement period separately

# This function scales the data based on the accepted values of the light and heavy reference gases
# The scaling is done for each measurement period separately
def scaleData(df, project):

    df["DateTimeMeasured"] = pd.to_datetime(df["DateTimeMeasured"])

    # Perform the scaling for each meausrement period separately
    df_samples = pd.DataFrame()
    grouped = df.groupby("MeasurementPeriod")
    if grouped.ngroups == 1:
        SuppFig = [""]
    else:
        SuppFig = ["a","b","c","d"]
    FigNum = 0
    for period, group in grouped:

        print(f"Measurement period {period}:")

        heavy_d18O_SD = group[group["SampleName"].str.contains("heavy")]["d18O"].std()
        heavy_Dp17O_SD = group[group["SampleName"].str.contains("heavy")]["Dp17O"].std()
        N_heavy = len(group[group["SampleName"].str.contains("heavy")])
        print(f"Heavy reference gas N = {N_heavy}, d18O SD = ±{heavy_d18O_SD:.3f}‰, ∆'17O SD = ±{heavy_Dp17O_SD:.0f} ppm")

        light_d18O_SD = group[group["SampleName"].str.contains("light")]["d18O"].std()
        light_Dp17O_SD = group[group["SampleName"].str.contains("light")]["Dp17O"].std()
        N_light = len(group[group["SampleName"].str.contains("light")])
        print(f"Light reference gas N = {N_light}, d18O SD: ±{light_d18O_SD:.3f}‰, ∆'17O SD: ±{light_Dp17O_SD:.0f} ppm")

        # Do the scaling here, based on the accepted values of the light and heavy reference gases
        heavy_d18O_measured = group[group["SampleName"].str.contains("heavy")]["d18O"].mean()
        heavy_d17O_measured = group[group["SampleName"].str.contains("heavy")]["d17O"].mean()

        light_d18O_measured = group[group["SampleName"].str.contains("light")]["d18O"].mean()
        light_d17O_measured = group[group["SampleName"].str.contains("light")]["d17O"].mean()

        heavy_d18O_accepted = 76.820
        heavy_Dp17O_accepted = -213
        heavy_d17O_accepted = unprime(heavy_Dp17O_accepted/1000 + 0.528 * prime(heavy_d18O_accepted))

        light_d18O_accepted = -1.509
        light_Dp17O_accepted = -141
        light_d17O_accepted = unprime(light_Dp17O_accepted/1000 + 0.528 * prime(light_d18O_accepted))

        slope_d18O = (light_d18O_accepted - heavy_d18O_accepted) / (light_d18O_measured - heavy_d18O_measured)
        intercept_d18O = heavy_d18O_accepted - slope_d18O * heavy_d18O_measured

        slope_d17O = (light_d17O_accepted - heavy_d17O_accepted) / (light_d17O_measured - heavy_d17O_measured)
        intercept_d17O = heavy_d17O_accepted - slope_d17O * heavy_d17O_measured

        group["d18O_scaled"] = slope_d18O*group['d18O']+intercept_d18O
        group["d17O_scaled"] = slope_d17O*group['d17O']+intercept_d17O
        group["Dp17O_scaled"] = (prime(group["d17O_scaled"]) - 0.528 * prime(group["d18O_scaled"]))*1000

        # Assign colors and markers to samples
        categories = group["SampleName"].unique()
        markers = dict(zip(categories, ["o", "s", "D", "v", "^",
                                        "<", ">", "p", "P", "*",
                                        "o", "s", "D", "v", "^",
                                        "<", ">", "p", "P", "*"]))
        colors = dict(zip(categories, plt.cm.tab20(np.linspace(0, 1, 20))))

        # Figure: unscaled Dp17O vs time

        fig, ax = plt.subplots()

        for cat in categories:
            data = group[group["SampleName"] == cat]
            ax.scatter(data["DateTimeMeasured"], data["Dp17O"],
                       marker=markers[cat], fc=colors[cat], label=cat)
            if np.isnan(data["Dp17OError"]).any() == False:
                plt.errorbar(group["DateTimeMeasured"], group["Dp17O"],
                             yerr=group["Dp17OError"],
                             fmt="none", color="#cacaca", zorder=0)
        
        plt.title(f"Measurement period: {period}")
        plt.ylabel("$\Delta^{\prime 17}$O (ppm, unscaled CO$_2$)")
        plt.xlabel("Measurement date")
        plt.legend(loc='upper right', bbox_to_anchor=(1.15, 1))

        plt.text(0.98, 0.98, SuppFig[FigNum], size=14, ha="right", va="top",
                 transform=ax.transAxes, fontweight="bold")

        plt.savefig(sys.path[0] + "/" + f"{project} Figure S1{SuppFig[FigNum]}.png")
        plt.close()

        # Exclude light and heavy CO2 from the exported dataframe
        group = group[~group["SampleName"].str.contains("heavy")]
        group = group[~group["SampleName"].str.contains("light")]

        FigNum += 1
        df_samples = pd.concat([df_samples, group])

        print("\n")

    return df_samples


# This function averages the scaled data from the different measurement periods
# and applies and acid fractionation factor
def applyAFF(df, mineral):
    # Calculate the mean values from the replicate measurements
    df = df.loc[:, ["SampleName", "d18O_scaled", "d17O_scaled", "Dp17O_scaled"]]
    df_mean = df.groupby('SampleName').mean().reset_index()
    df_mean = df_mean.rename(columns={'d18O_scaled': 'd18O_CO2', 'd17O_scaled': 'd17O_CO2', 'Dp17O_scaled': 'Dp17O_CO2'})

    # Calculate the standard deviation from the replicate measurements
    df_std = df.groupby('SampleName').std().reset_index()
    df_std = df_std.rename(columns={'d18O_scaled': 'd18O_error', 'd17O_scaled': 'd17O_error', 'Dp17O_scaled': 'Dp17O_error'})

    dfMerged = df_mean.merge(df_std, on='SampleName')
    dfMerged['Replicates'] = df.groupby('SampleName').size().reset_index(name='counts')['counts']

    # Acid fractionation correction
    if mineral == "calcite":
        dfMerged["d18O_AC"] = (dfMerged["d18O_CO2"] + 1000) / 1.01025 - 1000
        dfMerged["d17O_AC"] = (dfMerged["d17O_CO2"] + 1000) / (1.01025 ** 0.523) - 1000
        dfMerged["Dp17O_AC"] = Dp17O(dfMerged["d17O_AC"], dfMerged["d18O_AC"])
    
    elif mineral == "aragonite":
        dfMerged["d18O_AC"] = (dfMerged["d18O_CO2"] + 1000) / 1.01063 - 1000
        dfMerged["d17O_AC"] = (dfMerged["d17O_CO2"] + 1000) / (1.01063 ** 0.523) - 1000
        dfMerged["Dp17O_AC"] = Dp17O(dfMerged["d17O_AC"], dfMerged["d18O_AC"])
    
    print("All sample replicates averaged:")
    print(df.round({"Dp17O_CO2": 0, "Dp17O_error": 0, "Dp17O_AC": 0}).round(3))

    return dfMerged


# Here we go!

# Scale the data
df = scaleData(pd.read_csv(sys.path[0] + "/SK Table S-2.csv"), "SK")

df_wAFF = applyAFF(df, "aragonite")

# Exclude carbonate references
df_wAFF = df_wAFF[~df_wAFF["SampleName"].str.contains("NBS|DH11|IAEA")].reset_index(drop=True)
df = df[~df["SampleName"].str.contains("NBS|DH11|IAEA")].reset_index(drop=True)

# Export CSV
df_wAFF.to_csv(sys.path[0] + "/SK Table S-3 part-1.csv", index=False)


# Create Figure S2
fig, (ax1, ax2) = plt.subplots(1, 2)

# Assign colors and markers to samples
df.sort_values(by="SampleName", inplace=True)
categories = df["SampleName"].unique()
markers = dict(zip(categories, ["o", "s", "D", "^", "v", "X", "P", "*", "o",
               "s", "D", "^", "v", "X", "P", "*", "o", "s", "D", "^", "v", "X", "P", "*"]))
colors = dict(zip(categories, plt.cm.tab20(
    np.linspace(0, 1, len(categories)))))


# Subplot A

for cat in categories:
    data = df[df["SampleName"] == cat]
    ax1.scatter(prime(data["d18O_scaled"]), data["Dp17O_scaled"],
                marker=markers[cat], fc=colors[cat], label=cat)
    ax1.errorbar(prime(df["d18O_scaled"]), df["Dp17O_scaled"],
                 yerr=df["Dp17OError"], xerr=df["d18OError"],
                 fmt="none", color="#cacaca", zorder=0)

ax1.text(0.98, 0.98, "a", size=14, ha="right", va="top",
         transform=ax1.transAxes, fontweight="bold")

ax1.set_ylabel("$\Delta^{\prime 17}$O (ppm, CO$_2$)")
ax1.set_xlabel("$\delta^{\prime 18}$O (‰, VSMOW, CO$_2$)")

ylim = ax1.get_ylim()
xlim = ax1.get_xlim()

# Subplot B

for cat in categories:
    data = df_wAFF[df_wAFF["SampleName"] == cat]
    ax2.scatter(prime(data["d18O_CO2"]), data["Dp17O_CO2"],
                marker=markers[cat], fc=colors[cat], label=cat)
    ax2.errorbar(prime(df_wAFF["d18O_CO2"]), df_wAFF["Dp17O_CO2"],
                 yerr=df_wAFF["Dp17O_error"], xerr=df_wAFF["d18O_error"],
                 fmt="none", color="#cacaca", zorder=0)

ax2.text(0.98, 0.98, "b", size=14, ha="right", va="top",
         transform=ax2.transAxes, fontweight="bold")

ax2.set_ylabel("$\Delta^{\prime 17}$O (ppm, CO$_2$)")
ax2.set_xlabel("$\delta^{\prime 18}$O (‰, VSMOW, CO$_2$)")

ax2.set_ylim(ylim)
ax2.set_xlim(xlim)

ax2.legend(loc='upper right', bbox_to_anchor=(1.35, 1))

plt.savefig(sys.path[0] + "/SK Figure S2.png")
plt.close("all")
