# This code is partly based on isoForam by Daëron & Gray (2023)
# The code reads in the coral collection sites (location and depth) from "SK sample info.csv"
# and then assigns a seawater T and d18O value to each site based on the gridded model of Breitkreuz et al. (2018)

# INPUT: SK Table S-1.csv, SK Table S-3 part-1.csv, D18O_Breitkreuz_et_al_2018.nc
# OUTPUT: SK Figure S3.png, SK Table S-3 part-2.csv

# >>>>>>>>>

# Import libraries
from scipy.interpolate import interp1d
import pandas as pd
from csv import DictReader
import netCDF4 as nc
from pylab import *
import sys

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

d18Osw_model_sigma = 0.1
Tsw_model_sigma = 1.0

# This file is not included in the repository due to its size
# It can be downloaded from https://doi.pangaea.de/10.1594/PANGAEA.889922
ds = nc.Dataset(sys.path[0]+"/D18O_Breitkreuz_et_al_2018.nc")

lats = ds['lat_1deg_center'][:,0]
lons = ds['lon_1deg_center'][0,:]
zc = ds['depth_center'][:]
ze = ds['depth_edge'][:]
depths = -concatenate((ze[:1], zc[:]))

d18o = ds['D18O_1deg'][:,:,:,:] # month, depth, lat, lon
d18o = concatenate((d18o[:,:1,:,:], d18o[:,:,:,:]), axis = 1)

T = ds['THETA_1deg'][:,:,:,:] # month, depth, lat, lon
T = concatenate((T[:,:1,:,:], T[:,:,:,:]), axis = 1)

d18o = d18o.filled(fill_value = nan)
d18o_mask = (isnan(d18o)).astype(int)

T = T.filled(fill_value = nan)
T_mask = (isnan(T)).astype(int)

glon, glat = meshgrid(lons, lats)
gx = cos(glon * pi / 180) * cos(glat * pi / 180)
gy = sin(glon * pi / 180) * cos(glat * pi / 180)
gz = sin(glat * pi / 180)

with open(sys.path[0] + "/SK Table S-1.csv") as f:
	Samples = [{k: r[k] for k in r} for r in DictReader(f)]
	print(Samples[0].keys())
for r in Samples:
	for k in r:
		if k not in ["SampleName", "Type", "Species", "T_measured", "d13CDIC", "d18Osw_measured"]:
			r[k] = float(r[k])

print('Extracting seawater d18O for corals...')

df_model = pd.DataFrame(columns=['SampleName', 'T_modeled', 'T_modeled_err', 'd18Osw_modeled', 'd18Osw_modeled_err'])

for s in Samples:
	Sample, lon, lat, depth = s['SampleName'], s['Long'], s['Lat'], s['Depth']

	print(f'\tProcessing {Sample}')

	x = cos(lon * pi / 180) * cos(lat * pi / 180)
	y = sin(lon * pi / 180) * cos(lat * pi / 180)
	z = sin(lat * pi / 180)
	sqdistance = (gx-x)**2 + (gy-y)**2 + (gz-z)**2

	i = [i for i, _ in enumerate(depths) if _ >= depth][0]

	sqdistance += d18o_mask[0,i,:,:] * 10
	min_index = np.unravel_index(np.argmin(sqdistance, axis=None), sqdistance.shape)
	j, k = [int(_) for _ in min_index]

	fig = figure(figsize = (8,4))
	ax1, ax2 = subplot(121), subplot(122)
	subplots_adjust(.15, .15, .95, .9, .25)
	
	X, Y, Tloc, M = depths[:], d18o[:,:,j,k], T[:,:,j,k], d18o_mask[0,:,j,k]
	X, Y, Tloc = X[M<1], Y[:,M<1], Tloc[:,M<1]

	maxdepth = X[-1]

	d18values = []
	Tvalues = []
	for y in Y:
		sca(ax1)
		plot(y, -X, 'b-', alpha = .1, label = 'database $\delta^{18}$O$_{sw}$')
		f = interp1d(X,y)
		d18values += [f(depth)]

	for t in Tloc:
		sca(ax2)
		plot(t, -X, 'r-', alpha = .1, label = 'database T')
		f = interp1d(X,t)
		Tvalues += [f(depth)]

	kw = dict(elinewidth = 1.5, alpha = 1, capsize = 5, marker = '+', ls = 'None', capthick = 1.5)

	d18values = array(d18values)
	d18, sd18 = d18values.mean(), (d18Osw_model_sigma**2 + d18values.std(ddof = 1)**2)**.5
	sca(ax1)
	errorbar(d18, -depth, None, 1.96*sd18, ecolor='b', c='b',
				label="estimated $\delta^{18}$O$_{sw}$", **kw)
	xlabel("$\delta^{18}$O$_{seawater}$ (‰ VSMOW)")
	ylabel('depth (m)')
	plt.legend()

	handles, labels = ax1.get_legend_handles_labels()
	by_label = dict(zip(labels, handles))
	ax1.legend(by_label.values(), by_label.keys())

	text(.5, .97, f'{Sample}', va = 'top', ha = 'center', transform = fig.transFigure)

	Tvalues = array(Tvalues)
	t, st = Tvalues.mean(), (Tsw_model_sigma**2 + Tvalues.std(ddof = 1)**2)**.5
	sca(ax2)
	errorbar(t, -depth, None, 1.96*st, ecolor='r',
				c='r', label="estimated T", **kw)

	handles, labels = ax2.get_legend_handles_labels()
	by_label = dict(zip(labels, handles))
	ax2.legend(by_label.values(), by_label.keys())

	xlabel('Temperature ($^{\circ}$C)')

	new_data = {'SampleName': Sample, 'T_modeled': t, 'T_modeled_err': st, 'd18Osw_modeled': d18, 'd18Osw_modeled_err': sd18}
	df_model = df_model._append(new_data, ignore_index=True)

	# Save figures
	# savefig(sys.path[0] + "/isoForam models/" + f'{Sample} model.png', dpi=150)

	close(fig)


df_measurements = pd.read_csv(sys.path[0] + "/SK Table S-3 part-1.csv")
df_Info = pd.read_csv(sys.path[0] + "/SK Table S-1.csv")[
    ["SampleName", "Species", "Type", "T_measured", "d18Osw_measured"]]

df = df_measurements.merge(
	df_Info, on='SampleName').merge(df_model, on='SampleName')

# Save data to CSV
df.to_csv(sys.path[0] + "/SK Table S-3 part-2.csv", index=False)


# Create Figure S3

# Assign colors and markers
categories = df["SampleName"].unique()
markers = dict(zip(categories, [
    "o", "s", "D", "^", "v", "X", "P", "*", "o", "s", "D", "^", "v", "X", "P", "*", "o", "s"]))
colors = dict(zip(categories, plt.cm.tab20(
    np.linspace(0, 1, len(categories)))))


# Subplot A: Difference between measured and modelled d18Osw
fig, (ax1, ax2) = plt.subplots(1, 2)

for cat in categories:
    data = df[df["SampleName"] == cat]
    ax1.scatter(data['d18Osw_measured'], data['d18Osw_modeled'],
                marker=markers[cat], fc=colors[cat], label=cat)

# Calculate and annotate difference between measured and modeled d18Osw
D18Osw = df['d18Osw_measured']-df['d18Osw_modeled']
for i, difference in enumerate(D18Osw):
    ax1.annotate(f"{difference:.1f}", (df['d18Osw_measured'][i]-0.05, df['d18Osw_modeled'][i]),
                 ha='right', va='center')

# 1:1 line
ax1.plot([0, 1.75], [0, 1.75], c = "k", ls="dashed", zorder = -1)
xmin, xmax = ax1.get_xlim()
ymin, ymax = ax1.get_ylim()
angle = np.arctan((xmax-xmin)/(ymax-ymin)) * 180 / np.pi
ax1.text(1.5, 1.45, "1:1", ha="center", va="center", rotation=angle)

ax1.text(0.02, 0.98, "a", size=14, ha="left", va="top",
         transform=ax1.transAxes, fontweight="bold")

ax1.set_xlabel('Measured $\delta^{18}$O$_{sw}$ (‰, VSMOW)')
ax1.set_ylabel('Modelled $\delta^{18}$O$_{sw}$ (‰, VSMOW)')


# Subplot B: Difference between measured and modeled temperature

for cat in categories:
    data = df[df["SampleName"] == cat]
    ax2.scatter(data['T_measured'], data['T_modeled'],
                marker=markers[cat], fc=colors[cat], label=cat)

# Calculate and annotate difference between measured and modeled temperature
DT = df['T_measured']-df['T_modeled']
for i, difference in enumerate(DT):
    ax2.annotate(f"{difference:.1f}", (df['T_measured'][i]-1, df['T_modeled'][i]),
                 ha='right', va='center')

# 1:1 line
ax2.plot([0, 30], [0, 30], c="k", ls="dashed", zorder=-1)
xmin, xmax = ax1.get_xlim()
ymin, ymax = ax1.get_ylim()
angle = np.arctan((xmax-xmin)/(ymax-ymin)) * 180 / np.pi
ax2.text(20, 19, "1:1", ha="center", va="center", rotation=angle)

plt.text(0.02, 0.98, "b", size=14, ha="left", va="top",
         transform=ax2.transAxes, fontweight="bold")

ax2.set_xlabel('Measured temperature (°C)')
ax2.set_ylabel('Modelled temperature (°C)')

ax2.legend(loc='upper right', bbox_to_anchor=(1.35, 1))

plt.savefig(sys.path[0] + "/SK Figure S3.png")
plt.close("all")