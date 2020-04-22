# %%%
import geopandas as gpd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import mpl_toolkits.axes_grid1
import pandas as pd
sns.set(font_scale=2.0)
sns.set_style("whitegrid")
sns.set_palette("deep")

# %%
statname_exp = "nse_assimilated"
statname_ref = "nse_baseline"
gdf = gpd.read_file("./MissouriShapes/Missouri_rivers_MERITv0.0.shp").set_index("COMID")
distdf = pd.read_csv("./bamdist.csv", index_col=[0])

# %%
stats_exp = pd.read_csv("../improvementplot/stats/stats/dischargeAssim_cal10_log_stats.csv").loc[:, ["hrrid", statname_exp]].sort_values("hrrid")
stats_ref = pd.read_csv("../improvementplot/stats/stats/discharge_cal10_baseline_stats.csv").loc[:, ["hrrid", statname_ref]].sort_values("hrrid")
stats_exp["hrrid"] = stats_exp["hrrid"].astype(int)
stats_ref["hrrid"] = stats_ref["hrrid"].astype(int)
stats_exp.columns = ["hrrid", "nse"]
stats_ref.columns = ["hrrid", "nse"]
stats_imp = stats_exp - stats_ref
stats_imp["hrrid"] = stats_exp["hrrid"]

#%%
locinfo = pd.read_csv("./gaugeInfo.csv")
locinfo["COMID"] = locinfo["COMID"].astype(int)
locinfo = locinfo.set_index("COMID")
ref = pd.read_csv("./comid_hrrid.csv").astype(int)
ref = ref.set_index("HRRID")
stats = stats_imp.assign(COMID=stats_imp.hrrid.apply(lambda x: ref.loc[x].values[0]))
dwstreaches = distdf[distdf.upnums > 0].index.tolist()
dwstcomids = [ref.loc[r].values[0] for r in dwstreaches]

#%%
stats = locinfo.join(stats.set_index("COMID")).dropna()
stats = stats.assign(lat=stats.hrrid.apply(lambda x: locinfo.loc[ref.loc[x]].lat.values[0]))
stats = stats.assign(lon=stats.hrrid.apply(lambda x: locinfo.loc[ref.loc[x]].lon.values[0]))

# %%
# my_cmap = sns.cubehelix_palette(24, start=.5, rot=-.75, as_cmap=True)
my_cmap = colors.ListedColormap(["#C0392B", "#D4E6F1", "#5499C7", "#1A5276"])
bounds = np.array([-1, -0.1, 0.1, 1, 5])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=4)
ax = gdf.plot(figsize=(20,15), color="grey")
gdf.loc[dwstcomids].plot(color="k", ax=ax, linewidth=2.0)
plt.xlabel("Lon [deg.]", weight="bold")
plt.ylabel("Lat [deg.]", weight="bold")
plt.scatter(stats.lon, stats.lat, c=stats.nse, s=200, zorder=2, cmap=my_cmap, norm=norm, linewidths=1, edgecolors="k")
# plt.scatter(stats.lon, stats.lat, c=stats.nse, s=200, zorder=2, norm=norm)
divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
cax = divider.append_axes('right', '5%', pad='3%')
plt.colorbar(cax=cax, extend="max")
plt.savefig("nse_mean_map.png", dpi=300, bbox_inches="tight", pad_inches=0.1)

# %%


# %%
