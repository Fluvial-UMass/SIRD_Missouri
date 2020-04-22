# %%
import numpy as np
import xarray as xr
import h5py
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import datetime
import os
sns.set_style("whitegrid")


def read_dataset(ncpath, dsetname, key):
    """
    Args:
        ncpath (str)
        dsetname (str)
        key (str): "kind" coordinate name.
    Notes:
        netCDF dataset should have "kind" coordinate.
    """
    data = xr.open_dataset(ncpath)
    darray = data[dsetname].sel(kind=key)
    return darray


def read_dataset_ens(ncpath, dsetname, key):
    """
    Args:
        ncpath (str)
        dsetname (str)
        key (str): "ens" coordinate name.
    Notes:
        netCDF dataset should have "kind" coordinate.
    """
    data = xr.open_dataset(ncpath)
    darray = data[dsetname].sel(ens=key)
    return darray


def read_dataset_baseline(ncpath, dsetname, key):
    data = xr.open_dataset(ncpath)
    darray = data[dsetname]
    print(darray)
    return darray


def draw_hydrograph(ax, darrays, labels, colors, linestyles,
                     linewidths, timelabel="time"):
    for i in range(len(darrays)):
        darray = darrays[i]
        if np.isnan(darray.values).all():
            continue
        plt.plot(darray[timelabel],
                 darray.values,
                 label=labels[i],
                 color=colors[i],
                 linestyle=linestyles[i],
                 linewidth=linewidths[i])
    plt.xlabel("Date")
    plt.ylabel("Discharge")
    # plt.title(darray.reach.values, weight="bold")
    return ax


def add_shades_on_plot(ax, top_darrays, bottom_darrays, colors,
                       linestyles, linewidths, timelabel="time"):
    for i in range(len(top_darrays)):
        top = top_darrays[i]
        bottom = bottom_darrays[i]
        if np.isnan(top.values).all():
            continue
        plt.plot(top[timelabel],
                 top.values,
                 color=colors[i],
                 linestyle=linestyles[i],
                 linewidth=linewidths[i])
        plt.plot(bottom[timelabel],
                 bottom.values,
                 color=colors[i],
                 linestyle=linestyles[i],
                 linewidth=linewidths[i])
        plt.fill_between(top[timelabel],
                         bottom.values,
                         top.values,
                         color=colors[i],
                         alpha=0.25)
    return ax


def add_obs_on_plot(ax, darrays_point, darrays_err, labels,
                    colors, styles, timelabel="time"):
    for i in range(len(darrays_point)):
        if np.isnan(darrays_point[i].values).all():
            continue
        points = darrays_point[i].values
        errs = darrays_err[i].values
        plt.errorbar(darrays_point[i][timelabel].values,
                     points,
                     yerr=errs,
                     color=colors[i],
                     marker=styles[i],
                     label=labels[i])
    return ax


def add_obstime_on_plot_core(ax, patch, obsdarray, ymax,
                             timelabel="time"):
    obsreaches = obsdarray["reach"].values.tolist()
    obsreaches_patch = list(set(patch) & set(obsreaches))
    if len(obsreaches_patch) == 0:
        return ax
    else:
        obsdarray_patch = obsdarray.sel(reach=obsreaches_patch).values
        obsflags = np.where(np.isnan(obsdarray_patch), 0, 1).sum(axis=1)
        dates = obsdarray[timelabel].values
        dates_plot = dates[obsflags > 0]
        obsflags_plot = obsflags[obsflags > 0]
        plt.axvline(dates_plot[0], ymin=0, ymax=ymax, color="#F39C12",
                    linewidth=2, alpha=0.75, label="overpass")
        for idx, date_plot in enumerate(dates_plot[1::]):
            plt.axvline(date_plot, ymin=0, ymax=ymax, color="#F39C12",
                        linewidth=2, alpha=0.75)
        return ax


def add_obstime_on_plot(ax, reach, patches, upstreams, obsdarray, ymax,
                        reachstart=1):
    patch = np.array(patches[reach-reachstart]) + 1  # from 0-start to 1-start
    upstream = np.array(upstreams[reach-reachstart]) + 1 # from 0-start to 1-start
    combined_patch = patch.tolist()
    combined_patch.extend(upstream.tolist())
    combined_patch = list(set(combined_patch))
    ax = add_obstime_on_plot_core(ax, combined_patch, obsdarray, ymax)
    return ax


def tester():
    # define path
    ncpath_baseline = "../../out/baseline_cal10/discharge.nc"
    ncpath_assim = "../../out/old/MSR_cal10_log/dischargeAssim.nc"
    # ncpath_assim = "../../out/MSR_cal10_insert/dischargeInsert.nc"
    ncpath_gauge = "../../data/obs/gauge.nc"
    ncpath_bam = "../../data/obs/bam_wgauge_1984_2010.nc"
    path_upreaches = "../../data/hrr/upReaches.hdf5"
    path_localpatch = "../../data/hrr/localPatch_MERIT.hdf5"
    outdir = "./paper"

    # define date
    sdate = datetime.datetime(2004, 1, 1)
    edate = datetime.datetime(2005, 1, 1)

    # read data
    baseline = read_dataset_baseline(ncpath_baseline, "discharge", "discharge")
    assim = read_dataset(ncpath_assim, "rpr", "mean")
    # assim = read_dataset_baseline(ncpath_assim, "discharge", "dicharge")
    assim_top = read_dataset(ncpath_assim, "rpr", "p75")
    assim_bottom = read_dataset(ncpath_assim, "rpr", "p25")
    gauge = read_dataset(ncpath_gauge, "gauge", "discharge")
    bam_mean = read_dataset(ncpath_bam, "bam", "mean")
    bam_std = read_dataset(ncpath_bam, "bam", "std")
    upstreams = h5py.File(path_upreaches, mode="r")["upstream_reaches"][:]
    patches = h5py.File(path_localpatch, mode="r")["network"][:]

    # create dummy data array
    dummy = baseline.copy()
    dummy.values = np.ones_like(baseline.values)*np.nan

    # mkdir
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # decide which reaches we draw
    valreaches = gauge["reach"].values.tolist()
    # bamreaches = bam_mean["reach"].values.tolist()
    # valreaches.extend(bamreaches)
    # valreaches = list(set(valreaches))

    # for valreach in valreaches:
    for valreach in [28691, 28835]:
        print(valreach)
        # select data and set dummy where no data
        b_darray = baseline.sel(reach=valreach, time=slice(sdate, edate))
        a_darray = assim.sel(reach=valreach, time=slice(sdate, edate))
        a_top_darray = assim_top.sel(reach=valreach, time=slice(sdate, edate))
        a_bottom_darray = assim_bottom.sel(reach=valreach, time=slice(sdate, edate))
        try:
            g_darray = gauge.sel(reach=valreach, time=slice(sdate, edate))
        except(KeyError):
            g_darray = dummy.sel(reach=valreach, time=slice(sdate, edate))
        try:
            bm_darray = bam_mean.sel(reach=valreach, time=slice(sdate, edate))
            bs_darray = bam_std.sel(reach=valreach, time=slice(sdate, edate))
        except(KeyError):
            bm_darray = dummy.sel(reach=valreach, time=slice(sdate, edate))
            bs_darray = dummy.sel(reach=valreach, time=slice(sdate, edate))
        # main lines
        darrays = [b_darray, a_darray, g_darray]
        labels = ["baseline", "assim", "gauge"]
        colors = ["#5499C7", "#E74C3C", "k"]
        linestyles = ["-", "-", "--"]
        linewidths = [2, 2, 2]
        # shades (uncertainty)
        top_darrays = [a_top_darray]
        bottom_darrays = [a_bottom_darray]
        colors_shade = ["#E74C3C"]
        linestyles_shade = ["-"]
        linewidths_shade = [0.5]
        # draw
        ax = plt.figure(figsize=(10, 5))
        ymax_ = max(a_darray.max().values,
                    b_darray.max().values,
                    g_darray.max().values)
        ymax = ymax_ + ymax_/10
        # add shades
        x1 = datetime.datetime(2004,3,1)
        x2 = datetime.datetime(2004,4,1)
        y = [0, ymax*2]
        plt.fill_betweenx(y, x1, x2, color="gray", alpha=0.75, label="peak-cut period")
        # observaiton time
        obsdarray = bam_mean.sel(time=slice(sdate, edate))
        ax = add_obstime_on_plot(ax, valreach, patches, upstreams,
                                 obsdarray, ymax)
        ax = draw_hydrograph(ax,darrays, labels, colors, linestyles, linewidths)
        ax = add_shades_on_plot(ax, top_darrays, bottom_darrays,
                                colors_shade, linestyles_shade, linewidths_shade)
        ax = add_obs_on_plot(ax, [bm_darray], [bs_darray], ["bam"],
                             ["#1F618D"], ["o"])
        plt.legend(loc="upper right")
        plt.ylim(0, ymax)
        plt.xlabel("Date", weight="bold")
        plt.ylabel("Discharge [cms]", weight="bold")
        plt.savefig(os.path.join(outdir, "{0}.png".format(valreach)),
                    dpi=300)
        plt.close()



if __name__ == "__main__":
    tester()
