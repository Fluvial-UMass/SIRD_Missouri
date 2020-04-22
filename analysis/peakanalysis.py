import numpy as np
import pandas as pd
import xarray as xr
import datetime

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


def read_dataset_baseline(ncpath, dsetname, key):
    data = xr.open_dataset(ncpath)
    darray = data[dsetname]
    print(darray)
    return darray


def main(exp, ref, gauge, reaches, sdate, edate, case):
    statsall = []
    for reach in reaches:
        expvec = exp.sel(reach=reach, time=slice(sdate, edate)).values
        refvec = ref.sel(reach=reach, time=slice(sdate, edate)).values
        gaugevec = gauge.sel(reach=reach, time=slice(sdate, edate)).values
        stats = analysis(case, expvec, refvec, gaugevec)
        statsall.append(stats)
    statsall_df = pd.DataFrame(statsall)
    statsall_df.columns = ["nrmse_ref", "variation_ref", "bias_ref", "nrmse_exp", "variation_exp", "bias_exp"]
    # print(statsall_df)
    statsall_df.to_csv("./stats/peakstats_{0}.csv".format(case),
                       index=False, header=True)
    # print(statsall_df.describe())
    return statsall_df


def analysis(case, expvec, refvec, gaugevec, threshold=0.75):
    # peak flow in ref.
    thsld_val = np.quantile(refvec, threshold)
    if case == "above":
        peak_index = np.where(refvec >= thsld_val)
    elif case == "below":
        peak_index = np.where(refvec < thsld_val)
    else:
        raise KeyError(case)
    # extract for other case
    ref_peak = refvec[peak_index]
    exp_peak = expvec[peak_index]
    gauge_peak = gaugevec[peak_index]
    nrmse_ref = calc_nrmse(ref_peak, gauge_peak)
    nrmse_exp = calc_nrmse(exp_peak, gauge_peak)
    variation_ref = calc_variation(ref_peak, gauge_peak)
    variation_exp = calc_variation(exp_peak, gauge_peak)
    bias_ref = calc_volumetricbias(ref_peak, gauge_peak)
    bias_exp = calc_volumetricbias(exp_peak, gauge_peak)
    return [nrmse_ref, variation_ref, bias_ref, nrmse_exp, variation_exp, bias_exp]


def calc_nrmse(sim, obs):
    mse = np.nanmean(((sim-obs)**2))
    rmse = np.sqrt(mse)
    nrmse = rmse/np.nanmean(obs)
    return nrmse


def calc_variation(sim, obs):
    return np.nanstd(sim)/np.nanstd(obs) if np.nanstd(obs) != 0 else np.nan


def calc_volumetricbias(sim, obs):
    return np.nanmean(sim)/np.nanmean(obs)


if __name__ == "__main__":
    ncpath_assim = "../../out/MSR_cal10_logdelta/dischargeAssim.nc"
    ncpath_baseline = "../../out/baseline_cal10/discharge.nc"
    ncpath_gauge = "../../data/obs/gauge.nc"
    sdate = datetime.datetime(2003,1,1)
    edate = datetime.datetime(2010,12,31)
    exp = read_dataset(ncpath_assim,
                       "rpr", "mean")
    ref = read_dataset_baseline(ncpath_baseline,
                                "discharge", "discharge")
    gauge = read_dataset(ncpath_gauge, "gauge", "discharge")

    reaches = gauge["reach"].values.tolist()
    statsdf = main(exp, ref, gauge, reaches, sdate, edate, "above")
    statsdf = main(exp, ref, gauge, reaches, sdate, edate, "below")
