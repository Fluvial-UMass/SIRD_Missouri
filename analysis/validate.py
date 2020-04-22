import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import datetime
import calendar

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


def main(ncpath_gauge, ncpath_sim, sDate, eDate, case, stype="mean", ens=None):
    print(stype, ens)
    gauge = read_dataset(ncpath_gauge, "gauge", "discharge")
    gaugeReaches = gauge["reach"].values.tolist()
    if ens is None:
        simulation = read_dataset(ncpath_sim, "rpr", stype)
    elif isinstance(ens, int):
        simulation = read_dataset_ens(ncpath_sim, "all", ens)
    elif ens == "baseline":
        simulation = read_dataset_baseline(ncpath_sim, "discharge", "discharge")
    elif ens == "insert":
        simulation = read_dataset_ens(ncpath_sim, "all", 0)
    else:
        KeyError("invalid data type: {}".format(isinstance(ens)))
    sList = []
    for reach in gaugeReaches:
        sim = simulation.sel(reach=reach, time=slice(sDate, eDate))
        gag = gauge.sel(reach=reach, time=slice(sDate, eDate))
        df = pd.DataFrame([sim.values, gag.values]).T
        df = df.dropna(how="any")
        sts = stats(df)
        sts.append(reach)
        sList.append(sts)
    sArray = np.array(sList)
    df = pd.DataFrame(sArray,columns=["r","rmse","nrmse","rrmse","rbias","nse","kge","hrrid"])
    outname = ncpath_sim.split(".")[-2].split("/")[-1]
    df.to_csv("./stats/{0}_{1}_stats.csv".format(outname, case), index=False)


def stats(df):

    obs = df.iloc[:,1].values
    sim = df.iloc[:,0].values

    # valid_index = filterout_outliers(sim)
    # obs = obs[valid_index]
    # sim = sim[valid_index]

    r = np.corrcoef(obs,sim)[0][1]
    rmse = RMSE(obs,sim)
    nrmse = nRMSE(obs,sim)
    rrmse = rRMSE(obs,sim)
    rbias = rBIAS(obs,sim)
    nse = NSE(obs,sim)
    KGE = kge(obs,sim)

    return [r,rmse,nrmse,rrmse,rbias,nse,KGE]


def filterout_outliers(array, thsld=50):
    index = np.where(array<array.mean()*thsld)
    return index


def RMSE(obs,sim):
    rmse = (np.nanmean(((obs-sim)**2)))**0.5
    return rmse


def nRMSE(obs,sim):
    rmse = (np.nanmean(((obs-sim)**2)))**0.5
    nRmse = rmse/np.nanmean(obs)
    return nRmse


def rRMSE(obs,sim):
    obs_f = obs[obs!=0]
    sim_f = sim[obs!=0]
    rRMSE = (np.nanmean(((obs_f-sim_f)/obs_f)**2))**0.5
    return rRMSE


def rBIAS(obs, sim):
    rbias = (sim - obs).sum() / obs.sum()
    return rbias


def NSE(obs,sim):
    right = np.nansum((obs-sim)**2)/np.nansum((obs-np.nanmean(obs))**2)
    nse = 1 - right
    return nse


def kge(obs, sim):
    pcc = np.corrcoef(sim, obs)[0, 1]
    alpha = (np.nanstd(sim)/np.nanmean(sim))/(np.nanstd(obs)/np.nanmean(obs))
    beta = np.nanmean(sim)/np.nanmean(obs)
    kgesq = (pcc-1)**2 + (alpha-1)**2 + (beta-1)**2
    kge = 1 - kgesq**0.5
    if np.isinf(kge):
        kge = np.nan
    return kge


if __name__ == "__main__":
    ncpath_gauge = "../../data/obs/gauge.nc"
    ncpath_sims = ["../../out/old/MSR_cal10_log/dischargeAssim.nc",
                "../../out/MSR_uncal_insert/dischargeInsert_summarized.nc"]
    case = ["cal10_logdelta",
            "uncal_insert"]
    stype = ["mean",
            "mean"]
    enum = [None,
            "insert"]
    sDate = datetime.datetime(2003,1,1)
    eDate = datetime.datetime(2011,1,1)
    for idx, ncpath_sim in enumerate(ncpath_sims):
        main(ncpath_gauge, ncpath_sim, sDate, eDate, case[idx],
             stype=stype[idx], ens=enum[idx])
