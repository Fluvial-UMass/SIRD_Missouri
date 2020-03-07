import pandas as pd
import numpy as np
import datetime
import argparse
import xarray as xr
import os

parser = \
    argparse.ArgumentParser(description="summarize output as netCDF")
parser.add_argument("-n", "--name", type=str, help="experiment name",
                    required=True)
parser.add_argument("-c", "--case", type=str, help="case [open/assim]",
                    required=True)
parser.add_argument("-s", "--start", type=str,
                    required=True,
                    help="start date in YYYYMMDD")
parser.add_argument("-e", "--end", type=str,
                    required=True,
                    help="end date in YYYYMMDD")
parser.add_argument("-t", "--total", type=int,
                    required=True,
                    help="total ensembles you have.")
parser.add_argument("-o", "--outdir", type=str,
                    required=False,
                    default="../../out/",
                    help="output directory")
parser.add_argument("--log", default=False,
                    required=False,
                    action="store_true",
                    help="take log before take mean or percentiles")
parser.add_argument("--otype", type=str,
                    required=False,
                    default="nc",
                    help="output format; csv/nc")
parser.add_argument("-v", "--verbose", default=False,
                    required=False,
                    action="store_true",
                    help="show every dataframe")


def main(args):
    sdate = datetime.datetime.strptime(args.start, "%Y%m%d")
    edate = datetime.datetime.strptime(args.end, "%Y%m%d")
    if args.case == "assim":
        path = os.path.join(args.outdir, args.name,
                            "{0:02d}", "dischargeAssim.csv")
        if args.otype == "nc":
            outpath = os.path.join(args.outdir, args.name,
                                   "dischargeAssim.nc")
        elif args.otype == "csv":
            outpath = os.path.join(args.outdir, args.name,
                                   "dischargeAssim_mean.csv")
        else:
            raise(KeyError(args.otype))
    elif args.case == "open":
        path = os.path.join(args.outdir, args.name,
                            "{0:02d}", "discharge.csv")
        if args.otype == "nc":
            outpath = os.path.join(args.outdir, args.name,
                                   "discharge.nc")
        elif args.otype == "csv":
            outpath = os.path.join(args.outdir, args.name,
                                   "discharge_mean.csv")
        else:
            raise(KeyError(args.otype))
    else:
        KeyError("case {0} is not a valid argument.".format(args.name))
    all = []
    for idx, ens in enumerate(range(args.total)):
        print(ens)
        df = pd.read_csv(path.format(ens),
                         header=0,
                         parse_dates=[0]).set_index("Date")
        df_sel = df[sdate:edate]
        df_sel = df_sel.loc[~df_sel.index.duplicated(keep="last")]
        if args.verbose:
            print(df_sel)
        values = df_sel.values
        cols = df_sel.columns.astype(int)
        inds = df_sel[sdate:edate].index
        dim0 = values.shape[0]
        dim1 = values.shape[1]
        values = values.reshape(1, dim0, dim1)
        all.append(values)
    allarray = np.vstack(all)
    ensMean = (take_mean(allarray, args.log)).reshape(1, dim0, dim1)
    if args.otype == "csv":
        outdf = pd.DataFrame(ensMean[0], columns=cols, index=inds)
        outdf.to_csv(outpath, header=True, index=True)
    else:
        allDArray = xr.DataArray(allarray,
                                 dims=["ens", "time", "reach"],
                                 coords=[np.arange(0, args.total),
                                         inds,
                                         cols])
        ens75th = (take_quantile(allarray,
                                 0.75,
                                 args.log)).reshape(1, dim0, dim1)
        ens25th = (take_quantile(allarray,
                                 0.25,
                                 args.log)).reshape(1, dim0, dim1)
        ensMax = allarray.max(axis=0).reshape(1, dim0, dim1)
        ensMin = allarray.min(axis=0).reshape(1, dim0, dim1)
        reprarray = np.vstack([ensMean,
                               ens75th,
                               ens25th,
                               ensMax,
                               ensMin])
        reprDArray = xr.DataArray(reprarray,
                                  dims=["kind", "time", "reach"],
                                  coords=[["mean", "p75", "p25", "max", "min"],
                                          inds,
                                          cols])
        outDset = xr.Dataset({"all": allDArray, "rpr": reprDArray})
        outDset.to_netcdf(outpath)


def take_mean(array3d, log):
    if log:
        array3d = np.where(array3d == 0, 1e-8, array3d)
        logMean = (np.log(array3d)).mean(axis=0)
        ensMean = np.exp(logMean)
    else:
        ensMean = array3d.mean(axis=0)
    return ensMean


def take_quantile(array3d, qvalue, log):
    if log:
        array3d = np.where(array3d == 0, 1e-8, array3d)
        logQuant = np.quantile(np.log(array3d), qvalue, axis=0)
        ensQuant = np.exp(logQuant)
    else:
        ensQuant = np.quantile(array3d, qvalue, axis=0)
    return ensQuant


if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
