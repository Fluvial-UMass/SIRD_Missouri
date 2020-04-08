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
parser.add_argument("--itype", type=str,
                    required=False,
                    default="nc",
                    help="input format; csv/nc")
parser.add_argument("--otype", type=str,
                    required=False,
                    default="nc",
                    help="output format; csv/nc")
parser.add_argument("--skiprows", type=int,
                    required=False,
                    default=0,
                    help="skip rows from the top")
parser.add_argument("--drop_duplicates", default=False,
                    required=False,
                    action="store_true",
                    help="[beta] remove duplicated indices")
parser.add_argument("-v", "--verbose", default=False,
                    required=False,
                    action="store_true",
                    help="show every dataframe")


def main(args):
    sdate = datetime.datetime.strptime(args.start, "%Y%m%d")
    edate = datetime.datetime.strptime(args.end, "%Y%m%d")
    path, outpath = define_paths(args)
    allarray, \
    cols, inds, dim_time, dim_reach = read_all(args, path, sdate, edate)
    print(allarray.shape)
    ensMean = (take_mean(allarray, args.log)).reshape(1, dim_time, dim_reach)
    if args.otype == "csv":
        outdf = pd.DataFrame(ensMean[0], columns=cols, index=inds)
        outdf.to_csv(outpath, header=True, index=True)
    elif args.otype == "nc":
        outDset = \
            summarize(args, ensMean, allarray, cols, inds, dim_time, dim_reach)
        outDset.to_netcdf(outpath)


def define_paths(args):
    if not args.itype in ["nc", "csv"]:
        raise(KeyError("invalid itype:", args.itype))
    if not args.otype in ["nc", "csv"]:
        raise(KeyError("invalid otype:", args.otype))
    if args.case == "assim":
        path = os.path.join(args.outdir, args.name,
                            "{0:02d}", "dischargeAssim.{0}".format(args.itype))
        outpath = os.path.join(args.outdir, args.name,
                               "dischargeAssim.{0}".format(args.itype))
    elif args.case == "open":
        path = os.path.join(args.outdir, args.name,
                            "{0:02d}", "discharge.{0}".format(args.itype))
        outpath = os.path.join(args.outdir, args.name,
                               "discharge.{0}".format(args.otype))
    else:
        KeyError("case {0} is not a valid argument.".format(args.name))
    return path, outpath


def read_all(args, path, sdate, edate, drop_duplicates=False):
    all = []
    for idx, ens in enumerate(range(args.total)):
        srcpath = path.format(ens)
        print(srcpath)
        values, cols, inds = read_input(srcpath, sdate, edate,
                                        args.drop_duplicates)
        dim0 = values.shape[0]
        dim1 = values.shape[1]
        values = values.reshape(1, dim0, dim1)
        all.append(values)
    allarray = np.vstack(all)
    return allarray, cols, inds, dim0, dim1


def read_input_csv(path, sdate, edate, drop_duplicates=False):
    df = pd.read_csv(path,
                     header=0,
                     parse_dates=[0]).set_index("Date")
    df = df.iloc[args.skiprows::]
    df_sel = df[sdate:edate]
    if drop_duplicates:
        df_sel = df_sel.loc[~df_sel.index.duplicated(keep="last")]
    if args.verbose:
        print(df_sel)
    values = df_sel.values
    cols = df_sel.columns.astype(int)
    inds = df_sel[sdate:edate].index
    return values, cols, inds


def read_input_nc(path, sdate, edate, drop_duplicates=False,
                  dsetname="discharge", vkey="discharge",
                  kkey="reach"):
    darr = xr.open_dataset(path)["discharge"]
    if args.verbose:
        print(darr)
    darr_sel = darr.sel(time=slice(sdate, edate))
    values = darr_sel.values
    cols = darr_sel[kkey].values
    inds = darr_sel["time"].values
    if drop_duplicates:
        print(
            RuntimeWarning(
                "drop_duplicates is not fully examined yet."+
                " Use this for your own risk.")
            )
        _, newidx = np.unique(inds[::-1], return_index=True)
        values = values[newidx]
        inds = inds[newidx]
    return values, cols, inds


def read_input(path, sdate, edate, drop_duplicates=False):
    if args.itype == "csv":
        values, cols, inds = read_input_csv(path, sdate, edate,
                                            drop_duplicates)
    elif args.itype == "nc":
        values, cols, inds = read_input_nc(path, sdate, edate,
                                           drop_duplicates)
    return values, cols, inds


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


def summarize(args, ensMean, allarray, cols, inds, dim0, dim1):
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
    return outDset


if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
