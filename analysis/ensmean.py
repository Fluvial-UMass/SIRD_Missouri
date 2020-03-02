import pandas as pd
import numpy as np
import datetime
import argparse
import os

parser = \
    argparse.ArgumentParser(description="take ensemble mean")
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
        outpath = os.path.join(args.outdir, args.name,
                               "dischargeAssim_mean.csv")
    elif args.case == "open":
        path = os.path.join(args.outdir, args.name,
                            "{0:02d}", "discharge.csv")
        outpath = os.path.join(args.outdir, args.name,
                               "discharge_mean.csv")
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
        cols = df_sel.columns
        inds = df_sel[sdate:edate].index
        dim0 = values.shape[0]
        dim1 = values.shape[1]
        all.append(values)
    ensmean = np.array(all).reshape(-1, dim0, dim1).mean(axis=0)
    print("output shape:", ensmean.shape)
    enDf = pd.DataFrame(ensmean, columns=cols, index=inds)
    enDf.to_csv(outpath)


if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
