import pyHRR as pyHRR
import pandas as pd
import configparser
from tqdm import tqdm
import datetime
import argparse
import os

def run(config, sdate, edate):
    runner = pyHRR.HRR(config)
    initFile = configparser.ConfigParser()
    initFile.read(config)
    runoffdir = initFile.get("IO", "runoffDir")
    outdir = os.path.join(initFile.get("IO", "outDir"),
                          initFile.get("experiment", "expName"))
    outfile = initFile.get("IO", "oFile")
    outpath = os.path.join(outdir, outfile)
    dates = pd.date_range(sdate, edate, freq="D")

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    # initial run
    df, ndate = \
        runner.main_day(dates[0], flag="initial", restart="restart.bin",
                        runoffDir=runoffdir, outDir=outdir)
    runner.output(df, outpath, mode="w")
    for date in tqdm(dates[1::]):
        df, _ = \
            runner.main_day(date, flag="restart", restart="restart.bin",
                            runoffDir=runoffdir, outDir=outdir)
        runner.output(df, outpath, mode="a")


def backup_restart(date, restartpath):
    outpath = "{0}.{1}".format(restartpath,
                               date.strftime("%Y%m%d%H"))
    subprocess.check_call(["cp", restartpath, outpath])


if __name__ == "__main__":
    parser = \
        argparse.ArgumentParser(description="HRR single run experiment")
    parser.add_argument("-c", "--config", type=str, help="config file",
                        required=True)
    parser.add_argument("-s", "--start", type=str,
                        required=True,
                        help="start date in YYYYMMDD")
    parser.add_argument("-e", "--end", type=str,
                        required=True,
                        help="end date in YYYYMMDD")
    args = parser.parse_args()
    if not os.path.exists(args.config):
        raise(FileNotFoundError("config file does not exists."))
    sdate = datetime.datetime.strptime(args.start,
                                       "%Y%m%d")
    edate = datetime.datetime.strptime(args.end,
                                       "%Y%m%d")
    run(args.config, sdate, edate)
