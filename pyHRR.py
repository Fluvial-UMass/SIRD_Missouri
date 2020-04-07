import configparser
import datetime
import os
import subprocess
import sys
try:
    import netCDF4 as nc
except(ImportError):
    pass
import numpy as np
import pandas as pd


class HRR(object):
    def __init__(self, config, compile_=False):
        """"
            THis is a python-wrapper to handle HRR model in conjunction with a data assimilation module.
            Small edition of original source code was performed.
            v.1.0.0: Daily simulation.
            Coding: Yuta Ishitsuka
        """
        # read input info. for HRR.
        self.config = configparser.ConfigParser()
        self.config.read(config)
        self.rootDir = self.config.get("model", "rootDir")

        # file name configuration
        self.srcDir = self.config.get("model", "srcDir")
        self.channelPath = self.config.get("model", "channelPath")
        self.planesPath = self.config.get("model", "planesPath")
        self.exe = os.path.join(self.srcDir, "run")
        self.output_calibration = os.path.join(self.srcDir,
                                               "output_calibration.txt")
        self.backend = "netCDF4"  # output file format
        # for netCDF4 outputs
        self.outputids = self.read_OutID(self.output_calibration)
        self.numpyEpoch = np.datetime64('1970-01-01')

        # compile or not
        if compile_:
            self.__compile()

    def __compile(self):

        print("compilation is activated. Make clean/all.")
        os.chdir(self.srcDir)
        print("make clean")
        subprocess.check_call(["make", "clean"])
        print("make all")
        subprocess.check_call(["make", "all"])
        os.chdir(self.rootDir)

    def __createInputFile(self, date, runoffDir, restart, mode, outDir):

        assimUpdate = mode
        roffData = os.path.join(runoffDir,
                                "%s.txt" % (date.strftime("%Y%m%d")))
        pfafunits = self.config.get("input", "pfafunits")
        ndx = self.config.get("input", "ndx")
        ndt = 24  # for future improvement
        dtis = 3600  # for future improvement
        self.outerDt = dtis * ndt
        iyear = date.year
        imonth = date.month
        iday = date.day
        doy = (date - datetime.datetime(date.year, 1, 1, 0)).days + 1
        ks_all = self.config.get("input", "ks_all")
        sbRate = self.config.get("input", "sbRate")
        n_ch_all = self.config.get("input", "n_ch_all")
        qlat_all = self.config.get("input", "qlat_all")

        VARS = [
            "pfafunits", "ndx", "ndt", "dtis", "iyear", "imonth", "iday",
            "Julian Day", "ks_all", "setfsub_rate", "n_ch_all", "qlat_all"
        ]
        VALS = {
            "pfafunits": pfafunits,
            "ndx": ndx,
            "ndt": ndt,
            "dtis": dtis,
            "iyear": iyear,
            "imonth": imonth,
            "iday": iday,
            "Julian Day": doy,
            "ks_all": ks_all,
            "setfsub_rate": sbRate,
            "n_ch_all": n_ch_all,
            "qlat_all": qlat_all
        }

        self.infoPath = os.path.join(outDir, "input.txt")
        with open(self.infoPath, mode="w") as f:
            f.write("%s\n" % assimUpdate)
            f.write("%s\n" % self.srcDir)
            f.write("%s\n" % roffData)
            f.write("%s\n" % self.channelPath)
            f.write("%s\n" % self.planesPath)
            f.write("%s\n" % restart)
            f.write("%s\n" % outDir)
            [f.write("%s    %s\n" % (str(VALS[var]), var)) for var in VARS]

    def __create_output(outpath):
        if self.backend == "netCDF4":
            self.__create_netcdf(outpath)
        elif self.backend == "csv":
            pass
        else:
            raise KeyError("format {0} is not supported.".format(backend))

    def __create_netcdf(self, outpath):
        outnc = nc.Dataset(outpath, mode="w", format="NETCDF4")
        # dimensions
        time = outnc.createDimension("time", None)
        reach = outnc.createDimension("reach",
                                      len(self.outputids))
        # variables
        times = outnc.createVariable("time", "f8", ("time",))
        reaches = outnc.createVariable("reach", "i4", ("reach",))
        discharge = outnc.createVariable("discharge", "f4", ("time",
                                                             "reach",))
        # attributes
        outnc.description = "HRR discharge output"
        now = datetime.datetime.utcnow().strftime("%Y-%m:%d:%H")
        outnc.history = "Created {0}:UTC".format(now)
        reaches.units = "HRRID"
        discharge.units = "cms"
        times.units = "hours since 0001-01-01 00:00:00.0"
        times.calendar = "gregorian"
        # add variables
        reaches[:] = self.outputids
        # close and flush changes
        outnc.close()

    def output_nc(self, df, outpath, mode="a"):
        if mode == "w":
            self.__create_netcdf(outpath)
        outnc = nc.Dataset(outpath, mode="a")
        timeidx = outnc["/discharge"].shape[0]
        # get datetime object from np.datetime64
        npdate = df.index.values[0]
        tstamp = (npdate - self.numpyEpoch) / np.timedelta64(1, 's')
        date = datetime.datetime.utcfromtimestamp(tstamp)
        outputs = df.values
        outnc["/discharge"][timeidx, :] = outputs
        outnc["/time"][timeidx] = \
            nc.date2num(date,
                        units=outnc["/time"].units,
                        calendar=outnc["/time"].calendar)
        outnc.close()

    def output_csv(self, df, outpath, mode="a"):
        if mode == "w":
            with open(outpath, "w") as f:
                df = df.reset_index().rename({"index": "Date"}, axis=1)
                df.to_csv(f, index=False)
        elif mode == "a":
            with open(outpath, "a") as f:
                df.to_csv(f, header=False)
        else:
            raise KeyError("mode %s is unsupported." % mode)

    def output(self, df, outpath, mode="a"):
        if self.backend == "netCDF4":
            self.output_nc(df, outpath, mode=mode)
        elif self.backend == "csv":
            self.output_csv(df, outpath, mode=mode)
        else:
            raise KeyError("format {0} is not supported.".format(backend))

    def main_day(self,
                 date,
                 flag="restart",
                 restart="restart.txt",
                 runoffDir="../data/case1/",
                 mode="normal",
                 outDir="./out"):
        """
            main API to handle HRR
        """

        # check operation mode
        if not flag == "initial" and not flag == "restart":
            raise IOError("Undefined flag mode %s" % flag)

        # create input information file for HRR
        self.__createInputFile(date, runoffDir, restart, mode, outDir)

        # activate HRR
        subprocess.check_call([self.exe, flag, self.infoPath])
        nDate = date + datetime.timedelta(seconds=self.outerDt)
        #print('hrr output is:', os.path.join(outDir, "discharge_cms.txt"))
        out = pd.read_csv(os.path.join(outDir, "discharge_cms.txt"),
                          header=0,
                          sep="\s+").rename(index={0: date})

        return out, nDate

    def read_OutID(self, file="output_calibration.txt"):
        path = os.path.join(self.srcDir, file)
        with open(path, "r") as f:
            lines = f.readlines()
            numout = int(lines[0])
            outids = []
            for i in range(1, numout+1):
                outids.append(int(lines[i].split()[0]))
        return np.array(outids)


if __name__ == "__main__":
    import time
    stime = time.time()
    test = HRR("./config/config_Mackenzie_v2_case1.ini")
    runoffdir = "/project/uma_colin_gleason/Dongmei/DataAssim/Mackenzie/data/case1/00"
    sdate = datetime.datetime(1984, 10, 1)
    edate = datetime.datetime(1984, 10, 5)
    df, ndate = \
        test.main_day(sdate, flag="initial", restart="restart.bin",
                      runoffDir=runoffdir)
    test.output(df, "test.nc", mode="w")
    date = ndate
    while date < edate:
        print(date)
        df, date = \
            test.main_day(date, flag="restart", restart="restart.bin",
                          runoffDir=runoffdir)
        test.output(df, "test.nc", mode="a")
    print("elapsed:", time.time()-stime)
