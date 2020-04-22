import os
import argparse
import numpy as np
import pandas as pd
import datetime
import configparser
import subprocess
from tqdm import tqdm
import pyHRR as pyHRR


class DA_pyHRR(object):
    def __init__(self, config):

        self.config = config
        self.initFile = configparser.ConfigParser()
        self.initFile.read(config)

        # general settings
        self.expName = self.initFile.get("experiment", "expName")
        self.undef = self.initFile.get("observation", "undef")
        self.rootDir = self.initFile.get("model", "rootDir")
        self.runoffDir = self.initFile.get("IO", "runoffDir")
        self.outDir = os.path.join(self.initFile.get("IO", "outDir"),
                                   self.expName)
        self.oFile = self.initFile.get("IO", "oFile")
        self.oFileInsert = self.initFile.get("IO", "oFileInsert")
        self.oName = os.path.join(self.outDir, self.oFile)
        self.oNameInsert = os.path.join(self.outDir, self.oFileInsert)
        self.compile_ = False

        # modeling variables
        self.nReach = int(self.initFile.get("model", "nReach"))
        self.reachStart = 1
        self.ndx = int(self.initFile.get("input", "ndx"))
        self.channelPath = self.initFile.get("model", "channelPath")
        self.channelInfo = pd.read_csv(self.channelPath,
                                       sep="\s+",
                                       header=None)
        self.spinupPath = os.path.join(self.rootDir,
                                       self.outDir,
                                       "00/restart.bin")
        self.undef = -9999
        # initialize pyHRR
        self.model = pyHRR.HRR(config, compile_=self.compile_)
        self.outReachID = self.model.read_OutID()  # HRRID
        self.outReachIdx = self.outReachID - self.reachStart  # idx
        # read obs
        self.obsDf = self.readObs()
        self.insertDates = self.obsDf.index
        # load upstream reach info.
        nupa = [3, 4, 5, 6]  # columns containing upstream reach info.
        upas = []
        for i in nupa:
            upa = self.channelInfo.iloc[:, i]
            upa = upa.astype(np.int64).values.reshape(1, -1)
            upas.append(upa - self.reachStart)
        self.upas = np.concatenate(upas, axis=0)

    def spinup(self, sDate, eDate, model):
        runoffDir = self.runoffDir
        storageDir = \
            os.path.join(self.outDir,
                         "spinup/{0}".format(sDate.strftime("%Y%m%d%H")))
        if not os.path.exists(storageDir):
            print("makedirs")
            os.makedirs(storageDir)
        out, nDate = model.main_day(sDate,
                                    flag="initial",
                                    restart="restart.bin",
                                    runoffDir=runoffDir,
                                    outDir=storageDir)
        date = nDate
        print("start spinning up...")
        while date < eDate:
            print(date)
            out, nDate = model.main_day(date,
                                        flag="restart",
                                        restart="restart.bin",
                                        runoffDir=runoffDir,
                                        outDir=storageDir)
            date = nDate
        print("end spinup")

        return os.path.join(storageDir, "restart.bin")

    def startFromInit(self, sDate, eDate,
                      spinup=False, sDate_spin=None, eDate_spin=None,
                      verbose=False):
        if spinup:
            assert (isinstance(sDate_spin, datetime.datetime))
            assert (isinstance(eDate_spin, datetime.datetime))
            spinupPath = self.spinup(sDate_spin,
                                     eDate_spin,
                                     self.model)
        else:
            spinupPath = self.spinupPath
        sdate = self.initializeFromRestart(self.model,
                                           sDate,
                                           spinupPath,
                                           verbose)

        dates = pd.date_range(sdate, edate, freq="D")
        print("start main run.")
        if verbose:
            itr = tqdm(dates)
            for date in itr:
                self.forward(date)
        else:
            itr = dates
            for date in itr:
                print(date)
                self.forward(date)

    def restartFromRecent(self, sdate, edate, verbose=False):
        dates = pd.date_range(sdate, edate, freq="D")
        print("start main run.")
        if verbose:
            itr = tqdm(dates)
            for date in itr:
                self.forward(date)
        else:
            itr = dates
            for date in itr:
                print(date)
                self.forward(date)

    def forward(self, date):
        # forwarding step
        self.simulation_core(date)
        # assimilation if applicable
        if (date == self.insertDates).any():
            xa = self.dataInsert(date)
            self.insertOut(xa[self.outReachIdx], date)
        if date.month == 12 and date.day == 31:
            self.backupRestart(date)

    def initializeFromRestart(self, model, sDate, spinupPath, verbose):
        restartFile = spinupPath
        print("initializing...")
        outDir = self.outDir
        if not os.path.exists(outDir):
            os.makedirs(outDir)
        subprocess.check_call(["cp", restartFile, outDir])
        df, nDate = model.main_day(sDate,
                                   flag="restart",
                                   restart="restart.bin",
                                   runoffDir=self.runoffDir,
                                   outDir=self.outDir)
        subprocess.check_call(["cp",
                               restartFile,
                               os.path.join(outDir,
                                            "restartInsert.bin")])
        # model.output(df, self.oName, mode="w")
        model.output(df, self.oNameInsert, mode="w")
        print("end initializing")
        return nDate

    def dataInsert(self, date, cache=False):
        obsMean = self.__constObs(self.obsDf, date)
        print(obsMean[obsMean!=self.undef])
        self.restarts = []
        dschg_cfs = self.__readRestart("restartInsert.bin")
        xa = np.where(obsMean==self.undef,
                      dschg_cfs,
                      obsMean)
        innovation = xa[0] - dschg_cfs[0]
        print(innovation[innovation!=0])
        self.updateChannel(innovation)
        return xa[0]

    def __readRestart(self, restartFile):
        restPath = os.path.join(self.outDir, restartFile)
        rest = np.memmap(restPath, mode="r+",
                         shape=(8, self.ndx, self.nReach), dtype=np.float32)
        self.restart = rest
        old_qs = rest[-1, -1, :].copy()  # cfs
        return old_qs.reshape(1, -1)

    def __constObs(self, obsDf, date, infl=1):

        obs = obsDf[obsDf.index == date]
        # python index starts from 0
        reaches = obs.reach.values - self.reachStart
        obsMean = np.ones([self.nReach]) * self.undef
        obsMean[reaches] = self.cms2cfs(obs["mean"].values)
        obsMean = obsMean.reshape(1, -1)
        return obsMean

    def cfs2cms(self, cfs):
        cms = cfs*(0.3048**3)
        return cms

    def cms2cfs(self, cms):
        cfs = cms/(0.3048**3)
        cfs_undef = self.undef/(0.3048**3)
        cfs[np.where(cfs-cfs_undef < 0.0001)] = self.undef
        return cfs

    def readObs(self):
        obsPath = self.initFile.get("observation", "obsPath")
        df = pd.read_csv(obsPath, parse_dates=[0], index_col=[0])
        return df

    def backupRestart(self, date):
        restFile = os.path.join(self.outDir,
                                "restart.bin")
        outRestFile = restFile + ".{0}".format(date.strftime("%Y%m%d"))
        restInsertFile = os.path.join(self.outDir,
                                      "restartInsert.bin")
        outRestInsertFile = restInsertFile + \
            ".{0}".format(date.strftime("%Y%m%d"))
        # subprocess.call(["cp", restFile, outRestFile])
        subprocess.call(["cp", restInsertFile, outRestInsertFile])

    def submit(self, sdate, edate, restart, days_spinup=5, verbose=False):
        if restart:
            self.restartFromRecent(sdate, edate, verbose)
        else:
            sDate_spin = sdate
            eDate_spin = sdate + datetime.timedelta(days=days_spinup)
            sDate = eDate_spin
            eDate = edate
            self.startFromInit(sDate, eDate,
                               spinup=True,
                               sDate_spin=sDate_spin,
                               eDate_spin=eDate_spin,
                               verbose=verbose)

    def simulation_core(self, date):
        # df, nDate = self.model.main_day(date,
        #                            flag="restart",
        #                            restart="restart.bin",
        #                            runoffDir=self.runoffDir,
        #                            outDir=self.outDir)
        adf, nDate = self.model.main_day(date,
                                         flag="restart",
                                         restart="restartInsert.bin",
                                         runoffDir=self.runoffDir,
                                         outDir=self.outDir)
        # self.model.output(df, oName, mode="a")
        self.model.output(adf, self.oNameInsert, mode="a")
        return nDate

    def updateChannel(self, aArray1d):
        rest = self.restart
        # [1,2,3,4] => [1,1,1,...,1,2,2,2,...,2,...]
        update_innovation = np.repeat(aArray1d, self.ndx).reshape(self.nReach,
                                                                  self.ndx).T
        print(update_innovation)
        # before = rest.copy()
        rest[0, :, :] = rest[0, :, :] + update_innovation
        rest[1, :, :] = rest[1, :, :] + update_innovation
        rest[2, :, :] = rest[2, :, :] + update_innovation
        # assert before[0, :, :].sum()!=rest[0, :, :].sum(), "not updated"
        del rest

    def insertOut(self, xa, date):
        df = pd.DataFrame(np.vectorize(self.cfs2cms)(xa))
        df = df.T
        df.index = [date]
        df.reset_index().rename({"index": "Date"}, axis=1).set_index("Date")
        self.model.output(df, self.oNameInsert, mode="a")


if __name__ == "__main__":
    parser = \
        argparse.ArgumentParser(description="HRR Data Assimilation experiment")
    parser.add_argument("-c", "--config", type=str, help="config file",
                        required=True)
    parser.add_argument("-s", "--start", type=str,
                        required=True,
                        help="start date in YYYYMMDD")
    parser.add_argument("-e", "--end", type=str,
                        required=True,
                        help="end date in YYYYMMDD")
    parser.add_argument("--restart", action="store_true",
                        default=False,
                        help="restarting from existing restart.bin?")
    parser.add_argument("--spdays", type=int, default=5,
                        help="days for spinup")
    parser.add_argument("-v", "--verbose", default=False,
                        action="store_true",
                        help="show progress bar")
    args = parser.parse_args()
    if not os.path.exists(args.config):
        raise(FileNotFoundError("config file does not exists."))
    sdate = datetime.datetime.strptime(args.start,
                                       "%Y%m%d")
    edate = datetime.datetime.strptime(args.end,
                                       "%Y%m%d")
    runner = DA_pyHRR(args.config)
    if args.verbose:
        from tqdm import tqdm
    runner.submit(sdate, edate, args.restart, args.spdays, args.verbose)
