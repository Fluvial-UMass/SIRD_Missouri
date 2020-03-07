import os
import argparse
import numpy as np
import pandas as pd
import datetime
import configparser
import subprocess
import pickle
from multiprocessing import Pool
import pyletkf
import pyHRR as pyHRR
import daTools


class DA_pyHRR(object):
    def __init__(self, config):

        self.config = config
        self.initFile = configparser.ConfigParser()
        self.initFile.read(config)

        # general settings
        self.expName = self.initFile.get("experiment", "expName")
        self.undef = self.initFile.get("observation", "undef")
        self.rootDir = self.initFile.get("model", "rootDir")
        self.runoffDir = os.path.join(self.initFile.get("IO", "runoffDir"),
                                      "{0:02d}")
        self.oDir = os.path.join(self.initFile.get("IO", "outDir"),
                                 self.expName, "{0:02d}")
        self.assimCacheDir = self.initFile.get("assimilation", "assimCacheDir")
        self.oFile = self.initFile.get("IO", "oFile")
        self.oFileAssim = self.initFile.get("IO", "oFileAssim")
        self.oName = os.path.join(self.oDir, self.oFile)
        self.oNameAssim = os.path.join(self.oDir, self.oFileAssim)
        self.compile_ = False
        self.eTot = int(self.initFile.get("assimilation", "ensMem"))

        # modeling variables
        self.nReach = int(self.initFile.get("model", "nReach"))
        self.reachStart = 1
        self.ndx = int(self.initFile.get("input", "ndx"))
        self.channelPath = self.initFile.get("model", "channelPath")
        self.channelInfo = pd.read_csv(self.channelPath,
                                       sep="\s+",
                                       header=None)
        self.daCore = pyletkf.LETKF_core(config,
                                         mode="vector",
                                         use_cache=True)
        self.spinupPath = os.path.join(self.rootDir,
                                       self.oDir,
                                       "00/restart.bin")
        self.nCpus = int(self.initFile.get("model", "nCpus"))
        # initialize pyletkf
        self.daCore.initialize(backend="h5py")
        # initialize pyHRR
        self.model = pyHRR.HRR(config, compile_=self.compile_)
        self.outReachID = self.model.read_OutID()  # HRRID
        self.outReachIdx = self.outReachID - self.reachStart  # idx
        # read obs
        self.obsDf = self.readObs()
        self.assimDates = self.obsDf.index
        # load upstream reach info.
        nupa = [3, 4, 5, 6]  # columns containing upstream reach info.
        upas = []
        for i in nupa:
            upa = self.channelInfo.iloc[:, i]
            upa = upa.astype(np.int64).values.reshape(1, -1)
            upas.append(upa - self.reachStart)
        self.upas = np.concatenate(upas, axis=0)

    def spinup(self, sDate, eDate, model, eNum=0):
        runoffDir = self.runoffDir.format(eNum)
        storageDir = \
            os.path.join(self.oDir.format(eNum),
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
                                     self.model,
                                     eNum=0)
        else:
            spinupPath = self.spinupPath
        sdate = self.initializeFromRestart(self.eTot,
                                           self.model,
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
        simulation_parallel(date, self.eTot, self.model,
                            self.runoffDir, self.oDir, self.oName,
                            self.oNameAssim, ncpus=self.nCpus)
        # assimilation if applicable
        if (date == self.assimDates).any():
            xa = self.dataAssim(date, self.obsDf, self.eTot)
            assimOut_parallel(xa[:, :, self.outReachIdx], date, self.model,
                              self.cfs2cms, self.oNameAssim,
                              self.eTot, ncpus=self.nCpus)
        if date.month == 12 and date.day == 31:
            self.backupRestart(date)

    def initializeFromRestart(self, eTot, model, sDate, spinupPath, verbose):
        restartFile = spinupPath
        print("initializing...")
        if verbose:
            itr = tqdm(range(eTot))
        else:
            itr = range(eTot)
        for eNum in itr:
            outDir = os.path.join(self.rootDir, self.oDir.format(eNum))
            if not os.path.exists(outDir):
                os.makedirs(outDir)
            subprocess.check_call(["cp", restartFile, outDir])
            runoffDir = self.runoffDir.format(eNum)
            df, nDate = model.main_day(sDate,
                                       flag="restart",
                                       restart="restart.bin",
                                       runoffDir=runoffDir,
                                       outDir=outDir)
            subprocess.check_call(["cp",
                                   restartFile,
                                   os.path.join(outDir,
                                                "restartAssim.bin")])
            model.output(df, self.oName.format(eNum), mode="w")
            model.output(df, self.oNameAssim.format(eNum), mode="w")
        print("end initializing")
        return nDate

    def dataAssim(self, date, obsDf, eTot, cache=False, limitval=1e+10):
        obsMean, obsStd = self.__constObs(obsDf, date)
        self.restarts = []
        # take logs
        dschg_cfs_ens = \
            np.concatenate([self.__readRestart(eNum, "restartAssim.bin")
                            for eNum in range(eTot)], axis=0)
        dschg_cfs_ens = \
            dschg_cfs_ens.astype(np.float64).reshape(1, eTot, 1, self.nReach)
        dschg_cfs_ens = np.where(dschg_cfs_ens == 0, 1e-8, dschg_cfs_ens)
        dschg_cfs_ens = np.log(dschg_cfs_ens)
        obsvars = [1]
        # As an API requirements,
        # the input simulation array should be (nvars, eTot, nT, nReach)
        xa, Ws = self.daCore.letkf_vector(dschg_cfs_ens,
                                          obsMean,
                                          obsStd.astype(np.float64),
                                          obsvars,
                                          guess="prior",
                                          nCPUs=self.nCpus)
        xa = xa[0, :, :, :]
        if cache:
            outname = date.strftime("%Y%m%d%h_Ws.pkl")
            outPath = os.path.join(self.assimCacheDir, outname)
            with open(outPath, "wb") as f:
                pickle.dump(Ws, f)
        # limiter: to avoid diversion.
        xa[xa > np.log(limitval)] = dschg_cfs_ens[0][xa > np.log(limitval)]
        xa = np.exp(xa)  # convert from log
        updateChannel_parallel(xa, self.ndx, self.upas,
                               self.nReach, self.restarts,
                               eTot, ncpus=self.nCpus)
        return xa

    def take_nLog(self, array):
        """
        For observation data which contains missing values
        """
        array = np.where(array == 0, 1e-8, array)
        outArray = np.log(array)
        outArray[np.where(array == self.daCore.undef)] = self.daCore.undef
        return outArray

    def __readRestart(self, eNum, restartFile):
        restPath = os.path.join(self.oDir.format(eNum), restartFile)
        rest = np.memmap(restPath, mode="r+",
                         shape=(4, self.ndx, self.nReach), dtype=np.float32)
        self.restarts.append(rest)
        old_qs = rest[0, -1, :].copy()  # cfs
        return old_qs.reshape(1, -1)

    def __constObs(self, obsDf, date, infl=10):

        obs = obsDf[obsDf.index == date]
        # python index starts from 0
        reaches = obs.reach.values - self.reachStart
        obsMean = np.ones([self.nReach]) * self.daCore.undef
        obsStd = np.ones([self.nReach]) * 0.01

        # take logs
        obsMean[reaches] = self.take_nLog(self.cms2cfs(obs["mean"].values))
        obsConfLow = self.take_nLog(self.cms2cfs(obs["conf.low"].values))
        obsConfUpp = self.take_nLog(self.cms2cfs(obs["conf.high"].values))
        obsStd[reaches] = infl*(obsConfUpp - obsConfLow) / (2*1.96)
        print(obsMean[reaches])
        print(obsStd[reaches])
        obsMean = obsMean.reshape(1, -1)
        obsStd = obsStd.reshape(1, -1)
        return obsMean, obsStd

    def cfs2cms(self, cfs):
        cms = cfs*(0.3048**3)
        return cms

    def cms2cfs(self, cms):
        cfs = cms/(0.3048**3)
        cfs_undef = self.daCore.undef/(0.3048**3)
        cfs[np.where(cfs-cfs_undef < 0.0001)] = self.daCore.undef
        return cfs

    def readObs(self):
        obsPath = self.initFile.get("observation", "obsPath")
        df = pd.read_csv(obsPath, parse_dates=[0], index_col=[0])
        return df

    def backupRestart(self, date):
        for eNum in range(self.eTot):
            restFile = os.path.join(self.oDir.format(eNum),
                                    "restart.bin")
            outRestFile = restFile + ".{0}".format(date.strftime("%Y%m%d"))
            restAssimFile = os.path.join(self.rootDir,
                                         self.oDir.format(eNum),
                                         "restartAssim.bin")
            outRestAssimFile = restAssimFile + \
                ".{0}".format(date.strftime("%Y%m%d"))
            subprocess.call(["cp", restFile, outRestFile])
            subprocess.call(["cp", restAssimFile, outRestAssimFile])

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


# For multi-processing
def simulation_core(argList):
    date = argList[0]
    eNum = argList[1]
    model = argList[2]
    runoffDir = os.path.join(argList[3].format(eNum))
    outDir = os.path.join(argList[4].format(eNum))
    oName = argList[5].format(eNum)
    oNameAssim = argList[6].format(eNum)
    df, nDate = model.main_day(date,
                               flag="restart",
                               restart="restart.bin",
                               runoffDir=runoffDir,
                               outDir=outDir)
    adf, nDate = model.main_day(date,
                                flag="restart",
                                restart="restartAssim.bin",
                                runoffDir=runoffDir,
                                outDir=outDir)
    model.output(df, oName, mode="a")
    model.output(adf, oNameAssim, mode="a")
    return nDate


def simulation_parallel(date, eTot, model,
                        runoffRoot, outRoot, oName,
                        oNameAssim, ncpus=20):
    args = [[date, eNum, model, runoffRoot, outRoot, oName, oNameAssim]
            for eNum in range(eTot)]
    results = []
    pool = Pool(processes=ncpus)
    results = pool.map(simulation_core, args)
    pool.close()
    pool.join()
    return results[0]


def updateChannel(argsList):
    aArray1d = argsList[0].astype(np.float32)
    ndx = argsList[1]
    upas = argsList[2].astype(np.int32)
    nReach = argsList[3]
    rest = argsList[4]
    # [1,2,3,4] => [1,1,1,...,1,2,2,2,...,2,...]
    update_qs = np.repeat(aArray1d, ndx).reshape(nReach, ndx).T
    update_qin, update_qout = \
        daTools.updateQinout_pa(aArray1d, upas, nReach, ndx)
    rest[0, :, :] = update_qs
    rest[1, :, :] = update_qin
    rest[2, :, :] = update_qout
    del rest


def updateChannel_parallel(xa, ndx, upas,
                           nReach, restarts,
                           eTot, ncpus=20):
    args = [[xa[eNum, -1, :], ndx, upas,
             nReach, restarts[eNum]] for eNum in range(eTot)]
    for i in range(eTot):
        updateChannel(args[i])
    # slow due to IO conflict?
    # pool = Pool(processes=ncpus)
    # pool.map(updateChannel, args)
    # pool.close()
    # pool.join()


def assimOut(argsList):
    xa_each = argsList[0]
    eNum = argsList[1]
    date = argsList[2]
    modelInstance = argsList[3]
    cfs2cms = argsList[4]
    oNameAssim = argsList[5]
    df = pd.DataFrame(np.vectorize(cfs2cms)(xa_each))
    df = df.T
    df.index = [date]
    df.reset_index().rename({"index": "Date"}, axis=1).set_index("Date")
    modelInstance.output(df, oNameAssim.format(eNum), mode="a")


def assimOut_parallel(xa, date, modelInstance,
                      cfs2cms, oNameAssim, eTot, ncpus=20):
    args = [[xa[eNum, -1], eNum, date, modelInstance,
             cfs2cms, oNameAssim] for eNum in range(eTot)]
    for i in range(eTot):
        assimOut(args[i])
    # slow due to IO conflict?
    # pool = Pool(processes=ncpus)
    # pool.map(assimOut, args)
    # pool.close()
    # pool.join()


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
