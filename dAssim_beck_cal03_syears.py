import os
import sys

import numpy as np
import pandas as pd
import datetime
import configparser
import subprocess
import pickle
from multiprocessing import Pool

import parallel.pyletkf.pyletkf.pyletkf as pyletkf
import pytHRR as pyHRR
import daTools

class DA_pyHRR(object):
    def __init__(self,config):

        self.config = config
        self.initFile = configparser.ConfigParser()
        self.initFile.read(config)

        self.undef = self.initFile.get("observation","undef")
        self.rootDir = self.initFile.get("model","rootDir")
        self.runoffDir = "/nl/uma_colin_gleason/yuta/gridAssign/src/daily/meritvic_shyears/corrected/%02d/"
        self.oDir = "/project/uma_colin_gleason/yuta/DA/missouri/pyHRR/out/MERIT_VIC_beck_cal03_syears/%02d"
        self.assimCacheDir = "/project/uma_colin_gleason/yuta/DA/missouri/pyHRR/cache/"
        self.oFile = "discharge.csv"
        self.oFileAssim = "dischargeAssim.csv"
        self.oFileAssimSmooth = "dischargeAssim_smooth.csv"
        self.oName = os.path.join(self.oDir,self.oFile)
        self.oNameAssim = os.path.join(self.oDir,self.oFileAssim)
        self.oNameAssimSmooth = os.path.join(self.oDir,self.oFileAssimSmooth)

        self.compile_ = False

        self.eTot = int(self.initFile.get("assimilation","ensMem"))
        self.smoothMax = 365*2 #depends on your memory size and your observation frequency

        # modeling variables
        self.nReach = int(self.initFile.get("model","nReach"))
        self.reachStart = 1
        self.ndx = 20
        self.channelPath =os.path.join(self.rootDir,"hillslope_src","channels_cal_1984_2002.txt")
        self.channelInfo = pd.read_csv(self.channelPath,sep="\s+",header=None)
        self.daCore = pyletkf.LETKF_core(config,mode="vector",use_cache=True)
        self.spinupPath = os.path.join(self.rootDir,self.oDir,"00/restart.txt")
        self.nCpus = int(self.initFile.get("model","nCpus"))
        #initialize pyletkf
        self.daCore.initialize(backend="pickle")
        #initialize pyHRR
        self.model = pyHRR.HRR(config, compile_=self.compile_)
        #read obs
        self.obsDf = self.readObs()
        self.assimDates = self.obsDf.index
        #load upstream reach info.
        nupa = [3,4,5,6] # number of columns that contains upstream reach info.
        upas = []
        for i in nupa:
            upa = self.channelInfo.iloc[:,i].astype(np.int64).values.reshape(1,-1)
            upas.append(upa-1) # reach number starts from 1
        self.upas = np.concatenate(upas,axis=0)

    def spinup(self, sDate, eDate, model, eNum=0):
        runoffDir = os.path.join(self.rootDir,self.runoffDir%eNum)
        storageDir = os.path.join(self.rootDir,self.oDir%eNum,"spinup/%s"%sDate.strftime("%Y%m%d%H"))
        if os.path.exists(storageDir) == False:
            print("makedirs")
            os.makedirs(storageDir)
        out, nDate = model.main_day(sDate, flag="initial", restart="restart.txt", runoffDir=runoffDir, outDir=storageDir)
        date = nDate
        while date < eDate:
            out, nDate = model.main_day(date, flag="restart", restart="restart.txt", runoffDir=runoffDir, outDir=storageDir)
            date = nDate

        return os.path.join(storageDir, "restart.txt")

    def mainRun(self, sDate, eDate, spinup=False, sDate_spin=None, eDate_spin=None):
        if spinup:
            assert (isinstance(sDate_spin,datetime.datetime))
            assert (isinstance(eDate_spin,datetime.datetime))
            spinupPath = self.spinup(sDate_spin, eDate_spin, self.model, eNum=0)
        else:
            spinupPath = self.spinupPath
        date = self.initializeFromRestart(self.eTot, self.model, sDate, spinupPath)
        rnofRt = os.path.join(self.rootDir,self.runoffDir)
        outRt = os.path.join(self.rootDir,self.oDir)
        while date < eDate:
            nDate = simulation_parallel(date, self.eTot, self.model, rnofRt, outRt, self.oName, self.oNameAssim, ncpus=self.nCpus)
            if (date == self.assimDates).any():
                xa = self.dataAssim(date, self.obsDf, self.eTot)
                print(xa[:, -1, 20809])
                # [self.__assimOut(xa[eNum, -1],eNum,date,self.model) for eNum in range(0, self.eTot)]
                assimOut_parallel(xa, date, self.model, self.cfs2cms, self.oNameAssim, self.eTot, ncpus=self.nCpus)
            if date.month == 12 and date.day == 31:
                self.backupRestart(date)
            date = nDate

    def initializeFromRestart(self,eTot,model,sDate,spinupPath):
        restartFile = spinupPath
        for eNum in range(0,eTot):
            outDir = os.path.join(self.rootDir,self.oDir%eNum)
            if os.path.exists(outDir) == False:
                os.makedirs(outDir)
            subprocess.check_call(["cp",restartFile,outDir])
            runoffDir = os.path.join(self.rootDir,self.runoffDir%eNum)
            df, nDate = model.main_day(sDate, flag="restart", restart="restart.txt", runoffDir=runoffDir, outDir=outDir) # This will create the most recent restart.tx, and
            subprocess.check_call(["cp",restartFile,os.path.join(outDir,"restartAssim.txt")]) # copy that for the next simulation (assim).
            model.output(df,self.oName%eNum,mode="w")
            model.output(df,self.oNameAssim%eNum,mode="w")
        return nDate

    def dataAssim(self, date, obsDf, eTot):
        obsMean, obsStd = self.__constObs(obsDf,date)
        self.restarts = []
        dschg_cfs_ens = np.concatenate([self.__readRestart(eNum) for eNum in range(0,eTot)], axis=0)
        dschg_cfs_ens = dschg_cfs_ens.astype(np.float64).reshape(1,eTot,1,self.nReach)
        dschg_cfs_ens = np.where(dschg_cfs_ens == 0, 1e-8, dschg_cfs_ens)
        dschg_cfs_ens = np.log(dschg_cfs_ens)
        obsvars = [1]
        # As an API requirements, the input simulation array should be (nvars, eTot, nT, nReach)
        xa, Ws = self.daCore.letkf_vector(dschg_cfs_ens,obsMean,obsStd.astype(np.float64),obsvars,nCPUs=self.nCpus)
        xa = xa[0, :, :, :]
        outname = date.strftime("%Y%m%d%h_Ws.pkl")
        # outPath = os.path.join(self.assimCacheDir, outname)
        # with open(outPath,"wb") as f:
        #    pickle.dump(Ws,f)
        #xa[np.where(xa<0.)] = np.absolute(xa[np.where(xa<0.)])
        xa[xa > np.log(1e+10)] = dschg_cfs_ens[0][xa > np.log(1e+10)]  # limiter: to avoid diversion.
        xa = np.exp(xa) # convert from log
        # [self.__updateChannel(xa[eNum,-1,:],self.restarts[eNum],eNum) for eNum in range(0,eTot)]
        updateChannel_parallel(xa, self.ndx, self.upas, self.nReach, self.restarts, self.rootDir, self.oDir, eTot, ncpus=self.nCpus)
        return xa

    def take_nLog(self, array):
        """
        For observation data which contains missing values
        """
        array = np.where(array==0, 1e-8, array)
        outArray = np.log(array)
        outArray[np.where(array==self.daCore.undef)] = self.daCore.undef
        return outArray

    def __assimOut(self, xa_each, eNum, date, modelInstance):
        df = pd.DataFrame(xa_each).apply(self.cfs2cms, axis=1)
        df = df.T
        df.index = [date]
        df.reset_index().rename({"index":"Date"}, axis=1).set_index("Date")
        modelInstance.output(df, self.oNameAssim%eNum, mode="a")

    def __readRestart(self,eNum):
        dfPath = os.path.join(self.rootDir,self.oDir%eNum,"restartAssim.txt")
        edf = pd.read_csv(dfPath)
        self.restarts.append(edf)
        old_qs = (edf.groupby("i").max())["old_q"].values #cfs
        return old_qs.reshape(1,-1)

    def __constObs(self,obsDf,date):

        obs = obsDf[obsDf.index==date]
        print(obs)
        reaches = obs.reach.values - 1 #for python which starts from 0
        print(reaches)
        obsMean = np.ones([self.nReach]) * self.daCore.undef
        obsStd  = np.ones([self.nReach]) * 0.01

        obsMean[reaches] = self.take_nLog(self.cms2cfs(obs["mean"].values)) # log
        obsConfLow = self.take_nLog(self.cms2cfs(obs["conf.low"].values))
        obsConfUpp = self.take_nLog(self.cms2cfs(obs["conf.high"].values))
        obsStd[reaches] = (obsConfUpp - obsConfLow) / (2*1.96)
        
        obsMean = obsMean.reshape(1, -1)
        obsStd = obsStd.reshape(1, -1)
        return obsMean, obsStd

    def cfs2cms(self,cfs):
        cms = cfs*0.3048**3
        return cms

    def cms2cfs(self,cms):
        cfs = cms/(0.3048**3)
        cfs_undef = self.daCore.undef/(0.3048**3)
        cfs[np.where(cfs-cfs_undef < 0.0001)] = self.daCore.undef
        return cfs

    def readObs(self):

        obsPath = self.initFile.get("observation","obsPath")
        df = pd.read_csv(obsPath,parse_dates=[0],index_col=[0])

        return df

    def backupRestart(self, date):
        for eNum in range(self.eTot):
            restFile = os.path.join(self.rootDir,self.oDir%eNum,"restart.txt")
            outRestFile = restFile + ".%s" % date.strftime("%Y%m%d")
            restAssimFile = os.path.join(self.rootDir,self.oDir%eNum,"restartAssim.txt")
            outRestAssimFile = restAssimFile + ".%s" % date.strftime("%Y%m%d")
            subprocess.call(["cp", restFile, outRestFile])
            subprocess.call(["cp", restAssimFile, outRestAssimFile])

    def tester(self):
        sDate_spin = datetime.datetime(2002,10,1)
        eDate_spin = datetime.datetime(2002,10,2)
        sDate = eDate_spin
        eDate = datetime.datetime(2011,1,1)
        self.mainRun(sDate, eDate, spinup=True, sDate_spin=sDate_spin, eDate_spin=eDate_spin)

    def restart(self):
        restartDate = datetime.datetime(2002,11,9)
        eDate = datetime.datetime(2011,1,1)
        rnofRt = os.path.join(self.rootDir,self.runoffDir)
        outRt = os.path.join(self.rootDir,self.oDir)
        date = restartDate
        while date < eDate:
            nDate = simulation_parallel(date, self.eTot, self.model, rnofRt, outRt, self.oName, self.oNameAssim, ncpus=self.nCpus)
            if (date == self.assimDates).any():
                xa = self.dataAssim(date, self.obsDf, self.eTot)
                [self.__assimOut(xa[eNum, -1],eNum,date,self.model) for eNum in range(0, self.eTot)]
            if date.month == 12 and date.day == 31:
                self.backupRestart(date)
            date = nDate


# For multi-processing
def simulation_core(argList):
    date = argList[0]
    eNum = argList[1]
    model = argList[2]
    runoffDir = os.path.join(argList[3]%eNum)
    outDir = os.path.join(argList[4]%eNum)
    oName = argList[5]%eNum
    oNameAssim = argList[6]%eNum
    df, nDate = model.main_day(date, flag="restart", restart="restart.txt", runoffDir=runoffDir, outDir=outDir)
    adf, nDate = model.main_day(date, flag="restart", restart="restartAssim.txt", runoffDir=runoffDir, outDir=outDir)
    model.output(df, oName, mode="a")
    model.output(adf, oNameAssim, mode="a")
    return nDate


def simulation_parallel(date, eTot, model, runoffRoot, outRoot, oName, oNameAssim, ncpus=2):
    args = [[date,eNum,model,runoffRoot,outRoot,oName,oNameAssim] for eNum in range(0,eTot)]
    pool = Pool(processes=ncpus)
    results = pool.map(simulation_core,args)
    pool.close()
    pool.join()
    return results[0]


def updateChannel(argsList):
        aArray1d = argsList[0]
        ndx = argsList[1]
        upas = argsList[2]
        nReach = argsList[3]
        resDf = argsList[4]
        rootDir = argsList[5]
        oDir = argsList[6]
        eNum = argsList[7]
        nupa = argsList[8]
        update_qs = np.repeat(aArray1d, ndx) # [1,2,3,4] => [1,1,1,...,1,2,2,2,...,2,...]
        print(eNum)
        update_qin,\
        update_qout = daTools.updateQinout(aArray1d, upas, nReach, ndx)
        resDf["old_q"] = update_qs
        resDf["old_q_ch_in"] = update_qin
        resDf["old_q_ch_out"] = update_qout
        resOut = os.path.join(rootDir, oDir%eNum,"restartAssim.txt")
        resDf.to_csv(resOut,index=False)


def updateChannel_parallel(xa, ndx, upas, nReach, restarts, rootDir, oDir, eTot, nupa=[3,4,5], ncpus=2):
    args = [[xa[eNum,-1,:], ndx, upas, nReach, restarts[eNum], rootDir, oDir, eNum, nupa] for eNum in range(eTot)]
    pool = Pool(processes=ncpus)
    results = pool.map(updateChannel, args)
    pool.close()
    pool.join()


def assimOut(argsList):
        xa_each = argsList[0]
        eNum = argsList[1]
        date = argsList[2]
        modelInstance = argsList[3]
        cfs2cms = argsList[4]
        oNameAssim = argsList[5]
        df = pd.DataFrame(xa_each).apply(cfs2cms, axis=1)
        df = df.T
        df.index = [date]
        df.reset_index().rename({"index":"Date"}, axis=1).set_index("Date")
        modelInstance.output(df, oNameAssim%eNum, mode="a")


def assimOut_parallel(xa, date, modelInstance, cfs2cms, oNameAssim, eTot, ncpus=2):
    args = [[xa[eNum, -1], eNum, date, modelInstance, cfs2cms, oNameAssim] for eNum in range(eTot)]
    pool = Pool(processes=ncpus)
    results = pool.map(assimOut, args)
    pool.close()
    pool.join()


if __name__ == "__main__":
    test = DA_pyHRR("/project/uma_colin_gleason/yuta/DA/missouri/pyHRR/config_MERIT_cal03.ini")
    test.tester()
    #test.restart()
