import os
import sys

import numpy as np
import pandas as pd
import datetime
import configparser
import subprocess
from multiprocessing import Pool

import pyletkf.src.pyletkf as pyletkf
import pyHRR
import daTools

class DA_pyHRR(object):
    def __init__(self,config):

        self.config = config
        self.initFile = configparser.ConfigParser()
        self.initFile.read(config)

        self.undef = self.initFile.get("observation","undef")
        self.rootDir = self.initFile.get("model","rootDir")
        self.runoffDir = "data/case3/%02d"
        self.oDir = "out/case3/%02d"
        self.oFile = "discharge.csv"
        self.oFileAssim = "dischargeAssim.csv"
        self.oFileAssimSmooth = "dischargeAssim_smooth.csv"
        self.oName = os.path.join(self.oDir,self.oFile)
        self.oNameAssim = os.path.join(self.oDir,self.oFileAssim)
        self.oNameAssimSmooth = os.path.join(self.oDir,self.oFileAssimSmooth)
        if not os.path.exists(self.oDir):
            os.makedirs(self.oDir)

        self.compile_ = False

        self.eTot = int(self.initFile.get("assimilation","ensMem"))
        self.smoothMax = 365*2 #depends on your memory size and your observation frequency

        # modeling variables
        self.nReach = int(self.initFile.get("model","nReach"))
        self.reachStart = 1
        self.ndx = 20
        self.channelPath =os.path.join(self.rootDir,"src","channels.txt")
        self.channelInfo = pd.read_csv(self.channelPath,sep="\s+",header=None)
        self.daCore = pyletkf.LETKF_core(config,mode="vector",use_cache=False)
        self.spinupPath = os.path.join(self.rootDir,self.oDir,"00/restart.txt")
        #initialize pyletkf
        self.daCore.initialize()
        #initialize pyHRR
        self.model = pyHRR.HRR(config, compile_=self.compile_)
        #read obs
        self.obsDf = self.readObs()
        self.assimDates = self.obsDf.index
        #load upstream reach info.
        nupa = [3,4,5] # number of columns that contains upstream reach info.
        upas = []
        for i in nupa:
            upa = self.channelInfo.iloc[:,i].astype(np.int64).values.reshape(1,-1)
            upas.append(upa-1) # reach number starts from 1
        self.upas = np.concatenate(upas,axis=1)

    def spinup(self, sDate, eDate, model, eNum=0):
        runoffDir = os.path.join(self.rootDir,self.runoffDir%eNum)
        storageDir = os.path.join(self.rootDir,self.oDir,"spinup/%s"%sDate.strftime("%Y%m%d%H"))
        if os.path.exists(storageDir) == False:
            os.makedirs(storageDir)
        out, nDate = model.main_day(sDate, flag="initial", restart="restart.txt", runoffDir=runoffDir, outDir=storageDir)
        date = nDate
        while date < eDate:
            out, nDate = model.main_day(date, flag="restart", restart="restart.txt", runoffDir=runoffDir, outDir=storageDir)
            date = nDate

        return os.path.join(storageDir, "restart.txt")

    def mainRun(self, sDate, eDate, spinup=False, sDate_spin=None, eDate_spin=None, ncpus=2):
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
            nDate = simulation_parallel(date, self.eTot, self.model, rnofRt, outRt, self.oName, self.oNameAssim, ncpus=ncpus)
            if (date == self.assimDates).any():
                xa = self.dataAssim(date, self.obsDf, self.eTot)
                [self.__assimOut(xa[eNum],eNum,date,self.model) for eNum in range(0, self.eTot)]
        
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
            model.output(df,self.oName,mode="w")
            model.output(df,self.oNameAssim,mode="w")
        return nDate

    def dataAssim(self, date, obsDf, eTot):
        obsMean, obsStd = self.__constObs(obsDf,date)
        self.restarts = []
        dschg_cfs_ens = np.concatenate([self.__readRestart(eNum) for eNum in range(0,eTot)], axis=0)
        dschg_cfs_ens = dschg_cfs_ens.astype(np.float64).reshape(eTot,1,self.nReach) 
        # As an API requirements, thd input simulation array should be (eTot, nT, nReach)

        xa = self.daCore.letkf_vector(dschg_cfs_ens,obsMean,obsStd.astype(np.float64))
        xa[np.where(xa<0.)] = np.absolute(xa[np.where(xa<0.)])
        [self.__updateChannel(xa[eNum],self.restarts[eNum],eNum) for eNum in range(0,eTot)]
        return xa

    def __assimOut(self, xa_each, eNum, date, modelInstance):
        df = pd.DataFrame(xa_each).apply(self.cfs2cms, axis=1)
        df.index = [date]
        df.reset_index().rename({"index":"Date"}, axis=1).set_index("Date")
        modelInstance.output(df, self.oNameAssim, mode="a")

    def __readRestart(self,eNum):
        dfPath = os.path.join(self.rootDir,self.oDir%eNum,"restart.txt")
        edf = pd.read_csv(dfPath)
        self.restarts.append(edf)
        old_qs = (edf.groupby("i").max())["old_q"].values #cfs
        return old_qs.reshape(1,-1)

    def __updateChannel(self, aArray1d, resDf, eNum, nupa=[3,4,5]):
        update_qs = np.repeat(aArray1d,self.ndx) # [1,2,3,4] => [1,1,1,...,1,2,2,2,...,2,...]
        update_qin,\
        update_qout = daTools.updateQinout(update_qs,self.upas,self.nReach,self.ndx)
        resDf["old_q"] = update_qs
        resDf["old_q_ch_in"] = update_qin
        resDf["old_q_ch_out"] = update_qout
        resOut = os.path.join(self.rootDir,self.oDir%eNum,"restartAssim.txt")
        resDf.to_csv(resOut,index=False)
    
    def __constObs(self,obsDf,date):

        obs = obsDf[obsDf.index==date]

        reaches = obs.Reach.values - 1 #for python which starts from 0
        obsMean = np.ones([self.nReach]) * self.daCore.undef
        obsStd  = np.ones([self.nReach]) * 0.01

        obsMean[reaches] = self.cms2cfs(obs.Mean.values)
        obsStd[reaches] = self.cms2cfs(obs.Std.values)

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

    def tester(self):
        sDate_spin = datetime.datetime(2002,9,25)
        eDate_spin = datetime.datetime(2002,10,1)
        sDate = eDate_spin
        eDate = datetime.datetime(2002,10,5)
        self.mainRun(sDate, eDate, spinup=True, sDate_spin=sDate_spin, eDate_spin=eDate_spin)

# For multi-processing
def simulation_core(argList):
    date = argList[0]
    eNum = argList[1]
    model = argList[2]
    runoffDir = os.path.join(argList[3]%eNum)
    outDir = os.path.join(argList[4]%eNum)
    oName = argList[5]
    oNameAssim = argList[6]
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

if __name__ == "__main__":
    test = DA_pyHRR("config.ini")
    test.tester()

