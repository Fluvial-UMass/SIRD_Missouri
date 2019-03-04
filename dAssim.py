import numpy as np
import pandas as pd
import datetime
import configparser
import subprocess
import os
import sys
import pyletkf.src.pyletkf as pyletkf
import pyHRR


class DataAssim(object):

    def __init__(self,config):
        
        self.config = config
        self.initFile = configparser.ConfigParser()
        self.initFile.read(config)
        self.sDate = datetime.datetime(2000,1,1,0)
        self.eDate = datetime.datetime(2011,1,1,0)

        self.undef = float(self.initFile.get("observation","undef"))
        self.runoffDir = "../data/case3/%02d"
        self.oDir = "./out/case3/"
        self.oFile = "discharge_%02d.csv"
        self.oFileAssim = "dischargeAssim_%02d.csv"
        self.oName = os.path.join(self.oDir,self.oFile)
        self.oNameAssim = os.path.join(self.oDir,self.oFileAssim)
        if not os.path.exists(self.oDir):
            os.makedirs(self.oDir)

        self.compile_ = False

        self.eNum = int(self.initFile.get("assimilation","ensMem"))

        self.nReach = 58
        self.reachStart = 1
        self.ndx = 20
        self.channelInfo = pd.read_csv("./src/channels.txt",sep="\s+",header=None)
        self.daCore = pyletkf.LETKF_core("config.ini",mode="vector",use_cache=False)
        #initialize pyletkf
        self.daCore.initialize()

    
    def main(self):

        obsDf = self.readObs()
        assimDate = obsDf.index

        # initial HRR run
        model = pyHRR.HRR(self.config,compile_=self.compile_)
        for ens in range(0,self.eNum):
            subprocess.call(["cp","./spinup/20000101/restart_%02d.txt"%ens,"./src/restart_%02d.txt"%ens])
            df,nDate = model.main_day(self.sDate,flag="restart",restart="restart_%02d.txt"%ens,runoffDir=self.runoffDir%ens)
            subprocess.call(["cp","./src/restart_%02d.txt"%ens,"./src/restartAssim_%02d.txt"%ens])
            model.output(df,self.oName%ens,mode="w")
            model.output(df,self.oNameAssim%ens,mode="w")
        date = nDate
        #date = self.sDate
        flag = True
        while date < self.eDate:
            # ensemble loop
            for ens in range(0,self.eNum):
                # open loop
                df,nDate = model.main_day(date,flag="restart",restart="restart_%02d.txt"%ens,runoffDir=self.runoffDir%ens)
                model.output(df,self.oName%ens,mode="a")
                # assim loop
                dfAssim,nDate = model.main_day(date,flag="restart",restart="restartAssim_%02d.txt"%ens,runoffDir=self.runoffDir%ens)
                if not (date == assimDate).any():
                    # no assimilation. just update a state.
                    model.output(dfAssim,self.oNameAssim%ens,mode="a")
            # assimilation if obs. is available
            if (date == assimDate).any():
                xa = self.dataAssim(date,obsDf)
                [self.__assimOut(xa[ens],ens,date,model) for ens in range(0,self.eNum)]
                
            date = nDate

    
    def __assimOut(self,xa_each,ens,date,model):

        data = [date]
        [data.append(self.__cfs2cms(xa_each[reach-self.reachStart])) for reach in range(1,self.nReach+1)] # unit conversion cfs cms
        df = pd.DataFrame([np.array(data).T]).rename(columns={0:"Date"})
        df = df.set_index("Date")
        model.output(df,self.oNameAssim%ens,mode="a")


    def __cfs2cms(self,cfs):

        cms = cfs*0.3048**3

        return cms


    def __cms2cfs(self,cms):

        cfs = cms/(0.3048**3)
        cfs_undef = self.undef/(0.3048**3)
        cfs[np.where(cfs-cfs_undef < 0.0001)] = self.undef

        return cfs


    def readObs(self):

        obsPath = self.initFile.get("observation","obsPath")
        df = pd.read_csv(obsPath,parse_dates=[0],index_col=[0])

        return df


    def dataAssim(self,date,obsDf):

        obsMean, obsStd = self.__constObs(obsDf,date)
        sim = np.zeros([self.eNum,self.nReach])
        RES = []
        for ens in range(0,self.eNum):
            e = pd.read_csv("./src/restartAssim_%02d.txt"%ens)
            RES.append(e)
            for reach in range(1,self.nReach+1):
                qout = e[e.i == reach]["old_q"].values[0]
                sim[ens,reach-self.reachStart] = qout
        xa = self.daCore.letkf_vector(sim.astype(np.float64),obsMean,obsStd.astype(np.float64))
        # never allow 0
        xa[np.where(xa<0.)] = np.absolute(xa[np.where(xa<0.)])

        [self.__updateChannel(xa,RES[ens],ens) for ens in range(0,self.eNum)]

        return xa


    def __constObs(self,obsDf,date):

        obs = obsDf[obsDf.index==date]

        reaches = obs.Reach.values - 1 #for python which starts from 0
        obsMean = np.ones([self.nReach]) * self.daCore.undef
        obsStd  = np.ones([self.nReach]) * 0.01

        obsMean[reaches] = self.__cms2cfs(obs.Mean.values)
        obsStd[reaches] = self.__cms2cfs(obs.Std.values)

        return obsMean, obsStd


    def __updateChannel(self,xa,e,ens):

        for reach in range(1,self.nReach+1):
            # update old_q,self.nReach
            e.loc[e.i == reach, "old_q"] = xa[ens][reach-self.reachStart]
            n = pd.read_csv("./src/restart_%02d.txt"%ens)
        for reach in range(1,self.nReach+1):
            # update old_q_ch
            UPSTRM_raw = self.channelInfo.iloc[reach-1,3:6].values.astype(np.int64).tolist()
            UPSTRM = []
            [UPSTRM.append(i) for i in UPSTRM_raw if not i == -1]
            q_ch_int = 0
            for upstrm in UPSTRM:
                old_q = e.loc[e.i == upstrm, "old_q"].values[0]
                q_ch_int = q_ch_int + old_q
            q_ch_last = e.loc[e.i == reach, "old_q"].values[0]
            if q_ch_last < 0.00001:
                q_ch_in = np.zeros([self.ndx])
                q_ch_out = np.zeros([self.ndx])
            elif q_ch_int - q_ch_last < 0.00001:
                q_ch_in = np.ones([self.ndx]) * q_ch_int
                q_ch_out = np.ones([self.ndx]) * q_ch_int
            else:
                dq = (q_ch_last - q_ch_int)/self.ndx
                q_ch_in = np.arange(q_ch_int,q_ch_last,dq)[0:self.ndx]
                q_ch_out = np.arange(q_ch_int+dq,q_ch_last+dq,dq)[0:self.ndx]
            e.loc[e.i == reach, "old_q_ch_in"] = q_ch_in
            e.loc[e.i == reach, "old_q_ch_out"] = q_ch_out
        e.to_csv("./src/restartAssim_%02d.txt"%ens,index=False)    


if __name__ == "__main__":
    da = DataAssim("./config.ini")
    da.main()
