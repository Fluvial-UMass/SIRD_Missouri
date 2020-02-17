import pandas as pd
import datetime
import numpy as np
import gc

sdate = datetime.datetime(2002,10,10)
edate = datetime.datetime(2011,1,1)

# dfiles = ["/home/yi79a/DA/missouri/pyHRR/out/MERIT_VIC_beck_cal03/%02d/discharge.csv", "/home/yi79a/DA/missouri/pyHRR/out/MERIT_VIC_beck_cal03/%02d/dischargeAssim.csv"]
dfiles = ["/home/yi79a/DA/missouri/pyHRR/out/MERIT_VIC_beck_cal03_syears/%02d/dischargeAssim.csv"]
# outs = ["./meritvic_beck_open.csv", "./meritvic_beck_assim.csv"]
# outs = ["./dischargeInsert_uncal.csv"]
outs = ["meritvic_syears_assim_cal03.csv"]

validation_reaches = pd.read_csv("./validation_reaches.csv")["reach"].values
nReach = 28998
eNum = 20

product = []
for idx,fname in enumerate(dfiles):
    for ens in range(0,eNum):
        print(ens)
        out = []
        with open(fname%ens, mode="r") as f:
            lines = f.readlines()
        for l in lines:
            elements = l.strip("\n").split(",")
            if len(elements) == 404:
                out.append(elements)
            else:
                print(len(elements))
                date = elements[0]
                part = np.array(elements[1::])[(validation_reaches - 1)]
                print(date)
                part = np.insert(part, 0, date)
                print(len(part))
                out.append(part)
        df = pd.DataFrame(out[1::], columns=out[0]).set_index("Date").astype(float)
        df.index = pd.to_datetime(df.index)
        df = df[sdate:edate]
        print(df)
        if ens == 0:
            DF = df
        else:
            DF = df + DF
        del df
    enDf = DF/eNum
    print(enDf)
    enDf = enDf[sdate:edate]
    enDf.to_csv(outs[idx])
