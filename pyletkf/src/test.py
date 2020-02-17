import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
import pyletkf

nT = 25
noises = np.random.normal(0,0.5,[20,58])

def forward(x):
    xNext = np.abs(x*np.cos(5*x))
    return xNext

obs = np.ones([nT,58])*10 + np.abs(np.random.normal(10,5,[nT,58]))
for t in range(1,nT):
    obs[t] = forward(obs[t-1]) + np.random.normal(0,0.5)
sim = np.ones([20,nT,58])*10 + np.abs(np.random.normal(10,5,[20,nT,58])) + 10
for e in range(0,20):
    for t in range(1,nT):
        sim[e,t,:] = forward(sim[e,t-1,:] + sim[e,t-1,:]*noises[e]) + 10

t = 24
reach = 25
obs_t = obs[t]
obs_e = np.ones([58])*1
sim_t = sim[:,0:t,:]
am = pyletkf.LETKF_core("./config.ini",mode="vector",use_cache=False)
am.initialize()
data1 = am.letkf_vector(sim_t.astype(np.float64),obs_t.astype(np.float64),obs_e.astype(np.float64),smoother=True)
data2 = am.letkf_vector(sim_t.astype(np.float64),obs_t.astype(np.float64),obs_e.astype(np.float64),smoother=False)
assim1 = sim.copy()
assim2 = sim.copy()
assim1[:,0:t,:] = data1
assim2[:,0:t,:] = data2
plt.figure()
plt.plot(obs[:,reach],color="k",linewidth=2,linestyle="--")
for e in range(0,20):
    plt.plot(sim[e,:,reach],color="grey",alpha=0.5)
for e in range(0,20):
    plt.plot(assim1[e,:,reach],alpha=0.5)
plt.plot(sim[:,:,reach].mean(axis=0),color="b",marker="^",linewidth=2,label="sim")
plt.plot(assim1[:,:,reach].mean(axis=0),color="r",linewidth=2,label="assim_smoothed")
plt.plot(assim2[:,:,reach].mean(axis=0),color="orange",marker="o",linewidth=2,label="assim_original")
plt.xlabel("Time step")
plt.ylabel("Forecasted value")
plt.ylim(0,250)
plt.title("The effect of Ensemble Kalman Smoother")
plt.legend()
plt.show()
