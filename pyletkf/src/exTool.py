#!/usr/bin/env python
# -*- conding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import datetime
import configparser
import pickle
from numba import jit
from . import letkf
#import letkf

@jit
def constLocalPatch_vector(networkFile,patchArea,nReach,localPatchPath,reach_start=1):

    PATCHES = []
    if reach_start > 1:
        sys.stderr.write("RuntimeWarning: reach_start > 1 might cause serious error. It should be set to 0 (C-style) or 1 (Fortran-style).")
    print("reading river network...")
    network = pd.read_csv(networkFile)
    print("constructing local patches...")
    print("Maximum Local patch size: %s" %"{:.3f}".format(patchArea))
    for reach in range(reach_start,nReach+reach_start):
        PATCH = [reach]
        area = 0 # initialize
        uArea = network.loc[network.reach == reach]["upArea"].values[0]
        area = uArea
        flag = "up"
        upperBoundary = False
        bottomBoundary = False
        #while area < patchArea:
        for i in range(0,100):
            if upperBoundary == True and bottomBoundary == True:
                # no up/downstream catchments. To avoid infinite loop.
                break
            elif flag == "up":
                upReaches,area = concatUparea(network,reach,area,patchArea)
                flag = "down"
                for r in upReaches:
                    PATCH.append(r)
            elif flag == "down":
                downReaches,area = concatDownarea(network,reach,area,patchArea)
                flag = "up_iter"
                for r in downReaches:
                    PATCH.append(r)
            elif flag == "up_iter":
                UP = []
                for upReach in upReaches:
                    upReaches_each, area = concatUparea(network,upReach,area,patchArea)
                    for r in upReaches_each:
                        UP.append(r)
                    if area > patchArea:
                        break
                upReaches = UP
                for r in upReaches:
                    PATCH.append(r)
                flag = "down_iter"
                if len(upReaches) == 0:
                    # no upstream catchments.
                    upperBoundary = True
            elif flag == "down_iter":
                DOWN = []
                for downReach in downReaches:
                    downReaches_each, area = concatDownarea(network,downReach,area,patchArea)
                    for r in downReaches_each:
                        PATCH.append(r)
                    if area > patchArea:
                        break
                downReaches = DOWN
                for r in downReaches:
                    PATCH.append(r)
                flag = "up_iter"
                if len(downReaches) == 0:
                    # no downstream catchments.
                    bottomBoundary = True
            else:
                raise IOError("iteration goes wrong conditional branch. Fatal Error, Abort.")
        PATCH = (np.array(PATCH) - reach_start).tolist() #pythonize
        print("reach: %d"%reach,"\n\t",PATCH)
        PATCHES.append(PATCH)

    return PATCHES
 

@jit
def concatUparea(network,reach,area,patchArea,max_upReachNum=5):
    
    REACHES = []
    for upNum in range(0,max_upReachNum):
        columnName = "u%d"%upNum
        upReach = int(network.loc[network.reach == reach][columnName])
        if upReach == -1:
            continue
        uArea = network.loc[network.reach == upReach]["upArea"].values[0]
        area = area + uArea
        if area > patchArea:
            break
        REACHES.append(upReach)
    
    return REACHES, area


@jit
def concatDownarea(network,reach,area,patchArea,max_downReachNum=5):

    REACHES = []
    for downNum in range(0,max_downReachNum):
        columnName = "d%d"%downNum
        downReach = int(network.loc[network.reach == reach][columnName])
        if downReach == -1:
            continue
        uArea = network.loc[network.reach == downReach]["upArea"].values[0]
        area = area + uArea
        if area > patchArea:
            break
        REACHES.append(downReach)

    return REACHES, area
