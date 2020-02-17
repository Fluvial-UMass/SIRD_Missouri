from __future__ import division
import numpy as np
import scipy.linalg as sp_linalg
cimport numpy as np

DTYPE = np.float64
DTYPE_int = np.int64
ctypedef np.float64_t DTYPE_t
ctypedef np.int64_t DTYPE_t_int

def letkf(np.ndarray[DTYPE_t,ndim=3] allx, np.ndarray[DTYPE_t,ndim=2] observation, np.ndarray[DTYPE_t,ndim=2] obserr, np.ndarray[DTYPE_t,ndim=2] ocean, np.ndarray[DTYPE_t,ndim=2] excGrids, int patch, int eNum, float assimE, float assimW, float assimN, float assimS, float east, float west, float north, float south, float res, float undef):
    
    """
    Data Assimilation with Local Ensemble Transformed Kalman Filter
        inputs:
            allx: numpy.ndarray([nLat,nLon,eNum]): ensemble simulation
            observation: numpy.ndarray([nLat,nLon]): gridded observation with observed or undef values
            ocean: numpy.ndarray([nLat,nLon]): gridded ocean mask
            excGrids: numpy.ndarray([nLat,nLon]): grids to exclude
            patch: int: patch size
            eNum: int: number of ensemble members
            assimE,assimW,assimN,assimS: float: assimilation boundary
            east,west,north,south: float: simulation boundary
            res: float: simulation horizontal resolution in degree
            undef: float: undef value for the observation
            errfix: float: observation error to constract Rinv. Only float value is supported at V-1.0.0.
    """

    # c type declaration
    assert allx.dtype==DTYPE and observation.dtype==DTYPE and obserr.dtype==DTYPE and ocean.dtype==DTYPE and excGrids.dtype==DTYPE

    cdef int patch_side = patch*2+1
    cdef int patch_nums = patch_side**2
    cdef int nLat = allx.shape[1]
    cdef int nLon = allx.shape[2]
    cdef np.ndarray[DTYPE_t,ndim=3] globalxa = np.zeros([eNum,nLat,nLon],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=1] xt = np.ones([patch_nums],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] xf = np.ones([patch_nums,eNum],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=1] xf_mean = np.ones([patch_nums],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] xf_m = np.ones([patch_nums,1],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] xf_me = np.ones([patch_nums,eNum],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] xa = np.ones([patch_nums,eNum],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] local_obs = np.ones([patch_side,patch_side],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=1] local_obs_line = np.ones([patch_nums],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=1] local_obsErr_line = np.ones([patch_nums],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=1] local_ocean_line = np.ones([patch_nums],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=1] local_exc_line = np.ones([patch_nums],dtype=DTYPE)
    cdef int ovs
    cdef int center
    cdef int lon_start = int((assimW+west)/res)
    cdef int lon_end = int((assimE+west)/res)
    cdef int lat_start = int((north-assimN)/res)
    cdef int lat_end = int((north-assimS)/res)
    cdef int lon_cent
    cdef int lat_cent
    cdef int pLon_start 
    cdef int pLon_end
    cdef int pLat_start
    cdef int pLat_end
    cdef int diff
    cdef int i
    cdef np.ndarray[DTYPE_t,ndim=1] dummy = np.zeros([patch_side],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t_int,ndim=1] xRange = np.zeros([patch_side],dtype=DTYPE_int)
    cdef np.ndarray[DTYPE_t_int,ndim=2] xRange_ = np.zeros([patch_side,patch_side],dtype=DTYPE_int)
    cdef np.ndarray[DTYPE_t_int,ndim=1] yRange = np.zeros([patch_side],dtype=DTYPE_int)
    cdef np.ndarray[DTYPE_t_int,ndim=2] yRange_ = np.zeros([patch_side,patch_side],dtype=DTYPE_int)

    cdef np.ndarray[DTYPE_t,ndim=2] H = np.zeros([patch_nums,patch_nums],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] I = np.identity(eNum,dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] R = np.identity(patch_nums,dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] Rinv = np.identity(patch_nums,dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] Ef = np.zeros([patch_nums,eNum],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] Ea = np.zeros([patch_nums,eNum],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] Pa = np.zeros([eNum,eNum],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] Pa_sqr = np.zeros([eNum,eNum],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] d = np.zeros([patch_nums,1],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] Warr = np.zeros([eNum,eNum],dtype=DTYPE)
    # type declaration ends

    # main loop
    for lon_cent in range(lon_start,lon_end+1):
        for lat_cent in range(lat_start,lat_end+1):
            # local patch
            pLon_start = lon_cent-patch
            pLon_end = lon_cent+patch
            xRange = np.arange(pLon_start,pLon_end+1,1)
            pLat_start = lat_cent-patch
            pLat_end = lat_cent+patch
            yRange = np.arange(pLat_start,pLat_end+1,1)
            xRange_, yRange_ = np.meshgrid(xRange,yRange)
            xRange_[np.where(xRange_>lon_end)] = 0 #whatever !circulation look up
            yRange_[np.where(yRange_>lat_end)] = 0 #whatever !circulation look up

            xt = observation[yRange_,xRange_].flatten()
            local_obs = np.ones([patch_side,patch_side],dtype=DTYPE)
            local_obsErr_line = obserr[yRange_,xRange_].flatten()
            local_ocean_line = ocean[yRange_,xRange_].flatten()
            local_exc_line = excGrids[yRange_,xRange_].flatten()
            xf = allx[:,yRange_,xRange_].reshape(-1,eNum)
            xf_mean = xf.mean(axis=1)
            for i in range(0,eNum):
                xf_me[:,i] = xf_mean
            xf_m = xf_mean.reshape(-1,1)
            
            if pLat_start < lat_start:
                diff = lat_start - pLat_start
                local_obs[0:diff] = 0
            if pLat_end > lat_end:
                diff = lat_end - pLat_start
                local_obs[diff+1::] = 0
            if pLon_start < lon_start:
                diff = lon_start - pLon_start
                local_obs[:,0:diff] = 0
            if pLon_end > lon_end:
                diff = lon_end - pLon_start
                local_obs[:,diff+1::] = 0
            local_obs_line = local_obs.flatten()


            local_obs_line[np.where(xt-undef < 1.)] = 0
            local_obs_line[np.where(local_ocean_line == 1.)] = 0
            local_obs_line[np.where(local_exc_line == 0.)] = 0
            ovs = local_obs_line.sum()
            center = int(patch_nums/2)
            if ovs > 0 and local_ocean_line[center] == 0:
                """
                    observation is available in a local patch and the center of the patch is not ocean.
                    LETKF is activated.
                """
                # initialize
                H = np.zeros([patch_nums,patch_nums],dtype=DTYPE)
                #

                H[np.where(local_obs_line == 1.)[0],np.where(local_obs_line == 1.)[0]] = 1
                Ef = xf - xf_m
                Rinv = np.diag((local_obsErr_line**2)**(-1)).astype(DTYPE)
                
                HEft = np.dot(H,Ef).T # (patch_num x eNum)T = eNum x patch_num
                HEf = np.dot(H,Ef) # patch_num x eNum
                HEftRinvHEf = np.dot(np.dot(HEft,Rinv),HEf) # (eNum x patch_num) * (patch_num x patch_num) *(patch_num x eNum) = (eNum x eNum)

                VDVt = I + HEftRinvHEf
                w,v = np.linalg.eigh(VDVt,UPLO="U")
                
                Dinv = np.diag((w+1e-20)**(-1))
                Dsqr = sp_linalg.sqrtm(Dinv)

                Pa = np.dot(np.dot(v,Dinv),v.T)
                Pa_sqr = np.dot(np.dot(v,Dsqr),v.T)
                
                d = (np.dot(H,xt) - np.dot(H,xf_m).T).reshape(-1,1)
                Wvec = np.dot(np.dot(np.dot(Pa,np.dot(H,Ef).T),Rinv),d)
                for i in range(0,eNum):
                    Warr[:,i] = Wvec.reshape(eNum)
                W = Pa_sqr + Warr
                Ea = np.dot(Ef,W)
                
                xa = xf_me + Ea

            else:
                """
                    No observation is available in a local patch or the center of the patch is ocean.
                    No-assimilation. Return prior ensemble mean as a best guess.
                """
                xa = xf_me


            globalxa[:,lat_cent,lon_cent] = xa[center,:]

    return globalxa

