import numpy as np
cimport numpy as np
DTYPE_float = np.float64
DTYPE_int = np.int64
ctypedef np.float_t DTYPE_float_t
ctypedef np.int_t DTYPE_int_t

def updateQinout(np.ndarray[DTYPE_float_t, ndim=1] update_qs, np.ndarray[DTYPE_int_t, ndim=2] upas, int nReach, int ndx):
    cdef list qinAll
    cdef list qoutAll
    cdef list upa_reach
    cdef int reach
    cdef float q_init
    cdef float q_tail
    cdef np.ndarray[DTYPE_int_t, ndim=1] upa_reach_array
    cdef np.ndarray[DTYPE_float_t, ndim=1] q_ch_in
    cdef np.ndarray[DTYPE_float_t, ndim=1] q_ch_out
    qinAll = []
    qoutAll = []
    for reach in range(0, nReach):
        upa_reach_array = upas[:, reach]
        if (upa_reach_array == -1).any():
            upa_reach = upa_reach_array.tolist()
            upa_reach.remove(-1)
        else:
            upa_reach = upa_reach_array.tolist()
        q_init = update_qs[upa_reach].sum() # Fancy indexing
        q_tail = update_qs[reach]
        if q_tail < 0.00001:
            q_ch_in = np.zeros(ndx)
            q_ch_out = np.zeors(ndx)
        if (q_tail - q_init) < 0.00001:
            q_ch_in = np.ones(ndx)*q_init
            q_ch_out = q_ch_in
        else:
            dq = (q_tail - q_init)/ndx
            q_ch_in = np.arange(q_init,q_tail,dq)[0:ndx]
            q_ch_out = np.arange(q_init+dq,q_tail+dq,dq)[0:ndx]
        qinAll.append(q_ch_in)
        qoutAll.append(q_ch_out)
    q_ch_in = np.concatenate(qinAll)
    q_ch_out = np.concatenate(qoutAll)
    return q_ch_in, q_ch_out
