import numpy as np
from cython.parallel import prange
cimport numpy as np
cimport cython
DTYPE_float = np.float32
DTYPE_int = np.int32
ctypedef np.float_t DTYPE_float_t
ctypedef np.int_t DTYPE_int_t


@cython.boundscheck(False)
@cython.wraparound(False)
def updateQinout_pa(const np.float32_t[:] update_qs, const np.int32_t[:, :] upas, int nReach, int ndx):
    cdef list qinAll
    cdef list qoutAll
    cdef list upa_reach
    cdef int reach
    cdef int r
    cdef int i
    cdef int numupas = upas.shape[0]
    cdef float sqrt
    cdef float dq
    cdef float q_init
    cdef float q_tail
    q_ch_in = np.zeros([ndx, nReach], dtype=np.float32)
    cdef np.float32_t [:, :] q_ch_in_view = q_ch_in
    q_ch_out = np.zeros([ndx, nReach], dtype=np.float32)
    cdef np.float32_t [:, :] q_ch_out_view = q_ch_out
    for reach in prange(nReach, nogil=True):
        q_init = 0.
        for i in prange(numupas):
            r = upas[i, reach]
            if r < 0:
                # !shole be less than 0,
                # !as reach 1 is 0 in index.
                continue
            else:
                q_init += update_qs[r]
        q_tail = update_qs[reach]
        sqrt = ((q_tail - q_init)**2)**0.5
        if sqrt < 0.00001:
            for i in prange(ndx):
                q_ch_in_view[i, reach] = q_init
                q_ch_out_view[i, reach] = q_init
        else:
            dq = (q_tail - q_init)/ndx
            for i in prange(ndx):
                q_ch_in_view[i, reach] = q_init + dq*i
                q_ch_out_view[i, reach] = q_init + dq*(i+1)
    return q_ch_in, q_ch_out


def updateQinout(np.ndarray[DTYPE_float_t, ndim=1] update_qs, np.ndarray[DTYPE_int_t, ndim=2] upas, int nReach, int ndx):
    cdef list qinAll
    cdef list qoutAll
    cdef list upa_reach
    cdef int reach
    cdef float dq
    cdef float q_init
    cdef float q_tail
    cdef np.ndarray[DTYPE_float_t, ndim=1] q_ch_in
    cdef np.ndarray[DTYPE_float_t, ndim=1] q_ch_out
    qinAll = []
    qoutAll = []
    for reach in range(0, nReach):
        upa_reach = upas[:, reach].tolist()
        upa_reach = [r for r in upa_reach if r >= 0]
        q_init = update_qs[upa_reach].sum() # Fancy indexing
        q_tail = update_qs[reach]
        if q_tail < 0.00001:
            q_ch_in = np.zeros(ndx)
            q_ch_out = np.zeros(ndx)
        if abs(q_tail - q_init) < 0.00001:
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
