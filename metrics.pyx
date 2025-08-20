#cython: initializedcheck=False, boundscheck=False, wraparound=False, cdivision=True, profile=False
cimport cython

# cimport the Cython declarations for numpy
cimport numpy as np
from libc.string cimport memcpy
import numpy as np
import time

# if you want to use the Numpy-C-API from Cython
# (not strictly necessary for this example, but good practice)
np.import_array()

#cdef extern from "correntropy.h":
#    ctypedef double dtype_t
#    dtype_t crossCorrentropy(dtype_t*, dtype_t*, dtype_t*, int, int, dtype_t, int, int)

cdef extern from "metrics.h":
    void cont(double*, double*, int, double*, long)
    void trust(double*, double*, int, double*, long)
    void qm(double*, double*, int, int, double*, long)

ctypedef np.float_t DTYPE_t
DTYPE = np.float


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function 
def cont_c(np.ndarray[DTYPE_t, ndim=2, mode="c"] x not None,  np.ndarray[DTYPE_t, ndim=2, mode="c"] y not None, int k):
    cdef unsigned long N = y.shape[1]
    cdef np.ndarray[DTYPE_t, ndim=1] V = np.zeros((1), dtype=DTYPE)
    cont(<double*> np.PyArray_DATA(x), <double*> np.PyArray_DATA(y), k, <double*> np.PyArray_DATA(V), N)
    return V

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function 
def trust_c(np.ndarray[DTYPE_t, ndim=2, mode="c"] x not None,  np.ndarray[DTYPE_t, ndim=2, mode="c"] y not None, int k):
    cdef unsigned long N = y.shape[1]
    cdef np.ndarray[DTYPE_t, ndim=1] V = np.zeros((1), dtype=DTYPE)
    trust(<double*> np.PyArray_DATA(x), <double*> np.PyArray_DATA(y), k, <double*> np.PyArray_DATA(V), N)
    return V

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function 
def qm_c(np.ndarray[DTYPE_t, ndim=2, mode="c"] x not None,  np.ndarray[DTYPE_t, ndim=2, mode="c"] y not None, int n, int k):
    cdef unsigned long N = y.shape[1]
    cdef np.ndarray[DTYPE_t, ndim=1] V = np.zeros((1), dtype=DTYPE)
    qm(<double*> np.PyArray_DATA(x), <double*> np.PyArray_DATA(y), n, k, <double*> np.PyArray_DATA(V), N)
    return V

