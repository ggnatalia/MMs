import numpy as np 
cimport numpy as np 
cimport cython
from cpython.array cimport array
 

# Use Numpy-C-API
#np.import_array()


DTYPE = np.uint8
ctypedef np.uint8_t DTYPE_t

DTYPE2 = np.float32
ctypedef np.float32_t DTYPE2_t


#calculate_distances_cython(multiprocessing_globals_seqs[0], multiprocessing_globals_seqs[1])
@cython.boundscheck(False)
#if boundcheck is true, it will check to see if the index is inside the range of the vector before reading or writing to memory. And if not it will throw an error.
#If boundcheck is false, it will read or write to the pointer even if the index is out of bounds, given out false data by reading and writing to memory corrupting data by writing.
@cython.wraparound(False)
def calculate_distances_cython(DTYPE_t [:] arr1, DTYPE_t [:] arr2):
    cdef Py_ssize_t size1 = arr1.shape[0]
    cdef Py_ssize_t size2 = arr2.shape[0]
    cdef np.float32_t diffs = 0
    cdef Py_ssize_t invalid_pos_1 = 0
    cdef Py_ssize_t invalid_pos_2 = 0
    cdef Py_ssize_t i
    #cdef int size0 
    #assert ( size1 == size2)
    for i in range(size1):
        if not arr1[i] == 0 or arr2[i] == 0: # If i ==0 => i = '.', ignore that position in the distance calc
            if arr1[i] != arr2[i]:
                diffs += 1 # If not are the same nt, sum 1 to the counter diff
                    #print(diffs)
        if arr1[i] == 0 or arr1[i] == 5: #is arr1 it's '.'-0 or gap-5
            invalid_pos_1 +=1
        if arr2[i] == 0 or arr2[i] == 5: #is arr1 it's '.'-0 or gap-5
            invalid_pos_2 +=1
    cdef np.float32_t d
    #print(diffs)
    #print(size1)
    #print(size2)
    #print(invalid_pos_1)
    #print(invalid_pos_2)
    #print(size1-invalid_pos_1)
    #print(size2-invalid_pos_2)
    cdef Py_ssize_t maxl
    if size1-invalid_pos_1 > size2-invalid_pos_2:
        maxl = size1-invalid_pos_1
    else:
        maxl = size2-invalid_pos_2
    d = diffs/maxl
    return d
    








#def dict(arr, i1, i2)

#@cython.boundscheck(False)
#@cython.wraparound(False) 
#cdef inline Py_ssize_t _getIdx(DTYPE_t[:] p, DTYPE_t v) nogil:
#    cdef Py_ssize_t i
#    cdef Py_ssize_t size = p.shape[0]
#    for i in range(size):
#        if p[i] == v:
#            return i
#    else:
#        return PY_SSIZE_T_MAX # use PY_SSIZE_T_MAX to denote v not being inside p

#[ : ]

