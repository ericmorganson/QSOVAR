import numpy as np
from numpy.ctypeslib import load_library
from numpyctypes import c_ndarray

mylib = load_library('libmyextension', 'lib')       # '.' is the directory of the C++ lib  

def ln_like_fast(array1, array2):
    arg1 = c_ndarray(array1, dtype=np.double, ndim = 3, shape = (4,3,2))
    arg2 = c_ndarray(array2, dtype=np.double, ndim = 3, shape = (4,3,2))
    return mylib.myfunc(arg1, arg2)


