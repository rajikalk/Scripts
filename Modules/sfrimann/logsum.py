#!/usr/bin/env python

# so that python 2.x and 3.x will work the same
from __future__ import print_function, absolute_import, division

import numpy as np

def logadd(logx, logy):
    """
    Return log(x+y), avoiding arithmetic underflow/overflow.
    logx: log(x)
    logy: log(y)
    Rationale:
    x + y    = e^logx + e^logy
             = e^logx (1 + e^(logy-logx))
    log(x+y) = logx + log(1 + e^(logy-logx)) (1)
    Likewise, 
    log(x+y) = logy + log(1 + e^(logx-logy)) (2)
    The computation of the exponential overflows earlier and is less precise
    for big values than for small values. Due to the presence of logy-logx
    (resp. logx-logy), (1) is preferred when logx > logy and (2) is preferred
    otherwise.
    """

    if (logx == np.inf) & (logy == np.inf):
        return np.inf
    if (logx == -np.inf) & (logy == -np.inf):
        return -np.inf

    if logx > logy:
        return logx + np.log(1 + np.exp(logy-logx))
    else:
        return logy + np.log(1 + np.exp(logx-logy))

def log10add(logx, logy):
    """
    Return log(x+y), avoiding arithmetic underflow/overflow.
    logx: log(x)
    logy: log(y)
    Rationale:
    x + y    = e^logx + e^logy
             = e^logx (1 + e^(logy-logx))
    log(x+y) = logx + log(1 + e^(logy-logx)) (1)
    Likewise, 
    log(x+y) = logy + log(1 + e^(logx-logy)) (2)
    The computation of the exponential overflows earlier and is less precise
    for big values than for small values. Due to the presence of logy-logx
    (resp. logx-logy), (1) is preferred when logx > logy and (2) is preferred
    otherwise.
    """
    if (logx == np.inf) & (logy == np.inf):
        return np.inf
    if (logx == -np.inf) & (logy == -np.inf):
        return -np.inf
    
    if logx > logy:
        return logx + np.log10(1 + 10**(logy-logx))
    else:
        return logy + np.log10(1 + 10**(logx-logy))

"""
logsum_ufunc is a numpy ufunc (universal function) and as a result contains
reduce, accumulate, reduceat and outer.
"""
logsum_ufunc = np.frompyfunc(logadd, 2, 1)
log10sum_ufunc = np.frompyfunc(log10add, 2, 1)

"""
logsum(loga, axis=0, dtype=None, out=None)
Take the log of the sum of array elements over a given axis.
For example, for an array a=[a_1,...,a_N], it returns \log \sum_n a_n.
loga: numpy.log(a)
axis: The axis along which to apply the log sum.
dtype: The type used to represent the intermediate results.
out: A location into which the result is stored.
Examples:
logsum(1darray) => scalar
logsum(2darray, axis=0) => 1darray
"""
logsum = logsum_ufunc.reduce
log10sum = log10sum_ufunc.reduce

# Unit-tests...
if __name__ == "__main__":
    import unittest

    class Test(unittest.TestCase):
        def test_1d(self):
            inp = np.arange(1,10)
            out = logsum(np.log(inp))
            expected = np.log(np.sum(inp))
            self.assertAlmostEquals(out, expected)


        def test_2d(self):
            for a in (0,1):
                inp = np.arange(1,10).reshape(3,3)
                out = logsum(np.log(inp), axis=a)
                expected = np.log(np.sum(inp, axis=a))

                #self.assertTrue(np.allclose(out, expected))
                for i in range(3):
                    self.assertAlmostEquals(out[i], expected[i])

    unittest.main()