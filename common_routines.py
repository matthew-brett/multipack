from Numeric import *

def myasarray(a):
    if type(a) in [type(1.0),type(1L),type(1),type(1j)]:
        return asarray([a])
    else:
        return asarray(a)

def check_func(thefunc, x0, args, numinputs, output_shape=None):
    args = (x0[:numinputs],) + args
    res = myasarray(apply(thefunc,args))
    if (output_shape != None) and (res.shape != output_shape):
        if (output_shape[0] != 1) and (output_shape[1] != 1):
           raise TypeError, "There is a mismatch between the input and output shape of %s." % thefunc.func_name
    return res.shape



