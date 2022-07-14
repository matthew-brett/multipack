
from Multipack import *
from Numeric import *
from RandomArray import random

def fun1(x,p):
    y = zeros(len(x),'d')
    y[0] = x[0]*x[0] + x[1]*x[1] - p
    y[1] = x[0] - x[1]
    return y

def Dfun1(x,p):
    dy = zeros((2,2),'d')
    dy[0] = [2*x[0], 2*x[1]]
    dy[1] = [1.0, -1.0]
    return dy

def quad_fun(x,w):
    return sin(w*x)

def quad_fun2(x,w):
    return exp(w*x)

def ode_fun(y,t):
    yp = [0]*3
    yp[0] = -0.04*y[0] + 1e4*y[1]*y[2]
    yp[1] = 0.04*y[0] - 1e4*y[1]*y[2] - 3e7*y[1]**2
    yp[2] = 3e7*y[1]**2
    return yp

# for full-bandwidth case, needs to be given so that df[i]/dy[j] is in
#     element [i,j]    (function down row, derivative across columns) unless
#     col_deriv=1 is given as argument to odeint (faster since no internal
#     transpose.)
def Dode_fun(y,t):  
    dyp = zeros((3,3))
    dyp[0] = [-0.04, 1e4*y[2], 1e4*y[1]]
    dyp[1] = [0.04, -1e4*y[2]-6e7*y[1], -1e4*y[1]]
    dyp[2] = [0, 6e7*y[1], 0]
    return dyp;

# for banded case, need the diagonals of Dfun placed in the columns of
# the return matrix  (so return matrix is num_equations * (ml + mu + 1))
# start with the lowest diagonal.
# If ml or mu is given but the other one is not the one not given is assumed
#  to be zero.  If neither is given the matrix is assumed full.


def thefunc(a,data,t):
    return exp(a[0]*t) + exp(a[1]*t) - data


if __name__ == '__main__':
    x0 = [0.3,0.3];
    for p in range(1,10):
        g = fsolve(fun1,x0,args=(p,))
        print "p = %d, sqrt(p/2.0) = %f, g[0] = %f, g[1] = %f" % (p,sqrt(p/2.0),g[0],g[1])

    x0 = [1.0,1.0];
    for p in range(1,10):
        g = fsolve(fun1,x0,args=(p,),Dfun=Dfun1)
        print "p = %d, sqrt(p/2.0) = %f, g[0] = %f, g[1] = %f" % (p,sqrt(p/2.0),g[0],g[1])

    y = quad(quad_fun,0,3,(2,))
    print "Error is " + `y[0]-1/2.0*(1-cos(2*3))`

    y = quad(quad_fun2,-Inf,0,(2,))
    print y[0]
    
    an = odeint(ode_fun,[1,0,0],[0,0.4,4.0,40])
    print an
    an = odeint(ode_fun,[1,0,0],[0,0.4,4.0,40],Dfun=Dode_fun)
    print an
    an = odeint(ode_fun,[1,0,0],[0,0.4,4.0,40],Dfun=Dode_fun,full_output=1)
    print an
    
    # To test leastsq we will fit some data to a sum of exponentials.
    # First the data:
    a1 = -3.0; a2 = -2.0;
    noise_width = 0.05
    t = arange(0,4,0.1)
    data = exp(a1*t) + exp(a2*t) + noise_width*random(t.shape)

    ret = leastsq(thefunc,[-1,-1],args=(data,t))
    print ret[-1]
    print ret[0]

    

    

    

