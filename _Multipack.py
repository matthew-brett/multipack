"""Multipack module

   This module defines some functions that are higher-level interfaces to
   some low level wrappers in the _multipackmodule.c file.  All of these
   functions take a python function as the first argument.

   Contained in this module.

   fsolve    -  find the roots of N nonlinear equations in N unknowns
   leastsq   -  minimize the sum of squares of M(>N) equations in N unknowns.
   odeint    -  integrate a set of ODE's
   quad      -  integrate a function of one variable.
"""

