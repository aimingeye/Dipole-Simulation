# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 23:30:04 2021

@author: Mihir Gadgi
"""
import numpy as np
import scipy.integrate
import numdifftools as nd
d = 1
q_0 = 1
w = np.pi
e = 1
c = 1
epsilon_0 = 1
mu_0 = 1
def q(t,d,w,e,q_0):
    return q_0*np.cos(w*t)
def i(t,d,w,e,q_0):
    return -2*q_0*w*np.sin(w*t)*np.array([0,0,1], float)
def A(r,t):
    k = scipy.integrate.quad_vec(
        lambda r_p: i(t-np.linalg.norm(r-r_p)/c,d,w,e,q_0)/np.linalg.norm(r-r_p),
        -d/2,
        d/2
    )
    return mu_0/(4*np.pi)*k[0]
def phi(r,t):
    r_p = d/2*np.array([[0,0,1],[0,0,-1]],float)
    return sum([((-1)**i)*(1/(4*np.pi*epsilon_0))*
                q(t-np.linalg.norm(r-r_p[i])/c,d,w,e,q_0)/np.linalg.norm(r-r_p[i]) for i in range(2)])
def E(x, t):
    return -nd.Gradient(lambda x: phi(x,t))(x) - nd.Gradient(lambda t: A(x,t))(t)
def curl(f,x):
    jac = nd.Jacobian(f)(x)
    return np.array([jac[2,1]-jac[1,2],jac[0,2]-jac[2,0],jac[1,0]-jac[0,1]])
def B(x,t):
    return curl(lambda x: A(x,t), x)