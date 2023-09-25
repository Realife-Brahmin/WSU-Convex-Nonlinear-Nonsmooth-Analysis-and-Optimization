#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 17 10:32:24 2023

@author: asaki
"""

def TestFunction1(x,p,nargout):

    import numpy as np

    a=4-2.1*x[0]**2+(1/3)*x[0]**4
    c=4*x[1]**2-4
    f=a*x[0]**2+x[0]*x[1]+c*x[1]**2+p

    if nargout>1:
        g=np.zeros((2,1))
        g[0]=2*x[0]*a+x[0]**2*((4/3)*x[0]**3-4.2*x[0])+x[1]
        g[1]=x[0]+2*x[1]*c+8*x[1]**3;
        return f,g
    else:
        return f
    
####################################################################

def TestFunction2(x,p,nargout):
    
    a=p[0]
    b=p[1]
    c=p[2]

    den=2*len(x)
    f=sum(a*x**4+b*x**2+c*x)/den+40

    if nargout>1:
        g=(4*a*x**3+2*b*x+c)/den
        return f,g
    else:
        return f
    
####################################################################

def TestFunction3(x,p,nargout):
    
    beta=p
    n=len(x)
    t=np.zeros((n,1))
    s=np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            s[i,j]=(j+1+beta)*(x[j]**(i+1)-(j+1)**(-i-1))
        t[i]=sum(s[i,:])

    f=sum(t**2)

    if nargout>1:
        g=np.zeros((n,1))
        for j in range(n):
            for i in range(n):
                g[j]+=(t[i]*(i+1)*(j+1+beta))*x[j]**i
        return f,g
    else:
        return f

