#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
jupyter_tools.py - provides simple redundant routines for the jupyter files
"""
from random import random

# Class tools...

def print_multiply(a,b,s,r):
    print("""
    a = {}
    b = {}
a {} b = {}
""".format(a,b,s,r))

def random_vector(obj):
    return obj([random() for i in range(obj.dim(obj))])

def random_int_vector(obj,n=20):
    return obj([ n*int(random()-1/2) for i in range(obj.dim(obj)) ])

def random_imaginary_vector(obj):
    return obj([0] + [random() for i in range(obj.dim(obj)-1)])

def unit_list (obj):
    d = obj.dim(obj)
    return [ obj( ( [0]*i + [1] + [0]*(d-i-1) )) for i in range(d) ]

