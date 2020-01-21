#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
involution.albegra Unit Tests

general class wide, product table, algebraic, identity and quantum emulation tests

  file: test.py
source: https://github.com/peawormsworth/PyInvolution
author: Jeffrey B Anderson - truejeffanderson at gmail.com
"""

from involution.algebra import *
import unittest
import pandas as pd
from math import sqrt
from random import random, uniform

VERBOSE = 1
DEBUG   = 0

# Test Data...

complex_table =\
""" 1  i
    i -1  """

dual_table =\
""" 1  i
    i  0  """

split_table =\
""" 1  i
    i  1  """

quaternion_table =\
""" 1  i  j  k
    i -1  k -j
    j -k -1  i
    k  j -i -1  """

split_quaternion_table =\
""" 1  i  j  k
    i -1  k -j
    j -k  1 -i
    k  j  i  1  """

dual_complex_table =\
""" 1  i  j  k
    i -1  k -j
    j -k  0  0
    k  j  0  0  """

hyperbolic_quaternion_table =\
""" 1  i  j  k
    i  1  k  j
    j -k  1 -i
    k -j  i -1  """

# not currently used...
bicomplex_table =\
""" 1  i  j  k
    i -1  k -j
    j  k  1  i
    k -j  i -1  """

# not currently used...
split_quaternion_table =\
""" 1  i  j  k
    i -1  k -j
    j -k  1 -i
    k  j  i  1  """

octonion_table =\
""" 1  i  j  k  l  m  n  o  
    i -1  k -j  m -l -o  n 
    j -k -1  i  n  o -l -m 
    k  j -i -1  o -n  m -l 
    l -m -n -o -1  i  j  k 
    m  l -o  n -i -1 -k  j
    n  o  l -m -j  k -1 -i 
    o -n  m  l -k -j  i -1  """

split_octonion_table =\
""" 1  i  j  k  l  m  n  o  
    i -1  k -j -m  l -o  n 
    j -k -1  i -n  o  l -m 
    k  j -i -1 -o -n  m  l 
    l  m  n  o  1  i  j  k 
    m -l -o  n -i  1  k -j 
    n  o -l -m -j -k  1  i 
    o -n  m -l -k  j -i  1  """

dual_quaternion_table =\
""" 1  i  j  k  l  m  n  o
    i -1  k -j  m -l  o -n
    j -k -1  i  n -o -l  m
    k  j -i -1  o  n -m -l
    l -m -n -o  0  0  0  0
    m  l  o -n  0  0  0  0
    n -o  l  m  0  0  0  0
    o  n -m  l  0  0  0  0  """

sedenion_table =\
""" 1  i  j  k  l  m  n  o  p  q  r  s  t  u  v  w
    i -1  k -j  m -l -o  n  q -p -s  r -u  t  w -v
    j -k -1  i  n  o -l -m  r  s -p -q -v -w  t  u
    k  j -i -1  o -n  m -l  s -r  q -p -w  v -u  t
    l -m -n -o -1  i  j  k  t  u  v  w -p -q -r -s
    m  l -o  n -i -1 -k  j  u -t  w -v  q -p  s -r
    n  o  l -m -j  k -1 -i  v -w -t  u  r -s -p  q
    o -n  m  l -k -j  i -1  w  v -u -t  s  r -q -p
    p -q -r -s -t -u -v -w -1  i  j  k  l  m  n  o
    q  p -s  r -u  t  w -v -i -1 -k  j -m  l  o -n
    r  s  p -q -v -w  t  u -j  k -1 -i -n -o  l  m
    s -r  q  p -w  v -u  t -k -j  i -1 -o  n -m  l
    t  u  v  w  p -q -r -s -l  m  n  o -1 -i -j -k
    u -t  w -v  q  p  s -r -m -l  o -n  i -1  k -j
    v -w -t  u  r -s  p  q -n -o -l  m  j -k -1  i
    w  v -u -t  s  r -q  p -o  n -m -l  k  j -i -1  """

complex_abs_data = [
    [ [ 1, 2], sqrt(5 ) ],
    [ [ 3,-2], sqrt(13) ]
]

complex_add_data = [
    [ [ 1, 2], [ 3,-3], [ 4,-1] ],
    [ [ 3,-3], [ 1, 2], [ 4,-1] ]
]

complex_sub_data = [
    [ [ 1, 2], [ 3,-3], [-2, 5] ],
    [ [ 3,-3], [ 1, 2], [ 2,-5] ]
]

complex_product_data = [
    [ [ 1, 2], [ 3,-3], [ 9, 3] ],
    [ [ 3,-3], [ 1, 2], [ 9, 3] ]
]

complex_division_data = [
    [ [ 1, 2], [ 3,-3], [-1/6, 1/2] ],
    [ [ 3,-3], [ 1, 2], [-3/5,-9/5] ]
]



dual_abs_data = [
    [ [ 1, 2], 1 ],
    [ [ 3,-2], 3 ]
]

dual_add_data = [
    [ [ 1, 2], [ 3,-3], [ 4,-1] ],
    [ [ 3,-3], [ 1, 2], [ 4,-1] ]
]

dual_sub_data = [
    [ [ 1, 2], [ 3,-3], [-2, 5] ],
    [ [ 3,-3], [ 1, 2], [ 2,-5] ]
]

dual_product_data = [
    [ [ 1, 2], [ 3,-3], [ 3, 3] ],
    [ [ 3,-3], [ 1, 2], [ 3, 3] ]
]

dual_division_data = [
    [ [ 1, 2], [ 3,-3], [1/3, 1] ],
    [ [ 3,-3], [ 1, 2], [  3,-9] ]
]


# Data Generators...

complex_abs_record      = [ (Complex(a),expect)                for a,expect in complex_abs_data      ]
complex_add_record      = [ (Complex(a),Complex(b),Complex(c)) for a,b,c    in complex_add_data      ]
complex_sub_record      = [ (Complex(a),Complex(b),Complex(c)) for a,b,c    in complex_sub_data      ]
complex_product_record  = [ (Complex(a),Complex(b),Complex(c)) for a,b,c    in complex_product_data  ]
complex_division_record = [ (Complex(a),Complex(b),Complex(c)) for a,b,c    in complex_division_data ]

dual_abs_record      = [ (Dual(a),expect)          for a,expect in dual_abs_data      ]
dual_add_record      = [ (Dual(a),Dual(b),Dual(c)) for a,b,c    in dual_add_data      ]
dual_sub_record      = [ (Dual(a),Dual(b),Dual(c)) for a,b,c    in dual_sub_data      ]
dual_product_record  = [ (Dual(a),Dual(b),Dual(c)) for a,b,c    in dual_product_data  ]
dual_division_record = [ (Dual(a),Dual(b),Dual(c)) for a,b,c    in dual_division_data ]


# Class tools...

def random_vector(obj):
    return obj([random() for i in range(dim(obj))])

def random_imaginary_vector(obj):
    return obj([0] + [random() for i in range(dim(obj)-1)])

def unit_list (obj):
    d = 2 ** len(obj.dp)
    return [ obj( ( [0]*i + [1] + [0]*(d-i-1) )) for i in range(d) ]

def generate_table (obj):
    """create a multiplication table for a given Algebra object in n×n matrix format (n=dimensions)"""
    units  = unit_list(obj)
    return [ [j*i for i in units] for j in units]

def generate_str (obj):
    """create a multiplication table for a given Algebra object and return the elements in string format"""
    units  = unit_list(obj)
    return [ [str(j*i) for i in units] for j in units]

def dim (obj):
    """the expected number of dimensions given the number of doubling products of this object"""
    return 2**len(obj.dp)

def two_square_identity (x,y):
    a,b = x
    c,d = y
    return [a*c - d*b, d*a + b*c]


def four_square_identity (x,y):
    a,b,c,d = x
    e,f,g,h = y

    r = a*e - b*f - c*g - d*h
    s = a*f + b*e + c*h - d*g
    t = a*g - b*h + c*e + d*f
    u = a*h + b*g - c*f + d*e
    return [r,s,t,u]

def eight_square_identity (x,y):
    a,b,c,d,e,f,g,h = x
    i,j,k,l,m,n,o,p = y

    # assuming i^2 = -1 ...
    z1 = a*i - b*j - c*k - d*l - m*e - n*f - o*g - p*h
    z2 = a*j + b*i + c*l - d*k - m*f + n*e + o*h - p*g
    z3 = a*k - b*l + c*i + d*j - m*g - n*h + o*e + p*f
    z4 = a*l + b*k - c*j + d*i - m*h + n*g - o*f + p*e
    z5 = m*a - n*b - o*c - p*d + e*i + f*j + g*k + h*l
    z6 = m*b + n*a + o*d - p*c - e*j + f*i - g*l + h*k
    z7 = m*c - n*d + o*a + p*b - e*k + f*l + g*i - h*j
    z8 = m*d + n*c - o*b + p*a - e*l - f*k + g*j + h*i
    return [z1,z2,z3,z4,z5,z6,z7,z8]


# this was calculated by hand. There is a fair chance that it is wrong.
def sixteen_square_identity (x,y):
    a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16 = x
    b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16 = y
    # assuming i^2 = -1 ...
    z1  = a1 * b1  - a2  * b2  - a3  * b3  - a4  * b4  - b5  * a5  - b6  * a6  - b7  * a7  - b8  * a8  - b9  * a9  - b10 * a10 - b11 * a11 - b12 * a12 - a13 * b13 - a14 * b14 - a15 * b15 - a16 * b16
    z2  = a1 * b2  + a2  * b1  + a3  * b4  - a4  * b3  - b5  * a6  + b6  * a5  + b7  * a8  - b8  * a7  - b9  * a10 + b10 * a9  + b11 * a12 - b12 * a11 - a13 * b14 + a14 * b13 + a15 * b16 - a16 * b15
    z3  = a1 * b3  - a2  * b4  + a3  * b1  + a4  * b2  - b5  * a7  - b6  * a8  + b7  * a5  + b8  * a6  - b9  * a11 - b10 * a12 + b11 * a9  + b12 * a10 - a13 * b15 - a14 * b16 + a15 * b13 + a16 * b14
    z4  = a1 * b4  + a2  * b3  - a3  * b2  + a4  * b1  - b5  * a8  + b6  * a7  - b7  * a6  + b8  * a5  - b9  * a12 + b10 * a11 - b11 * a10 + b12 * a9  - a13 * b16 + a14 * b15 - a15 * b14 + a16 * b13
    z5  = b5 * a1  - b6  * a2  - b7  * a3  - b8  * a4  + a5  * b1  + a6  * b2  + a7  * b3  + a8  * b4  - a13 * b9  - a14 * b10 - a15 * b11 - a16 * b12 + b13 * a9  + b14 * a10 + b15 * a11 + b16 * a12
    z6  = b5 * a2  + b6  * a1  + b7  * a4  - b8  * a3  - a5  * b2  + a6  * b1  - a7  * b4  + a8  * b3  + a13 * b10 - a14 * b9  + a15 * b12 - a16 * b11 - b13 * a10 + b14 * a9  - b15 * a12 + b16 * a11
    z7  = b5 * a3  - b6  * a4  + b7  * a1  + b8  * a2  - a5  * b3  + a6  * b4  + a7  * b1  - a8  * b2  + a13 * b11 - a14 * b12 - a15 * b9  + a16 * b10 - b13 * a11 + b14 * a12 + b15 * a9  - b16 * a10
    z8  = b5 * a4  + b6  * a3  - b7  * a2  + b8  * a1  - a5  * b4  - a6  * b3  + a7  * b2  + a8  * b1  + a13 * b12 + a14 * b11 - a15 * b10 - a16 * b9  - b13 * a12 - b14 * a11 + b15 * a10 + b16 * a9
    z9  = b9 * a1  - b10 * a2  - b11 * a3  - b12 * a4  - a5  * b13 - a6  * b14 - a7  * b15 - a8  * b16 + a9  * b1  + a10 * b2  + a11 * b3  + a12 * b4  + b5  * a13 + b6  * a14 + b7  * a15 + b8  * a16
    z10 = b9 * a2  + b10 * a1  + b11 * a4  - b12 * a3  - a5  * b14 + a6  * b13 + a7  * b16 - a8  * b15 - a9  * b2  + a10 * b1  - a11 * b4  + a12 * b3  + b5  * a14 - b6  * a13 - b7  * a16 + b8  * a15
    z11 = b9 * a3  - b10 * a4  + b11 * a1  + b12 * a2  - a5  * b15 - a6  * b16 + a7  * b13 + a8  * b14 - a9  * b3  + a10 * b4  + a11 * b1  - a12 * b2  + b5  * a15 + b6  * a16 - b7  * a13 - b8  * a14
    z12 = b9 * a4  + b10 * a3  - b11 * a2  + b12 * a1  - a5  * b16 + a6  * b15 - a7  * b14 + a8  * b13 - a9  * b4  - a10 * b3  + a11 * b2  + a12 * b1  + b5  * a16 - b6  * a15 + b7  * a14 - b8  * a13
    z13 = a5 * b9  - a6  * b10 - a7  * b11 - a8  * b12 + b13 * a1  + b14 * a2  + b15 * a3  + b16 * a4  - b5  * a9  + b6  * a10 + b7  * a11 + b8  * a12 + a13 * b1  - a14 * b2  - a15 * b3  - a16 * b4
    z14 = a5 * b10 + a6  * b9  + a7  * b12 - a8  * b11 - b13 * a2  + b14 * a1  - b15 * a4  + b16 * a3  - b5  * a10 - b6  * a9  - b7  * a12 + b8  * a11 + a13 * b2  + a14 * b1  + a15 * b4  - a16 * b3
    z15 = a5 * b11 - a6  * b12 + a7  * b9  + a8  * b10 - b13 * a3  + b14 * a4  + b15 * a1  - b16 * a2  - b5  * a11 + b6  * a12 - b7  * a9  - b8  * a10 + a13 * b3  - a14 * b4  + a15 * b1  + a16 * b2
    z16 = a5 * b12 + a6  * b11 - a7  * b10 + a8  * b9  - b13 * a4  - b14 * a3  + b15 * a2  + b16 * a1  - b5  * a12 - b6  * a11 + b7  * a10 - b8  * a9  + a13 * b4  + a14 * b3  - a15 * b2  + a16 * b1
    return [z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15,z16]


# Class begins

class TestComplex(unittest.TestCase):
    obj = Complex

    def test_unit_multiplication (self):
        "Complex number unit product table"
        # ref: https://en.wikipedia.org/wiki/Complex_number
        expect = [ a.split() for a in complex_table.split("\n") ]
        calc = generate_str(self.obj)
        if DEBUG or VERBOSE:
            print()
        if DEBUG: print("\ncalc:   %s\nexpect: %s" % (calc, expect))
        if VERBOSE:
            imaginaries = '1ijklmnopqrstuvw'
            df = pd.DataFrame (expect, index = list(imaginaries[0:2]), columns = list(imaginaries[0:2]))
            print ("\n=== Expected %s table ===\n%s" % (self.obj.__name__,df))
            df2 = pd.DataFrame (calc, index = list(imaginaries[0:2]), columns = list(imaginaries[0:2]))
            print ("\n=== Calculated %s table ===\n%s" % (self.obj.__name__,df2))
        if DEBUG or VERBOSE: print('...')
        self.assertListEqual(calc, expect)


    def test_abs_results (self):
        "various absolute results"
        if DEBUG or VERBOSE:
            print()
        for row in complex_abs_record:
            input, expect = row
            calc = abs(input)
            if DEBUG or VERBOSE:
                print("|%s| =\nexpect = %s\n  calc = %s" % (input, expect, calc))
            self.assertEqual(calc, expect)


    def test_add_results (self):
        "various addition results"
        if DEBUG or VERBOSE:
            print()
        for row in complex_add_record:
            a,b,expect = row
            calc = a+b
            if DEBUG or VERBOSE:
                print("\n%s + %s =\nexpect = %s\n  calc = %s" % (a,b, expect, calc))
            self.assertEqual(expect,calc)
            self.assertIsInstance(expect, Complex)


    def test_sub_results (self):
        "various subtraction results"
        if DEBUG or VERBOSE:
            print()
        for row in complex_sub_record:
            a,b,expect = row
            calc = a-b
            if DEBUG or VERBOSE:
                print("\n%s - %s =\nexpect = %s\n  calc = %s" % (a,b, expect, calc))
            self.assertEqual(expect,calc)
            self.assertIsInstance(expect, Complex)


    def test_product_results (self):
        "various product results"
        if DEBUG or VERBOSE:
            print()
        for row in complex_product_record:
            a,b,expect = row
            calc = a*b
            if DEBUG or VERBOSE:
                print("\n%s × %s =\nexpect = %s\n  calc = %s" % (a,b, expect, calc))
            self.assertEqual(expect,calc)
            self.assertIsInstance(expect, Complex)


    def test_division_results (self):
        "various division results"
        if DEBUG or VERBOSE:
            print()
        for row in complex_division_record:
            a,b,expect = row
            calc = a/b
            if DEBUG or VERBOSE:
                print("\n%s / %s =\nexpect = %s\n  calc = %s" % (a,b, expect, calc))
            self.assertEqual(expect,calc)
            self.assertIsInstance(expect, Complex)


class TestDual(unittest.TestCase):
    """
Warning: there is a problem with this test.
The dual numbers do not all have absolute values.
So equality testing does not work.
Equality is based on distance between points.
So the distance between a & b = |a-b|
But since the absolute value of a dual number ignores the imaginary...
two dual numbers that have the same real value will have the same magnitude
regardless of the imaginary quantity.
It means that the dual class thinks (1,0) == (1,1000000)... which can not be right.
Should equality be forced to use pythagorus for distance calculations OR 
is that not allowed?
How can I tell if two dual numbers are so close they're probably be equal?
Until then, these tests are probably garbage.
"""

    obj = Dual

    def test_unit_multiplication (self):
        "Dual number unit product table"
        # ref: https://en.wikipedia.org/wiki/Dual_number
        expect = [ a.split() for a in dual_table.split("\n") ]
        calc = generate_str(self.obj)
        if DEBUG or VERBOSE:
            print()
        if DEBUG: print("\ncalc:   %s\nexpect: %s" % (calc, expect))
        if VERBOSE:
            imaginaries = '1ijklmnopqrstuvw'
            df = pd.DataFrame (expect, index = list(imaginaries[0:2]), columns = list(imaginaries[0:2]))
            print ("\n=== Expected %s table ===\n%s" % (self.obj.__name__,df))
            df2 = pd.DataFrame (calc, index = list(imaginaries[0:2]), columns = list(imaginaries[0:2]))
            print ("\n=== Calculated %s table ===\n%s" % (self.obj.__name__,df2))
        if DEBUG or VERBOSE: print('...')
        self.assertListEqual(calc, expect)


    def test_abs_results (self):
        "various absolute results"
        if DEBUG or VERBOSE:
            print()
        for row in dual_abs_record:
            input, expect = row
            calc = abs(input)
            if DEBUG or VERBOSE:
                print("|%s| =\nexpect = %s\n  calc = %s" % (input, expect, calc))
            self.assertEqual(calc, expect)


    def test_add_results (self):
        "various addition results"
        if DEBUG or VERBOSE:
            print()
        for row in dual_add_record:
            a,b,expect = row
            calc = a+b
            if DEBUG or VERBOSE:
                print("\n%s + %s =\nexpect = %s\n  calc = %s" % (a,b, expect, calc))
            self.assertEqual(expect,calc)
            self.assertIsInstance(expect, self.obj)


    def test_sub_results (self):
        "various subtraction results"
        if DEBUG or VERBOSE:
            print()
        for row in dual_sub_record:
            a,b,expect = row
            calc = a-b
            if DEBUG or VERBOSE:
                print("\n%s - %s =\nexpect = %s\n  calc = %s" % (a,b, expect, calc))
            self.assertEqual(expect,calc)
            self.assertIsInstance(expect, self.obj)


    def test_product_results (self):
        "various product results"
        if DEBUG or VERBOSE:
            print()
        for row in dual_product_record:
            a,b,expect = row
            calc = a*b
            if DEBUG or VERBOSE:
                print("\n%s × %s =\nexpect = %s\n  calc = %s" % (a,b, expect, calc))
            self.assertEqual(expect,calc)
            self.assertIsInstance(expect, self.obj)


    def test_division_results (self):
        "various division results"
        if DEBUG or VERBOSE:
            print()
        for row in dual_division_record:
            a,b,expect = row
            calc = a/b
            if DEBUG or VERBOSE:
                print("\n%s / %s =\nexpect = %s\n  calc = %s" % (a,b, expect, calc))
            self.assertEqual(expect,calc)
            self.assertIsInstance(expect, self.obj)


class TestSplit(unittest.TestCase):
    obj = Split

    def test_unit_multiplication (self):
        "Split number unit product table"
        # ref: https://en.wikipedia.org/wiki/Split_number
        expect = [ a.split() for a in split_table.split("\n") ]
        calc = generate_str(self.obj)
        if DEBUG or VERBOSE:
            print()
        if DEBUG: print("\ncalc:   %s\nexpect: %s" % (calc, expect))
        if VERBOSE:
            imaginaries = '1ijklmnopqrstuvw'
            df = pd.DataFrame (expect, index = list(imaginaries[0:2]), columns = list(imaginaries[0:2]))
            print ("\n=== Expected %s table ===\n%s" % (self.obj.__name__,df))
            df2 = pd.DataFrame (calc, index = list(imaginaries[0:2]), columns = list(imaginaries[0:2]))
            print ("\n=== Calculated %s table ===\n%s" % (self.obj.__name__,df2))
        if DEBUG or VERBOSE: print('...')
        self.assertListEqual(calc, expect)


class TestQuaternion(unittest.TestCase):
    obj = Quaternion

    def test_unit_multiplication (self):
        "Quaternion unit product table"
        # ref: https://en.wikipedia.org/wiki/Quaternion
        expect = [ a.split() for a in quaternion_table.split("\n") ]
        calc = generate_str(self.obj)
        if DEBUG or VERBOSE:
            print()
        if DEBUG: print("\ncalc:   %s\nexpect: %s" % (calc, expect))
        if VERBOSE:
            imaginaries = '1ijklmnopqrstuvw'
            df = pd.DataFrame (expect, index = list(imaginaries[0:4]), columns = list(imaginaries[0:4]))
            print ("\n=== Expected %s table ===\n%s" % (self.obj.__name__,df))
            df2 = pd.DataFrame (calc, index = list(imaginaries[0:4]), columns = list(imaginaries[0:4]))
            print ("\n=== Calculated %s table ===\n%s" % (self.obj.__name__,df2))
        if DEBUG or VERBOSE: print('...')
        self.assertListEqual(calc, expect)


    def test_conjugation (self):
        "Long form congugation"
        # ref: https://en.wikipedia.org/wiki/Quaternion#Conjugation,_the_norm,_and_reciprocal
        loops = 100
        if VERBOSE or DEBUG: 
            print(loops, 'loops')
        for n in range(loops):
            x = random_vector(self.obj)
            calc = x.conj()
            o,i,j,k = unit_list(self.obj)
            expect = -1/2*( x + (i*x)*i + (j*x)*j + (k*x)*k )
            if DEBUG:
                print("   calc = %r" % (calc))
                print(" expect = %r" % (expect))
            self.assertEqual(calc,expect)
        if DEBUG or VERBOSE: 
            print()
            print("\nComparing conjugate against this formula:\nconjugate(x) = \
                 -1/2*( x + (i*x)*i + (j*x)*j + (k*x)*k")
        if VERBOSE:
            print("These should be equal...\n")
            print("   calc = %s" % (calc))
            print(" expect = %s" % (expect))


    def test_conjugate_product (self):
        "the product of a vector with its conjugate is the multidimensional Pythagarus formula"
        # ref: https://en.wikipedia.org/wiki/Octonion#Conjugate,_norm,_and_inverse
        loops = 100
        if VERBOSE or DEBUG: 
            print(loops, 'loops')
        for n in range(loops):
            x = random_vector(self.obj)
            calc = x*x.conj()
            expect = sum([a ** 2 for a in x])
            self.assertEqual(calc,expect)
        if DEBUG or VERBOSE: 
            print()
        if DEBUG: 
            print("   calc = %r" % (calc))
            print(" expect = %r" % (expect))
        if VERBOSE:
            print("These should be equal...\n")
            print("   calc = %s" % (calc))
            print(" expect = %s" % (expect))


    def test_addition_metric_space (self):
        "Test addition is continuous in Quaternion metric topology"
        # ref: https://en.wikipedia.org/wiki/Quaternion#Conjugation,_the_norm,_and_reciprocal
        object = Quaternion
        loops = 100
        if VERBOSE or DEBUG: 
            print(loops, 'loops')
        for n in range(loops):
            p1 = random_vector(self.obj)
            p2 = random_vector(self.obj)
            q1 = random_vector(self.obj)
            q2 = random_vector(self.obj)
            a  = uniform(0,10)
            calc   = abs((p1 + a*p2 + q1 + a*q2) - (p1+q1))
            expect = a*abs(p2+q2)
            if DEBUG:
                print('p1:', p1)
                print('q1:', q1)
                print('p2:', p2)
                print('q2:', q2)
                print(' a:', a)
                print("   calc = %r" % (calc))
                print(" expect = %r" % (expect))
            self.assertAlmostEqual(calc,expect)
        if VERBOSE:
            print()
            print("These should be equal...\n")
            print("   calc = %s" % (calc))
            print(" expect = %s" % (expect))


    def test_dot_product (self):
        "Component vs component free dot product"
        # ref: https://en.wikipedia.org/wiki/Quaternion#Quaternions_and_the_geometry_of_R3
        object = Quaternion
        loops = 100
        if VERBOSE or DEBUG: 
            print(loops, 'loops')
        for n in range(loops):
            p = random_imaginary_vector(self.obj)
            q = random_imaginary_vector(self.obj)
            b1,c1,d1 = p[1:]
            b2,c2,d2 = q[1:]
            calc1 = 1/2*(p.conj()*q + q.conj()*p)
            calc2 = 1/2*(p*q.conj() + q*p.conj())
            expect = b1*b2 + c1*c2 + d1*d2
            if DEBUG:
                print("\nGiven:")
                print('      p:', p)
                print('      q:', q)
                print("  calc1 = %r" % (calc1))
                print("  calc2 = %r" % (calc2))
                print(" expect = %r" % (expect))
            self.assertAlmostEqual(calc1,expect)
            self.assertAlmostEqual(calc2,expect)
            self.assertAlmostEqual(calc1,calc2)
        if VERBOSE:
            print()
            print("These should be equal...\n")
            print("  calc1 = %s" % (calc1))
            print("  calc2 = %s" % (calc2))
            print(" expect = %s" % (expect))


class TestOctonion(unittest.TestCase):
    obj = Octonion

    def test_unit_multiplication (self):
        "Octonion unit product table"
        # ref: https://en.wikipedia.org/wiki/Octonion
        expect = [ a.split() for a in octonion_table.split("\n") ]
        calc = generate_str(self.obj)
        if DEBUG or VERBOSE:
            print()
        if DEBUG: print("\ncalc:   %s\nexpect: %s" % (calc, expect))
        if VERBOSE:
            imaginaries = '1ijklmnopqrstuvw'
            df = pd.DataFrame (expect, index = list(imaginaries[0:8]), columns = list(imaginaries[0:8]))
            print ("\n=== Expected %s table ===\n%s" % (self.obj.__name__,df))
            df2 = pd.DataFrame (calc, index = list(imaginaries[0:8]), columns = list(imaginaries[0:8]))
            print ("\n=== Calculated %s table ===\n%s" % (self.obj.__name__,df2))
        if DEBUG or VERBOSE: print('...')
        self.assertListEqual(calc, expect)

    def test_conjugation (self):
        "Long form congugation"
        # ref: https://en.wikipedia.org/wiki/Octonion#Conjugate,_norm,_and_inverse
        loops = 100
        if VERBOSE or DEBUG: 
            print(loops, 'loops')
        for n in range(loops):
            x = random_vector(self.obj)
            calc = x.conj()
            o,i,j,k,l,m,n,o = unit_list(self.obj)
            expect = -1/6*( x + (i*x)*i + (j*x)*j + (k*x)*k + (l*x)*l + (m*x)*m + (n*x)*n + (o*x)*o )
            if DEBUG:
                print("   calc = %r" % (calc))
                print(" expect = %r" % (expect))
            self.assertEqual(calc,expect)
        if DEBUG or VERBOSE: 
            print()
            print("\nComparing conjugate against this formula:\nconjugate(x) = \
                -1/6*( x + (i*x)*i + (j*x)*j + (k*x)*k + (l*x)*l + (m*x)*m + (n*x)*n + (o*x)*o )\n\n")
        if VERBOSE:
            print("These should be equal...\n")
            print("   calc = %s" % (calc))
            print(" expect = %s" % (expect))


    def test_conjugate_product (self):
        "the product of a vector with its conjugate is the multidimensional Pythagarus formula"
        # ref: https://en.wikipedia.org/wiki/Octonion#Conjugate,_norm,_and_inverse
        loops = 100
        if VERBOSE or DEBUG: 
            print(loops, 'loops')
        for n in range(loops):
            x = random_vector(self.obj)
            calc = x*x.conj()
            expect = sum([a ** 2 for a in x])
            print('calc: %r' % calc)
            print('expect: %r' % expect)
            self.assertEqual(calc,expect)
        if DEBUG or VERBOSE: 
            print()
        if DEBUG: 
            print("   calc = %r" % (calc))
            print(" expect = %r" % (expect))
        if VERBOSE:
            print("These should be equal...\n")
            print("   calc = %s" % (calc))
            print(" expect = %s" % (expect))


class TestSedenion(unittest.TestCase):
    obj = Sedenion

    def test_unit_multiplication (self):
        "Sedenion unit product table"
        # ref: https://en.wikipedia.org/wiki/Sedenion
        expect = [ a.split() for a in sedenion_table.split("\n") ]
        calc = generate_str(self.obj)
        if DEBUG or VERBOSE:
            print()
        if DEBUG: print("\ncalc:   %s\nexpect: %s" % (calc, expect))
        if VERBOSE:
            imaginaries = '1ijklmnopqrstuvw'
            df = pd.DataFrame (expect, index = list(imaginaries[0:16]), columns = list(imaginaries[0:16]))
            print ("\n=== Expected %s table ===\n%s" % (self.obj.__name__,df))
            df2 = pd.DataFrame (calc, index = list(imaginaries[0:16]), columns = list(imaginaries[0:16]))
            print ("\n=== Calculated %s table ===\n%s" % (self.obj.__name__,df2))
        if DEBUG or VERBOSE: print('...')
        self.assertListEqual(calc, expect)

# Exotic number tests

class TestSplitQuaternion(unittest.TestCase):
    obj = SplitQuaternion

    def test_unit_multiplication (self):
        "Split Quaternion unit product table"
        expect = [ a.split() for a in split_quaternion_table.split("\n") ]
        calc = generate_str(self.obj)
        if DEBUG or VERBOSE:
            print()
        if DEBUG: print("\ncalc:   %s\nexpect: %s" % (calc, expect))
        if VERBOSE:
            imaginaries = '1ijklmnopqrstuvw'
            df = pd.DataFrame (expect, index = list(imaginaries[0:4]), columns = list(imaginaries[0:4]))
            print ("\n=== Expected %s table ===\n%s" % (self.obj.__name__,df))
            df2 = pd.DataFrame (calc, index = list(imaginaries[0:4]), columns = list(imaginaries[0:4]))
            print ("\n=== Calculated %s table ===\n%s" % (self.obj.__name__,df2))
        if DEBUG or VERBOSE: print('...')
        self.assertListEqual(calc, expect)

    def test_split_complex_generation (self):
        "Generation Split Quaternion from split-complex numbers"
        # ref: https://en.wikipedia.org/wiki/Split-quaternion#Generation_from_split-complex_numbers
        loops = 100
        if VERBOSE or DEBUG: 
            print(loops, 'loops')
        if DEBUG or VERBOSE: 
            print()
        for n in range(loops):
            q = random_vector(self.obj)
            a = Split(q[:2])
            b = Split(q[2:])
            w,z = a
            y,x = b
            pre_calc = w**2 + x**2 - y**2 - z**2
            calc = Split([pre_calc,0])
            expect = a*a.conj() - b*b.conj()
            if DEBUG:
                print("   calc = %r" % (calc))
                print(" expect = %r" % (expect))
            self.assertEqual(calc,expect)
        if VERBOSE:
            print("Given:")
            print("    q = (a,b) = ((w+zi),(y+xi))")
            print("  conj(x) × x = w² + x² - y² - z²")
            print("\nThese should be equal...\n")
            print("   calc = %s" % (calc))
            print(" expect = %s" % (expect))

    def test_nilpotent(self):
        "Test Split Quaternion nilpotent - a number whose square is zero"
        # ref: https://en.wikipedia.org/wiki/Split-quaternion
        o,i,j,k = unit_list(self.obj)
        for unit in [-k, -j, k, j]:
            q = i-unit
            calc = q*q
            expect = self.obj([0,0,0,0])
            if DEBUG:
                print("      q = %r" % q)
                print("    q×q = %r" % calc)
            if VERBOSE:
                print("\nGiven the nilpotent:")
                print("      q = %s" %  q)
                print("It squares to zero...")
                print("    q×q = %s" % (calc))
                print(" expect = %s" % (expect))
            self.assertAlmostEqual(calc,expect)

    def test_idempotent(self):
        "Test Split Quaternion nilpotent - a number whose square is zero"
        # ref: https://en.wikipedia.org/wiki/Split-quaternion
        o,i,j,k = unit_list(self.obj)
        q = 1/2*(o+j)
        calc = q*q
        expect = q
        if DEBUG:
            print("      q = %r" % q)
            print("    q×q = %r" % calc)
        if VERBOSE:
            print("\nGiven the idempotent:")
            print("      q = %s" %  q)
            print("\nIt's square should be itself...\n")
            print("    q×q = %s" % (calc))
            print(" expect = %s" % (expect))
        self.assertAlmostEqual(calc,expect)


class TestSplitOctonion(unittest.TestCase):
    obj = SplitOctonion

    def test_unit_multiplication (self):
        "SplitOctonion unit product table"
        # ref: https://en.wikipedia.org/wiki/SplitOctonion
        expect = [ a.split() for a in split_octonion_table.split("\n") ]
        calc = generate_str(self.obj)
        if DEBUG or VERBOSE:
            print()
        if DEBUG: print("\ncalc:   %s\nexpect: %s" % (calc, expect))
        if VERBOSE:
            imaginaries = '1ijklmnopqrstuvw'
            df = pd.DataFrame (expect, index = list(imaginaries[0:8]), columns = list(imaginaries[0:8]))
            print ("\n=== Expected %s table ===\n%s" % (self.obj.__name__,df))
            df2 = pd.DataFrame (calc, index = list(imaginaries[0:8]), columns = list(imaginaries[0:8]))
            print ("\n=== Calculated %s table ===\n%s" % (self.obj.__name__,df2))
        if DEBUG or VERBOSE: print('...')
        self.assertListEqual(calc, expect)

    def test_conjugation (self):
        "Long form congugation"
        # ref: https://en.wikipedia.org/wiki/Split-octonion#Conjugate,_norm_and_inverse
        loops = 100
        if VERBOSE or DEBUG: 
            print(loops, 'loops')
        for _ in range(loops):
            z = random_vector(self.obj)
            calc = z.conj()
            x,i,j,k,l,m,n,o = unit_list(self.obj)
            expect = x*z[0] - i*z[1] - j*z[2] - k*z[3]- l*z[4] - m*z[5] - n*z[6] - o*z[7]
            if DEBUG:
                print("   calc = %r" % (calc))
                print(" expect = %r" % (expect))
            self.assertEqual(calc,expect)
        if DEBUG or VERBOSE: 
            print()
            print("""
Comparing conjugate against this formula:
  conjugate(x) = -1/6*( x + (i*x)*i + (j*x)*j + (k*x)*k + (l*x)*l + (m*x)*m + (n*x)*n + (o*x)*o )
            """)
        if VERBOSE:
            print("""
These should be equal...
   calc = %s 
 expect = %s 
            """ % (calc,expect))


    def test_conjugate_product (self):
        "the product of a vector with its conjugate is the multidimensional Pythagarus formula"
        # ref: https://en.wikipedia.org/wiki/Octonion#Conjugate,_norm,_and_inverse
        loops = 100
        if VERBOSE or DEBUG: 
            print(loops, 'loops')
        for n in range(loops):
            x = random_vector(self.obj)
            calc = x*x.conj()
            expect = sum([a ** 2 for a in x[:4]]) - sum([a ** 2 for a in x[4:]])
            self.assertEqual(calc,expect)
        if DEBUG or VERBOSE: 
            print()
        if DEBUG: 
            print("   calc = %r" % (calc))
            print(" expect = %r" % (expect))
        if VERBOSE:
            print("These should be equal...\n")
            print("   calc = %s" % (calc))
            print(" expect = %s" % (expect))



class TestDualComplex(unittest.TestCase):
    obj = DualComplex

    def test_unit_multiplication (self):
        "DualComplex unit product table"
        # ref: https://en.wikipedia.org/wiki/DualComplex
        expect = [ a.split() for a in dual_complex_table.split("\n") ]
        calc = generate_str(self.obj)
        if DEBUG or VERBOSE:
            print()
        if DEBUG: print("\ncalc:   %s\nexpect: %s" % (calc, expect))
        if VERBOSE:
            imaginaries = '1ijklmnopqrstuvw'
            df = pd.DataFrame (expect, index = list(imaginaries[0:4]), columns = list(imaginaries[0:4]))
            print ("\n=== Expected %s table ===\n%s" % (self.obj.__name__,df))
            df2 = pd.DataFrame (calc, index = list(imaginaries[0:4]), columns = list(imaginaries[0:4]))
            print ("\n=== Calculated %s table ===\n%s" % (self.obj.__name__,df2))
        if DEBUG or VERBOSE: print('...')
        self.assertListEqual(calc, expect)


class TestDualQuaternion(unittest.TestCase):
    obj = DualQuaternion

    def test_unit_multiplication (self):
        "DualQuaternion unit product table"
        # ref: https://en.wikipedia.org/wiki/DualQuaternion
        expect = [ a.split() for a in dual_quaternion_table.split("\n") ]
        calc = generate_str(self.obj)
        if DEBUG or VERBOSE:
            print()
        if DEBUG: print("\ncalc:   %s\nexpect: %s" % (calc, expect))
        if VERBOSE:
            imaginaries = '1ijklmnopqrstuvw'
            df = pd.DataFrame (expect, index = list(imaginaries[0:8]), columns = list(imaginaries[0:8]))
            print ("\n=== Expected %s table ===\n%s" % (self.obj.__name__,df))
            df2 = pd.DataFrame (calc, index = list(imaginaries[0:8]), columns = list(imaginaries[0:8]))
            print ("\n=== Calculated %s table ===\n%s" % (self.obj.__name__,df2))
        if DEBUG or VERBOSE: print('...')
        self.assertListEqual(calc, expect)


class TestHyperbolicQuaternion(unittest.TestCase):
    obj = HyperbolicQuaternion

    def test_unit_multiplication (self):
        "HyperbolicQuaternion unit product table"
        # ref: https://en.wikipedia.org/wiki/HyperbolicQuaternion
        expect = [ a.split() for a in hyperbolic_quaternion_table.split("\n") ]
        calc = generate_str(self.obj)
        if DEBUG or VERBOSE:
            print()
        if DEBUG: print("\ncalc:   %s\nexpect: %s" % (calc, expect))
        if VERBOSE:
            imaginaries = '1ijklmnopqrstuvw'
            df = pd.DataFrame (expect, index = list(imaginaries[0:4]), columns = list(imaginaries[0:4]))
            print ("\n=== Expected %s table ===\n%s" % (self.obj.__name__,df))
            df2 = pd.DataFrame (calc, index = list(imaginaries[0:4]), columns = list(imaginaries[0:4]))
            print ("\n=== Calculated %s table ===\n%s" % (self.obj.__name__,df2))
        if DEBUG or VERBOSE: print('...')
        self.assertListEqual(calc, expect)


class Commutative(unittest.TestCase):
    """Commutative Product tests"""

    def _is_commutative(self, obj, loops=100):
        d = dim(obj)
        if VERBOSE or DEBUG: 
            print(loops, 'loops')
        for n in range(loops):
            x = random_vector(obj)
            y = random_vector(obj)
            x * y
            if d > 2:
                self.assertNotEqual(x*y,y*x)
            else:
                self.assertEqual(x*y,y*x)
        if DEBUG or VERBOSE: 
            print()
        if DEBUG: 
            print("x × y = %r" % (x*y))
            print("y × x = %r" % (y*x))
        if VERBOSE:
            print("These should%s be equal...\n" %(' NOT' if d > 2 else ''))
            print("x × y = %s" % (x*y))
            print("y × x = %s" % (y*x))

    def test_complex(self):
        self._is_commutative(Complex, loops=5000)

    def test_quaternion(self):
        self._is_commutative(Quaternion, loops=1000)

    def test_octonion(self):
        self._is_commutative(Octonion, loops=300)

    def test_sedenion(self):
        self._is_commutative(Sedenion)

    def test_cd32(self):
        self._is_commutative(Cd32,loops=50)

    # this did not work out for now...
    #def test_split_octonion(self):
        #self._is_commutative(SplitOctonion, loops=3)


class WeakAlternativeCondition(unittest.TestCase):

    def _is_weak_alternative(self, obj, loops=100):
        "Weak Alternative Condition tests"
        d = dim(obj)
        if VERBOSE or DEBUG: 
            print(loops, 'loops')
        for n in range(loops):
            x = random_vector(obj)
            y = random_vector(obj)
            #print('x,y: ',x,y)
            if d < 16:
                self.assertEqual((y*x)*x, y*(x*x))
                self.assertEqual((x*y)*x, x*(y*x))
                self.assertEqual((x*x)*y, x*(x*y))
                self.assertEqual((x*y)*y, x*(y*y))
                self.assertEqual((y*x)*y, y*(x*y))
                self.assertEqual((y*y)*x, y*(y*x))
            else:
                self.assertEqual((x*y)*x, x*(y*x))
                self.assertEqual((y*x)*y, y*(x*y))
                self.assertNotEqual((y*x)*x, y*(x*x))
                self.assertNotEqual((x*x)*y, x*(x*y))
                self.assertNotEqual((x*y)*y, x*(y*y))
                self.assertNotEqual((y*y)*x, y*(y*x))
        if VERBOSE:
            print("\nGiven:\n")
            print("%12s : %s\n" % ('          x ',  x         ))
            print("%12s : %s\n" % ('          y ',  y         ))
            print("These should%s be equal...\n" %(' NOT' if d > 8 else ''))
            print("%12s : %s\n" % ('(y × x) × x ', (y* x) *x  ))
            print("%12s : %s\n" % (' y ×(x  × x)',  y*(x  *x) ))
            print("These should%s be equal...\n" %(' NOT' if d > 8 else ''))
            print("%12s : %s\n" % ('(x × y) × x ', (x* y) *x  ))
            print("%12s : %s\n" % (' x ×(y  × x ',  x*(y  *x) ))
            print("These should%s be equal...\n" %(' NOT' if d > 8 else ''))
            print("%12s : %s\n" % ('(x × x) × y ', (x* x) *y  ))
            print("%12s : %s\n" % (' x ×(x  × y)',  x*(x  *y) ))
            print("These should%s be equal...\n" %(' NOT' if d > 8 else ''))
            print("%12s : %s\n" % ('(x × y) × y ', (x* y) *y  ))
            print("%12s : %s\n" % (' x ×(y  × y)',  x*(y  *y) ))
            print("These should%s be equal...\n" %(' NOT' if d > 8 else ''))
            print("%12s : %s\n" % ('(y × x) × y ', (y* x) *y  ))
            print("%12s : %s\n" % (' y ×(x  × y)',  y*(x  *y) ))
            print("These should%s be equal...\n" %(' NOT' if d > 8 else ''))
            print("%12s : %s\n" % ('(y × y) × x ', (y* y) *x  ))
            print("%12s : %s\n" % (' y ×(y  × x)',  y*(y  *x) ))

    def test_complex(self):
        self._is_weak_alternative(Complex,loops=1000)

    def test_quaternion(self):
        self._is_weak_alternative(Quaternion, loops=300)

    def test_octonion(self):
        self._is_weak_alternative(Octonion)

    def test_sedenion(self):
        self._is_weak_alternative(Sedenion,loops=15)

    def test_32ion(self):
        self._is_weak_alternative(Cd32,loops=10)

    #def test_split_octonion(self):
        #self._is_moufang_condition(SplitOctonion, loops=300)


class DiophantusIdentity(unittest.TestCase):
    """
    Diophantus identity test
    Brahmagupta–Fibonacci / Diophantus identity
    """

    def _is_diophantus_identity(self, obj, loops=100):
        "Diophantus identity test"
        PRECISION = 10**-9
        d = dim(obj)
        if VERBOSE or DEBUG: 
            print(loops, 'loops')
        for n in range(loops):
            x = random_vector(obj)
            y = random_vector(obj)
            if d > 8:
                self.assertNotAlmostEqual(abs(x) * abs(y), abs(x*y), delta=PRECISION)
            else:
                self.assertAlmostEqual(abs(x) * abs(y), abs(x*y), delta=PRECISION)
        if VERBOSE or DEBUG:
            print()
        if DEBUG:
            print("%8s: %r\n" % (            'x',  x            ))
            print("%8s: %r\n" % (            'y',  y            ))
            print("%8s: %r\n" % (        'abs(x)',abs(x)        ))
        if VERBOSE:
            print("\nGiven:\n")
            print("%8s: %s\n" % (            'x',  x            ))
            print("%8s: %s\n" % (            'y',  y            ))
            print("\nHaving absolute values:\n")
            print("%8s: %s\n" % (        'abs(x)',abs(x)        ))
            print("%8s: %s\n" % (        'abs(y)',abs(y)        ))
            print("\nThese should%s be equal...\n" %(' NOT' if d > 8 else ''))
            print("%8s: %s\n" % ( 'abs(x)×abs(y)',abs(x)*abs(y) ))
            print("%8s: %s\n" % (      'abs(x×y)',abs(x*y)      ))


    def test_complex(self):
        self._is_diophantus_identity(Complex,loops=20000)

    def test_quaternion(self):
        self._is_diophantus_identity(Quaternion, loops=6000)

    def test_octonion(self):
        self._is_diophantus_identity(Octonion, loops=2000)

    def test_sedenion(self):
        self._is_diophantus_identity(Sedenion,loops=300)

    def test_32ion(self):
        self._is_diophantus_identity(Cd32,loops=200)

    #def test_split_octonion(self):
        #self._is_diophantus_identity(SplitOctonion, loops=300)


class MoufangCondition(unittest.TestCase):
    """Moufang condition tests"""

    def _is_moufang_condition(self, obj, loops=100):
        "Moufang condition tests"
        dm = dim(obj)
        if VERBOSE or DEBUG: 
            print(loops, 'loops')
        for n in range(loops):
            x = random_vector(obj)
            y = random_vector(obj)
            z = random_vector(obj)
            a = z*(x*(z*y))
            b = ((z*x)*z)*y
            c = x*(z*(y*z))
            d = ((x*z)*y)*z
            e = (z*x)*(y*z)
            f = (z*(x*y))*z
            g = (z*x)*(y*z)
            h = z*((x*y)*z)
            if dm > 8:
                self.assertNotEqual(a,b)
                self.assertNotEqual(c,d)
                self.assertNotEqual(e,f)
                self.assertNotEqual(g,h)
            else:
                self.assertEqual(a,b)
                self.assertEqual(c,d)
                self.assertEqual(e,f)
                self.assertEqual(g,h)
        if DEBUG or VERBOSE:
            print()
        if DEBUG:
            print("%8s: %r\n" % ('x',x))
            print("%8s: %r\n" % ('y',y))
            print("%8s: %r\n" % ('z',z))
        if VERBOSE:
            print("\n")
            a = z*(x*(z*y))
            b = ((z*x)*z)*y
            c = x*(z*(y*z))
            d = ((x*z)*y)*z
            e = (z*x)*(y*z)
            f = (z*(x*y))*z
            g = (z*x)*(y*z)
            h = z*((x*y)*z)
            print("\nGiven:\n")
            print("%8s: %s\n" % ('x',x))
            print("%8s: %s\n" % ('y',y))
            print("%8s: %s\n" % ('z',z))
            print("\nThese should%s be equal...\n" %(' NOT' if dm > 8 else ''))
            print("%8s: %s\n" % ('z × (x × (z × y))' , a) )
            print("%8s: %s\n" % ('((z  × x) × z) × y', b) )
            print("\nThese should%s be equal...\n" %(' NOT' if dm > 8 else ''))
            print("%8s: %s\n" % ('x × (z × (y × z ))', c) )
            print("%8s: %s\n" % ('(( x × z) × y) × z', d) )
            print("\nThese should%s be equal...\n" %(' NOT' if dm > 8 else ''))
            print("%8s: %s\n" % ('(z × x) × (y × z )', e) )
            print("%8s: %s\n" % ('(z × (x × y))× z'  , f) )
            print("\nThese should%s be equal...\n" %(' NOT' if dm > 8 else ''))
            print("%8s: %s\n" % ('(z × x) × (y × z )', g) )
            print("%8s: %s\n" % ('z ×((x × y) × z )' , h) )
            print("\nThese should%s be equal...\n" %(' NOT' if dm > 8 else ''))


    def test_complex(self):
        self._is_moufang_condition(Complex,loops=500)

    def test_quaternion(self):
        self._is_moufang_condition(Quaternion, loops=150)

    def test_octonion(self):
        self._is_moufang_condition(Octonion, loops=50)

    def test_sedenion(self):
        self._is_moufang_condition(Sedenion,loops=15)

    def test_32ion(self):
        self._is_moufang_condition(Cd32,loops=10)

    #def test_split_octonion(self):
        #self._is_moufang_condition(SplitOctonion, loops=300)


class PowerAssociative(unittest.TestCase):
    """Power Associative tests"""

    def _is_power_associative(self, obj, loops=100):
        PRECISION = 10**-9
        d = dim(obj)
        if VERBOSE or DEBUG: 
            print(loops, 'loops')
        for n in range(loops):
            x = random_vector(obj).normalize()
            y = random_vector(obj).normalize()
            z = x * y
            self.assertAlmostEqual(abs(x),1,delta=PRECISION)
            self.assertAlmostEqual(abs(y),1,delta=PRECISION)
            if d > 8:
                # this logic is questionable!!!
                self.assertAlmostEqual(abs(z),1,delta=10**-1)
            else:
                self.assertAlmostEqual(abs(z),1,delta=PRECISION)
        if DEBUG or VERBOSE:
            print()
        if DEBUG:
            print("%8s: %r\n" % ('x',  x))
            print("%8s: %r\n" % ('y',  y))
            print("%8s: %r\n" % ('x*y',z))
        if VERBOSE:
            print("\n")
            print("\nGiven random unit vectors...\n")
            print("%8s: %s\n" % ('x',  x))
            print("%8s: %s\n" % ('y',  y))
            print("%8s: %s\n" % ('x*y',z))
            print("\nHaving a magnitude (abs) of 1...\n")
            print("%8s: %s\n" % ('abs(x)',abs(x)))
            print("%8s: %s\n" % ('abs(y)',abs(y)))
            print("\nProduce a product of magnitude 1...\n")
            print("%8s: %s\n" % ('abs(x*y)',abs(z)))
            print("\nThey should be in the same spot...\n")
            print("%8s: %s\n" % ('distance to 1',abs(z)-1))
        if DEBUG or VERBOSE:
            print()

    def test_complex(self):
        self._is_power_associative(Complex,loops=8000)

    def test_quaternion(self):
        self._is_power_associative(Quaternion, loops=2400)

    def test_octonion(self):
        self._is_power_associative(Octonion, loops=800)

    def test_sedenion(self):
        self._is_power_associative(Sedenion,loops=120)

    def test_32ion(self):
        self._is_power_associative(Cd32,loops=80)

    #def test_split_octonion(self):
        #self._is_power_associative(SplitOctonion, loops=300)


class TwoSquareIdentity(unittest.TestCase):
    """Two square identity"""

    def test_two_square_identity(self):
        obj = Complex
        x = random_vector(obj)
        y = random_vector(obj)
        formula = obj(two_square_identity(x,y))
        if VERBOSE:
            print("\nGiven:\n")
            print('           x: ', x)
            print('           y: ', y)
            print("\nThese should be the same...\n")
            print("Brahmagupta-Fibonacci's Two-Square Identity test product is:\n")
            print('formula(y*x): ', formula)
            print("\nThe Involution product is:\n")
            print('   code(y*x): ', y*x)
        self.assertEqual(formula, x*y)


class FourSquareIdentity(unittest.TestCase):
    """Four square identity"""

    def test_four_square_identity(self):
        obj = Quaternion
        x = random_vector(obj)
        y = random_vector(obj)
        formula = obj(four_square_identity(x,y))
        if VERBOSE:
            print("\nGiven:\n")
            print('           x: ', x)
            print('           y: ', y)
            print("\nThese should be the same...\n")
            print("Euler's Four-Square Identity test product is:\n")
            print('formula(x*y): ', formula)
            print("\nThe Involution product is:\n")
            print('   code(x*y): ', x*y)
        self.assertEqual(formula, x*y)

class EightSquareIdentity(unittest.TestCase):
    """Eight square identity"""

    def test_eight_square_identity(self):
        obj = Octonion
        x = random_vector(obj)
        y = random_vector(obj)
        formula = obj(eight_square_identity(x,y))
        if VERBOSE:
            print("\nGiven:\n")
            print('           x: ', x)
            print('           y: ', y)
            print("\nThese should be the same...\n")
            print("Degen's Eight-Square Identity test product is:\n")
            print('formula(x*y): ', formula)
            print("\nThe Involution product is:\n")
            print('   code(x*y): ', x*y)
        self.assertEqual(formula, x*y)


class SixteenSquareIdentity(unittest.TestCase):
    """Sixteen square identity"""

    def test_sixteen_square_identity(self):
        obj = Sedenion
        x = random_vector(obj)
        y = random_vector(obj)
        formula = obj(sixteen_square_identity(x,y))
        if VERBOSE:
            print("\nGiven:\n")
            print('           x: ', x)
            print('           y: ', y)
            print("\nThese should be the same...\n")
            print("Pfister's Sixteen-Square Identity test product is:\n")
            print('formula(x*y): ', formula)
            print("\nThe Involution product is:\n")
            print('   code(x*y): ', x*y)
        self.assertAlmostEqual(formula, x*y)


# run all tests on execute...

if __name__ == '__main__':

    unittest.main()
