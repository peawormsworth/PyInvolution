#!/usr/bin/python3
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
from random import random

VERBOSE = 0
DEBUG   = 0

def random_vector(object):
    return object([random() for i in range(dim(object))])

def unit_list (object):
    d = 2 ** len(object.dp)
    return [ object( ( [0]*i + [1] + [0]*(d-i-1) )) for i in range(d) ]

def generate_table (object):
    units  = unit_list(object)
    return [ [j*i for i in units] for j in units]

def generate_str (object):
    units  = unit_list(object)
    return [ [str(j*i) for i in units] for j in units]

def dim (object):
    return 2**len(object.dp)


class TestComplex(unittest.TestCase):

    object = Complex

    def test_unit_multiplication (self):
        "Complex number unit product table"
        # ref: https://en.wikipedia.org/wiki/Complex_number
        expect = [
            ' 1  i '.split(), 
            ' i -1 '.split()
        ]
        calc = generate_str(self.object)
        if DEBUG or VERBOSE: print()
        if DEBUG: print("\ncalc:   %s\nexpect: %s" % (calc, expect))
        if VERBOSE:
            imaginaries = '1ijklmnopqrstuvw'
            df = pd.DataFrame (expect, index = list(imaginaries[0:2]), columns = list(imaginaries[0:2]))
            print ("\n=== Expected %s table ===\n%s" % (self.object.__name__,df))
            df2 = pd.DataFrame (calc, index = list(imaginaries[0:2]), columns = list(imaginaries[0:2]))
            print ("\n=== Calculated %s table ===\n%s" % (self.object.__name__,df2))
        if DEBUG or VERBOSE: print('...')
        self.assertListEqual(calc, expect)

    def test_abs_results (self):
        "various absolute results"
        a = Complex([ 1, 2])
        b = Complex([ 3,-2])
        if DEBUG or VERBOSE:
            print()
            print("abs(%s):\nexpect = %s\n  calc = %s" % (a, sqrt(5), abs(a)))
            print("abs(%s):\nexpect = %s\n  calc = %s" % (b, sqrt(13), abs(b)))
        self.assertEqual(abs(a), sqrt(5))
        self.assertEqual(abs(b), sqrt(13))

    def test_add_results (self):
        "various addition results"
        a = Complex([1, 2])
        b = Complex([3,-3])
        expect_a = Complex([4,-1])
        expect_b = Complex([4,-1])
        if DEBUG or VERBOSE:
            print()
            print("%s + %s:\nexpect = %s\n  calc = %s" % (a,b, expect_a, a+b))
            print("%s + %s:\nexpect = %s\n  calc = %s" % (b,a, expect_b, b+a))
        self.assertEqual(a+b, expect_a)
        self.assertEqual(b+a, expect_b)
        self.assertIsInstance(a+b, Complex)

    def test_sub_results (self):
        "various subtraction results"
        a = Complex([1, 2])
        b = Complex([3,-3])
        expect_a = Complex([-2, 5])
        expect_b = Complex([ 2,-5])
        if DEBUG or VERBOSE:
            print()
            print("%s - %s:\nexpect = %s\n  calc = %s" % (a,b, expect_a, a-b))
            print("%s - %s:\nexpect = %s\n  calc = %s" % (b,a, expect_b, b-a))
        self.assertEqual(a-b, expect_a)
        self.assertEqual(b-a, expect_b)
        self.assertIsInstance(a-b, Complex)

    def test_product_results (self):
        "various product results"
        a = Complex([1, 2])
        b = Complex([3,-3])
        expect_a = Complex([ 9, 3])
        expect_b = Complex([ 9, 3])
        if DEBUG or VERBOSE:
            print()
            print("%s * %s:\nexpect = %s\n  calc = %s" % (a,b, expect_a, a*b))
            print("%s * %s:\nexpect = %s\n  calc = %s" % (b,a, expect_b, b*a))
        self.assertEqual(a*b, expect_a)
        self.assertEqual(b*a, expect_b)
        self.assertIsInstance(a*b, Complex)

    def test_division_results (self):
        "various division results"
        a = Complex([1, 2])
        b = Complex([3,-3])
        expect_a = Complex([-1/6, 1/2])
        expect_b = Complex([-3/5,-9/5])
        if DEBUG or VERBOSE:
            print()
            print("%s / %s:\nexpect = %s\n  calc = %s" % (a,b, expect_a, a/b))
            print("%s / %s:\nexpect = %s\n  calc = %s" % (b,a, expect_b, b/a))
        self.assertEqual(a/b, expect_a)
        self.assertEqual(b/a, expect_b)
        self.assertIsInstance(a/b, Complex)


class TestDual(unittest.TestCase):
    def test_unit_multiplication (self):
        "Dual number unit product table"
        object = Dual
        # ref: https://en.wikipedia.org/wiki/Dual_number
        expect = [
            ' 1  i '.split(), 
            ' i  0 '.split()
        ]
        calc = generate_str(object)
        if DEBUG or VERBOSE: print()
        if DEBUG: print("\ncalc:   %s\nexpect: %s" % (calc, expect))
        if VERBOSE:
            imaginaries = '1ijklmnopqrstuvw'
            df = pd.DataFrame (expect, index = list(imaginaries[0:2]), columns = list(imaginaries[0:2]))
            print ("\n=== Expected %s table ===\n%s" % (object.__name__,df))
            df2 = pd.DataFrame (calc, index = list(imaginaries[0:2]), columns = list(imaginaries[0:2]))
            print ("\n=== Calculated %s table ===\n%s" % (object.__name__,df2))
        if DEBUG or VERBOSE: print('...')
        self.assertListEqual(calc, expect)

    def test_abs_results (self):
        "various absolute results"
        a = Dual([1,2])
        b = Dual([3,-2])
        expect_a = 1
        expect_b = 3
        if DEBUG or VERBOSE:
            print()
            print("abs(%s):\nexpect = %s\n  calc = %s" % (a, expect_a, abs(a)))
            print("abs(%s):\nexpect = %s\n  calc = %s" % (b, expect_b, abs(b)))
        self.assertEqual(abs(a), expect_a)
        self.assertEqual(abs(b), expect_b)

    def test_add_results (self):
        "various addition results"
        a = Dual([1, 2])
        b = Dual([3,-3])
        expect_a = Dual([4,-1])
        expect_b = Dual([4,-1])
        if DEBUG or VERBOSE:
            print()
            print("%s + %s:\nexpect = %s\n  calc = %s" % (a,b, expect_a, a+b))
            print("%s + %s:\nexpect = %s\n  calc = %s" % (b,a, expect_b, b+a))
        self.assertEqual(a+b, expect_a)
        self.assertEqual(b+a, expect_b)
        self.assertIsInstance(a+b, Dual)


    def test_sub_results (self):
        "various subtraction results"
        a = Dual([1, 2])
        b = Dual([3,-3])
        expect_a = Dual([-2, 5])
        expect_b = Dual([ 2,-5])
        if DEBUG or VERBOSE:
            print()
            print("%s - %s:\nexpect = %s\n  calc = %s" % (a,b, expect_a, a-b))
            print("%s - %s:\nexpect = %s\n  calc = %s" % (b,a, expect_b, b-a))
        self.assertEqual(a-b, expect_a)
        self.assertEqual(b-a, expect_b)
        self.assertIsInstance(a-b, Dual)


    def test_product_results (self):
        "various product results"
        a = Dual([1, 2])
        b = Dual([3,-3])
        expect_a = Dual([ 3, 3])
        expect_b = Dual([ 3, 3])
        if DEBUG or VERBOSE:
            print()
            print("%s * %s:\nexpect = %s\n  calc = %s" % (a,b, expect_a, a*b))
            print("%s * %s:\nexpect = %s\n  calc = %s" % (b,a, expect_b, b*a))
        self.assertEqual(a*b, expect_a)
        self.assertEqual(b*a, expect_b)
        self.assertIsInstance(a*b, Dual)


    def test_division_results (self):
        "various division results"
        a = Dual([1, 2])
        b = Dual([3,-3])
        expect_a = Dual([1/3,1])
        #expect_a = Dual([-1/6, 1/2])
        expect_b = Dual([3,3])
        if DEBUG or VERBOSE:
            print()
            print("%s / %s:\nexpect = %s\n  calc = %s" % (a,b, expect_a, a/b))
            print("%s / %s:\nexpect = %s\n  calc = %s" % (b,a, expect_b, b/a))
        self.assertEqual(a/b, expect_a)
        self.assertEqual(b/a, expect_b)
        self.assertIsInstance(a/b, Dual)


class TestSplit(unittest.TestCase):

    def test_unit_multiplication (self):
        "Split number unit product table"
        object = Split
        # ref: https://en.wikipedia.org/wiki/Split-complex_number
        expect = [
            ' 1  i '.split(), 
            ' i  1 '.split()
        ]
        calc = generate_str(object)
        if DEBUG or VERBOSE: print()
        if DEBUG: print("\ncalc:   %s\nexpect: %s" % (calc, expect))
        if VERBOSE:
            imaginaries = '1ijklmnopqrstuvw'
            df = pd.DataFrame (expect, index = list(imaginaries[0:2]), columns = list(imaginaries[0:2]))
            print ("\n=== Expected %s table ===\n%s" % (object.__name__,df))
            df2 = pd.DataFrame (calc, index = list(imaginaries[0:2]), columns = list(imaginaries[0:2]))
            print ("\n=== Calculated %s table ===\n%s" % (object.__name__,df2))
        if DEBUG or VERBOSE: print('...')
        self.assertListEqual(calc, expect)


class TestQuaternion(unittest.TestCase):

    def test_unit_multiplication (self):
        "Quaternion unit product table"
        object = Quaternion
        # ref: https://en.wikipedia.org/wiki/Quaternion
        expect = [
            ' 1  i  j  k '.split(), 
            ' i -1  k -j '.split(),
            ' j -k -1  i '.split(),
            ' k  j -i -1 '.split()
        ]
        calc = generate_str(object)
        if DEBUG or VERBOSE: print()
        if DEBUG: print("\ncalc:   %s\nexpect: %s" % (calc, expect))
        if VERBOSE:
            imaginaries = '1ijklmnopqrstuvw'
            df = pd.DataFrame (expect, index = list(imaginaries[0:4]), columns = list(imaginaries[0:4]))
            print ("\n=== Expected %s table ===\n%s" % (object.__name__,df))
            df2 = pd.DataFrame (calc, index = list(imaginaries[0:4]), columns = list(imaginaries[0:4]))
            print ("\n=== Calculated %s table ===\n%s" % (object.__name__,df2))
        if DEBUG or VERBOSE: print('...')
        self.assertListEqual(calc, expect)


class TestOctonion(unittest.TestCase):

    def test_unit_multiplication (self):
        "Octonion unit product table"
        object = Octonion 
        # ref: https://en.wikipedia.org/wiki/Octonion
        expect = [
            ' 1  i  j  k  l  m  n  o '.split(), 
            ' i -1  k -j  m -l -o  n '.split(),
            ' j -k -1  i  n  o -l -m '.split(),
            ' k  j -i -1  o -n  m -l '.split(),
            ' l -m -n -o -1  i  j  k '.split(),
            ' m  l -o  n -i -1 -k  j '.split(),
            ' n  o  l -m -j  k -1 -i '.split(),
            ' o -n  m  l -k -j  i -1 '.split()
        ]
        calc = generate_str(object)
        if DEBUG or VERBOSE: print()
        if DEBUG: print("\ncalc:   %s\nexpect: %s" % (calc, expect))
        if VERBOSE:
            imaginaries = '1ijklmnopqrstuvw'
            df = pd.DataFrame (expect, index = list(imaginaries[0:8]), columns = list(imaginaries[0:8]))
            print ("\n=== Expected %s table ===\n%s" % (object.__name__,df))
            df2 = pd.DataFrame (calc, index = list(imaginaries[0:8]), columns = list(imaginaries[0:8]))
            print ("\n=== Calculated %s table ===\n%s" % (object.__name__,df2))
        if DEBUG or VERBOSE: print('...')
        self.assertListEqual(calc, expect)



class TestSedenion(unittest.TestCase):

    def test_unit_multiplication (self):
        "Sedenion unit product table"
        object = Sedenion 
        # ref: https://en.wikipedia.org/wiki/Sedenion
        expect = [
            ' 1  i  j  k  l  m  n  o  p  q  r  s  t  u  v  w'.split(),
            ' i -1  k -j  m -l -o  n  q -p -s  r -u  t  w -v'.split(),
            ' j -k -1  i  n  o -l -m  r  s -p -q -v -w  t  u'.split(),
            ' k  j -i -1  o -n  m -l  s -r  q -p -w  v -u  t'.split(),
            ' l -m -n -o -1  i  j  k  t  u  v  w -p -q -r -s'.split(),
            ' m  l -o  n -i -1 -k  j  u -t  w -v  q -p  s -r'.split(),
            ' n  o  l -m -j  k -1 -i  v -w -t  u  r -s -p  q'.split(),
            ' o -n  m  l -k -j  i -1  w  v -u -t  s  r -q -p'.split(),
            ' p -q -r -s -t -u -v -w -1  i  j  k  l  m  n  o'.split(),
            ' q  p -s  r -u  t  w -v -i -1 -k  j -m  l  o -n'.split(),
            ' r  s  p -q -v -w  t  u -j  k -1 -i -n -o  l  m'.split(),
            ' s -r  q  p -w  v -u  t -k -j  i -1 -o  n -m  l'.split(),
            ' t  u  v  w  p -q -r -s -l  m  n  o -1 -i -j -k'.split(),
            ' u -t  w -v  q  p  s -r -m -l  o -n  i -1  k -j'.split(),
            ' v -w -t  u  r -s  p  q -n -o -l  m  j -k -1  i'.split(),
            ' w  v -u -t  s  r -q  p -o  n -m -l  k  j -i -1'.split()
        ]
        calc = generate_str(object)
        if DEBUG or VERBOSE: print()
        if DEBUG: print("\ncalc:   %s\nexpect: %s" % (calc, expect))
        if VERBOSE:
            imaginaries = '1ijklmnopqrstuvw'
            df = pd.DataFrame (expect, index = list(imaginaries[0:16]), columns = list(imaginaries[0:16]))
            print ("\n=== Expected %s table ===\n%s" % (object.__name__,df))
            df2 = pd.DataFrame (calc, index = list(imaginaries[0:16]), columns = list(imaginaries[0:16]))
            print ("\n=== Calculated %s table ===\n%s" % (object.__name__,df2))
        if DEBUG or VERBOSE: print('...')
        self.assertListEqual(calc, expect)


class TestSplitQuaternion(unittest.TestCase):

    def test_unit_multiplication (self):
        "SplitQuaternion unit product table"
        object = SplitQuaternion 
        # ref: https://en.wikipedia.org/wiki/Split-quaternion
        expect = [
            ' 1  i  j  k '.split(), 
            ' i -1  k -j '.split(),
            ' j -k  1 -i '.split(),
            ' k  j  i  1 '.split()
        ]
        calc = generate_str(object)
        if DEBUG or VERBOSE: print()
        if DEBUG: print("\ncalc:   %s\nexpect: %s" % (calc, expect))
        if VERBOSE:
            imaginaries = '1ijklmnopqrstuvw'
            df = pd.DataFrame (expect, index = list(imaginaries[0:4]), columns = list(imaginaries[0:4]))
            print ("\n=== Expected %s table ===\n%s" % (object.__name__,df))
            df2 = pd.DataFrame (calc, index = list(imaginaries[0:4]), columns = list(imaginaries[0:4]))
            print ("\n=== Calculated %s table ===\n%s" % (object.__name__,df2))
        if DEBUG or VERBOSE: print('...')
        self.assertListEqual(calc, expect)


class TestSplitOctonion(unittest.TestCase):

    def test_unit_multiplication (self):
        "SplitOctonion unit product table"
        object = SplitOctonion 
        # ref: https://en.wikipedia.org/wiki/Split-octonion
        expect = [
            ' 1  i  j  k  l  m  n  o '.split(), 
            ' i -1  k -j -m  l -o  n '.split(),
            ' j -k -1  i -n  o  l -m '.split(),
            ' k  j -i -1 -o -n  m  l '.split(),
            ' l  m  n  o  1  i  j  k '.split(),
            ' m -l -o  n -i  1  k -j '.split(),
            ' n  o -l -m -j -k  1  i '.split(),
            ' o -n  m -l -k  j -i  1 '.split()
        ]
        calc = generate_str(object)
        if DEBUG or VERBOSE: print()
        if DEBUG: print("\ncalc:   %s\nexpect: %s" % (calc, expect))
        if VERBOSE:
            imaginaries = '1ijklmnopqrstuvw'
            df = pd.DataFrame (expect, index = list(imaginaries[0:8]), columns = list(imaginaries[0:8]))
            print ("\n=== Expected %s table ===\n%s" % (object.__name__,df))
            df2 = pd.DataFrame (calc, index = list(imaginaries[0:8]), columns = list(imaginaries[0:8]))
            print ("\n=== Calculated %s table ===\n%s" % (object.__name__,df2))
        if DEBUG or VERBOSE: print('...')
        self.assertListEqual(calc, expect)


class TestDualComplex(unittest.TestCase):

    def test_unit_multiplication (self):
        "DualComplex unit product table"
        object = DualComplex
        # ref: https://en.wikipedia.org/wiki/Dual-complex_number
        expect = [
            ' 1  i  j  k '.split(), 
            ' i -1  k -j '.split(),
            ' j -k  0  0 '.split(),
            ' k  j  0  0 '.split()
        ]
        calc = generate_str(object)
        if DEBUG or VERBOSE: print()
        if DEBUG: print("\ncalc:   %s\nexpect: %s" % (calc, expect))
        if VERBOSE:
            imaginaries = '1ijklmnopqrstuvw'
            df = pd.DataFrame (expect, index = list(imaginaries[0:4]), columns = list(imaginaries[0:4]))
            print ("\n=== Expected %s table ===\n%s" % (object.__name__,df))
            df2 = pd.DataFrame (calc, index = list(imaginaries[0:4]), columns = list(imaginaries[0:4]))
            print ("\n=== Calculated %s table ===\n%s" % (object.__name__,df2))
        if DEBUG or VERBOSE: print('...')
        self.assertListEqual(calc, expect)


class TestDualQuaternion(unittest.TestCase):

    def test_unit_multiplication (self):
        "DualQuaternion unit product table"
        object = DualQuaternion 
        # ref: https://en.wikipedia.org/wiki/Dual_quaternion
        UNexpectEDfromWIKI = [
            ' 1  i  j  k  l  m  n  o '.split(), 
            ' i -1  k -j  m -l  o -n '.split(),
            ' j -k -1  i  n -o -l  m '.split(),
            ' k  j -i -1  o  n -m -l '.split(),
            ' l  m  n  o  0  0  0  0 '.split(),
            ' m -l  o -n  0  0  0  0 '.split(),
            ' n -o -l  m  0  0  0  0 '.split(),
            ' o  n -m -l  0  0  0  0 '.split()
        ]
        # but I think it should be this because otherwise:
        # lk = kl and that cannot be right???
        expect = [
            ' 1  i  j  k  l  m  n  o '.split(), 
            ' i -1  k -j  m -l  o -n '.split(),
            ' j -k -1  i  n -o -l  m '.split(),
            ' k  j -i -1  o  n -m -l '.split(),
            ' l -m -n -o  0  0  0  0 '.split(),
            ' m  l  o -n  0  0  0  0 '.split(),
            ' n -o  l  m  0  0  0  0 '.split(),
            ' o  n -m  l  0  0  0  0 '.split()
        ]
        calc = generate_str(object)
        if DEBUG or VERBOSE:
            print("\nWarning: this class is known to be incorrect. Do not use in production")
        if DEBUG: print("\ncalc:   %s\nexpect: %s" % (calc, expect))
        if VERBOSE:
            imaginaries = '1ijklmnopqrstuvw'
            df = pd.DataFrame (expect, index = list(imaginaries[0:8]), columns = list(imaginaries[0:8]))
            print ("\n=== Expected %s table ===\n%s" % (object.__name__,df))
            df2 = pd.DataFrame (calc, index = list(imaginaries[0:8]), columns = list(imaginaries[0:8]))
            print ("\n=== Calculated %s table ===\n%s" % (object.__name__,df2))
        if DEBUG or VERBOSE: print('...')
        self.assertListEqual(calc, expect)


class TestHyperbolicQuaternion(unittest.TestCase):

    def test_unit_multiplication (self):
        "Hyperbolic Quaternion unit product table"
        object = HyperbolicQuaternion
        # ref: https://en.wikipedia.org/wiki/Hyperbolic_quaternion
        expectEDfromWIKI = [
            ' 1  i  j  k '.split(), 
            ' i  1  k -j '.split(),
            ' j -k  1  i '.split(),
            ' k  j -i  1 '.split()
        ]
        expect = [
            ' 1  i  j  k '.split(), 
            ' i  1  k  j '.split(),
            ' j -k  1 -i '.split(),
            ' k -j  i -1 '.split()
        ]
        calc = generate_str(object)
        if DEBUG or VERBOSE:
            print("\nWarning: this class is known to be incorrect. Do not use in production")
        if DEBUG: print("\ncalc:   %s\nexpect: %s" % (calc, expect))
        if VERBOSE:
            imaginaries = '1ijklmnopqrstuvw'
            df = pd.DataFrame (expect, index = list(imaginaries[0:4]), columns = list(imaginaries[0:4]))
            print ("\n=== Expected %s table ===\n%s" % (object.__name__,df))
            df2 = pd.DataFrame (calc, index = list(imaginaries[0:4]), columns = list(imaginaries[0:4]))
            print ("\n=== Calculated %s table ===\n%s" % (object.__name__,df2))
        if DEBUG or VERBOSE: print('...')
        self.assertListEqual(calc, expect)


class TestCommutative(unittest.TestCase):
    """Commutative Product tests"""

    def _is_commutative(self, obj, LOOPS=100):
        d = dim(obj)
        if VERBOSE or DEBUG: 
            print(LOOPS, 'loops')
        for n in range(LOOPS):
            x = random_vector(obj)
            y = random_vector(obj)
            if d > 2:
                self.assertNotEqual(x*y,y*x)
            else:
                self.assertEqual(x*y,y*x)
        if DEBUG: 
            print('x: %r, y: %r' % (x,y))

    def test_complex_commutative(self):
        self._is_commutative(Complex)

    def test_quaternion(self):
        self._is_commutative(Quaternion)

    def test_octonion(self):
        self._is_commutative(Octonion)

    def test_sedenion(self):
        self._is_commutative(Sedenion)

    def test_cd32(self):
        self._is_commutative(Cd32)


class TestWeakAlternativeCondition(unittest.TestCase):

    def _is_weak_alternative(self, object, LOOPS=100):
        "Weak Alternative Condition tests"
        d = dim(object)
        if VERBOSE or DEBUG: 
            print(LOOPS, 'loops')
        for n in range(LOOPS):
            x = random_vector(object)
            y = random_vector(object)
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
            print("\n")
            print("%12s : %s\n" % ('          x ',  x         ))
            print("%12s : %s\n" % ('          y ',  y         ))
            print("%12s : %s\n" % ('(y * x) * x ', (y* x) *x  ))
            print("%12s : %s\n" % (' y *(x  * x ',  y*(x  *x) ))
            print("%12s : %s\n" % ('(x * y) * x ', (x* y) *x  ))
            print("%12s : %s\n" % (' x *(y  * x ',  x*(y  *x) ))
            print("%12s : %s\n" % ('(x * x) * y ', (x* x) *y  ))
            print("%12s : %s\n" % (' x *(x  * y)',  x*(x  *y) ))
            print("%12s : %s\n" % ('(x * y) * y ', (x* y) *y  ))
            print("%12s : %s\n" % (' x *(y  * y)',  x*(y  *y) ))
            print("%12s : %s\n" % ('(y * x) * y ', (y* x) *y  ))
            print("%12s : %s\n" % (' y *(x  * y)',  y*(x  *y) ))
            print("%12s : %s\n" % ('(y * y) * x ', (y* y) *x  ))
            print("%12s : %s\n" % (' y *(y  * x)',  y*(y  *x) ))

    def test_complex(self):
        self._is_weak_alternative(Complex)

    def test_quaternion(self):
        self._is_weak_alternative(Quaternion)

    def test_octonion(self):
        self._is_weak_alternative(Octonion)

    def test_sedenion(self):
        self._is_weak_alternative(Sedenion)

    def test_32ion(self):
        self._is_weak_alternative(Cd32)


class TestDiophantusIdentity(unittest.TestCase):
    """
    Diophantus identity test
    Brahmagupta–Fibonacci / Diophantus identity
    """

    def _is_diophantus_identity(self, object, LOOPS=100):
        "Diophantus identity test"
        PRECISION = 10**-9
        d = dim(object)
        if VERBOSE or DEBUG: 
            print(LOOPS, 'loops')
        for n in range(LOOPS):
            x = random_vector(object)
            y = random_vector(object)
            if d > 8:
                self.assertNotAlmostEqual(abs(x) * abs(y), abs(x*y), delta=PRECISION)
            else:
                self.assertAlmostEqual(abs(x) * abs(y), abs(x*y), delta=PRECISION)
        if VERBOSE:
            print("\n")
            print("%8s: %s\n" % (            'x',  x            ))
            print("%8s: %s\n" % (            'y',  y            ))
            print("%8s: %s\n" % (        'abs(x)',abs(x)        ))
            print("%8s: %s\n" % (        'abs(y)',abs(y)        ))
            print("%8s: %s\n" % ( 'abs(x)*abs(y)',abs(x)*abs(y) ))
            print("%8s: %s\n" % (      'abs(x*y)',abs(x*y)      ))


    def test_complex(self):
        self._is_diophantus_identity(Complex)

    def test_quaternion(self):
        self._is_diophantus_identity(Quaternion)

    def test_octonion(self):
        self._is_diophantus_identity(Octonion)

    def test_sedenion(self):
        self._is_diophantus_identity(Sedenion)

    def test_32ion(self):
        self._is_diophantus_identity(Cd32)


class TestMoufangCondition(unittest.TestCase):
    """Moufang condition tests"""

    def _is_moufang_condition(self, object, LOOPS=100):
        "Moufang condition tests"
        dm = dim(object)
        if VERBOSE or DEBUG: 
            print(LOOPS, 'loops')
        for n in range(LOOPS):
            x = random_vector(object)
            y = random_vector(object)
            z = random_vector(object)
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
            print("%8s: %s\n" % ('a',a))
            print("%8s: %s\n" % ('b',b))
            print("%8s: %s\n" % ('c',c))
            print("%8s: %s\n" % ('d',d))
            print("%8s: %s\n" % ('e',e))
            print("%8s: %s\n" % ('f',f))
            print("%8s: %s\n" % ('g',g))
            print("%8s: %s\n" % ('h',h))
            print("%8s: %s\n" % ('z*(x*(z*y))',a))
            print("%8s: %s\n" % ('((z*x)*z)*y',b))
            print("%8s: %s\n" % ('x*(z*(y*z))',c))
            print("%8s: %s\n" % ('((x*z)*y)*z',d))
            print("%8s: %s\n" % ('(z*x)*(y*z)',e))
            print("%8s: %s\n" % ('(z*(x*y))*z',f))
            print("%8s: %s\n" % ('(z*x)*(y*z)',g))
            print("%8s: %s\n" % ('z*((x*y)*z)',h))

    def test_complex(self):
        self._is_moufang_condition(Complex)

    def test_quaternion(self):
        self._is_moufang_condition(Quaternion)

    def test_octonion(self):
        self._is_moufang_condition(Octonion)

    def test_sedenion(self):
        self._is_moufang_condition(Sedenion)

    def atest_32ion(self):
        self._is_moufang_condition(Cd32)


class TestPowerAssociative(unittest.TestCase):
    """Power Associative tests"""

    def _is_power_associative(self, object, LOOPS=100):
        PRECISION = 10**-9
        d = dim(object)
        for n in range(LOOPS):
            x = random_vector(object).normalize()
            y = random_vector(object).normalize()
            z = x * y
            self.assertAlmostEqual(abs(x),1,delta=PRECISION)
            self.assertAlmostEqual(abs(y),1,delta=PRECISION)
            if d > 8:
                self.assertAlmostEqual(abs(z),1,delta=10**-1)
            else:
                self.assertAlmostEqual(abs(z),1,delta=PRECISION)
        if VERBOSE:
            print("\n")
            print("%8s: %s\n" % ('x',  x))
            print("%8s: %s\n" % ('y',  y))
            print("%8s: %s\n" % ('x*y',z))
            print("%8s: %s\n" % ('abs(x)',abs(x)))
            print("%8s: %s\n" % ('abs(y)',abs(y)))
            print("%8s: %s\n" % ('abs(x*y)',abs(z)))
            print("%8s: %s\n" % ('distance to 1',abs(z)-1))


    def test_complex(self):
        self._is_diophantus_identity(Complex)

    def test_complex(self):
        self._is_power_associative(Complex)

    def test_quaternion(self):
        self._is_power_associative(Quaternion)

    def test_octonion(self):
        self._is_power_associative(Octonion)

    def test_sedenion(self):
        self._is_power_associative(Sedenion)

    def atest_32ion(self):
        self._is_power_associative(Cd32ion)


class TestTwoSquareIdentity(unittest.TestCase):
    """Two square identity"""

    def test_two_square_identity(self):
        object = Complex
        x = random_vector(object)
        y = random_vector(object)
        a, b = x
        c, d = y
        formula = object([a*c - d*b, d*a + b*c])
        if VERBOSE:
            print('           x: ', x)
            print('           y: ', y)
            print('formula(y*x): ', formula)
            print('   code(y*x): ', y*x)
        self.assertEqual(formula, x*y)


class TestFourSquareIdentity(unittest.TestCase):
    """Four square identity"""

    def test_four_square_identity(self):
        object = Quaternion
        x = random_vector(object)
        y = random_vector(object)
        a,b,c,d = x
        e,f,g,h = y
        r = a*e - b*f - c*g - d*h
        s = a*f + b*e + c*h - d*g
        t = a*g - b*h + c*e + d*f
        u = a*h + b*g - c*f + d*e
        formula = object([r,s,t,u])
        if VERBOSE:
            print('           x: ', x)
            print('           y: ', y)
            print('formula(y*x): ', formula)
            print('   code(y*x): ', y*x)
        self.assertEqual(formula, x*y)


class TestEightSquareIdentity(unittest.TestCase):
    """Eight square identity"""

    def test_eight_square_identity(self):
        object = Octonion
        x = random_vector(object)
        y = random_vector(object)
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
        formula = object([z1,z2,z3,z4,z5,z6,z7,z8])
        if VERBOSE:
            print('           x: ', x)
            print('           y: ', y)
            print('formula(x*y): ', formula)
            print('   code(x*y): ', x*y)
        self.assertEqual(formula, x*y)


class TestSixteenSquareIdentity(unittest.TestCase):
    """Sixteen square identity"""

    def test_sixteen_square_identity(self):
        object = Sedenion
        x = random_vector(object)
        y = random_vector(object)
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
        formula = object([z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15,z16])
        if VERBOSE:
            print('           x: ', x)
            print('           y: ', y)
            print('formula(x*y): ', formula)
            print('   code(x*y): ', x*y)
        self.assertEqual(formula, x*y)


if __name__ == '__main__':

    unittest.main()
