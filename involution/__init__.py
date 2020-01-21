# -*- coding: utf-8 -*-
# $Id: involution/__init__.py $
# Author: Jeff Anderson <truejeffanderson@gmail.com>
# Copyright: AGPL
"""

involution

python class for Cayley Dickson construction and operation

source: https://github.com/peawormsworth
author: Jeffrey B Anderson - truejeffanderson at gmail.com
"""

import numpy as np
from math import sqrt, log, ceil

class Algebra():
    """
    Generate involution algebras of Cayley Dickson constructions.

    provides a multi-dimensional number calculator for vaious dimensions and algebraic forms and relationships.
    
    see: involution.albegra for common algebraic constructions.
    """

    #ii = None
    #dp = None
    precision = 10 ** -9
    str_func  = None


    def conj (m):
        """conjugate this object"""
        conj = m._conj(m[:])
        return m.__class__(conj, dp=m.dp, ii=m.ii)

        
    # look into np.conjugate() and remove this routine...
    def _conj(m,x):
        """conjugate a given list according to the construct of this object"""
        conj = -x
        conj[0] *= -1
        return conj


    def _curse_mul (m,x,y):
        """recursively multiply the list according to the construction of this object"""
        n = len(x)
        h = n // 2
        if h:
            a,b = x[:h],x[h:]
            c,d = y[:h],y[h:]
            z = np.zeros(n)
            level = int(log(h,2))
            ii = m.ii[level]
            dp = m.dp[level]
            mult = m._curse_mul
            if dp == 'p0':
                z[:h] = mult(c,a) + ii * mult(m._conj(b),d)
                z[h:] = mult(d,m._conj(a)) + mult(b,c)
            if dp == 'p1':
                z[:h] = mult(c,a) + ii * mult(d,m._conj(b))
                z[h:] = mult(m._conj(a),d) + mult(c,b)
            if dp == 'p2':
                z[:h] = mult(a,c) + ii * mult(m._conj(b),d)
                z[h:] = mult(d,m._conj(a)) + mult(b,c)
            if dp == 'p3':
                z[:h] = mult(a,c) + ii * mult(d,m._conj(b))
                z[h:] = mult(m._conj(a),d) + mult(c,b)
            if dp == 'pt0':
                z[:h] = mult(c,a) + ii * mult(b,m._conj(d))
                z[h:] = mult(a,d) + mult(m._conj(c),b)
            if dp == 'pt1':
                z[:h] = mult(c,a) + ii * mult(m._conj(d),b)
                z[h:] = mult(d,a) + mult(b,m._conj(c))
            if dp == 'pt2':
                z[:h] = mult(a,c) + ii * mult(b,m._conj(d))
                z[h:] = mult(a,d) + mult(m._conj(c),b)
            if dp == 'pt3':
                z[:h] = mult(a,c) + ii * mult(m._conj(d),b)
                z[h:] = mult(d,a) + mult(b,m._conj(c))

        else:
            z = x * y
        return z


    def __mul__ (m,o):
        """multiply two objects of similar type"""
        try: 
            o_state = o.state
        except:
            o_state = np.zeros(len(m))
            o_state[0] = o
        return m.__class__(m._curse_mul(m.state,o_state), dp=m.dp, ii=m.ii)


    def __abs__ (m,state=None):
        """absolute value of this object"""
        if state is None:
            return m.__abs__(m.state)
        h = len(state) // 2
        if h:
            a,b = state[:h],state[h:]
            # obtain imaginary squared value based on list size...
            level = ceil(log(h,2))
            ii = m.ii[level]
            return sqrt(m.__abs__(a) ** 2 - ii * m.__abs__(b) ** 2)
        return state[0]


    def __getitem__ (m, index=None):
        """the coefficient of the basis for the provided index"""
        return m.state[index]


    def __add__ (m, z):
        """Addition: z1+z2 = (a,b)+(c,d) = (a+c,b+d) auto called for: a+b"""
        try:
            sum = m.state + z.state
        except:
            sum = m[:].tolist()
            sum[0] = sum[0] + z
        return m.__class__(sum, dp=m.dp, ii=m.ii)

        return int(log(len(m),2))


    def level (m):
        """object level is 1 for list size of two and incrments as list size doubles"""
        return int(log(len(m),2))


    def __len__ (m):
        """object length as a count of dimensions or list size"""
        return len(m.state)


    def __str__ (m):
        """generic string function caller"""
        return m.str_func()

    def str_ijk (m):
        """output the object in a readable form (i,j,k notation)"""
        string = ''
        if len(m) > 16:
            symbols = 'abcdefghijklmnopqrstuvwxyzABCDEF'
        else:
            symbols = ' ijklmnopqrstuvwx'
        for i in range(len(m)):
            try:
                value = abs(m[i])
                if value:
                    sign = ''
                    if m[i] > 0:
                        if len(string):
                            sign = sign + '+'
                    else:
                        sign = sign + '-'
                        pass
                    if value % 1 < m.precision:
                        value = int(value)
                    if value == 1 and i:
                        value = ''
                    if i:
                        symbol = symbols[i]
                    else:
                        symbol = ''
                    string += sign + str(value) + symbol
            except:
                if string:
                    string = string + ',' + str(m[i])
                else:
                    string = str(m[i])
        if string == '':
            string = '0'
        return string


    def __repr__ (m):
        """replicate: output this number as a string suitable for evaluation"""
        return "%r([%s], dp=%r, ii=%r)" % (str(type(m).__name__), 
            ','.join(map(str,m[:])), m.dp, m.ii)


    def __radd__ (m, z):
        """Called for a+b, when a is a number and b is this class"""
        return m + z


    def __rsub__ (m, z):
        """Called for a-b, where a is a number and b is of this class"""
        return -m + z


    def __sub__ (m, z):
        """Subtraction: z1-z2 = (a,b)-(c,d) = (a-c,b-d) auto called for a-b"""
        return m + -z


    def __rmul__ (m, z):
        """Called for a*b, where a is a number and b is of this class"""
        return m * z


    def __rtruediv__ (m, z):
        """Called for a/b, where a is a number and b is of this class"""
        return ~m * z


    def __truediv__ (m, z):
        """Division: z1/z2 = (a,b) × (c,d)⁻¹ = (a,b) × inverse(c,d)"""
        if isinstance(z, Algebra):
            return  ~z * m
        else:
            return 1/z * m


    def __invert__ (m):
        """Invert: z⁻¹. called with ~ object"""
        return m.conj()  * (1/abs(m) ** 2)


    def __neg__ (m):
        """Negate: -z = -1 × z. called automatically with -obj"""
        return -1 * m


    def __pos__ (m):
        """Positive: +z = z"""
        return 1 * m


    def replace (m, z):
        """Replace: the existing coefficients with those of the given one"""
        print('new state: ', z)
        m.state = z.state.copy()
        return m


    def __iadd__ (m, z):
        """Addition with assignment: z += x"""
        return m.replace(m + z)


    def __isub__ (m, z):
        """Subtraction with assignment: z -= x"""
        return m.replace(m / z)


    def __imul__ (m, z):
        """Multiplication with assignment: z *= x"""
        return m.replace(m * z)


    def __idiv__ (m, z):
        """Division with assignment: z /= x"""
        return m.replace(m / z)


    def __eq__ (m, z):
        """Equality condition: true if z = x"""
        return abs(m - z) <= m.precision


    def __ne__ (m, z):
        """Inequality: true is z ≠ x"""
        return not m == z


    def normalize (m):
        """Normalize: z/|z| = zn, where norm of zn = 1"""
        return m / abs(m)


    def __init__ (m, state, dp=None, ii=None, str_func=None):
        """object constructor"""
        try:
            # check list input is an even 2^n
            assert(len(state)), 'input list required'
            assert(log(len(state),2).is_integer()), 'input list must be a power of 2 (list size = 2**n)'
            if dp:
                assert(2 ** len(dp) == len(state)), 'dp list must be 2 to the power of input list size'
            if ii:
                assert(2 ** len(ii) == len(state)), 'ii list must be 2 to the power of input list size'

            m.state = np.asarray(state, dtype = np.float64)
            if dp:
                m.dp = dp
            elif m.dp is None:
                m.dp = ('pt3',) * m.level()

            if ii:
                m.ii = ii
            elif m.ii is None:
                m.ii = (-1,) * m.level()

            if str_func:
                m.str_func = str_func
            elif m.str_func is None:
                m.str_func = m.str_ijk

        except:
            print("unknown exception here")
            raise
        
