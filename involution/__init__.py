# $Id: involution/__init__.py $
# Author: Jeff Anderson <truejeffanderson@gmail.com>
# Copyright: AGPL
"""

involution

python class for Cayley Dickson construction and operation

use:

import involution Algebra, dp

class EightDimField (Albegra):
    dim  = 8             # doubling product function from involution/product.py
    ii   = 1             # imaginary square: 1, 0 or -1
    dp   = dp.PT0        # number of dimensions
    base = FourDimField  # ojbect capable of multiplication such as this one.

field = EightDimField(1,2,3,4,3,2,1,0)

# field will function as every the eight dimensional number of the specified algebra 

see:
    involution.product - 8 doubling products for dp
    involution.algebra - common named coonstructions: Quaternion, Octonion, etc.

source: https://github.com/peawormsworth
author: Jeffrey B Anderson - truejeffanderson at gmail.com
"""

from math import sqrt, log, cos, sin, acos, floor

class Algebra(object):
    """
    Generat involution algebras of Cayley Dickson constructions.

    provides a multi-dimensional number calculator for vaious dimensions and algebraic forms and relationships.
    
    see: involution.albegra for common algebraic constructions.
    """

    precision = 10 ** -9

    def conj (m):
        """
        Conjugate: z* = (a,b)* = (a*,-b)
        It is the negation of imaginary values
        x*=x when x is a real number
        """

        try:    ac = m.a.conj()
        except: ac = m.a
        return  m.__class__(ac, -m.b)

    def __mul__ (m,z):
        return m.__class__(*m.dp(m,z))

    def __getitem__ (m, index=None):
        """the coefficient of the basis for the provided index"""
        try:    return (m.a[:], m.b[:])[index]
        except: return (m.a, m.b)[index]

        return m._flatX()[index]

    def __init__ (m,*ls):
        try:
            a,b = ls
            if (type(a) == type(b) == m.base):
                m.a(a)
                m.b(b)
            else:
                raise TypeError
        except:
            # divide ls in half as inputs to two new base objects.
            # put the resulting objects in a and b...
            # todo: there is a better way to do this...
            if len(ls) == 2:
                a,b = ls
                m.a(m.base(a))
                m.b(m.base(b))
            if len(ls) == 4:
                a,b,c,d = ls
                m.a(m.base(a,b))
                m.b(m.base(c,d))
            if len(ls) == 8:
                a,b,c,d,e,f,g,h = ls
                m.a(m.base(a,b,c,d))
                m.b(m.base(e,f,g,h))
            if len(ls) == 16:
                a,b,c,d,e,f,g,h,i,j,k,l,M,n,o,p = ls
                m.a(m.base(a,b,c,d,e,f,g,h))
                m.b(m.base(i,j,k,l,M,n,o,p))

    def __str__ (m):
        return '(' + ','.join(map(str,m[:])) + ')'

    def __repr__ (m):
        """replicate: output this number as a string suitable for evaluation"""
        return "%r(%r)" % ( str(type(m).__name__), m[:] )

    def __add__ (m, z):
        """Addition: z1+z2 = (a,b)+(c,d) = (a+c,b+d) auto called for: a+b"""
        try:
            a, b = m.a + z.a, m.b + z.b
        except:
            a, b = m.a + z, m.b
        return m.__class__(a, b)

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

    def __abs__ (m):
        """Absolute Value = Norm: √(norm(a)²+norm(b)²). called with abs(obj)"""
        return sqrt(abs(m.a) ** 2 - m.ii * abs(m.b) ** 2)

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
        m.a(z.a)
        m.b(z.b)
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

    def a (m, a=None):
        """a: left half of number"""
        if a is not None:
            m.a = a
        return m.a

    def b (m, b=None):
        """b: right half of number"""
        if b is not None:
            m.b = b
        return m.b


