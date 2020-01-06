# $Id: involution/algebra.py $
# Author: Jeff Anderson <truejeffanderson@gmail.com>
# Copyright: AGPL

"""
involution.algebra

python class of algebriac types of involution.Algebra

source: https://github.com/peawormsworth
author: Jeffrey B Anderson - truejeffanderson at gmail.com
"""

from involution         import Algebra
from involution.product import *


class Complex         (Algebra):
    ii   = -1
    dp   = P3
    dim  = 2
    base = float


class Dual            (Algebra):
    ii   = 0
    dp   = PT3
    dim  = 2
    base = float


class Split           (Algebra):
    ii   = 1
    dp   = PT3
    dim  = 2
    base = float


class Quaternion      (Complex):
    ii   = -1
    dp   = P3
    dim  = 4
    base = Complex


class Octonion        (Quaternion):
    ii   = -1
    dp   = PT0
    dim  = 8
    base = Quaternion


class Sedenion        (Octonion):
    ii   = -1
    dp   = PT3
    dim  = 16
    base = Octonion


class Cd32            (Sedenion):
    ii   = -1
    dp   = PT3
    dim  = 32
    base = Sedenion


class Cd64            (Cd32):
    ii   = -1
    dp   = PT3
    dim  = 64
    base = Cd32


class SplitQuaternion (Complex):
    ii   = 1
    dp   = P0
    dim  = 4
    base = Complex


class DualComplex     (Complex):
    ii   = 0
    dp   = P3
    dim  = 4
    base = Complex


class DualQuaternion  (Quaternion):
    ii   = 0
    dp   = PT1
    dim  = 8
    base = Quaternion


# Remaining classes incomplete or do not match provided tables
# placeholders for testing/development.

# this is not right...
class BiComplex       (Split):
    ii   = -1
    dp   = PT3
    dim  = 4
    base = Split

    # weird conjugate...
    def conj (m):
        return  m.__class__(m.a, -m.b)


# this is not right...
class SplitOctonion   (Quaternion):
    ii   = 1
    dp   = P3
    dim  = 8
    base = Quaternion


class SplitBiQuaternion (Algebra):
    """
    SplitBiQuaternions are also known as:

    elliptic biquaternions, Clifford biquaternion, dyquaternions, 𝔻 ⊗ ℍ where 𝔻 is a dual number, ℍ ⊕ ℍ
    """
    pass


class BiQuaternion         (Algebra): pass
class BiOctonion           (Algebra): pass
class HyperbolicQuaternion (Algebra): pass
class MulticomplexNumber   (Algebra): pass
class Spacetime            (Algebra): pass


