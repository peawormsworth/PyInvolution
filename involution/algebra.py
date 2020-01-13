# $Id: involution/algebra.py $
# Author: Jeff Anderson <truejeffanderson@gmail.com>
# Copyright: AGPL

"""
involution.algebra

python class of algebriac types of involution.Algebra

source: https://github.com/peawormsworth
author: Jeffrey B Anderson - truejeffanderson at gmail.com
"""

from involution import Algebra

class Complex (Algebra):
    ii = [ -1  ]
    dp = ['pt3' ]

class Dual (Algebra):
    ii = [  0  ]
    dp = ['pt3']

class Split (Algebra):
    ii = [  1  ]
    dp = ['pt3']

class Quaternion (Algebra):
    ii = [  -1 ] * 2
    dp = ['pt3'] * 2

class Octonion (Algebra):
    ii = [  -1 ] * 3
    dp = ['pt3'] * 3

class Sedenion (Algebra):
    ii = [  -1 ] * 4
    dp = ['pt3'] * 4

class SplitQuaternion (Algebra):
    ii = [  -1 ,   1 ]
    dp = ['pt3','pt3']

class SplitOctonion (Algebra):
    ii = [  -1 ,  -1 ,   1 ]
    dp = ['pt3','pt3', 'p3']

class DualComplex (Algebra):
    ii = [  -1 ,   0 ]
    dp = [ 'p3','pt2']


### Half-tested: I think I am right but a reference says otherwise...

class HyperbolicQuaternion (Algebra):
    ## close...
    ii = [   1 ,    1 ]
    dp = ['pt3',  'pt3']

class DualQuaternion (Algebra):
    ## close...
    ii = [  -1 ,  -1 ,   0 ]
    dp = ['pt3','pt3','pt2']

#########################################
# Untested ...
#########################################

class Cd32 (Algebra):
    ii = [  -1 ] * 5
    dp = ['pt3'] * 5

class Cd64 (Algebra):
    ii = [  -1 ] * 6
    dp = ['pt3'] * 6


class BiComplex            (Algebra): pass
class BiQuaternion         (Algebra): pass
class BiOctonion           (Algebra): pass
class SplitBiQuaternion    (Algebra): pass
class MulticomplexNumber   (Algebra): pass
class Spacetime            (Algebra): pass

