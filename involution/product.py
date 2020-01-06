# $Id: involution/product.py $
# Author: Jeff Anderson <truejeffanderson@gmail.com>
# Copyright: AGPL
"""

involution.product

Cayley Dickson doubling products for involution.Algebra

name selection based on: https://arxiv.org/abs/1707.07318
author: Jeffrey B Anderson - truejeffanderson at gmail.com
"""

def _seperate(z):
    try:
        a = z.a
        b = z.b
    except: 
        a = z
        b = 0
    def conjugate (a):
        try:    return a.conj()
        except: return a
    return a, conjugate(a), b, conjugate(b)

           
def P0 (m,z1,z2):
    a,ac,b,bc = _seperate(z1)
    c,cc,d,dc = _seperate(z2)
    return c*a + z1.ii * bc*d, d*ac  + b*c


def P1 (m,z1,z2):
    a,ac,b,bc = _seperate(z1)
    c,cc,d,dc = _seperate(z2)
    return c*a + z1.ii * d*bc, ac*d  + c*b


def P2 (m,z1,z2):
    a,ac,b,bc = __seperate(z1)
    c,cc,d,dc = __seperate(z2)
    return a*c + z1.ii * bc*d, d*ac  + b*c


def P3 (m,z1,z2):
    a,ac,b,bc = _seperate(z1)
    c,cc,d,dc = _seperate(z2)
    return a*c + z1.ii * d*bc, ac*d  + c*b


def PT0 (m,z1,z2):
    a,ac,b,bc = _seperate(z1)
    c,cc,d,dc = _seperate(z2)
    return c*a + z1.ii * b*dc, a*d  + cc*b


def PT1 (m,z1,z2):
    a,ac,b,bc = _seperate(z1)
    c,cc,d,dc = _seperate(z2)
    return c*a + z1.ii * dc*b, d*a  + b*cc


def PT2 (m,z1,z2):
    a,ac,b,bc = _seperate(z1)
    c,cc,d,dc = _seperate(z2)
    return a*c + z1.ii * b*dc, a*d  + cc*b


def PT3 (m,z1,z2):
    a,ac,b,bc = _seperate(z1)
    c,cc,d,dc = _seperate(z2)
    return a*c + z1.ii * dc*b, d*a  + b*cc


