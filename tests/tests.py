from sympy import *
from pprint import pprint
import itertools

from context import operators

from operators.cubic_rotations import *
from operators.operators import QuarkField, DiracIdx, ColorIdx, OperatorRepresentation, Operator
from operators.tensors import Gamma
import operators.grassmann as gr
from operators.grassmann import perform_contractions

g = Gamma()

C5p = Array(g.chargeConj * g.five * g.parityPlus)
C1p = Array(g.chargeConj * g.one * g.parityPlus)
C2p = Array(g.chargeConj * g.two * g.parityPlus)
C3p = Array(g.chargeConj * g.three * g.parityPlus)

u = QuarkField.create('u')
d = QuarkField.create('d')
s = QuarkField.create('s')

a1 = ColorIdx('a1')
b1 = ColorIdx('b1')
c1 = ColorIdx('c1')

a2 = ColorIdx('a2')
b2 = ColorIdx('b2')
c2 = ColorIdx('c2')


i1 = DiracIdx('i1')
j1 = DiracIdx('j1')
k1 = DiracIdx('k1')

i2 = DiracIdx('i2')
j2 = DiracIdx('j2')
k2 = DiracIdx('k2')

sud_sud = Eijk(a1,b1,c1) * u[a1,i1] * C5p[i1,j1] * d[b1,j1] * s[c1,k1] * C5p[k1,k2]  * Eijk(a2,b2,c2) * u[a2,i2] * C5p[i2,j2] * d[b2,j2] * s[c2,k2]

op = Operator(sud_sud)
print(op.simplified)

'''
op_rot = op.rotate(C4y)
for term in op_rot.getTerms():
  print(term)
  print('\n')

'''

#op_rep = OperatorRepresentation(op)
#print(op_rep.littleGroupContents(True))

