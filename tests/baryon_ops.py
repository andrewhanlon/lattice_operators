from sympy import *
from pprint import pprint
import itertools

from context import operators

from operators.cubic_rotations import *
from operators.cubic_rotations import _POINT_GROUP
from operators.operators import QuarkField, DiracIdx, ColorIdx, OperatorRepresentation, Operator
from operators.gamma import Gamma
import operators.grassmann as gr
from operators.grassmann import perform_contractions

g = Gamma()

C5p = Array(g.chargeConj * g.five * g.parityPlus)
'''
C5 = Array(g.chargeConj * g.five)
pp = Array(g.parityPlus)
'''

u = QuarkField.create('u')
d = QuarkField.create('d')
s = QuarkField.create('s')

i = ColorIdx('i')
j = ColorIdx('j')
k = ColorIdx('k')

a = DiracIdx('a')
b = DiracIdx('b')
c = DiracIdx('c')

baryon_uds = Eijk(i,j,k) * u[i,a] * C5p[a,b] * d[j,b] * s[k,c]
baryon_uds_0 = Operator(baryon_uds.subs(c,0))
baryon_uds_1 = Operator(baryon_uds.subs(c,1))
baryon_uds_2 = Operator(baryon_uds.subs(c,2))
baryon_uds_3 = Operator(baryon_uds.subs(c,3))
baryon_uds_ops = [baryon_uds_0, baryon_uds_1, baryon_uds_2, baryon_uds_3]
bar_uds_op_rep = OperatorRepresentation(*baryon_uds_ops)

baryon_uud = Eijk(i,j,k) * u[i,a] * C5p[a,b] * d[j,b] * u[k,c]
baryon_uud_0 = Operator(baryon_uud.subs(c,0))
baryon_uud_1 = Operator(baryon_uud.subs(c,1))
baryon_uud_2 = Operator(baryon_uud.subs(c,2))
baryon_uud_3 = Operator(baryon_uud.subs(c,3))
baryon_uud_ops = [baryon_uud_0, baryon_uud_1, baryon_uud_2, baryon_uud_3]
bar_uud_op_rep = OperatorRepresentation(*baryon_uud_ops)



for element in _POINT_GROUP:
  uds_mat = bar_uds_op_rep.getRepresentationMatrix(element, False).applyfunc(simplify)
  uud_mat = bar_uud_op_rep.getRepresentationMatrix(element, False).applyfunc(simplify)
  q_mat = spinor_representation.rotation(element, False)
  if uds_mat != uud_mat:
    print()
    pprint(uds_mat)
    pprint(uud_mat)

  print()
  print()
  print()
  print(element)
  pprint(uds_mat)
  print()
  pprint(q_mat)

