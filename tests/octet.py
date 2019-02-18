'''
import sys
sys.path.insert(0, '/home/ahanlon/git/sympy')
'''

from sympy import *
from pprint import pprint
import itertools
import time
import cProfile

from context import operators

from operators.gamma import Gamma
from operators.quarks import DiracIdx
from operators.octets import OctetBaryon, OctetBaryonElement
from operators.operators import rotate_to, coefficients, project_momentum, momentum
from operators.cubic_rotations import *
from operators.cubic_rotations import _POINT_GROUP
from operators.operator_rep import OperatorRepresentation
from operators.grassmann import GrassmannSymbol


g = Gamma()
C5p = Array(g.chargeConj * g.five * g.parityPlus)
C1p = Array(g.chargeConj * g.one * g.parityPlus)
C2p = Array(g.chargeConj * g.two * g.parityPlus)
C3p = Array(g.chargeConj * g.three * g.parityPlus)

a = DiracIdx('a')
b = DiracIdx('b')

L = OctetBaryon('L', P0)
L2 = OctetBaryon('L', P(1,0,0))
#L2 = OctetBaryon('L', P0)
print(L2(0).args)
'''
op = L[a] * C5p[a,b] * L[b]
op_p = project_momentum(op, [1,0,0], [-1,0,0])
print(op_p.args)
print()
op_p += project_momentum(op, [0,0,1], [0,0,-1])
print(op_p)
print()
'''

exit()
op = project_momentum(op, [0,0,0], [0,0,0])

op_rep = OperatorRepresentation(op)
lg_occs = op_rep.lgIrrepOccurences()
print(lg_occs)
exit()

p = OctetBaryon('p')
n = OctetBaryon('n')
xi_m = OctetBaryon('X-')
xi_0 = OctetBaryon('X0')


op_1 = p[a]*C1p[a,b]*xi_m[b] - n[a]*C1p[a,b]*xi_0[b] - xi_m[a]*C1p[a,b]*p[b] + xi_0[a]*C1p[a,b]*n[b]
op_2 = p[a]*C2p[a,b]*xi_m[b] - n[a]*C2p[a,b]*xi_0[b] - xi_m[a]*C2p[a,b]*p[b] + xi_0[a]*C2p[a,b]*n[b]
op_3 = p[a]*C3p[a,b]*xi_m[b] - n[a]*C3p[a,b]*xi_0[b] - xi_m[a]*C3p[a,b]*p[b] + xi_0[a]*C3p[a,b]*n[b]

op_rep = OperatorRepresentation(op_1, op_2, op_3)
lg_occs = op_rep.lgIrrepOccurences()

'''
pr = cProfile.Profile()
pr.enable()
lg_occs = op_rep.lgIrrepOccurences()
pr.disable()
pr.dump_stats('profile.log')
print(lg_occs)
'''
