from sympy import *
from pprint import pprint
import sys

from context import operators

from operators.operators import *
from operators.cubic_rotations import *
from operators.tensors import Gamma

g = Gamma()


u = QuarkField.create('u')
d = QuarkField.create('d')
ubar = AntiQuarkField.create('u')
dbar = AntiQuarkField.create('d')


a = ColorIdx('a')
b = ColorIdx('b')
c = ColorIdx('c')
i = DiracIdx('i')
j = DiracIdx('j')
k = DiracIdx('k')


pion = KroneckerDelta(a, b)*dbar[a,i]*Array(g.five)[i,j]*u[b,j]

pion_op_1 = Operator(pion, P0)*Operator(pion, P([0,1,0]))*Operator(pion, P([0,-1,0]))
pion_op_2 = Operator(pion, P0)*Operator(pion, P([1,0,0]))*Operator(pion, P([-1,0,0]))
pion_op_3 = Operator(pion, P0)*Operator(pion, P([0,-1,0]))*Operator(pion, P([0,1,0]))
pion_op_4 = Operator(pion, P0)*Operator(pion, P([-1,0,0]))*Operator(pion, P([1,0,0]))

pion_op_r1 = pion_op_1 - pion_op_2 + pion_op_3 - pion_op_4

pion_op_1 = Operator(pion, P0)*Operator(pion, P([0,0,1]))*Operator(pion, P([0,0,-1]))
pion_op_2 = Operator(pion, P0)*Operator(pion, P([0,1,0]))*Operator(pion, P([0,-1,0]))
pion_op_3 = Operator(pion, P0)*Operator(pion, P([1,0,0]))*Operator(pion, P([-1,0,0]))
pion_op_4 = Operator(pion, P0)*Operator(pion, P([0,0,-1]))*Operator(pion, P([0,0,1]))
pion_op_5 = Operator(pion, P0)*Operator(pion, P([0,-1,0]))*Operator(pion, P([0,1,0]))
pion_op_6 = Operator(pion, P0)*Operator(pion, P([-1,0,0]))*Operator(pion, P([1,0,0]))

pion_op_r2 = -2/sqrt(3)*pion_op_1 + 1/sqrt(3)*pion_op_2 + 1/sqrt(3)*pion_op_3 - 2/sqrt(3)*pion_op_4 + 1/sqrt(3)*pion_op_5 + 1/sqrt(3)*pion_op_6

op_rep = OperatorRepresentation(pion_op_r1, pion_op_r2)

print("\nPion:")
print(op_rep.littleGroupContents(True, True))
