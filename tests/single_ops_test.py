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
dbar = AntiQuarkField.create('d')


a = ColorIdx('a')
b = ColorIdx('b')
c = ColorIdx('c')
i = DiracIdx('i')
j = DiracIdx('j')
k = DiracIdx('k')

Cg5 = Array(g.chargeConj*g.five)

nucleon_1 = Eijk(a,b,c)*u[a,0]*u[b,j]*Cg5[j,k]*d[c,k]
nucleon_2 = Eijk(a,b,c)*u[a,1]*u[b,j]*Cg5[j,k]*d[c,k]

nucleon = [nucleon_1, nucleon_2]

op_rep = OperatorRepresentation(*nucleon)
print(op_rep.littleGroupContents(True, False))

exit()


delta = Eijk(a,b,c)*u[a,i]*u[b,j]*u[c,k]
ops = list()
for i_int in range(4):
  for j_int in range(i_int,4):
    for k_int in range(j_int,4):
      op = delta.subs({i:i_int,j:j_int,k:k_int})
      ops.append(Operator(op, Momentum([0,0,1])))

op_rep = OperatorRepresentation(*ops)

print("\nMoving Delta:")
print(op_rep.littleGroupContents(True, False))

ops_new = list()
for op in ops:
  ops_new.append(op.projectMomentum(P0))

op_rep = OperatorRepresentation(*ops_new)

print("\nAt rest Delta:")
print(op_rep.littleGroupContents(True))

nucleon = Eijk(a,b,c)*(u[a,i]*u[b,j]*d[c,k] - d[a,i]*u[b,j]*u[c,k])
ops = list()
ops.append(Operator(nucleon.subs({i:1, j:0, k:0}), P0))
ops.append(Operator(nucleon.subs({i:1, j:1, k:0}), P0))
ops.append(Operator(nucleon.subs({i:2, j:0, k:0}), P0))
ops.append(Operator(nucleon.subs({i:2, j:0, k:1}), P0))
ops.append(Operator(nucleon.subs({i:2, j:1, k:0}), P0))
ops.append(Operator(nucleon.subs({i:2, j:1, k:1}), P0))
ops.append(Operator(nucleon.subs({i:2, j:2, k:0}), P0))
ops.append(Operator(nucleon.subs({i:2, j:2, k:1}), P0))
ops.append(Operator(nucleon.subs({i:3, j:0, k:0}), P0))
ops.append(Operator(nucleon.subs({i:3, j:0, k:1}), P0))
ops.append(Operator(nucleon.subs({i:3, j:0, k:2}), P0))
ops.append(Operator(nucleon.subs({i:3, j:1, k:0}), P0))
ops.append(Operator(nucleon.subs({i:3, j:1, k:1}), P0))
ops.append(Operator(nucleon.subs({i:3, j:1, k:2}), P0))
ops.append(Operator(nucleon.subs({i:3, j:2, k:0}), P0))
ops.append(Operator(nucleon.subs({i:3, j:2, k:1}), P0))
ops.append(Operator(nucleon.subs({i:3, j:2, k:2}), P0))
ops.append(Operator(nucleon.subs({i:3, j:3, k:0}), P0))
ops.append(Operator(nucleon.subs({i:3, j:3, k:1}), P0))
ops.append(Operator(nucleon.subs({i:3, j:3, k:2}), P0))


op_rep = OperatorRepresentation(*ops)

print("\nNucleon:")
print(op_rep.littleGroupContents(True))

pion = KroneckerDelta(a, b) * dbar[a,i] * u[b,j]
ops = list()
for i_int in range(4):
  for j_int in range(4):
    op = pion.subs({i:i_int,j:j_int})
    ops.append(Operator(op, P0))

op_rep = OperatorRepresentation(*ops)

print("\nPion:")
print(op_rep.littleGroupContents(True))
