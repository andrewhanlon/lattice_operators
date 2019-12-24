#!/usr/bin/env python

from math import isclose
from sympy import Eijk, Array, S

from context import operators

from operators.operators import QuarkField, AntiQuarkField, DiracIdx, ColorIdx, Operator, \
    OperatorRepresentation
from operators.tensors import Gamma
from operators.cubic_rotations import P, P0

g = Gamma()

u = QuarkField.create('u')
d = QuarkField.create('d')

a = ColorIdx('a')
b = ColorIdx('b')
c = ColorIdx('c')

nucleon_1 = 0.5*Eijk(a,b,c)*u[a,0]*u[b,1]*d[c,0] - Eijk(a,b,c)*u[a,0]*u[b,0]*d[c,1] + 0.5*Eijk(a,b,c)*u[a,1]*u[b,0]*d[c,0]
nucleon_2 = -0.5*Eijk(a,b,c)*u[a,1]*u[b,0]*d[c,1] + Eijk(a,b,c)*u[a,1]*u[b,1]*d[c,0] - 0.5*Eijk(a,b,c)*u[a,0]*u[b,1]*d[c,1]

#pion = KroneckerDelta(a, b)*dbar[a,i]*Array(g.five)[i,j]*u[b,j]


def main():

  nucleon_op_1 = Operator(nucleon_1, P([1,1,1]))
  nucleon_op_2 = Operator(nucleon_2, P([1,1,1]))
  op_rep = OperatorRepresentation(nucleon_op_1, nucleon_op_2)

  rep_contents = op_rep.littleGroupContents(False, False)

  nice_str = ""
  for irrep, occurences in rep_contents.items():
    if isclose(occurences, 1., rel_tol=1e-09, abs_tol=1e-08):
      nice_str += "{} + ".format(irrep)
    elif not isclose(occurences, 0., rel_tol=1e-09, abs_tol=1e-08):
      nice_str += "{} {} + ".format(occurences, irrep)

  nice_str = nice_str[:-3]
  print(nice_str)



if __name__ == "__main__":
  main()
