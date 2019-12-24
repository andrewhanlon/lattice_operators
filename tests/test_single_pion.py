#!/usr/bin/env python

from math import isclose
from sympy import KroneckerDelta, Array, S

from context import operators

from operators.operators import QuarkField, AntiQuarkField, DiracIdx, ColorIdx, Operator, \
    OperatorRepresentation
from operators.tensors import Gamma
from operators.cubic_rotations import P, P0

g = Gamma()

u = QuarkField.create('u')
dbar = AntiQuarkField.create('d')

a = ColorIdx('a')
b = ColorIdx('b')
i = DiracIdx('i')
j = DiracIdx('j')

pion = KroneckerDelta(a, b)*dbar[a,i]*Array(g.five)[i,j]*u[b,j]


def main():

  pion_op = Operator(pion, P([0,-1,1]))
  op_rep = OperatorRepresentation(pion_op)

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
