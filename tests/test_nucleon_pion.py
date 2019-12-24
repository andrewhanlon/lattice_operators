#!/usr/bin/env python

from math import isclose
from sympy import KroneckerDelta, Eijk, Array, S
from sympy import im, re, N
import argparse

from context import operators

from operators.operators import QuarkField, AntiQuarkField, DiracIdx, ColorIdx, Operator, \
    OperatorRepresentation
from operators.tensors import Gamma
from operators.cubic_rotations import P, P0

g = Gamma()

u = QuarkField.create('u')
d = QuarkField.create('d')
dbar = AntiQuarkField.create('d')

a = ColorIdx('a')
b = ColorIdx('b')
c = ColorIdx('c')
i = DiracIdx('i')
j = DiracIdx('j')

pion = KroneckerDelta(a, b)*dbar[a,i]*Array(g.five)[i,j]*u[b,j]

C5 = Array(g.chargeConj * g.five)

nucleon_1 = 0.5*Eijk(a,b,c)*u[a,0]*u[b,1]*d[c,0] - Eijk(a,b,c)*u[a,0]*u[b,0]*d[c,1] + 0.5*Eijk(a,b,c)*u[a,1]*u[b,0]*d[c,0]
nucleon_2 = -0.5*Eijk(a,b,c)*u[a,1]*u[b,0]*d[c,1] + Eijk(a,b,c)*u[a,1]*u[b,1]*d[c,0] - 0.5*Eijk(a,b,c)*u[a,0]*u[b,1]*d[c,1]
nucleon = [nucleon_1, nucleon_2]

def main():
  parser = argparse.ArgumentParser(description="Check a pi-N operator")
  parser.add_argument('op_files', metavar='OP_FILES', type=str, nargs='+',
                      help="Specify operator files")

  args = parser.parse_args()

  check_operator(args.op_files)


def check_operator(op_files):
  ops = list()
  for op_file in op_files:
    f_handler = open(op_file, 'r')
    num_ops = int(f_handler.readline().strip())
    operator = S.Zero
    for op_line in f_handler:
      operator += get_operator(op_line)

    ops.append(operator)

  op_rep = OperatorRepresentation(*ops)
  use_generators = op_rep.momentum == P0
  rep_contents = op_rep.littleGroupContents(False, use_generators)

  nice_str = ""
  for irrep, occurences in rep_contents.items():
    if isclose(im(occurences), 0., rel_tol=1e-09, abs_tol=1e-08):
      occurences = N(re(occurences))

    if isclose(occurences, 1., rel_tol=1e-09, abs_tol=1e-08):
      nice_str += "{} + ".format(irrep)
    elif not isclose(occurences, 0., rel_tol=1e-09, abs_tol=1e-08):
      nice_str += "{} {} + ".format(occurences, irrep)

  nice_str = nice_str[:-3]
  print(nice_str)


def get_operator(op_line):
  op_line = op_line.split()
  pi_mom = P([int(op_line[0]), int(op_line[1]), int(op_line[2])])
  n_mom = P([int(op_line[4]), int(op_line[5]), int(op_line[6])])

  n_row = int(op_line[7])

  coeff_re = float(op_line[8])
  coeff_im = float(op_line[9])

  coeff = coeff_re + coeff_im*1j

  op = coeff*Operator(pion, pi_mom)*Operator(nucleon[n_row-1], n_mom)

  return op


if __name__ == "__main__":
  main()
