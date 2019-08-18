#!/usr/bin/env python

from math import isclose
from sympy import KroneckerDelta, Array, S
import argparse

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
  parser = argparse.ArgumentParser(description="Check a pi-pi-pi operator")
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
    if isclose(occurences, 1., rel_tol=1e-09, abs_tol=1e-08):
      nice_str += "{} + ".format(irrep)
    elif not isclose(occurences, 0., rel_tol=1e-09, abs_tol=1e-08):
      nice_str += "{} {} + ".format(occurences, irrep)

  nice_str = nice_str[:-3]
  print(nice_str)


def get_operator(op_line):
  op_line = op_line.split()
  mom_1 = P([int(op_line[0]), int(op_line[1]), int(op_line[2])])
  mom_2 = P([int(op_line[5]), int(op_line[6]), int(op_line[7])])
  mom_3 = P([int(op_line[10]), int(op_line[11]), int(op_line[12])])

  coeff_re = float(op_line[15])
  coeff_im = float(op_line[16])

  coeff = coeff_re + coeff_im*1j

  op = coeff*Operator(pion, mom_1)*Operator(pion, mom_2)*Operator(pion, mom_3)

  return op


if __name__ == "__main__":
  main()
