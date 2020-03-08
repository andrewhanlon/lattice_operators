#!/usr/bin/env python

from math import isclose
from sympy import KroneckerDelta, Array, S, N, re, im
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
  parser = argparse.ArgumentParser(description="Check if I=3 pi-pi-pi operator is zero")
  parser.add_argument('op_file', metavar='OP_FILE', type=str,
                      help="Specify operator file")

  args = parser.parse_args()

  print(zero(args.op_file))


def zero(op_file):
  f_handler = open(op_file, 'r')
  num_ops = int(f_handler.readline().strip())
  operator = S.Zero
  for op_line in f_handler:
    operator += get_operator(op_line)
  f_handler.close()

  for val in operator.coefficients.values():
    sim_val = N(val)
    re_close = isclose(re(sim_val), 0., rel_tol=1e-09, abs_tol=1e-08)
    im_close = isclose(im(sim_val), 0., rel_tol=1e-09, abs_tol=1e-08)
    if not re_close or not im_close:
      return False

  return True



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
