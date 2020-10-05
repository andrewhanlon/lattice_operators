#!/usr/bin/env python

from math import isclose
from sympy import S, im, re
import argparse

import single_hadrons as sh

from context import operators

from operators.operators import Operator, OperatorRepresentation
from operators.cubic_rotations import P, P0


def main():
  parser = argparse.ArgumentParser(description="Check a N-N operator")
  parser.add_argument("--anti", action="store_true", required=False,
                      help="Use anti-symmetric NN")
  parser.add_argument("--spin", action="store_true", required=False,
                      help="Determine spin")
  parser.add_argument("--orib", action="store_true", required=False,
                      help="Determine oribitatl angular momentum")
  parser.add_argument('op_files', metavar='OP_FILES', type=str, nargs='+',
                      help="Specify operator files")

  args = parser.parse_args()

  check_operator(args.op_files, args.anti, args.spin, args.orib)


def check_operator(op_files, anti, spin, orib):
  ops = list()
  for op_file in op_files:
    f_handler = open(op_file, 'r')
    num_ops = int(f_handler.readline().strip())
    operator = S.Zero
    for op_line in f_handler:
      operator += get_operator(op_line, anti, spin, orib)

    ops.append(operator)

  op_rep = OperatorRepresentation(*ops)
  use_generators = op_rep.momentum == P0
  rep_contents = op_rep.littleGroupContents(False, use_generators)

  nice_str = ""
  for irrep, occurences in rep_contents.items():
    if not isclose(im(occurences), 0., rel_tol=1e-09, abs_tol=1e-08):
      print("imaginary occurrent")
      exit()
    occurences = re(occurences)

    if isclose(occurences, 1., rel_tol=1e-09, abs_tol=1e-08):
      nice_str += "{} + ".format(irrep)
    elif not isclose(occurences, 0., rel_tol=1e-09, abs_tol=1e-08):
      nice_str += "{} {} + ".format(occurences, irrep)

  nice_str = nice_str[:-3]
  print(nice_str)


def get_operator(op_line, anti, spin, orib):
  op_line = op_line.split()
  if spin:
    mom_1 = P([0,0,0])
  else:
    mom_1 = P([int(op_line[0]), int(op_line[1]), int(op_line[2])])

  row_1 = int(op_line[3])

  if spin:
    mom_2 = P([0,0,0])
  else:
    mom_2 = P([int(op_line[4]), int(op_line[5]), int(op_line[6])])

  row_2 = int(op_line[7])

  coeff_re = float(op_line[8])
  coeff_im = float(op_line[9])

  coeff = coeff_re + coeff_im*1j

  if anti:
    op = coeff*Operator(sh.nucleon_p[row_1-1], mom_1)*Operator(sh.nucleon_m[row_2-1], mom_2)
  else:
    op = coeff*Operator(sh.nucleon_p[row_1-1], mom_1)*Operator(sh.nucleon_p[row_2-1], mom_2)

  return op


if __name__ == "__main__":
  main()
