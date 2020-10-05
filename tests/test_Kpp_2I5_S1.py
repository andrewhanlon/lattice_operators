#!/usr/bin/env python

from math import isclose
from sympy import S
import argparse

import single_hadrons as sh

from context import operators

from operators.operators import Operator, OperatorRepresentation
from operators.cubic_rotations import P, P0


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
  strange_1 = int(op_line[3])
  row_1 = int(op_line[4])
  mom_2 = P([int(op_line[5]), int(op_line[6]), int(op_line[7])])
  strange_2 = int(op_line[8])
  row_2 = int(op_line[9])
  mom_3 = P([int(op_line[10]), int(op_line[11]), int(op_line[12])])
  strange_3 = int(op_line[13])
  row_3 = int(op_line[14])

  if row_1 != 1 or row_2 != 1 or row_3 != 1:
    print("operator with irrep row > 1")
    exit()

  if strange_1 != 1 and strange_2 != 0 and strange_3 != 0:
    print("operator has wrong strangeness")
    exit()


  coeff_re = float(op_line[15])
  coeff_im = float(op_line[16])

  coeff = coeff_re + coeff_im*1j

  op = coeff*Operator(sh.kaon_p, mom_1)*Operator(sh.pion_p, mom_2)*Operator(sh.pion_p, mom_3)

  return op


if __name__ == "__main__":
  main()
