#!/usr/bin/env python

from math import isclose
from sympy import KroneckerDelta, Array, S, simplify, expand
import argparse

from context import operators

from operators.operators import QuarkField, AntiQuarkField, DiracIdx, ColorIdx, Operator, \
    OperatorRepresentation
from operators.tensors import Gamma
from operators.cubic_rotations import P, P0, _POINT_GROUP

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
  parser.add_argument('--operator', metavar='OP_FILES', type=str, nargs='+', action='append',
                      help="Specify operator files")

  args = parser.parse_args()

  check_operator(args.operator)


def check_operator(operator_files):
  op_reps = list()
  for operator_set_files in operator_files:
    ops = list()
    for operator_file in operator_set_files:
      f_handler = open(operator_file, 'r')
      num_ops = int(f_handler.readline().strip())
      operator = S.Zero
      for op_line in f_handler:
        operator += get_operator(op_line)

      f_handler.close()

      ops.append(operator)

    op_rep = OperatorRepresentation(*ops)
    op_reps.append(op_rep)

  for el in _POINT_GROUP:
    mats = list()
    for mom_rep in op_reps:
      try:
        mats.append(mom_rep.getRepresentationMatrix(el).applyfunc(expand).applyfunc(simplify))
      except Exception as e:
        print(f"\tFAILED: {e}\n")
        return

    if not all(mat == mats[0] for mat in mats):
      print("\tFAILED: not all operators transform in the same way")
      print("\t\tELEMENT: {}".format(el))
      for mat in mats:
        print()
        pprint(mat)

      return

  print("PASSED")


def get_operator(op_line):
  op_line = op_line.split()
  mom_1 = P([int(op_line[0]), int(op_line[1]), int(op_line[2])])
  mom_2 = P([int(op_line[5]), int(op_line[6]), int(op_line[7])])

  coeff_re = float(op_line[10])
  coeff_im = float(op_line[11])

  coeff = coeff_re + coeff_im*1j

  op = coeff*Operator(pion, mom_1)*Operator(pion, mom_2)

  return op


if __name__ == "__main__":
  main()
