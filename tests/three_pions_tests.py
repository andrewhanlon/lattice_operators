from sympy import *
import os
from math import isclose

from collections import namedtuple

from context import operators

from operators.operators import *
from operators.cubic_rotations import *
from operators.tensors import Gamma

base_dir = "/home/ahanlon/research/"
ops_file = "mmmOps/README.selOps"

operator = namedtuple('operator', ['opfiles', 'atrest', 'irrep'])

g = Gamma()

u = QuarkField.create('u')
dbar = AntiQuarkField.create('d')

a = ColorIdx('a')
b = ColorIdx('b')
i = DiracIdx('i')
j = DiracIdx('j')

pion = KroneckerDelta(a, b)*dbar[a,i]*Array(g.five)[i,j]*u[b,j]

def main():

  ops_fullfile = os.path.join(base_dir, ops_file)
  f_handler = open(ops_fullfile, 'r')
  op_files = list()
  for op_file in f_handler:
    op_fullfile = os.path.join(base_dir, op_file.strip())
    tokens = op_fullfile.split('/')
    irrep, irreprow = tokens[-2].split('_')
    irrep = irrep[:-1]
    mom_ray = tokens[-3]
    atrest = False
    if mom_ray == "mom_ray_000":
      atrest = True

    if irreprow == "1":
      op = operator([op_fullfile], atrest, irrep)
      op_files.append(op)
    else:
      op_files[-1].opfiles.append(op_fullfile)

  for op_file_list in op_files:
    check_operator(op_file_list)

  f_handler.close()


def check_operator(op):
  ops = list()
  for op_file in op.opfiles:
    f_handler = open(op_file, 'r')
    num_ops = int(f_handler.readline().strip())
    operator = S.Zero
    for op_line in f_handler:
      operator += get_operator(op_line)

    ops.append(operator)

  print(f"TESTING: {op.opfiles}")
  op_rep = OperatorRepresentation(*ops)
  rep_contents = op_rep.littleGroupContents(False, op.atrest)

  nice_str = ""
  for irrep, occurences in rep_contents.items():
    if isclose(occurences, 1., rel_tol=1e-09, abs_tol=1e-08):
      nice_str += "{} + ".format(irrep)
    elif not isclose(occurences, 0., rel_tol=1e-09, abs_tol=1e-08):
      nice_str += "{} {} + ".format(occurences, irrep)

  nice_str = nice_str[:-3]

  if nice_str != op.irrep:
    print(f"FAILED: Transforms as {nice_str}!")
  else:
    print("PASSED")

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
