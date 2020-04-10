#!/usr/bin/env python

import argparse

from sympy import S

from context import operators

import isospin_states
from operators.cubic_rotations import P


operator_map = {
    'pi_I3': isospin_states.isoseptet_pion_pion_pion,
    'pi_I0': isospin_states.isosinglet_pion_pion_pion,
    'K_3I2': isospin_states.isoquartet_kaon_kaon_kaon,
}

particle_num_map = {
    'pi_I3': 3,
    'pi_I0': 3,
    'K_3I2': 3,
}


def main():
  parser = argparse.ArgumentParser(description="Isospin project an operator")
  parser.add_argument('--op-file', type=str, required=True,
                      help="Specify operator files")
  parser.add_argument('--proj', type=str, required=True,
                      help="What to project to")
             
  args = parser.parse_args()

  operator = get_operator(args.op_file, args.proj)

  if operator.is_nearly_zero():
    print("zero")
  else:
    print("not zero")


def get_operator(op_file, operator_string):
  f_handler = open(op_file, 'r')
  num_ops = int(f_handler.readline().strip())
  operator = S.Zero
  for op_line in f_handler:
    operator += get_elemental(op_line, operator_string)

  f_handler.close()

  return operator

def get_elemental(op_line, operator_string):
  num_particles = particle_num_map[operator_string]

  op_line = op_line.split()
  moms = list()
  for particle_num in range(num_particles):
    col = 5*particle_num
    moms.append(P([int(op_line[col]), int(op_line[col+1]), int(op_line[col+2])]))

  coeff_re = float(op_line[5*num_particles])
  coeff_im = float(op_line[5*num_particles+1])

  coeff = coeff_re + coeff_im*1j

  op = coeff*operator_map[operator_string](*moms)

  return op


if __name__ == "__main__":
  main()
