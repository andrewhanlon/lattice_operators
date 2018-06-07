import sys

from sympy import Eijk
from sympy import get_contraction_structure, get_indices
from sympy import Array
from sympy import sqrt
from sympy import sympify
from sympy import Integer
from sympy import I
from sympy import Matrix

from context import operators

from operators.operators import QuarkField, DiracIdx, ColorIdx, OperatorRepresentation, Operator, OperatorBasis
from operators.cubic_rotations import *
from operators.tensors import Gamma
import operators.grassmann as gr

from pprint import pprint

g = Gamma()

C5p = Array(g.chargeConj * g.five * g.parityPlus)
C1p = Array(g.chargeConj * g.one * g.parityPlus)
C2p = Array(g.chargeConj * g.two * g.parityPlus)
C3p = Array(g.chargeConj * g.three * g.parityPlus)

i = DiracIdx('i')
j = DiracIdx('j')


def main():
  u = QuarkField.create('u')
  d = QuarkField.create('d')
  s = QuarkField.create('s')

  sud_i = Baryon(s,u,d, i)

  lam_lam_ops = flavor_term(sud_i, sud_i)

  test_spin_zero_op(lam_lam_ops[0])

  uud_i = Baryon(u,u,d, i)
  ssd_i = Baryon(s,s,d, i)
  dud_i = Baryon(d,u,d, i)
  ssu_i = Baryon(s,s,u, i)
  ssd_i = Baryon(s,s,d, i)

  uud_ssd = flavor_term(uud_i, ssd_i)
  dud_ssu = flavor_term(dud_i, ssu_i)
  ssd_uud = flavor_term(ssd_i, uud_i)
  ssu_dud = flavor_term(ssu_i, dud_i)

  N_Xi_a_I0_S1 = uud_ssd[1] - dud_ssu[1] - ssd_uud[1] + ssu_dud[1]
  N_Xi_a_I0_S2 = uud_ssd[2] - dud_ssu[2] - ssd_uud[2] + ssu_dud[2]
  N_Xi_a_I0_S3 = uud_ssd[3] - dud_ssu[3] - ssd_uud[3] + ssu_dud[3]

  test_spin_one_ops(N_Xi_a_I0_S1, N_Xi_a_I0_S2, N_Xi_a_I0_S3)

  '''
  uus_i = Baryon(u,u,s, i)
  dds_i = Baryon(d,d,s, i)
  dus_i = Baryon(d,u,s, i)
  uds_i = Baryon(u,d,s, i)

  uus_dds_ops = flavor_term(uus_i, dds_i)
  dus_uds_ops = flavor_term(dus_i, uds_i)
  dus_dus_ops = flavor_term(dus_i, dus_i)
  uds_dus_ops = flavor_term(uds_i, dus_i)
  uds_uds_ops = flavor_term(uds_i, uds_i)
  dds_uus_ops = flavor_term(dds_i, uus_i)

  sig_sig_0 = 2*uus_dds_ops[0] - dus_uds_ops[0] - dus_dus_ops[0] - uds_dus_ops[0] - uds_uds_ops[0] + 2*dds_uus_ops[0]
  print(type(sig_sig_0))
  '''




def test_spin_zero_op(op0):

  ''' Tested
  print("A1+ (p_i^2 = 0)")
  A1p_op = op0.projectMomentum(P0, P0)

  op_rep = OperatorRepresentation(A1p_op)
  print("\t{}\n".format(op_rep.littleGroupContents(True)))
  '''

  ''' Tested
  print("A1+ (p_i^2 = 1)")
  A1p_op = op0.projectMomentum(P([0,0,1]), P([0,0,-1])) \
      + op0.projectMomentum(P([0,1,0]), P([0,-1,0])) \
      + op0.projectMomentum(P([1,0,0]), P([-1,0,0])) \
      + op0.projectMomentum(P([0,0,-1]), P([0,0,1])) \
      + op0.projectMomentum(P([0,-1,0]), P([0,1,0])) \
      + op0.projectMomentum(P([-1,0,0]), P([1,0,0]))

  op_rep = OperatorRepresentation(A1p_op)
  print("\t{}\n".format(op_rep.littleGroupContents(True)))
  '''

  
def test_spin_one_ops(op1, op2, op3):

  print("T1p 1")
  T1p_1_op = op1.projectMomentum(P([0,0,1]), P([0,0,-1])) \
      + op1.projectMomentum(P([0,1,0]), P([0,-1,0])) \
      + op1.projectMomentum(P([1,0,0]), P([-1,0,0])) \
      + op1.projectMomentum(P([0,0,-1]), P([0,0,1])) \
      + op1.projectMomentum(P([0,-1,0]), P([0,1,0])) \
      + op1.projectMomentum(P([-1,0,0]), P([1,0,0]))

  T1p_2_op = op2.projectMomentum(P([0,0,1]), P([0,0,-1])) \
      + op2.projectMomentum(P([0,1,0]), P([0,-1,0])) \
      + op2.projectMomentum(P([1,0,0]), P([-1,0,0])) \
      + op2.projectMomentum(P([0,0,-1]), P([0,0,1])) \
      + op2.projectMomentum(P([0,-1,0]), P([0,1,0])) \
      + op2.projectMomentum(P([-1,0,0]), P([1,0,0]))

  T1p_3_op = op3.projectMomentum(P([0,0,1]), P([0,0,-1])) \
      + op3.projectMomentum(P([0,1,0]), P([0,-1,0])) \
      + op3.projectMomentum(P([1,0,0]), P([-1,0,0])) \
      + op3.projectMomentum(P([0,0,-1]), P([0,0,1])) \
      + op3.projectMomentum(P([0,-1,0]), P([0,1,0])) \
      + op3.projectMomentum(P([-1,0,0]), P([1,0,0]))

  op_rep = OperatorRepresentation(T1p_1_op, T1p_2_op, T1p_3_op)
  print("\t{}\n".format(op_rep.littleGroupContents(True)))


  print("T1p 2")
  T1p_1_op = 3*op_1.projectMomentum(P([1,0,0]), P([-1,0,0])) \
      - op_1.projectMomentum(P([1,0,0]), P([-1,0,0])) \
      - op_1.projectMomentum(P([0,1,0]), P([0,-1,0])) \
      - op_1.projectMomentum(P([0,0,1]), P([0,0,-1]))

  T1p_2_op = 3*op_2.projectMomentum(P([0,1,0]), P([0,-1,0])) \
      - op_2.projectMomentum(P([1,0,0]), P([-1,0,0])) \
      - op_2.projectMomentum(P([0,1,0]), P([0,-1,0])) \
      - op_2.projectMomentum(P([0,0,1]), P([0,0,-1]))

  T1p_3_op = 3*op_3.projectMomentum(P([0,0,1]), P([0,0,-1])) \
      - op_3.projectMomentum(P([1,0,0]), P([-1,0,0])) \
      - op_3.projectMomentum(P([0,1,0]), P([0,-1,0])) \
      - op_3.projectMomentum(P([0,0,1]), P([0,0,-1]))


  op_rep = OperatorRepresentation(T1p_1_op, T1p_2_op, T1p_3_op)
  print("\t{}\n".format(op_rep.littleGroupContents(True)))


def flavor_term(fl1_i, fl2_i):

  fl1_0_j = fl1_i * C5p[i,j]
  fl1_1_j = fl1_i * C1p[i,j]
  fl1_2_j = fl1_i * C2p[i,j]
  fl1_3_j = fl1_i * C3p[i,j]

  fl1_0_0_op = Operator(fl1_0_j.subs(j,0))
  fl1_0_1_op = Operator(fl1_0_j.subs(j,1))
  fl1_0_2_op = Operator(fl1_0_j.subs(j,2))
  fl1_0_3_op = Operator(fl1_0_j.subs(j,3))

  fl1_1_0_op = Operator(fl1_1_j.subs(j,0))
  fl1_1_1_op = Operator(fl1_1_j.subs(j,1))
  fl1_1_2_op = Operator(fl1_1_j.subs(j,2))
  fl1_1_3_op = Operator(fl1_1_j.subs(j,3))

  fl1_2_0_op = Operator(fl1_2_j.subs(j,0))
  fl1_2_1_op = Operator(fl1_2_j.subs(j,1))
  fl1_2_2_op = Operator(fl1_2_j.subs(j,2))
  fl1_2_3_op = Operator(fl1_2_j.subs(j,3))

  fl1_3_0_op = Operator(fl1_3_j.subs(j,0))
  fl1_3_1_op = Operator(fl1_3_j.subs(j,1))
  fl1_3_2_op = Operator(fl1_3_j.subs(j,2))
  fl1_3_3_op = Operator(fl1_3_j.subs(j,3))


  fl2_0_op = Operator(fl2_i.subs(i,0))
  fl2_1_op = Operator(fl2_i.subs(i,1))
  fl2_2_op = Operator(fl2_i.subs(i,2))
  fl2_3_op = Operator(fl2_i.subs(i,3))


  op_0 = fl1_0_0_op*fl2_0_op + fl1_0_1_op*fl2_1_op + fl1_0_2_op*fl2_2_op + fl1_0_3_op*fl2_3_op
  op_1 = fl1_1_0_op*fl2_0_op + fl1_1_1_op*fl2_1_op + fl1_1_2_op*fl2_2_op + fl1_1_3_op*fl2_3_op
  op_2 = fl1_2_0_op*fl2_0_op + fl1_2_1_op*fl2_1_op + fl1_2_2_op*fl2_2_op + fl1_2_3_op*fl2_3_op
  op_3 = fl1_3_0_op*fl2_0_op + fl1_3_1_op*fl2_1_op + fl1_3_2_op*fl2_2_op + fl1_3_3_op*fl2_3_op

  return [op_0, op_1, op_2, op_3]


def Baryon(q1, q2, q3, alpha):
  aB = ColorIdx('aB')
  bB = ColorIdx('bB')
  cB = ColorIdx('cB')

  iB = DiracIdx('iB')
  jB = DiracIdx('jB')

  baryon = Eijk(aB, bB, cB) * q2[aB,iB] * C5p[iB,jB] * q3[bB,jB] * q1[cB,alpha]
  return baryon


if __name__ == "__main__":
  main()
