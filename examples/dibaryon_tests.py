from sympy import Eijk
from sympy import Array
from sympy import sqrt
from sympy import S
from sympy import *

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

  print("Forming single baryon operators")
  uud_i = Baryon(u,u,d, i)
  dud_i = Baryon(d,u,d, i)
  #uus_i = Baryon(u,u,s, i)
  #dus_i = Baryon(d,u,s, i)
  #uds_i = Baryon(u,d,s, i)
  #dds_i = Baryon(d,d,s, i)
  sud_i = Baryon(s,u,d, i)
  ssd_i = Baryon(s,s,d, i)
  ssu_i = Baryon(s,s,u, i)
  ssd_i = Baryon(s,s,d, i)

  print("Forming baryon-baryon flavor terms")
  sud_sud = flavor_term(sud_i, sud_i)
  #uus_dds = flavor_term(uus_i, dds_i)
  #dus_uds = flavor_term(dus_i, uds_i)
  #dus_dus = flavor_term(dus_i, dus_i)
  #uds_dus = flavor_term(uds_i, dus_i)
  #uds_uds = flavor_term(uds_i, uds_i)
  #dds_uus = flavor_term(dds_i, uus_i)
  uud_ssd = flavor_term(uud_i, ssd_i)
  dud_ssu = flavor_term(dud_i, ssu_i)
  ssd_uud = flavor_term(ssd_i, uud_i)
  ssu_dud = flavor_term(ssu_i, dud_i)


  print("Forming the lambda-lambda flavor operators")
  L_L_s_I0_S0 = sud_sud[0]
  L_L_s_I0_S1 = sud_sud[1]
  L_L_s_I0_S2 = sud_sud[2]
  L_L_s_I0_S3 = sud_sud[3]
  L_L_s_I0 = [L_L_s_I0_S0, L_L_s_I0_S1, L_L_s_I0_S2, L_L_s_I0_S3]

  print("Forming the nucleon-xi flavor anti-symmetric operators")
  N_X_a_I0_S0 = uud_ssd[0] - dud_ssu[0] - ssd_uud[0] + ssu_dud[0]
  N_X_a_I0_S1 = uud_ssd[1] - dud_ssu[1] - ssd_uud[1] + ssu_dud[1]
  N_X_a_I0_S2 = uud_ssd[2] - dud_ssu[2] - ssd_uud[2] + ssu_dud[2]
  N_X_a_I0_S3 = uud_ssd[3] - dud_ssu[3] - ssd_uud[3] + ssu_dud[3]
  N_X_a_I0 = [N_X_a_I0_S0, N_X_a_I0_S1, N_X_a_I0_S2, N_X_a_I0_S3]

  #test_P0_A1p(L_L_s_I0, 0)  # PASSED
  #test_P0_A1p(L_L_s_I0, 1)  # PASSED
  #test_P0_A1p(L_L_s_I0, 2)  # PASSED
  #test_P0_A1p(L_L_s_I0, 3)  # PASSED

  #test_P0_T1p_Ls(N_X_a_I0, 0)  # PASSED
  #test_P0_T1p_Ls(N_X_a_I0, 1)  # PASSED
  #test_P0_T1p_Ls(N_X_a_I0, 2)  # PASSED
  #test_P0_T1p_Ls(N_X_a_I0, 3)  # PASSED

  #test_P0_T1p_Ld_1(N_X_a_I0)  # PASSED

  #test_P0_T1p_LE_2(N_X_a_I0)  # PASSED

  #test_P0_T1p_LT2_2(N_X_a_I0)  # PASSED

  #test_P1_A1_1_0_S0(L_L_s_I0, P([0,0,1]))  # PASSED
  #test_P1_A1_1_0_S0(L_L_s_I0, P([0,-1,0]))  # PASSED
  #test_P1_A1_1_0_S0(N_X_a_I0, P([0,0,1]))  # PASSED
  #test_P1_A1_1_0_S0(N_X_a_I0, P([0,-1,0]))  # PASSED

  #test_P1_A1_2_1_Ls(L_L_s_I0, P([0,0,1]))  # PASSED
  #test_P1_A1_2_1_Ls(L_L_s_I0, P([0,-1,0]))  # PASSED
  #test_P1_A1_2_1_Ls(N_X_a_I0, P([0,0,1]))  # PASSED
  #test_P1_A1_2_1_Ls(N_X_a_I0, P([0,-1,0]))  # PASSED

  #test_P1_A1_2_1_Lp(L_L_s_I0, P([0,0,1]))  # PASSED
  #test_P1_A1_2_1_Lp(N_X_a_I0, P([0,0,1]))  # PASSED


  #test_P1_A2_1_0_S1(L_L_s_I0, P([0,0,1]))  # PASSED
  #test_P1_A2_1_0_S1(L_L_s_I0, P([-1,0,0]))  # PASSED
  #test_P1_A2_1_0_S1(N_X_a_I0, P([0,0,1]))  # PASSED
  #test_P1_A2_1_0_S1(N_X_a_I0, P([-1,0,0]))  # PASSED

  #test_P1_A2_2_1_Ls(L_L_s_I0, P([0,0,1]))  # PASSED
  #test_P1_A2_2_1_Ls(L_L_s_I0, P([-1,0,0]))  # PASSED
  #test_P1_A2_2_1_Ls(N_X_a_I0, P([0,0,1]))  # PASSED
  #test_P1_A2_2_1_Ls(N_X_a_I0, P([-1,0,0]))  # PASSED

  #test_P1_A2_2_1_Lp(L_L_s_I0, P([0,0,1]))  # PASSED
  #test_P1_A2_2_1_Lp(L_L_s_I0, P([-1,0,0]))  # PASSED
  #test_P1_A2_2_1_Lp(N_X_a_I0, P([0,0,1]))  # PASSED
  #test_P1_A2_2_1_Lp(N_X_a_I0, P([-1,0,0]))  # PASSED

  #test_P1_B1_2_1_S0(L_L_s_I0, P([0,0,1]))  # PASSED B1
  #test_P1_B1_2_1_S0(L_L_s_I0, P([-1,0,0]))  # PASSED B1
  #test_P1_B1_2_1_S0(N_X_a_I0, P([0,0,1]))  # PASSED B1
  #test_P1_B1_2_1_S0(N_X_a_I0, P([-1,0,0]))  # PASSED B1

  #test_P1_B1_2_1_S1(L_L_s_I0, P([0,0,1]))  # PASSED B1
  #test_P1_B1_2_1_S1(L_L_s_I0, P([-1,0,0]))  # PASSED B1
  #test_P1_B1_2_1_S1(N_X_a_I0, P([0,0,1]))  # PASSED B1
  #test_P1_B1_2_1_S1(N_X_a_I0, P([-1,0,0]))  # PASSED B1

  #test_P1_B2_2_1_LB1(L_L_s_I0, P([0,0,1]))  # PASSED B2
  #test_P1_B2_2_1_LB1(L_L_s_I0, P([-1,0,0]))  # PASSED B2
  #test_P1_B2_2_1_LB1(N_X_a_I0, P([0,0,1]))  # PASSED B2
  #test_P1_B2_2_1_LB1(N_X_a_I0, P([-1,0,0]))  # PASSED B2

  #test_P1_B2_2_1_LE2(L_L_s_I0, P([0,0,1]))  # PASSED B2
  #test_P1_B2_2_1_LE2(L_L_s_I0, P([-1,0,0]))  # PASSED B2
  #test_P1_B2_2_1_LE2(N_X_a_I0, P([0,0,1]))  # PASSED B2
  #test_P1_B2_2_1_LE2(N_X_a_I0, P([-1,0,0]))  # PASSED B2

  #test_P1_E2_1_0_S0(L_L_s_I0, P([0,0,1]))  # PASSED
  #test_P1_E2_1_0_S0(L_L_s_I0, P([-1,0,0]))  # PASSED
  #test_P1_E2_1_0_S0(N_X_a_I0, P([0,0,1]))  # PASSED
  #test_P1_E2_1_0_S0(N_X_a_I0, P([-1,0,0]))  # PASSED

  #test_P1_E2_2_1_S0(L_L_s_I0, P([0,0,1]))  # PASSED
  #test_P1_E2_2_1_S0(L_L_s_I0, P([-1,0,0]))  # PASSED
  #test_P1_E2_2_1_S0(N_X_a_I0, P([0,0,1]))  # PASSED
  #test_P1_E2_2_1_S0(N_X_a_I0, P([-1,0,0]))  # PASSED

  #test_P1_E2_2_1_S1_A1(L_L_s_I0, P([0,0,1]))  # PASSED
  #test_P1_E2_2_1_S1_A1(L_L_s_I0, P([-1,0,0]))  # PASSED
  #test_P1_E2_2_1_S1_A1(N_X_a_I0, P([0,0,1]))  # PASSED
  #test_P1_E2_2_1_S1_A1(N_X_a_I0, P([-1,0,0]))  # PASSED

  #test_P1_E2_2_1_S1_B1(L_L_s_I0, P([0,0,1]))  # PASSED
  #test_P1_E2_2_1_S1_B1(L_L_s_I0, P([-1,0,0]))  # PASSED
  #test_P1_E2_2_1_S1_B1(N_X_a_I0, P([0,0,1]))  # PASSED
  #test_P1_E2_2_1_S1_B1(N_X_a_I0, P([-1,0,0]))  # PASSED

  #test_P1_E2_2_1_S1_E2(L_L_s_I0, P([0,0,1]))  # PASSED
  #test_P1_E2_2_1_S1_E2(L_L_s_I0, P([-1,0,0]))  # PASSED
  #test_P1_E2_2_1_S1_E2(N_X_a_I0, P([0,0,1]))  # PASSED
  #test_P1_E2_2_1_S1_E2(N_X_a_I0, P([-1,0,0]))  # PASSED

  #test_P2_A1_2_0_S0(L_L_s_I0, P([0,1,1]))  # PASSED
  #test_P2_A1_2_0_S0(L_L_s_I0, P([-1,1,0]))  # PASSED
  #test_P2_A1_2_0_S0(N_X_a_I0, P([0,1,1]))  # PASSED
  #test_P2_A1_2_0_S0(N_X_a_I0, P([-1,1,0]))  # PASSED

  #test_P2_A1_1_1_S0(L_L_s_I0, P([0,1,0]), P([0,0,1]))  # PASSED
  #test_P2_A1_1_1_S0(L_L_s_I0, P([0,1,0]), P([-1,0,0]))  # PASSED

  #test_P2_A1_1_1_S1(L_L_s_I0, P([0,1,0]), P([0,0,1]))  # PASSED
  #test_P2_A1_1_1_S1(L_L_s_I0, P([0,1,0]), P([-1,0,0]))  # PASSED


  #test_P2_A2_2_0_S1(L_L_s_I0, P([0,1,1]))  # PASSED
  #test_P2_A2_2_0_S1(L_L_s_I0, P([-1,0,1]))  # PASSED
  #test_P2_A2_2_0_S1(N_X_a_I0, P([0,1,1]))  # PASSED
  #test_P2_A2_2_0_S1(N_X_a_I0, P([-1,0,1]))  # PASSED

  #test_P2_A2_1_1_Fs(L_L_s_I0, P([0,0,1]), P([0,1,0]))  # PASSED
  #test_P2_A2_1_1_Fs(L_L_s_I0, P([0,0,1]), P([-1,0,0]))  # PASSED

  #test_P2_A2_1_1_Fa(N_X_a_I0, P([0,0,1]), P([0,1,0]))  # PASSED
  #test_P2_A2_1_1_Fa(N_X_a_I0, P([0,0,1]), P([-1,0,0]))  # PASSED

  #test_P2_B1_2_0_S1(L_L_s_I0, P([0,0,1]), P([0,1,0]))  # PASSED B2
  #test_P2_B1_2_0_S1(L_L_s_I0, P([0,0,1]), P([0,1,0]))  # PASSED B2
  #test_P2_B1_2_0_S1(N_X_a_I0, P([0,0,1]), P([0,1,0]))  # PASSED B2
  #test_P2_B1_2_0_S1(N_X_a_I0, P([0,0,1]), P([0,1,0]))  # PASSED B2

  #test_P2_B1_1_1_S1(N_X_a_I0, P([0,0,1]), P([0,1,0]))  # PASSED B2
  #test_P2_B1_1_1_S1(N_X_a_I0, P([0,0,1]), P([0,-1,0]))  # PASSED B2

  #test_P2_B1_1_1_S0(N_X_a_I0, P([0,0,1]), P([0,1,0]))  # PASSED B2
  #test_P2_B1_1_1_S0(N_X_a_I0, P([0,0,1]), P([0,-1,0]))  # PASSED B2
  

  #test_P2_B2_2_0_S1(L_L_s_I0, P([0,0,1]), P([0,1,0]))  # PASSED B1
  #test_P2_B2_2_0_S1(N_X_a_I0, P([0,0,1]), P([0,1,0]))  # PASSED B1

  #test_P2_B2_1_1_Fs(L_L_s_I0, P([0,0,1]), P([0,1,0]))  # PASSED B1
  #test_P2_B2_1_1_Fs(L_L_s_I0, P([0,0,1]), P([0,-1,0]))  # PASSED B1

  #test_P2_B2_1_1_Fa(N_X_a_I0, P([0,0,1]), P([0,1,0]))  # PASSED B1
  #test_P2_B2_1_1_Fa(N_X_a_I0, P([0,0,1]), P([0,-1,0]))  # PASSED B1

  #test_P3_A1_3_0_S0(L_L_s_I0, P([1,1,1]))  # PASSED
  #test_P3_A1_3_0_S0(L_L_s_I0, P([1,-1,1]))  # PASSED
  #test_P3_A1_3_0_S0(N_X_a_I0, P([1,1,1]))  # PASSED
  #test_P3_A1_3_0_S0(N_X_a_I0, P([1,-1,1]))  # PASSED

  #test_P3_A1_2_1_S0(L_L_s_I0, P([0,1,0]), P([1,0,0]), P([0,0,1]))  # PASSED
  #test_P3_A1_2_1_S0(L_L_s_I0, P([-1,0,0]), P([0,0,-1]), P([0,1,0]))  # PASSED
  #test_P3_A1_2_1_S0(N_X_a_I0, P([0,1,0]), P([1,0,0]), P([0,0,1]))  # PASSED
  #test_P3_A1_2_1_S0(N_X_a_I0, P([-1,0,0]), P([0,0,-1]), P([0,1,0]))  # PASSED

  #test_P3_A1_2_1_S1(L_L_s_I0, P([1,1,1]))  # PASSED
  #test_P3_A1_2_1_S1(L_L_s_I0, P([-1,1,-1]))  # PASSED
  #test_P3_A1_2_1_S1(N_X_a_I0, P([1,1,1]))  # PASSED
  #test_P3_A1_2_1_S1(N_X_a_I0, P([-1,1,-1]))  # PASSED

  #test_P3_A2_3_0_S1(L_L_s_I0, P([1,1,1]))  # PASSED
  #test_P3_A2_3_0_S1(L_L_s_I0, P([-1,1,1]))  # PASSED
  #test_P3_A2_3_0_S1(N_X_a_I0, P([1,1,1]))  # PASSED
  #test_P3_A2_3_0_S1(N_X_a_I0, P([-1,1,1])) # PASSED

  #test_P3_A2_2_1_Ls(L_L_s_I0, P([1,1,1]))  # PASSED
  #test_P3_A2_2_1_Ls(L_L_s_I0, P([-1,1,1]))  # PASSED
  #test_P3_A2_2_1_Ls(N_X_a_I0, P([1,1,1]))  # PASSED
  #test_P3_A2_2_1_Ls(N_X_a_I0, P([-1,1,1]))  # PASSED

  #test_P3_A2_2_1_Lp(L_L_s_I0, P([0,1,0]), P([1,0,0]), P([0,0,1]))  # PASSED
  #test_P3_A2_2_1_Lp(L_L_s_I0, P([-1,0,0]), P([0,0,-1]), P([0,1,0]))  # PASSED
  #test_P3_A2_2_1_Lp(N_X_a_I0, P([0,1,0]), P([1,0,0]), P([0,0,1]))  # PASSED
  #test_P3_A2_2_1_Lp(N_X_a_I0, P([-1,0,0]), P([0,0,-1]), P([0,1,0]))  # PASSED

  #test_P4_A1_Fs(L_L_s_I0, P([0,0,2]))  # PASSED
  #test_P4_A1_Fs(L_L_s_I0, P([2,0,0]))  # PASSED

  #test_P4_A2_Fa(N_X_a_I0, P([0,0,2]))  # PASSED
  #test_P4_A2_Fa(N_X_a_I0, P([2,0,0]))  # PASSED

  #test_P4_E1_Fa(N_X_a_I0, P([0,0,2]))  # PASSED
  #test_P4_E1_Fa(N_X_a_I0, P([2,0,0]))  # PASSED


  #test_equiv_P1_E2_1_0_S0(L_L_s_I0)
  #test_equiv_P1_A1_1_0_S0(L_L_s_I0)


# P=000 tests
def test_P0_A1p(ops, n=0):
  print("P=000 A1+, (pi^2 = {})".format(n))

  if n == 0:
    A1p_op = ops[0].projectMomentum(P0, P0)

  elif n == 1:
    A1p_op = ops[0].projectMomentum(P([0,0,1]), P([0,0,-1])) \
        + ops[0].projectMomentum(P([0,1,0]), P([0,-1,0])) \
        + ops[0].projectMomentum(P([1,0,0]), P([-1,0,0])) \
        + ops[0].projectMomentum(P([0,0,-1]), P([0,0,1])) \
        + ops[0].projectMomentum(P([0,-1,0]), P([0,1,0])) \
        + ops[0].projectMomentum(P([-1,0,0]), P([1,0,0]))

  elif n == 2:
    A1p_op = ops[0].projectMomentum(P([0,1,1]), P([0,-1,-1])) \
        + ops[0].projectMomentum(P([0,-1,-1]), P([0,1,1])) \
        + ops[0].projectMomentum(P([0,1,-1]), P([0,-1,1])) \
        + ops[0].projectMomentum(P([0,-1,1]), P([0,1,-1])) \
        + ops[0].projectMomentum(P([1,0,1]), P([-1,0,-1])) \
        + ops[0].projectMomentum(P([-1,0,-1]), P([1,0,1])) \
        + ops[0].projectMomentum(P([1,0,-1]), P([-1,0,1])) \
        + ops[0].projectMomentum(P([-1,0,1]), P([1,0,-1])) \
        + ops[0].projectMomentum(P([1,1,0]), P([-1,-1,0])) \
        + ops[0].projectMomentum(P([-1,-1,0]), P([1,1,0])) \
        + ops[0].projectMomentum(P([1,-1,0]), P([-1,1,0])) \
        + ops[0].projectMomentum(P([-1,1,0]), P([1,-1,0]))

  elif n == 3:
    A1p_op = ops[0].projectMomentum(P([1,1,1]), P([-1,-1,-1])) \
        + ops[0].projectMomentum(P([1,1,-1]), P([-1,-1,1])) \
        + ops[0].projectMomentum(P([1,-1,1]), P([-1,1,-1])) \
        + ops[0].projectMomentum(P([1,-1,-1]), P([-1,1,1])) \
        + ops[0].projectMomentum(P([-1,1,1]), P([1,-1,-1])) \
        + ops[0].projectMomentum(P([-1,1,-1]), P([1,-1,1])) \
        + ops[0].projectMomentum(P([-1,-1,1]), P([1,1,-1])) \
        + ops[0].projectMomentum(P([-1,-1,-1]), P([1,1,1]))

  else:
    print("pi^2 = {} not implemented".format(n))
    return

  op_rep = OperatorRepresentation(A1p_op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, True)))


def test_P0_T1p_Ls(ops, n=0):
  print("P=000 T1+ (S-wave), (pi^2 = {})".format(n))

  if n == 0:
    T1p_1_op = ops[1].projectMomentum(P([0,0,0]), P([0,0,0]))
    T1p_2_op = ops[2].projectMomentum(P([0,0,0]), P([0,0,0]))
    T1p_3_op = ops[3].projectMomentum(P([0,0,0]), P([0,0,0]))

  elif n == 1:
    T1p_1_op = ops[1].projectMomentum(P([0,0,1]), P([0,0,-1])) \
        + ops[1].projectMomentum(P([0,1,0]), P([0,-1,0])) \
        + ops[1].projectMomentum(P([1,0,0]), P([-1,0,0])) \
        + ops[1].projectMomentum(P([0,0,-1]), P([0,0,1])) \
        + ops[1].projectMomentum(P([0,-1,0]), P([0,1,0])) \
        + ops[1].projectMomentum(P([-1,0,0]), P([1,0,0]))

    T1p_2_op = ops[2].projectMomentum(P([0,0,1]), P([0,0,-1])) \
        + ops[2].projectMomentum(P([0,1,0]), P([0,-1,0])) \
        + ops[2].projectMomentum(P([1,0,0]), P([-1,0,0])) \
        + ops[2].projectMomentum(P([0,0,-1]), P([0,0,1])) \
        + ops[2].projectMomentum(P([0,-1,0]), P([0,1,0])) \
        + ops[2].projectMomentum(P([-1,0,0]), P([1,0,0]))

    T1p_3_op = ops[3].projectMomentum(P([0,0,1]), P([0,0,-1])) \
        + ops[3].projectMomentum(P([0,1,0]), P([0,-1,0])) \
        + ops[3].projectMomentum(P([1,0,0]), P([-1,0,0])) \
        + ops[3].projectMomentum(P([0,0,-1]), P([0,0,1])) \
        + ops[3].projectMomentum(P([0,-1,0]), P([0,1,0])) \
        + ops[3].projectMomentum(P([-1,0,0]), P([1,0,0]))

  elif n == 2:
    T1p_1_op = ops[1].projectMomentum(P([0,1,1]), P([0,-1,-1])) \
        + ops[1].projectMomentum(P([0,-1,-1]), P([0,1,1])) \
        + ops[1].projectMomentum(P([0,1,-1]), P([0,-1,1])) \
        + ops[1].projectMomentum(P([0,-1,1]), P([0,1,-1])) \
        + ops[1].projectMomentum(P([1,0,1]), P([-1,0,-1])) \
        + ops[1].projectMomentum(P([-1,0,-1]), P([1,0,1])) \
        + ops[1].projectMomentum(P([1,0,-1]), P([-1,0,1])) \
        + ops[1].projectMomentum(P([-1,0,1]), P([1,0,-1])) \
        + ops[1].projectMomentum(P([1,1,0]), P([-1,-1,0])) \
        + ops[1].projectMomentum(P([-1,-1,0]), P([1,1,0])) \
        + ops[1].projectMomentum(P([1,-1,0]), P([-1,1,0])) \
        + ops[1].projectMomentum(P([-1,1,0]), P([1,-1,0]))

    T1p_2_op = ops[2].projectMomentum(P([0,1,1]), P([0,-1,-1])) \
        + ops[2].projectMomentum(P([0,-1,-1]), P([0,1,1])) \
        + ops[2].projectMomentum(P([0,1,-1]), P([0,-1,1])) \
        + ops[2].projectMomentum(P([0,-1,1]), P([0,1,-1])) \
        + ops[2].projectMomentum(P([1,0,1]), P([-1,0,-1])) \
        + ops[2].projectMomentum(P([-1,0,-1]), P([1,0,1])) \
        + ops[2].projectMomentum(P([1,0,-1]), P([-1,0,1])) \
        + ops[2].projectMomentum(P([-1,0,1]), P([1,0,-1])) \
        + ops[2].projectMomentum(P([1,1,0]), P([-1,-1,0])) \
        + ops[2].projectMomentum(P([-1,-1,0]), P([1,1,0])) \
        + ops[2].projectMomentum(P([1,-1,0]), P([-1,1,0])) \
        + ops[2].projectMomentum(P([-1,1,0]), P([1,-1,0]))

    T1p_3_op = ops[3].projectMomentum(P([0,1,1]), P([0,-1,-1])) \
        + ops[3].projectMomentum(P([0,-1,-1]), P([0,1,1])) \
        + ops[3].projectMomentum(P([0,1,-1]), P([0,-1,1])) \
        + ops[3].projectMomentum(P([0,-1,1]), P([0,1,-1])) \
        + ops[3].projectMomentum(P([1,0,1]), P([-1,0,-1])) \
        + ops[3].projectMomentum(P([-1,0,-1]), P([1,0,1])) \
        + ops[3].projectMomentum(P([1,0,-1]), P([-1,0,1])) \
        + ops[3].projectMomentum(P([-1,0,1]), P([1,0,-1])) \
        + ops[3].projectMomentum(P([1,1,0]), P([-1,-1,0])) \
        + ops[3].projectMomentum(P([-1,-1,0]), P([1,1,0])) \
        + ops[3].projectMomentum(P([1,-1,0]), P([-1,1,0])) \
        + ops[3].projectMomentum(P([-1,1,0]), P([1,-1,0]))

  elif n == 3:
    T1p_1_op = ops[1].projectMomentum(P([1,1,1]), P([-1,-1,-1])) \
        + ops[1].projectMomentum(P([1,1,-1]), P([-1,-1,1])) \
        + ops[1].projectMomentum(P([1,-1,1]), P([-1,1,-1])) \
        + ops[1].projectMomentum(P([1,-1,-1]), P([-1,1,1])) \
        + ops[1].projectMomentum(P([-1,1,1]), P([1,-1,-1])) \
        + ops[1].projectMomentum(P([-1,1,-1]), P([1,-1,1])) \
        + ops[1].projectMomentum(P([-1,-1,1]), P([1,1,-1])) \
        + ops[1].projectMomentum(P([-1,-1,-1]), P([1,1,1]))

    T1p_2_op = ops[2].projectMomentum(P([1,1,1]), P([-1,-1,-1])) \
        + ops[2].projectMomentum(P([1,1,-1]), P([-1,-1,1])) \
        + ops[2].projectMomentum(P([1,-1,1]), P([-1,1,-1])) \
        + ops[2].projectMomentum(P([1,-1,-1]), P([-1,1,1])) \
        + ops[2].projectMomentum(P([-1,1,1]), P([1,-1,-1])) \
        + ops[2].projectMomentum(P([-1,1,-1]), P([1,-1,1])) \
        + ops[2].projectMomentum(P([-1,-1,1]), P([1,1,-1])) \
        + ops[2].projectMomentum(P([-1,-1,-1]), P([1,1,1]))

    T1p_3_op = ops[3].projectMomentum(P([1,1,1]), P([-1,-1,-1])) \
        + ops[3].projectMomentum(P([1,1,-1]), P([-1,-1,1])) \
        + ops[3].projectMomentum(P([1,-1,1]), P([-1,1,-1])) \
        + ops[3].projectMomentum(P([1,-1,-1]), P([-1,1,1])) \
        + ops[3].projectMomentum(P([-1,1,1]), P([1,-1,-1])) \
        + ops[3].projectMomentum(P([-1,1,-1]), P([1,-1,1])) \
        + ops[3].projectMomentum(P([-1,-1,1]), P([1,1,-1])) \
        + ops[3].projectMomentum(P([-1,-1,-1]), P([1,1,1]))

  else:
    print("pi^2 = {} not implemented".format(n))
    return

  op_rep = OperatorRepresentation(T1p_1_op, T1p_2_op, T1p_3_op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, True)))


def test_P0_T1p_Ld_1(ops):
  print("P=000 T1+ (D-wave), (pi^2 = 1)")
  
  T1p_1_op = 3*ops[1].projectMomentum(P([1,0,0]), P([-1,0,0])) \
      - ops[1].projectMomentum(P([1,0,0]), P([-1,0,0])) \
      - ops[1].projectMomentum(P([0,1,0]), P([0,-1,0])) \
      - ops[1].projectMomentum(P([0,0,1]), P([0,0,-1]))

  T1p_2_op = 3*ops[2].projectMomentum(P([0,1,0]), P([0,-1,0])) \
      - ops[2].projectMomentum(P([1,0,0]), P([-1,0,0])) \
      - ops[2].projectMomentum(P([0,1,0]), P([0,-1,0])) \
      - ops[2].projectMomentum(P([0,0,1]), P([0,0,-1]))

  T1p_3_op = 3*ops[3].projectMomentum(P([0,0,1]), P([0,0,-1])) \
      - ops[3].projectMomentum(P([1,0,0]), P([-1,0,0])) \
      - ops[3].projectMomentum(P([0,1,0]), P([0,-1,0])) \
      - ops[3].projectMomentum(P([0,0,1]), P([0,0,-1]))


  op_rep = OperatorRepresentation(T1p_1_op, T1p_2_op, T1p_3_op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, True)))


def test_P0_T1p_LE_2(ops):
  print("P=000 T1+ (E), (pi^2 = 2)")

  T1p_1_op = 3*(ops[1].projectMomentum(P([0,1,1]), P([0,-1,-1])) \
      + ops[1].projectMomentum(P([0,1,-1]), P([0,-1,1]))) \
      - (ops[1].projectMomentum(P([1,1,0]), P([-1,-1,0])) \
         + ops[1].projectMomentum(P([1,-1,0]), P([-1,1,0]))) \
      - (ops[1].projectMomentum(P([1,0,1]), P([-1,0,-1])) \
         + ops[1].projectMomentum(P([1,0,-1]), P([-1,0,1]))) \
      - (ops[1].projectMomentum(P([0,1,1]), P([0,-1,-1])) \
         + ops[1].projectMomentum(P([0,1,-1]), P([0,-1,1])))

  T1p_2_op = 3*(ops[2].projectMomentum(P([1,0,1]), P([-1,0,-1])) \
      + ops[2].projectMomentum(P([1,0,-1]), P([-1,0,1]))) \
      - (ops[2].projectMomentum(P([1,1,0]), P([-1,-1,0])) \
         + ops[2].projectMomentum(P([1,-1,0]), P([-1,1,0]))) \
      - (ops[2].projectMomentum(P([1,0,1]), P([-1,0,-1])) \
         + ops[2].projectMomentum(P([1,0,-1]), P([-1,0,1]))) \
      - (ops[2].projectMomentum(P([0,1,1]), P([0,-1,-1])) \
         + ops[2].projectMomentum(P([0,1,-1]), P([0,-1,1])))

  T1p_3_op = 3*(ops[3].projectMomentum(P([1,1,0]), P([-1,-1,0])) \
      + ops[3].projectMomentum(P([-1,1,0]), P([1,-1,0]))) \
      - (ops[3].projectMomentum(P([1,1,0]), P([-1,-1,0])) \
         + ops[3].projectMomentum(P([1,-1,0]), P([-1,1,0]))) \
      - (ops[3].projectMomentum(P([1,0,1]), P([-1,0,-1])) \
         + ops[3].projectMomentum(P([1,0,-1]), P([-1,0,1]))) \
      - (ops[3].projectMomentum(P([0,1,1]), P([0,-1,-1])) \
         + ops[3].projectMomentum(P([0,1,-1]), P([0,-1,1])))

  op_rep = OperatorRepresentation(T1p_1_op, T1p_2_op, T1p_3_op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, True)))


def test_P0_T1p_LT2_2(ops):
  print("P=000 T1+ (T2), (pi^2 = 2)")

  T1p_1_op = ops[2].projectMomentum(P([1,1,0]), P([-1,-1,0])) \
      - ops[2].projectMomentum(P([1,-1,0]), P([-1,1,0])) \
      + ops[3].projectMomentum(P([1,0,1]), P([-1,0,-1])) \
      - ops[3].projectMomentum(P([1,0,-1]), P([-1,0,1]))

  T1p_2_op = ops[1].projectMomentum(P([1,1,0]), P([-1,-1,0])) \
      - ops[1].projectMomentum(P([-1,1,0]), P([1,-1,0])) \
      + ops[3].projectMomentum(P([0,1,1]), P([0,-1,-1])) \
      - ops[3].projectMomentum(P([0,1,-1]), P([0,-1,1]))

  T1p_3_op = ops[1].projectMomentum(P([1,0,1]), P([-1,0,-1])) \
      - ops[1].projectMomentum(P([-1,0,1]), P([1,0,-1])) \
      + ops[2].projectMomentum(P([0,1,1]), P([0,-1,-1])) \
      - ops[2].projectMomentum(P([0,-1,1]), P([0,1,-1]))

  op_rep = OperatorRepresentation(T1p_1_op, T1p_2_op, T1p_3_op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, True)))


# P^2=1 A1 tests
def test_P1_A1_1_0_S0(ops, Ptot):

  if Ptot.psq != 1:
    print("Total P != 1! Skipping test")
    return

  print("P^2 = 1 A1 S0 (p1^2 = 1, p2^2 = 0)")
  op = ops[0].projectMomentum(Ptot, P0)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))
  
def test_P1_A1_2_1_Ls(ops, Ptot):

  if Ptot.psq != 1:
    print("Total P != 1! Skipping test")
    return

  print("P^2 = 1 A1 S-wave (p1^2 = 2, p2^2 = 1)")
  op = S.Zero
  if Ptot.x == 0:
    px = P([1,0,0])
    op += ops[0].projectMomentum(Ptot+px, -px) + ops[0].projectMomentum(Ptot-px, px)

  if Ptot.y == 0:
    py = P([0,1,0])
    op += ops[0].projectMomentum(Ptot+py, -py) + ops[0].projectMomentum(Ptot-py, py)

  if Ptot.z == 0:
    pz = P([0,0,1])
    op += ops[0].projectMomentum(Ptot+pz, -pz) + ops[0].projectMomentum(Ptot-pz, pz)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))


def test_P1_A1_2_1_Lp(ops, Ptot):

  if Ptot.psq != 1:
    print("Total P != 1! Skipping test")
    return

  print("P^2 = 1 A1 P-wave (p1^2 = 2, p2^2 = 1)")
  op = S.Zero
  if Ptot.x == 0:
    p = P([1,0,0]) * Ptot
    op += ops[1].projectMomentum(Ptot+p, -p) - ops[1].projectMomentum(Ptot-p, p)

  if Ptot.y == 0:
    p = P([0,1,0]) * Ptot
    op += ops[2].projectMomentum(Ptot+p, -p) - ops[2].projectMomentum(Ptot-p, p)

  if Ptot.z == 0:
    p = P([0,0,1]) * Ptot
    op += ops[3].projectMomentum(Ptot+p, -p) - ops[3].projectMomentum(Ptot-p, p)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))


# P=001 A2 tests
def test_P1_A2_1_0_S1(ops, Ptot):
  if Ptot.psq != 1:
    print("Total P != 1! Skipping test")
    return

  print("P^2 = 1 A2 S1 (p1^2 = 1, p2^2 = 0)")
  if Ptot.x != 0:
    op = ops[1].projectMomentum(Ptot, P0)
  elif Ptot.y != 0:
    op = ops[2].projectMomentum(Ptot, P0)
  elif Ptot.z != 0:
    op = ops[3].projectMomentum(Ptot, P0)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))

def test_P1_A2_2_1_Ls(ops, Ptot):
  if Ptot.psq != 1:
    print("Total P != 1! Skipping test")
    return

  px = P([1,0,0])
  py = P([0,1,0])
  pz = P([0,0,1])

  print("P^2 = 1 A2 (S-wave), (p1^2 = 2, p2^2 = 1)")
  if Ptot.x != 0:
    op = ops[1].projectMomentum(Ptot+py, -py) + ops[1].projectMomentum(Ptot-py, py) \
        + ops[1].projectMomentum(Ptot+pz, -pz) + ops[1].projectMomentum(Ptot-pz, pz)

  elif Ptot.y != 0:
    op = ops[2].projectMomentum(Ptot+px, -px) + ops[2].projectMomentum(Ptot-px, px) \
        + ops[2].projectMomentum(Ptot+pz, -pz) + ops[2].projectMomentum(Ptot-pz, pz)

  elif Ptot.z != 0:
    op = ops[3].projectMomentum(Ptot+px, -px) + ops[3].projectMomentum(Ptot-px, px) \
        + ops[3].projectMomentum(Ptot+py, -py) + ops[3].projectMomentum(Ptot-py, py)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))


def test_P1_A2_2_1_Lp(ops, Ptot):
  if Ptot.psq != 1:
    print("Total P != 1! Skipping test")
    return

  px = P([1,0,0])
  py = P([0,1,0])
  pz = P([0,0,1])

  print("P^2 = 1 A2 (P-wave), (p1^2 = 2, p2^2 = 1)")
  if Ptot.x != 0:
    op = ops[2].projectMomentum(Ptot+py, -py) - ops[2].projectMomentum(Ptot-py, py) \
        + ops[3].projectMomentum(Ptot+pz, -pz) - ops[3].projectMomentum(Ptot-pz, pz)

  elif Ptot.y != 0:
    op = ops[1].projectMomentum(Ptot+px, -px) - ops[1].projectMomentum(Ptot-px, px) \
        + ops[3].projectMomentum(Ptot+pz, -pz) - ops[3].projectMomentum(Ptot-pz, pz)

  elif Ptot.z != 0:
    op = ops[1].projectMomentum(Ptot+px, -px) - ops[1].projectMomentum(Ptot-px, px) \
        + ops[2].projectMomentum(Ptot+py, -py) - ops[2].projectMomentum(Ptot-py, py)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))

def test_P1_B1_2_1_S0(ops, Ptot):
  if Ptot.psq != 1:
    print("Total P != 1! Skipping test")
    return

  if Ptot.x != 0:
    k1 = P([0,1,0])
    k2 = P([0,0,1])

  elif Ptot.y != 0:
    k1 = P([1,0,0])
    k2 = P([0,0,1])

  elif Ptot.z != 0:
    k1 = P([0,1,0])
    k2 = P([1,0,0])


  print("P^2 = 1 B1 (S0), (p1^2 = 2, p2^2 = 1)")

  op = ops[0].projectMomentum(Ptot-k1, k1) + ops[0].projectMomentum(Ptot+k1, -k1) \
      - ops[0].projectMomentum(Ptot-k2, k2) - ops[0].projectMomentum(Ptot+k2, -k2)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))


def test_P1_B1_2_1_S1(ops, Ptot):
  if Ptot.psq != 1:
    print("Total P != 1! Skipping test")
    return

  if Ptot.x != 0:
    k1 = P([0,1,0])
    k2 = P([0,0,1])

    k1s = 2
    k2s = 3

  elif Ptot.y != 0:
    k1 = P([1,0,0])
    k2 = P([0,0,1])

    k1s = 1
    k2s = 3

  elif Ptot.z != 0:
    k1 = P([0,1,0])
    k2 = P([1,0,0])

    k1s = 2
    k2s = 1

  print("P^2 = 1 B1 (S1), (p1^2 = 2, p2^2 = 1)")

  op = ops[k1s].projectMomentum(Ptot + k1*Ptot, -k1*Ptot) - ops[k1s].projectMomentum(Ptot - k1*Ptot, k1*Ptot) \
      - ops[k2s].projectMomentum(Ptot + k2*Ptot, -k2*Ptot) + ops[k2s].projectMomentum(Ptot - k2*Ptot, k2*Ptot)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))


def test_P1_B2_2_1_LB1(ops, Ptot):
  if Ptot.psq != 1:
    print("Total P != 1! Skipping test")
    return

  if Ptot.x != 0:
    k1 = P([0,1,0])
    k2 = P([0,0,1])

    j = 1

  elif Ptot.y != 0:
    k1 = P([1,0,0])
    k2 = P([0,0,1])

    j = 2

  elif Ptot.z != 0:
    k1 = P([0,1,0])
    k2 = P([1,0,0])

    j = 3

  print("P^2 = 1 B2 (L-B1), (p1^2 = 2, p2^2 = 1)")

  op = ops[j].projectMomentum(Ptot - k1, k1) + ops[j].projectMomentum(Ptot + k1, -k1) \
      - ops[j].projectMomentum(Ptot - k2, k2) - ops[j].projectMomentum(Ptot + k2, -k2)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))



def test_P1_B2_2_1_LE2(ops, Ptot):
  if Ptot.psq != 1:
    print("Total P != 1! Skipping test")
    return

  if Ptot.x != 0:
    k1 = P([0,1,0])
    k2 = P([0,0,1])

    k1s = 2
    k2s = 3

  elif Ptot.y != 0:
    k1 = P([1,0,0])
    k2 = P([0,0,1])

    k1s = 1
    k2s = 3

  elif Ptot.z != 0:
    k1 = P([0,1,0])
    k2 = P([1,0,0])

    k1s = 2
    k2s = 1

  print("P^2 = 1 B2 (L-E2), (p1^2 = 2, p2^2 = 1)")

  op = ops[k1s].projectMomentum(Ptot + k1, -k1) - ops[k1s].projectMomentum(Ptot - k1, k1) \
      - ops[k2s].projectMomentum(Ptot + k2, -k2) + ops[k2s].projectMomentum(Ptot - k2, k2)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))



# P=001 E tests
def test_P1_E2_1_0_S0(ops, Ptot):
  if Ptot.psq != 1:
    print("Total P != 1! Skipping test")
    return

  print("P^2 = 1 E2 (p1^2 = 1, p2^2 = 0")

  if Ptot.x != 0:
    op1 = ops[2].projectMomentum(Ptot, P0)
    op2 = ops[3].projectMomentum(Ptot, P0)

  elif Ptot.y != 0:
    op1 = ops[1].projectMomentum(Ptot, P0)
    op2 = ops[3].projectMomentum(Ptot, P0)

  elif Ptot.z != 0:
    op1 = ops[1].projectMomentum(Ptot, P0)
    op2 = ops[2].projectMomentum(Ptot, P0)

  op_rep = OperatorRepresentation(op1, op2)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))


def test_P1_E2_2_1_S0(ops, Ptot):
  if Ptot.psq != 1:
    print("Total P != 1! Skipping test")
    return

  print("P^2 = 1 E2 (S 0) (p1^2 = 2, p2^2 = 1")


  if Ptot.x != 0:
    k1 = P([0,1,0])
    k2 = P([0,0,1])

  elif Ptot.y != 0:
    k1 = P([1,0,0])
    k2 = P([0,0,1])

  elif Ptot.z != 0:
    k1 = P([0,1,0])
    k2 = P([1,0,0])


  op_1 = ops[0].projectMomentum(Ptot + k1*Ptot, -k1*Ptot) - ops[0].projectMomentum(Ptot - k1*Ptot, k1*Ptot)
  op_2 = ops[0].projectMomentum(Ptot + k2*Ptot, -k2*Ptot) - ops[0].projectMomentum(Ptot - k2*Ptot, k2*Ptot)

  op_rep = OperatorRepresentation(op_1, op_2)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))


def test_P1_E2_2_1_S1_A1(ops, Ptot):
  if Ptot.psq != 1:
    print("Total P != 1! Skipping test")
    return

  print("P^2 = 1 E2 (L A1) (p1^2 = 2, p2^2 = 1")


  if Ptot.x != 0:
    jk1_1 = P([1,0,0])
    jk1_2 = P([0,0,1])

    k1s = 2

    jk2_1 = P([1,0,0])
    jk2_2 = P([0,1,0])

    k2s = 3

  elif Ptot.y != 0:
    jk1_1 = P([0,1,0])
    jk1_2 = P([0,0,1])

    k1s = 1

    jk2_1 = P([1,0,0])
    jk2_2 = P([0,1,0])

    k2s = 3

  elif Ptot.z != 0:
    jk1_1 = P([1,0,0])
    jk1_2 = P([0,0,1])

    k1s = 2

    jk2_1 = P([0,0,1])
    jk2_2 = P([0,1,0])

    k2s = 1


  op_1 = ops[k1s].projectMomentum(Ptot + jk1_1, -jk1_1) + ops[k1s].projectMomentum(Ptot - jk1_1, jk1_1) \
      + ops[k1s].projectMomentum(Ptot + jk1_2, -jk1_2) + ops[k1s].projectMomentum(Ptot - jk1_2, jk1_2)

  op_2 = ops[k2s].projectMomentum(Ptot + jk2_1, -jk2_1) + ops[k2s].projectMomentum(Ptot - jk2_1, jk2_1) \
      + ops[k2s].projectMomentum(Ptot + jk2_2, -jk2_2) + ops[k2s].projectMomentum(Ptot - jk2_2, jk2_2)

  op_rep = OperatorRepresentation(op_1, op_2)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))


def test_P1_E2_2_1_S1_B1(ops, Ptot):
  if Ptot.psq != 1:
    print("Total P != 1! Skipping test")
    return

  if Ptot.x != 0:
    k1 = P([0,1,0])
    k2 = P([0,0,1])

    k1s = 2
    k2s = 3

  elif Ptot.y != 0:
    k1 = P([1,0,0])
    k2 = P([0,0,1])

    k1s = 1
    k2s = 3

  elif Ptot.z != 0:
    k1 = P([0,1,0])
    k2 = P([1,0,0])

    k1s = 2
    k2s = 1

  print("P^2 = 1 E2 (L B1), (p1^2 = 2, p2^2 = 1)")

  op_1 = ops[k1s].projectMomentum(Ptot - k1, k1) + ops[k1s].projectMomentum(Ptot + k1, -k1) \
      - ops[k1s].projectMomentum(Ptot - k2, k2) - ops[k1s].projectMomentum(Ptot + k2, -k2)

  op_2 = ops[k2s].projectMomentum(Ptot - k1, k1) + ops[k2s].projectMomentum(Ptot + k1, -k1) \
      - ops[k2s].projectMomentum(Ptot - k2, k2) - ops[k2s].projectMomentum(Ptot + k2, -k2)

  op_rep = OperatorRepresentation(op_1, op_2)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))


def test_P1_E2_2_1_S1_E2(ops, Ptot):
  if Ptot.psq != 1:
    print("Total P != 1! Skipping test")
    return

  print("P^2 = 1 E2 (L E2) (p1^2 = 2, p2^2 = 1")

  if Ptot.x != 0:
    k1 = P([0,1,0])
    k2 = P([0,0,1])

    op1 = ops[1].projectMomentum(Ptot + k1, -k1) - ops[1].projectMomentum(Ptot - k1, k1)
    op2 = -ops[1].projectMomentum(Ptot + k2, -k2) + ops[1].projectMomentum(Ptot - k2, k2)

  elif Ptot.y != 0:
    k1 = P([1,0,0])
    k2 = P([0,0,1])

    op1 = ops[2].projectMomentum(Ptot + k1, -k1) - ops[2].projectMomentum(Ptot - k1, k1)
    op2 = -ops[2].projectMomentum(Ptot + k2, -k2) + ops[2].projectMomentum(Ptot - k2, k2)

  elif Ptot.z != 0:
    k1 = P([0,1,0])
    k2 = P([1,0,0])

    op1 = ops[3].projectMomentum(Ptot + k1, -k1) - ops[3].projectMomentum(Ptot - k1, k1)
    op2 = -ops[3].projectMomentum(Ptot + k2, -k2) + ops[3].projectMomentum(Ptot - k2, k2)

  op_rep = OperatorRepresentation(op1, op2)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))


# P=011 A1 tests
def test_P2_A1_2_0_S0(ops, Ptot):
  if Ptot.psq != 2:
    print("Total P != 2! Skipping test")
    return

  print("P^2 = 2 A1 S0 (p1^2 = 2, p2^2 = 0)")
  op = ops[0].projectMomentum(Ptot, P0)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))


def test_P2_A1_1_1_S0(ops, p1, p2):
  Ptot = p1 + p2
  if Ptot.psq != 2:
    print("Total P != 2! Skipping test")
    return

  print("P^2 = 2 A1 S0 (p1^2 = 1, p2^2 = 1)")
  op = ops[0].projectMomentum(p1, p2)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))


def test_P2_A1_1_1_S1(ops, p1, p2):
  Ptot = p1 + p2
  if Ptot.psq != 2:
    print("Total P != 2! Skipping test")
    return

  print("P^2 = 2 A1 S1 (p1^2 = 1, p2^2 = 1)")

  p = p1*p2
  op = S.Zero
  if p.x != 0:
    op += p.x * ops[1]

  if p.y != 0:
    op += p.y * ops[2]

  if p.z != 0:
    op += p.z * ops[3]

  op = op.projectMomentum(p1, p2)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))

# P=011 A2 tests
def test_P2_A2_2_0_S1(ops, Ptot):
  if Ptot.psq != 2:
    print("Total P != 2! Skipping test")
    return

  print("P^2 = 2 A2 S1 (p1^2 = 2, p2^2 = 0)")

  op = S.Zero
  if Ptot.x != 0:
    op += Ptot.x * ops[1]

  if Ptot.y != 0:
    op += Ptot.y * ops[2]

  if Ptot.z != 0:
    op += Ptot.z * ops[3]

  op = op.projectMomentum(Ptot, P0)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))


def test_P2_A2_1_1_Fs(ops, p1, p2):
  Ptot = p1 + p2
  if Ptot.psq != 2:
    print("Total P != 2! Skipping test")
    return

  print("P^2 = 2 A2 Fs (p1^2 = 1, p2^2 = 1)")

  p = p1 - p2
  op = S.Zero
  if p.x != 0:
    op += p.x * ops[1]

  if p.y != 0:
    op += p.y * ops[2]

  if p.z != 0:
    op += p.z * ops[3]

  op = op.projectMomentum(p1, p2)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))


def test_P2_A2_1_1_Fa(ops, p1, p2):
  Ptot = p1 + p2
  if Ptot.psq != 2:
    print("Total P != 2! Skipping test")
    return

  print("P^2 = 2 A2 Fa (p1^2 = 1, p2^2 = 1)")

  p = p1 + p2
  op = S.Zero
  if p.x != 0:
    op += p.x * ops[1]

  if p.y != 0:
    op += p.y * ops[2]

  if p.z != 0:
    op += p.z * ops[3]

  op = op.projectMomentum(p1, p2)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))

def test_P2_B1_2_0_S1(ops, p1, p2):
  Ptot = p1 + p2
  if Ptot.psq != 2:
    print("Total P != 2! Skipping test")
    return

  print("P^2 = 2 B1 S1 (p1^2 = 2, p2^2 = 0)")

  p = p1*p2
  op = S.Zero
  if p.x != 0:
    op += p.x * ops[1]

  if p.y != 0:
    op += p.y * ops[2]

  if p.z != 0:
    op += p.z * ops[3]

  op = op.projectMomentum(Ptot, P0)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))

def test_P2_B1_1_1_S1(ops, p1, p2):
  Ptot = p1 + p2
  if Ptot.psq != 2:
    print("Total P != 2! Skipping test")
    return

  print("P^2 = 2 B1 S1 (p1^2 = 1, p2^2 = 1)")

  p = p1*p2
  op = S.Zero
  if p.x != 0:
    op += p.x * ops[1]

  if p.y != 0:
    op += p.y * ops[2]

  if p.z != 0:
    op += p.z * ops[3]

  op = op.projectMomentum(p1, p2)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))

def test_P2_B1_1_1_S0(ops, p1, p2):

  Ptot = p1 + p2
  if Ptot.psq != 2:
    print("Total P != 2! Skipping test")
    return

  print("P^2 = 2 B1 S1 (p1^2 = 1, p2^2 = 1)")

  op = ops[0].projectMomentum(p1,p2)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))


def test_P2_B2_2_0_S1(ops, p1, p2):
  Ptot = p1 + p2
  if Ptot.psq != 2:
    print("Total P != 2! Skipping test")
    return

  print("P^2 = 2 B2 S1 (p1^2 = 2, p2^2 = 0)")

  p = p1 - p2
  op = S.Zero
  if p.x != 0:
    op += p.x * ops[1]

  if p.y != 0:
    op += p.y * ops[2]

  if p.z != 0:
    op += p.z * ops[3]

  op = op.projectMomentum(Ptot, P0)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))

def test_P2_B2_1_1_Fs(ops, p1, p2):
  Ptot = p1 + p2
  if Ptot.psq != 2:
    print("Total P != 2! Skipping test")
    return

  print("P^2 = 2 B2 Fs (p1^2 = 1, p2^2 = 1)")

  p = p1 + p2
  op = S.Zero
  if p.x != 0:
    op += p.x * ops[1]

  if p.y != 0:
    op += p.y * ops[2]

  if p.z != 0:
    op += p.z * ops[3]

  op = op.projectMomentum(p1, p2)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))


def test_P2_B2_1_1_Fa(ops, p1, p2):
  Ptot = p1 + p2
  if Ptot.psq != 2:
    print("Total P != 2! Skipping test")
    return

  print("P^2 = 2 B2 Fa (p1^2 = 1, p2^2 = 1)")

  p = p1 - p2
  op = S.Zero
  if p.x != 0:
    op += p.x * ops[1]

  if p.y != 0:
    op += p.y * ops[2]

  if p.z != 0:
    op += p.z * ops[3]

  op = op.projectMomentum(Ptot, P0)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))


def test_P3_A1_3_0_S0(ops, Ptot):
  if Ptot.psq != 3:
    print("Total P != 3! Skipping test")
    return

  print("P^2 = 3 A1 S0 (p1^2 = 3, p2^2 = 0)")
  op = ops[0].projectMomentum(Ptot, P0)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))


def test_P3_A1_2_1_S0(ops, p1, p2, p3):
  Ptot = p1 + p2 + p3
  if Ptot.psq != 3:
    print("Total P != 3! Skipping test")
    return

  print("P^2 = 3 A1 S0 (p1^2 = 2, p2^2 = 1)")

  op = ops[0].projectMomentum(Ptot - p1, p1) \
      + ops[0].projectMomentum(Ptot - p2, p2) \
      + ops[0].projectMomentum(Ptot - p3, p3)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))


def test_P3_A1_2_1_S1(ops, Ptot):
  if Ptot.psq != 3:
    print("Total P != 3! Skipping test")
    return

  p1 = P([Ptot.x,0,0])
  p2 = P([0,Ptot.y,0])
  p3 = P([0,0,Ptot.z])
  print("P^2 = 3 A1 S1 (p1^2 = 2, p2^2 = 1)")
  op = p2.y*ops[2].projectMomentum(Ptot-p3,p3) \
      - p3.z*ops[3].projectMomentum(Ptot-p2,p2) \
      - p1.x*ops[1].projectMomentum(Ptot-p3,p3) \
      + p3.z*ops[3].projectMomentum(Ptot-p1,p1) \
      + p1.x*ops[1].projectMomentum(Ptot-p2,p2) \
      - p2.y*ops[2].projectMomentum(Ptot-p1,p1)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))

def test_P3_A2_3_0_S1(ops, Ptot):
  if Ptot.psq != 3:
    print("Total P != 3! Skipping test")
    return

  op = Ptot.x * ops[1].projectMomentum(Ptot, P0) \
      + Ptot.y * ops[2].projectMomentum(Ptot, P0) \
      + Ptot.z * ops[3].projectMomentum(Ptot, P0)

  print("P^2 = 3 A2 S1 (p1^2 = 3, p2^2 = 0)")

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))

def test_P3_A2_2_1_Ls(ops, Ptot):
  if Ptot.psq != 3:
    print("Total P != 3! Skipping test")
    return

  print("P^2 = 3 A2 Ls (p1^2 = 2, p2^2 = 1)")
  p1 = Momentum([Ptot.x, 0, 0])
  p2 = Momentum([0, Ptot.y, 0])
  p3 = Momentum([0, 0, Ptot.z])

  op = Ptot.x * ops[1].projectMomentum(Ptot - p1, p1) \
      + Ptot.x * ops[1].projectMomentum(Ptot - p2, p2) \
      + Ptot.x * ops[1].projectMomentum(Ptot - p3, p3) \
      + Ptot.y * ops[2].projectMomentum(Ptot - p1, p1) \
      + Ptot.y * ops[2].projectMomentum(Ptot - p2, p2) \
      + Ptot.y * ops[2].projectMomentum(Ptot - p3, p3) \
      + Ptot.z * ops[3].projectMomentum(Ptot - p1, p1) \
      + Ptot.z * ops[3].projectMomentum(Ptot - p2, p2) \
      + Ptot.z * ops[3].projectMomentum(Ptot - p3, p3)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))



def test_P3_A2_2_1_Lp(ops, p1, p2, p3):
  Ptot = p1 + p2 + p3
  if Ptot.psq != 3:
    print("Total P != 3! Skipping test")
    return

  print("P^2 = 3 A2 Lp (p1^2 = 2, p2^2 = 1)")

  op = (p1*3 - p1 - p2 - p3).x * ops[1].projectMomentum(Ptot-p1, p1) \
      + (p1*3 - p1 - p2 - p3).y * ops[2].projectMomentum(Ptot-p1, p1) \
      + (p1*3 - p1 - p2 - p3).z * ops[3].projectMomentum(Ptot-p1, p1) \
      + (p2*3 - p1 - p2 - p3).x * ops[1].projectMomentum(Ptot-p2, p2) \
      + (p2*3 - p1 - p2 - p3).y * ops[2].projectMomentum(Ptot-p2, p2) \
      + (p2*3 - p1 - p2 - p3).z * ops[3].projectMomentum(Ptot-p2, p2) \
      + (p3*3 - p1 - p2 - p3).x * ops[1].projectMomentum(Ptot-p3, p3) \
      + (p3*3 - p1 - p2 - p3).y * ops[2].projectMomentum(Ptot-p3, p3) \
      + (p3*3 - p1 - p2 - p3).z * ops[3].projectMomentum(Ptot-p3, p3)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))

def test_P4_A1_Fs(ops, Ptot):
  if Ptot.psq != 4:
    print("Total P != 4! Skipping test")
    return

  print("P^2 = 4 A1 Fs (p1^2 = 1, p2^2 = 1)")

  op = ops[0].projectMomentum(Ptot/2, Ptot/2)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))

def test_P4_A2_Fa(ops, Ptot):
  if Ptot.psq != 4:
    print("Total P != 4! Skipping test")
    return

  print("P^2 = 4 A2 Fa (p1^2 = 1, p2^2 = 1)")

  if Ptot.x != 0:
    op = Ptot.x * ops[1]

  elif Ptot.y != 0:
    op = Ptot.y * ops[2]

  elif Ptot.z != 0:
    op = Ptot.z * ops[3]

  op = op.projectMomentum(Ptot/2, Ptot/2)

  op_rep = OperatorRepresentation(op)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))

def test_P4_E1_Fa(ops, Ptot):
  if Ptot.psq != 4:
    print("Total P != 4! Skipping test")
    return

  print("P^2 = 4 E1 Fa (p1^2 = 1, p2^2 = 1)")

  if Ptot.x != 0:
    op1 = ops[2].projectMomentum(Ptot/2, Ptot/2)
    op2 = ops[3].projectMomentum(Ptot/2, Ptot/2)

  elif Ptot.y != 0:
    op1 = ops[1].projectMomentum(Ptot/2, Ptot/2)
    op2 = ops[3].projectMomentum(Ptot/2, Ptot/2)

  elif Ptot.z != 0:
    op1 = ops[1].projectMomentum(Ptot/2, Ptot/2)
    op2 = ops[2].projectMomentum(Ptot/2, Ptot/2)

  op_rep = OperatorRepresentation(op1, op2)
  print("\t{}\n".format(op_rep.littleGroupContents(True, False)))

# Equivalent Momentum checks
def test_equiv_P1_E2_1_0_S0(ops):

  print("Testing Equivalent: P^2 = 1 E2 (p1^2 = 1, p2^2 = 0")

  pz = P([0,0,1])
  py = P([0,1,0])

  op1_z = ops[2].projectMomentum(pz, P0)
  op1_y = ops[1].projectMomentum(py, P0)

  op1_z_to_y = op1_z.rotate(C4xi)

  op1_y_coeffs = op1_y.coefficients
  op1_z_to_y_coeffs = op1_z_to_y.coefficients

  for term, coeff in op1_y_coeffs.items():
    print("{} == {}".format(coeff, op1_z_to_y_coeffs[term]))


def test_equiv_P1_A1_1_0_S0(ops):


  print("Testing Equivalent: P^2 = 1 A1 S0 (p1^2 = 1, p2^2 = 0)")

  pz = P([0,0,1])
  py = P([0,1,0])

  op_z = ops[0].projectMomentum(pz, P0)
  op_y = ops[0].projectMomentum(py, P0)

  op_z_to_y = op_z.rotate(C4xi)

  op_y_coeffs = op_y.coefficients
  op_z_coeffs = op_z.coefficients
  op_z_to_y_coeffs = op_z_to_y.coefficients

  print(op_y_coeffs == op_z_to_y_coeffs)





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
