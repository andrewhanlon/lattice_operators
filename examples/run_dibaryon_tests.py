from context import operators

from operators.operators import QuarkField
from dibaryon_tests import *

def main():
  u = QuarkField.create('u')
  d = QuarkField.create('d')
  s = QuarkField.create('s')

  print("Forming single baryon creation operators")
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

  print("Forming baryon-baryon creation flavor terms")
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


  print("Forming the lambda-lambda flavor creation operators")
  L_L_s_I0_S0 = sud_sud[0]
  L_L_s_I0_S1 = sud_sud[1]
  L_L_s_I0_S2 = sud_sud[2]
  L_L_s_I0_S3 = sud_sud[3]
  L_L_s_I0 = [L_L_s_I0_S0, L_L_s_I0_S1, L_L_s_I0_S2, L_L_s_I0_S3]

  print("Forming the nucleon-xi flavor anti-symmetric creation operators")
  N_X_a_I0_S0 = uud_ssd[0] - dud_ssu[0] - ssd_uud[0] + ssu_dud[0]
  N_X_a_I0_S1 = uud_ssd[1] - dud_ssu[1] - ssd_uud[1] + ssu_dud[1]
  N_X_a_I0_S2 = uud_ssd[2] - dud_ssu[2] - ssd_uud[2] + ssu_dud[2]
  N_X_a_I0_S3 = uud_ssd[3] - dud_ssu[3] - ssd_uud[3] + ssu_dud[3]
  N_X_a_I0 = [N_X_a_I0_S0, N_X_a_I0_S1, N_X_a_I0_S2, N_X_a_I0_S3]


  # Start irrep tests
  passed = 0
  total = 0

  # P = 0 (A1p)
  #passed += test_irrep(P0_A1p(L_L_s_I0, 0), 'A1g', 'P0 A1p (p_i^2 = 0); LL', True); total += 1  # PASSED
  #passed += test_irrep(P0_A1p(L_L_s_I0, 1), 'A1g', 'P0 A1p (p_i^2 = 1); LL', True); total += 1  # PASSED
  #passed += test_irrep(P0_A1p(L_L_s_I0, 2), 'A1g', 'P0 A1p (p_i^2 = 2); LL', True); total += 1  # PASSED
  #passed += test_irrep(P0_A1p(L_L_s_I0, 3), 'A1g', 'P0 A1p (p_i^2 = 3); LL', True); total += 1  # PASSED

  # P = 0 (T1g)
  #passed += test_irrep(P0_T1p_A1p(N_X_a_I0, 0), 'T1g', 'P0 T1p (A1p) (p_i^2 = 0); NXa I=0', True); total += 1  # PASSED
  #passed += test_irrep(P0_T1p_A1p(N_X_a_I0, 1), 'T1g', 'P0 T1p (A1p) (p_i^2 = 1); NXa I=0', True); total += 1  # PASSED
  #passed += test_irrep(P0_T1p_A1p(N_X_a_I0, 2), 'T1g', 'P0 T1p (A1p) (p_i^2 = 2); NXa I=0', True); total += 1  # PASSED
  #passed += test_irrep(P0_T1p_A1p(N_X_a_I0, 3), 'T1g', 'P0 T1p (A1p) (p_i^2 = 3); NXa I=0', True); total += 1  # PASSED

  #passed += test_irrep(P0_T1p_Ep_1(N_X_a_I0), 'T1g', 'P0 T1p (Ep) (p_i^2 = 1); NXa I=0', True); total += 1  # PASSED

  #passed += test_irrep(P0_T1p_Ep_2(N_X_a_I0), 'T1g', 'P0 T1p (Ep) (p_i^2 = 2); NXa I=0', True); total += 1  # PASSED

  #passed += test_irrep(P0_T1p_T2p_2(N_X_a_I0), 'T1g', 'P0 T1p (T2p) (p_i^2 = 2); NXa I=0', True); total += 1  # PASSED


  # Psq = 1 (A1)
  #passed += test_irrep(P1_A1_A1_1_0(L_L_s_I0, P([0,0,1])), 'A1', 'P=(0,0,1) A1 (A1) (1,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_A1_1_0(L_L_s_I0, P([0,0,-1])), 'A1', 'P=(0,0,-1) A1 (A1) (1,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_A1_1_0(L_L_s_I0, P([0,1,0])), 'A1', 'P=(0,1,0) A1 (A1) (1,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_A1_1_0(L_L_s_I0, P([0,-1,0])), 'A1', 'P=(0,-1,0) A1 (A1) (1,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_A1_1_0(L_L_s_I0, P([1,0,0])), 'A1', 'P=(1,0,0) A1 (A1) (1,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_A1_1_0(L_L_s_I0, P([-1,0,0])), 'A1', 'P=(-1,0,0) A1 (A1) (1,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_A1_1_0(N_X_a_I0, P([0,0,1])), 'A1', 'P=(0,0,1) A1 (A1) (1,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_A1_1_0(N_X_a_I0, P([0,0,-1])), 'A1', 'P=(0,0,-1) A1 (A1) (1,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_A1_1_0(N_X_a_I0, P([0,1,0])), 'A1', 'P=(0,1,0) A1 (A1) (1,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_A1_1_0(N_X_a_I0, P([0,-1,0])), 'A1', 'P=(0,-1,0) A1 (A1) (1,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_A1_1_0(N_X_a_I0, P([1,0,0])), 'A1', 'P=(1,0,0) A1 (A1) (1,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_A1_1_0(N_X_a_I0, P([-1,0,0])), 'A1', 'P=(-1,0,0) A1 (A1) (1,0); NXa I=0', False); total += 1  # PASSED

  #passed += test_irrep(P1_A1_A1_2_1(L_L_s_I0, P([0,0,1])), 'A1', 'P=(0,0,1) A1 (A1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_A1_2_1(L_L_s_I0, P([0,0,-1])), 'A1', 'P=(0,0,-1) A1 (A1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_A1_2_1(L_L_s_I0, P([0,1,0])), 'A1', 'P=(0,1,0) A1 (A1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_A1_2_1(L_L_s_I0, P([0,-1,0])), 'A1', 'P=(0,-1,0) A1 (A1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_A1_2_1(L_L_s_I0, P([1,0,0])), 'A1', 'P=(1,0,0) A1 (A1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_A1_2_1(L_L_s_I0, P([-1,0,0])), 'A1', 'P=(-1,0,0) A1 (A1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_A1_2_1(N_X_a_I0, P([0,0,1])), 'A1', 'P=(0,0,1) A1 (A1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_A1_2_1(N_X_a_I0, P([0,0,-1])), 'A1', 'P=(0,0,-1) A1 (A1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_A1_2_1(N_X_a_I0, P([0,1,0])), 'A1', 'P=(0,1,0) A1 (A1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_A1_2_1(N_X_a_I0, P([0,-1,0])), 'A1', 'P=(0,-1,0) A1 (A1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_A1_2_1(N_X_a_I0, P([1,0,0])), 'A1', 'P=(1,0,0) A1 (A1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_A1_2_1(N_X_a_I0, P([-1,0,0])), 'A1', 'P=(-1,0,0) A1 (A1) (2,1); NXa I=0', False); total += 1  # PASSED

  #passed += test_irrep(P1_A1_E2_2_1(L_L_s_I0, P([0,0,1])), 'A1', 'P=(0,0,1) A1 (E2) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_E2_2_1(L_L_s_I0, P([0,0,-1])), 'A1', 'P=(0,0,-1) A1 (E2) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_E2_2_1(L_L_s_I0, P([0,1,0])), 'A1', 'P=(0,1,0) A1 (E2) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_E2_2_1(L_L_s_I0, P([0,-1,0])), 'A1', 'P=(0,-1,0) A1 (E2) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_E2_2_1(L_L_s_I0, P([1,0,0])), 'A1', 'P=(1,0,0) A1 (E2) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_E2_2_1(L_L_s_I0, P([-1,0,0])), 'A1', 'P=(-1,0,0) A1 (E2) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_E2_2_1(N_X_a_I0, P([0,0,1])), 'A1', 'P=(0,0,1) A1 (E2) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_E2_2_1(N_X_a_I0, P([0,0,-1])), 'A1', 'P=(0,0,-1) A1 (E2) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_E2_2_1(N_X_a_I0, P([0,1,0])), 'A1', 'P=(0,1,0) A1 (E2) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_E2_2_1(N_X_a_I0, P([0,-1,0])), 'A1', 'P=(0,-1,0) A1 (E2) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_E2_2_1(N_X_a_I0, P([1,0,0])), 'A1', 'P=(1,0,0) A1 (E2) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A1_E2_2_1(N_X_a_I0, P([-1,0,0])), 'A1', 'P=(-1,0,0) A1 (E2) (2,1); NXa I=0', False); total += 1  # PASSED


  # Psq = 1 (A2)
  #passed += test_irrep(P1_A2_A1_1_0(L_L_s_I0, P([0,0,1])), 'A2', 'P=(0,0,1) A2 (A1) (1,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_A1_1_0(L_L_s_I0, P([0,0,-1])), 'A2', 'P=(0,0,-1) A2 (A1) (1,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_A1_1_0(L_L_s_I0, P([0,1,0])), 'A2', 'P=(0,1,0) A2 (A1) (1,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_A1_1_0(L_L_s_I0, P([0,-1,0])), 'A2', 'P=(0,-1,0) A2 (A1) (1,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_A1_1_0(L_L_s_I0, P([1,0,0])), 'A2', 'P=(1,0,0) A2 (A1) (1,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_A1_1_0(L_L_s_I0, P([-1,0,0])), 'A2', 'P=(-1,0,0) A2 (A1) (1,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_A1_1_0(N_X_a_I0, P([0,0,1])), 'A2', 'P=(0,0,1) A2 (A1) (1,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_A1_1_0(N_X_a_I0, P([0,0,-1])), 'A2', 'P=(0,0,-1) A2 (A1) (1,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_A1_1_0(N_X_a_I0, P([0,1,0])), 'A2', 'P=(0,1,0) A2 (A1) (1,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_A1_1_0(N_X_a_I0, P([0,-1,0])), 'A2', 'P=(0,-1,0) A2 (A1) (1,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_A1_1_0(N_X_a_I0, P([1,0,0])), 'A2', 'P=(1,0,0) A2 (A1) (1,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_A1_1_0(N_X_a_I0, P([-1,0,0])), 'A2', 'P=(-1,0,0) A2 (A1) (1,0); NXa I=0', False); total += 1  # PASSED

  #passed += test_irrep(P1_A2_A1_2_1(L_L_s_I0, P([0,0,1])), 'A2', 'P=(0,0,1) A2 (A1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_A1_2_1(L_L_s_I0, P([0,0,-1])), 'A2', 'P=(0,0,-1) A2 (A1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_A1_2_1(L_L_s_I0, P([0,1,0])), 'A2', 'P=(0,1,0) A2 (A1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_A1_2_1(L_L_s_I0, P([0,-1,0])), 'A2', 'P=(0,-1,0) A2 (A1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_A1_2_1(L_L_s_I0, P([1,0,0])), 'A2', 'P=(1,0,0) A2 (A1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_A1_2_1(L_L_s_I0, P([-1,0,0])), 'A2', 'P=(-1,0,0) A2 (A1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_A1_2_1(N_X_a_I0, P([0,0,1])), 'A2', 'P=(0,0,1) A2 (A1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_A1_2_1(N_X_a_I0, P([0,0,-1])), 'A2', 'P=(0,0,-1) A2 (A1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_A1_2_1(N_X_a_I0, P([0,1,0])), 'A2', 'P=(0,1,0) A2 (A1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_A1_2_1(N_X_a_I0, P([0,-1,0])), 'A2', 'P=(0,-1,0) A2 (A1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_A1_2_1(N_X_a_I0, P([1,0,0])), 'A2', 'P=(1,0,0) A2 (A1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_A1_2_1(N_X_a_I0, P([-1,0,0])), 'A2', 'P=(-1,0,0) A2 (A1) (2,1); NXa I=0', False); total += 1  # PASSED

  #passed += test_irrep(P1_A2_E2_2_1(L_L_s_I0, P([0,0,1])), 'A2', 'P=(0,0,1) A2 (E2) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_E2_2_1(L_L_s_I0, P([0,0,-1])), 'A2', 'P=(0,0,-1) A2 (E2) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_E2_2_1(L_L_s_I0, P([0,1,0])), 'A2', 'P=(0,1,0) A2 (E2) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_E2_2_1(L_L_s_I0, P([0,-1,0])), 'A2', 'P=(0,-1,0) A2 (E2) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_E2_2_1(L_L_s_I0, P([1,0,0])), 'A2', 'P=(1,0,0) A2 (E2) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_E2_2_1(L_L_s_I0, P([-1,0,0])), 'A2', 'P=(-1,0,0) A2 (E2) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_E2_2_1(N_X_a_I0, P([0,0,1])), 'A2', 'P=(0,0,1) A2 (E2) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_E2_2_1(N_X_a_I0, P([0,0,-1])), 'A2', 'P=(0,0,-1) A2 (E2) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_E2_2_1(N_X_a_I0, P([0,1,0])), 'A2', 'P=(0,1,0) A2 (E2) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_E2_2_1(N_X_a_I0, P([0,-1,0])), 'A2', 'P=(0,-1,0) A2 (E2) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_E2_2_1(N_X_a_I0, P([1,0,0])), 'A2', 'P=(1,0,0) A2 (E2) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_A2_E2_2_1(N_X_a_I0, P([-1,0,0])), 'A2', 'P=(-1,0,0) A2 (E2) (2,1); NXa I=0', False); total += 1  # PASSED

  # Psq = 1 (B1)
  #passed += test_irrep(P1_B1_B1_2_1(L_L_s_I0, P([0,0,1])), 'B1', 'P=(0,0,1) B1 (B1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_B1_B1_2_1(L_L_s_I0, P([0,0,-1])), 'B1', 'P=(0,0,-1) B1 (B1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_B1_B1_2_1(L_L_s_I0, P([0,1,0])), 'B1', 'P=(0,1,0) B1 (B1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_B1_B1_2_1(L_L_s_I0, P([0,-1,0])), 'B1', 'P=(0,-1,0) B1 (B1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_B1_B1_2_1(L_L_s_I0, P([1,0,0])), 'B1', 'P=(1,0,0) B1 (B1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_B1_B1_2_1(L_L_s_I0, P([-1,0,0])), 'B1', 'P=(-1,0,0) B1 (B1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_B1_B1_2_1(N_X_a_I0, P([0,0,1])), 'B1', 'P=(0,0,1) B1 (B1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_B1_B1_2_1(N_X_a_I0, P([0,0,-1])), 'B1', 'P=(0,0,-1) B1 (B1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_B1_B1_2_1(N_X_a_I0, P([0,1,0])), 'B1', 'P=(0,1,0) B1 (B1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_B1_B1_2_1(N_X_a_I0, P([0,-1,0])), 'B1', 'P=(0,-1,0) B1 (B1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_B1_B1_2_1(N_X_a_I0, P([1,0,0])), 'B1', 'P=(1,0,0) B1 (B1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_B1_B1_2_1(N_X_a_I0, P([-1,0,0])), 'B1', 'P=(-1,0,0) B1 (B1) (2,1); NXa I=0', False); total += 1  # PASSED

  #passed += test_irrep(P1_B1_E2_2_1(L_L_s_I0, P([0,0,1])), 'B1', 'P=(0,0,1) B1 (E2) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_B1_E2_2_1(L_L_s_I0, P([0,0,-1])), 'B1', 'P=(0,0,-1) B1 (E2) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_B1_E2_2_1(L_L_s_I0, P([0,1,0])), 'B1', 'P=(0,1,0) B1 (E2) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_B1_E2_2_1(L_L_s_I0, P([0,-1,0])), 'B1', 'P=(0,-1,0) B1 (E2) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_B1_E2_2_1(L_L_s_I0, P([1,0,0])), 'B1', 'P=(1,0,0) B1 (E2) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_B1_E2_2_1(L_L_s_I0, P([-1,0,0])), 'B1', 'P=(-1,0,0) B1 (E2) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_B1_E2_2_1(N_X_a_I0, P([0,0,1])), 'B1', 'P=(0,0,1) B1 (E2) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_B1_E2_2_1(N_X_a_I0, P([0,0,-1])), 'B1', 'P=(0,0,-1) B1 (E2) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_B1_E2_2_1(N_X_a_I0, P([0,1,0])), 'B1', 'P=(0,1,0) B1 (E2) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_B1_E2_2_1(N_X_a_I0, P([0,-1,0])), 'B1', 'P=(0,-1,0) B1 (E2) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_B1_E2_2_1(N_X_a_I0, P([1,0,0])), 'B1', 'P=(1,0,0) B1 (E2) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_B1_E2_2_1(N_X_a_I0, P([-1,0,0])), 'B1', 'P=(-1,0,0) B1 (E2) (2,1); NXa I=0', False); total += 1  # PASSED

  # Psq = 1 (B2)
  #passed += test_irrep(P1_B2_B1_2_1(L_L_s_I0, P([0,0,1])), 'B2', 'P=(0,0,1) B2 (B1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_B2_B1_2_1(L_L_s_I0, P([0,0,-1])), 'B2', 'P=(0,0,-1) B2 (B1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_B2_B1_2_1(L_L_s_I0, P([0,1,0])), 'B2', 'P=(0,1,0) B2 (B1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_B2_B1_2_1(L_L_s_I0, P([0,-1,0])), 'B2', 'P=(0,-1,0) B2 (B1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_B2_B1_2_1(L_L_s_I0, P([1,0,0])), 'B2', 'P=(1,0,0) B2 (B1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_B2_B1_2_1(L_L_s_I0, P([-1,0,0])), 'B2', 'P=(-1,0,0) B2 (B1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_B2_B1_2_1(N_X_a_I0, P([0,0,1])), 'B2', 'P=(0,0,1) B2 (B1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_B2_B1_2_1(N_X_a_I0, P([0,0,-1])), 'B2', 'P=(0,0,-1) B2 (B1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_B2_B1_2_1(N_X_a_I0, P([0,1,0])), 'B2', 'P=(0,1,0) B2 (B1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_B2_B1_2_1(N_X_a_I0, P([0,-1,0])), 'B2', 'P=(0,-1,0) B2 (B1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_B2_B1_2_1(N_X_a_I0, P([1,0,0])), 'B2', 'P=(1,0,0) B2 (B1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_B2_B1_2_1(N_X_a_I0, P([-1,0,0])), 'B2', 'P=(-1,0,0) B2 (B1) (2,1); NXa I=0', False); total += 1  # PASSED

  #passed += test_irrep(P1_B2_E2_2_1(L_L_s_I0, P([0,0,1])), 'B2', 'P=(0,0,1) B2 (E2) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_B2_E2_2_1(L_L_s_I0, P([0,0,-1])), 'B2', 'P=(0,0,-1) B2 (E2) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_B2_E2_2_1(L_L_s_I0, P([0,1,0])), 'B2', 'P=(0,1,0) B2 (E2) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_B2_E2_2_1(L_L_s_I0, P([0,-1,0])), 'B2', 'P=(0,-1,0) B2 (E2) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_B2_E2_2_1(L_L_s_I0, P([1,0,0])), 'B2', 'P=(1,0,0) B2 (E2) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_B2_E2_2_1(L_L_s_I0, P([-1,0,0])), 'B2', 'P=(-1,0,0) B2 (E2) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_B2_E2_2_1(N_X_a_I0, P([0,0,1])), 'B2', 'P=(0,0,1) B2 (E2) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_B2_E2_2_1(N_X_a_I0, P([0,0,-1])), 'B2', 'P=(0,0,-1) B2 (E2) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_B2_E2_2_1(N_X_a_I0, P([0,1,0])), 'B2', 'P=(0,1,0) B2 (E2) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_B2_E2_2_1(N_X_a_I0, P([0,-1,0])), 'B2', 'P=(0,-1,0) B2 (E2) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_B2_E2_2_1(N_X_a_I0, P([1,0,0])), 'B2', 'P=(1,0,0) B2 (E2) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_B2_E2_2_1(N_X_a_I0, P([-1,0,0])), 'B2', 'P=(-1,0,0) B2 (E2) (2,1); NXa I=0', False); total += 1  # PASSED

  # Psq = 1 (E2)
  #passed += test_irrep(P1_E2_A1_1_0(L_L_s_I0, P([0,0,1])), 'E', 'P=(0,0,1) E2 (A1) (1,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_A1_1_0(L_L_s_I0, P([0,0,-1])), 'E', 'P=(0,0,-1) E2 (A1) (1,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_A1_1_0(L_L_s_I0, P([0,1,0])), 'E', 'P=(0,1,0) E2 (A1) (1,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_A1_1_0(L_L_s_I0, P([0,-1,0])), 'E', 'P=(0,-1,0) E2 (A1) (1,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_A1_1_0(L_L_s_I0, P([1,0,0])), 'E', 'P=(1,0,0) E2 (A1) (1,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_A1_1_0(L_L_s_I0, P([-1,0,0])), 'E', 'P=(-1,0,0) E2 (A1) (1,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_A1_1_0(N_X_a_I0, P([0,0,1])), 'E', 'P=(0,0,1) E2 (A1) (1,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_A1_1_0(N_X_a_I0, P([0,0,-1])), 'E', 'P=(0,0,-1) E2 (A1) (1,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_A1_1_0(N_X_a_I0, P([0,1,0])), 'E', 'P=(0,1,0) E2 (A1) (1,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_A1_1_0(N_X_a_I0, P([0,-1,0])), 'E', 'P=(0,-1,0) E2 (A1) (1,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_A1_1_0(N_X_a_I0, P([1,0,0])), 'E', 'P=(1,0,0) E2 (A1) (1,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_A1_1_0(N_X_a_I0, P([-1,0,0])), 'E', 'P=(-1,0,0) E2 (A1) (1,0); NXa I=0', False); total += 1  # PASSED

  #passed += test_irrep(P1_E2_E2_2_1_S0(L_L_s_I0, P([0,0,1])), 'E', 'P=(0,0,1) E2 (E2) (2,1) S0; LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_E2_2_1_S0(L_L_s_I0, P([0,0,-1])), 'E', 'P=(0,0,-1) E2 (E2) (2,1) S0; LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_E2_2_1_S0(L_L_s_I0, P([0,1,0])), 'E', 'P=(0,1,0) E2 (E2) (2,1) S0; LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_E2_2_1_S0(L_L_s_I0, P([0,-1,0])), 'E', 'P=(0,-1,0) E2 (E2) (2,1) S0; LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_E2_2_1_S0(L_L_s_I0, P([1,0,0])), 'E', 'P=(1,0,0) E2 (E2) (2,1) S0; LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_E2_2_1_S0(L_L_s_I0, P([-1,0,0])), 'E', 'P=(-1,0,0) E2 (E2) (2,1) S0; LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_E2_2_1_S0(N_X_a_I0, P([0,0,1])), 'E', 'P=(0,0,1) E2 (E2) (2,1) S0; NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_E2_2_1_S0(N_X_a_I0, P([0,0,-1])), 'E', 'P=(0,0,-1) E2 (E2) (2,1) S0; NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_E2_2_1_S0(N_X_a_I0, P([0,1,0])), 'E', 'P=(0,1,0) E2 (E2) (2,1) S0; NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_E2_2_1_S0(N_X_a_I0, P([0,-1,0])), 'E', 'P=(0,-1,0) E2 (E2) (2,1) S0; NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_E2_2_1_S0(N_X_a_I0, P([1,0,0])), 'E', 'P=(1,0,0) E2 (E2) (2,1) S0; NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_E2_2_1_S0(N_X_a_I0, P([-1,0,0])), 'E', 'P=(-1,0,0) E2 (E2) (2,1) S0; NXa I=0', False); total += 1  # PASSED

  #passed += test_irrep(P1_E2_A1_2_1(L_L_s_I0, P([0,0,1])), 'E', 'P=(0,0,1) E2 (A1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_A1_2_1(L_L_s_I0, P([0,0,-1])), 'E', 'P=(0,0,-1) E2 (A1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_A1_2_1(L_L_s_I0, P([0,1,0])), 'E', 'P=(0,1,0) E2 (A1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_A1_2_1(L_L_s_I0, P([0,-1,0])), 'E', 'P=(0,-1,0) E2 (A1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_A1_2_1(L_L_s_I0, P([1,0,0])), 'E', 'P=(1,0,0) E2 (A1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_A1_2_1(L_L_s_I0, P([-1,0,0])), 'E', 'P=(-1,0,0) E2 (A1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_A1_2_1(N_X_a_I0, P([0,0,1])), 'E', 'P=(0,0,1) E2 (A1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_A1_2_1(N_X_a_I0, P([0,0,-1])), 'E', 'P=(0,0,-1) E2 (A1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_A1_2_1(N_X_a_I0, P([0,1,0])), 'E', 'P=(0,1,0) E2 (A1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_A1_2_1(N_X_a_I0, P([0,-1,0])), 'E', 'P=(0,-1,0) E2 (A1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_A1_2_1(N_X_a_I0, P([1,0,0])), 'E', 'P=(1,0,0) E2 (A1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_A1_2_1(N_X_a_I0, P([-1,0,0])), 'E', 'P=(-1,0,0) E2 (A1) (2,1); NXa I=0', False); total += 1  # PASSED

  #passed += test_irrep(P1_E2_B1_2_1(L_L_s_I0, P([0,0,1])), 'E', 'P=(0,0,1) E2 (B1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_B1_2_1(L_L_s_I0, P([0,0,-1])), 'E', 'P=(0,0,-1) E2 (B1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_B1_2_1(L_L_s_I0, P([0,1,0])), 'E', 'P=(0,1,0) E2 (B1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_B1_2_1(L_L_s_I0, P([0,-1,0])), 'E', 'P=(0,-1,0) E2 (B1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_B1_2_1(L_L_s_I0, P([1,0,0])), 'E', 'P=(1,0,0) E2 (B1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_B1_2_1(L_L_s_I0, P([-1,0,0])), 'E', 'P=(-1,0,0) E2 (B1) (2,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_B1_2_1(N_X_a_I0, P([0,0,1])), 'E', 'P=(0,0,1) E2 (B1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_B1_2_1(N_X_a_I0, P([0,0,-1])), 'E', 'P=(0,0,-1) E2 (B1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_B1_2_1(N_X_a_I0, P([0,1,0])), 'E', 'P=(0,1,0) E2 (B1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_B1_2_1(N_X_a_I0, P([0,-1,0])), 'E', 'P=(0,-1,0) E2 (B1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_B1_2_1(N_X_a_I0, P([1,0,0])), 'E', 'P=(1,0,0) E2 (B1) (2,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_B1_2_1(N_X_a_I0, P([-1,0,0])), 'E', 'P=(-1,0,0) E2 (B1) (2,1); NXa I=0', False); total += 1  # PASSED

  #passed += test_irrep(P1_E2_E2_2_1_S1(L_L_s_I0, P([0,0,1])), 'E', 'P=(0,0,1) E2 (E2) (2,1) S1; LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_E2_2_1_S1(L_L_s_I0, P([0,0,-1])), 'E', 'P=(0,0,-1) E2 (E2) (2,1) S1; LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_E2_2_1_S1(L_L_s_I0, P([0,1,0])), 'E', 'P=(0,1,0) E2 (E2) (2,1) S1; LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_E2_2_1_S1(L_L_s_I0, P([0,-1,0])), 'E', 'P=(0,-1,0) E2 (E2) (2,1) S1; LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_E2_2_1_S1(L_L_s_I0, P([1,0,0])), 'E', 'P=(1,0,0) E2 (E2) (2,1) S1; LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_E2_2_1_S1(L_L_s_I0, P([-1,0,0])), 'E', 'P=(-1,0,0) E2 (E2) (2,1) S1; LL', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_E2_2_1_S1(N_X_a_I0, P([0,0,1])), 'E', 'P=(0,0,1) E2 (E2) (2,1) S1; NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_E2_2_1_S1(N_X_a_I0, P([0,0,-1])), 'E', 'P=(0,0,-1) E2 (E2) (2,1) S1; NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_E2_2_1_S1(N_X_a_I0, P([0,1,0])), 'E', 'P=(0,1,0) E2 (E2) (2,1) S1; NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_E2_2_1_S1(N_X_a_I0, P([0,-1,0])), 'E', 'P=(0,-1,0) E2 (E2) (2,1) S1; NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_E2_2_1_S1(N_X_a_I0, P([1,0,0])), 'E', 'P=(1,0,0) E2 (E2) (2,1) S1; NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P1_E2_E2_2_1_S1(N_X_a_I0, P([-1,0,0])), 'E', 'P=(-1,0,0) E2 (E2) (2,1) S1; NXa I=0', False); total += 1  # PASSED

  # Psq = 2 (A1)
  #passed += test_irrep(P2_A1_A1_2_0(L_L_s_I0, P([0,1,1])), 'A1', 'P=(0,1,1) A1 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_2_0(L_L_s_I0, P([0,1,-1])), 'A1', 'P=(0,1,-1) A1 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_2_0(L_L_s_I0, P([0,-1,1])), 'A1', 'P=(0,-1,1) A1 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_2_0(L_L_s_I0, P([0,-1,-1])), 'A1', 'P=(0,-1,-1) A1 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_2_0(L_L_s_I0, P([1,0,1])), 'A1', 'P=(1,0,1) A1 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_2_0(L_L_s_I0, P([1,0,-1])), 'A1', 'P=(1,0,-1) A1 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_2_0(L_L_s_I0, P([-1,0,1])), 'A1', 'P=(-1,0,1) A1 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_2_0(L_L_s_I0, P([-1,0,-1])), 'A1', 'P=(-1,0,-1) A1 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_2_0(L_L_s_I0, P([1,1,0])), 'A1', 'P=(1,1,0) A1 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_2_0(L_L_s_I0, P([1,-1,0])), 'A1', 'P=(1,-1,0) A1 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_2_0(L_L_s_I0, P([-1,1,0])), 'A1', 'P=(-1,1,0) A1 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_2_0(L_L_s_I0, P([-1,-1,0])), 'A1', 'P=(-1,-1,0) A1 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_2_0(N_X_a_I0, P([0,1,1])), 'A1', 'P=(0,1,1) A1 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_2_0(N_X_a_I0, P([0,1,-1])), 'A1', 'P=(0,1,-1) A1 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_2_0(N_X_a_I0, P([0,-1,1])), 'A1', 'P=(0,-1,1) A1 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_2_0(N_X_a_I0, P([0,-1,-1])), 'A1', 'P=(0,-1,-1) A1 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_2_0(N_X_a_I0, P([1,0,1])), 'A1', 'P=(1,0,1) A1 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_2_0(N_X_a_I0, P([1,0,-1])), 'A1', 'P=(1,0,-1) A1 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_2_0(N_X_a_I0, P([-1,0,1])), 'A1', 'P=(-1,0,1) A1 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_2_0(N_X_a_I0, P([-1,0,-1])), 'A1', 'P=(-1,0,-1) A1 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_2_0(N_X_a_I0, P([1,1,0])), 'A1', 'P=(1,1,0) A1 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_2_0(N_X_a_I0, P([1,-1,0])), 'A1', 'P=(1,-1,0) A1 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_2_0(N_X_a_I0, P([-1,1,0])), 'A1', 'P=(-1,1,0) A1 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_2_0(N_X_a_I0, P([-1,-1,0])), 'A1', 'P=(-1,-1,0) A1 (A1) (2,0); NXa I=0', False); total += 1  # PASSED

  #passed += test_irrep(P2_A1_A1_1_1(L_L_s_I0, P([0,1,1])), 'A1', 'P=(0,1,1) A1 (A1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_1_1(L_L_s_I0, P([0,1,-1])), 'A1', 'P=(0,1,-1) A1 (A1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_1_1(L_L_s_I0, P([0,-1,1])), 'A1', 'P=(0,-1,1) A1 (A1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_1_1(L_L_s_I0, P([0,-1,-1])), 'A1', 'P=(0,-1,-1) A1 (A1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_1_1(L_L_s_I0, P([1,0,1])), 'A1', 'P=(1,0,1) A1 (A1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_1_1(L_L_s_I0, P([1,0,-1])), 'A1', 'P=(1,0,-1) A1 (A1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_1_1(L_L_s_I0, P([-1,0,1])), 'A1', 'P=(-1,0,1) A1 (A1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_1_1(L_L_s_I0, P([-1,0,-1])), 'A1', 'P=(-1,0,-1) A1 (A1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_1_1(L_L_s_I0, P([1,1,0])), 'A1', 'P=(1,1,0) A1 (A1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_1_1(L_L_s_I0, P([1,-1,0])), 'A1', 'P=(1,-1,0) A1 (A1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_1_1(L_L_s_I0, P([-1,1,0])), 'A1', 'P=(-1,1,0) A1 (A1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_A1_1_1(L_L_s_I0, P([-1,-1,0])), 'A1', 'P=(-1,-1,0) A1 (A1) (1,1); LL', False); total += 1  # PASSED

  #passed += test_irrep(P2_A1_B1_1_1(L_L_s_I0, P([0,1,1])), 'A1', 'P=(0,1,1) A1 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_B1_1_1(L_L_s_I0, P([0,1,-1])), 'A1', 'P=(0,1,-1) A1 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_B1_1_1(L_L_s_I0, P([0,-1,1])), 'A1', 'P=(0,-1,1) A1 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_B1_1_1(L_L_s_I0, P([0,-1,-1])), 'A1', 'P=(0,-1,-1) A1 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_B1_1_1(L_L_s_I0, P([1,0,1])), 'A1', 'P=(1,0,1) A1 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_B1_1_1(L_L_s_I0, P([1,0,-1])), 'A1', 'P=(1,0,-1) A1 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_B1_1_1(L_L_s_I0, P([-1,0,1])), 'A1', 'P=(-1,0,1) A1 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_B1_1_1(L_L_s_I0, P([-1,0,-1])), 'A1', 'P=(-1,0,-1) A1 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_B1_1_1(L_L_s_I0, P([1,1,0])), 'A1', 'P=(1,1,0) A1 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_B1_1_1(L_L_s_I0, P([1,-1,0])), 'A1', 'P=(1,-1,0) A1 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_B1_1_1(L_L_s_I0, P([-1,1,0])), 'A1', 'P=(-1,1,0) A1 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A1_B1_1_1(L_L_s_I0, P([-1,-1,0])), 'A1', 'P=(-1,-1,0) A1 (B1) (1,1); LL', False); total += 1  # PASSED


  # Psq = 2 (A2)
  #passed += test_irrep(P2_A2_A1_2_0(L_L_s_I0, P([0,1,1])), 'A2', 'P=(0,1,1) A2 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_A1_2_0(L_L_s_I0, P([0,1,-1])), 'A2', 'P=(0,1,-1) A2 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_A1_2_0(L_L_s_I0, P([0,-1,1])), 'A2', 'P=(0,-1,1) A2 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_A1_2_0(L_L_s_I0, P([0,-1,-1])), 'A2', 'P=(0,-1,-1) A2 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_A1_2_0(L_L_s_I0, P([1,0,1])), 'A2', 'P=(1,0,1) A2 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_A1_2_0(L_L_s_I0, P([1,0,-1])), 'A2', 'P=(1,0,-1) A2 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_A1_2_0(L_L_s_I0, P([-1,0,1])), 'A2', 'P=(-1,0,1) A2 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_A1_2_0(L_L_s_I0, P([-1,0,-1])), 'A2', 'P=(-1,0,-1) A2 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_A1_2_0(L_L_s_I0, P([1,1,0])), 'A2', 'P=(1,1,0) A2 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_A1_2_0(L_L_s_I0, P([1,-1,0])), 'A2', 'P=(1,-1,0) A2 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_A1_2_0(L_L_s_I0, P([-1,1,0])), 'A2', 'P=(-1,1,0) A2 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_A1_2_0(L_L_s_I0, P([-1,-1,0])), 'A2', 'P=(-1,-1,0) A2 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_A1_2_0(N_X_a_I0, P([0,1,1])), 'A2', 'P=(0,1,1) A2 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_A1_2_0(N_X_a_I0, P([0,1,-1])), 'A2', 'P=(0,1,-1) A2 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_A1_2_0(N_X_a_I0, P([0,-1,1])), 'A2', 'P=(0,-1,1) A2 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_A1_2_0(N_X_a_I0, P([0,-1,-1])), 'A2', 'P=(0,-1,-1) A2 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_A1_2_0(N_X_a_I0, P([1,0,1])), 'A2', 'P=(1,0,1) A2 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_A1_2_0(N_X_a_I0, P([1,0,-1])), 'A2', 'P=(1,0,-1) A2 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_A1_2_0(N_X_a_I0, P([-1,0,1])), 'A2', 'P=(-1,0,1) A2 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_A1_2_0(N_X_a_I0, P([-1,0,-1])), 'A2', 'P=(-1,0,-1) A2 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_A1_2_0(N_X_a_I0, P([1,1,0])), 'A2', 'P=(1,1,0) A2 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_A1_2_0(N_X_a_I0, P([1,-1,0])), 'A2', 'P=(1,-1,0) A2 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_A1_2_0(N_X_a_I0, P([-1,1,0])), 'A2', 'P=(-1,1,0) A2 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_A1_2_0(N_X_a_I0, P([-1,-1,0])), 'A2', 'P=(-1,-1,0) A2 (A1) (2,0); NXa I=0', False); total += 1  # PASSED

  #passed += test_irrep(P2_A2_B1_1_1_Fs(L_L_s_I0, P([0,1,1])), 'A2', 'P=(0,1,1) A2 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_B1_1_1_Fs(L_L_s_I0, P([0,1,-1])), 'A2', 'P=(0,1,-1) A2 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_B1_1_1_Fs(L_L_s_I0, P([0,-1,1])), 'A2', 'P=(0,-1,1) A2 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_B1_1_1_Fs(L_L_s_I0, P([0,-1,-1])), 'A2', 'P=(0,-1,-1) A2 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_B1_1_1_Fs(L_L_s_I0, P([1,0,1])), 'A2', 'P=(1,0,1) A2 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_B1_1_1_Fs(L_L_s_I0, P([1,0,-1])), 'A2', 'P=(1,0,-1) A2 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_B1_1_1_Fs(L_L_s_I0, P([-1,0,1])), 'A2', 'P=(-1,0,1) A2 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_B1_1_1_Fs(L_L_s_I0, P([-1,0,-1])), 'A2', 'P=(-1,0,-1) A2 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_B1_1_1_Fs(L_L_s_I0, P([1,1,0])), 'A2', 'P=(1,1,0) A2 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_B1_1_1_Fs(L_L_s_I0, P([1,-1,0])), 'A2', 'P=(1,-1,0) A2 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_B1_1_1_Fs(L_L_s_I0, P([-1,1,0])), 'A2', 'P=(-1,1,0) A2 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_B1_1_1_Fs(L_L_s_I0, P([-1,-1,0])), 'A2', 'P=(-1,-1,0) A2 (B1) (1,1); LL', False); total += 1  # PASSED

  #passed += test_irrep(P2_A2_B1_1_1_Fa(N_X_a_I0, P([0,1,1])), 'A2', 'P=(0,1,1) A2 (B1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_B1_1_1_Fa(N_X_a_I0, P([0,1,-1])), 'A2', 'P=(0,1,-1) A2 (B1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_B1_1_1_Fa(N_X_a_I0, P([0,-1,1])), 'A2', 'P=(0,-1,1) A2 (B1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_B1_1_1_Fa(N_X_a_I0, P([0,-1,-1])), 'A2', 'P=(0,-1,-1) A2 (B1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_B1_1_1_Fa(N_X_a_I0, P([1,0,1])), 'A2', 'P=(1,0,1) A2 (B1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_B1_1_1_Fa(N_X_a_I0, P([1,0,-1])), 'A2', 'P=(1,0,-1) A2 (B1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_B1_1_1_Fa(N_X_a_I0, P([-1,0,1])), 'A2', 'P=(-1,0,1) A2 (B1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_B1_1_1_Fa(N_X_a_I0, P([-1,0,-1])), 'A2', 'P=(-1,0,-1) A2 (B1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_B1_1_1_Fa(N_X_a_I0, P([1,1,0])), 'A2', 'P=(1,1,0) A2 (B1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_B1_1_1_Fa(N_X_a_I0, P([1,-1,0])), 'A2', 'P=(1,-1,0) A2 (B1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_B1_1_1_Fa(N_X_a_I0, P([-1,1,0])), 'A2', 'P=(-1,1,0) A2 (B1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_A2_B1_1_1_Fa(N_X_a_I0, P([-1,-1,0])), 'A2', 'P=(-1,-1,0) A2 (B1) (1,1); NXa I=0', False); total += 1  # PASSED


  # Psq = 2 (B1)
  #passed += test_irrep(P2_B1_A1_2_0(L_L_s_I0, P([0,1,1])), 'B2', 'P=(0,1,1) B1 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_2_0(L_L_s_I0, P([0,1,-1])), 'B2', 'P=(0,1,-1) B1 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_2_0(L_L_s_I0, P([0,-1,1])), 'B2', 'P=(0,-1,1) B1 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_2_0(L_L_s_I0, P([0,-1,-1])), 'B2', 'P=(0,-1,-1) B1 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_2_0(L_L_s_I0, P([1,0,1])), 'B2', 'P=(1,0,1) B1 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_2_0(L_L_s_I0, P([1,0,-1])), 'B2', 'P=(1,0,-1) B1 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_2_0(L_L_s_I0, P([-1,0,1])), 'B2', 'P=(-1,0,1) B1 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_2_0(L_L_s_I0, P([-1,0,-1])), 'B2', 'P=(-1,0,-1) B1 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_2_0(L_L_s_I0, P([1,1,0])), 'B2', 'P=(1,1,0) B1 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_2_0(L_L_s_I0, P([1,-1,0])), 'B2', 'P=(1,-1,0) B1 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_2_0(L_L_s_I0, P([-1,1,0])), 'B2', 'P=(-1,1,0) B1 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_2_0(L_L_s_I0, P([-1,-1,0])), 'B2', 'P=(-1,-1,0) B1 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_2_0(N_X_a_I0, P([0,1,1])), 'B2', 'P=(0,1,1) B1 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_2_0(N_X_a_I0, P([0,1,-1])), 'B2', 'P=(0,1,-1) B1 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_2_0(N_X_a_I0, P([0,-1,1])), 'B2', 'P=(0,-1,1) B1 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_2_0(N_X_a_I0, P([0,-1,-1])), 'B2', 'P=(0,-1,-1) B1 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_2_0(N_X_a_I0, P([1,0,1])), 'B2', 'P=(1,0,1) B1 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_2_0(N_X_a_I0, P([1,0,-1])), 'B2', 'P=(1,0,-1) B1 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_2_0(N_X_a_I0, P([-1,0,1])), 'B2', 'P=(-1,0,1) B1 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_2_0(N_X_a_I0, P([-1,0,-1])), 'B2', 'P=(-1,0,-1) B1 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_2_0(N_X_a_I0, P([1,1,0])), 'B2', 'P=(1,1,0) B1 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_2_0(N_X_a_I0, P([1,-1,0])), 'B2', 'P=(1,-1,0) B1 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_2_0(N_X_a_I0, P([-1,1,0])), 'B2', 'P=(-1,1,0) B1 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_2_0(N_X_a_I0, P([-1,-1,0])), 'B2', 'P=(-1,-1,0) B1 (A1) (2,0); NXa I=0', False); total += 1  # PASSED

  #passed += test_irrep(P2_B1_A1_1_1(N_X_a_I0, P([0,1,1])), 'B2', 'P=(0,1,1) B1 (A1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_1_1(N_X_a_I0, P([0,1,-1])), 'B2', 'P=(0,1,-1) B1 (A1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_1_1(N_X_a_I0, P([0,-1,1])), 'B2', 'P=(0,-1,1) B1 (A1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_1_1(N_X_a_I0, P([0,-1,-1])), 'B2', 'P=(0,-1,-1) B1 (A1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_1_1(N_X_a_I0, P([1,0,1])), 'B2', 'P=(1,0,1) B1 (A1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_1_1(N_X_a_I0, P([1,0,-1])), 'B2', 'P=(1,0,-1) B1 (A1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_1_1(N_X_a_I0, P([-1,0,1])), 'B2', 'P=(-1,0,1) B1 (A1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_1_1(N_X_a_I0, P([-1,0,-1])), 'B2', 'P=(-1,0,-1) B1 (A1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_1_1(N_X_a_I0, P([1,1,0])), 'B2', 'P=(1,1,0) B1 (A1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_1_1(N_X_a_I0, P([1,-1,0])), 'B2', 'P=(1,-1,0) B1 (A1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_1_1(N_X_a_I0, P([-1,1,0])), 'B2', 'P=(-1,1,0) B1 (A1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_A1_1_1(N_X_a_I0, P([-1,-1,0])), 'B2', 'P=(-1,-1,0) B1 (A1) (1,1); NXa I=0', False); total += 1  # PASSED

  #passed += test_irrep(P2_B1_B1_1_1(N_X_a_I0, P([0,1,1])), 'B2', 'P=(0,1,1) B1 (B1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_B1_1_1(N_X_a_I0, P([0,1,-1])), 'B2', 'P=(0,1,-1) B1 (B1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_B1_1_1(N_X_a_I0, P([0,-1,1])), 'B2', 'P=(0,-1,1) B1 (B1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_B1_1_1(N_X_a_I0, P([0,-1,-1])), 'B2', 'P=(0,-1,-1) B1 (B1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_B1_1_1(N_X_a_I0, P([1,0,1])), 'B2', 'P=(1,0,1) B1 (B1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_B1_1_1(N_X_a_I0, P([1,0,-1])), 'B2', 'P=(1,0,-1) B1 (B1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_B1_1_1(N_X_a_I0, P([-1,0,1])), 'B2', 'P=(-1,0,1) B1 (B1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_B1_1_1(N_X_a_I0, P([-1,0,-1])), 'B2', 'P=(-1,0,-1) B1 (B1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_B1_1_1(N_X_a_I0, P([1,1,0])), 'B2', 'P=(1,1,0) B1 (B1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_B1_1_1(N_X_a_I0, P([1,-1,0])), 'B2', 'P=(1,-1,0) B1 (B1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_B1_1_1(N_X_a_I0, P([-1,1,0])), 'B2', 'P=(-1,1,0) B1 (B1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B1_B1_1_1(N_X_a_I0, P([-1,-1,0])), 'B2', 'P=(-1,-1,0) B1 (B1) (1,1); NXa I=0', False); total += 1  # PASSED


  # Psq = 2 (B2)
  #passed += test_irrep(P2_B2_A1_2_0(L_L_s_I0, P([0,1,1])), 'B1', 'P=(0,1,1) B2 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_2_0(L_L_s_I0, P([0,1,-1])), 'B1', 'P=(0,1,-1) B2 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_2_0(L_L_s_I0, P([0,-1,1])), 'B1', 'P=(0,-1,1) B2 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_2_0(L_L_s_I0, P([0,-1,-1])), 'B1', 'P=(0,-1,-1) B2 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_2_0(L_L_s_I0, P([1,0,1])), 'B1', 'P=(1,0,1) B2 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_2_0(L_L_s_I0, P([1,0,-1])), 'B1', 'P=(1,0,-1) B2 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_2_0(L_L_s_I0, P([-1,0,1])), 'B1', 'P=(-1,0,1) B2 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_2_0(L_L_s_I0, P([-1,0,-1])), 'B1', 'P=(-1,0,-1) B2 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_2_0(L_L_s_I0, P([1,1,0])), 'B1', 'P=(1,1,0) B2 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_2_0(L_L_s_I0, P([1,-1,0])), 'B1', 'P=(1,-1,0) B2 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_2_0(L_L_s_I0, P([-1,1,0])), 'B1', 'P=(-1,1,0) B2 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_2_0(L_L_s_I0, P([-1,-1,0])), 'B1', 'P=(-1,-1,0) B2 (A1) (2,0); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_2_0(N_X_a_I0, P([0,1,1])), 'B1', 'P=(0,1,1) B2 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_2_0(N_X_a_I0, P([0,1,-1])), 'B1', 'P=(0,1,-1) B2 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_2_0(N_X_a_I0, P([0,-1,1])), 'B1', 'P=(0,-1,1) B2 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_2_0(N_X_a_I0, P([0,-1,-1])), 'B1', 'P=(0,-1,-1) B2 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_2_0(N_X_a_I0, P([1,0,1])), 'B1', 'P=(1,0,1) B2 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_2_0(N_X_a_I0, P([1,0,-1])), 'B1', 'P=(1,0,-1) B2 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_2_0(N_X_a_I0, P([-1,0,1])), 'B1', 'P=(-1,0,1) B2 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_2_0(N_X_a_I0, P([-1,0,-1])), 'B1', 'P=(-1,0,-1) B2 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_2_0(N_X_a_I0, P([1,1,0])), 'B1', 'P=(1,1,0) B2 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_2_0(N_X_a_I0, P([1,-1,0])), 'B1', 'P=(1,-1,0) B2 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_2_0(N_X_a_I0, P([-1,1,0])), 'B1', 'P=(-1,1,0) B2 (A1) (2,0); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_2_0(N_X_a_I0, P([-1,-1,0])), 'B1', 'P=(-1,-1,0) B2 (A1) (2,0); NXa I=0', False); total += 1  # PASSED

  #passed += test_irrep(P2_B2_B1_1_1(L_L_s_I0, P([0,1,1])), 'B1', 'P=(0,1,1) B2 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_B1_1_1(L_L_s_I0, P([0,1,-1])), 'B1', 'P=(0,1,-1) B2 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_B1_1_1(L_L_s_I0, P([0,-1,1])), 'B1', 'P=(0,-1,1) B2 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_B1_1_1(L_L_s_I0, P([0,-1,-1])), 'B1', 'P=(0,-1,-1) B2 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_B1_1_1(L_L_s_I0, P([1,0,1])), 'B1', 'P=(1,0,1) B2 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_B1_1_1(L_L_s_I0, P([1,0,-1])), 'B1', 'P=(1,0,-1) B2 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_B1_1_1(L_L_s_I0, P([-1,0,1])), 'B1', 'P=(-1,0,1) B2 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_B1_1_1(L_L_s_I0, P([-1,0,-1])), 'B1', 'P=(-1,0,-1) B2 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_B1_1_1(L_L_s_I0, P([1,1,0])), 'B1', 'P=(1,1,0) B2 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_B1_1_1(L_L_s_I0, P([1,-1,0])), 'B1', 'P=(1,-1,0) B2 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_B1_1_1(L_L_s_I0, P([-1,1,0])), 'B1', 'P=(-1,1,0) B2 (B1) (1,1); LL', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_B1_1_1(L_L_s_I0, P([-1,-1,0])), 'B1', 'P=(-1,-1,0) B2 (B1) (1,1); LL', False); total += 1  # PASSED

  #passed += test_irrep(P2_B2_A1_1_1(N_X_a_I0, P([0,1,1])), 'B1', 'P=(0,1,1) B2 (A1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_1_1(N_X_a_I0, P([0,1,-1])), 'B1', 'P=(0,1,-1) B2 (A1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_1_1(N_X_a_I0, P([0,-1,1])), 'B1', 'P=(0,-1,1) B2 (A1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_1_1(N_X_a_I0, P([0,-1,-1])), 'B1', 'P=(0,-1,-1) B2 (A1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_1_1(N_X_a_I0, P([1,0,1])), 'B1', 'P=(1,0,1) B2 (A1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_1_1(N_X_a_I0, P([1,0,-1])), 'B1', 'P=(1,0,-1) B2 (A1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_1_1(N_X_a_I0, P([-1,0,1])), 'B1', 'P=(-1,0,1) B2 (A1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_1_1(N_X_a_I0, P([-1,0,-1])), 'B1', 'P=(-1,0,-1) B2 (A1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_1_1(N_X_a_I0, P([1,1,0])), 'B1', 'P=(1,1,0) B2 (A1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_1_1(N_X_a_I0, P([1,-1,0])), 'B1', 'P=(1,-1,0) B2 (A1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_1_1(N_X_a_I0, P([-1,1,0])), 'B1', 'P=(-1,1,0) B2 (A1) (1,1); NXa I=0', False); total += 1  # PASSED
  #passed += test_irrep(P2_B2_A1_1_1(N_X_a_I0, P([-1,-1,0])), 'B1', 'P=(-1,-1,0) B2 (A1) (1,1); NXa I=0', False); total += 1  # PASSED


  # Psq = 3 (A1)
  #passed += test_irrep(P3_A1_A1_3_0(L_L_s_I0, P([1,1,1])), 'A1', 'P=(1,1,1) A1 (A1) (3,0); LL', False); total += 1
  #passed += test_irrep(P3_A1_A1_3_0(L_L_s_I0, P([1,1,-1])), 'A1', 'P=(1,1,-1) A1 (A1) (3,0); LL', False); total += 1
  #passed += test_irrep(P3_A1_A1_3_0(L_L_s_I0, P([1,-1,1])), 'A1', 'P=(1,-1,1) A1 (A1) (3,0); LL', False); total += 1
  #passed += test_irrep(P3_A1_A1_3_0(L_L_s_I0, P([-1,1,1])), 'A1', 'P=(-1,1,1) A1 (A1) (3,0); LL', False); total += 1
  #passed += test_irrep(P3_A1_A1_3_0(L_L_s_I0, P([1,-1,-1])), 'A1', 'P=(1,-1,-1) A1 (A1) (3,0); LL', False); total += 1
  #passed += test_irrep(P3_A1_A1_3_0(L_L_s_I0, P([-1,1,-1])), 'A1', 'P=(-1,1,-1) A1 (A1) (3,0); LL', False); total += 1
  #passed += test_irrep(P3_A1_A1_3_0(L_L_s_I0, P([-1,-1,1])), 'A1', 'P=(-1,-1,1) A1 (A1) (3,0); LL', False); total += 1
  #passed += test_irrep(P3_A1_A1_3_0(L_L_s_I0, P([-1,-1,-1])), 'A1', 'P=(-1,-1,-1) A1 (A1) (3,0); LL', False); total += 1
  #passed += test_irrep(P3_A1_A1_3_0(N_X_a_I0, P([1,1,1])), 'A1', 'P=(1,1,1) A1 (A1) (3,0); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A1_A1_3_0(N_X_a_I0, P([1,1,-1])), 'A1', 'P=(1,1,-1) A1 (A1) (3,0); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A1_A1_3_0(N_X_a_I0, P([1,-1,1])), 'A1', 'P=(1,-1,1) A1 (A1) (3,0); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A1_A1_3_0(N_X_a_I0, P([-1,1,1])), 'A1', 'P=(-1,1,1) A1 (A1) (3,0); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A1_A1_3_0(N_X_a_I0, P([1,-1,-1])), 'A1', 'P=(1,-1,-1) A1 (A1) (3,0); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A1_A1_3_0(N_X_a_I0, P([-1,1,-1])), 'A1', 'P=(-1,1,-1) A1 (A1) (3,0); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A1_A1_3_0(N_X_a_I0, P([-1,-1,1])), 'A1', 'P=(-1,-1,1) A1 (A1) (3,0); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A1_A1_3_0(N_X_a_I0, P([-1,-1,-1])), 'A1', 'P=(-1,-1,-1) A1 (A1) (3,0); NXa I=0', False); total += 1

  #passed += test_irrep(P3_A1_A1_2_1(L_L_s_I0, P([1,1,1])), 'A1', 'P=(1,1,1) A1 (A1) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A1_A1_2_1(L_L_s_I0, P([1,1,-1])), 'A1', 'P=(1,1,-1) A1 (A1) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A1_A1_2_1(L_L_s_I0, P([1,-1,1])), 'A1', 'P=(1,-1,1) A1 (A1) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A1_A1_2_1(L_L_s_I0, P([-1,1,1])), 'A1', 'P=(-1,1,1) A1 (A1) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A1_A1_2_1(L_L_s_I0, P([1,-1,-1])), 'A1', 'P=(1,-1,-1) A1 (A1) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A1_A1_2_1(L_L_s_I0, P([-1,1,-1])), 'A1', 'P=(-1,1,-1) A1 (A1) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A1_A1_2_1(L_L_s_I0, P([-1,-1,1])), 'A1', 'P=(-1,-1,1) A1 (A1) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A1_A1_2_1(L_L_s_I0, P([-1,-1,-1])), 'A1', 'P=(-1,-1,-1) A1 (A1) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A1_A1_2_1(N_X_a_I0, P([1,1,1])), 'A1', 'P=(1,1,1) A1 (A1) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A1_A1_2_1(N_X_a_I0, P([1,1,-1])), 'A1', 'P=(1,1,-1) A1 (A1) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A1_A1_2_1(N_X_a_I0, P([1,-1,1])), 'A1', 'P=(1,-1,1) A1 (A1) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A1_A1_2_1(N_X_a_I0, P([-1,1,1])), 'A1', 'P=(-1,1,1) A1 (A1) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A1_A1_2_1(N_X_a_I0, P([1,-1,-1])), 'A1', 'P=(1,-1,-1) A1 (A1) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A1_A1_2_1(N_X_a_I0, P([-1,1,-1])), 'A1', 'P=(-1,1,-1) A1 (A1) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A1_A1_2_1(N_X_a_I0, P([-1,-1,1])), 'A1', 'P=(-1,-1,1) A1 (A1) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A1_A1_2_1(N_X_a_I0, P([-1,-1,-1])), 'A1', 'P=(-1,-1,-1) A1 (A1) (2,1); NXa I=0', False); total += 1

  #passed += test_irrep(P3_A1_E2_2_1(L_L_s_I0, P([1,1,1])), 'A1', 'P=(1,1,1) A1 (E2) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A1_E2_2_1(L_L_s_I0, P([1,1,-1])), 'A1', 'P=(1,1,-1) A1 (E2) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A1_E2_2_1(L_L_s_I0, P([1,-1,1])), 'A1', 'P=(1,-1,1) A1 (E2) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A1_E2_2_1(L_L_s_I0, P([-1,1,1])), 'A1', 'P=(-1,1,1) A1 (E2) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A1_E2_2_1(L_L_s_I0, P([1,-1,-1])), 'A1', 'P=(1,-1,-1) A1 (E2) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A1_E2_2_1(L_L_s_I0, P([-1,1,-1])), 'A1', 'P=(-1,1,-1) A1 (E2) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A1_E2_2_1(L_L_s_I0, P([-1,-1,1])), 'A1', 'P=(-1,-1,1) A1 (E2) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A1_E2_2_1(L_L_s_I0, P([-1,-1,-1])), 'A1', 'P=(-1,-1,-1) A1 (E2) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A1_E2_2_1(N_X_a_I0, P([1,1,1])), 'A1', 'P=(1,1,1) A1 (E2) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A1_E2_2_1(N_X_a_I0, P([1,1,-1])), 'A1', 'P=(1,1,-1) A1 (E2) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A1_E2_2_1(N_X_a_I0, P([1,-1,1])), 'A1', 'P=(1,-1,1) A1 (E2) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A1_E2_2_1(N_X_a_I0, P([-1,1,1])), 'A1', 'P=(-1,1,1) A1 (E2) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A1_E2_2_1(N_X_a_I0, P([1,-1,-1])), 'A1', 'P=(1,-1,-1) A1 (E2) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A1_E2_2_1(N_X_a_I0, P([-1,1,-1])), 'A1', 'P=(-1,1,-1) A1 (E2) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A1_E2_2_1(N_X_a_I0, P([-1,-1,1])), 'A1', 'P=(-1,-1,1) A1 (E2) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A1_E2_2_1(N_X_a_I0, P([-1,-1,-1])), 'A1', 'P=(-1,-1,-1) A1 (E2) (2,1); NXa I=0', False); total += 1


  # Psq = 3 (A2)
  #passed += test_irrep(P3_A2_A1_3_0(L_L_s_I0, P([1,1,1])), 'A2', 'P=(1,1,1) A2 (A1) (3,0); LL', False); total += 1
  #passed += test_irrep(P3_A2_A1_3_0(L_L_s_I0, P([1,1,-1])), 'A2', 'P=(1,1,-1) A2 (A1) (3,0); LL', False); total += 1
  #passed += test_irrep(P3_A2_A1_3_0(L_L_s_I0, P([1,-1,1])), 'A2', 'P=(1,-1,1) A2 (A1) (3,0); LL', False); total += 1
  #passed += test_irrep(P3_A2_A1_3_0(L_L_s_I0, P([-1,1,1])), 'A2', 'P=(-1,1,1) A2 (A1) (3,0); LL', False); total += 1
  #passed += test_irrep(P3_A2_A1_3_0(L_L_s_I0, P([1,-1,-1])), 'A2', 'P=(1,-1,-1) A2 (A1) (3,0); LL', False); total += 1
  #passed += test_irrep(P3_A2_A1_3_0(L_L_s_I0, P([-1,1,-1])), 'A2', 'P=(-1,1,-1) A2 (A1) (3,0); LL', False); total += 1
  #passed += test_irrep(P3_A2_A1_3_0(L_L_s_I0, P([-1,-1,1])), 'A2', 'P=(-1,-1,1) A2 (A1) (3,0); LL', False); total += 1
  #passed += test_irrep(P3_A2_A1_3_0(L_L_s_I0, P([-1,-1,-1])), 'A2', 'P=(-1,-1,-1) A2 (A1) (3,0); LL', False); total += 1
  #passed += test_irrep(P3_A2_A1_3_0(N_X_a_I0, P([1,1,1])), 'A2', 'P=(1,1,1) A2 (A1) (3,0); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A2_A1_3_0(N_X_a_I0, P([1,1,-1])), 'A2', 'P=(1,1,-1) A2 (A1) (3,0); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A2_A1_3_0(N_X_a_I0, P([1,-1,1])), 'A2', 'P=(1,-1,1) A2 (A1) (3,0); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A2_A1_3_0(N_X_a_I0, P([-1,1,1])), 'A2', 'P=(-1,1,1) A2 (A1) (3,0); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A2_A1_3_0(N_X_a_I0, P([1,-1,-1])), 'A2', 'P=(1,-1,-1) A2 (A1) (3,0); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A2_A1_3_0(N_X_a_I0, P([-1,1,-1])), 'A2', 'P=(-1,1,-1) A2 (A1) (3,0); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A2_A1_3_0(N_X_a_I0, P([-1,-1,1])), 'A2', 'P=(-1,-1,1) A2 (A1) (3,0); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A2_A1_3_0(N_X_a_I0, P([-1,-1,-1])), 'A2', 'P=(-1,-1,-1) A2 (A1) (3,0); NXa I=0', False); total += 1

  #passed += test_irrep(P3_A2_A1_2_1(L_L_s_I0, P([1,1,1])), 'A2', 'P=(1,1,1) A2 (A1) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A2_A1_2_1(L_L_s_I0, P([1,1,-1])), 'A2', 'P=(1,1,-1) A2 (A1) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A2_A1_2_1(L_L_s_I0, P([1,-1,1])), 'A2', 'P=(1,-1,1) A2 (A1) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A2_A1_2_1(L_L_s_I0, P([-1,1,1])), 'A2', 'P=(-1,1,1) A2 (A1) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A2_A1_2_1(L_L_s_I0, P([1,-1,-1])), 'A2', 'P=(1,-1,-1) A2 (A1) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A2_A1_2_1(L_L_s_I0, P([-1,1,-1])), 'A2', 'P=(-1,1,-1) A2 (A1) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A2_A1_2_1(L_L_s_I0, P([-1,-1,1])), 'A2', 'P=(-1,-1,1) A2 (A1) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A2_A1_2_1(L_L_s_I0, P([-1,-1,-1])), 'A2', 'P=(-1,-1,-1) A2 (A1) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A2_A1_2_1(N_X_a_I0, P([1,1,1])), 'A2', 'P=(1,1,1) A2 (A1) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A2_A1_2_1(N_X_a_I0, P([1,1,-1])), 'A2', 'P=(1,1,-1) A2 (A1) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A2_A1_2_1(N_X_a_I0, P([1,-1,1])), 'A2', 'P=(1,-1,1) A2 (A1) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A2_A1_2_1(N_X_a_I0, P([-1,1,1])), 'A2', 'P=(-1,1,1) A2 (A1) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A2_A1_2_1(N_X_a_I0, P([1,-1,-1])), 'A2', 'P=(1,-1,-1) A2 (A1) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A2_A1_2_1(N_X_a_I0, P([-1,1,-1])), 'A2', 'P=(-1,1,-1) A2 (A1) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A2_A1_2_1(N_X_a_I0, P([-1,-1,1])), 'A2', 'P=(-1,-1,1) A2 (A1) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A2_A1_2_1(N_X_a_I0, P([-1,-1,-1])), 'A2', 'P=(-1,-1,-1) A2 (A1) (2,1); NXa I=0', False); total += 1

  #passed += test_irrep(P3_A2_E2_2_1(L_L_s_I0, P([1,1,1])), 'A2', 'P=(1,1,1) A2 (E2) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A2_E2_2_1(L_L_s_I0, P([1,1,-1])), 'A2', 'P=(1,1,-1) A2 (E2) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A2_E2_2_1(L_L_s_I0, P([1,-1,1])), 'A2', 'P=(1,-1,1) A2 (E2) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A2_E2_2_1(L_L_s_I0, P([-1,1,1])), 'A2', 'P=(-1,1,1) A2 (E2) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A2_E2_2_1(L_L_s_I0, P([1,-1,-1])), 'A2', 'P=(1,-1,-1) A2 (E2) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A2_E2_2_1(L_L_s_I0, P([-1,1,-1])), 'A2', 'P=(-1,1,-1) A2 (E2) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A2_E2_2_1(L_L_s_I0, P([-1,-1,1])), 'A2', 'P=(-1,-1,1) A2 (E2) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A2_E2_2_1(L_L_s_I0, P([-1,-1,-1])), 'A2', 'P=(-1,-1,-1) A2 (E2) (2,1); LL', False); total += 1
  #passed += test_irrep(P3_A2_E2_2_1(N_X_a_I0, P([1,1,1])), 'A2', 'P=(1,1,1) A2 (E2) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A2_E2_2_1(N_X_a_I0, P([1,1,-1])), 'A2', 'P=(1,1,-1) A2 (E2) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A2_E2_2_1(N_X_a_I0, P([1,-1,1])), 'A2', 'P=(1,-1,1) A2 (E2) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A2_E2_2_1(N_X_a_I0, P([-1,1,1])), 'A2', 'P=(-1,1,1) A2 (E2) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A2_E2_2_1(N_X_a_I0, P([1,-1,-1])), 'A2', 'P=(1,-1,-1) A2 (E2) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A2_E2_2_1(N_X_a_I0, P([-1,1,-1])), 'A2', 'P=(-1,1,-1) A2 (E2) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A2_E2_2_1(N_X_a_I0, P([-1,-1,1])), 'A2', 'P=(-1,-1,1) A2 (E2) (2,1); NXa I=0', False); total += 1
  #passed += test_irrep(P3_A2_E2_2_1(N_X_a_I0, P([-1,-1,-1])), 'A2', 'P=(-1,-1,-1) A2 (E2) (2,1); NXa I=0', False); total += 1


  # Psq = 4 (A1)
  #passed += test_irrep(P4_A1_A1_1_1(L_L_s_I0, P([0,0,2])), 'A1', 'P=(0,0,2) A1 (A1) (1,1); LL', False); total += 1
  #passed += test_irrep(P4_A1_A1_1_1(L_L_s_I0, P([0,0,-2])), 'A1', 'P=(0,0,-2) A1 (A1) (1,1); LL', False); total += 1
  #passed += test_irrep(P4_A1_A1_1_1(L_L_s_I0, P([0,2,0])), 'A1', 'P=(0,2,0) A1 (A1) (1,1); LL', False); total += 1
  #passed += test_irrep(P4_A1_A1_1_1(L_L_s_I0, P([0,-2,0])), 'A1', 'P=(0,-2,0) A1 (A1) (1,1); LL', False); total += 1
  #passed += test_irrep(P4_A1_A1_1_1(L_L_s_I0, P([2,0,0])), 'A1', 'P=(2,0,0) A1 (A1) (1,1); LL', False); total += 1
  #passed += test_irrep(P4_A1_A1_1_1(L_L_s_I0, P([-2,0,0])), 'A1', 'P=(-2,0,0) A1 (A1) (1,1); LL', False); total += 1

  # Psq = 4 (A2)
  #passed += test_irrep(P4_A2_A1_1_1(N_X_a_I0, P([0,0,2])), 'A2', 'P=(0,0,2) A2 (A1) (1,1); NXa I=0', False); total += 1
  #passed += test_irrep(P4_A2_A1_1_1(N_X_a_I0, P([0,0,-2])), 'A2', 'P=(0,0,-2) A2 (A1) (1,1); NXa I=0', False); total += 1
  #passed += test_irrep(P4_A2_A1_1_1(N_X_a_I0, P([0,2,0])), 'A2', 'P=(0,2,0) A2 (A1) (1,1); NXa I=0', False); total += 1
  #passed += test_irrep(P4_A2_A1_1_1(N_X_a_I0, P([0,-2,0])), 'A2', 'P=(0,-2,0) A2 (A1) (1,1); NXa I=0', False); total += 1
  #passed += test_irrep(P4_A2_A1_1_1(N_X_a_I0, P([2,0,0])), 'A2', 'P=(2,0,0) A2 (A1) (1,1); NXa I=0', False); total += 1
  #passed += test_irrep(P4_A2_A1_1_1(N_X_a_I0, P([-2,0,0])), 'A2', 'P=(-2,0,0) A2 (A1) (1,1); NXa I=0', False); total += 1

  # Psq = 4 (E1)
  #passed += test_irrep(P4_E1_A1_1_1(N_X_a_I0, P([0,0,2])), 'E', 'P=(0,0,2) E1 (A1) (1,1); NXa I=0', False); total += 1
  #passed += test_irrep(P4_E1_A1_1_1(N_X_a_I0, P([0,0,-2])), 'E', 'P=(0,0,-2) E1 (A1) (1,1); NXa I=0', False); total += 1
  #passed += test_irrep(P4_E1_A1_1_1(N_X_a_I0, P([0,2,0])), 'E', 'P=(0,2,0) E1 (A1) (1,1); NXa I=0', False); total += 1
  #passed += test_irrep(P4_E1_A1_1_1(N_X_a_I0, P([0,-2,0])), 'E', 'P=(0,-2,0) E1 (A1) (1,1); NXa I=0', False); total += 1
  #passed += test_irrep(P4_E1_A1_1_1(N_X_a_I0, P([2,0,0])), 'E', 'P=(2,0,0) E1 (A1) (1,1); NXa I=0', False); total += 1
  #passed += test_irrep(P4_E1_A1_1_1(N_X_a_I0, P([-2,0,0])), 'E', 'P=(-2,0,0) E1 (A1) (1,1); NXa I=0', False); total += 1




  print("Irrep tests: ({}/{}) PASSED".format(passed, total))
  if total != passed:
    print("\t{} FAILED!!".format(total-passed))


  # Test momentum equivalence
  
  #test_P1_A1(L_L_s_I0)  # PASSED
  #test_P1_A1(N_X_a_I0)  # PASSED
  #test_P1_A1(L_L_s_I0, N_X_a_I0)  # PASSED
  #test_P1_A2(L_L_s_I0)  # PASSED
  #test_P1_A2(N_X_a_I0)  # PASSED
  #test_P1_A2(L_L_s_I0, N_X_a_I0)  # PASSED

  #test_P1_B1(L_L_s_I0)  # PASSED
  #test_P1_B1(N_X_a_I0)  # PASSED
  #test_P1_B1(L_L_s_I0, N_X_a_I0)  # PASSED
  #test_P1_B2(L_L_s_I0)  # PASSED
  #test_P1_B2(N_X_a_I0)  # PASSED
  #test_P1_B2(L_L_s_I0, N_X_a_I0)  # PASSED
  #test_P1_E2(L_L_s_I0)  # PASSED
  #test_P1_E2(N_X_a_I0)  # PASSED
  #test_P1_E2(L_L_s_I0, N_X_a_I0)  # PASSED

  #test_P2_A1(L_L_s_I0, N_X_a_I0)  # PASSED
  #test_P2_A2(L_L_s_I0, N_X_a_I0)  # PASSED
  #test_P2_B1(L_L_s_I0, N_X_a_I0)  # PASSED
  #test_P2_B2(L_L_s_I0, N_X_a_I0)  # PASSED

  #test_P3_A1(L_L_s_I0)  # PASSED
  #test_P3_A1(L_L_s_I0, N_X_a_I0)  # PASSED
  #test_P3_A1(N_X_a_I0)  # PASSED
  #test_P3_A2(L_L_s_I0)  # PASSED
  #test_P3_A2(L_L_s_I0, N_X_a_I0)  # PASSED

if __name__ == "__main__":
  main()
