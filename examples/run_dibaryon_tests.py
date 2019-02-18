import time
import argparse

from sympy import Array

from context import operators
from operators.quarks import DiracIdx
from operators.octets import OctetBaryon
from operators.gamma import Gamma

from dibaryon_tests import *

def main():
  parser = argparse.ArgumentParser(description="Run Dibaryon Operator Checks")
  parser.add_argument("-q", "--quark-level", action="store_true", required=False,
                      help="Create Baryon operators at the quark level")

  args = parser.parse_args()

  g = Gamma()

  C5p = Array(g.chargeConj * g.five * g.parityPlus)
  C1p = Array(g.chargeConj * g.one * g.parityPlus)
  C2p = Array(g.chargeConj * g.two * g.parityPlus)
  C3p = Array(g.chargeConj * g.three * g.parityPlus)

  if args.quark_level:
    return NotImplemented
  else:
    p = OctetBaryon('p')
    n = OctetBaryon('n')

    Sp = OctetBaryon('S+')
    S0 = OctetBaryon('S0')
    Sm = OctetBaryon('S-')

    L = OctetBaryon('L')

    X0 = OctetBaryon('X0')
    Xm = OctetBaryon('X-')

  a = DiracIdx('a')
  b = DiracIdx('b')

  # I = 0 operators
  LL_I0_0 = L[a]*C5p[a,b]*L[b]
  LL_I0_1 = L[a]*C1p[a,b]*L[b]
  LL_I0_2 = L[a]*C2p[a,b]*L[b]
  LL_I0_3 = L[a]*C3p[a,b]*L[b]
  LL_I0 = [LL_I0_0, LL_I0_1, LL_I0_2, LL_I0_3]

  SS_I0_0 = Sp[a]*C5p[a,b]*Sm[b] - S0[a]*C5p[a,b]*S0[b] + Sm[a]*C5p[a,b]*Sp[b]
  SS_I0_1 = Sp[a]*C1p[a,b]*Sm[b] - S0[a]*C1p[a,b]*S0[b] + Sm[a]*C1p[a,b]*Sp[b]
  SS_I0_2 = Sp[a]*C2p[a,b]*Sm[b] - S0[a]*C2p[a,b]*S0[b] + Sm[a]*C2p[a,b]*Sp[b]
  SS_I0_3 = Sp[a]*C3p[a,b]*Sm[b] - S0[a]*C3p[a,b]*S0[b] + Sm[a]*C3p[a,b]*Sp[b]
  SS_I0 = [SS_I0_0, SS_I0_1, SS_I0_2, SS_I0_3]

  NXs_I0_0 = p[a]*C5p[a,b]*Xm[b] - n[a]*C5p[a,b]*X0[b] + Xm[a]*C5p[a,b]*p[b] - X0[a]*C5p[a,b]*n[b]
  NXs_I0_1 = p[a]*C1p[a,b]*Xm[b] - n[a]*C1p[a,b]*X0[b] + Xm[a]*C1p[a,b]*p[b] - X0[a]*C1p[a,b]*n[b]
  NXs_I0_2 = p[a]*C2p[a,b]*Xm[b] - n[a]*C2p[a,b]*X0[b] + Xm[a]*C2p[a,b]*p[b] - X0[a]*C2p[a,b]*n[b]
  NXs_I0_3 = p[a]*C3p[a,b]*Xm[b] - n[a]*C3p[a,b]*X0[b] + Xm[a]*C3p[a,b]*p[b] - X0[a]*C3p[a,b]*n[b]
  NXs_I0 = [NXs_I0_0, NXs_I0_1, NXs_I0_2, NXs_I0_3]

  NXa_I0_0 = p[a]*C5p[a,b]*Xm[b] - n[a]*C5p[a,b]*X0[b] - Xm[a]*C5p[a,b]*p[b] + X0[a]*C5p[a,b]*n[b]
  NXa_I0_1 = p[a]*C1p[a,b]*Xm[b] - n[a]*C1p[a,b]*X0[b] - Xm[a]*C1p[a,b]*p[b] + X0[a]*C1p[a,b]*n[b]
  NXa_I0_2 = p[a]*C2p[a,b]*Xm[b] - n[a]*C2p[a,b]*X0[b] - Xm[a]*C2p[a,b]*p[b] + X0[a]*C2p[a,b]*n[b]
  NXa_I0_3 = p[a]*C3p[a,b]*Xm[b] - n[a]*C3p[a,b]*X0[b] - Xm[a]*C3p[a,b]*p[b] + X0[a]*C3p[a,b]*n[b]
  NXa_I0 = [NXa_I0_0, NXa_I0_1, NXa_I0_2, NXa_I0_3]

  # I = 1 operators
  SLa_I1_0 = S0[a]*C5p[a,b]*L[b] - L[b]*C5p[a,b]*S0[b]
  SLa_I1_1 = S0[a]*C1p[a,b]*L[b] - L[b]*C1p[a,b]*S0[b]
  SLa_I1_2 = S0[a]*C2p[a,b]*L[b] - L[b]*C2p[a,b]*S0[b]
  SLa_I1_3 = S0[a]*C3p[a,b]*L[b] - L[b]*C3p[a,b]*S0[b]
  SLa_I1 = [SLa_I1_0, SLa_I1_1, SLa_I1_2, SLa_I1_3]

  SS_I1_0 = Sp[a]*C5p[a,b]*Sm[b] - Sm[a]*C5p[a,b]*Sp[b]
  SS_I1_1 = Sp[a]*C1p[a,b]*Sm[b] - Sm[a]*C1p[a,b]*Sp[b]
  SS_I1_2 = Sp[a]*C2p[a,b]*Sm[b] - Sm[a]*C2p[a,b]*Sp[b]
  SS_I1_3 = Sp[a]*C3p[a,b]*Sm[b] - Sm[a]*C3p[a,b]*Sp[b]
  SS_I1 = [SS_I1_0, SS_I1_1, SS_I1_2, SS_I1_3]

  NXa_I1_0 = p[a]*C5p[a,b]*Xm[b] + n[a]*C5p[a,b]*X0[b] - Xm[a]*C5p[a,b]*p[b] - X0[a]*C5p[a,b]*n[b]
  NXa_I1_1 = p[a]*C1p[a,b]*Xm[b] + n[a]*C1p[a,b]*X0[b] - Xm[a]*C1p[a,b]*p[b] - X0[a]*C1p[a,b]*n[b]
  NXa_I1_2 = p[a]*C2p[a,b]*Xm[b] + n[a]*C2p[a,b]*X0[b] - Xm[a]*C2p[a,b]*p[b] - X0[a]*C2p[a,b]*n[b]
  NXa_I1_3 = p[a]*C3p[a,b]*Xm[b] + n[a]*C3p[a,b]*X0[b] - Xm[a]*C3p[a,b]*p[b] - X0[a]*C3p[a,b]*n[b]
  NXa_I1 = [NXa_I1_0, NXa_I1_1, NXa_I1_2, NXa_I1_3]

  SLs_I1_0 = S0[a]*C5p[a,b]*L[b] + L[b]*C5p[a,b]*S0[b]
  SLs_I1_1 = S0[a]*C1p[a,b]*L[b] + L[b]*C1p[a,b]*S0[b]
  SLs_I1_2 = S0[a]*C2p[a,b]*L[b] + L[b]*C2p[a,b]*S0[b]
  SLs_I1_3 = S0[a]*C3p[a,b]*L[b] + L[b]*C3p[a,b]*S0[b]
  SLs_I1 = [SLs_I1_0, SLs_I1_1, SLs_I1_2, SLs_I1_3]

  NXs_I1_0 = p[a]*C5p[a,b]*Xm[b] + n[a]*C5p[a,b]*X0[b] + Xm[a]*C5p[a,b]*p[b] + X0[a]*C5p[a,b]*n[b]
  NXs_I1_1 = p[a]*C1p[a,b]*Xm[b] + n[a]*C1p[a,b]*X0[b] + Xm[a]*C1p[a,b]*p[b] + X0[a]*C1p[a,b]*n[b]
  NXs_I1_2 = p[a]*C2p[a,b]*Xm[b] + n[a]*C2p[a,b]*X0[b] + Xm[a]*C2p[a,b]*p[b] + X0[a]*C2p[a,b]*n[b]
  NXs_I1_3 = p[a]*C3p[a,b]*Xm[b] + n[a]*C3p[a,b]*X0[b] + Xm[a]*C3p[a,b]*p[b] + X0[a]*C3p[a,b]*n[b]
  NXs_I1 = [NXs_I1_0, NXs_I1_1, NXs_I1_2, NXs_I1_3]

  # Start irrep tests

  # P = 0 (A1p)
  P0_A1p.test_irrep(LL_I0, momenta=[0,1,2,3], message='LL (I=0)')
  P0_A1p.test_irrep(NXa_I0, message='NXa (I=0) - should fail', count=False)

  # P = 0 (A2p)
  P0_A2p.test_irrep(NXa_I0, message='NXa (I=0)')
  P0_A2p.test_irrep(LL_I0, message='LL (I=0) - should fail', count=False)

  # P = 0 (Ep) - LL, I=0
  P0_Ep_Ep_1.test_irrep(LL_I0, message='LL (I=0)')
  P0_Ep_Ep_1.test_irrep(NXa_I0, message='NXa (I=0) - should fail', count=False)

  P0_Ep_Ep_2.test_irrep(LL_I0, message='LL (I=0)')
  P0_Ep_Ep_2.test_irrep(NXa_I0, message='NXa (I=0) - should fail', count=False)

  P0_Ep_T2p.test_irrep(NXa_I0, message='NXa (I=0)')
  P0_Ep_T2p.test_irrep(LL_I0, message='LL (I=0) - should fail', count=False)

  # P = 0 (T1g)
  P0_T1p_A1p.test_irrep(NXa_I0, momenta=[0,1,2,3], message='NXa (I=0)')
  P0_T1p_A1p.test_irrep(LL_I0, message='LL (I=0) - should fail') # FAILS due to Pow not being handled...!

  P0_T1p_Ep_1.test_irrep(NXa_I0, message='NXa (I=0)')
  P0_T1p_Ep_1.test_irrep(LL_I0, message='LL (I=0) - should fail', count=False)

  P0_T1p_Ep_2.test_irrep(NXa_I0, message='NXa (I=0)')
  P0_T1p_Ep_2.test_irrep(LL_I0, message='LL (I=0) - should fail', count=False)

  P0_T1p_T2p_2.test_irrep(NXa_I0, message='NXa (I=0)')
  P0_T1p_T2p_2.test_irrep(LL_I0, message='LL (I=0) - should fail', count=False)

  P0_T1p_1.test_irrep(NXa_I0, message='NXa (I=0)')
  P0_T1p_1.test_irrep(LL_I0, message='LL (I=0) - should fail', count=False)

  P0_T1p_3.test_irrep(NXa_I0, message='NXa (I=0)')
  P0_T1p_3.test_irrep(LL_I0, message='LL (I=0) - should fail', count=False)

  # P = 0 (T2g)
  P0_T2p_Ep_1.test_irrep(NXa_I0, message='NXa (I=0)')
  P0_T2p_Ep_1.test_irrep(LL_I0, message='LL (I=0) - should fail', count=False)

  P0_T2p_T2p_2_0.test_irrep(LL_I0, message='LL (I=0)')
  P0_T2p_T2p_2_0.test_irrep(NXa_I0, message='NXa (I=0) - should fail', count=False)

  P0_T2p_Ep_2.test_irrep(NXa_I0, message='NXa (I=0)')
  P0_T2p_Ep_2.test_irrep(LL_I0, message='LL (I=0) - should fail', count=False)

  P0_T2p_T2p_2_1.test_irrep(NXa_I0, message='NXa (I=0)')
  P0_T2p_T2p_2_1.test_irrep(LL_I0, message='LL (I=0) - should fail', count=False)

  P0_T2p_2.test_irrep(NXa_I0, message='NXa (I=0)')
  P0_T2p_2.test_irrep(LL_I0, message='LL (I=0) - should fail', count=False)

  P0_T2p_3.test_irrep(NXa_I0, message='NXa (I=0)')
  P0_T2p_3.test_irrep(LL_I0, message='LL (I=0) - should fail', count=False)

  # P = 0 (A1m)
  P0_A1m.test_irrep(LL_I0, momenta=[1,2,3], message='LL (I=0)')
  P0_A1m.test_irrep(NXa_I0, message='NXa (I=0) - should fail', count=False)

  # P = 0 (A2m)
  P0_A2m.test_irrep(LL_I0, message='LL (I=0)')
  P0_A2m.test_irrep(NXa_I0, message='NXa (I=0) - should fail', count=False)

  # P = 0 (Em)
  P0_Em_T1m.test_irrep(LL_I0, momenta=[1,2,3], message='LL (I=0)')
  P0_Em_T1m.test_irrep(NXa_I0, message='NXa (I=0) - should fail', count=False)

  P0_Em_T2m.test_irrep(LL_I0, message='LL (I=0)')
  P0_Em_T2m.test_irrep(NXa_I0, message='NXa (I=0) - should fail', count=False)


  # P = 0 (T1m)
  P0_T1m_T1m_0.test_irrep(NXa_I0, momenta=[1,2,3], message='NXa (I=0)')
  P0_T1m_T1m_0.test_irrep(LL_I0, message='LL (I=0) - should fail', count=False)

  P0_T1m_T1m_1.test_irrep(LL_I0, momenta=[1,2,3], message='LL (I=0)')
  P0_T1m_T1m_1.test_irrep(NXa_I0, message='NXa (I=0) - should fail', count=False)

  P0_T1m_T2m.test_irrep(LL_I0, message='LL (I=0)')
  P0_T1m_T2m.test_irrep(NXa_I0, message='NXa (I=0) - should fail', count=False)

  

  # P = 0 (T2m)
  P0_T2m_T1m.test_irrep(LL_I0, momenta=[1,2,3], message='LL (I=0)')
  P0_T2m_T1m.test_irrep(NXa_I0, message='NXa (I=0) - should fail', count=False)

  P0_T2m_T2m_1.test_irrep(LL_I0, message='LL (I=0)')
  P0_T2m_T2m_1.test_irrep(NXa_I0, message='Nxa (I=0) - should fail', count=False)

  P0_T2m_T2m_0.test_irrep(NXa_I0, message='Nxa (I=0)')
  P0_T2m_T2m_0.test_irrep(LL_I0, message='Nxa (I=0) - should fail', count=False)



  # Psq = 1 (A1)
  P1_A1_LA1_10.test_irrep(LL_I0, message='LL (I=0)')
  P1_A1_LA1_10.test_irrep(NXa_I0, message='NXa (I=0)')

  P1_A1_LA1_21.test_irrep(LL_I0, message='LL (I=0)')
  P1_A1_LA1_21.test_irrep(NXa_I0, message='NXa (I=0)')

  P1_A1_LE2_21.test_irrep(LL_I0, message='LL (I=0)')
  P1_A1_LE2_21.test_irrep(NXa_I0, message='NXa (I=0)')

  # Psq = 1 (A2)
  P1_A2_LA1_10.test_irrep(LL_I0, message='LL (I=0)')
  P1_A2_LA1_10.test_irrep(NXa_I0, message='NXa (I=0)')

  P1_A2_LA1_21.test_irrep(LL_I0, message='LL (I=0)')
  P1_A2_LA1_21.test_irrep(NXa_I0, message='NXa (I=0)')

  P1_A2_LE2_21.test_irrep(LL_I0, message='LL (I=0)')
  P1_A2_LE2_21.test_irrep(NXa_I0, message='NXa (I=0)')


  # Psq = 1 (B1)
  P1_B1_LB1_21.test_irrep(LL_I0, message='LL (I=0)')
  P1_B1_LB1_21.test_irrep(NXa_I0, message='NXa (I=0)')

  P1_B1_LE2_21.test_irrep(LL_I0, message='LL (I=0)')
  P1_B1_LE2_21.test_irrep(NXa_I0, message='NXa (I=0)')

  # Psq = 1 (B2)
  P1_B2_LB1_21.test_irrep(LL_I0, message='LL (I=0)')
  P1_B2_LB1_21.test_irrep(NXa_I0, message='NXa (I=0)')

  P1_B2_LE2_21.test_irrep(LL_I0, message='LL (I=0)')
  P1_B2_LE2_21.test_irrep(NXa_I0, message='NXa (I=0)')

  # Psq = 1 (E2)
  P1_E2_LA1_10.test_irrep(LL_I0, message='LL (I=0)')
  P1_E2_LA1_10.test_irrep(NXa_I0, message='NXa (I=0)')

  P1_E2_LE2_21_0.test_irrep(LL_I0, message='LL (I=0)')
  P1_E2_LE2_21_0.test_irrep(NXa_I0, message='NXa (I=0)')

  P1_E2_LA1_21.test_irrep(LL_I0, message='LL (I=0)')
  P1_E2_LA1_21.test_irrep(NXa_I0, message='NXa (I=0)')

  P1_E2_LB1_21.test_irrep(LL_I0, message='LL (I=0)')
  P1_E2_LB1_21.test_irrep(NXa_I0, message='NXa (I=0)')

  P1_E2_LE2_21_1.test_irrep(LL_I0, message='LL (I=0)')
  P1_E2_LE2_21_1.test_irrep(NXa_I0, message='NXa (I=0)')

  # Psq = 2 (A1)
  P2_A1_LA1_20.test_irrep(LL_I0, message='LL (I=0)')
  P2_A1_LA1_20.test_irrep(NXa_I0, message='NXa (I=0)')

  P2_A1_LA1_11.test_irrep(LL_I0, message='LL (I=0)')
  P2_A1_LA1_11.test_irrep(NXa_I0, message='NXa (I=0) - should fail', count=False)

  P2_A1_LB1_11.test_irrep(LL_I0, message='LL (I=0)')
  P2_A1_LB1_11.test_irrep(NXa_I0, message='NXa (I=0) - should fail', count=False)

  # Psq = 2 (A2)
  P2_A2_LA1_20.test_irrep(LL_I0, message='LL (I=0)')
  P2_A2_LA1_20.test_irrep(NXa_I0, message='NXa (I=0)')

  P2_A2_LB1_11.test_irrep(LL_I0, message='LL (I=0)')
  P2_A2_LB1_11.test_irrep(NXa_I0, message='NXa (I=0) - should fail', count=False)

  P2_A2_LA1_11.test_irrep(LL_I0, message='LL (I=0) - should fail', count=False)
  P2_A2_LA1_11.test_irrep(NXa_I0, message='NXa (I=0)')

  # Psq = 2 (B1)
  P2_B1_LA1_20.test_irrep(LL_I0, message='LL (I=0)')
  P2_B1_LA1_20.test_irrep(LL_I0, message='LL (I=0)')

  P2_B1_LA1_11.test_irrep(LL_I0, message='LL (I=0) - should fail', count=False)
  P2_B1_LA1_11.test_irrep(NXa_I0, message='NXa (I=0)')

  P2_B1_LB1_11.test_irrep(LL_I0, message='LL (I=0) - should fail', count=False)
  P2_B1_LB1_11.test_irrep(NXa_I0, message='NXa (I=0)')

  # Psq = 2 (B2)
  P2_B2_LA1_20.test_irrep(LL_I0, message='LL (I=0)')
  P2_B2_LA1_20.test_irrep(NXa_I0, message='NXa (I=0)')

  P2_B2_LB1_11.test_irrep(LL_I0, message='LL (I=0)')
  P2_B2_LB1_11.test_irrep(NXa_I0, message='NXa (I=0) - should fail', count=False)

  P2_B2_LA1_11.test_irrep(LL_I0, message='LL (I=0) - should fail', count=False)
  P2_B2_LA1_11.test_irrep(NXa_I0, message='NXa (I=0)')

  # Psq = 3 (A1)
  P3_A1_LA1_30.test_irrep(LL_I0, message='LL (I=0)')
  P3_A1_LA1_30.test_irrep(NXa_I0, message='NXa (I=0)')
  
  P3_A1_LA1_21.test_irrep(LL_I0, message='LL (I=0)')
  P3_A1_LA1_21.test_irrep(NXa_I0, message='NXa (I=0)')

  P3_A1_LE2_21.test_irrep(LL_I0, message='LL (I=0)')
  P3_A1_LE2_21.test_irrep(NXa_I0, message='NXa (I=0)')

  # Psq = 3 (A2)
  P3_A2_LA1_30.test_irrep(LL_I0, message='LL (I=0)')
  P3_A2_LA1_30.test_irrep(NXa_I0, message='NXa (I=0)')

  P3_A2_LA1_21.test_irrep(LL_I0, message='LL (I=0)')
  P3_A2_LA1_21.test_irrep(NXa_I0, message='NXa (I=0)')

  P3_A2_LE2_21.test_irrep(LL_I0, message='LL (I=0)')
  P3_A2_LE2_21.test_irrep(NXa_I0, message='NXa (I=0)')

  # Psq = 3 (E2)
  P3_E2_LA1_30.test_irrep(LL_I0, message='LL (I=0)')
  P3_E2_LA1_30.test_irrep(NXa_I0, message='NXa (I=0)')

  P3_E2_LA1_21.test_irrep(LL_I0, message='LL (I=0)')
  P3_E2_LA1_21.test_irrep(NXa_I0, message='NXa (I=0)')

  P3_E2_LE2_21_0.test_irrep(LL_I0, message='LL (I=0)')
  P3_E2_LE2_21_0.test_irrep(NXa_I0, message='NXa (I=0)')

  P3_E2_LE2_21_1L.test_irrep(LL_I0, message='LL (I=0)')
  P3_E2_LE2_21_1L.test_irrep(NXa_I0, message='NXa (I=0)')

  P3_E2_LE2_21_1T.test_irrep(LL_I0, message='LL (I=0)')
  P3_E2_LE2_21_1T.test_irrep(NXa_I0, message='NXa (I=0)')

  # Psq = 4 (A1)
  P4_A1_LA1_11.test_irrep(LL_I0, message='LL (I=0)')
  P4_A1_LA1_11.test_irrep(NXa_I0, message='NXa (I=0) - should fail', count=False)

  # Psq = 4 (A2)
  P4_A2_LA1_11.test_irrep(LL_I0, message='LL (I=0) - should fail', count=False)
  P4_A2_LA1_11.test_irrep(NXa_I0, message='NXa (I=0)')

  # Psq = 4 (E1)
  P4_E2_LA1_11.test_irrep(LL_I0, message='LL (I=0) - should fail', count=False)
  P4_E2_LA1_11.test_irrep(NXa_I0, message='NXa (I=0)')


  # Test momentum equivalence
  OperatorTest.test_equivalent_transforms(P1_A1_LA1_10.getMomRep(LL_I0),
                                          P1_A1_LA1_10.getMomRep(NXa_I0),
                                          P1_A1_LA1_21.getMomRep(LL_I0),
                                          P1_A1_LA1_21.getMomRep(NXa_I0),
                                          P1_A1_LE2_21.getMomRep(LL_I0),
                                          P1_A1_LE2_21.getMomRep(NXa_I0),
                                          message='P1 A1')

  OperatorTest.test_equivalent_transforms(P1_A2_LA1_10.getMomRep(LL_I0),
                                          P1_A2_LA1_10.getMomRep(NXa_I0),
                                          P1_A2_LA1_21.getMomRep(LL_I0),
                                          P1_A2_LA1_21.getMomRep(NXa_I0),
                                          P1_A2_LE2_21.getMomRep(LL_I0),
                                          P1_A2_LE2_21.getMomRep(NXa_I0),
                                          message='P1 A2')

  OperatorTest.test_equivalent_transforms(P1_B1_LB1_21.getMomRep(LL_I0),
                                          P1_B1_LB1_21.getMomRep(NXa_I0),
                                          P1_B1_LE2_21.getMomRep(LL_I0),
                                          P1_B1_LE2_21.getMomRep(NXa_I0),
                                          message='P1 B1')

  OperatorTest.test_equivalent_transforms(P1_B2_LB1_21.getMomRep(LL_I0),
                                          P1_B2_LB1_21.getMomRep(NXa_I0),
                                          P1_B2_LE2_21.getMomRep(LL_I0),
                                          P1_B2_LE2_21.getMomRep(NXa_I0),
                                          message='P1 B2')

  OperatorTest.test_equivalent_transforms(P1_E2_LA1_10.getMomRep(LL_I0),
                                          P1_E2_LA1_10.getMomRep(NXa_I0),
                                          P1_E2_LE2_21_0.getMomRep(LL_I0),
                                          P1_E2_LE2_21_0.getMomRep(NXa_I0),
                                          P1_E2_LA1_21.getMomRep(LL_I0),
                                          P1_E2_LA1_21.getMomRep(NXa_I0),
                                          P1_E2_LB1_21.getMomRep(LL_I0),
                                          P1_E2_LB1_21.getMomRep(NXa_I0),
                                          P1_E2_LE2_21_1.getMomRep(LL_I0),
                                          P1_E2_LE2_21_1.getMomRep(NXa_I0),
                                          message='P1 E2')

  OperatorTest.test_equivalent_transforms(P2_A1_LA1_20.getMomRep(LL_I0),
                                          P2_A1_LA1_20.getMomRep(NXa_I0),
                                          P2_A1_LA1_11.getMomRep(LL_I0),
                                          P2_A1_LB1_11.getMomRep(LL_I0),
                                          message='P2 A1')


  OperatorTest.test_equivalent_transforms(P2_A2_LA1_20.getMomRep(LL_I0),
                                          P2_A2_LA1_20.getMomRep(NXa_I0),
                                          P2_A2_LB1_11.getMomRep(LL_I0),
                                          P2_A2_LA1_11.getMomRep(NXa_I0),
                                          message='P2 A2')


  OperatorTest.test_equivalent_transforms(P2_B1_LA1_20.getMomRep(LL_I0),
                                          P2_B1_LA1_20.getMomRep(LL_I0),
                                          P2_B1_LA1_11.getMomRep(NXa_I0),
                                          P2_B1_LB1_11.getMomRep(NXa_I0),
                                          message='P2 B1')

  OperatorTest.test_equivalent_transforms(P2_B2_LA1_20.getMomRep(LL_I0),
                                          P2_B2_LA1_20.getMomRep(NXa_I0),
                                          P2_B2_LB1_11.getMomRep(LL_I0),
                                          P2_B2_LA1_11.getMomRep(NXa_I0),
                                          message='P2 B2')

  OperatorTest.test_equivalent_transforms(P3_A1_LA1_30.getMomRep(LL_I0),
                                          P3_A1_LA1_30.getMomRep(NXa_I0),
                                          P3_A1_LA1_21.getMomRep(LL_I0),
                                          P3_A1_LA1_21.getMomRep(NXa_I0),
                                          P3_A1_LE2_21.getMomRep(LL_I0),
                                          P3_A1_LE2_21.getMomRep(NXa_I0),
                                          message='P3 A1')

  OperatorTest.test_equivalent_transforms(P3_A2_LA1_30.getMomRep(LL_I0),
                                          P3_A2_LA1_30.getMomRep(NXa_I0),
                                          P3_A2_LA1_21.getMomRep(LL_I0),
                                          P3_A2_LA1_21.getMomRep(NXa_I0),
                                          P3_A2_LE2_21.getMomRep(LL_I0),
                                          P3_A2_LE2_21.getMomRep(NXa_I0),
                                          message='P3 A2')


  OperatorTest.test_equivalent_transforms(P3_E2_LA1_30.getMomRep(LL_I0),
                                          P3_E2_LA1_30.getMomRep(NXa_I0),
                                          P3_E2_LA1_21.getMomRep(LL_I0),
                                          P3_E2_LA1_21.getMomRep(NXa_I0),
                                          P3_E2_LE2_21_0.getMomRep(LL_I0),
                                          P3_E2_LE2_21_0.getMomRep(NXa_I0),
                                          P3_E2_LE2_21_1L.getMomRep(LL_I0),
                                          P3_E2_LE2_21_1L.getMomRep(NXa_I0),
                                          P3_E2_LE2_21_1T.getMomRep(LL_I0),
                                          P3_E2_LE2_21_1T.getMomRep(NXa_I0),
                                          message='P3 E2')

  OperatorTest.print_results()
  

if __name__ == "__main__":
  start = time.time()
  main()
  total_time = time.time() - start

  print("TOTAL TIME: {} sec".format(total_time))
