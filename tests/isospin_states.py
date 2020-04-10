from sympy import Rational, sqrt

import single_hadrons as sh

from context import operators

from operators.operators import Operator

def isoseptet_pion_pion_pion(mom1, mom2, mom3):
  return Operator(sh.pion_p, mom1)*Operator(sh.pion_p, mom2)*Operator(sh.pion_p, mom3)

def isoquartet_kaon_kaon_kaon(mom1, mom2, mom3):
  return Operator(sh.kaon_p, mom1)*Operator(sh.kaon_p, mom2)*Operator(sh.kaon_p, mom3)

def isotriplet_kaon_kaon_pion_one(mom1, mom2, mom3):
  K0_Kp_pip = Operator(sh.kaon_0, mom1)*Operator(sh.kaon_p, mom2)*Operator(sh.pion_p, mom3)
  Kp_K0_pip = Operator(sh.kaon_p, mom1)*Operator(sh.kaon_0, mom2)*Operator(sh.pion_p, mom3)
  Kp_Kp_pi0 = Operator(sh.kaon_p, mom1)*Operator(sh.kaon_p, mom2)*Operator(sh.pion_0, mom3)

  return Rational(1,2)*K0_Kp_pip + Rational(1,2)*Kp_K0_pip - Rational(1,sqrt(2))*Kp_Kp_pi0

def isotriplet_kaon_kaon_pion_two(mom1, mom2, mom3):
  K0_Kp_pip = Operator(sh.kaon_0, mom1)*Operator(sh.kaon_p, mom2)*Operator(sh.pion_p, mom3)
  Kp_K0_pip = Operator(sh.kaon_p, mom1)*Operator(sh.kaon_0, mom2)*Operator(sh.pion_p, mom3)

  return sqrt(Rational(1,2))*(K0_Kp_pip - Kp_K0_pip)

def isoquartet_kaon_kaon_kaon_pion_one(mom1, mom2, mom3, mom4):
  K0_Kp_Kp_pip = Operator(sh.kaon_0, mom1)*Operator(sh.kaon_p, mom2)*Operator(sh.kaon_p, mom3)*Operator(sh.pion_p, mom4)
  Kp_K0_Kp_pip = Operator(sh.kaon_p, mom1)*Operator(sh.kaon_0, mom2)*Operator(sh.kaon_p, mom3)*Operator(sh.pion_p, mom4)
  Kp_Kp_K0_pip = Operator(sh.kaon_p, mom1)*Operator(sh.kaon_p, mom2)*Operator(sh.kaon_0, mom3)*Operator(sh.pion_p, mom4)
  Kp_Kp_Kp_pi0 = Operator(sh.kaon_p, mom1)*Operator(sh.kaon_p, mom2)*Operator(sh.kaon_p, mom3)*Operator(sh.pion_0, mom4)

  return -sqrt(Rational(2,15))*K0_Kp_Kp_pip - sqrt(Rational(2,15))*Kp_K0_Kp_pip - sqrt(Rational(2,15))*Kp_Kp_K0_pip + sqrt(Rational(3,5))*Kp_Kp_Kp_pi0

def isoquartet_kaon_kaon_kaon_pion_two(mom1, mom2, mom3, mom4):
  K0_Kp_Kp_pip = Operator(sh.kaon_0, mom1)*Operator(sh.kaon_p, mom2)*Operator(sh.kaon_p, mom3)*Operator(sh.pion_p, mom4)
  Kp_K0_Kp_pip = Operator(sh.kaon_p, mom1)*Operator(sh.kaon_0, mom2)*Operator(sh.kaon_p, mom3)*Operator(sh.pion_p, mom4)
  Kp_Kp_K0_pip = Operator(sh.kaon_p, mom1)*Operator(sh.kaon_p, mom2)*Operator(sh.kaon_0, mom3)*Operator(sh.pion_p, mom4)

  return sqrt(Rational(1,6))*(-2*K0_Kp_Kp_pip + Kp_K0_Kp_pip + Kp_Kp_K0_pip)

def isoquartet_kaon_kaon_kaon_pion_three(mom1, mom2, mom3, mom4):
  Kp_K0_Kp_pip = Operator(sh.kaon_p, mom1)*Operator(sh.kaon_0, mom2)*Operator(sh.kaon_p, mom3)*Operator(sh.pion_p, mom4)
  Kp_Kp_K0_pip = Operator(sh.kaon_p, mom1)*Operator(sh.kaon_p, mom2)*Operator(sh.kaon_0, mom3)*Operator(sh.pion_p, mom4)

  return Rational(1,2)*Kp_K0_Kp_pip - Rational(1,2)*Kp_Kp_K0_pip

def isosinglet_pion_pion_pion(mom1, mom2, mom3):
  pim_pi0_pip = Operator(sh.pion_m, mom1)*Operator(sh.pion_0, mom2)*Operator(sh.pion_p, mom3)
  pi0_pim_pip = Operator(sh.pion_0, mom1)*Operator(sh.pion_m, mom2)*Operator(sh.pion_p, mom3)
  pim_pip_pi0 = Operator(sh.pion_m, mom1)*Operator(sh.pion_p, mom2)*Operator(sh.pion_0, mom3)
  pip_pim_pi0 = Operator(sh.pion_p, mom1)*Operator(sh.pion_m, mom2)*Operator(sh.pion_0, mom3)
  pi0_pip_pim = Operator(sh.pion_0, mom1)*Operator(sh.pion_p, mom2)*Operator(sh.pion_m, mom3)
  pip_pi0_pim = Operator(sh.pion_p, mom1)*Operator(sh.pion_0, mom2)*Operator(sh.pion_m, mom3)

  return sqrt(Rational(1,6))*(-pim_pi0_pip + pi0_pim_pip + pim_pip_pi0 - pip_pim_pi0 - pi0_pip_pim + pip_pi0_pim)
