import time
import cProfile
from typing import NamedTuple
from types import FunctionType
from pprint import pprint

from sympy import S, sqrt, Eijk, KroneckerDelta, Rational, simplify, expand

from context import operators
from operators.operators import project_momentum
from operators.operator_rep import OperatorRepresentation, OperatorRepresentationError
from operators.cubic_rotations import _POINT_GROUP, P, P0, MOMENTA

####################################################################################################
#                                                                                                  #
#               OperatorTest Class                                                                 #
#                                                                                                  #
####################################################################################################


class OperatorTest(NamedTuple):
  name: str
  irrep: str
  psq: int
  operator_def: FunctionType

  PASSED = 0
  TOTAL_TESTS = 0
  TOTAL_TIME = 0

  def test_irrep(self, f_op, momenta=None, message='', count=True):
    start = time.time()
    print("TESTING: {}".format(self.name))
    if message:
      print("\tINFO: {}".format(message))

    if momenta is not None:
      print("\tMOMENTA: {}".format(momenta))
    elif momenta is None and self.psq:
      momenta = MOMENTA[self.psq]
      print("\tMOMENTA: {}".format(momenta))
    else:
      print("\tMOMENTA: none given, assuming one momentum")

    print()

    t_irreps = set()
    try:
      if momenta is not None:
        for momentum in momenta:
          print("\tTESTING MOMENTUM: {}".format(momentum))
          op_rep = OperatorRepresentation(*self.operator_def(f_op, momentum))
          t_irrep = op_rep.lgIrrepOccurences(True)
          t_irreps.add(t_irrep)
          print("\t\tRep: {}\n".format(t_irrep))
      else:
        op_rep = OperatorRepresentation(*self.operator_def(f_op))
        t_irrep = op_rep.lgIrrepOccurences(True)
        t_irreps.add(t_irrep)

    except Exception as e:
      print("\tFAILED: {}\n".format(e))
      OperatorTest.TOTAL_TESTS += count
      return

    if len(t_irreps) != 1:
      print("\tFAILED: Not all irreps in set transformed identically")
    else:
      t_irrep = t_irreps.pop()

      if t_irrep == self.irrep:
        OperatorTest.PASSED += count
        print("\tPASSED:")
        print("\t\tExpected Irrep: {}".format(self.irrep))
        print("\t\tActual Rep:     {}".format(t_irrep))

      else:
        print("\tFAILED:")
        print("\t\tExpected Irrep: {}".format(self.irrep))
        print("\t\tActual Rep:     {}".format(t_irrep))

    total = time.time() - start
    OperatorTest.TOTAL_TIME += total
    print("\n\tTIME: {}\n".format(total))

    OperatorTest.TOTAL_TESTS += count

  def getMomRep(self, f_op):
    if self.psq == 0:
      raise ValueError("Mometum Representation requires moving operator")

    ops = []
    for p in MOMENTA[self.psq]:
      ops.extend(self.operator_def(f_op, p))

    return OperatorRepresentation(*ops)

  @staticmethod
  def test_equivalent_transforms(*mom_reps, message=''):
    start = time.time()
    print("EQUIVALENT MOMENTUM TESTING:")
    if message:
      print("\tINFO: {}".format(message))

    failed = False
    # @ADH - Could just choose a representative element for each equivalent frame rather than the whole _POINT_GROUP
    for el in _POINT_GROUP:
      mats = list()
      for mom_rep in mom_reps:
        try:
          # expand doesn't seem to work correctly, and this appears faster than just calling simplify
          mats.append(mom_rep.getRepresentationMatrix(el).applyfunc(expand).applyfunc(simplify))
        except Exception as e:
          print("\tFAILED: {}\n".format(e))
          OperatorTest.TOTAL_TESTS += 1
          return

      if not all(mat == mats[0] for mat in mats):
        print("\tFAILED: not all operators transform in the same way")
        print("\t\tELEMENT: {}".format(el))
        for mat in mats:
          print()
          pprint(mat)

        failed = True
        break

    if not failed:
      OperatorTest.PASSED += 1
      print("\tPASSED")

    total = time.time() - start
    OperatorTest.TOTAL_TIME += total
    print("\n\tTIME: {}\n".format(total))

    OperatorTest.TOTAL_TESTS += 1

  @staticmethod
  def print_results():
    print("{}/{} PASSED".format(OperatorTest.PASSED, OperatorTest.TOTAL_TESTS))
    failed = OperatorTest.TOTAL_TESTS - OperatorTest.PASSED
    if failed:
      print("{} FAILED!!!".format(failed))

    print("TOTAL TEST TIME: {}".format(OperatorTest.TOTAL_TIME))

    OperatorTest.PASSED = 0
    OperatorTest.TOTAL_TESTS = 0
    OperatorTest.TOTAL_TIME = 0





####################################################################################################
#                                                                                                  #
#               Operator Definitions                                                               #
#                                                                                                  #
####################################################################################################

############ HELPERS ###############################################################################

L = {
    1: { 1: 1/sqrt(2), 2: -1/sqrt(2), 3:         0 },
    2: { 1: 1/sqrt(6), 2:  1/sqrt(6), 3: -2/sqrt(6)}
}

Pi = {
    1: P(1,0,0),
    2: P(0,1,0),
    3: P(0,0,1)
}


def mod3(i):
  return ((i - 1)%3) + 1

def Pks(p):
  pks = {}
  n = 1
  for i, pi in Pi.items():
    if p*pi != P0:
      pks[n] = pi
      n += 1

  return pks

def Pksi(p):
  pks = {}
  n = 1
  for i, pi in Pi.items():
    if p*pi != P0:
      pks[n] = i
      n += 1

  return pks

def Pc(p):
  ds = {}
  n = 1

  if p.x != 0:
    ds[n] = P(p.x,0,0)
    n += 1

  if p.y != 0:
    ds[n] = P(0, p.y,0)
    n += 1

  if p.z != 0:
    ds[n] = P(0,0,p.z)
    n += 1

  return ds

#############  P^2 = 0  ############################################################################
#############    A1p    ############################################################################

def _P0_A1p(f_ops, n=0):

  op = S.Zero
  N = 0
  for p in MOMENTA[n]:
    op += project_momentum(f_ops[0], p, -p)
    N += 1

  return [op/N]

#############    A2p    ############################################################################
def _P0_A2p_LT2p(f_ops):

  op = S.Zero
  for i in range(1,4):
    p1 = Pi[mod3(i+1)] + Pi[mod3(i+2)]
    p2 = Pi[mod3(i+1)] - Pi[mod3(i+2)]

    op += project_momentum(f_ops[i], p1, -p1) - project_momentum(f_ops[i], p2, -p2)

  return [op]

def _P0_LT2p_A2p(f_ops):

  op = S.Zero
  for i in range(1,4):
    p1 = Pi[mod3(i+1)] + Pi[mod3(i+2)]
    p2 = Pi[mod3(i+1)] - Pi[mod3(i+2)]

    op += project_momentum(f_ops[0], p1, -p1) - project_momentum(f_ops[0], p2, -p2)

  return [op]

#############    Ep     ############################################################################
def _P0_Ep_LEp_1(f_ops):

  ops = [S.Zero, S.Zero]
  for i in range(1,3):
    for j in range(1,4):
      ops[i-1] += L[i][j]*project_momentum(f_ops[0], Pi[j], -Pi[j])

  return ops

def _P0_Ep_LEp_2(f_ops):

  ops = [S.Zero, S.Zero]
  for i in range(1,3):
    for j in range(1,4):
      p1 = Pi[mod3(j+1)] + Pi[mod3(j+2)]
      p2 = Pi[mod3(j+1)] - Pi[mod3(j+2)]

      ops[i-1] += L[i][j]*(project_momentum(f_ops[0], p1, -p1) \
                         + project_momentum(f_ops[0], p2, -p2))

  return ops

def _P0_Ep_LT2p(f_ops):

  ops = [S.Zero, S.Zero]
  for i in range(1,3):
    for j in range(1,4):
      for k in range(1,3):
        p1 = Pi[mod3(j+1)] + Pi[mod3(j+2)]
        p2 = Pi[mod3(j+1)] - Pi[mod3(j+2)]

        ops[i-1] += Eijk(i,k)*L[k][j]*(project_momentum(f_ops[j], p1, -p1) \
                                     - project_momentum(f_ops[j], p2, -p2))

  return ops

def _P0_LT2p_Ep(f_ops):

  ops = [S.Zero, S.Zero]
  for i in range(1,3):
    for j in range(1,4):
      for k in range(1,3):
        p1 = Pi[mod3(j+1)] + Pi[mod3(j+2)]
        p2 = Pi[mod3(j+1)] - Pi[mod3(j+2)]

        ops[i-1] += Eijk(i,k)*L[k][j]*(project_momentum(f_ops[0], p1, -p1) \
                                     - project_momentum(f_ops[0], p2, -p2))

  return ops

#############    T1p    ############################################################################
def _P0_T1p_LA1p(f_ops, n=0):

  ops = [S.Zero, S.Zero, S.Zero]
  for i in range(1,4):
    for p in MOMENTA[n]:
      ops[i-1] += project_momentum(f_ops[i], p, -p)

  return ops

def _P0_LA1p_T1p(f_ops, n=0):

  op = S.Zero
  for p in MOMENTA[n]:
    op += project_momentum(f_ops[0], p, -p)

  return [op]

def _P0_T1p_LEp_1(f_ops):

  ops = [S.Zero, S.Zero, S.Zero]
  for i in range(1,4):
    ops[i-1] = project_momentum(f_ops[i], Pi[i], -Pi[i])

    for j in range(1,4):
      ops[i-1] -= Rational(1,3)*project_momentum(f_ops[i], Pi[j], -Pi[j])

  return ops

def _P0_LEp_T1p_1(f_ops):

  ops = [S.Zero, S.Zero, S.Zero]
  for i in range(1,4):
    ops[i-1] = project_momentum(f_ops[0], Pi[i], -Pi[i])

    for j in range(1,4):
      ops[i-1] -= Rational(1,3)*project_momentum(f_ops[0], Pi[j], -Pi[j])

  return ops

def _P0_T1p_LEp_2(f_ops):

  ops = [S.Zero, S.Zero, S.Zero]
  for i in range(1,4):
    p1 = Pi[mod3(i+1)] + Pi[mod3(i+2)]
    p2 = Pi[mod3(i+1)] - Pi[mod3(i+2)]
    ops[i-1] = Rational(1,2)*(project_momentum(f_ops[i], p1, -p1) \
                            + project_momentum(f_ops[i], p2, -p2))

    for j in range(1,4):
      for k in range(j+1,4):
        ops[i-1] -= Rational(1,6)*(project_momentum(f_ops[i], Pi[j]+Pi[k], -Pi[j]-Pi[k]) \
                                 + project_momentum(f_ops[i], Pi[j]-Pi[k], -Pi[j]+Pi[k]))

  return ops

def _P0_LEp_T1p_2(f_ops):

  ops = [S.Zero, S.Zero, S.Zero]
  for i in range(1,4):
    p1 = Pi[mod3(i+1)] + Pi[mod3(i+2)]
    p2 = Pi[mod3(i+1)] - Pi[mod3(i+2)]
    ops[i-1] = Rational(1,2)*(project_momentum(f_ops[0], p1, -p1) \
                            + project_momentum(f_ops[0], p2, -p2))

    for j in range(1,4):
      for k in range(j+1,4):
        ops[i-1] -= Rational(1,6)*(project_momentum(f_ops[0], Pi[j]+Pi[k], -Pi[j]-Pi[k]) \
                                 + project_momentum(f_ops[0], Pi[j]-Pi[k], -Pi[j]+Pi[k]))

  return ops

def _P0_T1p_LT2p(f_ops):

  ops = [S.Zero, S.Zero, S.Zero]
  for i in range(1,4):
    for j in range(1,3):
      i_j = mod3(i+j)
      ops[i-1] += Rational(1,4)*(project_momentum(f_ops[i_j], Pi[i]+Pi[i_j], -Pi[i]-Pi[i_j]) \
                               - project_momentum(f_ops[i_j], Pi[i]-Pi[i_j], -Pi[i]+Pi[i_j]))

  return ops

def _P0_LT2p_T1p(f_ops):

  ops = [S.Zero, S.Zero, S.Zero]
  for i in range(1,4):
    for j in range(1,3):
      i_j = mod3(i+j)
      ops[i-1] += Rational(1,4)*(project_momentum(f_ops[0], Pi[i]+Pi[i_j], -Pi[i]-Pi[i_j]) \
                               - project_momentum(f_ops[0], Pi[i]-Pi[i_j], -Pi[i]+Pi[i_j]))

  return ops

def _P0_T1p_1(f_ops):

  ops = [S.Zero, S.Zero, S.Zero]
  for i in range(1,4):
    for p in MOMENTA[2]:
      for j in range(1,4):
        ops[i-1] += (p[i]*p[j] - KroneckerDelta(i,j)*Rational(p.sq,3))*project_momentum(f_ops[j], p, -p)

  return ops

def _P0_T1p_3(f_ops):

  ops = [S.Zero, S.Zero, S.Zero]
  for i in range(1,4):
    for p in MOMENTA[2]:
      ops[i-1] += (p[i]**2 - Rational(p.sq,3))*project_momentum(f_ops[i], p, -p)
      for j in range(1,4):
        ops[i-1] -= Rational(2,5)*(p[i]*p[j] - KroneckerDelta(i,j)*Rational(p.sq,3))*project_momentum(f_ops[j], p, -p)

  return ops

#############    T2p    ############################################################################
def _P0_T2p_LEp_1(f_ops):

  ops = [S.Zero, S.Zero, S.Zero]
  for i in range(1,4):
    for j in range(1,3):
      i_j = mod3(i+j)
      ops[i-1] += (-1)**j*project_momentum(f_ops[i], Pi[i_j], -Pi[i_j])

  return ops

def _P0_LEp_T2p_1(f_ops):

  ops = [S.Zero, S.Zero, S.Zero]
  for i in range(1,4):
    for j in range(1,3):
      i_j = mod3(i+j)
      ops[i-1] += (-1)**j*project_momentum(f_ops[0], Pi[i_j], -Pi[i_j])

  return ops

def _P0_T2p_LT2p_0(f_ops):

  ops = [S.Zero, S.Zero, S.Zero]
  for i in range(1,4):
    p1 = Pi[mod3(i+1)] + Pi[mod3(i+2)]
    p2 = Pi[mod3(i+1)] - Pi[mod3(i+2)]
    ops[i-1] += project_momentum(f_ops[0], p1, -p1) - project_momentum(f_ops[0], p2, -p2)

  return ops

def _P0_T2p_LEp_2(f_ops):
  ops = [S.Zero, S.Zero, S.Zero]
  for i in range(1,4):
    for j in range(1,3):
      i_j = mod3(i+j)
      p1 = Pi[i] + Pi[i_j]
      p2 = Pi[i] - Pi[i_j]

      ops[i-1] += (-1)**j*(project_momentum(f_ops[i], p1, -p1) \
                         + project_momentum(f_ops[i], p2, -p2))

  return ops

def _P0_LEp_T2p_2(f_ops):
  ops = [S.Zero, S.Zero, S.Zero]
  for i in range(1,4):
    for j in range(1,3):
      i_j = mod3(i+j)
      p1 = Pi[i] + Pi[i_j]
      p2 = Pi[i] - Pi[i_j]

      ops[i-1] += (-1)**j*(project_momentum(f_ops[0], p1, -p1) \
                         + project_momentum(f_ops[0], p2, -p2))

  return ops

def _P0_T2p_LT2p_1(f_ops):

  ops = [S.Zero, S.Zero, S.Zero]
  for i in range(1,4):
    for j in range(1,4):
      for k in range(1,4):
        p1 = Pi[mod3(k+1)] + Pi[mod3(k+2)]
        p2 = Pi[mod3(k+1)] - Pi[mod3(k+2)]
        ops[i-1] += Eijk(i,j,k)*(project_momentum(f_ops[j], p1, -p1) \
                               - project_momentum(f_ops[j], p2, -p2))

  return ops

def _P0_LT2p_T2p_1(f_ops):

  ops = [S.Zero, S.Zero, S.Zero]
  for i in range(1,4):
    for j in range(1,4):
      for k in range(1,4):
        p1 = Pi[mod3(k+1)] + Pi[mod3(k+2)]
        p2 = Pi[mod3(k+1)] - Pi[mod3(k+2)]
        ops[i-1] += Eijk(i,j,k)*(project_momentum(f_ops[0], p1, -p1) \
                               - project_momentum(f_ops[0], p2, -p2))

  return ops

def _P0_T2p_2(f_ops):

  ops = [S.Zero, S.Zero, S.Zero]
  for i in range(1,4):
    for p in MOMENTA[2]:
      for j in range(1,3):
        i_j = mod3(i+j)
        ops[i-1] += (-1)**j*(p[i_j]**2*project_momentum(f_ops[i], p, -p) \
                         - p[i]*p[i_j]*project_momentum(f_ops[i_j], p, -p))

  return ops

def _P0_T2p_3(f_ops):

  ops = [S.Zero, S.Zero, S.Zero]
  for i in range(1,4):
    for p in MOMENTA[2]:
      for j in range(1,3):
        i_j = mod3(i+j)
        ops[i-1] += (-1)**j*(2*p[i]*p[i_j]*project_momentum(f_ops[i_j], p, -p) \
                               + p[i_j]**2*project_momentum(f_ops[i], p, -p))

  return ops

#############    A1m    ############################################################################

def _P0_A1m_LT1m(f_ops, n=1):
  op = S.Zero
  for p in MOMENTA[n]:
    for i in range(1,4):
      op += p[i]*project_momentum(f_ops[i], p, -p)

  return [op]

def _P0_LT1m_A1m(f_ops, n=1):
  op = S.Zero
  for p in MOMENTA[n]:
    for i in range(1,4):
      op += p[i]*project_momentum(f_ops[0], p, -p)

  return [op]

#############    A2m    ############################################################################

def _P0_A2m_LT2m(f_ops):
  op = S.Zero
  for i in range(1,4):
    for j in range(1,3):
      i_j = mod3(i+j)
      p1 = Pi[i] + Pi[i_j]
      p2 = Pi[i] - Pi[i_j]

      op += (-1)**j*(project_momentum(f_ops[i], p1, -p1) + project_momentum(f_ops[i], p2, -p2))

  return [op]

def _P0_LT2m_A2m(f_ops):
  op = S.Zero
  for i in range(1,4):
    for j in range(1,3):
      i_j = mod3(i+j)
      p1 = Pi[i] + Pi[i_j]
      p2 = Pi[i] - Pi[i_j]

      op += (-1)**j*(project_momentum(f_ops[0], p1, -p1) + project_momentum(f_ops[0], p2, -p2))

  return [op]

#############    Em     ############################################################################
def _P0_Em_LT1m(f_ops, n=1):

  ops = [S.Zero, S.Zero]
  for i in range(1,3):
    for p in MOMENTA[n]:
      for j in range(1,4):
        ops[i-1] += L[i][j]*p[j]*project_momentum(f_ops[j], p, -p)

  return ops

def _P0_LT1m_Em(f_ops, n=1):

  ops = [S.Zero, S.Zero]
  for i in range(1,3):
    for p in MOMENTA[n]:
      for j in range(1,4):
        ops[i-1] += L[i][j]*p[j]*project_momentum(f_ops[0], p, -p)

  return ops

def _P0_Em_LT2m(f_ops):

  ops = [S.Zero, S.Zero]
  for j in range(1,3):
    for l in range(1,3):
      for k in range(1,4):
        p1 = Pi[k] + Pi[mod3(k+j)]
        p2 = Pi[k] - Pi[mod3(k+j)]

        for i in range(1,3):
          ops[i-1] += Eijk(i,l)*L[l][k]*(-1)**j*(project_momentum(f_ops[k], p1, -p1) \
                                               + project_momentum(f_ops[k], p2, -p2))

  return ops

def _P0_LT2m_Em(f_ops):

  ops = [S.Zero, S.Zero]
  for j in range(1,3):
    for l in range(1,3):
      for k in range(1,4):
        p1 = Pi[k] + Pi[mod3(k+j)]
        p2 = Pi[k] - Pi[mod3(k+j)]

        for i in range(1,3):
          ops[i-1] += Eijk(i,l)*L[l][k]*(-1)**j*(project_momentum(f_ops[0], p1, -p1) \
                                               + project_momentum(f_ops[0], p2, -p2))

  return ops

#############    T1m    ############################################################################

def _P0_T1m_LT1m_0(f_ops, n=1):
  ops = [S.Zero, S.Zero, S.Zero]
  for i in range(1,4):
    for p in MOMENTA[n]:
      ops[i-1] += p[i]*project_momentum(f_ops[0], p, -p)

  return ops

def _P0_T1m_LT1m_1(f_ops, n=1):
  ops = [S.Zero, S.Zero, S.Zero]
  for i in range(1,4):
    for p in MOMENTA[n]:
      for j in range(1,4):
        for k in range(1,4):
          ops[i-1] += Eijk(i,j,k)*p[j]*project_momentum(f_ops[k], p, -p)

  return ops

def _P0_LT1m_T1m_1(f_ops, n=1):
  ops = [S.Zero, S.Zero, S.Zero]
  for i in range(1,4):
    for p in MOMENTA[n]:
      for j in range(1,4):
        for k in range(1,4):
          ops[i-1] += Eijk(i,j,k)*p[j]*project_momentum(f_ops[0], p, -p)

  return ops

def _P0_T1m_LT2m(f_ops):
  ops = [S.Zero, S.Zero, S.Zero]

  for i in range(1,4):
    for j in range(1,4):
      if i == j:
        continue
      for k in range(1,4):
        if i == k or j == k:
          continue

        for l in range(1,3):
          p1 = Pi[k] + Pi[mod3(k+l)]
          p2 = Pi[k] - Pi[mod3(k+l)]
          ops[i-1] += (-1)**l*(project_momentum(f_ops[j], p1, -p1) \
                             + project_momentum(f_ops[j], p2, -p2))

  return ops

def _P0_LT2m_T1m(f_ops):
  ops = [S.Zero, S.Zero, S.Zero]

  for i in range(1,4):
    for j in range(1,4):
      if i == j:
        continue
      for k in range(1,4):
        if i == k or j == k:
          continue

        for l in range(1,3):
          p1 = Pi[k] + Pi[mod3(k+l)]
          p2 = Pi[k] - Pi[mod3(k+l)]
          ops[i-1] += (-1)**l*(project_momentum(f_ops[0], p1, -p1) \
                             + project_momentum(f_ops[0], p2, -p2))

  return ops

#############    T2m    ############################################################################

def _P0_T2m_LT1m(f_ops, n=1):
  ops = [S.Zero, S.Zero, S.Zero]
  for i in range(1,4):
    for p in MOMENTA[n]:
      for j in range(1,4):
        if j == i:
          continue
        for k in range(1,4):
          if i == k or j == k:
            continue

          ops[i-1] += p[j]*project_momentum(f_ops[k], p, -p)

  return ops

def _P0_LT1m_T2m(f_ops, n=1):
  ops = [S.Zero, S.Zero, S.Zero]
  for i in range(1,4):
    for p in MOMENTA[n]:
      for j in range(1,4):
        if j == i:
          continue
        for k in range(1,4):
          if i == k or j == k:
            continue

          ops[i-1] += p[j]*project_momentum(f_ops[0], p, -p)

  return ops

def _P0_T2m_LT2m_1(f_ops):
  ops = [S.Zero, S.Zero, S.Zero]
  for i in range(1,4):
    for j in range(1,4):
      for k in range(1,4):
        for l in range(1,3):
          p1 = Pi[k] + Pi[mod3(k+l)]
          p2 = Pi[k] - Pi[mod3(k+l)]

          ops[i-1] += Eijk(i,j,k)*(-1)**l*(project_momentum(f_ops[j], p1, -p1) \
                                         + project_momentum(f_ops[j], p2, -p2))

  return ops

def _P0_LT2m_T2m_1(f_ops):
  ops = [S.Zero, S.Zero, S.Zero]
  for i in range(1,4):
    for j in range(1,4):
      for k in range(1,4):
        for l in range(1,3):
          p1 = Pi[k] + Pi[mod3(k+l)]
          p2 = Pi[k] - Pi[mod3(k+l)]

          ops[i-1] += Eijk(i,j,k)*(-1)**l*(project_momentum(f_ops[0], p1, -p1) \
                                         + project_momentum(f_ops[0], p2, -p2))

  return ops

def _P0_T2m_LT2m_0(f_ops):
  ops = [S.Zero, S.Zero, S.Zero]
  for i in range(1,4):
    for j in range(1,3):
      p1 = Pi[i] + Pi[mod3(i+j)]
      p2 = Pi[i] - Pi[mod3(i+j)]

      ops[i-1] += (-1)**j*(project_momentum(f_ops[0], p1, -p1) \
                         + project_momentum(f_ops[0], p2, -p2))
  
  return ops

#############  P^2 = 1  ############################################################################
#############    A1     ############################################################################

def _P1_A1_LA1_10(f_ops, Ptot):
  if Ptot.sq != 1:
    raise ValueError("P^2 != 1")

  op = project_momentum(f_ops[0], Ptot, P0)

  return [op]

def _P1_A1_LA1_21(f_ops, Ptot):
  if Ptot.sq != 1:
    raise ValueError("P^2 != 1")

  op = S.Zero
  for j in range(1,4):
    if Ptot[j] != 0:
      continue

    op += Rational(1,4)*(project_momentum(f_ops[0], Ptot + Pi[j], -Pi[j]) \
                       + project_momentum(f_ops[0], Ptot - Pi[j],  Pi[j]))

  return [op]

def _P1_A1_LE2_21(f_ops, Ptot):
  if Ptot.sq != 1:
    raise ValueError("P^2 != 1")

  op = S.Zero
  for j in range(1,4):
    if Ptot[j] != 0:
      continue

    op += Rational(1,4)*(project_momentum(f_ops[j], Ptot + Pi[j]*Ptot, -Pi[j]*Ptot) \
                       - project_momentum(f_ops[j], Ptot - Pi[j]*Ptot,  Pi[j]*Ptot))

  return [op]

def _P1_LE2_A1_21(f_ops, Ptot):
  if Ptot.sq != 1:
    raise ValueError("P^2 != 1")

  op = S.Zero
  for j in range(1,4):
    if Ptot[j] != 0:
      continue

    op += Rational(1,4)*(project_momentum(f_ops[0], Ptot+Pi[j]*Ptot, -Pi[j]*Ptot) \
                       - project_momentum(f_ops[0], Ptot-Pi[j]*Ptot,  Pi[j]*Ptot))

  return [op]

#############    A2     ############################################################################

def _P1_A2_LA1_10(f_ops, Ptot):
  if Ptot.sq != 1:
    raise ValueError("P^2 != 1")

  op = S.Zero
  for i in range(1,4):
    op += Ptot[i]*project_momentum(f_ops[i], Ptot, P0)

  return [op]

def _P1_LA1_A2_10(f_ops, Ptot):
  if Ptot.sq != 1:
    raise ValueError("P^2 != 1")

  op = S.Zero
  for i in range(1,4):
    op += Ptot[i]*project_momentum(f_ops[0], Ptot, P0)

  return [op]

def _P1_A2_LA1_21(f_ops, Ptot):
  if Ptot.sq != 1:
    raise ValueError("P^2 != 1")

  op = S.Zero
  for i in range(1,4):
    for j in range(1,4):
      if Ptot[j] != 0:
        continue

      op += Rational(Ptot[i],4)*(project_momentum(f_ops[i], Ptot + Pi[j], -Pi[j]) \
                               + project_momentum(f_ops[i], Ptot - Pi[j],  Pi[j]))
  return [op]

def _P1_LA1_A2_21(f_ops, Ptot):
  if Ptot.sq != 1:
    raise ValueError("P^2 != 1")

  op = S.Zero
  for i in range(1,4):
    for j in range(1,4):
      if Ptot[j] != 0:
        continue

      op += Rational(Ptot[i],4)*(project_momentum(f_ops[0], Ptot + Pi[j], -Pi[j]) \
                               + project_momentum(f_ops[0], Ptot - Pi[j],  Pi[j]))
  return [op]

def _P1_A2_LE2_21(f_ops, Ptot):
  if Ptot.sq != 1:
    raise ValueError("P^2 != 1")

  op = S.Zero
  for j in range(1,4):
    if Ptot[j] != 0:
      continue

    op += Rational(1,4)*(project_momentum(f_ops[j], Ptot + Pi[j], -Pi[j]) \
                       - project_momentum(f_ops[j], Ptot - Pi[j],  Pi[j]))

  return [op]

def _P1_LE2_A2_21(f_ops, Ptot):
  if Ptot.sq != 1:
    raise ValueError("P^2 != 1")

  op = S.Zero
  for j in range(1,4):
    if Ptot[j] != 0:
      continue

    op += Rational(1,4)*(project_momentum(f_ops[0], Ptot + Pi[j], -Pi[j]) \
                       - project_momentum(f_ops[0], Ptot - Pi[j],  Pi[j]))

  return [op]

#############    B1     ############################################################################

def _P1_B1_LB1_21(f_ops, Ptot):
  if Ptot.sq != 1:
    raise ValueError("P^2 != 1")

  k = Pks(Ptot)

  op = S.Zero
  for i in range(1,3):
    op += Rational(1,4)*(-1)**i*(project_momentum(f_ops[0], Ptot - k[i],  k[i]) \
                               + project_momentum(f_ops[0], Ptot + k[i], -k[i]))

  return [op]

def _P1_B1_LE2_21(f_ops, Ptot):
  if Ptot.sq != 1:
    raise ValueError("P^2 != 1")

  k = Pks(Ptot)
  ks = Pksi(Ptot)

  op = S.Zero
  for i in range(1,3):
    op += Rational(1,4)*(-1)**i*(project_momentum(f_ops[ks[i]], Ptot + k[i]*Ptot, -k[i]*Ptot) \
                               - project_momentum(f_ops[ks[i]], Ptot - k[i]*Ptot,  k[i]*Ptot))

  return [op]

def _P1_LE2_B1_21(f_ops, Ptot):
  if Ptot.sq != 1:
    raise ValueError("P^2 != 1")

  k = Pks(Ptot)

  op = S.Zero
  for i in range(1,3):
    op += Rational(1,4)*(-1)**i*(project_momentum(f_ops[0], Ptot + k[i]*Ptot, -k[i]*Ptot) \
                               - project_momentum(f_ops[0], Ptot - k[i]*Ptot,  k[i]*Ptot))

  return [op]

#############    B2     ############################################################################

def _P1_B2_LB1_21(f_ops, Ptot):
  if Ptot.sq != 1:
    raise ValueError("P^2 != 1")

  k = Pks(Ptot)
  
  op = S.Zero
  for i in range(1,3):
    for j in range(1,4):
      op += Rational(1,4)*(-1)**i*Ptot[j]*(project_momentum(f_ops[j], Ptot - k[i],  k[i]) \
                                         + project_momentum(f_ops[j], Ptot + k[i], -k[i]))

  return [op]

def _P1_LB1_B2_21(f_ops, Ptot):
  if Ptot.sq != 1:
    raise ValueError("P^2 != 1")

  k = Pks(Ptot)
  
  op = S.Zero
  for i in range(1,3):
    for j in range(1,4):
      op += Rational(1,4)*(-1)**i*Ptot[j]*(project_momentum(f_ops[0], Ptot - k[i],  k[i]) \
                                         + project_momentum(f_ops[0], Ptot + k[i], -k[i]))

  return [op]

def _P1_B2_LE2_21(f_ops, Ptot):
  if Ptot.sq != 1:
    raise ValueError("P^2 != 1")

  k = Pks(Ptot)
  ks = Pksi(Ptot)

  op = S.Zero
  for i in range(1,3):
    op += Rational(1,4)*(-1)**i*(project_momentum(f_ops[ks[i]], Ptot + k[i], -k[i]) \
                               - project_momentum(f_ops[ks[i]], Ptot - k[i],  k[i]))

  return [op]

def _P1_LE2_B2_21(f_ops, Ptot):
  if Ptot.sq != 1:
    raise ValueError("P^2 != 1")

  k = Pks(Ptot)

  op = S.Zero
  for i in range(1,3):
    op += Rational(1,4)*(-1)**i*(project_momentum(f_ops[0], Ptot + k[i], -k[i]) \
                               - project_momentum(f_ops[0], Ptot - k[i],  k[i]))

  return [op]

#############    E2     ############################################################################

def _P1_E2_LA1_10(f_ops, Ptot):
  if Ptot.sq != 1:
    raise ValueError("P^2 != 1")

  ks = Pksi(Ptot)

  op = [project_momentum(f_ops[ks[1]], Ptot, P0), project_momentum(f_ops[ks[2]], Ptot, P0)]

  return op

def _P1_LA1_E2_10(f_ops, Ptot):
  if Ptot.sq != 1:
    raise ValueError("P^2 != 1")

  op = project_momentum(f_ops[0], Ptot, P0)

  return [op]


def _P1_E2_LE2_21_0(f_ops, Ptot):
  if Ptot.sq != 1:
    raise ValueError("P^2 != 1")

  k = Pks(Ptot)

  op_1 = Rational(1,2)*(project_momentum(f_ops[0], Ptot + k[1]*Ptot, -k[1]*Ptot) \
                      - project_momentum(f_ops[0], Ptot - k[1]*Ptot, k[1]*Ptot))
  op_2 = Rational(1,2)*(project_momentum(f_ops[0], Ptot + k[2]*Ptot, -k[2]*Ptot) \
                      - project_momentum(f_ops[0], Ptot - k[2]*Ptot, k[2]*Ptot))

  return [op_1, op_2]


def _P1_E2_LA1_21(f_ops, Ptot):
  if Ptot.sq != 1:
    raise ValueError("P^2 != 1")

  ks = Pksi(Ptot)

  op1 = S.Zero
  op2 = S.Zero
  for j in range(1,4):
    if Ptot[j] != 0:
      continue

    op1 += Rational(1,4)*(project_momentum(f_ops[ks[1]], Ptot + Pi[j], -Pi[j]) \
                        + project_momentum(f_ops[ks[1]], Ptot - Pi[j],  Pi[j]))

    op2 += Rational(1,4)*(project_momentum(f_ops[ks[2]], Ptot + Pi[j], -Pi[j]) \
                        + project_momentum(f_ops[ks[2]], Ptot - Pi[j],  Pi[j]))

  return [op1, op2]

def _P1_LA1_E2_21(f_ops, Ptot):
  if Ptot.sq != 1:
    raise ValueError("P^2 != 1")

  ks = Pksi(Ptot)

  op1 = S.Zero
  op2 = S.Zero
  for j in range(1,4):
    if Ptot[j] != 0:
      continue

    op1 += Rational(1,4)*(project_momentum(f_ops[0], Ptot + Pi[j], -Pi[j]) \
                        + project_momentum(f_ops[0], Ptot - Pi[j],  Pi[j]))

    op2 += Rational(1,4)*(project_momentum(f_ops[0], Ptot + Pi[j], -Pi[j]) \
                        + project_momentum(f_ops[0], Ptot - Pi[j],  Pi[j]))

  return [op1, op2]

def _P1_E2_LB1_21(f_ops, Ptot):
  if Ptot.sq != 1:
    raise ValueError("P^2 != 1")

  k = Pks(Ptot)
  ks = Pksi(Ptot)

  op = [S.Zero, S.Zero]
  for i in range(1,3):
    for j in range(1,3):
      op[i-1] += Rational(1,4)*(-1)**(i+j)*(project_momentum(f_ops[ks[i]], Ptot - k[j],  k[j]) \
                                          + project_momentum(f_ops[ks[i]], Ptot + k[j], -k[j]))

  return op

def _P1_LB1_E2_21(f_ops, Ptot):
  if Ptot.sq != 1:
    raise ValueError("P^2 != 1")

  k = Pks(Ptot)

  op = [S.Zero, S.Zero]
  for i in range(1,3):
    for j in range(1,3):
      op[i-1] += Rational(1,4)*(-1)**(i+j)*(project_momentum(f_ops[0], Ptot - k[j],  k[j]) \
                                          + project_momentum(f_ops[0], Ptot + k[j], -k[j]))

  return op

def _P1_E2_LE2_21_1(f_ops, Ptot):
  if Ptot.sq != 1:
    raise ValueError("P^2 != 1")

  k = Pks(Ptot)

  op1 = S.Zero
  op2 = S.Zero
  for j in range(1,4):
    op1 += Rational(1,2)*Ptot[j]*(project_momentum(f_ops[j], Ptot + k[1], -k[1]) \
                                - project_momentum(f_ops[j], Ptot - k[1],  k[1]))

    op2 += Rational(1,2)*Ptot[j]*(project_momentum(f_ops[j], Ptot + k[2], -k[2]) \
                                - project_momentum(f_ops[j], Ptot - k[2],  k[2]))

  return [op1, op2]

def _P1_LE2_E2_21_1(f_ops, Ptot):
  if Ptot.sq != 1:
    raise ValueError("P^2 != 1")

  k = Pks(Ptot)

  op1 = S.Zero
  op2 = S.Zero
  for j in range(1,4):
    op1 += Rational(1,2)*Ptot[j]*(project_momentum(f_ops[0], Ptot + k[1], -k[1]) \
                                - project_momentum(f_ops[0], Ptot - k[1],  k[1]))

    op2 += Rational(1,2)*Ptot[j]*(project_momentum(f_ops[0], Ptot + k[2], -k[2]) \
                                - project_momentum(f_ops[0], Ptot - k[2],  k[2]))

  return [op1, op2]

#############  P^2 = 2  ############################################################################
#############    A1     ############################################################################

def _P2_A1_LA1_20(f_ops, Ptot):
  if Ptot.sq != 2:
    raise ValueError("P^2 != 2")

  op = project_momentum(f_ops[0], Ptot, P0)

  return [op]


def _P2_A1_LA1_11(f_ops, Ptot):
  if Ptot.sq != 2:
    raise ValueError("P^2 != 2")

  ds = Pc(Ptot)

  op = project_momentum(f_ops[0], ds[1], ds[2])

  return [op]


def _P2_A1_LB1_11(f_ops, Ptot):
  if Ptot.sq != 2:
    raise ValueError("P^2 != 2")

  ds = Pc(Ptot)

  p = ds[1]*ds[2]
  op = S.Zero

  for i in range(1,4):
    op += p[i]*project_momentum(f_ops[i], ds[1], ds[2])

  return [op]

def _P2_LB1_A1_11(f_ops, Ptot):
  if Ptot.sq != 2:
    raise ValueError("P^2 != 2")

  ds = Pc(Ptot)

  p = ds[1]*ds[2]
  op = S.Zero

  for i in range(1,4):
    op += p[i]*project_momentum(f_ops[0], ds[1], ds[2])

  return [op]

#############    A2     ############################################################################

def _P2_A2_LA1_20(f_ops, Ptot):
  if Ptot.sq != 2:
    raise ValueError("P^2 != 2")

  op = S.Zero

  for i in range(1,4):
    op += Ptot[i]*project_momentum(f_ops[i], Ptot, P0)

  return [op]

def _P2_LA1_A2_20(f_ops, Ptot):
  if Ptot.sq != 2:
    raise ValueError("P^2 != 2")

  op = S.Zero

  for i in range(1,4):
    op += Ptot[i]*project_momentum(f_ops[0], Ptot, P0)

  return [op]

def _P2_A2_LB1_11(f_ops, Ptot):
  if Ptot.sq != 2:
    raise ValueError("P^2 != 2")

  ds = Pc(Ptot)

  p = ds[1] - ds[2]

  op = S.Zero
  for i in range(1,4):
    op += p[i]*project_momentum(f_ops[i], ds[1], ds[2])

  return [op]

def _P2_LB1_A2_11(f_ops, Ptot):
  if Ptot.sq != 2:
    raise ValueError("P^2 != 2")

  ds = Pc(Ptot)

  p = ds[1] - ds[2]

  op = S.Zero
  for i in range(1,4):
    op += p[i]*project_momentum(f_ops[0], ds[1], ds[2])

  return [op]

def _P2_A2_LA1_11(f_ops, Ptot):
  if Ptot.sq != 2:
    raise ValueError("P^2 != 2")

  ds = Pc(Ptot)

  op = S.Zero
  for i in range(1,4):
    op += Ptot[i]*project_momentum(f_ops[i], ds[1], ds[2])

  return [op]

def _P2_LA1_A2_11(f_ops, Ptot):
  if Ptot.sq != 2:
    raise ValueError("P^2 != 2")

  ds = Pc(Ptot)

  op = S.Zero
  for i in range(1,4):
    op += Ptot[i]*project_momentum(f_ops[0], ds[1], ds[2])

  return [op]

#############    B1     ############################################################################

def _P2_B1_LA1_20(f_ops, Ptot):
  if Ptot.sq != 2:
    raise ValueError("P^2 != 2")

  ds = Pc(Ptot)

  p = ds[1]*ds[2]
  op = S.Zero
  for i in range(1,4):
    op += p[i]*project_momentum(f_ops[i], Ptot, P0)

  return [op]

def _P2_LA1_B1_20(f_ops, Ptot):
  if Ptot.sq != 2:
    raise ValueError("P^2 != 2")

  ds = Pc(Ptot)

  p = ds[1]*ds[2]
  op = S.Zero
  for i in range(1,4):
    op += p[i]*project_momentum(f_ops[0], Ptot, P0)

  return [op]

def _P2_B1_LA1_11(f_ops, Ptot):
  if Ptot.sq != 2:
    raise ValueError("P^2 != 2")

  ds = Pc(Ptot)

  p = ds[1]*ds[2]
  op = S.Zero
  for i in range(1,4):
    op += p[i]*project_momentum(f_ops[i], ds[1], ds[2])

  return [op]

def _P2_LA1_B1_11(f_ops, Ptot):
  if Ptot.sq != 2:
    raise ValueError("P^2 != 2")

  ds = Pc(Ptot)

  p = ds[1]*ds[2]
  op = S.Zero
  for i in range(1,4):
    op += p[i]*project_momentum(f_ops[0], ds[1], ds[2])

  return [op]

def _P2_B1_LB1_11(f_ops, Ptot):
  if Ptot.sq != 2:
    raise ValueError("P^2 != 2")

  ds = Pc(Ptot)

  op = project_momentum(f_ops[0], ds[1], ds[2])

  return [op]

#############    B2     ############################################################################

def _P2_B2_LA1_20(f_ops, Ptot):
  if Ptot.sq != 2:
    raise ValueError("P^2 != 2")

  ds = Pc(Ptot)

  p = ds[1] - ds[2]
  op = S.Zero
  for i in range(1,4):
    op += p[i]*project_momentum(f_ops[i], Ptot, P0)

  return [op]

def _P2_LA1_B2_20(f_ops, Ptot):
  if Ptot.sq != 2:
    raise ValueError("P^2 != 2")

  ds = Pc(Ptot)

  p = ds[1] - ds[2]
  op = S.Zero
  for i in range(1,4):
    op += p[i]*project_momentum(f_ops[0], Ptot, P0)

  return [op]

def _P2_B2_LB1_11(f_ops, Ptot):
  if Ptot.sq != 2:
    raise ValueError("P^2 != 2")

  ds = Pc(Ptot)

  op = S.Zero
  for i in range(1,4):
    op += Ptot[i]*project_momentum(f_ops[i], ds[1], ds[2])

  return [op]

def _P2_LB1_B2_11(f_ops, Ptot):
  if Ptot.sq != 2:
    raise ValueError("P^2 != 2")

  ds = Pc(Ptot)

  op = S.Zero
  for i in range(1,4):
    op += Ptot[i]*project_momentum(f_ops[0], ds[1], ds[2])

  return [op]

def _P2_B2_LA1_11(f_ops, Ptot):
  if Ptot.sq != 2:
    raise ValueError("P^2 != 2")

  ds = Pc(Ptot)

  p = ds[1] - ds[2]
  op = S.Zero
  for i in range(1,4):
    op += p[i]*project_momentum(f_ops[i], ds[1], ds[2])

  return [op]

def _P2_LA1_B2_11(f_ops, Ptot):
  if Ptot.sq != 2:
    raise ValueError("P^2 != 2")

  ds = Pc(Ptot)

  p = ds[1] - ds[2]
  op = S.Zero
  for i in range(1,4):
    op += p[i]*project_momentum(f_ops[0], ds[1], ds[2])

  return [op]


#############  P^2 = 3  ############################################################################
#############    A1     ############################################################################

def _P3_A1_LA1_30(f_ops, Ptot):
  if Ptot.sq != 3:
    raise ValueError("P^2 != 3")

  op = project_momentum(f_ops[0], Ptot, P0)

  return [op]


def _P3_A1_LA1_21(f_ops, Ptot):
  if Ptot.sq != 3:
    raise ValueError("P^2 != 3")

  ds = Pc(Ptot)

  op = S.Zero
  for i in range(1,4):
    op += Rational(1,3)*project_momentum(f_ops[0], Ptot - ds[i], ds[i])

  return [op]


def _P3_A1_LE2_21(f_ops, Ptot):
  if Ptot.sq != 3:
    raise ValueError("P^2 != 3")

  ds = Pc(Ptot)

  p1p2 = ds[1]*ds[2]
  coeff = p1p2.x*ds[3].x + p1p2.y*ds[3].y + p1p2.z*ds[3].z

  op = S.Zero
  for i in range(1,4):
    for j in range(1,4):
      for k in range(1,4):
        for l in range(1,4):
          op += coeff*Rational(1,6)*Eijk(i,j,k)*ds[j][l]*project_momentum(f_ops[l], Ptot - ds[k], ds[k])

  return [op]

def _P3_LE2_A1_21(f_ops, Ptot):
  if Ptot.sq != 3:
    raise ValueError("P^2 != 3")

  ds = Pc(Ptot)

  p1p2 = ds[1]*ds[2]
  coeff = p1p2.x*ds[3].x + p1p2.y*ds[3].y + p1p2.z*ds[3].z

  op = S.Zero
  for i in range(1,4):
    for j in range(1,4):
      for k in range(1,4):
        for l in range(1,4):
          op += coeff*Rational(1,6)*Eijk(i,j,k)*ds[j][l]*project_momentum(f_ops[0], Ptot - ds[k], ds[k])

  return [op]

#############    A2     ############################################################################

def _P3_A2_LA1_30(f_ops, Ptot):
  if Ptot.sq != 3:
    raise ValueError("P^2 != 3")

  op = S.Zero
  for i in range(1,4):
    op += Rational(Ptot[i],3)*project_momentum(f_ops[i], Ptot, P0)

  return [op]

def _P3_LA1_A2_30(f_ops, Ptot):
  if Ptot.sq != 3:
    raise ValueError("P^2 != 3")

  op = S.Zero
  for i in range(1,4):
    op += Rational(Ptot[i],3)*project_momentum(f_ops[0], Ptot, P0)

  return [op]

def _P3_A2_LA1_21(f_ops, Ptot):
  if Ptot.sq != 3:
    raise ValueError("P^2 != 3")

  ds = Pc(Ptot)

  op = S.Zero
  for i in range(1,4):
    for j in range(1,4):
      op += Rational(Ptot[i],3)*project_momentum(f_ops[i], Ptot - ds[j], ds[j])

  return [op]

def _P3_LA1_A2_21(f_ops, Ptot):
  if Ptot.sq != 3:
    raise ValueError("P^2 != 3")

  ds = Pc(Ptot)

  op = S.Zero
  for i in range(1,4):
    for j in range(1,4):
      op += Rational(Ptot[i],3)*project_momentum(f_ops[0], Ptot - ds[j], ds[j])

  return [op]

def _P3_A2_LE2_21(f_ops, Ptot):
  if Ptot.sq != 3:
    raise ValueError("P^2 != 3")

  ds = Pc(Ptot)

  op = S.Zero
  for i in range(1,4):
    for k in range(1,4):
      d = 3*ds[i]
      for j in range(1,4):
        d -= ds[j]

      op += Rational(1,9)*d[k]*project_momentum(f_ops[k], Ptot - ds[i], ds[i])

  return [op]

def _P3_LE2_A2_21(f_ops, Ptot):
  if Ptot.sq != 3:
    raise ValueError("P^2 != 3")

  ds = Pc(Ptot)

  op = S.Zero
  for i in range(1,4):
    for k in range(1,4):
      d = 3*ds[i]
      for j in range(1,4):
        d -= ds[j]

      op += Rational(1,9)*d[k]*project_momentum(f_ops[0], Ptot - ds[i], ds[i])

  return [op]

#############    E2     ############################################################################

def _P3_E2_LA1_30(f_ops, Ptot):
  if Ptot.sq != 3:
    raise ValueError("P^2 != 3")

  ops = [S.Zero, S.Zero]
  for i in range(1,3):
    for j in range(1,4):
      ops[i-1] += L[i][j]*Ptot[j]*project_momentum(f_ops[j], Ptot, P0)

  return ops

def _P3_LA1_E2_30(f_ops, Ptot):
  if Ptot.sq != 3:
    raise ValueError("P^2 != 3")

  ops = [S.Zero, S.Zero]
  for i in range(1,3):
    for j in range(1,4):
      ops[i-1] += L[i][j]*Ptot[j]*project_momentum(f_ops[0], Ptot, P0)

  return ops


def _P3_E2_LA1_21(f_ops, Ptot):
  if Ptot.sq != 3:
    raise ValueError("P^2 != 3")

  ds = Pc(Ptot)

  ops = [S.Zero, S.Zero]
  for i in range(1,3):
    for j in range(1,4):
      for k in range(1,4):
        ops[i-1] += Rational(1,3)*L[i][j]*Ptot[j]*project_momentum(f_ops[j], Ptot - ds[k], ds[k])

  return ops

def _P3_LA1_E2_21(f_ops, Ptot):
  if Ptot.sq != 3:
    raise ValueError("P^2 != 3")

  ds = Pc(Ptot)

  ops = [S.Zero, S.Zero]
  for i in range(1,3):
    for j in range(1,4):
      for k in range(1,4):
        ops[i-1] += Rational(1,3)*L[i][j]*Ptot[j]*project_momentum(f_ops[0], Ptot - ds[k], ds[k])

  return ops


def _P3_E2_LE2_21_0(f_ops, Ptot):
  if Ptot.sq != 3:
    raise ValueError("P^2 != 3")

  ds = Pc(Ptot)

  p1p2 = ds[1]*ds[2]
  coeff = p1p2.x*ds[3].x + p1p2.y*ds[3].y + p1p2.z*ds[3].z

  ops = [S.Zero, S.Zero]
  for i in range(1,3):
    for j in range(1,3):
      for k in range(1,4):
        ops[i-1] += coeff*Eijk(i,j)*L[j][k]*project_momentum(f_ops[0], Ptot - ds[k], ds[k])

  return ops


def _P3_E2_LE2_21_1L(f_ops, Ptot):
  if Ptot.sq != 3:
    raise ValueError("P^2 != 3")

  ds = Pc(Ptot)

  ops = [S.Zero, S.Zero]
  for i in range(1,3):
    for j in range(1,4):
      for k in range(1,4):
        ops[i-1] += L[i][j]*Ptot[k]*project_momentum(f_ops[k], Ptot - ds[j], ds[j])

  return ops

def _P3_LE2_E2_21_1L(f_ops, Ptot):
  if Ptot.sq != 3:
    raise ValueError("P^2 != 3")

  ds = Pc(Ptot)

  ops = [S.Zero, S.Zero]
  for i in range(1,3):
    for j in range(1,4):
      for k in range(1,4):
        ops[i-1] += L[i][j]*Ptot[k]*project_momentum(f_ops[0], Ptot - ds[j], ds[j])

  return ops

def _P3_E2_LE2_21_1T(f_ops, Ptot):
  if Ptot.sq != 3:
    raise ValueError("P^2 != 3")

  ds = Pc(Ptot)

  op1 = S.Zero
  op2 = S.Zero
  for j in range(1,4):
    for k in range(1,4):
      op1 += (L[1][j]*L[2][k] + L[2][j]*L[1][k])*Ptot[j]*project_momentum(f_ops[j], Ptot - ds[k], ds[k])
      op2 += (L[1][j]*L[1][k] - L[2][j]*L[2][k])*Ptot[j]*project_momentum(f_ops[j], Ptot - ds[k], ds[k])

  return [op1, op2]

def _P3_LE2_E2_21_1T(f_ops, Ptot):
  if Ptot.sq != 3:
    raise ValueError("P^2 != 3")

  ds = Pc(Ptot)

  op1 = S.Zero
  op2 = S.Zero
  for j in range(1,4):
    for k in range(1,4):
      op1 += (L[1][j]*L[2][k] + L[2][j]*L[1][k])*Ptot[j]*project_momentum(f_ops[0], Ptot - ds[k], ds[k])
      op2 += (L[1][j]*L[1][k] - L[2][j]*L[2][k])*Ptot[j]*project_momentum(f_ops[0], Ptot - ds[k], ds[k])

  return [op1, op2]



#############  P^2 = 4  ############################################################################
#############    A1     ############################################################################

def _P4_A1_LA1_11(f_ops, Ptot):
  if Ptot.sq != 4:
    raise ValueError("P^2 != 4")

  op = project_momentum(f_ops[0], Ptot/2, Ptot/2)

  return [op]

#############    A2     ############################################################################

def _P4_A2_LA1_11(f_ops, Ptot):
  if Ptot.sq != 4:
    raise ValueError("P^2 != 4")

  op = S.Zero
  for i in range(1,4):
    op += Ptot[i]/2*project_momentum(f_ops[i], Ptot/2, Ptot/2)

  return [op]

def _P4_LA1_A2_11(f_ops, Ptot):
  if Ptot.sq != 4:
    raise ValueError("P^2 != 4")

  op = S.Zero
  for i in range(1,4):
    op += Ptot[i]/2*project_momentum(f_ops[0], Ptot/2, Ptot/2)

  return [op]

#############    E2     ############################################################################

def _P4_E2_LA1_11(f_ops, Ptot):
  if Ptot.sq != 4:
    raise ValueError("P^2 != 4")

  if Ptot.x != 0:
    op1 = project_momentum(f_ops[2], Ptot/2, Ptot/2)
    op2 = project_momentum(f_ops[3], Ptot/2, Ptot/2)

  elif Ptot.y != 0:
    op1 = project_momentum(f_ops[1], Ptot/2, Ptot/2)
    op2 = project_momentum(f_ops[3], Ptot/2, Ptot/2)

  elif Ptot.z != 0:
    op1 = project_momentum(f_ops[1], Ptot/2, Ptot/2)
    op2 = project_momentum(f_ops[2], Ptot/2, Ptot/2)

  return [op1, op2]

def _P4_LA1_E2_11(f_ops, Ptot):
  if Ptot.sq != 4:
    raise ValueError("P^2 != 4")

  op = project_momentum(f_ops[0], Ptot/2, Ptot/2)

  return [op]


############### Operator Test Class Instances #####################################################

P0_A1p = OperatorTest("P0 A1p L=A1p", "A1g", 0, _P0_A1p)

P0_A2p_LT2p = OperatorTest("P0 A2p L=T2p", "A2g", 0, _P0_A2p_LT2p)
P0_LT2p_A2p = OperatorTest("P0 L=T2p A2p", "T2g", 0, _P0_LT2p_A2p)

P0_Ep_LEp_1 = OperatorTest("P0 Ep L=Ep (p_i^2 = 1)", "Eg", 0, _P0_Ep_LEp_1)
P0_Ep_LEp_2 = OperatorTest("P0 Ep L=Ep (p_i^2 = 2)", "Eg", 0, _P0_Ep_LEp_2)
P0_Ep_LT2p = OperatorTest("P0 Ep L=T2p (p_i^2 = 2)", "Eg", 0, _P0_Ep_LT2p)
P0_LT2p_Ep = OperatorTest("P0 L=T2p Ep (p_i^2 = 2)", "Eg", 0, _P0_LT2p_Ep)

P0_T1p_LA1p = OperatorTest("P0 T1p L=A1p", "T1g", 0, _P0_T1p_LA1p)
P0_LA1p_T1p = OperatorTest("P0 L=A1p T1p", "A1g", 0, _P0_LA1p_T1p)
P0_T1p_LEp_1 = OperatorTest("P0 T1p L=Ep (p_i^2 = 1)", "T1g", 0, _P0_T1p_LEp_1)
P0_LEp_T1p_1 = OperatorTest("P0 L=Ep T1p (p_i^2 = 1)", "Eg", 0, _P0_LEp_T1p_1)
P0_T1p_LEp_2 = OperatorTest("P0 T1p L=Ep (p_i^2 = 2)", "T1g", 0, _P0_T1p_LEp_2)
P0_LEp_T1p_2 = OperatorTest("P0 L=Ep T1p (p_i^2 = 2)", "Eg", 0, _P0_LEp_T1p_2)
P0_T1p_LT2p = OperatorTest("P0 T1p L=T2p (p_i^2 = 2)", "T1g", 0, _P0_T1p_LT2p)
P0_LT2p_T1p = OperatorTest("P0 L=T2p T1p (p_i^2 = 2)", "T2g", 0, _P0_LT2p_T1p)
P0_T1p_1 = OperatorTest("P0 T1p J=1 (p_i^2 = 2)", "T1g", 0, _P0_T1p_1)
P0_T1p_3 = OperatorTest("P0 T1p J=3 (p_i^2 = 2)", "T1g", 0, _P0_T1p_3)

P0_T2p_LEp_1 = OperatorTest("P0 T2p L=Ep (p_i^2 = 1)", "T2g", 0, _P0_T2p_LEp_1)
P0_LEp_T2p_1 = OperatorTest("P0 L=Ep T2p (p_i^2 = 1)", "Eg", 0, _P0_LEp_T2p_1)
P0_T2p_LT2p_0 = OperatorTest("P0 T2p L=T2p 0 (p_i^2 = 2)", "T2g", 0, _P0_T2p_LT2p_0)
P0_T2p_LEp_2 = OperatorTest("P0 T2p L=Ep (p_i^2 = 2)", "T2g", 0, _P0_T2p_LEp_2)
P0_LEp_T2p_2 = OperatorTest("P0 L=Ep T2p (p_i^2 = 2)", "Eg", 0, _P0_LEp_T2p_2)
P0_T2p_LT2p_1 = OperatorTest("P0 T2p L=T2p 1 (p_i^2 = 2)", "T2g", 0, _P0_T2p_LT2p_1)
P0_LT2p_T2p_1 = OperatorTest("P0 L=T2p T2p 1 (p_i^2 = 2)", "T2g", 0, _P0_LT2p_T2p_1)
P0_T2p_2 = OperatorTest("P0 T2p J=2 (p_i^2 = 2)", "T2g", 0, _P0_T2p_2)
P0_T2p_3 = OperatorTest("P0 T2p J=3 (p_i^2 = 2)", "T2g", 0, _P0_T2p_3)

P0_A1m_LT1m = OperatorTest("P0 A1m L=T1m", "A1u", 0, _P0_A1m_LT1m)
P0_LT1m_A1m = OperatorTest("P0 L=T1m A1m", "T1u", 0, _P0_LT1m_A1m)

P0_A2m_LT2m = OperatorTest("P0 A2m L=T2m (p_i^2 = 2)", "A2u", 0, _P0_A2m_LT2m)
P0_LT2m_A2m = OperatorTest("P0 L=T2m A2m (p_i^2 = 2)", "T2u", 0, _P0_LT2m_A2m)

P0_Em_LT1m = OperatorTest("P0 Em L=T1m", "Eu", 0, _P0_Em_LT1m)
P0_LT1m_Em = OperatorTest("P0 L=T1m Em", "T1u", 0, _P0_LT1m_Em)
P0_Em_LT2m = OperatorTest("P0 Em L=T2m (p_i^2 = 2)", "Eu", 0, _P0_Em_LT2m)
P0_LT2m_Em = OperatorTest("P0 L=T2m Em (p_i^2 = 2)", "T2u", 0, _P0_LT2m_Em)

P0_T1m_LT1m_0 = OperatorTest("P0 T1m L=T1m 0", "T1u", 0, _P0_T1m_LT1m_0)
P0_T1m_LT1m_1 = OperatorTest("P0 T1m L=T1m 1", "T1u", 0, _P0_T1m_LT1m_1)
P0_LT1m_T1m_1 = OperatorTest("P0 L=T1m T1m 1", "T1u", 0, _P0_LT1m_T1m_1)
P0_T1m_LT2m = OperatorTest("P0 T1m L=T2m (p_i^2 = 2)", "T1u", 0, _P0_T1m_LT2m)
P0_LT2m_T1m = OperatorTest("P0 L=T2m T1m (p_i^2 = 2)", "T2u", 0, _P0_LT2m_T1m)

P0_T2m_LT1m = OperatorTest("P0 T2m L=T1m", "T2u", 0, _P0_T2m_LT1m)
P0_LT1m_T2m = OperatorTest("P0 L=T1m T2m", "T1u", 0, _P0_LT1m_T2m)
P0_T2m_LT2m_1 = OperatorTest("P0 T2m L=T2m 1 (p_i^2 = 2)", "T2u", 0, _P0_T2m_LT2m_1)
P0_LT2m_T2m_1 = OperatorTest("P0 L=T2m T2m 1 (p_i^2 = 2)", "T2u", 0, _P0_LT2m_T2m_1)
P0_T2m_LT2m_0 = OperatorTest("P0 T2m L=T2m 0 (p_i^2 = 2)", "T2u", 0, _P0_T2m_LT2m_0)

P1_A1_LA1_10 = OperatorTest("P1 A1 L=A1 (1,0)", "A1", 1, _P1_A1_LA1_10)
P1_A1_LA1_21 = OperatorTest("P1 A1 L=A1 (2,1)", "A1", 1, _P1_A1_LA1_21)
P1_A1_LE2_21 = OperatorTest("P1 A1 L=E2 (2,1)", "A1", 1, _P1_A1_LE2_21)
P1_LE2_A1_21 = OperatorTest("P1 L=E2 A1 (2,1)", "E", 1, _P1_LE2_A1_21)

P1_A2_LA1_10 = OperatorTest("P1 A2 L=A1 (1,0)", "A2", 1, _P1_A2_LA1_10)
P1_LA1_A2_10 = OperatorTest("P1 L=A1 A2 (1,0)", "A1", 1, _P1_LA1_A2_10)
P1_A2_LA1_21 = OperatorTest("P1 A2 L=A1 (2,1)", "A2", 1, _P1_A2_LA1_21)
P1_LA1_A2_21 = OperatorTest("P1 L=A1 A2 (2,1)", "A1", 1, _P1_LA1_A2_21)
P1_A2_LE2_21 = OperatorTest("P1 A2 L=E2 (2,1)", "A2", 1, _P1_A2_LE2_21)
P1_LE2_A2_21 = OperatorTest("P1 L=E2 A2 (2,1)", "E", 1, _P1_LE2_A2_21)

P1_B1_LB1_21 = OperatorTest("P1 B1 L=B1 (2,1)", "B1", 1, _P1_B1_LB1_21)
P1_B1_LE2_21 = OperatorTest("P1 B1 L=E2 (2,1)", "B1", 1, _P1_B1_LE2_21)
P1_LE2_B1_21 = OperatorTest("P1 L=E2 B1 (2,1)", "E", 1, _P1_LE2_B1_21)

P1_B2_LB1_21 = OperatorTest("P1 B2 L=B1 (2,1)", "B2", 1, _P1_B2_LB1_21)
P1_LB1_B2_21 = OperatorTest("P1 L=B1 B2 (2,1)", "B1", 1, _P1_LB1_B2_21)
P1_B2_LE2_21 = OperatorTest("P1 B2 L=E2 (2,1)", "B2", 1, _P1_B2_LE2_21)
P1_LE2_B2_21 = OperatorTest("P1 L=E2 B2 (2,1)", "E", 1, _P1_LE2_B2_21)

P1_E2_LA1_10 = OperatorTest("P1 E2 L=A1 (1,0)", "E", 1, _P1_E2_LA1_10)
P1_LA1_E2_10 = OperatorTest("P1 L=A1 E2 (1,0)", "A1", 1, _P1_LA1_E2_10)
P1_E2_LE2_21_0 = OperatorTest("P1 E2 L=E2 0 (2,1)", "E", 1, _P1_E2_LE2_21_0)
P1_E2_LA1_21 = OperatorTest("P1 E2 L=A1 (2,1)", "E", 1, _P1_E2_LA1_21)
P1_LA1_E2_21 = OperatorTest("P1 L=A1 E2 (2,1)", "A1", 1, _P1_LA1_E2_21)
P1_E2_LB1_21 = OperatorTest("P1 E2 L=B1 (2,1)", "E", 1, _P1_E2_LB1_21)
P1_LB1_E2_21 = OperatorTest("P1 L=B1 E2 (2,1)", "B1", 1, _P1_LB1_E2_21)
P1_E2_LE2_21_1 = OperatorTest("P1 E2 L=E2 1 (2,1)", "E", 1, _P1_E2_LE2_21_1)
P1_LE2_E2_21_1 = OperatorTest("P1 L=E2 E2 1 (2,1)", "E", 1, _P1_LE2_E2_21_1)

P2_A1_LA1_20 = OperatorTest("P2 A1 L=A1 (2,0)", "A1", 2, _P2_A1_LA1_20)
P2_A1_LA1_11 = OperatorTest("P2 A1 L=A1 (1,1)", "A1", 2, _P2_A1_LA1_11)
P2_A1_LB1_11 = OperatorTest("P2 A1 L=B1 (1,1)", "A1", 2, _P2_A1_LB1_11)
P2_LB1_A1_11 = OperatorTest("P2 L=B1 A1 (1,1)", "B1", 2, _P2_LB1_A1_11)

P2_A2_LA1_20 = OperatorTest("P2 A2 L=A1 (2,0)", "A2", 2, _P2_A2_LA1_20)
P2_LA1_A2_20 = OperatorTest("P2 L=A1 A2 (2,0)", "A1", 2, _P2_LA1_A2_20)
P2_A2_LB1_11 = OperatorTest("P2 A2 L=B1 (1,1)", "A2", 2, _P2_A2_LB1_11)
P2_LB1_A2_11 = OperatorTest("P2 L=B1 A2 (1,1)", "B1", 2, _P2_LB1_A2_11)
P2_A2_LA1_11 = OperatorTest("P2 A2 L=A1 (1,1)", "A2", 2, _P2_A2_LA1_11)
P2_LA1_A2_11 = OperatorTest("P2 L=A1 A2 (1,1)", "A1", 2, _P2_LA1_A2_11)

P2_B1_LA1_20 = OperatorTest("P2 B1 L=A1 (2,0)", "B2", 2, _P2_B1_LA1_20)
P2_LA1_B1_20 = OperatorTest("P2 L=A1 B1 (2,0)", "A1", 2, _P2_LA1_B1_20)
P2_B1_LA1_11 = OperatorTest("P2 B1 L=A1 (1,1)", "B2", 2, _P2_B1_LA1_11)
P2_LA1_B1_11 = OperatorTest("P2 L=A1 B1 (1,1)", "A1", 2, _P2_LA1_B1_11)
P2_B1_LB1_11 = OperatorTest("P2 B1 L=B1 (1,1)", "B2", 2, _P2_B1_LB1_11)

P2_B2_LA1_20 = OperatorTest("P2 B2 L=A1 (2,0)", "B1", 2, _P2_B2_LA1_20)
P2_LA1_B2_20 = OperatorTest("P2 L=A1 B2 (2,0)", "A1", 2, _P2_LA1_B2_20)
P2_B2_LB1_11 = OperatorTest("P2 B2 L=B1 (1,1)", "B1", 2, _P2_B2_LB1_11)
P2_LB1_B2_11 = OperatorTest("P2 L=B1 B2 (1,1)", "B2", 2, _P2_LB1_B2_11)
P2_B2_LA1_11 = OperatorTest("P2 B2 L=A1 (1,1)", "B1", 2, _P2_B2_LA1_11)
P2_LA1_B2_11 = OperatorTest("P2 L=A1 B2 (1,1)", "A1", 2, _P2_LA1_B2_11)


P3_A1_LA1_30 = OperatorTest("P3 A1 L=A1 (3,0)", "A1", 3, _P3_A1_LA1_30)
P3_A1_LA1_21 = OperatorTest("P3 A1 L=A1 (2,1)", "A1", 3, _P3_A1_LA1_21)
P3_A1_LE2_21 = OperatorTest("P3 A1 L=E2 (2,1)", "A1", 3, _P3_A1_LE2_21)
P3_LE2_A1_21 = OperatorTest("P3 L=E2 A1 (2,1)", "E", 3, _P3_LE2_A1_21)

P3_A2_LA1_30 = OperatorTest("P3 A2 L=A1 (3,0)", "A2", 3, _P3_A2_LA1_30)
P3_LA1_A2_30 = OperatorTest("P3 L=A1 A2 (3,0)", "A1", 3, _P3_LA1_A2_30)
P3_A2_LA1_21 = OperatorTest("P3 A2 L=A1 (2,1)", "A2", 3, _P3_A2_LA1_21)
P3_LA1_A2_21 = OperatorTest("P3 L=A1 A2 (2,1)", "A1", 3, _P3_LA1_A2_21)
P3_A2_LE2_21 = OperatorTest("P3 A2 L=E2 (2,1)", "A2", 3, _P3_A2_LE2_21)
P3_LE2_A2_21 = OperatorTest("P3 L=E2 A2 (2,1)", "E", 3, _P3_LE2_A2_21)

P3_E2_LA1_30 = OperatorTest("P3 E2 L=A1 (3,0)", "E", 3, _P3_E2_LA1_30)
P3_LA1_E2_30 = OperatorTest("P3 L=A1 E2 (3,0)", "A1", 3, _P3_LA1_E2_30)
P3_E2_LA1_21 = OperatorTest("P3 E2 L=A1 (2,1)", "E", 3, _P3_E2_LA1_21)
P3_LA1_E2_21 = OperatorTest("P3 L=A1 E2 (2,1)", "A1", 3, _P3_LA1_E2_21)
P3_E2_LE2_21_0 = OperatorTest("P3 E2 L=E2 0 (2,1)", "E", 3, _P3_E2_LE2_21_0)
P3_E2_LE2_21_1L = OperatorTest("P3 E2 L=E2 1L (2,1)", "E", 3, _P3_E2_LE2_21_1L)
P3_LE2_E2_21_1L = OperatorTest("P3 L=E2 E2 1L (2,1)", "E", 3, _P3_LE2_E2_21_1L)
P3_E2_LE2_21_1T = OperatorTest("P3 E2 L=E2 1T (2,1)", "E", 3, _P3_E2_LE2_21_1T)
P3_LE2_E2_21_1T = OperatorTest("P3 L=E2 E2 1T (2,1)", "E", 3, _P3_LE2_E2_21_1T)

P4_A1_LA1_11 = OperatorTest("P4 A1 L=A1 (1,1)", "A1", 4, _P4_A1_LA1_11)
P4_A2_LA1_11 = OperatorTest("P4 A2 L=A1 (1,1)", "A2", 4, _P4_A2_LA1_11)
P4_LA1_A2_11 = OperatorTest("P4 L=A1 A2 (1,1)", "A1", 4, _P4_LA1_A2_11)
P4_E2_LA1_11 = OperatorTest("P4 E2 L=A1 (1,1)", "E", 4, _P4_E2_LA1_11)
P4_LA1_E2_11 = OperatorTest("P4 L=A1 E2 (1,1)", "A1", 4, _P4_LA1_E2_11)

