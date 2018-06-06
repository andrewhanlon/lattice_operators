from enum import Enum, auto

from sympy.tensor.array import Array
from sympy.matrices import Matrix, eye
from sympy.functions.special.tensor_functions import LeviCivita

from sympy import I

# Code for Gamma matrices
# Enumeration to specify different gamma matrix representations
class GammaRep(Enum):
  DIRAC_PAULI = auto()
  WEYL_CHIRAL = auto()
  DEGRAND_ROSSI = auto()

# explicit definition of gamma matrices
# The Dirac-Pauli Representation
_gamma_dp = [
    # gamma1
    Matrix([
        [ 0 , 0 , 0 ,-I ],
        [ 0 , 0 ,-I , 0 ],
        [ 0 , I , 0 , 0 ],
        [ I , 0 , 0 , 0 ]
    ]),
    # gamma2
    Matrix([
        [ 0 , 0 , 0 ,-1 ],
        [ 0 , 0 , 1 , 0 ],
        [ 0 , 1 , 0 , 0 ],
        [-1 , 0 , 0 , 0 ]
    ]),
    # gamma3
    Matrix([
        [ 0 , 0 ,-I , 0 ],
        [ 0 , 0 , 0 , I ],
        [ I , 0 , 0 , 0 ],
        [ 0 ,-I , 0 , 0 ]
    ]),
    # gamma4
    Matrix([
        [ 1 , 0 , 0 , 0 ],
        [ 0 , 1 , 0 , 0 ],
        [ 0 , 0 ,-1 , 0 ],
        [ 0 , 0 , 0 ,-1 ]
    ]),
    # gamma5
    Matrix([
        [ 0 , 0 , 1 , 0 ],
        [ 0 , 0 , 0 , 1 ],
        [ 1 , 0 , 0 , 0 ],
        [ 0 , 1 , 0 , 0 ]
    ])
]

# The Weyl chiral Representation
_gamma_wc = [
    # gamma1
    Matrix([
        [ 0 , 0 , 0 ,-I ],
        [ 0 , 0 ,-I , 0 ],
        [ 0 , I , 0 , 0 ],
        [ I , 0 , 0 , 0 ]
    ]),
    # gamma2
    Matrix([
        [ 0 , 0 , 0 ,-1 ],
        [ 0 , 0 , 1 , 0 ],
        [ 0 , 1 , 0 , 0 ],
        [-1 , 0 , 0 , 0 ]
    ]),
    # gamma3
    Matrix([
        [ 0 , 0 ,-I , 0 ],
        [ 0 , 0 , 0 , I ],
        [ I , 0 , 0 , 0 ],
        [ 0 ,-I , 0 , 0 ]
    ]),
    # gamma4
    Matrix([
        [ 0 , 0 , 1 , 0 ],
        [ 0 , 0 , 0 , 1 ],
        [ 1 , 0 , 0 , 0 ],
        [ 0 , 1 , 0 , 0 ]
    ]),
    # gamma5
    Matrix([
        [-1 , 0 , 0 , 0 ],
        [ 0 ,-1 , 0 , 0 ],
        [ 0 , 0 , 1 , 0 ],
        [ 0 , 0 , 0 , 1 ]
    ])
]

# The DeGrand-Rossi Representation
_gamma_dr = [
    # gamma1
    Matrix([
        [ 0 , 0 , 0 , I ],
        [ 0 , 0 , I , 0 ],
        [ 0 ,-I , 0 , 0 ],
        [-I , 0 , 0 , 0 ]
    ]),
    # gamma2
    Matrix([
        [ 0 , 0 , 0 ,-1 ],
        [ 0 , 0 , 1 , 0 ],
        [ 0 , 1 , 0 , 0 ],
        [-1 , 0 , 0 , 0 ]
    ]),
    # gamma3
    Matrix([
        [ 0 , 0 , I , 0 ],
        [ 0 , 0 , 0 ,-I ],
        [-I , 0 , 0 , 0 ],
        [ 0 , I , 0 , 0 ]
    ]),
    # gamma4
    Matrix([
        [ 0 , 0 , 1 , 0 ],
        [ 0 , 0 , 0 , 1 ],
        [ 1 , 0 , 0 , 0 ],
        [ 0 , 1 , 0 , 0 ]
    ]),
    # gamma5
    Matrix([
        [-1 , 0 , 0 , 0 ],
        [ 0 ,-1 , 0 , 0 ],
        [ 0 , 0 , 1 , 0 ],
        [ 0 , 0 , 0 , 1 ]
    ])
]


class Gamma:

  def __init__(self, rep=GammaRep.DIRAC_PAULI):
    if rep == GammaRep.DIRAC_PAULI:
      self._gamma = _gamma_dp
    elif rep == GammaRep.WEYL_CHIRAL:
      self._gamma = _gamma_wc
    elif rep == GammaRep.DEGRAND_ROSSI:
      self._gamma = _gamma_dr

    self._rep = rep

  @property
  def rep(self):
    return self._rep

  @property
  def one(self):
    return self._gamma[0]

  @property
  def two(self):
    return self._gamma[1]

  @property
  def three(self):
    return self._gamma[2]

  @property
  def four(self):
    return self._gamma[3]

  @property
  def five(self):
    return self._gamma[4]

  # @ADH - Is this how I want to do this?
  # @ADH - Also, is this expression representation dependent?
  @property
  def chargeConj(self):
    return self.four * self.two

  @property
  def parityPlus(self):
    return (eye(4) + self.four)/2

  @property
  def parityMinus(self):
    return (eye(4) - self.four)/2

