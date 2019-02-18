from collections import defaultdict

from .grassmann import GrassmannVector, grassmann_coefficients, GrassmannTerm
from .cubic_rotations import spinor_representation, E, P0, P
from .operators import OperatorVector, OperatorElement

class OctetBaryon(OperatorVector):

  is_commutative = False

  bosonic = False
  fermionic = True
  statistics = -1

  def __init__(self, name, momentum=P0):
    """
    TODO:
      set shape somehow?
    """
    self._name = name
    self._momentum = momentum
    self._orig_momentum = momentum
    self._grassmann_vector = GrassmannVector(name=name, shape=(4,))
    self._vectors = {E: self._grassmann_vector.applyfunc(lambda i: OctetBaryonElement(i, momentum))}

    self._rotated_to = E

  def rotateTo(self, element):
    if element not in self._vectors:
      rot_mat = spinor_representation.rotation(element, False)
      trans_momentum = element * self._orig_momentum
      self._momentum = trans_momentum
      octet_vector = self._grassmann_vector.applyfunc(lambda i: OctetBaryonElement(i, trans_momentum))
      self._vectors[element] = octet_vector.transformRight(rot_mat, 0)

    self._rotated_to = element

  def P(self, *momentum):
    momentum = P(*momentum)
    if momentum == self._orig_momentum:
      self.rotateTo(E)
      return self

    return OctetBaryon(self.name, momentum)

  @property
  def name(self):
    return self._name



class OctetBaryonElement(OperatorElement):

  is_commutative = False

  bosonic = False
  fermionic = True
  statistics = -1

  def __init__(self, expr, momentum=P0):
    self._expr = expr
    self._coefficients = None
    self._momentum = momentum

  @property
  def coefficients(self):
    if self._coefficients is None:
      self._coefficients = defaultdict(int)
      for grassmann_term, coeff in grassmann_coefficients(self._expr).items():
        self._coefficients[(GrassmannTerm(grassmann_term, self.momentum),)] = coeff

    return self._coefficients



