from sympy import Idx, Indexed, IndexedBase
from sympy.core.compatibility import is_sequence

from .grassmann import GrassmannVector
from .cubic_rotations import spinor_representation, E

class Quark(IndexedBase):

  is_commutative = False

  def __init__(self, name, shape=(3,4), antiquark=False):
    self._antiquark = antiquark
    if antiquark:
      name += "+"

    self._vectors = {E: GrassmannVector(name=name, shape=shape)}
    self._rotated_to = E

  @property
  def rotated_to(self):
    return self._rotated_to

  @property
  def vectors(self):
    return self._vectors

  @property
  def vector(self):
    return self.vectors[self.rotated_to]

  @property
  def antiquark(self):
    return self._antiquark

  def rotateTo(self, element):
    if element not in self.vectors:
      if self.antiquark:
        rot_mat = spinor_representation.rotation(element, True)
        self._vectors[element] = self.vectors[E].transformLeft(rot_mat, 1)
      else:
        rot_mat = spinor_representation.rotation(element, False)
        self._vectors[element] = self.vectors[E].transformRight(rot_mat, 1)

    self._rotated_to = element

  def __call__(self, *indices):
    return self.vectors[self.rotated_to][indices]

  def __getitem__(self, indices):
    if is_sequence(indices):
      # Special case needed because M[*my_tuple] is a syntax error.
      if self.shape and len(self.shape) != len(indices):
        raise IndexException("Rank mismatch.")
      return IndexedQuark(self, *indices)
    else:
      if self.shape and len(self.shape) != 1:
        raise IndexException("Rank mismatch.")
      return IndexedQuark(self, indices)

class IndexedQuark(Indexed):

  is_commutative = False


class DiracIdx(Idx):
  def __new__(cls, *args):
    try:
      return super().__new__(cls, args[0], 4)
    except IndexError:
      raise TypeError("DiracIdx requires a label")

class ColorIdx(Idx):
  def __new__(cls, *args):
    try:
      return super().__new__(cls, args[0], 3)
    except IndexError:
      raise TypeError("DiracIdx requires a label")



