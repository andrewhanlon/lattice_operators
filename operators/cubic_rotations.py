from enum import Enum, auto, unique
from math import gcd
from functools import reduce

from sortedcontainers import SortedSet

from sympy.functions.elementary.miscellaneous import sqrt
from sympy.matrices import eye, Identity, MatrixSymbol
from sympy import pi, Array, Matrix
from sympy import tensorcontraction, tensorproduct
from sympy import cos, sin, sqrt

from .tensors import GammaRep, Gamma


@unique
class Angle(Enum):
  E = 1
  HALF = 2
  THIRD = 3
  QUARTER = 4
  INV_THIRD = -3
  INV_QUARTER = -4

@unique
class Axis(Enum):
  X = auto()
  Y = auto()
  Z = auto()
  A = auto()
  B = auto()
  C = auto()
  D = auto()
  E = auto()
  F = auto()
  ALPHA = auto()
  BETA = auto()
  GAMMA = auto()
  DELTA = auto()

# @ADH - Force rotations to be in _POINT_GROUP ??
class CubicRotation:

  def __init__(self, rotation_angle, rotation_axis=Axis.Z, parity=False):
    self._rotation_angle = rotation_angle
    self._parity = parity
    self._matrix = None

    if rotation_angle == Angle.E:
      self._rotation_axis = Axis.Z
    else:
      self._rotation_axis = rotation_axis

  @property
  def angle(self):
    return self._rotation_angle

  @property
  def axis(self):
    return self._rotation_axis

  @property
  def parity(self):
    return self._parity

  @property
  def matrix(self):
    if self._matrix is None:
      a, b, c = _EULER_ANGLES[CubicRotation(self.angle, self.axis)]

      rot_mat = Matrix([
          [cos(a)*cos(b)*cos(c) - sin(a)*sin(c), -cos(a)*cos(b)*sin(c) - sin(a)*cos(c), cos(a)*sin(b)],
          [sin(a)*cos(b)*cos(c) + cos(a)*sin(c), -sin(a)*cos(b)*sin(c) + cos(a)*cos(c), sin(a)*sin(b)],
          [       -sin(b)*cos(c)               ,              sin(b)*sin(c)           ,    cos(b)    ]
      ])

      if self.parity:
        rot_mat *= -1

      self._matrix = rot_mat

    return self._matrix

  def inverse(self):
    mat_inv = self.matrix.inv()
    for rotation in _POINT_GROUP:
      if mat_inv == rotation.matrix:
        return rotation

    return None

  def __mul__(self, other):
    if isinstance(other, self.__class__):
      res_mat = self.matrix * other.matrix
      for rotation in _POINT_GROUP:
        if res_mat == rotation.matrix:
          return rotation

      return None

    elif isinstance(other, Momentum):
      return Momentum(tensorcontraction(tensorproduct(self.matrix, other), (1,2)))

    return NotImplemented
        
  def __str__(self):
    return _ROTATIONS[self]

  def __repr__(self):
    repr_str = "{}_{}".format(self.angle, self.axis)
    if self.parity:
      repr_str += "_I"

    return repr_str

  def __hash__(self):
    return hash(self.__repr__())

  def __eq__(self, other):
    if isinstance(other, self.__class__):
      return self.__repr__() == other.__repr__()
    return NotImplemented

  def __ne__(self, other):
    return not self.__eq__(other)

  def __lt__(self, other):
    if isinstance(other, self.__class__):
      return self.__repr__() < other.__repr__()
    return NotImplemented

  def __gt__(self, other):
    if isinstance(other, self.__class__):
      return self.__repr__() > other.__repr__()
    return NotImplemented

  def __le__(self, other):
    if isinstance(other, self.__class__):
      return self.__repr__() <= other.__repr__()
    return NotImplemented

  def __ge__(self, other):
    if isinstance(other, self.__class__):
      return self.__repr__() >= other.__repr__()
    return NotImplemented


E      = CubicRotation(Angle.E,           Axis.Z,     False)
C2x    = CubicRotation(Angle.HALF,        Axis.X,     False)
C2y    = CubicRotation(Angle.HALF,        Axis.Y,     False)
C2z    = CubicRotation(Angle.HALF,        Axis.Z,     False)
C2a    = CubicRotation(Angle.HALF,        Axis.A,     False)
C2b    = CubicRotation(Angle.HALF,        Axis.B,     False)
C2c    = CubicRotation(Angle.HALF,        Axis.C,     False)
C2d    = CubicRotation(Angle.HALF,        Axis.D,     False)
C2e    = CubicRotation(Angle.HALF,        Axis.E,     False)
C2f    = CubicRotation(Angle.HALF,        Axis.F,     False)
C3a    = CubicRotation(Angle.THIRD,       Axis.ALPHA, False)
C3b    = CubicRotation(Angle.THIRD,       Axis.BETA,  False)
C3c    = CubicRotation(Angle.THIRD,       Axis.GAMMA, False)
C3d    = CubicRotation(Angle.THIRD,       Axis.DELTA, False)
C3ai   = CubicRotation(Angle.INV_THIRD,   Axis.ALPHA, False)
C3bi   = CubicRotation(Angle.INV_THIRD,   Axis.BETA,  False)
C3ci   = CubicRotation(Angle.INV_THIRD,   Axis.GAMMA, False)
C3di   = CubicRotation(Angle.INV_THIRD,   Axis.DELTA, False)
C4x    = CubicRotation(Angle.QUARTER,     Axis.X,     False)
C4y    = CubicRotation(Angle.QUARTER,     Axis.Y,     False)
C4z    = CubicRotation(Angle.QUARTER,     Axis.Z,     False)
C4xi   = CubicRotation(Angle.INV_QUARTER, Axis.X,     False)
C4yi   = CubicRotation(Angle.INV_QUARTER, Axis.Y,     False)
C4zi   = CubicRotation(Angle.INV_QUARTER, Axis.Z,     False)
Is     = CubicRotation(Angle.E,           Axis.Z,     True)
I_C2x  = CubicRotation(Angle.HALF,        Axis.X,     True)
I_C2y  = CubicRotation(Angle.HALF,        Axis.Y,     True)
I_C2z  = CubicRotation(Angle.HALF,        Axis.Z,     True)
I_C2a  = CubicRotation(Angle.HALF,        Axis.A,     True)
I_C2b  = CubicRotation(Angle.HALF,        Axis.B,     True)
I_C2c  = CubicRotation(Angle.HALF,        Axis.C,     True)
I_C2d  = CubicRotation(Angle.HALF,        Axis.D,     True)
I_C2e  = CubicRotation(Angle.HALF,        Axis.E,     True)
I_C2f  = CubicRotation(Angle.HALF,        Axis.F,     True)
I_C3a  = CubicRotation(Angle.THIRD,       Axis.ALPHA, True)
I_C3b  = CubicRotation(Angle.THIRD,       Axis.BETA,  True)
I_C3c  = CubicRotation(Angle.THIRD,       Axis.GAMMA, True)
I_C3d  = CubicRotation(Angle.THIRD,       Axis.DELTA, True)
I_C3ai = CubicRotation(Angle.INV_THIRD,   Axis.ALPHA, True)
I_C3bi = CubicRotation(Angle.INV_THIRD,   Axis.BETA,  True)
I_C3ci = CubicRotation(Angle.INV_THIRD,   Axis.GAMMA, True)
I_C3di = CubicRotation(Angle.INV_THIRD,   Axis.DELTA, True)
I_C4x  = CubicRotation(Angle.QUARTER,     Axis.X,     True)
I_C4y  = CubicRotation(Angle.QUARTER,     Axis.Y,     True)
I_C4z  = CubicRotation(Angle.QUARTER,     Axis.Z,     True)
I_C4xi = CubicRotation(Angle.INV_QUARTER, Axis.X,     True)
I_C4yi = CubicRotation(Angle.INV_QUARTER, Axis.Y,     True)
I_C4zi = CubicRotation(Angle.INV_QUARTER, Axis.Z,     True)

# @ADH - These need to be done for each LittleGroup...
_GENERATORS = {
    E:      [],
    C2x:    [C4x,  C4x],
    C2y:    [C4y,  C4y],
    C2z:    [C4z,  C4z],
    C2a:    [C2y,  C4z],
    C2b:    [C2x,  C4z],
    C2c:    [C4y,  C2z],
    C2d:    [C2z,  C4y],
    C2e:    [C2z,  C4x],
    C2f:    [C2y,  C4x],
    C3a:    [C4yi, C4z],
    C3b:    [C4y,  C4zi],
    C3c:    [C4yi, C4zi],
    C3d:    [C4y,  C4z],
    C3ai:   [C4zi, C4y],
    C3bi:   [C4z,  C4yi],
    C3ci:   [C4z,  C4y],
    C3di:   [C4zi, C4yi],
    C4x:    [C4zi, C4y, C4z],
    C4y:    [],
    C4z:    [],
    C4xi:   "invert",
    C4yi:   "invert",
    C4zi:   "invert",
    Is:     [],
    I_C2x:  [Is, C2x],
    I_C2y:  [Is, C2y],
    I_C2z:  [Is, C2z],
    I_C2a:  [Is, C2a],
    I_C2b:  [Is, C2b],
    I_C2c:  [Is, C2c],
    I_C2d:  [Is, C2d],
    I_C2e:  [Is, C2e],
    I_C2f:  [Is, C2f],
    I_C3a:  [Is, C3a],
    I_C3b:  [Is, C3b],
    I_C3c:  [Is, C3c],
    I_C3d:  [Is, C3d],
    I_C3ai: [Is, C3ai],
    I_C3bi: [Is, C3bi],
    I_C3ci: [Is, C3ci],
    I_C3di: [Is, C3di],
    I_C4x:  [Is, C4x],
    I_C4y:  [Is, C4y],
    I_C4z:  [Is, C4z],
    I_C4xi: [Is, C4xi],
    I_C4yi: [Is, C4yi],
    I_C4zi: [Is, C4zi]
}


_ROTATIONS = {
    E:      "E",
    C2x:    "C_{2x}",
    C2y:    "C_{2y}",
    C2z:    "C_{2z}",
    C2a:    "C_{2a}",
    C2b:    "C_{2b}",
    C2c:    "C_{2c}",
    C2d:    "C_{2d}",
    C2e:    "C_{2e}",
    C2f:    "C_{2f}",
    C3a:    "C_{3A}",
    C3b:    "C_{3B}",
    C3c:    "C_{3Y}",
    C3d:    "C_{3D}",
    C3ai:   "C_{3A}^{-1}",
    C3bi:   "C_{3B}^{-1}",
    C3ci:   "C_{3Y}^{-1}",
    C3di:   "C_{3D}^{-1}",
    C4x:    "C_{4x}",
    C4y:    "C_{4y}",
    C4z:    "C_{4z}",
    C4xi:   "C_{4x}^{-1}",
    C4yi:   "C_{4y}^{-1}",
    C4zi:   "C_{4z}^{-1}",
    Is:      "I_S",
    I_C2x:  "I_S C_{2x}",
    I_C2y:  "I_S C_{2y}",
    I_C2z:  "I_S C_{2z}",
    I_C2a:  "I_S C_{2a}",
    I_C2b:  "I_S C_{2b}",
    I_C2c:  "I_S C_{2c}",
    I_C2d:  "I_S C_{2d}",
    I_C2e:  "I_S C_{2e}",
    I_C2f:  "I_S C_{2f}",
    I_C3a:  "I_S C_{3A}",
    I_C3b:  "I_S C_{3B}",
    I_C3c:  "I_S C_{3Y}",
    I_C3d:  "I_S C_{3D}",
    I_C3ai: "I_S C_{3A}^{-1}",
    I_C3bi: "I_S C_{3B}^{-1}",
    I_C3ci: "I_S C_{3Y}^{-1}",
    I_C3di: "I_S C_{3D}^{-1}",
    I_C4x:  "I_S C_{4x}",
    I_C4y:  "I_S C_{4y}",
    I_C4z:  "I_S C_{4z}",
    I_C4xi: "I_S C_{4x}^{-1}",
    I_C4yi: "I_S C_{4y}^{-1}",
    I_C4zi: "I_S C_{4z}^{-1}"
}


_EULER_ANGLES = {
    E:    (      0,      0,      0),
    C2x:  (  -pi/2,     pi,   pi/2),
    C2y:  (      0,     pi,      0),
    C2z:  (   pi/2,      0,   pi/2),
    C2a:  (  -pi/4,     pi,   pi/4),
    C2b:  (-3*pi/4,     pi, 3*pi/4),
    C2c:  (      0,   pi/2,     pi),
    C2d:  (     pi,   pi/2,      0),
    C2e:  (   pi/2,   pi/2,   pi/2),
    C2f:  (  -pi/2,   pi/2,  -pi/2),
    C3a:  (     pi,   pi/2,  -pi/2),
    C3b:  (      0,   pi/2,  -pi/2),
    C3c:  (    -pi,   pi/2,   pi/2),
    C3d:  (      0,   pi/2,   pi/2),
    C3ai: (  -pi/2,   pi/2,      0),
    C3bi: (  -pi/2,   pi/2,     pi),
    C3ci: (   pi/2,   pi/2,      0),
    C3di: (   pi/2,   pi/2,    -pi),
    C4x:  (  -pi/2,   pi/2,   pi/2),
    C4y:  (      0,   pi/2,      0),
    C4z:  (   pi/4,      0,   pi/4),
    C4xi: (   pi/2,   pi/2,  -pi/2),
    C4yi: (    -pi,   pi/2,     pi),
    C4zi: (  -pi/4,      0,  -pi/4)
}


_OCTAHEDRAL_GROUP = SortedSet([cubic_rotation for cubic_rotation in _ROTATIONS.keys() if not cubic_rotation.parity])
_POINT_GROUP = SortedSet([cubic_rotation for cubic_rotation in _ROTATIONS.keys()])


class Momentum(Array):

  def __init__(self, momentum=[0,0,0]):
    self._p_x = momentum[0]
    self._p_y = momentum[1]
    self._p_z = momentum[2]

  @property
  def x(self):
    return self._p_x

  @property
  def y(self):
    return self._p_y

  @property
  def z(self):
    return self._p_z

  @property
  def psq(self):
    return self._p_x**2 + self._p_y**2 + self._p_z**2

  @property
  def pref(self):
    return Momentum(sorted([abs(pi) for pi in self]))

  @property
  def reduced_pref(self):
    if self.psq == 0:
      return self

    factor = reduce(lambda x,y: gcd(x,y), self.pref)
    return Momentum([self.pref.x//factor, self.pref.y//factor, self.pref.z//factor])

  @property
  def reduced(self):
    if self.psq == 0:
      return self

    factor = reduce(lambda x,y: gcd(x,y), self)
    return Momentum([self.x//factor, self.y//factor, self.z//factor])

  def __str__(self):
    return "P=({},{},{})".format(self.x, self.y, self.z)

  def __repr__(self):
    return "P{{{}_{}_{}}}".format(self.x, self.y, self.z)

  def __hash__(self):
    return hash(self.__repr__())

  def __eq__(self, other):
    if isinstance(other, self.__class__):
      return self.__repr__() == other.__repr__()
    return NotImplemented

  def __ne__(self, other):
    return not self.__eq__(other)

  def __neg__(self):
    return Momentum([-self.x, -self.y, -self.z])

  def __add__(self, other):
    if isinstance(other, self.__class__):
      return Momentum([self.x + other.x, self.y + other.y, self.z + other.z])
    return NotImplemented
  
  def __sub__(self, other):
    return self.__add__(-other)

  def __mul__(self, other):
    if isinstance(other, self.__class__):
      px = self.y*other.z - self.z*other.y
      py = self.z*other.x - self.x*other.z
      pz = self.x*other.y - self.y*other.x
      return Momentum([px, py, pz])

    else:
      return Momentum([self.x * other, self.y * other, self.z * other])

  def __rmul(self, other):
    return Momentum([other * self.x, other * self.y, other * self.z])

  def __truediv__(self, other):
    return Momentum([self.x // other, self.y // other, self.z // other])

  def __rtruediv__(self, other):
    return Momentum([self.x // other, self.y // other, self.z // other])


P = Momentum
P0 = Momentum([0,0,0])

# @ADH - add the C_S refs
# @ADH - Is this cheating?
_REFERENCE_ROTATIONS = {
    P([ 0, 0, 0]): E,
    P([ 0, 0, 1]): E,
    P([ 0, 0,-1]): C2x,
    P([ 1, 0, 0]): C4y,
    P([-1, 0, 0]): C4yi,
    P([ 0,-1, 0]): C4x,
    P([ 0, 1, 0]): C4xi,
    P([ 0, 1, 1]): E,
    P([ 0,-1,-1]): C2x,
    P([ 0, 1,-1]): C4xi,
    P([ 0,-1, 1]): C4x,
    P([ 1, 0, 1]): C4zi,
    P([-1, 0,-1]): C2b,
    P([ 1, 0,-1]): C2a,
    P([-1, 0, 1]): C4z,
    P([ 1, 1, 0]): C4y,
    P([-1,-1, 0]): C2d,
    P([ 1,-1, 0]): C2c,
    P([-1, 1, 0]): C4yi,
    P([ 1, 1, 1]): E,
    P([ 1, 1,-1]): C4y,
    P([ 1,-1, 1]): C4x,
    P([ 1,-1,-1]): C2x,
    P([-1, 1, 1]): C4z,
    P([-1, 1,-1]): C2y,
    P([-1,-1, 1]): C2z,
    P([-1,-1,-1]): C2d,
}

_BOSONIC_LITTLE_GROUP_IRREPS = {
    P([0,0,0]): ["A1g", "A2g", "Eg", "T1g", "T2g", "A1u", "A2u", "Eu", "T1u", "T2u"],
    P([0,0,1]): ["A1", "A2", "B1", "B2", "E"],
    P([0,1,1]): ["A1", "A2", "B1", "B2"],
    P([1,1,1]): ["A1", "A2", "E"],
    P([0,1,2]): ["A1", "A2"],
    P([1,1,2]): ["A1", "A2"]
}

_BOSONIC_LITTLE_GROUPS = {
    P([0,0,0]): "O_h",
    P([0,0,1]): "C_{4v}",
    P([0,1,1]): "C_{2v}",
    P([1,1,1]): "C_{3v}",
    P([0,1,2]): "C_S",
    P([1,1,2]): "C_S"
}

_FERMIONIC_LITTLE_GROUP_IRREPS = {
    P([0,0,0]): ["G1g", "G2g", "Hg", "G1u", "G2u", "Hu"],
    P([0,0,1]): ["G1", "G2"],
    P([0,1,1]): ["G"],
    P([1,1,1]): ["F1", "F2", "G"],
    P([0,1,2]): ["F1", "F2"],
    P([1,1,2]): ["F1", "F2"]
}

_FERMIONIC_LITTLE_GROUPS = {
    P([0,0,0]): "O_h^D",
    P([0,0,1]): "C_{4v}^D",
    P([0,1,1]): "C_{2v}^D",
    P([1,1,1]): "C_{3v}^D",
    P([0,1,2]): "C_S^D",
    P([1,1,2]): "C_S^D"
}

_LITTLE_GROUPS = {
    P([0,0,0]): "Oh",
    P([0,0,1]): "C4v",
    P([0,1,1]): "C2v",
    P([1,1,1]): "C3v",
    P([0,1,2]): "CS",
    P([1,1,2]): "CS"
}


# Conjugacy Classes
Oh_1  = frozenset([E])
Oh_2  = frozenset([C3a, C3b, C3c, C3d, C3ai, C3bi, C3ci, C3di])
Oh_3  = frozenset([C2x, C2y, C2z])
Oh_4  = frozenset([C4x, C4y, C4z, C4xi, C4yi, C4zi])
Oh_5  = frozenset([C2a, C2b, C2c, C2d, C2e, C2f])
Oh_6  = frozenset([Is])
Oh_7  = frozenset([I_C3a, I_C3b, I_C3c, I_C3d, I_C3ai, I_C3bi, I_C3ci, I_C3di])
Oh_8  = frozenset([I_C2x, I_C2y, I_C2z])
Oh_9  = frozenset([I_C4x, I_C4y, I_C4z, I_C4xi, I_C4yi, I_C4zi])
Oh_10 = frozenset([I_C2a, I_C2b, I_C2c, I_C2d, I_C2e, I_C2f])

C4v_1 = frozenset([E])
C4v_2 = frozenset([C2z])
C4v_3 = frozenset([C4z, C4zi])
C4v_4 = frozenset([I_C2x, I_C2y])
C4v_5 = frozenset([I_C2a, I_C2b])

C2v_1 = frozenset([E])
C2v_2 = frozenset([C2e])
C2v_3 = frozenset([I_C2f])
C2v_4 = frozenset([I_C2x])

C3v_1 = frozenset([E])
C3v_2 = frozenset([C3d, C3di])
C3v_3 = frozenset([I_C2b, I_C2d, I_C2f])

Cs012_1 = frozenset([E])
Cs012_2 = frozenset([I_C2x])

Cs112_1 = frozenset([E])
Cs112_2 = frozenset([I_C2b])

_CHARACTERS = {
    ("Oh",  "A1g", Oh_1):  1,
    ("Oh",  "A1g", Oh_2):  1,
    ("Oh",  "A1g", Oh_3):  1,
    ("Oh",  "A1g", Oh_4):  1,
    ("Oh",  "A1g", Oh_5):  1,
    ("Oh",  "A1g", Oh_6):  1,
    ("Oh",  "A1g", Oh_7):  1,
    ("Oh",  "A1g", Oh_8):  1,
    ("Oh",  "A1g", Oh_9):  1,
    ("Oh",  "A1g", Oh_10): 1,
    ("Oh",  "A2g", Oh_1):  1,
    ("Oh",  "A2g", Oh_2):  1,
    ("Oh",  "A2g", Oh_3):  1,
    ("Oh",  "A2g", Oh_4):  -1,
    ("Oh",  "A2g", Oh_5):  -1,
    ("Oh",  "A2g", Oh_6):  1,
    ("Oh",  "A2g", Oh_7):  1,
    ("Oh",  "A2g", Oh_8):  1,
    ("Oh",  "A2g", Oh_9):  -1,
    ("Oh",  "A2g", Oh_10): -1,
    ("Oh",  "Eg",  Oh_1):  2,
    ("Oh",  "Eg",  Oh_2):  -1,
    ("Oh",  "Eg",  Oh_3):  2,
    ("Oh",  "Eg",  Oh_4):  0,
    ("Oh",  "Eg",  Oh_5):  0,
    ("Oh",  "Eg",  Oh_6):  2,
    ("Oh",  "Eg",  Oh_7):  -1,
    ("Oh",  "Eg",  Oh_8):  2,
    ("Oh",  "Eg",  Oh_9):  0,
    ("Oh",  "Eg",  Oh_10): 0,
    ("Oh",  "T1g", Oh_1):  3,
    ("Oh",  "T1g", Oh_2):  0,
    ("Oh",  "T1g", Oh_3):  -1,
    ("Oh",  "T1g", Oh_4):  1,
    ("Oh",  "T1g", Oh_5):  -1,
    ("Oh",  "T1g", Oh_6):  3,
    ("Oh",  "T1g", Oh_7):  0,
    ("Oh",  "T1g", Oh_8):  -1,
    ("Oh",  "T1g", Oh_9):  1,
    ("Oh",  "T1g", Oh_10): -1,
    ("Oh",  "T2g", Oh_1):  3,
    ("Oh",  "T2g", Oh_2):  0,
    ("Oh",  "T2g", Oh_3):  -1,
    ("Oh",  "T2g", Oh_4):  -1,
    ("Oh",  "T2g", Oh_5):  1,
    ("Oh",  "T2g", Oh_6):  3,
    ("Oh",  "T2g", Oh_7):  0,
    ("Oh",  "T2g", Oh_8):  -1,
    ("Oh",  "T2g", Oh_9):  -1,
    ("Oh",  "T2g", Oh_10): 1,
    ("Oh",  "G1g", Oh_1):  2,
    ("Oh",  "G1g", Oh_2):  1,
    ("Oh",  "G1g", Oh_3):  0,
    ("Oh",  "G1g", Oh_4):  sqrt(2),
    ("Oh",  "G1g", Oh_5):  0,
    ("Oh",  "G1g", Oh_6):  2,
    ("Oh",  "G1g", Oh_7):  1,
    ("Oh",  "G1g", Oh_8):  0,
    ("Oh",  "G1g", Oh_9):  sqrt(2),
    ("Oh",  "G1g", Oh_10): 0,
    ("Oh",  "G2g", Oh_1):  2,
    ("Oh",  "G2g", Oh_2):  1,
    ("Oh",  "G2g", Oh_3):  0,
    ("Oh",  "G2g", Oh_4):  -sqrt(2),
    ("Oh",  "G2g", Oh_5):  0,
    ("Oh",  "G2g", Oh_6):  2,
    ("Oh",  "G2g", Oh_7):  1,
    ("Oh",  "G2g", Oh_8):  0,
    ("Oh",  "G2g", Oh_9):  -sqrt(2),
    ("Oh",  "G2g", Oh_10): 0,
    ("Oh",  "Hg",  Oh_1):  4,
    ("Oh",  "Hg",  Oh_2):  -1,
    ("Oh",  "Hg",  Oh_3):  0,
    ("Oh",  "Hg",  Oh_4):  0,
    ("Oh",  "Hg",  Oh_5):  0,
    ("Oh",  "Hg",  Oh_6):  4,
    ("Oh",  "Hg",  Oh_7):  -1,
    ("Oh",  "Hg",  Oh_8):  0,
    ("Oh",  "Hg",  Oh_9):  0,
    ("Oh",  "Hg",  Oh_10): 0,
    ("Oh",  "A1u", Oh_1):  1,
    ("Oh",  "A1u", Oh_2):  1,
    ("Oh",  "A1u", Oh_3):  1,
    ("Oh",  "A1u", Oh_4):  1,
    ("Oh",  "A1u", Oh_5):  1,
    ("Oh",  "A1u", Oh_6):  -1,
    ("Oh",  "A1u", Oh_7):  -1,
    ("Oh",  "A1u", Oh_8):  -1,
    ("Oh",  "A1u", Oh_9):  -1,
    ("Oh",  "A1u", Oh_10): -1,
    ("Oh",  "A2u", Oh_1):  1,
    ("Oh",  "A2u", Oh_2):  1,
    ("Oh",  "A2u", Oh_3):  1,
    ("Oh",  "A2u", Oh_4):  -1,
    ("Oh",  "A2u", Oh_5):  -1,
    ("Oh",  "A2u", Oh_6):  -1,
    ("Oh",  "A2u", Oh_7):  -1,
    ("Oh",  "A2u", Oh_8):  -1,
    ("Oh",  "A2u", Oh_9):  1,
    ("Oh",  "A2u", Oh_10): 1,
    ("Oh",  "Eu",  Oh_1):  2,
    ("Oh",  "Eu",  Oh_2):  -1,
    ("Oh",  "Eu",  Oh_3):  2,
    ("Oh",  "Eu",  Oh_4):  0,
    ("Oh",  "Eu",  Oh_5):  0,
    ("Oh",  "Eu",  Oh_6):  -2,
    ("Oh",  "Eu",  Oh_7):  1,
    ("Oh",  "Eu",  Oh_8):  -2,
    ("Oh",  "Eu",  Oh_9):  0,
    ("Oh",  "Eu",  Oh_10): 0,
    ("Oh",  "T1u", Oh_1):  3,
    ("Oh",  "T1u", Oh_2):  0,
    ("Oh",  "T1u", Oh_3):  -1,
    ("Oh",  "T1u", Oh_4):  1,
    ("Oh",  "T1u", Oh_5):  -1,
    ("Oh",  "T1u", Oh_6):  -3,
    ("Oh",  "T1u", Oh_7):  0,
    ("Oh",  "T1u", Oh_8):  1,
    ("Oh",  "T1u", Oh_9):  -1,
    ("Oh",  "T1u", Oh_10): 1,
    ("Oh",  "T2u", Oh_1):  3,
    ("Oh",  "T2u", Oh_2):  0,
    ("Oh",  "T2u", Oh_3):  -1,
    ("Oh",  "T2u", Oh_4):  -1,
    ("Oh",  "T2u", Oh_5):  1,
    ("Oh",  "T2u", Oh_6):  -3,
    ("Oh",  "T2u", Oh_7):  0,
    ("Oh",  "T2u", Oh_8):  1,
    ("Oh",  "T2u", Oh_9):  1,
    ("Oh",  "T2u", Oh_10): -1,
    ("Oh",  "G1u", Oh_1):  2,
    ("Oh",  "G1u", Oh_2):  1,
    ("Oh",  "G1u", Oh_3):  0,
    ("Oh",  "G1u", Oh_4):  sqrt(2),
    ("Oh",  "G1u", Oh_5):  0,
    ("Oh",  "G1u", Oh_6):  -2,
    ("Oh",  "G1u", Oh_7):  -1,
    ("Oh",  "G1u", Oh_8):  0,
    ("Oh",  "G1u", Oh_9):  -sqrt(2),
    ("Oh",  "G1u", Oh_10): 0,
    ("Oh",  "G2u", Oh_1):  2,
    ("Oh",  "G2u", Oh_2):  1,
    ("Oh",  "G2u", Oh_3):  0,
    ("Oh",  "G2u", Oh_4):  -sqrt(2),
    ("Oh",  "G2u", Oh_5):  0,
    ("Oh",  "G2u", Oh_6):  -2,
    ("Oh",  "G2u", Oh_7):  -1,
    ("Oh",  "G2u", Oh_8):  0,
    ("Oh",  "G2u", Oh_9):  sqrt(2),
    ("Oh",  "G2u", Oh_10): 0,
    ("Oh",  "Hu",  Oh_1):  4,
    ("Oh",  "Hu",  Oh_2):  -1,
    ("Oh",  "Hu",  Oh_3):  0,
    ("Oh",  "Hu",  Oh_4):  0,
    ("Oh",  "Hu",  Oh_5):  0,
    ("Oh",  "Hu",  Oh_6):  -4,
    ("Oh",  "Hu",  Oh_7):  1,
    ("Oh",  "Hu",  Oh_8):  0,
    ("Oh",  "Hu",  Oh_9):  0,
    ("Oh",  "Hu",  Oh_10): 0,
    ("C4v", "A1",  C4v_1): 1,
    ("C4v", "A1",  C4v_2): 1,
    ("C4v", "A1",  C4v_3): 1,
    ("C4v", "A1",  C4v_4): 1,
    ("C4v", "A1",  C4v_5): 1,
    ("C4v", "A2",  C4v_1): 1,
    ("C4v", "A2",  C4v_2): 1,
    ("C4v", "A2",  C4v_3): 1,
    ("C4v", "A2",  C4v_4): -1,
    ("C4v", "A2",  C4v_5): -1,
    ("C4v", "B1",  C4v_1): 1,
    ("C4v", "B1",  C4v_2): 1,
    ("C4v", "B1",  C4v_3): -1,
    ("C4v", "B1",  C4v_4): 1,
    ("C4v", "B1",  C4v_5): -1,
    ("C4v", "B2",  C4v_1): 1,
    ("C4v", "B2",  C4v_2): 1,
    ("C4v", "B2",  C4v_3): -1,
    ("C4v", "B2",  C4v_4): -1,
    ("C4v", "B2",  C4v_5): 1,
    ("C4v", "E",   C4v_1): 2,
    ("C4v", "E",   C4v_2): -2,
    ("C4v", "E",   C4v_3): 0,
    ("C4v", "E",   C4v_4): 0,
    ("C4v", "E",   C4v_5): 0,
    ("C4v", "G1",  C4v_1): 2,
    ("C4v", "G1",  C4v_2): 0,
    ("C4v", "G1",  C4v_3): sqrt(2),
    ("C4v", "G1",  C4v_4): 0,
    ("C4v", "G1",  C4v_5): 0,
    ("C4v", "G2",  C4v_1): 2,
    ("C4v", "G2",  C4v_2): 0,
    ("C4v", "G2",  C4v_3): -sqrt(2),
    ("C4v", "G2",  C4v_4): 0,
    ("C4v", "G2",  C4v_5): 0,
    ("C2v", "A1",  C2v_1): 1,
    ("C2v", "A1",  C2v_2): 1,
    ("C2v", "A1",  C2v_3): 1,
    ("C2v", "A1",  C2v_4): 1,
    ("C2v", "A2",  C2v_1): 1,
    ("C2v", "A2",  C2v_2): 1,
    ("C2v", "A2",  C2v_3): -1,
    ("C2v", "A2",  C2v_4): -1,
    ("C2v", "B1",  C2v_1): 1,
    ("C2v", "B1",  C2v_2): -1,
    ("C2v", "B1",  C2v_3): 1,
    ("C2v", "B1",  C2v_4): -1,
    ("C2v", "B2",  C2v_1): 1,
    ("C2v", "B2",  C2v_2): -1,
    ("C2v", "B2",  C2v_3): -1,
    ("C2v", "B2",  C2v_4): 1,
    ("C2v", "G",   C2v_1): 2,
    ("C2v", "G",   C2v_2): 0,
    ("C2v", "G",   C2v_3): 0,
    ("C2v", "G",   C2v_4): 0,
    ("C3v", "A1",  C3v_1): 1,
    ("C3v", "A1",  C3v_2): 1,
    ("C3v", "A1",  C3v_3): 1,
    ("C3v", "A2",  C3v_1): 1,
    ("C3v", "A2",  C3v_2): 1,
    ("C3v", "A2",  C3v_3): -1,
    ("C3v", "E",   C3v_1): 2,
    ("C3v", "E",   C3v_2): -1,
    ("C3v", "E",   C3v_3): 0,
    ("C3v", "F1",  C3v_1): 1,
    ("C3v", "F1",  C3v_2): -1,
    ("C3v", "F1",  C3v_3): 1j,
    ("C3v", "F2",  C3v_1): 1,
    ("C3v", "F2",  C3v_2): -1,
    ("C3v", "F2",  C3v_3): -1j,
    ("C3v", "G",   C3v_1): 2,
    ("C3v", "G",   C3v_2): 1,
    ("C3v", "G",   C3v_3): 0
}

class LittleGroup:

  def __init__(self, bosonic, momentum=P0):
    self._momentum = momentum
    self._bosonic = bosonic
    self._elements = set()
    self._ref_elements = dict()
    self._conj_class = dict()

  @property
  def order(self):
    return len(self.elements)

  @property
  def momentum(self):
    return self._momentum

  @property
  def bosonic(self):
    return self._bosonic

  @property
  def fermionic(self):
    return not self._bosonic

  def getCharacter(self, irrep, element):
    if element not in self.elements:
      raise ValueError("Element not in Little Group")

    if irrep not in self.irreps:
      raise ValueError("Not a Little Group irrep")

    conj_class = self.getConjugacyClass(element)
    conj_class_ref = frozenset([self.reference_element(el) for el in conj_class])
    return _CHARACTERS[(self.little_group, irrep, conj_class_ref)]

  def getConjugacyClass(self, element):
    if element not in self.elements:
      raise ValueError("Element not in Little Group")

    if element in self._conj_class:
      return self._conj_class[element]

    conj_class = set()
    for el in self.elements:
      conj_class.add(el*element*el.inverse())

    conj_class = frozenset(conj_class)
    for el in conj_class:
      self._conj_class[el] = conj_class

    return conj_class
  

  @property
  def elements(self):
    if not self._elements:
      for rotation in _POINT_GROUP:
        mom_prime = rotation * self.momentum
        if mom_prime == self.momentum:
          self._elements.add(rotation)

    return self._elements

  def reference_element(self, element):
    if not self._ref_elements:
      self._make_ref_elements()

    return self._ref_elements[element]

  def _make_ref_elements(self):
    ref_rotation = _REFERENCE_ROTATIONS[self.momentum.reduced]
    for element in self.elements:
      ref_element = ref_rotation.inverse() * element * ref_rotation
      self._ref_elements[element] = ref_element

  @property
  def irreps(self):
    if self.bosonic:
      return _BOSONIC_LITTLE_GROUP_IRREPS[self.momentum.reduced_pref]
    else:
      return _FERMIONIC_LITTLE_GROUP_IRREPS[self.momentum.reduced_pref]

  @property
  def little_group(self):
    return _LITTLE_GROUPS[self.momentum.reduced_pref]


  def __str__(self):
    if self.bosonic:
      return _BOSONIC_LITTLE_GROUPS[self.momentum.reduced_pref]
    else:
      return _FERMIONIC_LITTLE_GROUPS[self.momentum.reduced_pref]



# Spinor Representations
I = Identity(4)
g1 = MatrixSymbol('g1', 4, 4)
g2 = MatrixSymbol('g2', 4, 4)
g3 = MatrixSymbol('g3', 4, 4)

_rotation_map = {
    Angle.HALF: {
        Axis.X: g2*g3,
        Axis.Y: g3*g1,
        Axis.Z: g1*g2,
        Axis.A:  1/sqrt(2)*(g2*g3 + g3*g1),
        Axis.B:  1/sqrt(2)*(g2*g3 - g3*g1),
        Axis.C:  1/sqrt(2)*(g2*g3 + g1*g2),
        Axis.D: -1/sqrt(2)*(g2*g3 - g1*g2),
        Axis.E:  1/sqrt(2)*(g3*g1 + g1*g2),
        Axis.F:  1/sqrt(2)*(g3*g1 - g1*g2)
    },
    Angle.THIRD: {
        Axis.ALPHA: (I - g2*g3 - g3*g1 + g1*g2)/2,
        Axis.BETA:  (I - g2*g3 + g3*g1 - g1*g2)/2,
        Axis.GAMMA: (I + g2*g3 - g3*g1 - g1*g2)/2,
        Axis.DELTA: (I + g2*g3 + g3*g1 + g1*g2)/2
    },
    Angle.INV_THIRD: {
        Axis.ALPHA: (I + g2*g3 + g3*g1 - g1*g2)/2,
        Axis.BETA:  (I + g2*g3 - g3*g1 + g1*g2)/2,
        Axis.GAMMA: (I - g2*g3 + g3*g1 + g1*g2)/2,
        Axis.DELTA: (I - g2*g3 - g3*g1 - g1*g2)/2
    },
    Angle.QUARTER: {
        Axis.X: 1/sqrt(2)*(I + g2*g3),
        Axis.Y: 1/sqrt(2)*(I + g3*g1),
        Axis.Z: 1/sqrt(2)*(I + g1*g2)
    },
    Angle.INV_QUARTER: {
        Axis.X: 1/sqrt(2)*(I - g2*g3),
        Axis.Y: 1/sqrt(2)*(I - g3*g1),
        Axis.Z: 1/sqrt(2)*(I - g1*g2)
    }
}

_conjugate_rotation_map = {
    Angle.HALF: {
        Axis.X: -g2*g3,
        Axis.Y: -g3*g1,
        Axis.Z: -g1*g2,
        Axis.A: -1/sqrt(2)*(g2*g3 + g3*g1),
        Axis.B: -1/sqrt(2)*(g2*g3 - g3*g1),
        Axis.C: -1/sqrt(2)*(g2*g3 + g1*g2),
        Axis.D:  1/sqrt(2)*(g2*g3 - g1*g2),
        Axis.E: -1/sqrt(2)*(g3*g1 + g1*g2),
        Axis.F: -1/sqrt(2)*(g3*g1 - g1*g2)
    },
    Angle.THIRD: {
        Axis.ALPHA: (I + g2*g3 + g3*g1 - g1*g2)/2,
        Axis.BETA:  (I + g2*g3 - g3*g1 + g1*g2)/2,
        Axis.GAMMA: (I - g2*g3 + g3*g1 + g1*g2)/2,
        Axis.DELTA: (I - g2*g3 - g3*g1 - g1*g2)/2
    },
    Angle.INV_THIRD: {
        Axis.ALPHA: (I - g2*g3 - g3*g1 + g1*g2)/2,
        Axis.BETA:  (I - g2*g3 + g3*g1 - g1*g2)/2,
        Axis.GAMMA: (I + g2*g3 - g3*g1 - g1*g2)/2,
        Axis.DELTA: (I + g2*g3 + g3*g1 + g1*g2)/2
    },
    Angle.QUARTER: {
        Axis.X: 1/sqrt(2)*(I - g2*g3),
        Axis.Y: 1/sqrt(2)*(I - g3*g1),
        Axis.Z: 1/sqrt(2)*(I - g1*g2)
    },
    Angle.INV_QUARTER: {
        Axis.X: 1/sqrt(2)*(I + g2*g3),
        Axis.Y: 1/sqrt(2)*(I + g3*g1),
        Axis.Z: 1/sqrt(2)*(I + g1*g2)
    }
}


class SpinorRepresentation:

  def __init__(self, gamma_rep=GammaRep.DIRAC_PAULI):
    self._gamma = Gamma(gamma_rep)

    self._setup_rotations()

  @property
  def gamma(self):
    return self._gamma

  @property
  def gammaRep(self):
    return self.gamma.rep

  @gammaRep.setter
  def gammaRep(self, gamma_rep):
    if self.gammaRep != gamma_rep:
      self._gamma = Gamma(gamma_rep)
      self._setup_rotations()

  def rotation(self, cubic_rotation, conjugate, double_element=False):

    if cubic_rotation.angle == Angle.E:
      repr_mat = eye(4)
    elif conjugate:
      repr_mat = self.conj_rotations[cubic_rotation.angle][cubic_rotation.axis]
    else:
      repr_mat = self.rotations[cubic_rotation.angle][cubic_rotation.axis]

    if cubic_rotation.parity:
      repr_mat = self.gamma.four * repr_mat

    if double_element:
      repr_mat = -repr_mat

    return repr_mat
  
  def _setup_rotations(self):
    substitutions = [(I,eye(4)), (g1,self.gamma.one), (g2,self.gamma.two), (g3,self.gamma.three)]

    self.rotations = dict()
    rotation_map = _rotation_map
    for angle, axes_map in rotation_map.items():
      self.rotations[angle] = dict()
      for axis, repr_map in axes_map.items():
        self.rotations[angle][axis] = Matrix(repr_map.subs(substitutions).doit())

    self.conj_rotations = dict()
    conj_rotation_map = _conjugate_rotation_map
    for angle, axes_map in conj_rotation_map.items():
      self.conj_rotations[angle] = dict()
      for axis, repr_map in axes_map.items():
        self.conj_rotations[angle][axis] = Matrix(repr_map.subs(substitutions).doit())


spinor_representation = SpinorRepresentation()
