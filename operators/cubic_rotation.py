from enum import Enum, auto, unique
from math import gcd
from functools import reduce

from sympy import sqrt, cos, sin, pi
from sympy import eye, Identity, MatrixSymbol, Array, Matrix, Rational
from sympy import Quaternion, S, I
from sympy import nan

from .gamma import GammaRep, Gamma


class SpinorRotation:
  """
  TODO: 
    - Maybe this should just extend Quaternion?
    - Better solution to _POINT_GROUP_NAMES?
  """

  def __init__(self, axis, angle, parity):
    """
    Args:
      axis (3-tuple): the rotation axis that defines the rotation
      angle (float): the rotation angle in [0, 4pi)
      parity (+1/-1): specifies if there is a parity
    """
    self._quaternion = Quaternion.from_axis_angle(axis, angle)
    self._parity = parity

  @property
  def axis(self):
    return self._quaternion.to_axis_angle()[0]

  @property
  def angle(self):
    return self._quaternion.to_axis_angle()[1]

  @property
  def parity(self):
    return self._parity

  @property
  def quaternion(self):
    return self._quaternion

  def __mul__(self, other):
    if isinstance(other, self.__class__):
      bare_spinor_rotation = SpinorRotation.__new__(SpinorRotation)
      bare_spinor_rotation._quaternion = self.quaternion*other.quaternion
      bare_spinor_rotation._parity = self.parity*other.parity
      return bare_spinor_rotation
    elif isinstance(other, Momentum):
      rotated = [self.parity*pi for pi in Quaternion.rotate_point((other.x, other.y, other.z), self.quaternion)]
      return Momentum(*rotated)
    return NotImplemented

  def __invert__(self):
    conj_quaternion = Quaternion(self.quaternion.a, -self.quaternion.b, -self.quaternion.c, -self.quaternion.d)
    bare_spinor_rotation = SpinorRotation.__new__(SpinorRotation)
    bare_spinor_rotation._quaternion = conj_quaternion
    bare_spinor_rotation._parity = self.parity
    return bare_spinor_rotation

  def __str__(self):
    if self in _POINT_GROUP_NAMES:
      return _POINT_GROUP_NAMES[self]

    parity_sign = "{0:+}".format(self.parity)[0]
    return f"{parity_sign}({self.quaternion})"
      
  def __repr__(self):
    return self.__str__()

  def __hash__(self):
    return hash((self.parity, self.quaternion))

  def __eq__(self, other):
    if isinstance(other, self.__class__):
      return (self.parity, self.quaternion) == (other.parity, other.quaternion)
    return NotImplemented

  def __ne__(self, other):
    return not self.__eq__(other)


point_group = set()
_POINT_GROUP_NAMES = dict()
E     = SpinorRotation(( 0, 0, 1),      0, +1,); point_group.add(E);     _POINT_GROUP_NAMES[E] = "E"
C2x   = SpinorRotation(( 1, 0, 0),     pi, +1,); point_group.add(C2x);   _POINT_GROUP_NAMES[C2x] = "C2x"
C2y   = SpinorRotation(( 0, 1, 0),     pi, +1,); point_group.add(C2y);   _POINT_GROUP_NAMES[C2y] = "C2y"
C2z   = SpinorRotation(( 0, 0, 1),     pi, +1,); point_group.add(C2z);   _POINT_GROUP_NAMES[C2z] = "C2z"
C2a   = SpinorRotation(( 1, 1, 0),     pi, +1,); point_group.add(C2a);   _POINT_GROUP_NAMES[C2a] = "C2a"
C2b   = SpinorRotation(( 1,-1, 0),     pi, +1,); point_group.add(C2b);   _POINT_GROUP_NAMES[C2b] = "C2b"
C2c   = SpinorRotation(( 1, 0, 1),     pi, +1,); point_group.add(C2c);   _POINT_GROUP_NAMES[C2c] = "C2c"
C2d   = SpinorRotation((-1, 0, 1),     pi, +1,); point_group.add(C2d);   _POINT_GROUP_NAMES[C2d] = "C2d"
C2e   = SpinorRotation(( 0, 1, 1),     pi, +1,); point_group.add(C2e);   _POINT_GROUP_NAMES[C2e] = "C2e"
C2f   = SpinorRotation(( 0, 1,-1),     pi, +1,); point_group.add(C2f);   _POINT_GROUP_NAMES[C2f] = "C2f"
C3A   = SpinorRotation((-1,-1, 1), 2*pi/3, +1,); point_group.add(C3A);   _POINT_GROUP_NAMES[C3A] = "C3A"
C3B   = SpinorRotation((-1, 1,-1), 2*pi/3, +1,); point_group.add(C3B);   _POINT_GROUP_NAMES[C3B] = "C3B"
C3C   = SpinorRotation(( 1,-1,-1), 2*pi/3, +1,); point_group.add(C3C);   _POINT_GROUP_NAMES[C3C] = "C3C"
C3D   = SpinorRotation(( 1, 1, 1), 2*pi/3, +1,); point_group.add(C3D);   _POINT_GROUP_NAMES[C3D] = "C3D"
C3Ai  = SpinorRotation(( 1, 1,-1), 2*pi/3, +1,); point_group.add(C3Ai);  _POINT_GROUP_NAMES[C3Ai] = "C3Ai"
C3Bi  = SpinorRotation(( 1,-1, 1), 2*pi/3, +1,); point_group.add(C3Bi);  _POINT_GROUP_NAMES[C3Bi] = "C3Bi"
C3Ci  = SpinorRotation((-1, 1, 1), 2*pi/3, +1,); point_group.add(C3Ci);  _POINT_GROUP_NAMES[C3Ci] = "C3Ci"
C3Di  = SpinorRotation((-1,-1,-1), 2*pi/3, +1,); point_group.add(C3Di);  _POINT_GROUP_NAMES[C3Di] = "C3Di"
C4x   = SpinorRotation(( 1, 0, 0),   pi/2, +1,); point_group.add(C4x);   _POINT_GROUP_NAMES[C4x] = "C4x"
C4y   = SpinorRotation(( 0, 1, 0),   pi/2, +1,); point_group.add(C4y);   _POINT_GROUP_NAMES[C4y] = "C4y"
C4z   = SpinorRotation(( 0, 0, 1),   pi/2, +1,); point_group.add(C4z);   _POINT_GROUP_NAMES[C4z] = "C4z"
C4xi  = SpinorRotation((-1, 0, 0),   pi/2, +1,); point_group.add(C4xi);  _POINT_GROUP_NAMES[C4xi] = "C4xi"
C4yi  = SpinorRotation(( 0,-1, 0),   pi/2, +1,); point_group.add(C4yi);  _POINT_GROUP_NAMES[C4yi] = "C4yi"
C4zi  = SpinorRotation(( 0, 0,-1),   pi/2, +1,); point_group.add(C4zi);  _POINT_GROUP_NAMES[C4zi] = "C4zi"
Is    = SpinorRotation(( 0, 0, 1),      0, -1,); point_group.add(Is);    _POINT_GROUP_NAMES[Is] = "Is"
I2x   = SpinorRotation(( 1, 0, 0),     pi, -1,); point_group.add(I2x);   _POINT_GROUP_NAMES[I2x] = "I2x"
I2y   = SpinorRotation(( 0, 1, 0),     pi, -1,); point_group.add(I2y);   _POINT_GROUP_NAMES[I2y] = "I2y"
I2z   = SpinorRotation(( 0, 0, 1),     pi, -1,); point_group.add(I2z);   _POINT_GROUP_NAMES[I2z] = "I2z"
I2a   = SpinorRotation(( 1, 1, 0),     pi, -1,); point_group.add(I2a);   _POINT_GROUP_NAMES[I2a] = "I2a"
I2b   = SpinorRotation(( 1,-1, 0),     pi, -1,); point_group.add(I2b);   _POINT_GROUP_NAMES[I2b] = "I2b"
I2c   = SpinorRotation(( 1, 0, 1),     pi, -1,); point_group.add(I2c);   _POINT_GROUP_NAMES[I2c] = "I2c"
I2d   = SpinorRotation((-1, 0, 1),     pi, -1,); point_group.add(I2d);   _POINT_GROUP_NAMES[I2d] = "I2d"
I2e   = SpinorRotation(( 0, 1, 1),     pi, -1,); point_group.add(I2e);   _POINT_GROUP_NAMES[I2e] = "I2e"
I2f   = SpinorRotation(( 0, 1,-1),     pi, -1,); point_group.add(I2f);   _POINT_GROUP_NAMES[I2f] = "I2f"
I3A   = SpinorRotation((-1,-1, 1), 2*pi/3, -1,); point_group.add(I3A);   _POINT_GROUP_NAMES[I3A] = "I3A"
I3B   = SpinorRotation((-1, 1,-1), 2*pi/3, -1,); point_group.add(I3B);   _POINT_GROUP_NAMES[I3B] = "I3B"
I3C   = SpinorRotation(( 1,-1,-1), 2*pi/3, -1,); point_group.add(I3C);   _POINT_GROUP_NAMES[I3C] = "I3C"
I3D   = SpinorRotation(( 1, 1, 1), 2*pi/3, -1,); point_group.add(I3D);   _POINT_GROUP_NAMES[I3D] = "I3D"
I3Ai  = SpinorRotation(( 1, 1,-1), 2*pi/3, -1,); point_group.add(I3Ai);  _POINT_GROUP_NAMES[I3Ai] = "I3Ai"
I3Bi  = SpinorRotation(( 1,-1, 1), 2*pi/3, -1,); point_group.add(I3Bi);  _POINT_GROUP_NAMES[I3Bi] = "I3Bi"
I3Ci  = SpinorRotation((-1, 1, 1), 2*pi/3, -1,); point_group.add(I3Ci);  _POINT_GROUP_NAMES[I3Ci] = "I3Ci"
I3Di  = SpinorRotation((-1,-1,-1), 2*pi/3, -1,); point_group.add(I3Di);  _POINT_GROUP_NAMES[I3Di] = "I3Di"
I4x   = SpinorRotation(( 1, 0, 0),   pi/2, -1,); point_group.add(I4x);   _POINT_GROUP_NAMES[I4x] = "I4x"
I4y   = SpinorRotation(( 0, 1, 0),   pi/2, -1,); point_group.add(I4y);   _POINT_GROUP_NAMES[I4y] = "I4y"
I4z   = SpinorRotation(( 0, 0, 1),   pi/2, -1,); point_group.add(I4z);   _POINT_GROUP_NAMES[I4z] = "I4z"
I4xi  = SpinorRotation((-1, 0, 0),   pi/2, -1,); point_group.add(I4xi);  _POINT_GROUP_NAMES[I4xi] = "I4xi"
I4yi  = SpinorRotation(( 0,-1, 0),   pi/2, -1,); point_group.add(I4yi);  _POINT_GROUP_NAMES[I4yi] = "I4yi"
I4zi  = SpinorRotation(( 0, 0,-1),   pi/2, -1,); point_group.add(I4zi);  _POINT_GROUP_NAMES[I4zi] = "I4zi"

EE    = SpinorRotation(( 0, 0, 1),   2*pi, +1,); point_group.add(EE);    _POINT_GROUP_NAMES[EE] = "EE"
EC2x  = SpinorRotation(( 1, 0, 0),   3*pi, +1,); point_group.add(EC2x);  _POINT_GROUP_NAMES[EC2x] = "EC2x"
EC2y  = SpinorRotation(( 0, 1, 0),   3*pi, +1,); point_group.add(EC2y);  _POINT_GROUP_NAMES[EC2y] = "EC2y"
EC2z  = SpinorRotation(( 0, 0, 1),   3*pi, +1,); point_group.add(EC2z);  _POINT_GROUP_NAMES[EC2z] = "EC2z"
EC2a  = SpinorRotation(( 1, 1, 0),   3*pi, +1,); point_group.add(EC2a);  _POINT_GROUP_NAMES[EC2a] = "EC2a"
EC2b  = SpinorRotation(( 1,-1, 0),   3*pi, +1,); point_group.add(EC2b);  _POINT_GROUP_NAMES[EC2b] = "EC2b"
EC2c  = SpinorRotation(( 1, 0, 1),   3*pi, +1,); point_group.add(EC2c);  _POINT_GROUP_NAMES[EC2c] = "EC2c"
EC2d  = SpinorRotation((-1, 0, 1),   3*pi, +1,); point_group.add(EC2d);  _POINT_GROUP_NAMES[EC2d] = "EC2d"
EC2e  = SpinorRotation(( 0, 1, 1),   3*pi, +1,); point_group.add(EC2e);  _POINT_GROUP_NAMES[EC2e] = "EC2e"
EC2f  = SpinorRotation(( 0, 1,-1),   3*pi, +1,); point_group.add(EC2f);  _POINT_GROUP_NAMES[EC2f] = "EC2f"
EC3A  = SpinorRotation((-1,-1, 1), 8*pi/3, +1,); point_group.add(EC3A);  _POINT_GROUP_NAMES[EC3A] = "EC3A"
EC3B  = SpinorRotation((-1, 1,-1), 8*pi/3, +1,); point_group.add(EC3B);  _POINT_GROUP_NAMES[EC3B] = "EC3B"
EC3C  = SpinorRotation(( 1,-1,-1), 8*pi/3, +1,); point_group.add(EC3C);  _POINT_GROUP_NAMES[EC3C] = "EC3C"
EC3D  = SpinorRotation(( 1, 1, 1), 8*pi/3, +1,); point_group.add(EC3D);  _POINT_GROUP_NAMES[EC3D] = "EC3D"
EC3Ai = SpinorRotation(( 1, 1,-1), 8*pi/3, +1,); point_group.add(EC3Ai); _POINT_GROUP_NAMES[EC3Ai] = "EC3Ai"
EC3Bi = SpinorRotation(( 1,-1, 1), 8*pi/3, +1,); point_group.add(EC3Bi); _POINT_GROUP_NAMES[EC3Bi] = "EC3Bi"
EC3Ci = SpinorRotation((-1, 1, 1), 8*pi/3, +1,); point_group.add(EC3Ci); _POINT_GROUP_NAMES[EC3Ci] = "EC3Ci"
EC3Di = SpinorRotation((-1,-1,-1), 8*pi/3, +1,); point_group.add(EC3Di); _POINT_GROUP_NAMES[EC3Di] = "EC3Di"
EC4x  = SpinorRotation(( 1, 0, 0), 5*pi/2, +1,); point_group.add(EC4x);  _POINT_GROUP_NAMES[EC4x] = "EC4x"
EC4y  = SpinorRotation(( 0, 1, 0), 5*pi/2, +1,); point_group.add(EC4y);  _POINT_GROUP_NAMES[EC4y] = "EC4y"
EC4z  = SpinorRotation(( 0, 0, 1), 5*pi/2, +1,); point_group.add(EC4z);  _POINT_GROUP_NAMES[EC4z] = "EC4z"
EC4xi = SpinorRotation((-1, 0, 0), 5*pi/2, +1,); point_group.add(EC4xi); _POINT_GROUP_NAMES[EC4xi] = "EC4xi"
EC4yi = SpinorRotation(( 0,-1, 0), 5*pi/2, +1,); point_group.add(EC4yi); _POINT_GROUP_NAMES[EC4yi] = "EC4yi"
EC4zi = SpinorRotation(( 0, 0,-1), 5*pi/2, +1,); point_group.add(EC4zi); _POINT_GROUP_NAMES[EC4zi] = "EC4zi"
EIs   = SpinorRotation(( 0, 0, 1),   2*pi, -1,); point_group.add(EIs);   _POINT_GROUP_NAMES[EIs] = "EIs"
EI2x  = SpinorRotation(( 1, 0, 0),   3*pi, -1,); point_group.add(EI2x);  _POINT_GROUP_NAMES[EI2x] = "EI2x"
EI2y  = SpinorRotation(( 0, 1, 0),   3*pi, -1,); point_group.add(EI2y);  _POINT_GROUP_NAMES[EI2y] = "EI2y"
EI2z  = SpinorRotation(( 0, 0, 1),   3*pi, -1,); point_group.add(EI2z);  _POINT_GROUP_NAMES[EI2z] = "EI2z"
EI2a  = SpinorRotation(( 1, 1, 0),   3*pi, -1,); point_group.add(EI2a);  _POINT_GROUP_NAMES[EI2a] = "EI2a"
EI2b  = SpinorRotation(( 1,-1, 0),   3*pi, -1,); point_group.add(EI2b);  _POINT_GROUP_NAMES[EI2b] = "EI2b"
EI2c  = SpinorRotation(( 1, 0, 1),   3*pi, -1,); point_group.add(EI2c);  _POINT_GROUP_NAMES[EI2c] = "EI2c"
EI2d  = SpinorRotation((-1, 0, 1),   3*pi, -1,); point_group.add(EI2d);  _POINT_GROUP_NAMES[EI2d] = "EI2d"
EI2e  = SpinorRotation(( 0, 1, 1),   3*pi, -1,); point_group.add(EI2e);  _POINT_GROUP_NAMES[EI2e] = "EI2e"
EI2f  = SpinorRotation(( 0, 1,-1),   3*pi, -1,); point_group.add(EI2f);  _POINT_GROUP_NAMES[EI2f] = "EI2f"
EI3A  = SpinorRotation((-1,-1, 1), 8*pi/3, -1,); point_group.add(EI3A);  _POINT_GROUP_NAMES[EI3A] = "EI3A"
EI3B  = SpinorRotation((-1, 1,-1), 8*pi/3, -1,); point_group.add(EI3B);  _POINT_GROUP_NAMES[EI3B] = "EI3B"
EI3C  = SpinorRotation(( 1,-1,-1), 8*pi/3, -1,); point_group.add(EI3C);  _POINT_GROUP_NAMES[EI3C] = "EI3C"
EI3D  = SpinorRotation(( 1, 1, 1), 8*pi/3, -1,); point_group.add(EI3D);  _POINT_GROUP_NAMES[EI3D] = "EI3D"
EI3Ai = SpinorRotation(( 1, 1,-1), 8*pi/3, -1,); point_group.add(EI3Ai); _POINT_GROUP_NAMES[EI3Ai] = "EI3Ai"
EI3Bi = SpinorRotation(( 1,-1, 1), 8*pi/3, -1,); point_group.add(EI3Bi); _POINT_GROUP_NAMES[EI3Bi] = "EI3Bi"
EI3Ci = SpinorRotation((-1, 1, 1), 8*pi/3, -1,); point_group.add(EI3Ci); _POINT_GROUP_NAMES[EI3Ci] = "EI3Ci"
EI3Di = SpinorRotation((-1,-1,-1), 8*pi/3, -1,); point_group.add(EI3Di); _POINT_GROUP_NAMES[EI3Di] = "EI3Di"
EI4x  = SpinorRotation(( 1, 0, 0), 5*pi/2, -1,); point_group.add(EI4x);  _POINT_GROUP_NAMES[EI4x] = "EI4x"
EI4y  = SpinorRotation(( 0, 1, 0), 5*pi/2, -1,); point_group.add(EI4y);  _POINT_GROUP_NAMES[EI4y] = "EI4y"
EI4z  = SpinorRotation(( 0, 0, 1), 5*pi/2, -1,); point_group.add(EI4z);  _POINT_GROUP_NAMES[EI4z] = "EI4z"
EI4xi = SpinorRotation((-1, 0, 0), 5*pi/2, -1,); point_group.add(EI4xi); _POINT_GROUP_NAMES[EI4xi] = "EI4xi"
EI4yi = SpinorRotation(( 0,-1, 0), 5*pi/2, -1,); point_group.add(EI4yi); _POINT_GROUP_NAMES[EI4yi] = "EI4yi"
EI4zi = SpinorRotation(( 0, 0,-1), 5*pi/2, -1,); point_group.add(EI4zi); _POINT_GROUP_NAMES[EI4zi] = "EI4zi"



class Momentum(Array):

  mom_map = {
       0: '0',
       1: '+',
      -1: '-',
       2: '#',
      -2: '=',
       3: 'T',
      -3: 't',
       4: 'V',
      -4: 'v',
       5: 'F',
      -5: 'f',
       6: 'S',
      -6: 's',
       8: 'A',
      -8: 'a',
       9: 'N',
      -9: 'n',
      10: 'Z',
     -10: 'z',
      11: 'E',
     -11: 'e',
  }

  mom_map_inv = {v: k for k, v in mom_map.items()}

  def __init__(self, *momentum):
    if not isinstance(momentum[0], Momentum):
      self._p_x = momentum[0]
      self._p_y = momentum[1]
      self._p_z = momentum[2]

  def __new__(cls, *momentum):
    if isinstance(momentum[0], Momentum):
      return momentum[0]

    return super().__new__(cls, momentum)

  @classmethod
  def createFromMomRay(cls, mom_ray):
    return cls(cls.mom_map_inv[mom_ray[0]], cls.mom_map_inv[mom_ray[1]], cls.mom_map_inv[mom_ray[2]])

  @property
  def mom_ray(self):
    return f"{self.mom_map[self.x]}{self.mom_map[self.y]}{self.mom_map[self.z]}"

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
  def sq(self):
    return self._p_x**2 + self._p_y**2 + self._p_z**2

  @property
  def pref(self):
    return Momentum(*sorted([abs(pi) for pi in self]))

  @property
  def reduced_pref(self):
    if self.sq == 0:
      return self

    factor = reduce(lambda x,y: gcd(x,y), self.pref)
    return Momentum(self.pref.x//factor, self.pref.y//factor, self.pref.z//factor)

  @property
  def reduced(self):
    if self.sq == 0:
      return self

    factor = reduce(lambda x,y: gcd(x,y), self)
    return Momentum(self.x//factor, self.y//factor, self.z//factor)

  '''
  def __getitem__(self, i):
    if i in range(1,4):
      return super().__getitem__(i-1)
    else:
      raise IndexError("Index out of range")
  '''
      

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

  def __lt__(self, other):
    if isinstance(other, self.__class__):
      return self.__repr__() < other.__repr__()
    return NotImplemented

  def __le__(self, other):
    if isinstance(other, self.__class__):
      return self.__repr__() <= other.__repr__()
    return NotImplemented

  def __gt__(self, other):
    if isinstance(other, self.__class__):
      return self.__repr__() > other.__repr__()
    return NotImplemented

  def __ge__(self, other):
    if isinstance(other, self.__class__):
      return self.__repr__() >= other.__repr__()
    return NotImplemented

  def __neg__(self):
    return Momentum(-self.x, -self.y, -self.z)

  def __add__(self, other):
    if isinstance(other, self.__class__):
      return Momentum(self.x + other.x, self.y + other.y, self.z + other.z)
    return NotImplemented
  
  def __sub__(self, other):
    return self.__add__(-other)

  def __mul__(self, other):
    if isinstance(other, self.__class__):
      px = self.y*other.z - self.z*other.y
      py = self.z*other.x - self.x*other.z
      pz = self.x*other.y - self.y*other.x
      return Momentum(px, py, pz)

    else:
      return Momentum(self.x * other, self.y * other, self.z * other)

  def __rmul__(self, other):
    return Momentum(other * self.x, other * self.y, other * self.z)

  def __truediv__(self, other):
    if isinstance(other, self.__class__):
      scales = set()
      if other.x != 0 and self.x != 0:
        scales.add(self.x/other.x)
      elif self.x != other.x:
        return None

      if other.y != 0 and self.y != 0:
        scales.add(self.y/other.y)
      elif self.y != other.y:
        return None

      if other.z != 0 and self.z != 0:
        scales.add(self.z/other.z)
      elif self.z != other.z:
        return None

      if len(scales) > 1:
        return None
      elif len(scales) == 0:
        return 1
      else:
        return scales.pop()

    return Momentum(self.x // other, self.y // other, self.z // other)

  def __rtruediv__(self, other):
    return Momentum(self.x // other, self.y // other, self.z // other)


P = Momentum
P0 = P(0,0,0)

_REFERENCE_ROTATIONS = {
    P( 0, 0, 0): E,

    P( 0, 0, 1): E,
    P( 0, 0,-1): C2x,
    P( 1, 0, 0): C4y,
    P(-1, 0, 0): C4yi,
    P( 0,-1, 0): C4x,
    P( 0, 1, 0): C4xi,

    P( 0, 1, 1): E,
    P( 0,-1,-1): C2x,
    P( 0, 1,-1): C4xi,
    P( 0,-1, 1): C4x,
    P( 1, 0, 1): C4zi,
    P(-1, 0,-1): C2b,
    P( 1, 0,-1): C2a,
    P(-1, 0, 1): C4z,
    P( 1, 1, 0): C4y,
    P(-1,-1, 0): C2d,
    P( 1,-1, 0): C2c,
    P(-1, 1, 0): C4yi,

    P( 1, 1, 1): E,
    P( 1, 1,-1): C4y,
    P( 1,-1, 1): C4x,
    P( 1,-1,-1): C2x,
    P(-1, 1, 1): C4z,
    P(-1, 1,-1): C2y,
    P(-1,-1, 1): C2z,
    P(-1,-1,-1): C2d,

    P( 0, 1, 2): E,
    P( 0, 2, 1): C2e,
    P( 1, 0, 2): C4zi,
    P( 2, 0, 1): C3D,
    P( 1, 2, 0): C3Di,
    P( 2, 1, 0): C4y,
    P( 0,-1, 2): C2z,
    P( 0, 2,-1): C4xi,
    P(-1, 0, 2): C4z,
    P( 2, 0,-1): C3B,
    P(-1, 2, 0): C3Ci,
    P( 2,-1, 0): C2c,
    P( 0, 1,-2): C2y,
    P( 0,-2, 1): C4x,
    P( 1, 0,-2): C2a,
    P(-2, 0, 1): C3C,
    P( 1,-2, 0): C3Ai,
    P(-2, 1, 0): C4yi,
    P( 0,-1,-2): C2x,
    P( 0,-2,-1): C2f,
    P(-1, 0,-2): C2b,
    P(-2, 0,-1): C3A,
    P(-1,-2, 0): C3Bi,
    P(-2,-1, 0): C2d,

    P( 1, 1, 2): E,
    P(-1, 1, 2): C4z,
    P( 1,-1, 2): C4zi,
    P(-1,-1, 2): C2z,
    P( 1, 1,-2): C2a,
    P(-1, 1,-2): C2y,
    P( 1,-1,-2): C2x,
    P(-1,-1,-2): C2b,
    P( 1, 2, 1): C3Di,
    P(-1, 2, 1): C2e,
    P( 1, 2,-1): C4xi,
    P(-1, 2,-1): C3Ci,
    P( 1,-2, 1): C4x,
    P(-1,-2, 1): C3Bi,
    P( 1,-2,-1): C3Ai,
    P(-1,-2,-1): C2f,
    P( 2, 1, 1): C3D,
    P( 2,-1, 1): C2c,
    P( 2, 1,-1): C4y,
    P( 2,-1,-1): C3B,
    P(-2, 1, 1): C4yi,
    P(-2,-1, 1): C3C,
    P(-2, 1,-1): C3A,
    P(-2,-1,-1): C2d,
}

_BOSONIC_LITTLE_GROUP_IRREPS = {
    "Oh":    ["A1g", "A2g", "Eg", "T1g", "T2g", "A1u", "A2u", "Eu", "T1u", "T2u"],
    "C4v":   ["A1", "A2", "B1", "B2", "E"],
    "C2v":   ["A1", "A2", "B1", "B2"],
    "C3v":   ["A1", "A2", "E"],
    "Cs012": ["A1", "A2"],
    "Cs112": ["A1", "A2"]
}

_FERMIONIC_LITTLE_GROUP_IRREPS = {
    "Oh":    ["G1g", "G2g", "Hg", "G1u", "G2u", "Hu"],
    "C4v":   ["G1", "G2"],
    "C2v":   ["G"],
    "C3v":   ["F1", "F2", "G"],
    "Cs012": ["F1", "F2"],
    "Cs112": ["F1", "F2"]
}

_BOSONIC_LITTLE_GROUP_IRREPS_ROWS = {
    ("Oh", "A1g"): 1, 
    ("Oh", "A2g"): 1,
    ("Oh", "Eg"):  2,
    ("Oh", "T1g"): 3,
    ("Oh", "T2g"): 3,
    ("Oh", "A1u"): 1,
    ("Oh", "A2u"): 1,
    ("Oh", "Eu"):  2,
    ("Oh", "T1u"): 3,
    ("Oh", "T2u"): 3,

    ("C4v", "A1"): 1,
    ("C4v", "A2"): 1,
    ("C4v", "B1"): 1,
    ("C4v", "B2"): 1,
    ("C4v", "E"):  2,

    ("C2v", "A1"): 1,
    ("C2v", "A2"): 1,
    ("C2v", "B1"): 1,
    ("C2v", "B2"): 1,

    ("C3v", "A1"): 1,
    ("C3v", "A2"): 1,
    ("C3v", "E"):  2,

    ("Cs012", "A1"): 1,
    ("Cs012", "A2"): 1,

    ("Cs112", "A1"): 1,
    ("Cs112", "A2"): 1,
}

_FERMIONIC_LITTLE_GROUP_IRREPS_ROWS = {
    ("Oh", "G1g"): 2,
    ("Oh", "G2g"): 2,
    ("Oh", "Hg"):  4,
    ("Oh", "G1u"): 2,
    ("Oh", "G2u"): 2,
    ("Oh", "Hu"):  4,

    ("C4v", "G1"): 2,
    ("C4v", "G2"): 2,

    ("C2v", "G"):  2,

    ("C3v", "F1"): 1,
    ("C3v", "F2"): 1,
    ("C3v", "G"):  2,

    ("Cs012", "F1"): 1,
    ("Cs012", "F2"): 1,

    ("Cs112", "F1"): 1,
    ("Cs112", "F2"): 1,
}

# Conjugacy Classes
Oh_1  = frozenset([E])
Oh_2  = frozenset([C3A, C3B, C3C, C3D, C3Ai, C3Bi, C3Ci, C3Di])
Oh_3  = frozenset([C2x, C2y, C2z, EC2x, EC2y, EC2z])
Oh_4  = frozenset([C4x, C4y, C4z, C4xi, C4yi, C4zi])
Oh_5  = frozenset([C2a, C2b, C2c, C2d, C2e, C2f, EC2a, EC2b, EC2c, EC2d, EC2e, EC2f])
Oh_6  = frozenset([EE])
Oh_7  = frozenset([EC3A, EC3B, EC3C, EC3D, EC3Ai, EC3Bi, EC3Ci, EC3Di])
Oh_8  = frozenset([EC4x, EC4y, EC4z, EC4xi, EC4yi, EC4zi])
Oh_9  = frozenset([Is])
Oh_10 = frozenset([I3A, I3B, I3C, I3D, I3Ai, I3Bi, I3Ci, I3Di])
Oh_11 = frozenset([I2x, I2y, I2z, EI2x, EI2y, EI2z])
Oh_12 = frozenset([I4x, I4y, I4z, I4xi, I4yi, I4zi])
Oh_13 = frozenset([I2a, I2b, I2c, I2d, I2e, I2f, EI2a, EI2b, EI2c, EI2d, EI2e, EI2f])
Oh_14 = frozenset([EIs])
Oh_15 = frozenset([EI3A, EI3B, EI3C, EI3D, EI3Ai, EI3Bi, EI3Ci, EI3Di])
Oh_16 = frozenset([EI4x, EI4y, EI4z, EI4xi, EI4yi, EI4zi])

C4v_1 = frozenset([E])
C4v_2 = frozenset([C2z, EC2z])
C4v_3 = frozenset([C4z, C4zi])
C4v_4 = frozenset([I2x, I2y, EI2x, EI2y])
C4v_5 = frozenset([I2a, I2b, EI2a, EI2b])
C4v_6 = frozenset([EE])
C4v_7 = frozenset([EC4z, EC4zi])

C2v_1 = frozenset([E])
C2v_2 = frozenset([C2e, EC2e])
C2v_3 = frozenset([I2f, EI2f])
C2v_4 = frozenset([I2x, EI2x])
C2v_5 = frozenset([EE])

C3v_1 = frozenset([E])
C3v_2 = frozenset([C3D, C3Di])
C3v_3 = frozenset([I2b, I2d, I2f])
C3v_4 = frozenset([EE])
C3v_5 = frozenset([EC3D, EC3Di])
C3v_6 = frozenset([EI2b, EI2d, EI2f])

Cs012_1 = frozenset([E])
Cs012_2 = frozenset([I2x])
Cs012_3 = frozenset([EI2x])
Cs012_4 = frozenset([EE])

Cs112_1 = frozenset([E])
Cs112_2 = frozenset([I2b])
Cs112_3 = frozenset([EI2b])
Cs112_4 = frozenset([EE])

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
    ("Oh",  "A1g", Oh_11): 1,
    ("Oh",  "A1g", Oh_12): 1,
    ("Oh",  "A1g", Oh_13): 1,
    ("Oh",  "A1g", Oh_14): 1,
    ("Oh",  "A1g", Oh_15): 1,
    ("Oh",  "A1g", Oh_16): 1,

    ("Oh",  "A2g", Oh_1):  1,
    ("Oh",  "A2g", Oh_2):  1,
    ("Oh",  "A2g", Oh_3):  1,
    ("Oh",  "A2g", Oh_4):  -1,
    ("Oh",  "A2g", Oh_5):  -1,
    ("Oh",  "A2g", Oh_6):  1,
    ("Oh",  "A2g", Oh_7):  1,
    ("Oh",  "A2g", Oh_8):  -1,
    ("Oh",  "A2g", Oh_9):  1,
    ("Oh",  "A2g", Oh_10): 1,
    ("Oh",  "A2g", Oh_11): 1,
    ("Oh",  "A2g", Oh_12): -1,
    ("Oh",  "A2g", Oh_13): -1,
    ("Oh",  "A2g", Oh_14): 1,
    ("Oh",  "A2g", Oh_15): 1,
    ("Oh",  "A2g", Oh_16): -1,

    ("Oh",  "Eg",  Oh_1):  2,
    ("Oh",  "Eg",  Oh_2):  -1,
    ("Oh",  "Eg",  Oh_3):  2,
    ("Oh",  "Eg",  Oh_4):  0,
    ("Oh",  "Eg",  Oh_5):  0,
    ("Oh",  "Eg",  Oh_6):  2,
    ("Oh",  "Eg",  Oh_7):  -1,
    ("Oh",  "Eg",  Oh_8):  0,
    ("Oh",  "Eg",  Oh_9):  2,
    ("Oh",  "Eg",  Oh_10): -1,
    ("Oh",  "Eg",  Oh_11): 2,
    ("Oh",  "Eg",  Oh_12): 0,
    ("Oh",  "Eg",  Oh_13): 0,
    ("Oh",  "Eg",  Oh_14): 2,
    ("Oh",  "Eg",  Oh_15): -1,
    ("Oh",  "Eg",  Oh_16): 0,

    ("Oh",  "T1g", Oh_1):  3,
    ("Oh",  "T1g", Oh_2):  0,
    ("Oh",  "T1g", Oh_3):  -1,
    ("Oh",  "T1g", Oh_4):  1,
    ("Oh",  "T1g", Oh_5):  -1,
    ("Oh",  "T1g", Oh_6):  3,
    ("Oh",  "T1g", Oh_7):  0,
    ("Oh",  "T1g", Oh_8):  1,
    ("Oh",  "T1g", Oh_9):  3,
    ("Oh",  "T1g", Oh_10): 0,
    ("Oh",  "T1g", Oh_11): -1,
    ("Oh",  "T1g", Oh_12): 1,
    ("Oh",  "T1g", Oh_13): -1,
    ("Oh",  "T1g", Oh_14): 3,
    ("Oh",  "T1g", Oh_15): 0,
    ("Oh",  "T1g", Oh_16): 1,

    ("Oh",  "T2g", Oh_1):  3,
    ("Oh",  "T2g", Oh_2):  0,
    ("Oh",  "T2g", Oh_3):  -1,
    ("Oh",  "T2g", Oh_4):  -1,
    ("Oh",  "T2g", Oh_5):  1,
    ("Oh",  "T2g", Oh_6):  3,
    ("Oh",  "T2g", Oh_7):  0,
    ("Oh",  "T2g", Oh_8):  -1,
    ("Oh",  "T2g", Oh_9):  3,
    ("Oh",  "T2g", Oh_10): 0,
    ("Oh",  "T2g", Oh_11): -1,
    ("Oh",  "T2g", Oh_12): -1,
    ("Oh",  "T2g", Oh_13): 1,
    ("Oh",  "T2g", Oh_14): 3,
    ("Oh",  "T2g", Oh_15): 0,
    ("Oh",  "T2g", Oh_16): -1,

    ("Oh",  "G1g", Oh_1):  2,
    ("Oh",  "G1g", Oh_2):  1,
    ("Oh",  "G1g", Oh_3):  0,
    ("Oh",  "G1g", Oh_4):  sqrt(2),
    ("Oh",  "G1g", Oh_5):  0,
    ("Oh",  "G1g", Oh_6):  -2,
    ("Oh",  "G1g", Oh_7):  -1,
    ("Oh",  "G1g", Oh_8):  -sqrt(2),
    ("Oh",  "G1g", Oh_9):  2,
    ("Oh",  "G1g", Oh_10): 1,
    ("Oh",  "G1g", Oh_11): 0,
    ("Oh",  "G1g", Oh_12): sqrt(2),
    ("Oh",  "G1g", Oh_13): 0,
    ("Oh",  "G1g", Oh_14): -2,
    ("Oh",  "G1g", Oh_15): -1,
    ("Oh",  "G1g", Oh_16): -sqrt(2),

    ("Oh",  "G2g", Oh_1):  2,
    ("Oh",  "G2g", Oh_2):  1,
    ("Oh",  "G2g", Oh_3):  0,
    ("Oh",  "G2g", Oh_4):  -sqrt(2),
    ("Oh",  "G2g", Oh_5):  0,
    ("Oh",  "G2g", Oh_6):  -2,
    ("Oh",  "G2g", Oh_7):  -1,
    ("Oh",  "G2g", Oh_8):  sqrt(2),
    ("Oh",  "G2g", Oh_9):  2,
    ("Oh",  "G2g", Oh_10): 1,
    ("Oh",  "G2g", Oh_11): 0,
    ("Oh",  "G2g", Oh_12): -sqrt(2),
    ("Oh",  "G2g", Oh_13): 0,
    ("Oh",  "G2g", Oh_14): -2,
    ("Oh",  "G2g", Oh_15): -1,
    ("Oh",  "G2g", Oh_16): sqrt(2),

    ("Oh",  "Hg",  Oh_1):  4,
    ("Oh",  "Hg",  Oh_2):  -1,
    ("Oh",  "Hg",  Oh_3):  0,
    ("Oh",  "Hg",  Oh_4):  0,
    ("Oh",  "Hg",  Oh_5):  0,
    ("Oh",  "Hg",  Oh_6):  -4,
    ("Oh",  "Hg",  Oh_7):  1,
    ("Oh",  "Hg",  Oh_8):  0,
    ("Oh",  "Hg",  Oh_9):  4,
    ("Oh",  "Hg",  Oh_10): -1,
    ("Oh",  "Hg",  Oh_11): 0,
    ("Oh",  "Hg",  Oh_12): 0,
    ("Oh",  "Hg",  Oh_13): 0,
    ("Oh",  "Hg",  Oh_14): -4,
    ("Oh",  "Hg",  Oh_15): 1,
    ("Oh",  "Hg",  Oh_16): 0,

    ("Oh",  "A1u", Oh_1):  1,
    ("Oh",  "A1u", Oh_2):  1,
    ("Oh",  "A1u", Oh_3):  1,
    ("Oh",  "A1u", Oh_4):  1,
    ("Oh",  "A1u", Oh_5):  1,
    ("Oh",  "A1u", Oh_6):  1,
    ("Oh",  "A1u", Oh_7):  1,
    ("Oh",  "A1u", Oh_8):  1,
    ("Oh",  "A1u", Oh_9):  -1,
    ("Oh",  "A1u", Oh_10): -1,
    ("Oh",  "A1u", Oh_11): -1,
    ("Oh",  "A1u", Oh_12): -1,
    ("Oh",  "A1u", Oh_13): -1,
    ("Oh",  "A1u", Oh_14): -1,
    ("Oh",  "A1u", Oh_15): -1,
    ("Oh",  "A1u", Oh_16): -1,

    ("Oh",  "A2u", Oh_1):  1,
    ("Oh",  "A2u", Oh_2):  1,
    ("Oh",  "A2u", Oh_3):  1,
    ("Oh",  "A2u", Oh_4):  -1,
    ("Oh",  "A2u", Oh_5):  -1,
    ("Oh",  "A2u", Oh_6):  1,
    ("Oh",  "A2u", Oh_7):  1,
    ("Oh",  "A2u", Oh_8):  -1,
    ("Oh",  "A2u", Oh_9):  -1,
    ("Oh",  "A2u", Oh_10): -1,
    ("Oh",  "A2u", Oh_11): -1,
    ("Oh",  "A2u", Oh_12): 1,
    ("Oh",  "A2u", Oh_13): 1,
    ("Oh",  "A2u", Oh_14): -1,
    ("Oh",  "A2u", Oh_15): -1,
    ("Oh",  "A2u", Oh_16): 1,

    ("Oh",  "Eu",  Oh_1):  2,
    ("Oh",  "Eu",  Oh_2):  -1,
    ("Oh",  "Eu",  Oh_3):  2,
    ("Oh",  "Eu",  Oh_4):  0,
    ("Oh",  "Eu",  Oh_5):  0,
    ("Oh",  "Eu",  Oh_6):  2,
    ("Oh",  "Eu",  Oh_7):  -1,
    ("Oh",  "Eu",  Oh_8):  0,
    ("Oh",  "Eu",  Oh_9):  -2,
    ("Oh",  "Eu",  Oh_10): 1,
    ("Oh",  "Eu",  Oh_11): -2,
    ("Oh",  "Eu",  Oh_12): 0,
    ("Oh",  "Eu",  Oh_13): 0,
    ("Oh",  "Eu",  Oh_14): -2,
    ("Oh",  "Eu",  Oh_15): 1,
    ("Oh",  "Eu",  Oh_16): 0,

    ("Oh",  "T1u", Oh_1):  3,
    ("Oh",  "T1u", Oh_2):  0,
    ("Oh",  "T1u", Oh_3):  -1,
    ("Oh",  "T1u", Oh_4):  1,
    ("Oh",  "T1u", Oh_5):  -1,
    ("Oh",  "T1u", Oh_6):  3,
    ("Oh",  "T1u", Oh_7):  0,
    ("Oh",  "T1u", Oh_8):  1,
    ("Oh",  "T1u", Oh_9):  -3,
    ("Oh",  "T1u", Oh_10): 0,
    ("Oh",  "T1u", Oh_11): 1,
    ("Oh",  "T1u", Oh_12): -1,
    ("Oh",  "T1u", Oh_13): 1,
    ("Oh",  "T1u", Oh_14): -3,
    ("Oh",  "T1u", Oh_15): 0,
    ("Oh",  "T1u", Oh_16): -1,

    ("Oh",  "T2u", Oh_1):  3,
    ("Oh",  "T2u", Oh_2):  0,
    ("Oh",  "T2u", Oh_3):  -1,
    ("Oh",  "T2u", Oh_4):  -1,
    ("Oh",  "T2u", Oh_5):  1,
    ("Oh",  "T2u", Oh_6):  3,
    ("Oh",  "T2u", Oh_7):  0,
    ("Oh",  "T2u", Oh_8):  -1,
    ("Oh",  "T2u", Oh_9):  -3,
    ("Oh",  "T2u", Oh_10): 0,
    ("Oh",  "T2u", Oh_11): 1,
    ("Oh",  "T2u", Oh_12): 1,
    ("Oh",  "T2u", Oh_13): -1,
    ("Oh",  "T2u", Oh_14): -3,
    ("Oh",  "T2u", Oh_15): 0,
    ("Oh",  "T2u", Oh_16): 1,

    ("Oh",  "G1u", Oh_1):  2,
    ("Oh",  "G1u", Oh_2):  1,
    ("Oh",  "G1u", Oh_3):  0,
    ("Oh",  "G1u", Oh_4):  sqrt(2),
    ("Oh",  "G1u", Oh_5):  0,
    ("Oh",  "G1u", Oh_6):  -2,
    ("Oh",  "G1u", Oh_7):  -1,
    ("Oh",  "G1u", Oh_8):  -sqrt(2),
    ("Oh",  "G1u", Oh_9):  -2,
    ("Oh",  "G1u", Oh_10): -1,
    ("Oh",  "G1u", Oh_11): 0,
    ("Oh",  "G1u", Oh_12): -sqrt(2),
    ("Oh",  "G1u", Oh_13): 0,
    ("Oh",  "G1u", Oh_14): 2,
    ("Oh",  "G1u", Oh_15): 1,
    ("Oh",  "G1u", Oh_16): sqrt(2),

    ("Oh",  "G2u", Oh_1):  2,
    ("Oh",  "G2u", Oh_2):  1,
    ("Oh",  "G2u", Oh_3):  0,
    ("Oh",  "G2u", Oh_4):  -sqrt(2),
    ("Oh",  "G2u", Oh_5):  0,
    ("Oh",  "G2u", Oh_6):  -2,
    ("Oh",  "G2u", Oh_7):  -1,
    ("Oh",  "G2u", Oh_8):  sqrt(2),
    ("Oh",  "G2u", Oh_9):  -2,
    ("Oh",  "G2u", Oh_10): -1,
    ("Oh",  "G2u", Oh_11): 0,
    ("Oh",  "G2u", Oh_12): sqrt(2),
    ("Oh",  "G2u", Oh_13): 0,
    ("Oh",  "G2u", Oh_14): 2,
    ("Oh",  "G2u", Oh_15): 1,
    ("Oh",  "G2u", Oh_16): -sqrt(2),

    ("Oh",  "Hu",  Oh_1):  4,
    ("Oh",  "Hu",  Oh_2):  -1,
    ("Oh",  "Hu",  Oh_3):  0,
    ("Oh",  "Hu",  Oh_4):  0,
    ("Oh",  "Hu",  Oh_5):  0,
    ("Oh",  "Hu",  Oh_6):  -4,
    ("Oh",  "Hu",  Oh_7):  1,
    ("Oh",  "Hu",  Oh_8):  0,
    ("Oh",  "Hu",  Oh_9):  -4,
    ("Oh",  "Hu",  Oh_10): 1,
    ("Oh",  "Hu",  Oh_11): 0,
    ("Oh",  "Hu",  Oh_12): 0,
    ("Oh",  "Hu",  Oh_13): 0,
    ("Oh",  "Hu",  Oh_14): 4,
    ("Oh",  "Hu",  Oh_15): -1,
    ("Oh",  "Hu",  Oh_16): 0,

    ("C4v", "A1",  C4v_1): 1,
    ("C4v", "A1",  C4v_2): 1,
    ("C4v", "A1",  C4v_3): 1,
    ("C4v", "A1",  C4v_4): 1,
    ("C4v", "A1",  C4v_5): 1,
    ("C4v", "A1",  C4v_6): 1,
    ("C4v", "A1",  C4v_7): 1,

    ("C4v", "A2",  C4v_1): 1,
    ("C4v", "A2",  C4v_2): 1,
    ("C4v", "A2",  C4v_3): 1,
    ("C4v", "A2",  C4v_4): -1,
    ("C4v", "A2",  C4v_5): -1,
    ("C4v", "A2",  C4v_6): 1,
    ("C4v", "A2",  C4v_7): 1,

    ("C4v", "B1",  C4v_1): 1,
    ("C4v", "B1",  C4v_2): 1,
    ("C4v", "B1",  C4v_3): -1,
    ("C4v", "B1",  C4v_4): 1,
    ("C4v", "B1",  C4v_5): -1,
    ("C4v", "B1",  C4v_6): 1,
    ("C4v", "B1",  C4v_7): -1,

    ("C4v", "B2",  C4v_1): 1,
    ("C4v", "B2",  C4v_2): 1,
    ("C4v", "B2",  C4v_3): -1,
    ("C4v", "B2",  C4v_4): -1,
    ("C4v", "B2",  C4v_5): 1,
    ("C4v", "B2",  C4v_6): 1,
    ("C4v", "B2",  C4v_7): -1,

    ("C4v", "E",   C4v_1): 2,
    ("C4v", "E",   C4v_2): -2,
    ("C4v", "E",   C4v_3): 0,
    ("C4v", "E",   C4v_4): 0,
    ("C4v", "E",   C4v_5): 0,
    ("C4v", "E",   C4v_6): 2,
    ("C4v", "E",   C4v_7): 0,

    ("C4v", "G1",  C4v_1): 2,
    ("C4v", "G1",  C4v_2): 0,
    ("C4v", "G1",  C4v_3): sqrt(2),
    ("C4v", "G1",  C4v_4): 0,
    ("C4v", "G1",  C4v_5): 0,
    ("C4v", "G1",  C4v_6): -2,
    ("C4v", "G1",  C4v_7): -sqrt(2),

    ("C4v", "G2",  C4v_1): 2,
    ("C4v", "G2",  C4v_2): 0,
    ("C4v", "G2",  C4v_3): -sqrt(2),
    ("C4v", "G2",  C4v_4): 0,
    ("C4v", "G2",  C4v_5): 0,
    ("C4v", "G2",  C4v_6): -2,
    ("C4v", "G2",  C4v_7): sqrt(2),

    ("C2v", "A1",  C2v_1): 1,
    ("C2v", "A1",  C2v_2): 1,
    ("C2v", "A1",  C2v_3): 1,
    ("C2v", "A1",  C2v_4): 1,
    ("C2v", "A1",  C2v_5): 1,

    ("C2v", "A2",  C2v_1): 1,
    ("C2v", "A2",  C2v_2): 1,
    ("C2v", "A2",  C2v_3): -1,
    ("C2v", "A2",  C2v_4): -1,
    ("C2v", "A2",  C2v_5): 1,

    ("C2v", "B1",  C2v_1): 1,
    ("C2v", "B1",  C2v_2): -1,
    ("C2v", "B1",  C2v_3): 1,
    ("C2v", "B1",  C2v_4): -1,
    ("C2v", "B1",  C2v_5): 1,

    ("C2v", "B2",  C2v_1): 1,
    ("C2v", "B2",  C2v_2): -1,
    ("C2v", "B2",  C2v_3): -1,
    ("C2v", "B2",  C2v_4): 1,
    ("C2v", "B2",  C2v_5): 1,

    ("C2v", "G",   C2v_1): 2,
    ("C2v", "G",   C2v_2): 0,
    ("C2v", "G",   C2v_3): 0,
    ("C2v", "G",   C2v_4): 0,
    ("C2v", "G",   C2v_5): -2,

    ("C3v", "A1",  C3v_1): 1,
    ("C3v", "A1",  C3v_2): 1,
    ("C3v", "A1",  C3v_3): 1,
    ("C3v", "A1",  C3v_4): 1,
    ("C3v", "A1",  C3v_5): 1,
    ("C3v", "A1",  C3v_6): 1,

    ("C3v", "A2",  C3v_1): 1,
    ("C3v", "A2",  C3v_2): 1,
    ("C3v", "A2",  C3v_3): -1,
    ("C3v", "A2",  C3v_4): 1,
    ("C3v", "A2",  C3v_5): 1,
    ("C3v", "A2",  C3v_6): -1,

    ("C3v", "E",   C3v_1): 2,
    ("C3v", "E",   C3v_2): -1,
    ("C3v", "E",   C3v_3): 0,
    ("C3v", "E",   C3v_4): 2,
    ("C3v", "E",   C3v_5): -1,
    ("C3v", "E",   C3v_6): 0,

    ("C3v", "F1",  C3v_1): 1,
    ("C3v", "F1",  C3v_2): -1,
    ("C3v", "F1",  C3v_3): I,
    ("C3v", "F1",  C3v_4): -1,
    ("C3v", "F1",  C3v_5): 1,
    ("C3v", "F1",  C3v_6): -I,

    ("C3v", "F2",  C3v_1): 1,
    ("C3v", "F2",  C3v_2): -1,
    ("C3v", "F2",  C3v_3): -I,
    ("C3v", "F2",  C3v_4): -1,
    ("C3v", "F2",  C3v_5): 1,
    ("C3v", "F2",  C3v_6): I,

    ("C3v", "G",   C3v_1): 2,
    ("C3v", "G",   C3v_2): 1,
    ("C3v", "G",   C3v_3): 0,
    ("C3v", "G",   C3v_4): -2,
    ("C3v", "G",   C3v_5): -1,
    ("C3v", "G",   C3v_6): 0,

    ("Cs012", "A1", Cs012_1): 1,
    ("Cs012", "A1", Cs012_2): 1,
    ("Cs012", "A1", Cs012_3): 1,
    ("Cs012", "A1", Cs012_4): 1,

    ("Cs012", "A2", Cs012_1): 1,
    ("Cs012", "A2", Cs012_2): -1,
    ("Cs012", "A2", Cs012_3): -1,
    ("Cs012", "A2", Cs012_4): 1,

    ("Cs012", "F1", Cs012_1): 1,
    ("Cs012", "F1", Cs012_2): I,
    ("Cs012", "F1", Cs012_3): -I,
    ("Cs012", "F1", Cs012_4): -1,

    ("Cs012", "F2", Cs012_1): 1,
    ("Cs012", "F2", Cs012_2): -I,
    ("Cs012", "F2", Cs012_3): I,
    ("Cs012", "F2", Cs012_4): -1,

    ("Cs112", "A1", Cs112_1): 1,
    ("Cs112", "A1", Cs112_2): 1,
    ("Cs112", "A1", Cs112_3): 1,
    ("Cs112", "A1", Cs112_4): 1,

    ("Cs112", "A2", Cs112_1): 1,
    ("Cs112", "A2", Cs112_2): -1,
    ("Cs112", "A2", Cs112_3): -1,
    ("Cs112", "A2", Cs112_4): 1,

    ("Cs112", "F1", Cs112_1): 1,
    ("Cs112", "F1", Cs112_2): I,
    ("Cs112", "F1", Cs112_3): -I,
    ("Cs112", "F1", Cs112_4): -1,

    ("Cs112", "F2", Cs112_1): 1,
    ("Cs112", "F2", Cs112_2): -I,
    ("Cs112", "F2", Cs112_3): I,
    ("Cs112", "F2", Cs112_4): -1,
}

_IRREP_MATRICES = {
    ("Oh",  "A1g", E):   eye(1),
    ("Oh",  "A1g", C4y): Matrix([1]),
    ("Oh",  "A1g", C4z): Matrix([1]),
    ("Oh",  "A1g", Is):  eye(1),

    ("Oh",  "A1u", E):   eye(1),
    ("Oh",  "A1u", C4y): Matrix([1]),
    ("Oh",  "A1u", C4z): Matrix([1]),
    ("Oh",  "A1u", Is):  -eye(1),

    ("Oh",  "A2g", E):   eye(1),
    ("Oh",  "A2g", C4y): Matrix([-1]),
    ("Oh",  "A2g", C4z): Matrix([-1]),
    ("Oh",  "A2g", Is):  eye(1),

    ("Oh",  "A2u", E):   eye(1),
    ("Oh",  "A2u", C4y): Matrix([-1]),
    ("Oh",  "A2u", C4z): Matrix([-1]),
    ("Oh",  "A2u", Is):  -eye(1),

    ("Oh",  "Eg",  E):   eye(2),

    ("Oh",  "Eg",  C4y): Rational(1,2)*Matrix([[   1   ,sqrt(3)],
                                               [sqrt(3),  -1   ]]),
    ("Oh",  "Eg",  C4z): Matrix([[-1, 0],
                                 [ 0, 1]]),
    ("Oh",  "Eg",  Is):  eye(2),


    ("Oh",  "Eu",  E):   eye(2),

    ("Oh",  "Eu",  C4y): Rational(1,2)*Matrix([[   1   ,sqrt(3)],
                                               [sqrt(3),  -1   ]]),
    ("Oh",  "Eu",  C4z): Matrix([[-1, 0],
                                 [ 0, 1]]),
    ("Oh",  "Eu",  Is):  -eye(2),


    ("Oh",  "T1g", E):   eye(3),


    ("Oh",  "T1g", C4y): Matrix([[ 0, 0, 1],
                                 [ 0, 1, 0],
                                 [-1, 0, 0]]),
    ("Oh",  "T1g", C4z): Matrix([[ 0,-1, 0],
                                 [ 1, 0, 0],
                                 [ 0, 0, 1]]),
    ("Oh",  "T1g", Is):  eye(3),



    ("Oh",  "T1u", E):  eye(3),


    ("Oh",  "T1u", C4y): Matrix([[ 0, 0, 1],
                                 [ 0, 1, 0],
                                 [-1, 0, 0]]),
    ("Oh",  "T1u", C4z): Matrix([[ 0,-1, 0],
                                 [ 1, 0, 0],
                                 [ 0, 0, 1]]),
    ("Oh",  "T1u", Is): -eye(3),



    ("Oh",  "T2g", E):   eye(3),


    ("Oh",  "T2g", C4y): Matrix([[ 0, 0,-1],
                                 [ 0,-1, 0],
                                 [ 1, 0, 0]]),
    ("Oh",  "T2g", C4z): Matrix([[ 0, 1, 0],
                                 [-1, 0, 0],
                                 [ 0, 0,-1]]),
    ("Oh",  "T2g", Is):  eye(3),



    ("Oh",  "T2u", E):   eye(3),


    ("Oh",  "T2u", C4y): Matrix([[ 0, 0,-1],
                                 [ 0,-1, 0],
                                 [ 1, 0, 0]]),
    ("Oh",  "T2u", C4z): Matrix([[ 0, 1, 0],
                                 [-1, 0, 0],
                                 [ 0, 0,-1]]),
    ("Oh",  "T2u", Is):  -eye(3),



    ("Oh",  "G1g", E):   eye(2),

    ("Oh",  "G1g", C4y): 1/sqrt(2)*Matrix([[ 1,-1],
                                           [ 1, 1]]),
    ("Oh",  "G1g", C4z): 1/sqrt(2)*Matrix([[1-I,  0  ],
                                           [ 0 , 1+I]]),
    ("Oh",  "G1g", Is):  eye(2),


    ("Oh",  "G1u", E):   eye(2),

    ("Oh",  "G1u", C4y): 1/sqrt(2)*Matrix([[ 1,-1],
                                           [ 1, 1]]),
    ("Oh",  "G1u", C4z): 1/sqrt(2)*Matrix([[1-I,  0  ],
                                           [ 0 , 1+I]]),
    ("Oh",  "G1u", Is):  -eye(2),


    ("Oh",  "G2g", E):   eye(2),

    ("Oh",  "G2g", C4y): -1/sqrt(2)*Matrix([[ 1,-1],
                                           [ 1, 1]]),
    ("Oh",  "G2g", C4z): -1/sqrt(2)*Matrix([[1-I,  0  ],
                                           [ 0 , 1+I]]),
    ("Oh",  "G2g", Is):  eye(2),


    ("Oh",  "G2u", E):   eye(2),

    ("Oh",  "G2u", C4y): -1/sqrt(2)*Matrix([[ 1,-1],
                                           [ 1, 1]]),
    ("Oh",  "G2u", C4z): -1/sqrt(2)*Matrix([[1-I,  0  ],
                                           [ 0 , 1+I]]),
    ("Oh",  "G2u", Is):  -eye(2),




    ("Oh",  "Hg",  E):   eye(4),
    


    ("Oh",  "Hg",  C4z): sqrt(2)/4*Matrix([[   1   ,-sqrt(3), sqrt(3),   -1   ],
                                          [ sqrt(3),   -1   ,   -1   , sqrt(3)],
                                          [ sqrt(3),    1   ,   -1   ,-sqrt(3)],
                                          [    1   , sqrt(3), sqrt(3),    1   ]]),
    ("Oh",  "Hg",  C4y): 1/sqrt(2)*Matrix([[-1-I,  0 ,  0 ,  0 ],
                                           [  0 , 1-I,  0 ,  0 ],
                                           [  0 ,  0 , 1+I,  0 ],
                                           [  0 ,  0 ,  0 ,-1+I]]),
    ("Oh",  "Hg",  Is):  eye(4),




    ("Oh",  "Hu",  E):   eye(4),
    


    ("Oh",  "Hu",  C4z): sqrt(2)/4*Matrix([[   1   ,-sqrt(3), sqrt(3),   -1   ],
                                          [ sqrt(3),   -1   ,   -1   , sqrt(3)],
                                          [ sqrt(3),    1   ,   -1   ,-sqrt(3)],
                                          [    1   , sqrt(3), sqrt(3),    1   ]]),
    ("Oh",  "Hu",  C4y): 1/sqrt(2)*Matrix([[-1-I,  0 ,  0 ,  0 ],
                                           [  0 , 1-I,  0 ,  0 ],
                                           [  0 ,  0 , 1+I,  0 ],
                                           [  0 ,  0 ,  0 ,-1+I]]),
    ("Oh",  "Hu",  Is):  -eye(4),





    ("C4v", "A1",  C4z): Matrix([1]),
    ("C4v", "A1",  I2y): Matrix([1]),

    ("C4v", "A2",  C4z): Matrix([ 1]),
    ("C4v", "A2",  I2y): Matrix([-1]),

    ("C4v", "B1",  C4z): Matrix([-1]),
    ("C4v", "B1",  I2y): Matrix([ 1]),

    ("C4v", "B2",  C4z): Matrix([-1]),
    ("C4v", "B2",  I2y): Matrix([-1]),


    ("C4v", "E",   C4z): Matrix([[ 0,-1],
                                 [ 1, 0]]),
    ("C4v", "E",   I2y): Matrix([[ 1, 0],
                                 [ 0,-1]]),


    ("C4v", "G1",  C4z): Matrix([[(1-I)/sqrt(2),      0      ],
                                 [      0      ,(1+I)/sqrt(2)]]),
    ("C4v", "G1",  I2y): Matrix([[ 0,-1],
                                 [ 1, 0]]),


    ("C4v", "G2",  C4z): Matrix([[ (I-1)/sqrt(2),       0      ],
                                 [       0      ,-(1+I)/sqrt(2)]]),
    ("C4v", "G2",  I2y): Matrix([[ 0,-1],
                                 [ 1, 0]]),



    ("C2v", "A1",  C2e): Matrix([1]),
    ("C2v", "A1",  I2f): Matrix([1]),

    ("C2v", "A2",  C2e): Matrix([ 1]),
    ("C2v", "A2",  I2f): Matrix([-1]),

    ("C2v", "B1",  C2e): Matrix([-1]),
    ("C2v", "B1",  I2f): Matrix([ 1]),

    ("C2v", "B2",  C2e): Matrix([-1]),
    ("C2v", "B2",  I2f): Matrix([-1]),


    ("C2v", "G",   C2e): Matrix([[-I/sqrt(2),-1/sqrt(2)],
                                 [ 1/sqrt(2), I/sqrt(2)]]),
    ("C2v", "G",   I2f): Matrix([[ I/sqrt(2),-1/sqrt(2)],
                                 [ 1/sqrt(2),-I/sqrt(2)]]),


    ("C3v", "A1",  C3D): Matrix([1]),
    ("C3v", "A1",  I2b): Matrix([1]),

    ("C3v", "A2",  C3D): Matrix([ 1]),
    ("C3v", "A2",  I2b): Matrix([-1]),


    ("C3v", "E",   C3D): Rational(1,2)*Matrix([[   -1   , sqrt(3)],
                                               [-sqrt(3),   -1   ]]),
    ("C3v", "E",   I2b): Matrix([[-1, 0],
                                 [ 0, 1]]),


    ("C3v", "F1",  C3D): Matrix([-1]),
    ("C3v", "F1",  I2b): Matrix([ I]),

    ("C3v", "F2",  C3D): Matrix([-1]),
    ("C3v", "F2",  I2b): Matrix([-I]),


    ("C3v", "G",   C3D): Rational(1,2)*Matrix([[ 1-I,-1-I],
                                               [ 1-I, 1+I]]),
    ("C3v", "G",   I2b): 1/sqrt(2)*Matrix([[  0 , 1-I],
                                           [-1-I,  0]]),


    ("Cs012", "A1", I2x): Matrix([ 1]),
    ("Cs012", "A2", I2x): Matrix([-1]),
    ("Cs012", "F1", I2x): Matrix([ I]),
    ("Cs012", "F2", I2x): Matrix([-I]),


    ("Cs112", "A1", I2b): Matrix([ 1]),
    ("Cs112", "A2", I2b): Matrix([-1]),
    ("Cs112", "F1", I2b): Matrix([ I]),
    ("Cs112", "F2", I2b): Matrix([-I]),
}

_GENERATORS = {
    "Oh": {
        C2y:   [C4y,   C4y],
        EC4yi: [C2y,   C4y],
        EE:    [EC4yi, C4y],
        EC4y:  [EE,    C4y],
        EC2y:  [EC4y,  C4y],
        C4yi:  [EC2y,  C4y],
        C2z:   [C4z,   C4z],
        EC4zi: [C2z,   C4z],
        EC4z:  [EE,    C4z],
        EC2z:  [EC4z,  C4z],
        C4zi:  [EC2z,  C4z],
        C3D:   [C4y,   C4z],
        C4x:   [C4zi,  C3D],
        C2x:   [C4x,   C4x],
        EC4xi: [C2x,   C4x],
        EC4x:  [EE,    C4x],
        EC2x:  [EC4x,  C4x],
        C4xi:  [EC2x,  C4x],
        C3A:   [C4yi,  C4z],
        C3B:   [C4y,   C4zi],
        C3C:   [C4yi,  C4zi],
        C3Ai:  [C4zi,  C4y],
        C3Bi:  [C4z,   C4yi],
        C3Ci:  [C4z,   C4y],
        C3Di:  [C4zi,  C4yi],
        C2a:   [C2y,   C4z],
        C2b:   [C2x,   C4z],
        C2c:   [C4y,   C2z],
        C2d:   [C2z,   C4y],
        C2e:   [C2z,   C4x],
        C2f:   [C2y,   C4x],
        EC3A:  [EE,    C3A],
        EC3B:  [EE,    C3B],
        EC3C:  [EE,    C3C],
        EC3D:  [EE,    C3D],
        EC3Ai: [EE,    C3Ai],
        EC3Bi: [EE,    C3Bi],
        EC3Ci: [EE,    C3Ci],
        EC3Di: [EE,    C3Di],
        EC2a:  [EE,    C2a],
        EC2b:  [EE,    C2b],
        EC2c:  [EE,    C2c],
        EC2d:  [EE,    C2d],
        EC2e:  [EE,    C2e],
        EC2f:  [EE,    C2f],
        I4y:   [Is,    C4y],
        I4z:   [Is,    C4z],
        I2y:   [Is,    C2y],
        EI4yi: [Is,    EC4yi],
        EIs:   [Is,    EE],
        EI4y:  [Is,    EC4y],
        EI2y:  [Is,    EC2y],
        I4yi:  [Is,    C4yi],
        I2z:   [Is,    C2z],
        EI4zi: [Is,    EC4zi],
        EI4z:  [Is,    EC4z],
        EI2z:  [Is,    EC2z],
        I4zi:  [Is,    C4zi],
        I3D:   [Is,    C3D],
        I4x:   [Is,    C4x],
        I2x:   [Is,    C2x],
        EI4xi: [Is,    EC4xi],
        EI4x:  [Is,    EC4x],
        EI2x:  [Is,    EC2x],
        I4xi:  [Is,    C4xi],
        I3A:   [Is,    C3A],
        I3B:   [Is,    C3B],
        I3C:   [Is,    C3C],
        I3Ai:  [Is,    C3Ai],
        I3Bi:  [Is,    C3Bi],
        I3Ci:  [Is,    C3Ci],
        I3Di:  [Is,    C3Di],
        I2a:   [Is,    C2a],
        I2b:   [Is,    C2b],
        I2c:   [Is,    C2c],
        I2d:   [Is,    C2d],
        I2e:   [Is,    C2e],
        I2f:   [Is,    C2f],
        EI3A:  [Is,    EC3A],
        EI3B:  [Is,    EC3B],
        EI3C:  [Is,    EC3C],
        EI3D:  [Is,    EC3D],
        EI3Ai: [Is,    EC3Ai],
        EI3Bi: [Is,    EC3Bi],
        EI3Ci: [Is,    EC3Ci],
        EI3Di: [Is,    EC3Di],
        EI2a:  [Is,    EC2a],
        EI2b:  [Is,    EC2b],
        EI2c:  [Is,    EC2c],
        EI2d:  [Is,    EC2d],
        EI2e:  [Is,    EC2e],
        EI2f:  [Is,    EC2f],
    },

    "C4v": {
        C2z:   [C4z,  C4z],
        EC4zi: [C4z,  C2z],
        EE:    [C4z,  EC4zi],
        EC4z:  [C4z,  EE],
        EC2z:  [C4z,  EC4z],
        C4zi:  [C4z,  EC2z],
        E:     [C4z,  C4zi],
        I2x:   [I2y,  C2z],
        I2a:   [I2y,  C4z],
        I2b:   [I2y,  EC4zi],
        EI2y:  [I2y,  EE],
        EI2x:  [I2y,  EC2z],
        EI2a:  [I2y,  EC4z],
        EI2b:  [I2y,  C4zi],
    },

    "C3v": {
        EC3Di: [C3D, C3D],
        EE:    [C3D, EC3Di],
        EC3D:  [C3D, EE],
        C3Di:  [C3D, EC3D],
        E:     [C3D, C3Di],
        I2d:   [I2b, EC3Di],
        I2f:   [I2b, EC3D],
        EI2b:  [I2b, EE],
        EI2d:  [I2b, C3Di],
        EI2f:  [I2b, C3D],
    },

    "C2v": {
        EE:   [C2e, C2e],
        EC2e: [C2e, EE],
        E:    [C2e, EC2e],
        I2x:  [I2f, C2e],
        EI2x: [I2f, EC2e],
        EI2f: [I2f, EE],
    },

    "Cs012": {
        EE:   [I2x, I2x],
        EI2x: [I2x, EE],
        E:    [EE,  EE],
    },

    "Cs112": {
        EE:   [I2b, I2b],
        EI2b: [I2b, EE],
        E:    [EE,  EE],
    },
}

# Construct the rest of the irrep matrices
for little_group, bos_lg_irreps in _BOSONIC_LITTLE_GROUP_IRREPS.items():
  for bos_lg_irrep in bos_lg_irreps:
    for lg_element, prod in _GENERATORS[little_group].items():
      _IRREP_MATRICES[(little_group, bos_lg_irrep, lg_element)] = \
          _IRREP_MATRICES[(little_group, bos_lg_irrep, prod[0])]*_IRREP_MATRICES[(little_group, bos_lg_irrep, prod[1])]


for little_group, fer_lg_irreps in _FERMIONIC_LITTLE_GROUP_IRREPS.items():
  for fer_lg_irrep in fer_lg_irreps:
    for lg_element, prod in _GENERATORS[little_group].items():
      _IRREP_MATRICES[(little_group, fer_lg_irrep, lg_element)] = \
          _IRREP_MATRICES[(little_group, fer_lg_irrep, prod[0])]*_IRREP_MATRICES[(little_group, fer_lg_irrep, prod[1])]


class LittleGroup:

  LITTLE_GROUPS = {
      P(0,0,0): "Oh",
      P(0,0,1): "C4v",
      P(0,1,1): "C2v",
      P(1,1,1): "C3v",
      P(0,1,2): "Cs012",
      P(1,1,2): "Cs112",
  }

  CHARACTERS = _CHARACTERS
  IRREP_MATRICES = _IRREP_MATRICES
  GENERATORS = _GENERATORS
  REFERENCE_ROTATIONS = _REFERENCE_ROTATIONS

  def __init__(self, momentum=P0):
    self._momentum = momentum

  _instances = {}
  def __new__(cls, momentum=P0):
    if momentum not in cls._instances:
      new_instance = super().__new__(cls)
      new_instance._setup(momentum)
      cls._instances[momentum] = new_instance

    return cls._instances[momentum]

  @property
  def order(self):
    return len(self.elements)

  @property
  def momentum(self):
    return self._momentum

  def getCharacter(self, irrep, element):
    conj_class = frozenset([self.wigner_rotation(el) for el in self.conjugacy_class(element)])
    return _CHARACTERS[(self.little_group, irrep, conj_class)]

  def conjugacy_class(self, element):
    if element not in self.elements:
      raise ValueError("Element not in Little Group")
    
    return self._conj_classes[element]

  @property
  def equivalent_momenta(self):
    return self._equivalent_momenta

  @property
  def elements(self):
    return self._elements

  def wigner_rotation(self, element):
    return self._wigner_rotations[element]

  def irrep_matrix(self, irrep, element):
    if irrep[-1] == 'm' or irrep[-1] == 'p':
      irrep = irrep[:-1]
    return self.IRREP_MATRICES[(self.little_group, irrep, self.wigner_rotation(element))]

  def reference_rotation(self, momentum):
    return self.REFERENCE_ROTATIONS[momentum.reduced]

  def _setup(self, momentum):
    # get elements of little group and equivalent momentum frames
    self._elements = set()
    self._equivalent_momenta = set()
    for rotation in point_group:
      mom_prime = rotation * momentum
      if mom_prime == momentum:
        self._elements.add(rotation)
      else:
        self._equivalent_momenta.add(mom_prime)

    # get the wigner rotations associated with each element
    self._wigner_rotations = dict()
    ref_rotation = _REFERENCE_ROTATIONS[momentum.reduced]
    for element in point_group:
      mom_rotated = element*momentum
      wigner_rotation = ~_REFERENCE_ROTATIONS[mom_rotated.reduced] * element * ref_rotation
      self._wigner_rotations[element] = wigner_rotation

    # setupt the conjugacy classes
    conj_classes = set()
    for el1 in self.elements:
      conj_class = set()
      for el2 in self.elements:
        conj_class.add(el2*el1*~el2)
      conj_class = frozenset(conj_class)
      conj_classes.add(conj_class)

    self._conj_classes = dict()
    for conj_class in conj_classes:
      for element in conj_class:
        self._conj_classes[element] = conj_class

    # setup irreps
    lg_str = self.LITTLE_GROUPS[momentum.reduced_pref]
    _bosonic_irreps = _BOSONIC_LITTLE_GROUP_IRREPS[lg_str]
    _fermionic_irreps = _FERMIONIC_LITTLE_GROUP_IRREPS[lg_str]

    self._bosonic_irreps = dict()
    self._fermionic_irreps = dict()
    for bosonic_irrep in _bosonic_irreps:
      self._bosonic_irreps[bosonic_irrep] = _BOSONIC_LITTLE_GROUP_IRREPS_ROWS[(lg_str, bosonic_irrep)]

    for fermionic_irrep in _fermionic_irreps:
      self._fermionic_irreps[fermionic_irrep] = _FERMIONIC_LITTLE_GROUP_IRREPS_ROWS[(lg_str, fermionic_irrep)]

  @property
  def irreps(self):
    return {**self._bosonic_irreps, **self._fermionic_irreps}

  @property
  def bosonic_irreps(self):
    return self._bosonic_irreps

  @property
  def fermionic_irreps(self):
    return self._fermionic_irreps

  def irrep_rows(self, irrep):
    if irrep[-1] == 'm' or irrep[-1] == 'p':
      irrep = irrep[:-1]
    return self.irreps[irrep]

  @property
  def little_group(self):
    return self.LITTLE_GROUPS[self.momentum.reduced_pref]

  @property
  def is_reference_momentum(self):
    return self.momentum.reduced == self.momentum.reduced_pref


'''
# Spinor Representations
I = Identity(4)
g1 = MatrixSymbol('g1', 4, 4)
g2 = MatrixSymbol('g2', 4, 4)
g3 = MatrixSymbol('g3', 4, 4)

_rotation_map = {
    Angle.HALF: {
        Axis.x: g2*g3,
        Axis.y: g3*g1,
        Axis.z: g1*g2,
        Axis.a:  1/sqrt(2)*(g2*g3 + g3*g1),
        Axis.b:  1/sqrt(2)*(g2*g3 - g3*g1),
        Axis.c:  1/sqrt(2)*(g2*g3 + g1*g2),
        Axis.d: -1/sqrt(2)*(g2*g3 - g1*g2),
        Axis.e:  1/sqrt(2)*(g3*g1 + g1*g2),
        Axis.f:  1/sqrt(2)*(g3*g1 - g1*g2)
    },
    Angle.THIRD: {
        Axis.A: (I - g2*g3 - g3*g1 + g1*g2)/2,
        Axis.B: (I - g2*g3 + g3*g1 - g1*g2)/2,
        Axis.C: (I + g2*g3 - g3*g1 - g1*g2)/2,
        Axis.D: (I + g2*g3 + g3*g1 + g1*g2)/2
    },
    Angle.INV_THIRD: {
        Axis.A: (I + g2*g3 + g3*g1 - g1*g2)/2,
        Axis.B: (I + g2*g3 - g3*g1 + g1*g2)/2,
        Axis.C: (I - g2*g3 + g3*g1 + g1*g2)/2,
        Axis.D: (I - g2*g3 - g3*g1 - g1*g2)/2
    },
    Angle.QUARTER: {
        Axis.x: 1/sqrt(2)*(I + g2*g3),
        Axis.y: 1/sqrt(2)*(I + g3*g1),
        Axis.z: 1/sqrt(2)*(I + g1*g2)
    },
    Angle.INV_QUARTER: {
        Axis.x: 1/sqrt(2)*(I - g2*g3),
        Axis.y: 1/sqrt(2)*(I - g3*g1),
        Axis.z: 1/sqrt(2)*(I - g1*g2)
    }
}

_conjugate_rotation_map = {
    Angle.HALF: {
        Axis.x: -g2*g3,
        Axis.y: -g3*g1,
        Axis.z: -g1*g2,
        Axis.a: -1/sqrt(2)*(g2*g3 + g3*g1),
        Axis.b: -1/sqrt(2)*(g2*g3 - g3*g1),
        Axis.c: -1/sqrt(2)*(g2*g3 + g1*g2),
        Axis.d:  1/sqrt(2)*(g2*g3 - g1*g2),
        Axis.e: -1/sqrt(2)*(g3*g1 + g1*g2),
        Axis.f: -1/sqrt(2)*(g3*g1 - g1*g2)
    },
    Angle.THIRD: {
        Axis.A: (I + g2*g3 + g3*g1 - g1*g2)/2,
        Axis.B: (I + g2*g3 - g3*g1 + g1*g2)/2,
        Axis.C: (I - g2*g3 + g3*g1 + g1*g2)/2,
        Axis.D: (I - g2*g3 - g3*g1 - g1*g2)/2
    },
    Angle.INV_THIRD: {
        Axis.A: (I - g2*g3 - g3*g1 + g1*g2)/2,
        Axis.B: (I - g2*g3 + g3*g1 - g1*g2)/2,
        Axis.C: (I + g2*g3 - g3*g1 - g1*g2)/2,
        Axis.D: (I + g2*g3 + g3*g1 + g1*g2)/2
    },
    Angle.QUARTER: {
        Axis.x: 1/sqrt(2)*(I - g2*g3),
        Axis.y: 1/sqrt(2)*(I - g3*g1),
        Axis.z: 1/sqrt(2)*(I - g1*g2)
    },
    Angle.INV_QUARTER: {
        Axis.x: 1/sqrt(2)*(I + g2*g3),
        Axis.y: 1/sqrt(2)*(I + g3*g1),
        Axis.z: 1/sqrt(2)*(I + g1*g2)
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
'''
