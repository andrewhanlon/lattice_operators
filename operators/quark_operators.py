from collections import defaultdict
import itertools
import functools
from operator import mul

from sympy import expand
from sympy import AtomicExpr, S
from sympy import Add, Mul
from sympy import get_indices, Idx

from .cubic_rotations import E, P0
from .grassmann import grassmann_simplify, grassmann_coefficients, GrassmannSymbol
from .vector import IndexedVector
from .quarks import Quark, IndexedQuark


class Operator(AtomicExpr):
  is_commutative = False

  @property
  def rotated_to(self):
    return self._rotated_to




class QuarkOperator(Operator):
  """QuarkOperator class

  A QuarkOperator corresponds to an object with one spatial index.
  """

  def __init__(self, op_expr):
    self._op_expr = op_expr
    self._rotated_to = E
    self._coefficients = dict()

  def __new__(cls, op_expr):
    outer_indices, index_syms = get_indices(op_expr)

    outer_indices = {index for index in outer_indices if isinstance(index, Idx)}

    if outer_indices:
      raise TypeError("Operators with free indices not yet supported")

    op_simp = grassmann_simplify(operator_replace(op_expr))
    if op_simp == S.Zero:
      return S.Zero

    obj = super().__new__(cls, op_expr)
    obj._simplified = {E: op_simp}

    return obj

  @property
  def expr(self):
    return self._op_expr

  @property
  def simplified(self):
    return self._simplified[self._rotated_to]

  @classmethod
  def _rotate(cls, element, expr):
    if isinstance(expr, IndexedQuark) and isinstance(expr.base, Quark):
      expr.base.rotateTo(element)
      return expr

    elif isinstance(expr, Quark):
      expr.rotateTo(element)
      return expr

    elif isinstance(expr, GrassmannSymbol):
      raise TypeError("Cannot rotate a GrassmannSymbol")

    elif expr.is_Atom:
      return expr

    else:
      args = [cls._rotate(element, arg) for arg in expr.args]
      return expr.func(*args)

  def rotateTo(self, element):
    self._rotated_to = element

    if element not in self._simplified:
      rot_expr = self.__class__._rotate(element, self.expr)
      self._simplified[element] = grassmann_simplify(operator_replace(rot_expr))

  @classmethod
  def _number_of_quarks(cls, expr):
    if isinstance(expr, IndexedQuark) and isinstance(expr.base, Quark):
      return 1

    elif isinstance(expr, Quark):
      return 1

    elif isinstance(expr, Add):
      quarks_per_term = set([cls._number_of_quarks(arg) for arg in expr.args])
      if len(quarks_per_term) != 1:
        raise ValueError("All terms should have the same number of Quarks")

      return quarks_per_term.pop()

    elif expr.is_Atom:
      return 0

    else:
      return sum([cls._number_of_quarks(arg) for arg in expr.args])

  @property
  def number_of_quarks(self):
    try:
      return self._quarks
    except AttributeError:
      self._quarks = self.__class__._number_of_quarks(self.expr)
      return self._quarks


  @property
  def bosonic(self):
    return self.number_of_quarks % 2 == 0

  @property
  def fermionic(self):
    return self.number_of_quarks % 2 == 1

  @property
  def coefficients(self):
    if self.rotated_to not in self._coefficients:

      self._coefficients[self.rotated_to] = defaultdict(int)
      for grassmann_term, coeff in grassmann_coefficients(self.simplified).items():
        self._coefficients[self.rotated_to][(grassmann_term.__repr__(),)] = coeff

    return self._coefficients[self.rotated_to]

  def __repr__(self):
    return self.simplified.__repr__()



class MomProjOperator(Operator):

  def __init__(self, operator, momentum):
    operator.rotateTo(E)
    self._operator = operator
    self._original_momentum = momentum
    self._current_momentum = momentum

  @property
  def rotated_to(self):
    _rotated_to = self.operator.rotated_to
    if self._current_momentum != (_rotated_to * self._original_momentum):
      raise ValueError("Momentum and operator differently rotated...")

    return _rotated_to

  @property
  def operator(self):
    return self._operator

  @property
  def expr(self):
    return self.operator.expr

  @property
  def simplified(self):
    return self.operator.simplified

  @property
  def momentum(self):
    return self._current_momentum

  def rotateTo(self, element):
    self._current_momentum = element * self._original_momentum
    self.operator.rotateTo(element)

  @property
  def number_of_quarks(self):
    return self.operator.number_of_quarks

  @property
  def bosonic(self):
    return self.operator.bosonic

  @property
  def fermionic(self):
    return self.operator.fermionic

  @property
  def coefficients(self):
    coefficients = self.operator.coefficients.copy()
    coefficients = dict((((self.momentum.__repr__(), grass_term[0]),), coeff) for (grass_term, coeff) in coefficients.items())
    return coefficients

  def __repr__(self):
    return "{} * ({})".format(self.momentum, self.operator.__repr__())


def operator_replace(op_expr):
  if isinstance(op_expr, IndexedQuark):
    return IndexedVector(op_expr.base.vector, *op_expr.indices)

  elif op_expr.is_Atom:
    return op_expr
  
  else:
    args = [operator_replace(arg) for arg in op_expr.args]
    return op_expr.func(*args)

def rotate_to(op_expr, element):
  if isinstance(op_expr, QuarkOperator) or isinstance(op_expr, MomProjOperator):
    op_expr.rotateTo(element)
    return op_expr

  elif op_expr.is_Atom:
    return op_expr

  else:
    rot_args = [rotate_to(arg, element) for arg in op_expr.args]
    return op_expr.func(*rot_args)

def operator_coefficients(op_expr):
  if isinstance(op_expr, QuarkOperator) or isinstance(op_expr, MomProjOperator):
    return op_expr.coefficients

  elif isinstance(op_expr, Mul):
    coeff_args = list()
    overall_coeff = S.One
    for arg in op_expr.args:
      if arg.is_number:
        overall_coeff *= arg
        continue

      coeff_arg = operator_coefficients(arg)
      coeff_args.append(coeff_arg.items())

    coeffs = dict()
    for prod in itertools.product(*coeff_args):
      final_mul = overall_coeff
      coeff = list(zip(*prod))
      term = sum(coeff[0], tuple())
      if term[0] == term[1]:
        continue

      # @ADH - TEMP FIX!!
      if term[1] < term[0]:
        final_mul *= -S.One
        term = (term[1], term[0])

      value = functools.reduce(mul, coeff[1])
      if final_mul != S.One:
        value *= final_mul

      if term in coeffs:
        coeffs[term] += value
      else:
        coeffs[term] = value

    coeffs = { k:v for k, v in coeffs.items() if expand(v) }

    return coeffs

  elif isinstance(op_expr, Add):
    coeffs = dict()
    for arg in op_expr.args:
      coeff_arg = operator_coefficients(arg)
      for term, coeff in coeff_arg.items():
        if term in coeffs:
          coeffs[term] += coeff
        else:
          coeffs[term] = coeff

    coeffs = { k:v for k, v in coeffs.items() if expand(v) }

    return coeffs

def operator_terms(op_expr):
  return operator_coefficients(op_expr).keys()

def project_momentum(op_expr, *momenta):
  if isinstance(op_expr, QuarkOperator):
    return MomProjOperator(op_expr, *momenta)
  
  elif isinstance(op_expr, MomProjOperator):
    return MomProjOperator(op_expr.operator, *momenta)

  elif isinstance(op_expr, Add):
    proj_args = (project_momentum(arg, *momenta) for arg in op_expr.args)
    return op_expr.func(*proj_args)

  elif isinstance(op_expr, Mul):
    proj_args = list()
    mom_ind = 0
    for arg in op_expr.args:
      num_ops = number_of_operators(arg)
      if num_ops == 0:
        proj_args.append(arg)
      else:
        reduced_momenta = momenta[mom_ind:mom_ind+num_ops]
        proj_args.append(project_momentum(arg, *reduced_momenta))
        mom_ind += num_ops

    if mom_ind != len(momenta):
      raise ValueError("Number of operators does not match with number of momenta provided")

    return op_expr.func(*proj_args)

  else:
    return op_expr

# @ADH - add note about how it is the number of operators in each expanded term...
def number_of_operators(op_expr):
  if isinstance(op_expr, QuarkOperator) or isinstance(op_expr, MomProjOperator):
    return 1

  elif isinstance(op_expr, Add):
    ops_per_term = set([number_of_operators(arg) for arg in op_expr.args])
    if len(ops_per_term) != 1:
      raise ValueError("All terms should have the same number of Operators")

    return ops_per_term.pop()

  elif isinstance(op_expr, Mul):
    return sum([number_of_operators(arg) for arg in op_expr.args])

  else:
    return 0

def bosonic(op_expr):
  return number_of_operators(op_expr) % 2 == 0

def fermionic(op_expr):
  return number_of_operators(op_expr) % 2 == 1

# Assumes non-momentum projected objects are at zero momentum
def momentum(op_expr):
  if isinstance(op_expr, MomProjOperator):
    return op_expr.momentum

  elif isinstance(op_expr, Add):
    mom_per_term = set([momentum(arg) for arg in op_expr.args])
    if len(mom_per_term) != 1:
      raise ValueError("All terms should have the same momentum")

    return mom_per_term.pop()

  elif isinstance(op_expr, Mul):
    return sum([momentum(arg) for arg in op_expr.args], P0)

  else:
    return P0

