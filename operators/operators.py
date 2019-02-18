from abc import ABCMeta, abstractmethod
from operator import mul
from functools import reduce
import itertools

from sympy import Indexed, IndexedBase
from sympy.tensor.indexed import IndexException
from sympy import AtomicExpr
from sympy import Add, Mul
from sympy import S, expand
from sympy.core.compatibility import is_sequence

from .grassmann import contract_indices, GrassmannSymbol
from .cubic_rotations import P0, P, E

"""
TODO: Doesn't handel Pow correctly...
"""

class Operator(AtomicExpr):

  is_commutative = False
  _momentum = P0

  @property
  def bosonic(self):
    pass

  @property
  def fermionic(self):
    pass

  @property
  def rotated_to(self):
    return self._rotated_to

  @abstractmethod
  def rotateTo(self, element):
    pass

  def rotateBy(self, element):
    new_element = element * self.rotated_to
    self.rotateTo(new_element)

  @property
  @abstractmethod
  def coefficients(self):
    pass

  @property
  @abstractmethod
  def statistics(self):
    pass

  @property
  def momentum(self):
    return self._momentum

  @abstractmethod
  def P(self, *momentum):
    pass


class OperatorVector(Operator, IndexedBase):

  is_commutative = False

  @property
  def vector(self):
    return self._vectors[self.rotated_to]

  def __getitem__(self, indices):
    """
    TODO:
      Check the indices?
    """
    if is_sequence(indices):
      # Special case needed because M[*my_tuple] is a syntax error.
      '''
      if self.shape and len(self.shape) != len(indices):
        raise IndexException("Rank mismatch.")
      '''
      return IndexedOperator(self, *indices)
    else:
      '''
      if self.shape and len(self.shape) != 1:
        raise IndexException("Rank mismatch.")
      '''
      return IndexedOperator(self, indices)

  def __call__(self, *indices):
    return self.vector[indices]


class IndexedOperator(Operator, Indexed):

  is_commutative = False

  def rotateTo(self, element):
    self.base.rotateTo(element)

  @property
  def statistics(self):
    return self.base.statistics

  @property
  def momentum(self):
    return self.base.momentum

  def P(self, *momentum):
    new_base = self.base.P(*momentum)
    return new_base[self.indices[0]]

  '''
  def __repr__(self):
    return "{}.{}".format(super().__repr__(), self.momentum.__repr__())

  def __str__(self):
    return "{}.{}".format(super().__str__(), self.momentum.__str__())
  '''

class OperatorElement(Operator):

  is_commutative = False

  def P(self, *momentum):
    self._momentum = P(*momentum)



def rotate_to(expr, element):
  if isinstance(expr, Operator):
    expr.rotateTo(element)
    return expr

  elif expr.is_Atom:
    return expr

  else:
    rotated_args = [rotate_to(arg, element) for arg in expr.args]
    return expr.func(*rotated_args)

def terms(expr):
  return coefficients(expr).keys()

def coefficients(expr):
  repl_expr = _replace(expr)
  cont_expr = contract_indices(repl_expr)
  coeffs = _coefficients(cont_expr)
  return coeffs

def _replace(expr):
  if isinstance(expr, IndexedOperator):
    return expr.base(*expr.indices)

  elif expr.is_Atom:
    return expr

  else:
    repl_args = [_replace(arg) for arg in expr.args]
    return expr.func(*repl_args)

def _coefficients(expr):

  if isinstance(expr, Operator):
    return expr.coefficients

  elif isinstance(expr, Mul):
    prod_coeffs = list()
    overall_coeff = S.One
    for prod in expr.args:
      if prod.is_number:
        overall_coeff *= prod
        continue

      prod_coeffs.append(_coefficients(prod).items())

    coeffs = dict()
    for prod_coeff in itertools.product(*prod_coeffs):
      prod = list(zip(*prod_coeff))
      term = list(sum(prod[0], tuple()))
      overall_term_coeff = _sort_term(term) * overall_coeff
      if overall_term_coeff == S.Zero:
        continue

      term = tuple(term)
      coeff = reduce(mul, prod[1])
      if overall_term_coeff != S.One:
        coeff *= overall_term_coeff

      if term in coeffs:
        coeffs[term] += coeff
      else:
        coeffs[term] = coeff

    simp_coeffs = dict()
    for k, v in coeffs.items():
      simp_v = expand(v)
      if simp_v != S.Zero:
        simp_coeffs[k] = simp_v

    return simp_coeffs

  elif isinstance(expr, Add):
    coeffs = dict()
    for op_term in expr.args:
      term_coeffs = _coefficients(op_term)
      for term, coeff in term_coeffs.items():
        if term in coeffs:
          coeffs[term] += coeff
        else:
          coeffs[term] = coeff

    simp_coeffs = dict()
    for k, v in coeffs.items():
      simp_v = expand(v)
      if simp_v != S.Zero:
        simp_coeffs[k] = simp_v

    return simp_coeffs

def _sort_term(terms):
  """
  sorts list of terms

  Args:
    terms: a list of GrassmannTerm's to sort

  Returns:
    S.One - if an even number of swaps of anti-commutative objects is performed
    -S.One - if an odd number of swaps of anti-commutative objects is performed
    S.Zero - if two identical anti-commutative objects exist in the list of terms

  Note:
    The list passed in is ordered in place. Unless S.Zero is returned, then
    the terms are not guaranteed to be ordered

  Todo:
    - This uses a bubble sort, which is very slow for anything other than a very small
      list size. However, we do have a very small size in this case, but it might be
      a good idea to improve this in the future

    - Also, this function should be tested.

    - may want to not use __repr__() for the comparison, and instead use the types comparison
      operator.
  """
  stats = S.One
  n = len(terms)
  while n > 1:
    newn = 0
    for i in range(1, n):
      if terms[i-1] > terms[i]:
        terms[i-1], terms[i] = terms[i], terms[i-1]
        newn = i
        if fermionic(terms[i-1].term) and fermionic(terms[i].term):
          stats = -stats

      if terms[i-1] == terms[i] and fermionic(terms[i-1].term) and fermionic(terms[i].term):
        return S.Zero

    n = newn

  return stats


def bosonic(expr):
  return statistics(expr) == 1

def fermionic(expr):
  return statistics(expr) == -1

def statistics(expr):
  """
  TODO:
    Find a more generic way to include GrassmannSymbol and other things...
  """

  if isinstance(expr, (Operator, GrassmannSymbol)):
    return expr.statistics

  elif isinstance(expr, Add):
    statistics_per_term = {statistics(term) for term in expr.args}
    if len(statistics_per_term) != 1:
      raise TypeError("All terms should have same statistics")

    return statistics_per_term.pop()

  elif isinstance(expr, Mul):
    return reduce(mul, [statistics(prod) for prod in expr.args], 1)

  else:
    return 1


def momentum(expr):
  if isinstance(expr, Operator):
    return expr.momentum

  elif isinstance(expr, Add):
    mom_per_term = {momentum(term) for term in expr.args}
    if len(mom_per_term) != 1:
      raise TypeError("All terms should have the same momentum")

    return mom_per_term.pop()

  elif isinstance(expr, Mul):
    return sum([momentum(prod) for prod in expr.args], P0)

  else:
    return P0

def project_momentum(expr, *momenta):
  """
  TODO:
    check that all momenta gets used up?
  """
  if isinstance(expr, Operator):
    return expr.P(*(momenta[0]))
  
  elif isinstance(expr, Add):
    proj_args = (project_momentum(arg, *momenta) for arg in expr.args)
    return expr.func(*proj_args)

  elif isinstance(expr, Mul):
    proj_args = list()
    mom_ind = 0
    for arg in expr.args:
      num_ops = number_of_operators(arg)
      if num_ops == 0:
        proj_args.append(arg)
      else:
        reduced_momenta = momenta[mom_ind:mom_ind+num_ops]
        proj_args.append(project_momentum(arg, *reduced_momenta))
        mom_ind += num_ops

    if mom_ind != len(momenta):
      raise ValueError("Number of operators does not match with number of momenta provided")

    return expr.func(*proj_args)

  else:
    return expr

def number_of_operators(expr):
  """
  TODO:
    add note about how it is the number of operators in each expanded term...
  """
  if isinstance(expr, Operator):
    return 1

  elif isinstance(expr, Add):
    ops_per_term = set([number_of_operators(arg) for arg in expr.args])
    if len(ops_per_term) != 1:
      raise ValueError("All terms should have the same number of Operators")

    return ops_per_term.pop()

  elif isinstance(expr, Mul):
    return sum([number_of_operators(arg) for arg in expr.args])

  else:
    return 0

'''
class MomProjOperator(Operator):

  is_commutative = False

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
  def momentum(self):
    return self._current_momentum

  def rotateTo(self, element):
    self._current_momentum = element * self._original_momentum
    self.operator.rotateTo(element)

  @property
  def coefficients(self):
    if self.momentum == P0:
      return self.operator.coefficients

    coefficients = { GrassmannTerm(k.term, self.momentum):v for k, v in self.operator.coefficients.items() }

    return coefficients

  @property
  def statistics(self):
    return self.operator.statistics
'''
