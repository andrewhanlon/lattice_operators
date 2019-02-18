""" grassmann module

This module extends functionality provided by SymPy in order to create
Grassmann variables
"""

import itertools
from collections import defaultdict
from functools import reduce

from sympy import S, Symbol, expand, symbols
from sympy import Mul, Add, Pow, Expr, Sum, Integer
from sympy import Array
from sympy import IndexedBase, Indexed, Idx
from sympy import get_contraction_structure, get_indices, tensorcontraction, tensorproduct, permutedims
import numpy as np

from .vector import Vector
from .cubic_rotations import P0


class GrassmannTerm:
  """
  Todo:
    - Should sorting be handled here?
    - Better name corresponding to the momentum projection?
      I.e. it is not a full term, it is just a term within a term that can have a definite momentum
  """

  def __init__(self, term, momentum=P0):
    self._term = term
    self._momentum = momentum

  @property
  def term(self):
    return self._term

  @property
  def momentum(self):
    return self._momentum

  def _cmp(self):
    return (self.momentum.__repr__(), self.term.__repr__())

  def __hash__(self):
    return hash(self._cmp())

  def __repr__(self):
    return "({}:{})".format(self.momentum.__str__(), self.term.__str__())

  def __str__(self):
    return "({}, {})".format(self.momentum.__str__(), self.term.__str__())

  def __eq__(self, other):
    if isinstance(other, self.__class__):
      return self._cmp() == other._cmp()
    return NotImplemented

  def __ne__(self, other):
    return not self.__eq__(other)

  def __lt__(self, other):
    if isinstance(other, self.__class__):
      return self._cmp() < other._smp()
    return NotImplemented

  def __gt__(self, other):
    if isinstance(other, self.__class__):
      return self._cmp() > other._cmp()
    return NotImplemented

  def __le__(self, other):
    if isinstance(other, self.__class__):
      return self._cmp() <= other._cmp()
    return NotImplemented

  def __ge__(self, other):
    if isinstance(other, self.__class__):
      return self._cmp() >= other._cmp()
    return NotImplemented



class GrassmannSymbol(Symbol):
  """GrassmannSymbol class

  This class extends the functionality of the Symbol class in sympy.
  Its purpose is to create a grassmann symbol.
  """

  statistics = -1

  def __new__(cls, *args, **kwargs):
    """the __new__ method constructs the GrassmannSymbol

    The args passed get passed directly to the Symbol __new__ method
    while also forcing commutative to false.
    """
    return super().__new__(cls, *args, **kwargs, commutative=False)

  def __mul__(self, other):
    if isinstance(other, GrassmannSymbol):
      if other == self:
        return S.Zero
      elif other.name < self.name:
        return -Symbol.__mul__(other,self)

    return super().__mul__(other)

  def __pow__(self, exponent):
    if exponent == 0:
      return S.One
    elif exponent == 1:
      return self
    else:
      return S.Zero


class GrassmannVector(Vector):

  is_commutative = False

  def __new__(cls, iterable=None, shape=None, name=None):
    if name is not None:
      if shape is None:
        return GrassmannSymbol(name)
      else:
        sym_str = "{}".format(name)
        for ind in shape:
          sym_str += ":{}".format(ind)

        return super().__new__(cls, symbols(sym_str, cls=GrassmannSymbol), shape)

    return super().__new__(cls, iterable, shape)


  def transformRight(self, matrix, index=0):
    return self.transform(matrix, index)

  def transformLeft(self, matrix, index=0):
    return self.transform(matrix.T, index)


  def transform(self, matrix, index=0):

    transformed = tensorcontraction(tensorproduct(matrix, self), (1, index+2))
    if index:
      perms = list(range(len(self.shape)))
      perms[0] = index
      perms[index] = 0
      transformed = permutedims(transformed, perms)

    return self.__class__(transformed)


def grassmann_simplify(expr):
  outer_indices, _ = get_indices(expr)
  outer_indices = {index for index in outer_indices if isinstance(index, Idx)}

  if outer_indices:
    raise TypeError("Free indices not yet supported")

  cont_expr = contract_indices(expr)
  simp_expr = expand(cont_expr)
  simp_expr = noncommutative_sort(simp_expr)
  return simp_expr

def contract_indices(expr):
  """
  TODO:
    Test this?
  """
  contractions = get_contraction_structure(expand(expr))
  if len(contractions.keys()) == 1 and None in contractions:
    return expr

  summation = S.Zero
  for indices in contractions:
    if isinstance(indices, Expr):
      continue

    partial_sum = S.Zero
    for term in contractions[indices]:
      partial_sum += term

    indices = {index for index in indices if isinstance(index, Idx)}
    summation += Sum(partial_sum, *indices).doit()

  return summation

def noncommutative_sort(expr):
  if expr.is_Atom:
    return expr
  else:
    sorted_args = (noncommutative_sort(arg) for arg in expr.args)
    if isinstance(expr, Mul):
      return _noncommutative_sort_product(Mul(*sorted_args))
    else:
      return expr.func(*sorted_args)

def _noncommutative_sort_product(product):
  """
  sorts a product of grassmann variables, paying attention to the anti-commutative
  nature of grassmann variables.

  Args:
    product: a product containing grassmann variables

  Returns:
    The product in canonical ordering

  Todo:
    This uses a bubble sort, which is very slow for anything other than a very small
    list size. However, we do have a very small size in this case, but it might be
    a good idea to improve this in the future
  """

  while True:
    if not isinstance(product, Mul):
      return product

    arglist = list(product.args)
    i = 0
    while i < len(arglist)-1:
      slice_prod = arglist[i]*arglist[i+1]
      is_mul = isinstance(slice_prod, Mul)
      arglist[i:i+2] = slice_prod.args if is_mul else [slice_prod]
      i += 1

    new_product = Mul(*arglist)
    if product == new_product:
      return new_product
    product = new_product


def grassmann_coefficients(expr):
  """
  TODO:
    Test this?
  """
  coeffs_dict = defaultdict(int)
  for term, coeff in expr.as_coefficients_dict().items():
    if term.is_Atom:
      coeffs_dict.update({term: coeff})
      continue
    products = [product for product in term.args if not product.is_complex]
    extra_coeffs = [coeffs for coeffs in term.args if coeffs.is_complex]
    coeffs = Integer(1)
    if extra_coeffs:
      coeffs = reduce(lambda x, y: x*y, extra_coeffs)
    coeffs *= coeff
    term = term.func(*products)
    if term in coeffs_dict:
      coeffs += coeffs_dict[term]

    coeffs_dict.update({term: coeffs})
  return coeffs_dict
