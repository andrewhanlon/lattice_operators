"""
This module extends functionality provided by SymPy in order to create
Grassmann variables
"""

import itertools
from collections import defaultdict
from functools import reduce

from sympy.core.symbol import Symbol
from sympy.core.singleton import S
from sympy.core.mul import Mul
from sympy import Add
from sympy.core.power import Pow
from sympy.core.expr import Expr
from sympy.core.function import expand
from sympy.tensor.array import ImmutableDenseNDimArray
from sympy.tensor.array import Array
from sympy.tensor.indexed import IndexedBase, Indexed, Idx
from sympy.tensor.index_methods import get_contraction_structure, get_indices
from sympy.concrete.summations import Sum
from sympy.tensor.array import tensorcontraction, tensorproduct, permutedims
from sympy import Integer
import numpy as np


class GrassmannSymbol(Symbol):
  def __new__(cls, *args, **kwargs):
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


class GrassmannField(ImmutableDenseNDimArray):

  is_commutative = False # @ADH - Is this necessary?

  def __new__(cls, iterable=None, shape=None, name=None):
    if name is not None: # create new object

      # @ADH - Fix this - I admit this is a strangeway to make this (i.e. making a numpy array that we then
      # use to make a Sympy array. But, I wasn't quite sure how else to do it, and this was
      # simple enough.

      grassmann_vars = np.empty(shape=shape, dtype=object)
      it = np.nditer(grassmann_vars, flags=["refs_ok", "multi_index"], op_flags=["writeonly"])
      itnext = not it.finished
      while itnext:
        multi_index = list(str(x) for x in it.multi_index)
        var_name = "{name}_{{{indices}}}".format(name=name, indices=','.join(multi_index))
        it[0] = GrassmannSymbol(var_name)
        itnext = it.iternext()

      return super().__new__(cls, grassmann_vars)

    return super().__new__(cls, iterable, shape)

  
  # @ADH - Maybe you can just call the super() and just change Indexed to GrassmannIndexed?
  def _check_symbolic_index(self, index):
    # Check if any index is symbolic:
    tuple_index = (index if isinstance(index, tuple) else (index,))
    if any([(isinstance(i, Expr) and (not i.is_number)) for i in tuple_index]):
      for i, nth_dim in zip(tuple_index, self.shape):
        if ((i < 0) == True) or ((i >= nth_dim) == True):
          raise ValueError("index out of range")
      return GrassmannIndexed(self, *tuple_index)
    return None

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


class GrassmannIndexed(Indexed):
  is_commutative = False



# @ADH - Do you think this is the fastest way to do this?
# @ADH - Unfortunately, if I do the contractions first (which seems wasteful not to),
#        then this fails due to 'Indices are not consistent'...
def grassmann_simplify(expr, contract=True):
  # @ADH - extremely annoying that get_indices returns a different
  #        ordering for the outer_indices on the same exact object
  #        on different runs
  #
  #        We should think of a way to remove this issue
  outer_indices, index_syms = get_indices(expr)

  # @ADH - Like in perform_contractions. I don't like removing the indices like this
  outer_indices = {index for index in outer_indices if isinstance(index, Idx)}

  if outer_indices:

    # Build list (all_subs) of all possible substitutions to make for
    # the free indices
    ranges = list()
    shape = list()
    for outer_index in outer_indices:
      ranges.append(range(outer_index.lower, outer_index.upper+1))
      shape.append(outer_index.upper - outer_index.lower + 1)

    shape = tuple(shape)

    all_subs = list()
    combinations = itertools.product(*ranges)
    for combination in combinations:
      subs = list()
      for sub in zip(outer_indices, combination):
        subs.append(sub)
      all_subs.append(subs)

    exprs = list()
    for subs in all_subs:
      subed_expr = expr.subs(subs)
      if contract:
        subed_expr = perform_contractions(subed_expr)

      subed_expr = expand(subed_expr)    # @ADH - Should expand be performed above?
      subed_expr = _grassmann_simplify(subed_expr)
      exprs.append(subed_expr)

    return Array(exprs, shape=shape)

  elif contract:
    cont_expr = perform_contractions(expr)
    simp_expr = expand(cont_expr)
    simp_expr = _grassmann_simplify(simp_expr)
    return simp_expr

  return expr
  

# @ADH - Make sure this is right, too
def perform_contractions(expr):
  contractions = get_contraction_structure(expand(expr))
  if len(contractions.keys()) == 1 and None in contractions:
    return expr

  summation = S.Zero
  for indices in contractions:
    partial_sum = S.Zero
    # @ADH - What does this if statement do?
    if isinstance(indices, Expr):
      continue
    for term in contractions[indices]:
      partial_sum += term

    # @ADH - I don't like this exactly, but it works for now
    #        I.e. making sure indices are of type Idx
    #        I don't like it, simply because I wish the get_contraction_structure
    #        wouldn't return non-indices (like integers)
    indices = [index for index in indices if isinstance(index, Idx)]
    summation += Sum(partial_sum, *indices).doit()

  return summation

'''
def grassmann_expand(expr):
  #cont_expr = perform_contractions(expr)
  array = make_array(expr, True)
  return array
'''

  
# @ADH - Need to make sure this function is right and deals with powers appropriately
def _grassmann_simplify(expr):

  if expr.is_Atom:
    return expr
  else:
    simplified_args = (_grassmann_simplify(arg) for arg in expr.args)
    if isinstance(expr, Mul):
      return _simplify_product(Mul(*simplified_args))
    else:
      return expr.func(*simplified_args) 


# @ADH - Need to make sure this function is right and deals with powers appropriately
def _simplify_product(product):
  while True:
    if not isinstance(product, Mul):
      return product

    arglist = list(product.args)
    i = 0
    while i < len(arglist)-1:
      slice_prod = arglist[i]*arglist[i+1]
      is_mul = isinstance(slice_prod,Mul)
      arglist[i:i+2] = slice_prod.args if is_mul else [slice_prod]
      i += 1

    new_product = Mul(*arglist)
    if product == new_product:
      return new_product
    product = new_product


# @ADH - This function sorta scares me. I would really like to make sure it is right
def coefficients(expr):
  coeffs_dict = defaultdict(int)
  for term, coeff in expr.as_coefficients_dict().items():
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

