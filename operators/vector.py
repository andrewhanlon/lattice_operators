import itertools
import collections

from sympy import Array, Indexed, Expr, Symbol, sympify, Tuple
from sympy.core.compatibility import string_types

class Vector(Array):

  is_commutative = False

  def _check_symbolic_index(self, index):
    # Check if any index is symbolic:
    tuple_index = (index if isinstance(index, tuple) else (index,))
    if any([(isinstance(i, Expr) and (not i.is_number)) for i in tuple_index]):
      for i, nth_dim in zip(tuple_index, self.shape):
        if ((i < 0) == True) or ((i >= nth_dim) == True):
          raise ValueError("index out of range")
      return IndexedVector(self, *tuple_index)
    return None


class IndexedVector(Indexed):
  is_commutative = False
