import itertools
from sortedcontainers import SortedSet
from collections import OrderedDict
from collections import defaultdict

from sympy import zeros, Expr, Idx, get_indices, Array, conjugate, Matrix, Transpose, IndexedBase, trace, Indexed
from sympy import eye, Integer, S
from sympy import simplify
from sympy import Add
from sympy import Sum

from .cubic_rotations import spinor_representation, LittleGroup, E
from .cubic_rotations import _GENERATORS as GENERATORS
from .cubic_rotations import P0
from .grassmann import grassmann_simplify, coefficients, GrassmannField
from .tensors import GammaRep

# @ADH - May want to not allow for contracted operators?
#        Maybe by requiring no grassmann variables, just grassmann fields
#        indexed.

# @ADH - Possible improvement: An operator with free_indices
#        will return multiple operators when a call is made
#        to sub_free_indices, but all of these returned operators
#        will have the same contraction structure, which makes it
#        wasteful to recompute the contractions. But, I want to have
#        each subbed operator as a separtge Operator object.


class QuarkField(GrassmannField):

  @classmethod
  def create(cls, name):
    return super().__new__(cls, shape=(3,4), name=name)

  def rotate(self, element):
    return self.transformRight(spinor_representation.rotation(element, False), 1)

  def colorRotate(self, matrix):
    return self.transformRight(matrix, 0)


class AntiQuarkField(GrassmannField):

  @classmethod
  def create(cls, name):
    return super().__new__(cls, shape=(3,4), name="{}+".format(name))

  def rotate(self, element):
    return self.transformLeft(spinor_representation.rotation(element, True), 1)

  def colorRotate(self, matrix):
    return self.transformLeft(matrix, 0)


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




class OperatorRepresentation:
  
  def __init__(self, *operators):

    if not operators:
      raise ValueError("Must provide at least one basis operator")

    self._basis = OperatorBasis(*operators)
    self._little_group = LittleGroup(self.basis.bosonic, self.basis.momentum)

    self._rep_matrices = dict()
    self._characters = dict()

  @property
  def momentum(self):
    return self.basis.momentum

  @property
  def bosonic(self):
    return self.basis.bosonic

  @property
  def fermionic(self):
    return self.basis.fermionic
  
  @property
  def dimension(self):
    return self.basis.dimension

  @property
  def basis(self):
    return self._basis

  @property
  def little_group(self):
    return self._little_group

  def littleGroupContents(self, nice=False, use_generators=True):
    contents = dict()
    for lgIrrep in self.little_group.irreps:
      #print(lgIrrep)
      occurences = S.Zero
      for element in self.little_group.elements:
        occurences += self.getCharacter(element, use_generators) * conjugate(self.little_group.getCharacter(lgIrrep, element))

      occurences = simplify(occurences)   # @ADH - I don't like that this was necessary
      occurences /= self.little_group.order
      contents[lgIrrep] = int(occurences)
      if occurences != contents[lgIrrep]:
        raise ValueError("Occurence is not an integer")

    if nice:
      nice_str = ""
      for irrep, occurences in contents.items():
        if occurences == 1:
          nice_str += "{} + ".format(irrep)
        elif occurences > 1:
          nice_str += "{} {} + ".format(occurences, irrep)

      nice_str = nice_str[:-3]
      return nice_str

    return contents

  def getCharacter(self, lg_element, use_generators=True):

    if lg_element in self._characters:
      return self._characters[lg_element]
    elif lg_element in self._rep_matrices:
      char = trace(self._rep_matrices[lg_element])
      self._characters[lg_element] = char
      return char
    else:
      char = trace(self.getRepresentationMatrix(lg_element, use_generators))
      self._characters[lg_element] = char
      return char

  def getRepresentationMatrix(self, lg_element, use_generators=True):

    if lg_element in self._rep_matrices:
      return self._rep_matrices[lg_element]
    elif use_generators and GENERATORS[lg_element] == "invert":
      self._rep_matrices[lg_element] = self.getRepresentationMatrix(lg_element.inverse(), True).inv()
      return self._rep_matrices[lg_element]
    elif use_generators and GENERATORS[lg_element]:
      rep_mat = Integer(1)
      for el in GENERATORS[lg_element]:
        rep_mat *= self.getRepresentationMatrix(el, True)

      self._rep_matrices[lg_element] = rep_mat
      return self._rep_matrices[lg_element]
    else:
      self._compute_rep_matrix(lg_element)
      return self._rep_matrices[lg_element]
      

    return self._rep_matrices[lg_element]

  def _compute_rep_matrix(self, lg_element):
    if lg_element == E:
      self._rep_matrices[lg_element] = eye(self.dimension)
      self._characters[lg_element] = self.dimension
    else:
      transformed_basis = self.basis.rotate(lg_element)
      self._rep_matrices[lg_element] = self.basis.matrix.pinv() * transformed_basis.matrix


  def irreducible(self, use_generators=True):
    inner_product = S.Zero
    for element in self.little_group.elements:
      inner_product += conjugate(self.getCharacter(element, use_generators)) * self.getCharacter(element, use_generators)

    inner_product = simplify(inner_product)   # @ADH - I don't like that this was necessary
    if inner_product == self.little_group.order:
      return True

    return False


class OperatorBasis:

  def __init__(self, *in_operators, grassmann_basis=None):

    if not in_operators:
      raise ValueError("Must provide at least one basis operator")

    same_type = all(in_op.bosonic == in_operators[0].bosonic for in_op in in_operators)
    same_mom = all(in_op.momentum == in_operators[0].momentum for in_op in in_operators)

    if not same_type or not same_mom:
      raise ValueError("All basis vectors must have same total momentum and be of same particle type")

    self._operators = in_operators
    self._bosonic = in_operators[0].bosonic
    self._fermionic = in_operators[0].fermionic
    self._momentum = in_operators[0].momentum

    self._grassmann_basis = grassmann_basis
    if grassmann_basis is None:
      self._create_grassmann_basis(in_operators)

    self._vectors = dict()
    self._matrix = None

  @property
  def momentum(self):
    return self._momentum

  @property
  def operators(self):
    return self._operators

  def _create_grassmann_basis(self, in_operators):

    terms = list()
    for in_operator in in_operators:
      terms.extend(list(in_operator.getTerms()))

    self._grassmann_basis = set(terms)

  def vector(self, operator):
    if operator not in self._vectors:
      _vector = list()
      coeffs_dict = operator.coefficients
      for grassmann_vector in self.grassmann_basis:
        if grassmann_vector in coeffs_dict:
          _vector.append(coeffs_dict[grassmann_vector])
          del coeffs_dict[grassmann_vector]
        else:
          _vector.append(0)

      if coeffs_dict:
        raise ValueError("Basis is not complete")

      self._vectors[operator] = _vector

    return self._vectors[operator]

  @property
  def matrix(self):
    if self._matrix is None:
      mat = list()
      for operator in self.operators:
        mat.append(self.vector(operator))

      self._matrix = Matrix(mat).T

    return self._matrix

  @property
  def bosonic(self):
    return self._bosonic

  @property
  def fermionic(self):
    return self._fermionic

  @property
  def grassmann_basis(self):
    return self._grassmann_basis

  @property
  def dimension(self):
    return len(self.operators)

  @property
  def operators(self):
    try:
      return self._operators
    except AttributeError:
      self._operators = SortedSet()
      return self._operators

  def rotate(self, lg_element):
    transformed_ops = list()
    for operator in self.operators:
      transformed_ops.append(operator.rotate(lg_element))

    return OperatorBasis(*transformed_ops, grassmann_basis=self.grassmann_basis)



# @ADH - In the future, I'd like support for indexed Operators
class Operator:

  # @ADH - I may want to check for any GrassmannSymbols in operator and raise a warning
  #        because this may correspond to a contracted object...
  def __init__(self, operator, momentum=None):
    outer_indices, index_syms = get_indices(operator)

    # @ADH - Again, I don't really like this hack of removing non-Idx indices
    outer_indices = {index for index in outer_indices if isinstance(index, Idx)}

    if outer_indices:
      raise TypeError("Operators with free indices not yet supported")

    self._operator = operator
    self._momentum = momentum

    self._simplified = None
    self._coefficients = None
    self._terms = None

  @property
  def momentum(self):
    return self._momentum

  def getTerms(self):
    if self._terms is None:
      self._terms = set(self.coefficients.keys())

    return self._terms

  @property
  def coefficients(self):
    if self._coefficients is None:
      self._coefficients = defaultdict(int)
      for grassmann_term, coeff in coefficients(self.simplified).items():
        self._coefficients[(self.momentum, grassmann_term)] = coeff

    return self._coefficients

  @property
  def zero(self):
    return self.simplified == S.Zero

  def projectMomentum(self, momentum):
    return self.__class__(self.operator, momentum)

  @staticmethod
  def _rotate(element, expr):
    if isinstance(expr, Indexed) and (isinstance(expr.base, QuarkField)
                                      or isinstance(expr.base, AntiQuarkField)):
      transformed_base = expr.base.rotate(element)
      return transformed_base[expr.indices]

    elif isinstance(expr, QuarkField) or isinstance(expr, AntiQuarkField):  # @ADH - Should this be allowed?
      return expr.rotate(element)

    elif expr.is_Atom:
      return expr

    else:
      args = [Operator._rotate(element, arg) for arg in expr.args]
      return expr.func(*args)

  def rotate(self, element):
    if self.momentum is None:
      trans_momentum = None
    else:
      trans_momentum = element * self.momentum

    trans_operator = Operator._rotate(element, self.operator)

    return self.__class__(trans_operator, trans_momentum)

  @staticmethod
  def _numberOfQuarks(expr):
    if isinstance(expr, Indexed) and (isinstance(expr.base, QuarkField)
                                      or isinstance(expr.base, AntiQuarkField)):
      return 1

    elif isinstance(expr, QuarkField) or isinstance(expr, AntiQuarkField):  # @ADH - Should this be allowed?
      return 1

    elif isinstance(expr, Add):
      quarks_per_term = set([Operator._numberOfQuarks(arg) for arg in expr.args])
      if len(quarks_per_term) != 1:
        raise ValueError("All terms should have the same number of Quarks")

      return quarks_per_term.pop()
    elif expr.is_Atom:
      return 0
    else:
      return sum([Operator._numberOfQuarks(arg) for arg in expr.args])

  @property
  def number_of_quarks(self):
    try:
      return self._number_of_quarks
    except AttributeError:
      self._number_of_quarks = Operator._numberOfQuarks(self.operator)
      return self._number_of_quarks


  @property
  def bosonic(self):
    return self.number_of_quarks % 2 == 0

  @property
  def fermionic(self):
    return self.number_of_quarks % 2 == 1

  @property
  def operator(self):
    return self._operator

  @property
  def simplified(self):
    if self._simplified is None:
      self._simplified = grassmann_simplify(self.operator, True)

    return self._simplified

  def __str__(self):
    if self.momentum is None:
      return self.operator.__str__()

    return "{}[{}]".format(self.operator.__str__(), self.momentum.__repr__())

  def __repr__(self):
    if self.momentum is None:
      return self.simplified.__repr__()

    return "{}_{}".format(self.simplified.__repr__(), self.momentum.__repr__())

  def __hash__(self):
    return hash(self.__repr__())

  def __eq__(self, other):
    if isinstance(other, self.__class__):
      return self.__repr__() == other.__repr__()
    return NotImplemented

  def __ne__(self, other):
    if isinstance(other, self.__class__):
      return not self.__eq__(other)
    return NotImplemented

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

  def __mul__(self, other):
    if self.zero:
      return S.Zero

    elif isinstance(other, self.__class__):
      if other.zero:
        return S.Zero

      return OperatorMul(self, other)

    # @ADH - Assumes all terms in OperatorMul are of type Operator (enforced by OperatorMul itself)
    elif isinstance(other, OperatorMul):
      return OperatorMul(self, *other.operators)

    elif isinstance(other, OperatorAdd):
      added_ops = list()
      for operator in other.operators:
        added_ops.append(self * operator)

      return OperatorAdd(*added_ops)

    else:
      op = self.operator * other
      return Operator(op, self.momentum)

  def __rmul__(self, other):
    op = other * self.operator
    return Operator(op, self.momentum)


  def __add__(self, other):
    if self.zero:
      return other

    elif isinstance(other, self.__class__) and self.momentum == other.momentum and self.number_of_quarks == other.number_of_quarks:
      if other.zero:
        return self

      return Operator(self.operator + other.operator, self.momentum)

    elif other == S.Zero:
      return self

    raise TypeError("Can only add another Operator to an Operator with same momentum and number of quarks")

  def __radd__(self, other):
    return self.__add__(other)

  def __sub__(self, other):
    return self.__add__(-other)

  def __neg__(self):
    return Operator(-self.operator, self.momentum)

'''
# @ADH - WARNING WARNING WARNING - ASSUMES products are of Two Baryons!!!
#      - FIX ME LATER!!!
class GrassmannProductBasis:
  
  def __init__(self, product_terms):
    self._terms = set()
    self._original = defaultdict(int)
    self._parity = defaultdict(int)

    for product_term in product_terms:
      if len(product_term) > 2:
        raise ValueError("three and higher particle operators not currently supported")

      # @ADH - ASSUMES Two Baryons!!
      if product_term[0] == product_term[1]:
        continue

      str_rep1 = "{}__{}".format(product_term[0][0].__repr__(), product_term[0][1].__repr__())
      str_rep2 = "{}__{}".format(product_term[1][0].__repr__(), product_term[1][1].__repr__())

      if str_rep1 < str_rep2:
        new_term = (product_term[1], product_term[0])
        self._parity[new_term] = -S.One
        self._original[new_term] = product_term

        if new_term in self._terms:
          raise ValueError('shit1')

        self._terms.add(new_term)
      else:
        self._parity[product_term] = S.One
        self._original[product_term] = product_term

        if product_term in self._terms:
          raise ValueError('shit2')

        self._terms.add(product_term)

  @property
  def terms(self):
    return self._terms

  def original(self, term):
    return self._original[term]

  def parity(self, term):
    return self._parity[term]
'''


class OperatorMul:

  def __init__(self, *operators):

    if len(operators) > 2:
      raise ValueError("3-particle and higher operators not currently supported")

    all_operators = all(isinstance(op, Operator) for op in operators)
    if not all_operators:
      raise ValueError("All objects passed to OperatorMul must be of type Operator")

    self._operators = operators

    self._coefficients = None

    self._raw_terms = None
    self._terms = None

  def __new__(self, *operators):
    if not operators:
      return S.Zero

    elif len(operators) == 1:
      return operators[0]

    # @ADH - Make sure all operators are of type Operator?

    return object.__new__(self)

  '''
  @property
  def basis(self):
    if self._basis is None:
      term_list = list()
      for op in self.operators:
        term_list.append(op.getTerms())

      self._basis = GrassmannProductBasis(set(itertools.product(*term_list)))

    return self._basis

  def getTerms(self):
    return self.basis.terms

  @property
  def coefficients(self):
    if self._coefficients is None:
      terms = self.getTerms()

      op_coeffs = list()
      for op in self.operators:
        op_coeffs.append(op.coefficients)

      coeffs_dict = defaultdict(int)
      for term in terms:
        coeffs_dict[term] = self.basis.parity(term)
        for op_term, op_coeff in zip(self.basis.original(term), op_coeffs):
          coeffs_dict[term] *= op_coeff[op_term]

      self._coefficients = coeffs_dict

    return self._coefficients
  '''

  
  @property
  def raw_terms(self):
    if self._raw_terms is None:
      term_list = list()
      for op in self.operators:
        term_list.append(op.getTerms())

      self._raw_terms = set(itertools.product(*term_list))

    return self._raw_terms

  

  # @ADH - ONLY WORKS ASSUMING TWO-BARYON OPERATORS
  @property
  def coefficients(self):
    if self._coefficients is None:
      op_coeffs = list()
      for op in self.operators:
        op_coeffs.append(op.coefficients)

      coeffs_dict = defaultdict(int)

      for term in self.raw_terms:
        if term[0] == term[1]:
          continue

        str_rep1 = "{}__{}".format(term[0][0].__repr__(), term[0][1].__repr__())
        str_rep2 = "{}__{}".format(term[1][0].__repr__(), term[1][1].__repr__())
        
        coeff = op_coeffs[0][term[0]] * op_coeffs[1][term[1]]

        new_term = term
        if str_rep1 < str_rep2:
          coeff = -coeff
          new_term = tuple([term[1], term[0]])

        if new_term in coeffs_dict:
          coeffs_dict[new_term] += coeff
        else:
          coeffs_dict[new_term] = coeff

      self._coefficients = { k:v for k, v in coeffs_dict.items() if v }

    return self._coefficients


  def getTerms(self):
    if self._terms is None:
      self._terms = set(self.coefficients.keys())

    return self._terms

        
  @property
  def bosonic(self):
    return self.number_of_quarks % 2 == 0

  @property
  def fermionic(self):
    return self.number_of_quarks % 2 == 1

  @property
  def operators(self):
    return self._operators

  @property
  def momenta(self):
    mom_list = list()
    for operator in self.operators:
      mom_list.append(operator.momentum)

    return tuple(mom_list)

  # @ADH - Check operators and momenta are same size?
  def projectMomentum(self, *momenta):
    operators = list()
    for operator, momentum in zip(self.operators, momenta):
      operators.append(Operator(operator.operator, momentum))

    return self.__class__(*operators)

  @property
  def number_of_operators(self):
    return len(self.operators)

  @property
  def number_of_quarks(self):
    return sum(op.number_of_quarks for op in self.operators)

  # @ADH - Could replace this with a sum function that doesn't start with zero
  @property
  def momentum(self):
    tot_mom = P0
    for operator in self.operators:
      if operator.momentum is not None:
        tot_mom += operator.momentum

    return tot_mom

  def rotate(self, element):
    trans_operators = list()
    for operator in self.operators:
      trans_operators.append(operator.rotate(element))

    return OperatorMul(*trans_operators)


  # @ADH - Anything else need to be the same between other and self?
  #        We can't add an Operator, because it has a different number of operators.
  #        Thus, I should probably check for the same number of operators (equivalently, the number of momentum projections).
  def __add__(self, other):
    if isinstance(other, self.__class__) and self.momentum == other.momentum and self.number_of_quarks == other.number_of_quarks:
      return OperatorAdd(self, other)

    elif isinstance(other, OperatorAdd) and self.momentum == other.momentum and self.number_of_quarks == other.number_of_quarks:
      return OperatorAdd(other, *other.operators)

    elif other == S.Zero:
      return self

    return NotImplemented

  def __radd__(self, other):
    return self.__add__(other)


  def __mul__(self, other):
    if isinstance(other, self.__class__):
      return OperatorMul(*self.operators, *other.operators)

    elif isinstance(other, Operator):
      if other.zero:
        return S.Zero

      return OperatorMul(*self.operators, other)

    elif isinstance(other, OperatorAdd):
      op_adds = list()
      for op in other.operators:
        op_adds.append(self*op)

      return OperatorAdd(*op_adds)

    elif other == S.Zero:
      return S.Zero

    operators = list(self.operators)
    operators[0] = operators[0] * other
    return OperaotrMul(*operators)

  def __rmul__(self, other):
    operators = list(self.operators)
    operators[0] = other * operators[0]
    return OperatorMul(*operators)

  def __sub__(self, other):
    return self.__add__(-other)

  def __neg__(self):
    operators = list(self.operators)
    operators[0] = -operators[0]
    return OperatorMul(*operators)


  def __str__(self):
    op_str = ""
    for op in self.operators:
      op_str += "({}) ".format(op.__str__())

    return op_str[:-1]

  def __repr__(self):
    op_repr = ""
    for op in self.operators:
      op_repr += "{}_*_".format(op.__repr__())

    return op_repr[:-3]

  def __hash__(self):
    return hash(self.__repr__())


  # @ADH - Do I want to include comparisons with OperatorAdd or Operator?

  def __eq__(self, other):
    if isinstance(other, self.__class__):
      return self.__repr__() == other.__repr__()
    return NotImplemented

  def __ne__(self, other):
    if isinstance(other, self.__class__):
      return not self.__eq__(other)
    return NotImplemented

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



class OperatorAdd:

  # @ADH - Add more checkes?
  def __init__(self, *operators):
    same_type = all(op.bosonic == operators[0].bosonic for op in operators)
    if not same_type:
      raise ValueError("Can not add a fermionic operator to a bosonic one")

    same_mom = all(op.momentum == operators[0].momentum for op in operators)
    if not same_mom:
      raise ValueError("Can not add operators with different total momentum")

    self._operators = SortedSet(operators)
    self._bosonic = operators[0].bosonic
    self._fermionic = operators[0].fermionic
    self._momentum = operators[0].momentum

    self._coefficients = None
    self._terms = None

  def __new__(self, *operators):
    if not operators:
      return S.Zero

    if len(operators) == 1:
      return operators[0]

    return object.__new__(self)

  @property
  def operators(self):
    return self._operators

  @property
  def momentum(self):
    return self._momentum

  @property
  def bosonic(self):
    return self._bosonic

  @property
  def fermionic(self):
    return self._fermionic

  def rotate(self, element):
    trans_operators = list()
    for operator in self.operators:
      trans_operators.append(operator.rotate(element))

    return OperatorAdd(*trans_operators)

  def projectMomentum(self, *momenta):
    operators = list()
    for operator in self.operators:
      operators.append(operator.projectMomentum(*momenta))

    return self.__class__(*operators)

  def getTerms(self):
    if self._terms is None:
      terms_list = list()
      for op in self.operators:
        for term in op.getTerms():
          terms_list.append(term)

      self._terms = set(terms_list)

    return self._terms

  @property
  def coefficients(self):
    if self._coefficients is None:
      coeffs = defaultdict(int)
      for term in self.getTerms():
        coeffs[term] = S.Zero
        for operator in self.operators:
          if term in operator.coefficients:
            coeffs[term] += operator.coefficients[term]

      self._coefficients = { k:v for k, v in coeffs.items() if v }

    return self._coefficients


  # @ADH - add more checks?
  def __add__(self, other):
    if isinstance(other, Operator) and self.momentum == other.momentum and self.bosonic == other.bosonic:
      return OperatorAdd(*self.operators, other)

    elif isinstance(other, self.__class__) and self.momentum == other.momentum and self.bosonic == other.bosonic:
      return OperatorAdd(*self.operators, *other.operators)

    elif isinstance(other, OperatorMul) and self.momentum == other.momentum and self.bosonic == other.bosonic:
      return OperatorAdd(*self.operators, other)

    elif other == S.Zero:
      return self

    return NotImplemented

  def __radd__(self, other):
    return self.__add__(other)


  def __mul__(self, other):
    add_ops = list()
    for op in self.operators:
      add_ops.append(op*other)

    return OperatorAdd(*add_ops)

  def __rmul__(self, other):
    add_ops = list()
    for op in self.operators:
      add_ops.append(other*op)

    return OperatorAdd(*add_ops)

  def __sub__(self, other):
    return self.__add__(-other)

  def __neg__(self):
    ops = list()
    for op in self.operators:
      ops.append(-op)

    return OperatorAdd(*ops)


  def __str__(self):
    op_str = ""
    for op in self.operators:
      op_str += "{} + ".format(op.__str__())

    return op_str[:-3]

  def __repr__(self):
    op_repr = ""
    for op in self.operators:
      op_repr += "{}_+_".format(op.__repr__())

    return op_repr[:-3]

  def __hash__(self):
    return hash(self.__repr__())

  # @ADH - Do I want to include comparisons with OperatorAdd or Operator?

  def __eq__(self, other):
    if isinstance(other, self.__class__):
      return self.__repr__() == other.__repr__()
    return NotImplemented

  def __ne__(self, other):
    if isinstance(other, self.__class__):
      return not self.__eq__(other)
    return NotImplemented

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

