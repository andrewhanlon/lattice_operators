from sympy import expand, S, Integer, trace, conjugate, Matrix, eye

from .operators import bosonic, fermionic, momentum, terms, coefficients, rotate_to
from .cubic_rotations import LittleGroup, E, P0
from .cubic_rotations import _GENERATORS as GENERATORS # @ADH - Maybe change the way this is named?

class OperatorRepresentationError(Exception):
  pass

class OperatorRepresentation:
  
  def __init__(self, *operators, use_generators=False):

    if not operators:
      raise ValueError("Must provide at least one basis operator")

    self._use_generators = use_generators
    self._basis = OperatorBasis(*operators)
    if self.basis.momentum is not None:
      self._little_group = LittleGroup(self.basis.bosonic, self.basis.momentum)

    full_rank = min(self.basis.matrix(E).cols, self.basis.matrix(E).rows)
    if self.basis.matrix(E).rank() < full_rank:
      raise OperatorRepresentationError("Rank defecient: Operators not linearly independent?")

    self._rep_matrices = dict()
    self._characters = dict()

  @property
  def basis(self):
    return self._basis

  @property
  def little_group(self):
    try:
      return self._little_group
    except AttributeError:
      return None

  def lgIrrepOccurences(self, nice=False):
    occurences = dict()
    for lgIrrep in self.little_group.irreps:
      occurence = S.Zero
      for element in self.little_group.elements:
        occurence += self.getCharacter(element) * conjugate(self.little_group.getCharacter(lgIrrep, element))

      occurence = expand(occurence)
      occurence /= self.little_group.order
      occurences[lgIrrep] = occurence

    if nice:
      nice_str = ""
      for irrep, occurence in occurences.items():
        if occurence == 1:
          nice_str += "{} + ".format(irrep)
        elif occurence > 0:
          nice_str += "{} {} + ".format(occurence, irrep)

      nice_str = nice_str[:-3]
      return nice_str

    return occurences

  def getCharacter(self, lg_element):

    if lg_element in self._characters:
      return self._characters[lg_element]
    elif lg_element in self._rep_matrices:
      char = trace(self._rep_matrices[lg_element])
      self._characters[lg_element] = expand(char)
      return char
    else:
      char = trace(self.getRepresentationMatrix(lg_element))
      self._characters[lg_element] = expand(char)
      return char

  def getRepresentationMatrix(self, lg_element):
    use_generators = (self.basis.momentum == P0) and self._use_generators

    # @ADH - tidy this up
    if lg_element in self._rep_matrices:
      return self._rep_matrices[lg_element]
    elif use_generators and GENERATORS[lg_element] == "invert":
      self._rep_matrices[lg_element] = self.getRepresentationMatrix(lg_element.inverse()).inv()
    elif use_generators and GENERATORS[lg_element]:
      rep_mat = Integer(1)
      for el in GENERATORS[lg_element]:
        rep_mat *= self.getRepresentationMatrix(el)

      self._rep_matrices[lg_element] = rep_mat
    else:
      self._compute_rep_matrix(lg_element)

    return self._rep_matrices[lg_element]

  def _compute_rep_matrix(self, lg_element):
    if lg_element == E:
      self._rep_matrices[lg_element] = eye(self.basis.dimension)
      self._characters[lg_element] = self.basis.dimension
    else:
      self._rep_matrices[lg_element] = self.basis.matrix(E).pinv() * self.basis.matrix(lg_element)


  @property
  def irreducible(self):
    inner_product = S.Zero
    for element in self.little_group.elements:
      inner_product += conjugate(self.getCharacter(element)) * self.getCharacter(element)

    inner_product = expand(inner_product)
    if inner_product == self.little_group.order:
      return True

    return False

  @property
  def reducible(self):
    return not self.irreducible


class OperatorBasis:

  def __init__(self, *in_operators):

    if not in_operators:
      raise ValueError("Must provide at least one basis operator")

    self._operators = in_operators

    self._bosonic = bosonic(in_operators[0])
    same_type = all(bosonic(in_op) == self._bosonic for in_op in in_operators)
    if not same_type:
      raise ValueError("Operators must all be bosonic or all be fermionic")

    self._momentum = momentum(in_operators[0])
    same_mom = all(momentum(in_op) == self._momentum for in_op in in_operators)

    '''
    Can't do momentum equivalent tests with this check...
    if not same_mom:
      raise ValueError("Operators must all have the same total momentum")
    '''

    self._create_grassmann_basis(in_operators)

    self._matrix = dict()

  @property
  def bosonic(self):
    return self._bosonic

  @property
  def fermionic(self):
    return self._fermionic

  @property
  def momentum(self):
    return self._momentum

  @property
  def grassmann_basis(self):
    return self._grassmann_basis

  @property
  def dimension(self):
    return len(self.operators)

  @property
  def operators(self):
    return self._operators

  def _create_grassmann_basis(self, in_operators):

    op_terms = list()
    for in_operator in in_operators:
      op_terms.extend(list(terms(in_operator)))

    self._grassmann_basis = set(op_terms)

  def vector(self, operator):
    _vector = list()
    coeffs_dict = coefficients(operator).copy()
    for grassmann_vector in self.grassmann_basis:
      if grassmann_vector in coeffs_dict:
        _vector.append(coeffs_dict[grassmann_vector])
        del coeffs_dict[grassmann_vector]
      else:
        _vector.append(0)

    if coeffs_dict:
      raise OperatorRepresentationError("Basis is not complete")

    return _vector

  def matrix(self, element):
    if element not in self._matrix:
      mat = list()
      self.rotate(element)
      for operator in self.operators:
        mat.append(self.vector(operator))

      self._matrix[element] = Matrix(mat).T

    return self._matrix[element]

  @property
  def operators(self):
    return self._operators

  def rotate(self, element):
    for operator in self.operators:
      rotate_to(operator, element)


