from sympy import IndexedBase

class BuildingBlock(IndexedBase):

  is_commutative = False

  def __getitem__(self, indices):
    if is_sequence(indices):
      return IndexedBuildingBlock(self, *indices)
    else:
      return IndexedBuildingBlock(self, indices)


class IndexedBuildingBlock(Indexed):

  is_commutative = False

