from sortedcontainers import SortedSet

from context import operators
import operators.cubic_rotations as cr

lg = cr.LittleGroup(True, cr.P([0,0,1]))

left_cosets = list()
right_cosets = list()
for element in cr._POINT_GROUP:
  left_coset = SortedSet()
  right_coset = SortedSet()
  for h in lg.elements:
    left_coset.add(element*h)
    right_coset.add(h*element)

  if left_coset not in left_cosets:
    left_cosets.append(left_coset)

  if right_coset not in right_cosets:
    right_cosets.append(right_coset)


print("Left Cosets")
for left_coset in left_cosets:
  print("\t{}\n".format(left_coset))

'''
pz = cr.P([0,0,1])
print(cr.I_C2c.inverse()*cr.C3a)
#print(cr.C3a*pz)
'''

