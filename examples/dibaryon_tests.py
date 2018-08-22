from sympy import Eijk
from sympy import Array
from sympy import S
from sympy import simplify

from context import operators

from operators.operators import DiracIdx, ColorIdx, OperatorRepresentation, Operator
from operators.cubic_rotations import _POINT_GROUP, P, P0
from operators.cubic_rotations import *
from operators.tensors import Gamma

from pprint import pprint

def test_P1_A1(ops, more_ops=[]):

  print("P1 A1 equivalent momentum test...")
  fails = 0

  ops_1 = list()
  ops_1.extend(P1_A1_A1_1_0(ops, P([0,0,1])))
  ops_1.extend(P1_A1_A1_1_0(ops, P([0,1,0])))
  ops_1.extend(P1_A1_A1_1_0(ops, P([1,0,0])))
  ops_1.extend(P1_A1_A1_1_0(ops, P([0,0,-1])))
  ops_1.extend(P1_A1_A1_1_0(ops, P([0,-1,0])))
  ops_1.extend(P1_A1_A1_1_0(ops, P([-1,0,0])))

  op1_rep = OperatorRepresentation(*ops_1)

  ops_2 = list()
  ops_2.extend(P1_A1_A1_2_1(ops, P([0,0,1])))
  ops_2.extend(P1_A1_A1_2_1(ops, P([0,1,0])))
  ops_2.extend(P1_A1_A1_2_1(ops, P([1,0,0])))
  ops_2.extend(P1_A1_A1_2_1(ops, P([0,0,-1])))
  ops_2.extend(P1_A1_A1_2_1(ops, P([0,-1,0])))
  ops_2.extend(P1_A1_A1_2_1(ops, P([-1,0,0])))

  op2_rep = OperatorRepresentation(*ops_2)

  ops_3 = list()
  ops_3.extend(P1_A1_E2_2_1(ops, P([0,0,1])))
  ops_3.extend(P1_A1_E2_2_1(ops, P([0,1,0])))
  ops_3.extend(P1_A1_E2_2_1(ops, P([1,0,0])))
  ops_3.extend(P1_A1_E2_2_1(ops, P([0,0,-1])))
  ops_3.extend(P1_A1_E2_2_1(ops, P([0,-1,0])))
  ops_3.extend(P1_A1_E2_2_1(ops, P([-1,0,0])))

  op3_rep = OperatorRepresentation(*ops_3)

  if more_ops:
    ops_4 = list()
    ops_4.extend(P1_A1_E2_2_1(more_ops, P([0,0,1])))
    ops_4.extend(P1_A1_E2_2_1(more_ops, P([0,1,0])))
    ops_4.extend(P1_A1_E2_2_1(more_ops, P([1,0,0])))
    ops_4.extend(P1_A1_E2_2_1(more_ops, P([0,0,-1])))
    ops_4.extend(P1_A1_E2_2_1(more_ops, P([0,-1,0])))
    ops_4.extend(P1_A1_E2_2_1(more_ops, P([-1,0,0])))

    op4_rep = OperatorRepresentation(*ops_4)

  for element in _POINT_GROUP:
    #print(element)
    op1_mat = op1_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op1_mat)
    op2_mat = op2_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op2_mat)
    op3_mat = op3_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op3_mat)
    if more_ops:
      op4_mat = op4_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
      #pprint(op4_mat)

    if op1_mat != op2_mat:
      print("{}: op1 != op2".format(element))
      fails += 1

    if op1_mat != op3_mat:
      print("{}: op1 != op3".format(element))
      fails += 1

    if op2_mat != op3_mat:
      print("{}: op2 != op3".format(element))
      fails += 1

    if more_ops:
      if op1_mat != op4_mat:
        print("{}: op1 != op4".format(element))
        fails += 1

  if fails:
    print("TEST FAILED")
  else:
    print("TEST PASSED")

def test_P1_A2(ops, more_ops=[]):

  print("P1 A2 equivalent momentum test...")
  fails = 0

  ops_1 = list()
  ops_1.extend(P1_A2_A1_1_0(ops, P([0,0,1])))
  ops_1.extend(P1_A2_A1_1_0(ops, P([0,1,0])))
  ops_1.extend(P1_A2_A1_1_0(ops, P([1,0,0])))
  ops_1.extend(P1_A2_A1_1_0(ops, P([0,0,-1])))
  ops_1.extend(P1_A2_A1_1_0(ops, P([0,-1,0])))
  ops_1.extend(P1_A2_A1_1_0(ops, P([-1,0,0])))

  op1_rep = OperatorRepresentation(*ops_1)

  ops_2 = list()
  ops_2.extend(P1_A2_A1_2_1(ops, P([0,0,1])))
  ops_2.extend(P1_A2_A1_2_1(ops, P([0,1,0])))
  ops_2.extend(P1_A2_A1_2_1(ops, P([1,0,0])))
  ops_2.extend(P1_A2_A1_2_1(ops, P([0,0,-1])))
  ops_2.extend(P1_A2_A1_2_1(ops, P([0,-1,0])))
  ops_2.extend(P1_A2_A1_2_1(ops, P([-1,0,0])))

  op2_rep = OperatorRepresentation(*ops_2)

  ops_3 = list()
  ops_3.extend(P1_A2_E2_2_1(ops, P([0,0,1])))
  ops_3.extend(P1_A2_E2_2_1(ops, P([0,1,0])))
  ops_3.extend(P1_A2_E2_2_1(ops, P([1,0,0])))
  ops_3.extend(P1_A2_E2_2_1(ops, P([0,0,-1])))
  ops_3.extend(P1_A2_E2_2_1(ops, P([0,-1,0])))
  ops_3.extend(P1_A2_E2_2_1(ops, P([-1,0,0])))

  op3_rep = OperatorRepresentation(*ops_3)

  if more_ops:
    ops_4 = list()
    ops_4.extend(P1_A2_E2_2_1(more_ops, P([0,0,1])))
    ops_4.extend(P1_A2_E2_2_1(more_ops, P([0,1,0])))
    ops_4.extend(P1_A2_E2_2_1(more_ops, P([1,0,0])))
    ops_4.extend(P1_A2_E2_2_1(more_ops, P([0,0,-1])))
    ops_4.extend(P1_A2_E2_2_1(more_ops, P([0,-1,0])))
    ops_4.extend(P1_A2_E2_2_1(more_ops, P([-1,0,0])))

    op4_rep = OperatorRepresentation(*ops_4)

  for element in _POINT_GROUP:
    #print(element)
    op1_mat = op1_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op1_mat)
    op2_mat = op2_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op2_mat)
    op3_mat = op3_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op3_mat)
    if more_ops:
      op4_mat = op4_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
      #pprint(op4_mat)


    if op1_mat != op2_mat:
      print("{}: op1 != op2".format(element))
      pprint(op1_mat)
      pprint(op2_mat)
      fails += 1

    if op1_mat != op3_mat:
      print("{}: op1 != op3".format(element))
      pprint(op1_mat)
      pprint(op3_mat)
      fails += 1

    if op2_mat != op3_mat:
      print("{}: op2 != op3".format(element))
      pprint(op2_mat)
      pprint(op3_mat)
      fails += 1

    if more_ops:
      if op1_mat != op4_mat:
        print("{}: op1 != op4".format(element))
        pprint(op1_mat)
        pprint(op4_mat)
        fails += 1

  if fails:
    print("TEST FAILED")
  else:
    print("TEST PASSED")

def test_P1_B1(ops, more_ops=[]):

  print("P1 B1 equivalent momentum test...")
  fails = 0

  ops_1 = list()
  ops_1.extend(P1_B1_B1_2_1(ops, P([0,0,1])))
  ops_1.extend(P1_B1_B1_2_1(ops, P([0,1,0])))
  ops_1.extend(P1_B1_B1_2_1(ops, P([1,0,0])))
  ops_1.extend(P1_B1_B1_2_1(ops, P([0,0,-1])))
  ops_1.extend(P1_B1_B1_2_1(ops, P([0,-1,0])))
  ops_1.extend(P1_B1_B1_2_1(ops, P([-1,0,0])))

  op1_rep = OperatorRepresentation(*ops_1)

  ops_2 = list()
  ops_2.extend(P1_B1_E2_2_1(ops, P([0,0,1])))
  ops_2.extend(P1_B1_E2_2_1(ops, P([0,1,0])))
  ops_2.extend(P1_B1_E2_2_1(ops, P([1,0,0])))
  ops_2.extend(P1_B1_E2_2_1(ops, P([0,0,-1])))
  ops_2.extend(P1_B1_E2_2_1(ops, P([0,-1,0])))
  ops_2.extend(P1_B1_E2_2_1(ops, P([-1,0,0])))

  op2_rep = OperatorRepresentation(*ops_2)

  if more_ops:
    ops_3 = list()
    ops_3.extend(P1_B1_E2_2_1(more_ops, P([0,0,1])))
    ops_3.extend(P1_B1_E2_2_1(more_ops, P([0,1,0])))
    ops_3.extend(P1_B1_E2_2_1(more_ops, P([1,0,0])))
    ops_3.extend(P1_B1_E2_2_1(more_ops, P([0,0,-1])))
    ops_3.extend(P1_B1_E2_2_1(more_ops, P([0,-1,0])))
    ops_3.extend(P1_B1_E2_2_1(more_ops, P([-1,0,0])))

    op3_rep = OperatorRepresentation(*ops_3)

  for element in _POINT_GROUP:
    #print(element)
    op1_mat = op1_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op1_mat)
    op2_mat = op2_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op2_mat)
    if more_ops:
      op3_mat = op3_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
      #pprint(op3_mat)

    if op1_mat != op2_mat:
      print("{}: op1 != op2".format(element))
      fails += 1

    if more_ops:
      if op1_mat != op3_mat:
        print("{}: op1 != op3".format(element))
        fails += 1

  if fails:
    print("TEST FAILED")
  else:
    print("TEST PASSED")

def test_P1_B2(ops, more_ops=[]):

  print("P1 B2 equivalent momentum test...")
  fails = 0

  ops_1 = list()
  ops_1.extend(P1_B2_B1_2_1(ops, P([0,0,1])))
  ops_1.extend(P1_B2_B1_2_1(ops, P([0,1,0])))
  ops_1.extend(P1_B2_B1_2_1(ops, P([1,0,0])))
  ops_1.extend(P1_B2_B1_2_1(ops, P([0,0,-1])))
  ops_1.extend(P1_B2_B1_2_1(ops, P([0,-1,0])))
  ops_1.extend(P1_B2_B1_2_1(ops, P([-1,0,0])))

  op1_rep = OperatorRepresentation(*ops_1)

  ops_2 = list()
  ops_2.extend(P1_B2_E2_2_1(ops, P([0,0,1])))
  ops_2.extend(P1_B2_E2_2_1(ops, P([0,1,0])))
  ops_2.extend(P1_B2_E2_2_1(ops, P([1,0,0])))
  ops_2.extend(P1_B2_E2_2_1(ops, P([0,0,-1])))
  ops_2.extend(P1_B2_E2_2_1(ops, P([0,-1,0])))
  ops_2.extend(P1_B2_E2_2_1(ops, P([-1,0,0])))

  op2_rep = OperatorRepresentation(*ops_2)

  if more_ops:
    ops_3 = list()
    ops_3.extend(P1_B2_E2_2_1(more_ops, P([0,0,1])))
    ops_3.extend(P1_B2_E2_2_1(more_ops, P([0,1,0])))
    ops_3.extend(P1_B2_E2_2_1(more_ops, P([1,0,0])))
    ops_3.extend(P1_B2_E2_2_1(more_ops, P([0,0,-1])))
    ops_3.extend(P1_B2_E2_2_1(more_ops, P([0,-1,0])))
    ops_3.extend(P1_B2_E2_2_1(more_ops, P([-1,0,0])))

    op3_rep = OperatorRepresentation(*ops_3)

  for element in _POINT_GROUP:
    #print(element)
    op1_mat = op1_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op1_mat)
    op2_mat = op2_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op2_mat)
    if more_ops:
      op3_mat = op3_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
      #pprint(op3_mat)

    if op1_mat != op2_mat:
      print("{}: op1 != op2".format(element))
      fails += 1

    if more_ops:
      if op1_mat != op3_mat:
        print("{}: op1 != op3".format(element))
        fails += 1

  if fails:
    print("TEST FAILED")
  else:
    print("TEST PASSED")

def test_P1_E2(ops, more_ops=[]):

  print("P1 E2 equivalent momentum test...")
  fails = 0

  ops_1 = list()
  ops_1.extend(P1_E2_A1_1_0(ops, P([0,0,1])))
  ops_1.extend(P1_E2_A1_1_0(ops, P([0,1,0])))
  ops_1.extend(P1_E2_A1_1_0(ops, P([1,0,0])))
  ops_1.extend(P1_E2_A1_1_0(ops, P([0,0,-1])))
  ops_1.extend(P1_E2_A1_1_0(ops, P([0,-1,0])))
  ops_1.extend(P1_E2_A1_1_0(ops, P([-1,0,0])))

  op1_rep = OperatorRepresentation(*ops_1)

  ops_2 = list()
  ops_2.extend(P1_E2_E2_2_1_S0(ops, P([0,0,1])))
  ops_2.extend(P1_E2_E2_2_1_S0(ops, P([0,1,0])))
  ops_2.extend(P1_E2_E2_2_1_S0(ops, P([1,0,0])))
  ops_2.extend(P1_E2_E2_2_1_S0(ops, P([0,0,-1])))
  ops_2.extend(P1_E2_E2_2_1_S0(ops, P([0,-1,0])))
  ops_2.extend(P1_E2_E2_2_1_S0(ops, P([-1,0,0])))

  op2_rep = OperatorRepresentation(*ops_2)

  ops_3 = list()
  ops_3.extend(P1_E2_A1_2_1(ops, P([0,0,1])))
  ops_3.extend(P1_E2_A1_2_1(ops, P([0,1,0])))
  ops_3.extend(P1_E2_A1_2_1(ops, P([1,0,0])))
  ops_3.extend(P1_E2_A1_2_1(ops, P([0,0,-1])))
  ops_3.extend(P1_E2_A1_2_1(ops, P([0,-1,0])))
  ops_3.extend(P1_E2_A1_2_1(ops, P([-1,0,0])))

  op3_rep = OperatorRepresentation(*ops_3)

  ops_4 = list()
  ops_4.extend(P1_E2_B1_2_1(ops, P([0,0,1])))
  ops_4.extend(P1_E2_B1_2_1(ops, P([0,1,0])))
  ops_4.extend(P1_E2_B1_2_1(ops, P([1,0,0])))
  ops_4.extend(P1_E2_B1_2_1(ops, P([0,0,-1])))
  ops_4.extend(P1_E2_B1_2_1(ops, P([0,-1,0])))
  ops_4.extend(P1_E2_B1_2_1(ops, P([-1,0,0])))

  op4_rep = OperatorRepresentation(*ops_4)

  ops_5 = list()
  ops_5.extend(P1_E2_E2_2_1_S1(ops, P([0,0,1])))
  ops_5.extend(P1_E2_E2_2_1_S1(ops, P([0,1,0])))
  ops_5.extend(P1_E2_E2_2_1_S1(ops, P([1,0,0])))
  ops_5.extend(P1_E2_E2_2_1_S1(ops, P([0,0,-1])))
  ops_5.extend(P1_E2_E2_2_1_S1(ops, P([0,-1,0])))
  ops_5.extend(P1_E2_E2_2_1_S1(ops, P([-1,0,0])))

  op5_rep = OperatorRepresentation(*ops_5)

  if more_ops:
    ops_6 = list()
    ops_6.extend(P1_E2_A1_1_0(more_ops, P([0,0,1])))
    ops_6.extend(P1_E2_A1_1_0(more_ops, P([0,1,0])))
    ops_6.extend(P1_E2_A1_1_0(more_ops, P([1,0,0])))
    ops_6.extend(P1_E2_A1_1_0(more_ops, P([0,0,-1])))
    ops_6.extend(P1_E2_A1_1_0(more_ops, P([0,-1,0])))
    ops_6.extend(P1_E2_A1_1_0(more_ops, P([-1,0,0])))

    op6_rep = OperatorRepresentation(*ops_6)

  for element in _POINT_GROUP:
    #print(element)
    op1_mat = op1_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op1_mat)
    op2_mat = op2_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op2_mat)
    op3_mat = op3_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op3_mat)
    op4_mat = op4_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op4_mat)
    op5_mat = op5_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op5_mat)
    if more_ops:
      op6_mat = op6_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
      #pprint(op6_mat)

    if op1_mat != op2_mat:
      print("{}: op1 != op2".format(element))
      pprint(op1_mat)
      pprint(op2_mat)
      fails += 1

    if op1_mat != op3_mat:
      print("{}: op1 != op3".format(element))
      fails += 1

    if op1_mat != op4_mat:
      print("{}: op1 != op4".format(element))
      fails += 1

    if op1_mat != op5_mat:
      print("{}: op1 != op5".format(element))
      fails += 1
      pprint(op1_mat)
      pprint(op5_mat)

    if op2_mat != op3_mat:
      print("{}: op2 != op3".format(element))
      fails += 1

    if op2_mat != op4_mat:
      print("{}: op2 != op4".format(element))
      fails += 1

    if op2_mat != op5_mat:
      print("{}: op2 != op5".format(element))
      fails += 1

    if op3_mat != op4_mat:
      print("{}: op3 != op4".format(element))
      fails += 1

    if op3_mat != op5_mat:
      print("{}: op3 != op5".format(element))
      fails += 1

    if op4_mat != op5_mat:
      print("{}: op4 != op5".format(element))
      fails += 1

    if more_ops:
      if op1_mat != op6_mat:
        print("{}: op1 != op6".format(element))
        fails += 1

  if fails:
    print("TEST FAILED")
  else:
    print("TEST PASSED")


def test_P2_A1(ops_sym, ops_asym):
  print("P2 A1 equivalent momentum test...")
  fails = 0

  ops_1 = list()
  ops_1.extend(P2_A1_A1_2_0(ops_sym, P([0,1,1])))
  ops_1.extend(P2_A1_A1_2_0(ops_sym, P([0,-1,-1])))
  ops_1.extend(P2_A1_A1_2_0(ops_sym, P([0,1,-1])))
  ops_1.extend(P2_A1_A1_2_0(ops_sym, P([0,-1,1])))
  ops_1.extend(P2_A1_A1_2_0(ops_sym, P([1,0,1])))
  ops_1.extend(P2_A1_A1_2_0(ops_sym, P([-1,0,-1])))
  ops_1.extend(P2_A1_A1_2_0(ops_sym, P([1,0,-1])))
  ops_1.extend(P2_A1_A1_2_0(ops_sym, P([-1,0,1])))
  ops_1.extend(P2_A1_A1_2_0(ops_sym, P([1,1,0])))
  ops_1.extend(P2_A1_A1_2_0(ops_sym, P([-1,-1,0])))
  ops_1.extend(P2_A1_A1_2_0(ops_sym, P([1,-1,0])))
  ops_1.extend(P2_A1_A1_2_0(ops_sym, P([-1,1,0])))

  op1_rep = OperatorRepresentation(*ops_1)

  ops_2 = list()
  ops_2.extend(P2_A1_A1_2_0(ops_asym, P([0,1,1])))
  ops_2.extend(P2_A1_A1_2_0(ops_asym, P([0,-1,-1])))
  ops_2.extend(P2_A1_A1_2_0(ops_asym, P([0,1,-1])))
  ops_2.extend(P2_A1_A1_2_0(ops_asym, P([0,-1,1])))
  ops_2.extend(P2_A1_A1_2_0(ops_asym, P([1,0,1])))
  ops_2.extend(P2_A1_A1_2_0(ops_asym, P([-1,0,-1])))
  ops_2.extend(P2_A1_A1_2_0(ops_asym, P([1,0,-1])))
  ops_2.extend(P2_A1_A1_2_0(ops_asym, P([-1,0,1])))
  ops_2.extend(P2_A1_A1_2_0(ops_asym, P([1,1,0])))
  ops_2.extend(P2_A1_A1_2_0(ops_asym, P([-1,-1,0])))
  ops_2.extend(P2_A1_A1_2_0(ops_asym, P([1,-1,0])))
  ops_2.extend(P2_A1_A1_2_0(ops_asym, P([-1,1,0])))

  op2_rep = OperatorRepresentation(*ops_2)

  ops_3 = list()
  ops_3.extend(P2_A1_A1_1_1(ops_sym, P([0,1,1])))
  ops_3.extend(P2_A1_A1_1_1(ops_sym, P([0,-1,-1])))
  ops_3.extend(P2_A1_A1_1_1(ops_sym, P([0,1,-1])))
  ops_3.extend(P2_A1_A1_1_1(ops_sym, P([0,-1,1])))
  ops_3.extend(P2_A1_A1_1_1(ops_sym, P([1,0,1])))
  ops_3.extend(P2_A1_A1_1_1(ops_sym, P([-1,0,-1])))
  ops_3.extend(P2_A1_A1_1_1(ops_sym, P([1,0,-1])))
  ops_3.extend(P2_A1_A1_1_1(ops_sym, P([-1,0,1])))
  ops_3.extend(P2_A1_A1_1_1(ops_sym, P([1,1,0])))
  ops_3.extend(P2_A1_A1_1_1(ops_sym, P([-1,-1,0])))
  ops_3.extend(P2_A1_A1_1_1(ops_sym, P([1,-1,0])))
  ops_3.extend(P2_A1_A1_1_1(ops_sym, P([-1,1,0])))

  op3_rep = OperatorRepresentation(*ops_3)

  ops_4 = list()
  ops_4.extend(P2_A1_B1_1_1(ops_sym, P([0,1,1])))
  ops_4.extend(P2_A1_B1_1_1(ops_sym, P([0,-1,-1])))
  ops_4.extend(P2_A1_B1_1_1(ops_sym, P([0,1,-1])))
  ops_4.extend(P2_A1_B1_1_1(ops_sym, P([0,-1,1])))
  ops_4.extend(P2_A1_B1_1_1(ops_sym, P([1,0,1])))
  ops_4.extend(P2_A1_B1_1_1(ops_sym, P([-1,0,-1])))
  ops_4.extend(P2_A1_B1_1_1(ops_sym, P([1,0,-1])))
  ops_4.extend(P2_A1_B1_1_1(ops_sym, P([-1,0,1])))
  ops_4.extend(P2_A1_B1_1_1(ops_sym, P([1,1,0])))
  ops_4.extend(P2_A1_B1_1_1(ops_sym, P([-1,-1,0])))
  ops_4.extend(P2_A1_B1_1_1(ops_sym, P([1,-1,0])))
  ops_4.extend(P2_A1_B1_1_1(ops_sym, P([-1,1,0])))

  op4_rep = OperatorRepresentation(*ops_4)

  for element in _POINT_GROUP:
    #print(element)
    op1_mat = op1_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op1_mat)
    op2_mat = op2_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op2_mat)
    op3_mat = op3_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op3_mat)
    op4_mat = op4_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op4_mat)

    if op1_mat != op2_mat:
      print("{}: op1 != op2".format(element))
      fails += 1

    if op1_mat != op3_mat:
      print("{}: op1 != op3".format(element))
      fails += 1

    if op1_mat != op4_mat:
      print("{}: op1 != op4".format(element))
      fails += 1

    if op2_mat != op3_mat:
      print("{}: op2 != op3".format(element))
      fails += 1

    if op2_mat != op4_mat:
      print("{}: op2 != op4".format(element))
      fails += 1

    if op3_mat != op4_mat:
      print("{}: op3 != op4".format(element))
      fails += 1

  if fails:
    print("TEST FAILED")
  else:
    print("TEST PASSED")


def test_P2_A2(ops_sym, ops_asym):
  print("P2 A2 equivalent momentum test...")
  fails = 0

  ops_1 = list()
  ops_1.extend(P2_A2_A1_2_0(ops_sym, P([0,1,1])))
  ops_1.extend(P2_A2_A1_2_0(ops_sym, P([0,-1,-1])))
  ops_1.extend(P2_A2_A1_2_0(ops_sym, P([0,1,-1])))
  ops_1.extend(P2_A2_A1_2_0(ops_sym, P([0,-1,1])))
  ops_1.extend(P2_A2_A1_2_0(ops_sym, P([1,0,1])))
  ops_1.extend(P2_A2_A1_2_0(ops_sym, P([-1,0,-1])))
  ops_1.extend(P2_A2_A1_2_0(ops_sym, P([1,0,-1])))
  ops_1.extend(P2_A2_A1_2_0(ops_sym, P([-1,0,1])))
  ops_1.extend(P2_A2_A1_2_0(ops_sym, P([1,1,0])))
  ops_1.extend(P2_A2_A1_2_0(ops_sym, P([-1,-1,0])))
  ops_1.extend(P2_A2_A1_2_0(ops_sym, P([1,-1,0])))
  ops_1.extend(P2_A2_A1_2_0(ops_sym, P([-1,1,0])))

  op1_rep = OperatorRepresentation(*ops_1)

  ops_2 = list()
  ops_2.extend(P2_A2_A1_2_0(ops_asym, P([0,1,1])))
  ops_2.extend(P2_A2_A1_2_0(ops_asym, P([0,-1,-1])))
  ops_2.extend(P2_A2_A1_2_0(ops_asym, P([0,1,-1])))
  ops_2.extend(P2_A2_A1_2_0(ops_asym, P([0,-1,1])))
  ops_2.extend(P2_A2_A1_2_0(ops_asym, P([1,0,1])))
  ops_2.extend(P2_A2_A1_2_0(ops_asym, P([-1,0,-1])))
  ops_2.extend(P2_A2_A1_2_0(ops_asym, P([1,0,-1])))
  ops_2.extend(P2_A2_A1_2_0(ops_asym, P([-1,0,1])))
  ops_2.extend(P2_A2_A1_2_0(ops_asym, P([1,1,0])))
  ops_2.extend(P2_A2_A1_2_0(ops_asym, P([-1,-1,0])))
  ops_2.extend(P2_A2_A1_2_0(ops_asym, P([1,-1,0])))
  ops_2.extend(P2_A2_A1_2_0(ops_asym, P([-1,1,0])))

  op2_rep = OperatorRepresentation(*ops_2)

  ops_3 = list()
  ops_3.extend(P2_A2_B1_1_1_Fs(ops_sym, P([0,1,1])))
  ops_3.extend(P2_A2_B1_1_1_Fs(ops_sym, P([0,-1,-1])))
  ops_3.extend(P2_A2_B1_1_1_Fs(ops_sym, P([0,1,-1])))
  ops_3.extend(P2_A2_B1_1_1_Fs(ops_sym, P([0,-1,1])))
  ops_3.extend(P2_A2_B1_1_1_Fs(ops_sym, P([1,0,1])))
  ops_3.extend(P2_A2_B1_1_1_Fs(ops_sym, P([-1,0,-1])))
  ops_3.extend(P2_A2_B1_1_1_Fs(ops_sym, P([1,0,-1])))
  ops_3.extend(P2_A2_B1_1_1_Fs(ops_sym, P([-1,0,1])))
  ops_3.extend(P2_A2_B1_1_1_Fs(ops_sym, P([1,1,0])))
  ops_3.extend(P2_A2_B1_1_1_Fs(ops_sym, P([-1,-1,0])))
  ops_3.extend(P2_A2_B1_1_1_Fs(ops_sym, P([1,-1,0])))
  ops_3.extend(P2_A2_B1_1_1_Fs(ops_sym, P([-1,1,0])))

  op3_rep = OperatorRepresentation(*ops_3)

  ops_4 = list()
  ops_4.extend(P2_A2_B1_1_1_Fa(ops_asym, P([0,1,1])))
  ops_4.extend(P2_A2_B1_1_1_Fa(ops_asym, P([0,-1,-1])))
  ops_4.extend(P2_A2_B1_1_1_Fa(ops_asym, P([0,1,-1])))
  ops_4.extend(P2_A2_B1_1_1_Fa(ops_asym, P([0,-1,1])))
  ops_4.extend(P2_A2_B1_1_1_Fa(ops_asym, P([1,0,1])))
  ops_4.extend(P2_A2_B1_1_1_Fa(ops_asym, P([-1,0,-1])))
  ops_4.extend(P2_A2_B1_1_1_Fa(ops_asym, P([1,0,-1])))
  ops_4.extend(P2_A2_B1_1_1_Fa(ops_asym, P([-1,0,1])))
  ops_4.extend(P2_A2_B1_1_1_Fa(ops_asym, P([1,1,0])))
  ops_4.extend(P2_A2_B1_1_1_Fa(ops_asym, P([-1,-1,0])))
  ops_4.extend(P2_A2_B1_1_1_Fa(ops_asym, P([1,-1,0])))
  ops_4.extend(P2_A2_B1_1_1_Fa(ops_asym, P([-1,1,0])))

  op4_rep = OperatorRepresentation(*ops_4)

  for element in _POINT_GROUP:
    #print(element)
    op1_mat = op1_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op1_mat)
    op2_mat = op2_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op2_mat)
    op3_mat = op3_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op3_mat)
    op4_mat = op4_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op4_mat)

    if op1_mat != op2_mat:
      print("{}: op1 != op2".format(element))
      fails += 1

    if op1_mat != op3_mat:
      print("{}: op1 != op3".format(element))
      fails += 1

    if op1_mat != op4_mat:
      print("{}: op1 != op4".format(element))
      fails += 1

    if op2_mat != op3_mat:
      print("{}: op2 != op3".format(element))
      fails += 1

    if op2_mat != op4_mat:
      print("{}: op2 != op4".format(element))
      fails += 1

    if op3_mat != op4_mat:
      print("{}: op3 != op4".format(element))
      fails += 1

  if fails:
    print("TEST FAILED")
  else:
    print("TEST PASSED")



def test_P2_B1(ops_sym, ops_asym):
  print("P2 B1 equivalent momentum test...")
  fails = 0

  ops_1 = list()
  ops_1.extend(P2_B1_A1_2_0(ops_sym, P([0,1,1])))
  ops_1.extend(P2_B1_A1_2_0(ops_sym, P([0,-1,-1])))
  ops_1.extend(P2_B1_A1_2_0(ops_sym, P([0,1,-1])))
  ops_1.extend(P2_B1_A1_2_0(ops_sym, P([0,-1,1])))
  ops_1.extend(P2_B1_A1_2_0(ops_sym, P([1,0,1])))
  ops_1.extend(P2_B1_A1_2_0(ops_sym, P([-1,0,-1])))
  ops_1.extend(P2_B1_A1_2_0(ops_sym, P([1,0,-1])))
  ops_1.extend(P2_B1_A1_2_0(ops_sym, P([-1,0,1])))
  ops_1.extend(P2_B1_A1_2_0(ops_sym, P([1,1,0])))
  ops_1.extend(P2_B1_A1_2_0(ops_sym, P([-1,-1,0])))
  ops_1.extend(P2_B1_A1_2_0(ops_sym, P([1,-1,0])))
  ops_1.extend(P2_B1_A1_2_0(ops_sym, P([-1,1,0])))

  op1_rep = OperatorRepresentation(*ops_1)

  ops_2 = list()
  ops_2.extend(P2_B1_A1_2_0(ops_asym, P([0,1,1])))
  ops_2.extend(P2_B1_A1_2_0(ops_asym, P([0,-1,-1])))
  ops_2.extend(P2_B1_A1_2_0(ops_asym, P([0,1,-1])))
  ops_2.extend(P2_B1_A1_2_0(ops_asym, P([0,-1,1])))
  ops_2.extend(P2_B1_A1_2_0(ops_asym, P([1,0,1])))
  ops_2.extend(P2_B1_A1_2_0(ops_asym, P([-1,0,-1])))
  ops_2.extend(P2_B1_A1_2_0(ops_asym, P([1,0,-1])))
  ops_2.extend(P2_B1_A1_2_0(ops_asym, P([-1,0,1])))
  ops_2.extend(P2_B1_A1_2_0(ops_asym, P([1,1,0])))
  ops_2.extend(P2_B1_A1_2_0(ops_asym, P([-1,-1,0])))
  ops_2.extend(P2_B1_A1_2_0(ops_asym, P([1,-1,0])))
  ops_2.extend(P2_B1_A1_2_0(ops_asym, P([-1,1,0])))

  op2_rep = OperatorRepresentation(*ops_2)

  ops_3 = list()
  ops_3.extend(P2_B1_A1_1_1(ops_asym, P([0,1,1])))
  ops_3.extend(P2_B1_A1_1_1(ops_asym, P([0,-1,-1])))
  ops_3.extend(P2_B1_A1_1_1(ops_asym, P([0,1,-1])))
  ops_3.extend(P2_B1_A1_1_1(ops_asym, P([0,-1,1])))
  ops_3.extend(P2_B1_A1_1_1(ops_asym, P([1,0,1])))
  ops_3.extend(P2_B1_A1_1_1(ops_asym, P([-1,0,-1])))
  ops_3.extend(P2_B1_A1_1_1(ops_asym, P([1,0,-1])))
  ops_3.extend(P2_B1_A1_1_1(ops_asym, P([-1,0,1])))
  ops_3.extend(P2_B1_A1_1_1(ops_asym, P([1,1,0])))
  ops_3.extend(P2_B1_A1_1_1(ops_asym, P([-1,-1,0])))
  ops_3.extend(P2_B1_A1_1_1(ops_asym, P([1,-1,0])))
  ops_3.extend(P2_B1_A1_1_1(ops_asym, P([-1,1,0])))

  op3_rep = OperatorRepresentation(*ops_3)

  ops_4 = list()
  ops_4.extend(P2_B1_B1_1_1(ops_asym, P([0,1,1])))
  ops_4.extend(P2_B1_B1_1_1(ops_asym, P([0,-1,-1])))
  ops_4.extend(P2_B1_B1_1_1(ops_asym, P([0,1,-1])))
  ops_4.extend(P2_B1_B1_1_1(ops_asym, P([0,-1,1])))
  ops_4.extend(P2_B1_B1_1_1(ops_asym, P([1,0,1])))
  ops_4.extend(P2_B1_B1_1_1(ops_asym, P([-1,0,-1])))
  ops_4.extend(P2_B1_B1_1_1(ops_asym, P([1,0,-1])))
  ops_4.extend(P2_B1_B1_1_1(ops_asym, P([-1,0,1])))
  ops_4.extend(P2_B1_B1_1_1(ops_asym, P([1,1,0])))
  ops_4.extend(P2_B1_B1_1_1(ops_asym, P([-1,-1,0])))
  ops_4.extend(P2_B1_B1_1_1(ops_asym, P([1,-1,0])))
  ops_4.extend(P2_B1_B1_1_1(ops_asym, P([-1,1,0])))

  op4_rep = OperatorRepresentation(*ops_4)

  for element in _POINT_GROUP:
    #print(element)
    op1_mat = op1_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op1_mat)
    op2_mat = op2_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op2_mat)
    op3_mat = op3_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op3_mat)
    op4_mat = op4_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op4_mat)

    if op1_mat != op2_mat:
      print("{}: op1 != op2".format(element))
      fails += 1

    if op1_mat != op3_mat:
      print("{}: op1 != op3".format(element))
      fails += 1

    if op1_mat != op4_mat:
      print("{}: op1 != op4".format(element))
      fails += 1

    if op2_mat != op3_mat:
      print("{}: op2 != op3".format(element))
      fails += 1

    if op2_mat != op4_mat:
      print("{}: op2 != op4".format(element))
      fails += 1

    if op3_mat != op4_mat:
      print("{}: op3 != op4".format(element))
      fails += 1

  if fails:
    print("TEST FAILED")
  else:
    print("TEST PASSED")



def test_P2_B2(ops_sym, ops_asym):

  print("P2 B2 equivalent momentum test...")
  fails = 0

  ops_1 = list()
  ops_1.extend(P2_B2_A1_2_0(ops_sym, P([0,1,1])))
  ops_1.extend(P2_B2_A1_2_0(ops_sym, P([0,-1,-1])))
  ops_1.extend(P2_B2_A1_2_0(ops_sym, P([0,1,-1])))
  ops_1.extend(P2_B2_A1_2_0(ops_sym, P([0,-1,1])))
  ops_1.extend(P2_B2_A1_2_0(ops_sym, P([1,0,1])))
  ops_1.extend(P2_B2_A1_2_0(ops_sym, P([-1,0,-1])))
  ops_1.extend(P2_B2_A1_2_0(ops_sym, P([1,0,-1])))
  ops_1.extend(P2_B2_A1_2_0(ops_sym, P([-1,0,1])))
  ops_1.extend(P2_B2_A1_2_0(ops_sym, P([1,1,0])))
  ops_1.extend(P2_B2_A1_2_0(ops_sym, P([-1,-1,0])))
  ops_1.extend(P2_B2_A1_2_0(ops_sym, P([1,-1,0])))
  ops_1.extend(P2_B2_A1_2_0(ops_sym, P([-1,1,0])))

  op1_rep = OperatorRepresentation(*ops_1)

  ops_2 = list()
  ops_2.extend(P2_B2_A1_2_0(ops_asym, P([0,1,1])))
  ops_2.extend(P2_B2_A1_2_0(ops_asym, P([0,-1,-1])))
  ops_2.extend(P2_B2_A1_2_0(ops_asym, P([0,1,-1])))
  ops_2.extend(P2_B2_A1_2_0(ops_asym, P([0,-1,1])))
  ops_2.extend(P2_B2_A1_2_0(ops_asym, P([1,0,1])))
  ops_2.extend(P2_B2_A1_2_0(ops_asym, P([-1,0,-1])))
  ops_2.extend(P2_B2_A1_2_0(ops_asym, P([1,0,-1])))
  ops_2.extend(P2_B2_A1_2_0(ops_asym, P([-1,0,1])))
  ops_2.extend(P2_B2_A1_2_0(ops_asym, P([1,1,0])))
  ops_2.extend(P2_B2_A1_2_0(ops_asym, P([-1,-1,0])))
  ops_2.extend(P2_B2_A1_2_0(ops_asym, P([1,-1,0])))
  ops_2.extend(P2_B2_A1_2_0(ops_asym, P([-1,1,0])))

  op2_rep = OperatorRepresentation(*ops_2)

  ops_3 = list()
  ops_3.extend(P2_B2_B1_1_1(ops_sym, P([0,1,1])))
  ops_3.extend(P2_B2_B1_1_1(ops_sym, P([0,-1,-1])))
  ops_3.extend(P2_B2_B1_1_1(ops_sym, P([0,1,-1])))
  ops_3.extend(P2_B2_B1_1_1(ops_sym, P([0,-1,1])))
  ops_3.extend(P2_B2_B1_1_1(ops_sym, P([1,0,1])))
  ops_3.extend(P2_B2_B1_1_1(ops_sym, P([-1,0,-1])))
  ops_3.extend(P2_B2_B1_1_1(ops_sym, P([1,0,-1])))
  ops_3.extend(P2_B2_B1_1_1(ops_sym, P([-1,0,1])))
  ops_3.extend(P2_B2_B1_1_1(ops_sym, P([1,1,0])))
  ops_3.extend(P2_B2_B1_1_1(ops_sym, P([-1,-1,0])))
  ops_3.extend(P2_B2_B1_1_1(ops_sym, P([1,-1,0])))
  ops_3.extend(P2_B2_B1_1_1(ops_sym, P([-1,1,0])))

  op3_rep = OperatorRepresentation(*ops_3)

  ops_4 = list()
  ops_4.extend(P2_B2_A1_1_1(ops_asym, P([0,1,1])))
  ops_4.extend(P2_B2_A1_1_1(ops_asym, P([0,-1,-1])))
  ops_4.extend(P2_B2_A1_1_1(ops_asym, P([0,1,-1])))
  ops_4.extend(P2_B2_A1_1_1(ops_asym, P([0,-1,1])))
  ops_4.extend(P2_B2_A1_1_1(ops_asym, P([1,0,1])))
  ops_4.extend(P2_B2_A1_1_1(ops_asym, P([-1,0,-1])))
  ops_4.extend(P2_B2_A1_1_1(ops_asym, P([1,0,-1])))
  ops_4.extend(P2_B2_A1_1_1(ops_asym, P([-1,0,1])))
  ops_4.extend(P2_B2_A1_1_1(ops_asym, P([1,1,0])))
  ops_4.extend(P2_B2_A1_1_1(ops_asym, P([-1,-1,0])))
  ops_4.extend(P2_B2_A1_1_1(ops_asym, P([1,-1,0])))
  ops_4.extend(P2_B2_A1_1_1(ops_asym, P([-1,1,0])))

  op4_rep = OperatorRepresentation(*ops_4)

  for element in _POINT_GROUP:
    #print(element)
    op1_mat = op1_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op1_mat)
    op2_mat = op2_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op2_mat)
    op3_mat = op3_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op3_mat)
    op4_mat = op4_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op4_mat)

    if op1_mat != op2_mat:
      print("{}: op1 != op2".format(element))
      fails += 1

    if op1_mat != op3_mat:
      print("{}: op1 != op3".format(element))
      fails += 1

    if op1_mat != op4_mat:
      print("{}: op1 != op4".format(element))
      fails += 1

    if op2_mat != op3_mat:
      print("{}: op2 != op3".format(element))
      fails += 1

    if op2_mat != op4_mat:
      print("{}: op2 != op4".format(element))
      fails += 1

    if op3_mat != op4_mat:
      print("{}: op3 != op4".format(element))
      fails += 1

  if fails:
    print("TEST FAILED")
  else:
    print("TEST PASSED")


def test_P3_A1(ops, more_ops=[]):
  print("P3 A1 equivalent momentum test...")
  fails = 0

  ops_1 = list()
  ops_1.extend(P3_A1_A1_3_0(ops, P([1,1,1])))
  ops_1.extend(P3_A1_A1_3_0(ops, P([1,1,-1])))
  ops_1.extend(P3_A1_A1_3_0(ops, P([1,-1,1])))
  ops_1.extend(P3_A1_A1_3_0(ops, P([-1,1,1])))
  ops_1.extend(P3_A1_A1_3_0(ops, P([1,-1,-1])))
  ops_1.extend(P3_A1_A1_3_0(ops, P([-1,1,-1])))
  ops_1.extend(P3_A1_A1_3_0(ops, P([-1,-1,1])))
  ops_1.extend(P3_A1_A1_3_0(ops, P([-1,-1,-1])))

  op1_rep = OperatorRepresentation(*ops_1)

  ops_2 = list()
  ops_2.extend(P3_A1_A1_2_1(ops, P([1,1,1])))
  ops_2.extend(P3_A1_A1_2_1(ops, P([1,1,-1])))
  ops_2.extend(P3_A1_A1_2_1(ops, P([1,-1,1])))
  ops_2.extend(P3_A1_A1_2_1(ops, P([-1,1,1])))
  ops_2.extend(P3_A1_A1_2_1(ops, P([1,-1,-1])))
  ops_2.extend(P3_A1_A1_2_1(ops, P([-1,1,-1])))
  ops_2.extend(P3_A1_A1_2_1(ops, P([-1,-1,1])))
  ops_2.extend(P3_A1_A1_2_1(ops, P([-1,-1,-1])))

  op2_rep = OperatorRepresentation(*ops_2)

  ops_3 = list()
  ops_3.extend(P3_A1_E2_2_1(ops, P([1,1,1])))
  ops_3.extend(P3_A1_E2_2_1(ops, P([1,1,-1])))
  ops_3.extend(P3_A1_E2_2_1(ops, P([1,-1,1])))
  ops_3.extend(P3_A1_E2_2_1(ops, P([-1,1,1])))
  ops_3.extend(P3_A1_E2_2_1(ops, P([1,-1,-1])))
  ops_3.extend(P3_A1_E2_2_1(ops, P([-1,1,-1])))
  ops_3.extend(P3_A1_E2_2_1(ops, P([-1,-1,1])))
  ops_3.extend(P3_A1_E2_2_1(ops, P([-1,-1,-1])))

  op3_rep = OperatorRepresentation(*ops_3)

  if more_ops:
    ops_4 = list()
    ops_4.extend(P3_A1_A1_3_0(more_ops, P([1,1,1])))
    ops_4.extend(P3_A1_A1_3_0(more_ops, P([1,1,-1])))
    ops_4.extend(P3_A1_A1_3_0(more_ops, P([1,-1,1])))
    ops_4.extend(P3_A1_A1_3_0(more_ops, P([-1,1,1])))
    ops_4.extend(P3_A1_A1_3_0(more_ops, P([1,-1,-1])))
    ops_4.extend(P3_A1_A1_3_0(more_ops, P([-1,1,-1])))
    ops_4.extend(P3_A1_A1_3_0(more_ops, P([-1,-1,1])))
    ops_4.extend(P3_A1_A1_3_0(more_ops, P([-1,-1,-1])))

    op4_rep = OperatorRepresentation(*ops_4)

  for element in _POINT_GROUP:
    #print(element)
    op1_mat = op1_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op1_mat)
    op2_mat = op2_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op2_mat)
    op3_mat = op3_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op3_mat)
    if more_ops:
      op4_mat = op4_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
      #pprint(op4_mat)


    if op1_mat != op2_mat:
      print("{}: op1 != op2".format(element))
      fails += 1

    if op1_mat != op3_mat:
      print("{}: op1 != op3".format(element))
      fails += 1

    if op2_mat != op3_mat:
      print("{}: op2 != op3".format(element))
      fails += 1

    if more_ops:
      if op1_mat != op4_mat:
        print("{}: op1 != op4".format(element))
        fails += 1

      if op2_mat != op4_mat:
        print("{}: op2 != op4".format(element))
        fails += 1

      if op3_mat != op4_mat:
        print("{}: op3 != op4".format(element))
        fails += 1

  if fails:
    print("TEST FAILED")
  else:
    print("TEST PASSED")


def test_P3_A2(ops, more_ops=[]):
  print("P3 A1 equivalent momentum test...")
  fails = 0

  ops_1 = list()
  ops_1.extend(P3_A2_A1_3_0(ops, P([1,1,1])))
  ops_1.extend(P3_A2_A1_3_0(ops, P([1,1,-1])))
  ops_1.extend(P3_A2_A1_3_0(ops, P([1,-1,1])))
  ops_1.extend(P3_A2_A1_3_0(ops, P([-1,1,1])))
  ops_1.extend(P3_A2_A1_3_0(ops, P([1,-1,-1])))
  ops_1.extend(P3_A2_A1_3_0(ops, P([-1,1,-1])))
  ops_1.extend(P3_A2_A1_3_0(ops, P([-1,-1,1])))
  ops_1.extend(P3_A2_A1_3_0(ops, P([-1,-1,-1])))

  op1_rep = OperatorRepresentation(*ops_1)

  ops_2 = list()
  ops_2.extend(P3_A2_A1_2_1(ops, P([1,1,1])))
  ops_2.extend(P3_A2_A1_2_1(ops, P([1,1,-1])))
  ops_2.extend(P3_A2_A1_2_1(ops, P([1,-1,1])))
  ops_2.extend(P3_A2_A1_2_1(ops, P([-1,1,1])))
  ops_2.extend(P3_A2_A1_2_1(ops, P([1,-1,-1])))
  ops_2.extend(P3_A2_A1_2_1(ops, P([-1,1,-1])))
  ops_2.extend(P3_A2_A1_2_1(ops, P([-1,-1,1])))
  ops_2.extend(P3_A2_A1_2_1(ops, P([-1,-1,-1])))

  op2_rep = OperatorRepresentation(*ops_2)

  ops_3 = list()
  ops_3.extend(P3_A2_E2_2_1(ops, P([1,1,1])))
  ops_3.extend(P3_A2_E2_2_1(ops, P([1,1,-1])))
  ops_3.extend(P3_A2_E2_2_1(ops, P([1,-1,1])))
  ops_3.extend(P3_A2_E2_2_1(ops, P([-1,1,1])))
  ops_3.extend(P3_A2_E2_2_1(ops, P([1,-1,-1])))
  ops_3.extend(P3_A2_E2_2_1(ops, P([-1,1,-1])))
  ops_3.extend(P3_A2_E2_2_1(ops, P([-1,-1,1])))
  ops_3.extend(P3_A2_E2_2_1(ops, P([-1,-1,-1])))

  op3_rep = OperatorRepresentation(*ops_3)

  if more_ops:
    ops_4 = list()
    ops_4.extend(P3_A2_A1_3_0(more_ops, P([1,1,1])))
    ops_4.extend(P3_A2_A1_3_0(more_ops, P([1,1,-1])))
    ops_4.extend(P3_A2_A1_3_0(more_ops, P([1,-1,1])))
    ops_4.extend(P3_A2_A1_3_0(more_ops, P([-1,1,1])))
    ops_4.extend(P3_A2_A1_3_0(more_ops, P([1,-1,-1])))
    ops_4.extend(P3_A2_A1_3_0(more_ops, P([-1,1,-1])))
    ops_4.extend(P3_A2_A1_3_0(more_ops, P([-1,-1,1])))
    ops_4.extend(P3_A2_A1_3_0(more_ops, P([-1,-1,-1])))

  op4_rep = OperatorRepresentation(*ops_4)

  for element in _POINT_GROUP:
    #print(element)
    op1_mat = op1_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op1_mat)
    op2_mat = op2_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op2_mat)
    op3_mat = op3_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
    #pprint(op3_mat)
    if more_ops:
      op4_mat = op4_rep.getRepresentationMatrix(element, True).applyfunc(simplify)
      #pprint(op4_mat)

    if op1_mat != op2_mat:
      print("{}: op1 != op2".format(element))
      fails += 1

    if op1_mat != op3_mat:
      print("{}: op1 != op3".format(element))
      fails += 1

    if op2_mat != op3_mat:
      print("{}: op2 != op3".format(element))
      fails += 1

    if more_ops:
      if op1_mat != op4_mat:
        print("{}: op1 != op4".format(element))
        fails += 1

  if fails:
    print("TEST FAILED")
  else:
    print("TEST PASSED")



def P0_A1p(ops, n=0):
  if n == 0:
    A1p_op = ops[0].projectMomentum(P0, P0)

  elif n == 1:
    A1p_op = ops[0].projectMomentum(P([0,0,1]), P([0,0,-1])) \
        + ops[0].projectMomentum(P([0,1,0]), P([0,-1,0])) \
        + ops[0].projectMomentum(P([1,0,0]), P([-1,0,0])) \
        + ops[0].projectMomentum(P([0,0,-1]), P([0,0,1])) \
        + ops[0].projectMomentum(P([0,-1,0]), P([0,1,0])) \
        + ops[0].projectMomentum(P([-1,0,0]), P([1,0,0]))

  elif n == 2:
    A1p_op = ops[0].projectMomentum(P([0,1,1]), P([0,-1,-1])) \
        + ops[0].projectMomentum(P([0,-1,-1]), P([0,1,1])) \
        + ops[0].projectMomentum(P([0,1,-1]), P([0,-1,1])) \
        + ops[0].projectMomentum(P([0,-1,1]), P([0,1,-1])) \
        + ops[0].projectMomentum(P([1,0,1]), P([-1,0,-1])) \
        + ops[0].projectMomentum(P([-1,0,-1]), P([1,0,1])) \
        + ops[0].projectMomentum(P([1,0,-1]), P([-1,0,1])) \
        + ops[0].projectMomentum(P([-1,0,1]), P([1,0,-1])) \
        + ops[0].projectMomentum(P([1,1,0]), P([-1,-1,0])) \
        + ops[0].projectMomentum(P([-1,-1,0]), P([1,1,0])) \
        + ops[0].projectMomentum(P([1,-1,0]), P([-1,1,0])) \
        + ops[0].projectMomentum(P([-1,1,0]), P([1,-1,0]))

  elif n == 3:
    A1p_op = ops[0].projectMomentum(P([1,1,1]), P([-1,-1,-1])) \
        + ops[0].projectMomentum(P([1,1,-1]), P([-1,-1,1])) \
        + ops[0].projectMomentum(P([1,-1,1]), P([-1,1,-1])) \
        + ops[0].projectMomentum(P([1,-1,-1]), P([-1,1,1])) \
        + ops[0].projectMomentum(P([-1,1,1]), P([1,-1,-1])) \
        + ops[0].projectMomentum(P([-1,1,-1]), P([1,-1,1])) \
        + ops[0].projectMomentum(P([-1,-1,1]), P([1,1,-1])) \
        + ops[0].projectMomentum(P([-1,-1,-1]), P([1,1,1]))

  else:
    print("P0_A1p: pi^2 = {} not implemented".format(n))
    return

  return [A1p_op]


def P0_T1p_A1p(ops, n=0):
  if n == 0:
    T1p_1_op = ops[1].projectMomentum(P([0,0,0]), P([0,0,0]))
    T1p_2_op = ops[2].projectMomentum(P([0,0,0]), P([0,0,0]))
    T1p_3_op = ops[3].projectMomentum(P([0,0,0]), P([0,0,0]))

  elif n == 1:
    T1p_1_op = ops[1].projectMomentum(P([0,0,1]), P([0,0,-1])) \
        + ops[1].projectMomentum(P([0,1,0]), P([0,-1,0])) \
        + ops[1].projectMomentum(P([1,0,0]), P([-1,0,0])) \
        + ops[1].projectMomentum(P([0,0,-1]), P([0,0,1])) \
        + ops[1].projectMomentum(P([0,-1,0]), P([0,1,0])) \
        + ops[1].projectMomentum(P([-1,0,0]), P([1,0,0]))

    T1p_2_op = ops[2].projectMomentum(P([0,0,1]), P([0,0,-1])) \
        + ops[2].projectMomentum(P([0,1,0]), P([0,-1,0])) \
        + ops[2].projectMomentum(P([1,0,0]), P([-1,0,0])) \
        + ops[2].projectMomentum(P([0,0,-1]), P([0,0,1])) \
        + ops[2].projectMomentum(P([0,-1,0]), P([0,1,0])) \
        + ops[2].projectMomentum(P([-1,0,0]), P([1,0,0]))

    T1p_3_op = ops[3].projectMomentum(P([0,0,1]), P([0,0,-1])) \
        + ops[3].projectMomentum(P([0,1,0]), P([0,-1,0])) \
        + ops[3].projectMomentum(P([1,0,0]), P([-1,0,0])) \
        + ops[3].projectMomentum(P([0,0,-1]), P([0,0,1])) \
        + ops[3].projectMomentum(P([0,-1,0]), P([0,1,0])) \
        + ops[3].projectMomentum(P([-1,0,0]), P([1,0,0]))

  elif n == 2:
    T1p_1_op = ops[1].projectMomentum(P([0,1,1]), P([0,-1,-1])) \
        + ops[1].projectMomentum(P([0,-1,-1]), P([0,1,1])) \
        + ops[1].projectMomentum(P([0,1,-1]), P([0,-1,1])) \
        + ops[1].projectMomentum(P([0,-1,1]), P([0,1,-1])) \
        + ops[1].projectMomentum(P([1,0,1]), P([-1,0,-1])) \
        + ops[1].projectMomentum(P([-1,0,-1]), P([1,0,1])) \
        + ops[1].projectMomentum(P([1,0,-1]), P([-1,0,1])) \
        + ops[1].projectMomentum(P([-1,0,1]), P([1,0,-1])) \
        + ops[1].projectMomentum(P([1,1,0]), P([-1,-1,0])) \
        + ops[1].projectMomentum(P([-1,-1,0]), P([1,1,0])) \
        + ops[1].projectMomentum(P([1,-1,0]), P([-1,1,0])) \
        + ops[1].projectMomentum(P([-1,1,0]), P([1,-1,0]))

    T1p_2_op = ops[2].projectMomentum(P([0,1,1]), P([0,-1,-1])) \
        + ops[2].projectMomentum(P([0,-1,-1]), P([0,1,1])) \
        + ops[2].projectMomentum(P([0,1,-1]), P([0,-1,1])) \
        + ops[2].projectMomentum(P([0,-1,1]), P([0,1,-1])) \
        + ops[2].projectMomentum(P([1,0,1]), P([-1,0,-1])) \
        + ops[2].projectMomentum(P([-1,0,-1]), P([1,0,1])) \
        + ops[2].projectMomentum(P([1,0,-1]), P([-1,0,1])) \
        + ops[2].projectMomentum(P([-1,0,1]), P([1,0,-1])) \
        + ops[2].projectMomentum(P([1,1,0]), P([-1,-1,0])) \
        + ops[2].projectMomentum(P([-1,-1,0]), P([1,1,0])) \
        + ops[2].projectMomentum(P([1,-1,0]), P([-1,1,0])) \
        + ops[2].projectMomentum(P([-1,1,0]), P([1,-1,0]))

    T1p_3_op = ops[3].projectMomentum(P([0,1,1]), P([0,-1,-1])) \
        + ops[3].projectMomentum(P([0,-1,-1]), P([0,1,1])) \
        + ops[3].projectMomentum(P([0,1,-1]), P([0,-1,1])) \
        + ops[3].projectMomentum(P([0,-1,1]), P([0,1,-1])) \
        + ops[3].projectMomentum(P([1,0,1]), P([-1,0,-1])) \
        + ops[3].projectMomentum(P([-1,0,-1]), P([1,0,1])) \
        + ops[3].projectMomentum(P([1,0,-1]), P([-1,0,1])) \
        + ops[3].projectMomentum(P([-1,0,1]), P([1,0,-1])) \
        + ops[3].projectMomentum(P([1,1,0]), P([-1,-1,0])) \
        + ops[3].projectMomentum(P([-1,-1,0]), P([1,1,0])) \
        + ops[3].projectMomentum(P([1,-1,0]), P([-1,1,0])) \
        + ops[3].projectMomentum(P([-1,1,0]), P([1,-1,0]))

  elif n == 3:
    T1p_1_op = ops[1].projectMomentum(P([1,1,1]), P([-1,-1,-1])) \
        + ops[1].projectMomentum(P([1,1,-1]), P([-1,-1,1])) \
        + ops[1].projectMomentum(P([1,-1,1]), P([-1,1,-1])) \
        + ops[1].projectMomentum(P([1,-1,-1]), P([-1,1,1])) \
        + ops[1].projectMomentum(P([-1,1,1]), P([1,-1,-1])) \
        + ops[1].projectMomentum(P([-1,1,-1]), P([1,-1,1])) \
        + ops[1].projectMomentum(P([-1,-1,1]), P([1,1,-1])) \
        + ops[1].projectMomentum(P([-1,-1,-1]), P([1,1,1]))

    T1p_2_op = ops[2].projectMomentum(P([1,1,1]), P([-1,-1,-1])) \
        + ops[2].projectMomentum(P([1,1,-1]), P([-1,-1,1])) \
        + ops[2].projectMomentum(P([1,-1,1]), P([-1,1,-1])) \
        + ops[2].projectMomentum(P([1,-1,-1]), P([-1,1,1])) \
        + ops[2].projectMomentum(P([-1,1,1]), P([1,-1,-1])) \
        + ops[2].projectMomentum(P([-1,1,-1]), P([1,-1,1])) \
        + ops[2].projectMomentum(P([-1,-1,1]), P([1,1,-1])) \
        + ops[2].projectMomentum(P([-1,-1,-1]), P([1,1,1]))

    T1p_3_op = ops[3].projectMomentum(P([1,1,1]), P([-1,-1,-1])) \
        + ops[3].projectMomentum(P([1,1,-1]), P([-1,-1,1])) \
        + ops[3].projectMomentum(P([1,-1,1]), P([-1,1,-1])) \
        + ops[3].projectMomentum(P([1,-1,-1]), P([-1,1,1])) \
        + ops[3].projectMomentum(P([-1,1,1]), P([1,-1,-1])) \
        + ops[3].projectMomentum(P([-1,1,-1]), P([1,-1,1])) \
        + ops[3].projectMomentum(P([-1,-1,1]), P([1,1,-1])) \
        + ops[3].projectMomentum(P([-1,-1,-1]), P([1,1,1]))

  else:
    print("P0_T1p_A1p: pi^2 = {} not implemented".format(n))
    return

  return [T1p_1_op, T1p_2_op, T1p_3_op]


def P0_T1p_Ep_1(ops):
  
  T1p_1_op = 3*ops[1].projectMomentum(P([1,0,0]), P([-1,0,0])) \
      - ops[1].projectMomentum(P([1,0,0]), P([-1,0,0])) \
      - ops[1].projectMomentum(P([0,1,0]), P([0,-1,0])) \
      - ops[1].projectMomentum(P([0,0,1]), P([0,0,-1]))

  T1p_2_op = 3*ops[2].projectMomentum(P([0,1,0]), P([0,-1,0])) \
      - ops[2].projectMomentum(P([1,0,0]), P([-1,0,0])) \
      - ops[2].projectMomentum(P([0,1,0]), P([0,-1,0])) \
      - ops[2].projectMomentum(P([0,0,1]), P([0,0,-1]))

  T1p_3_op = 3*ops[3].projectMomentum(P([0,0,1]), P([0,0,-1])) \
      - ops[3].projectMomentum(P([1,0,0]), P([-1,0,0])) \
      - ops[3].projectMomentum(P([0,1,0]), P([0,-1,0])) \
      - ops[3].projectMomentum(P([0,0,1]), P([0,0,-1]))

  return [T1p_1_op, T1p_2_op, T1p_3_op]


# @ADH - Consider other k_{i1} and k_{i2}?
def P0_T1p_Ep_2(ops):

  k11 = P([0,1,0]); k12 = P([0,0,1])
  k21 = P([1,0,0]); k22 = P([0,0,1])
  k31 = P([1,0,0]); k32 = P([0,1,0])

  T1p_1_op = 3*(ops[1].projectMomentum(k11 + k12, -k11 - k12) \
      + ops[1].projectMomentum(k11 - k12, -k11 + k12)) \
      - (ops[1].projectMomentum(P([1,1,0]), P([-1,-1,0])) \
         + ops[1].projectMomentum(P([1,-1,0]), P([-1,1,0]))) \
      - (ops[1].projectMomentum(P([1,0,1]), P([-1,0,-1])) \
         + ops[1].projectMomentum(P([1,0,-1]), P([-1,0,1]))) \
      - (ops[1].projectMomentum(P([0,1,1]), P([0,-1,-1])) \
         + ops[1].projectMomentum(P([0,1,-1]), P([0,-1,1])))

  T1p_2_op = 3*(ops[2].projectMomentum(k21 + k22, -k21 - k22) \
      + ops[2].projectMomentum(k21 - k22, -k21 + k22)) \
      - (ops[2].projectMomentum(P([1,1,0]), P([-1,-1,0])) \
         + ops[2].projectMomentum(P([1,-1,0]), P([-1,1,0]))) \
      - (ops[2].projectMomentum(P([1,0,1]), P([-1,0,-1])) \
         + ops[2].projectMomentum(P([1,0,-1]), P([-1,0,1]))) \
      - (ops[2].projectMomentum(P([0,1,1]), P([0,-1,-1])) \
         + ops[2].projectMomentum(P([0,1,-1]), P([0,-1,1])))

  T1p_3_op = 3*(ops[3].projectMomentum(k31 + k32, -k31 - k32) \
      + ops[3].projectMomentum(k31 - k32, -k31 + k32)) \
      - (ops[3].projectMomentum(P([1,1,0]), P([-1,-1,0])) \
         + ops[3].projectMomentum(P([1,-1,0]), P([-1,1,0]))) \
      - (ops[3].projectMomentum(P([1,0,1]), P([-1,0,-1])) \
         + ops[3].projectMomentum(P([1,0,-1]), P([-1,0,1]))) \
      - (ops[3].projectMomentum(P([0,1,1]), P([0,-1,-1])) \
         + ops[3].projectMomentum(P([0,1,-1]), P([0,-1,1])))

  return [T1p_1_op, T1p_2_op, T1p_3_op]

# @ADH - Consider other k_{i1} and k_{i2}?
def P0_T1p_T2p_2(ops):
  px = P([1,0,0]); py = P([0,1,0]); pz = P([0,0,1])

  k11 = P([0,1,0]); k12 = P([0,0,1])
  k21 = P([1,0,0]); k22 = P([0,0,1])
  k31 = P([1,0,0]); k32 = P([0,1,0])

  T1p_1_op = ops[2].projectMomentum(px + k11, -px - k11) \
      - ops[2].projectMomentum(px - k11, -px + k11) \
      + ops[3].projectMomentum(px + k12, -px - k12) \
      - ops[3].projectMomentum(px - k12, -px + k12)

  T1p_2_op = ops[1].projectMomentum(py + k21, -py - k21) \
      - ops[1].projectMomentum(py - k21, -py + k21) \
      + ops[3].projectMomentum(py + k22, -py - k22) \
      - ops[3].projectMomentum(py - k22, -py + k22)

  T1p_3_op = ops[1].projectMomentum(pz + k31, -pz - k31) \
      - ops[1].projectMomentum(pz - k31, -pz + k31) \
      + ops[2].projectMomentum(pz + k32, -pz - k32) \
      - ops[2].projectMomentum(pz - k32, -pz + k32)

  return [T1p_1_op, T1p_2_op, T1p_3_op]


def P1_A1_A1_1_0(ops, Ptot):
  if Ptot.psq != 1:
    print("Invalid P^2 given")
    return

  op = ops[0].projectMomentum(Ptot, P0)

  return [op]
  
def P1_A1_A1_2_1(ops, Ptot):

  if Ptot.psq != 1:
    print("Invalid P^2 given")
    return

  op = S.Zero
  if Ptot.x == 0:
    px = P([1,0,0])
    op += ops[0].projectMomentum(Ptot+px, -px) + ops[0].projectMomentum(Ptot-px, px)

  if Ptot.y == 0:
    py = P([0,1,0])
    op += ops[0].projectMomentum(Ptot+py, -py) + ops[0].projectMomentum(Ptot-py, py)

  if Ptot.z == 0:
    pz = P([0,0,1])
    op += ops[0].projectMomentum(Ptot+pz, -pz) + ops[0].projectMomentum(Ptot-pz, pz)

  return [op]


def P1_A1_E2_2_1(ops, Ptot):

  if Ptot.psq != 1:
    print("Invalid P^2 given")
    return

  op = S.Zero
  if Ptot.x == 0:
    p = P([1,0,0]) * Ptot
    op += ops[1].projectMomentum(Ptot+p, -p) - ops[1].projectMomentum(Ptot-p, p)

  if Ptot.y == 0:
    p = P([0,1,0]) * Ptot
    op += ops[2].projectMomentum(Ptot+p, -p) - ops[2].projectMomentum(Ptot-p, p)

  if Ptot.z == 0:
    p = P([0,0,1]) * Ptot
    op += ops[3].projectMomentum(Ptot+p, -p) - ops[3].projectMomentum(Ptot-p, p)

  return [op]


def P1_A2_A1_1_0(ops, Ptot):
  if Ptot.psq != 1:
    print("Invalid P^2 given")
    return

  if Ptot.x != 0:
    op = Ptot.x * ops[1].projectMomentum(Ptot, P0)
  elif Ptot.y != 0:
    op = Ptot.y * ops[2].projectMomentum(Ptot, P0)
  elif Ptot.z != 0:
    op = Ptot.z * ops[3].projectMomentum(Ptot, P0)

  return [op]

def P1_A2_A1_2_1(ops, Ptot):
  if Ptot.psq != 1:
    print("Invalid P^2 given")
    return

  px = P([1,0,0])
  py = P([0,1,0])
  pz = P([0,0,1])

  if Ptot.x != 0:
    op = Ptot.x * (ops[1].projectMomentum(Ptot+py, -py) + ops[1].projectMomentum(Ptot-py, py) \
        + ops[1].projectMomentum(Ptot+pz, -pz) + ops[1].projectMomentum(Ptot-pz, pz))

  elif Ptot.y != 0:
    op = Ptot.y * (ops[2].projectMomentum(Ptot+px, -px) + ops[2].projectMomentum(Ptot-px, px) \
        + ops[2].projectMomentum(Ptot+pz, -pz) + ops[2].projectMomentum(Ptot-pz, pz))

  elif Ptot.z != 0:
    op = Ptot.z * (ops[3].projectMomentum(Ptot+px, -px) + ops[3].projectMomentum(Ptot-px, px) \
        + ops[3].projectMomentum(Ptot+py, -py) + ops[3].projectMomentum(Ptot-py, py))

  return [op]


def P1_A2_E2_2_1(ops, Ptot):
  if Ptot.psq != 1:
    print("Invalid P^2 given")
    return

  px = P([1,0,0])
  py = P([0,1,0])
  pz = P([0,0,1])

  if Ptot.x != 0:
    op = ops[2].projectMomentum(Ptot+py, -py) - ops[2].projectMomentum(Ptot-py, py) \
        + ops[3].projectMomentum(Ptot+pz, -pz) - ops[3].projectMomentum(Ptot-pz, pz)

  elif Ptot.y != 0:
    op = ops[1].projectMomentum(Ptot+px, -px) - ops[1].projectMomentum(Ptot-px, px) \
        + ops[3].projectMomentum(Ptot+pz, -pz) - ops[3].projectMomentum(Ptot-pz, pz)

  elif Ptot.z != 0:
    op = ops[1].projectMomentum(Ptot+px, -px) - ops[1].projectMomentum(Ptot-px, px) \
        + ops[2].projectMomentum(Ptot+py, -py) - ops[2].projectMomentum(Ptot-py, py)

  return [op]

def P1_B1_B1_2_1(ops, Ptot):
  if Ptot.psq != 1:
    print("Invalid P^2 given")
    return

  if Ptot.x != 0:
    k1 = P([0,1,0])
    k2 = P([0,0,1])

  elif Ptot.y != 0:
    k1 = P([1,0,0])
    k2 = P([0,0,1])

  elif Ptot.z != 0:
    k1 = P([0,1,0])
    k2 = P([1,0,0])

  op = ops[0].projectMomentum(Ptot-k1, k1) + ops[0].projectMomentum(Ptot+k1, -k1) \
      - ops[0].projectMomentum(Ptot-k2, k2) - ops[0].projectMomentum(Ptot+k2, -k2)

  return [op]


def P1_B1_E2_2_1(ops, Ptot):
  if Ptot.psq != 1:
    print("Invalid P^2 given")
    return

  if Ptot.x != 0:
    k1 = P([0,1,0])
    k2 = P([0,0,1])

    k1s = 2
    k2s = 3

  elif Ptot.y != 0:
    k1 = P([1,0,0])
    k2 = P([0,0,1])

    k1s = 1
    k2s = 3

  elif Ptot.z != 0:
    k1 = P([0,1,0])
    k2 = P([1,0,0])

    k1s = 2
    k2s = 1

  op = ops[k1s].projectMomentum(Ptot + k1*Ptot, -k1*Ptot) - ops[k1s].projectMomentum(Ptot - k1*Ptot, k1*Ptot) \
      - ops[k2s].projectMomentum(Ptot + k2*Ptot, -k2*Ptot) + ops[k2s].projectMomentum(Ptot - k2*Ptot, k2*Ptot)

  return [op]


def P1_B2_B1_2_1(ops, Ptot):
  if Ptot.psq != 1:
    print("Invalid P^2 given")
    return

  if Ptot.x != 0:
    k1 = P([0,1,0])
    k2 = P([0,0,1])

    j = 1
    dj = Ptot.x

  elif Ptot.y != 0:
    k1 = P([1,0,0])
    k2 = P([0,0,1])

    j = 2
    dj = Ptot.y

  elif Ptot.z != 0:
    k1 = P([0,1,0])
    k2 = P([1,0,0])

    j = 3
    dj = Ptot.z

  op = dj * (ops[j].projectMomentum(Ptot - k1, k1) + ops[j].projectMomentum(Ptot + k1, -k1) \
      - ops[j].projectMomentum(Ptot - k2, k2) - ops[j].projectMomentum(Ptot + k2, -k2))

  return [op]

def P1_B2_E2_2_1(ops, Ptot):
  if Ptot.psq != 1:
    print("Invalid P^2 given")
    return

  if Ptot.x != 0:
    k1 = P([0,1,0])
    k2 = P([0,0,1])

    k1s = 2
    k2s = 3

  elif Ptot.y != 0:
    k1 = P([1,0,0])
    k2 = P([0,0,1])

    k1s = 1
    k2s = 3

  elif Ptot.z != 0:
    k1 = P([0,1,0])
    k2 = P([1,0,0])

    k1s = 2
    k2s = 1

  op = ops[k1s].projectMomentum(Ptot + k1, -k1) - ops[k1s].projectMomentum(Ptot - k1, k1) \
      - ops[k2s].projectMomentum(Ptot + k2, -k2) + ops[k2s].projectMomentum(Ptot - k2, k2)

  return [op]


def P1_E2_A1_1_0(ops, Ptot):
  if Ptot.psq != 1:
    print("Invalid P^2 given")
    return

  if Ptot.x != 0:
    op1 = ops[2].projectMomentum(Ptot, P0)
    op2 = ops[3].projectMomentum(Ptot, P0)

  elif Ptot.y != 0:
    op1 = ops[1].projectMomentum(Ptot, P0)
    op2 = ops[3].projectMomentum(Ptot, P0)

  elif Ptot.z != 0:
    op1 = ops[1].projectMomentum(Ptot, P0)
    op2 = ops[2].projectMomentum(Ptot, P0)

  return [op1, op2]


def P1_E2_E2_2_1_S0(ops, Ptot):
  if Ptot.psq != 1:
    print("Invalid P^2 given")
    return

  if Ptot.x != 0:
    k1 = P([0,1,0])
    k2 = P([0,0,1])

  elif Ptot.y != 0:
    k1 = P([1,0,0])
    k2 = P([0,0,1])

  elif Ptot.z != 0:
    k1 = P([1,0,0])
    k2 = P([0,1,0])

  op_1 = ops[0].projectMomentum(Ptot + k1*Ptot, -k1*Ptot) - ops[0].projectMomentum(Ptot - k1*Ptot, k1*Ptot)
  op_2 = ops[0].projectMomentum(Ptot + k2*Ptot, -k2*Ptot) - ops[0].projectMomentum(Ptot - k2*Ptot, k2*Ptot)

  return [op_1, op_2]


def P1_E2_A1_2_1(ops, Ptot):
  if Ptot.psq != 1:
    print("Invalid P^2 given")
    return

  if Ptot.x != 0:
    jk1_1 = P([1,0,0])
    jk1_2 = P([0,0,1])

    k1s = 2

    jk2_1 = P([1,0,0])
    jk2_2 = P([0,1,0])

    k2s = 3

  elif Ptot.y != 0:
    jk1_1 = P([0,1,0])
    jk1_2 = P([0,0,1])

    k1s = 1

    jk2_1 = P([1,0,0])
    jk2_2 = P([0,1,0])

    k2s = 3

  elif Ptot.z != 0:
    jk1_1 = P([0,0,1])
    jk1_2 = P([0,1,0])

    k1s = 1

    jk2_1 = P([1,0,0])
    jk2_2 = P([0,0,1])

    k2s = 2


  op_1 = ops[k1s].projectMomentum(Ptot + jk1_1, -jk1_1) + ops[k1s].projectMomentum(Ptot - jk1_1, jk1_1) \
      + ops[k1s].projectMomentum(Ptot + jk1_2, -jk1_2) + ops[k1s].projectMomentum(Ptot - jk1_2, jk1_2)

  op_2 = ops[k2s].projectMomentum(Ptot + jk2_1, -jk2_1) + ops[k2s].projectMomentum(Ptot - jk2_1, jk2_1) \
      + ops[k2s].projectMomentum(Ptot + jk2_2, -jk2_2) + ops[k2s].projectMomentum(Ptot - jk2_2, jk2_2)

  return [op_1, op_2]


def P1_E2_B1_2_1(ops, Ptot):
  if Ptot.psq != 1:
    print("Invalid P^2 given")
    return

  if Ptot.x != 0:
    k1 = P([0,1,0])
    k2 = P([0,0,1])

    k1s = 2
    k2s = 3

  elif Ptot.y != 0:
    k1 = P([1,0,0])
    k2 = P([0,0,1])

    k1s = 1
    k2s = 3

  elif Ptot.z != 0:
    k1 = P([1,0,0])
    k2 = P([0,1,0])

    k1s = 1
    k2s = 2

  op_1 = ops[k1s].projectMomentum(Ptot - k1, k1) + ops[k1s].projectMomentum(Ptot + k1, -k1) \
      - ops[k1s].projectMomentum(Ptot - k2, k2) - ops[k1s].projectMomentum(Ptot + k2, -k2)

  op_2 = -ops[k2s].projectMomentum(Ptot - k1, k1) - ops[k2s].projectMomentum(Ptot + k1, -k1) \
      + ops[k2s].projectMomentum(Ptot - k2, k2) + ops[k2s].projectMomentum(Ptot + k2, -k2)

  return [op_1, op_2]


def P1_E2_E2_2_1_S1(ops, Ptot):
  if Ptot.psq != 1:
    print("Invalid P^2 given")
    return

  if Ptot.x != 0:
    k1 = P([0,1,0])
    k2 = P([0,0,1])
    dj = Ptot.x

    op1 = dj * (ops[1].projectMomentum(Ptot + k1, -k1) - ops[1].projectMomentum(Ptot - k1, k1))
    op2 = dj * (ops[1].projectMomentum(Ptot + k2, -k2) - ops[1].projectMomentum(Ptot - k2, k2))

  elif Ptot.y != 0:
    k1 = P([1,0,0])
    k2 = P([0,0,1])
    dj = Ptot.y

    op1 = dj * (ops[2].projectMomentum(Ptot + k1, -k1) - ops[2].projectMomentum(Ptot - k1, k1))
    op2 = dj * (ops[2].projectMomentum(Ptot + k2, -k2) - ops[2].projectMomentum(Ptot - k2, k2))

  elif Ptot.z != 0:
    k1 = P([1,0,0])
    k2 = P([0,1,0])
    dj = Ptot.z

    op1 = dj * (ops[3].projectMomentum(Ptot + k1, -k1) - ops[3].projectMomentum(Ptot - k1, k1))
    op2 = dj * (ops[3].projectMomentum(Ptot + k2, -k2) - ops[3].projectMomentum(Ptot - k2, k2))

  return [op1, op2]


def P2_A1_A1_2_0(ops, Ptot):
  if Ptot.psq != 2:
    print("Invalid P^2 given")
    return

  op = ops[0].projectMomentum(Ptot, P0)

  return [op]


def P2_A1_A1_1_1(ops, Ptot):
  if Ptot.psq != 2:
    print("Invalid P^2 given")
    return

  if Ptot.x == 0:
    p1 = P([0,Ptot.y,0])
    p2 = P([0,0,Ptot.z])

  elif Ptot.y == 0:
    p1 = P([Ptot.x,0,0])
    p2 = P([0,0,Ptot.z])

  elif Ptot.z == 0:
    p1 = P([Ptot.x,0,0])
    p2 = P([0,Ptot.y,0])

  op = ops[0].projectMomentum(p1, p2)

  return [op]


def P2_A1_B1_1_1(ops, Ptot):
  if Ptot.psq != 2:
    print("Invalid P^2 given")
    return

  if Ptot.x == 0:
    p1 = P([0,Ptot.y,0])
    p2 = P([0,0,Ptot.z])

  elif Ptot.y == 0:
    p1 = P([Ptot.x,0,0])
    p2 = P([0,0,Ptot.z])

  elif Ptot.z == 0:
    p1 = P([Ptot.x,0,0])
    p2 = P([0,Ptot.y,0])

  p = p1*p2
  op = S.Zero
  if p.x != 0:
    op += p.x * ops[1]

  if p.y != 0:
    op += p.y * ops[2]

  if p.z != 0:
    op += p.z * ops[3]

  op = op.projectMomentum(p1, p2)

  return [op]


def P2_A2_A1_2_0(ops, Ptot):
  if Ptot.psq != 2:
    print("Invalid P^2 given")
    return

  op = S.Zero
  if Ptot.x != 0:
    op += Ptot.x * ops[1]

  if Ptot.y != 0:
    op += Ptot.y * ops[2]

  if Ptot.z != 0:
    op += Ptot.z * ops[3]

  op = op.projectMomentum(Ptot, P0)

  return [op]


def P2_A2_B1_1_1_Fs(ops, Ptot):
  if Ptot.psq != 2:
    print("Invalid P^2 given")
    return

  if Ptot.x == 0:
    p1 = P([0,Ptot.y,0])
    p2 = P([0,0,Ptot.z])

  elif Ptot.y == 0:
    p1 = P([Ptot.x,0,0])
    p2 = P([0,0,Ptot.z])

  elif Ptot.z == 0:
    p1 = P([Ptot.x,0,0])
    p2 = P([0,Ptot.y,0])

  p = p1 - p2
  op = S.Zero
  if p.x != 0:
    op += p.x * ops[1]

  if p.y != 0:
    op += p.y * ops[2]

  if p.z != 0:
    op += p.z * ops[3]

  op = op.projectMomentum(p1, p2)

  return [op]


def P2_A2_B1_1_1_Fa(ops, Ptot):
  if Ptot.psq != 2:
    print("Invalid P^2 given")
    return

  if Ptot.x == 0:
    p1 = P([0,Ptot.y,0])
    p2 = P([0,0,Ptot.z])

  elif Ptot.y == 0:
    p1 = P([Ptot.x,0,0])
    p2 = P([0,0,Ptot.z])

  elif Ptot.z == 0:
    p1 = P([Ptot.x,0,0])
    p2 = P([0,Ptot.y,0])

  p = p1 + p2
  op = S.Zero
  if p.x != 0:
    op += p.x * ops[1]

  if p.y != 0:
    op += p.y * ops[2]

  if p.z != 0:
    op += p.z * ops[3]

  op = op.projectMomentum(p1, p2)

  return [op]


def P2_B1_A1_2_0(ops, Ptot):
  if Ptot.psq != 2:
    print("Invalid P^2 given")
    return

  if Ptot.x == 0:
    p1 = P([0,Ptot.y,0])
    p2 = P([0,0,Ptot.z])

  elif Ptot.y == 0:
    p1 = P([Ptot.x,0,0])
    p2 = P([0,0,Ptot.z])

  elif Ptot.z == 0:
    p1 = P([Ptot.x,0,0])
    p2 = P([0,Ptot.y,0])

  p = p1*p2
  op = S.Zero
  if p.x != 0:
    op += p.x * ops[1]

  if p.y != 0:
    op += p.y * ops[2]

  if p.z != 0:
    op += p.z * ops[3]

  op = op.projectMomentum(Ptot, P0)

  return [op]

def P2_B1_A1_1_1(ops, Ptot):
  if Ptot.psq != 2:
    print("Invalid P^2 given")
    return

  if Ptot.x == 0:
    p1 = P([0,Ptot.y,0])
    p2 = P([0,0,Ptot.z])

  elif Ptot.y == 0:
    p1 = P([Ptot.x,0,0])
    p2 = P([0,0,Ptot.z])

  elif Ptot.z == 0:
    p1 = P([Ptot.x,0,0])
    p2 = P([0,Ptot.y,0])

  p = p1*p2
  op = S.Zero
  if p.x != 0:
    op += p.x * ops[1]

  if p.y != 0:
    op += p.y * ops[2]

  if p.z != 0:
    op += p.z * ops[3]

  op = op.projectMomentum(p1, p2)

  return [op]

def P2_B1_B1_1_1(ops, Ptot):
  if Ptot.psq != 2:
    print("Invalid P^2 given")
    return

  if Ptot.x == 0:
    p1 = P([0,Ptot.y,0])
    p2 = P([0,0,Ptot.z])

  elif Ptot.y == 0:
    p1 = P([Ptot.x,0,0])
    p2 = P([0,0,Ptot.z])

  elif Ptot.z == 0:
    p1 = P([Ptot.x,0,0])
    p2 = P([0,Ptot.y,0])

  op = ops[0].projectMomentum(p1,p2)

  return [op]


def P2_B2_A1_2_0(ops, Ptot):
  if Ptot.psq != 2:
    print("Invalid P^2 given")
    return

  if Ptot.x == 0:
    p1 = P([0,Ptot.y,0])
    p2 = P([0,0,Ptot.z])

  elif Ptot.y == 0:
    p1 = P([Ptot.x,0,0])
    p2 = P([0,0,Ptot.z])

  elif Ptot.z == 0:
    p1 = P([Ptot.x,0,0])
    p2 = P([0,Ptot.y,0])

  p = p1 - p2
  op = S.Zero
  if p.x != 0:
    op += p.x * ops[1]

  if p.y != 0:
    op += p.y * ops[2]

  if p.z != 0:
    op += p.z * ops[3]

  op = op.projectMomentum(Ptot, P0)

  return [op]

def P2_B2_B1_1_1(ops, Ptot):
  if Ptot.psq != 2:
    print("Invalid P^2 given")
    return

  if Ptot.x == 0:
    p1 = P([0,Ptot.y,0])
    p2 = P([0,0,Ptot.z])

  elif Ptot.y == 0:
    p1 = P([Ptot.x,0,0])
    p2 = P([0,0,Ptot.z])

  elif Ptot.z == 0:
    p1 = P([Ptot.x,0,0])
    p2 = P([0,Ptot.y,0])

  p = p1 + p2
  op = S.Zero
  if p.x != 0:
    op += p.x * ops[1]

  if p.y != 0:
    op += p.y * ops[2]

  if p.z != 0:
    op += p.z * ops[3]

  op = op.projectMomentum(p1, p2)

  return [op]


def P2_B2_A1_1_1(ops, Ptot):
  if Ptot.psq != 2:
    print("Invalid P^2 given")
    return

  if Ptot.x == 0:
    p1 = P([0,Ptot.y,0])
    p2 = P([0,0,Ptot.z])

  elif Ptot.y == 0:
    p1 = P([Ptot.x,0,0])
    p2 = P([0,0,Ptot.z])

  elif Ptot.z == 0:
    p1 = P([Ptot.x,0,0])
    p2 = P([0,Ptot.y,0])

  p = p1 - p2
  op = S.Zero
  if p.x != 0:
    op += p.x * ops[1]

  if p.y != 0:
    op += p.y * ops[2]

  if p.z != 0:
    op += p.z * ops[3]

  op = op.projectMomentum(p1, p2)

  return [op]


def P3_A1_A1_3_0(ops, Ptot):
  if Ptot.psq != 3:
    print("Invalid P^2 given")
    return

  op = ops[0].projectMomentum(Ptot, P0)

  return [op]


def P3_A1_A1_2_1(ops, Ptot):
  if Ptot.psq != 3:
    print("Invalid P^2 given")
    return

  p1 = Momentum([Ptot.x, 0, 0])
  p2 = Momentum([0, Ptot.y, 0])
  p3 = Momentum([0, 0, Ptot.z])

  op = ops[0].projectMomentum(Ptot - p1, p1) \
      + ops[0].projectMomentum(Ptot - p2, p2) \
      + ops[0].projectMomentum(Ptot - p3, p3)

  return [op]


def P3_A1_E2_2_1_old(ops, Ptot):
  if Ptot.psq != 3:
    print("Invalid P^2 given")
    return

  p1 = P([Ptot.x,0,0])
  p2 = P([0,Ptot.y,0])
  p3 = P([0,0,Ptot.z])

  op = p2.y*ops[2].projectMomentum(Ptot-p3,p3) \
      - p3.z*ops[3].projectMomentum(Ptot-p2,p2) \
      - p1.x*ops[1].projectMomentum(Ptot-p3,p3) \
      + p3.z*ops[3].projectMomentum(Ptot-p1,p1) \
      + p1.x*ops[1].projectMomentum(Ptot-p2,p2) \
      - p2.y*ops[2].projectMomentum(Ptot-p1,p1)

  return [op]

def P3_A1_E2_2_1(ops, Ptot):
  if Ptot.psq != 3:
    print("Invalid P^2 given")
    return

  p1 = P([Ptot.x,0,0])
  p2 = P([0,Ptot.y,0])
  p3 = P([0,0,Ptot.z])

  p1p2 = p1*p2
  coeff = p1p2.x * p3.x + p1p2.y * p3.y + p1p2.z * p3.z

  op = coeff * (p2.y*ops[2].projectMomentum(Ptot-p3,p3) \
      - p3.z*ops[3].projectMomentum(Ptot-p2,p2) \
      - p1.x*ops[1].projectMomentum(Ptot-p3,p3) \
      + p3.z*ops[3].projectMomentum(Ptot-p1,p1) \
      + p1.x*ops[1].projectMomentum(Ptot-p2,p2) \
      - p2.y*ops[2].projectMomentum(Ptot-p1,p1))

  return [op]

def P3_A2_A1_3_0(ops, Ptot):
  if Ptot.psq != 3:
    print("Invalid P^2 given")
    return

  op = Ptot.x * ops[1].projectMomentum(Ptot, P0) \
      + Ptot.y * ops[2].projectMomentum(Ptot, P0) \
      + Ptot.z * ops[3].projectMomentum(Ptot, P0)

  return [op]


def P3_A2_A1_2_1(ops, Ptot):
  if Ptot.psq != 3:
    print("Invalid P^2 given")
    return

  p1 = Momentum([Ptot.x, 0, 0])
  p2 = Momentum([0, Ptot.y, 0])
  p3 = Momentum([0, 0, Ptot.z])

  op = Ptot.x * ops[1].projectMomentum(Ptot - p1, p1) \
      + Ptot.x * ops[1].projectMomentum(Ptot - p2, p2) \
      + Ptot.x * ops[1].projectMomentum(Ptot - p3, p3) \
      + Ptot.y * ops[2].projectMomentum(Ptot - p1, p1) \
      + Ptot.y * ops[2].projectMomentum(Ptot - p2, p2) \
      + Ptot.y * ops[2].projectMomentum(Ptot - p3, p3) \
      + Ptot.z * ops[3].projectMomentum(Ptot - p1, p1) \
      + Ptot.z * ops[3].projectMomentum(Ptot - p2, p2) \
      + Ptot.z * ops[3].projectMomentum(Ptot - p3, p3)

  return [op]


def P3_A2_E2_2_1(ops, Ptot):
  if Ptot.psq != 3:
    print("Invalid P^2 given")
    return

  p1 = Momentum([Ptot.x, 0, 0])
  p2 = Momentum([0, Ptot.y, 0])
  p3 = Momentum([0, 0, Ptot.z])


  op = (p1*3 - p1 - p2 - p3).x * ops[1].projectMomentum(Ptot-p1, p1) \
      + (p1*3 - p1 - p2 - p3).y * ops[2].projectMomentum(Ptot-p1, p1) \
      + (p1*3 - p1 - p2 - p3).z * ops[3].projectMomentum(Ptot-p1, p1) \
      + (p2*3 - p1 - p2 - p3).x * ops[1].projectMomentum(Ptot-p2, p2) \
      + (p2*3 - p1 - p2 - p3).y * ops[2].projectMomentum(Ptot-p2, p2) \
      + (p2*3 - p1 - p2 - p3).z * ops[3].projectMomentum(Ptot-p2, p2) \
      + (p3*3 - p1 - p2 - p3).x * ops[1].projectMomentum(Ptot-p3, p3) \
      + (p3*3 - p1 - p2 - p3).y * ops[2].projectMomentum(Ptot-p3, p3) \
      + (p3*3 - p1 - p2 - p3).z * ops[3].projectMomentum(Ptot-p3, p3)

  return [op]

def P4_A1_A1_1_1(ops, Ptot):
  if Ptot.psq != 4:
    print("Invalid P^2 given")
    return

  op = ops[0].projectMomentum(Ptot/2, Ptot/2)

  return [op]

def P4_A2_A1_1_1(ops, Ptot):
  if Ptot.psq != 4:
    print("Invalid P^2 given")
    return

  if Ptot.x != 0:
    op = Ptot.x * ops[1]

  elif Ptot.y != 0:
    op = Ptot.y * ops[2]

  elif Ptot.z != 0:
    op = Ptot.z * ops[3]

  op = op.projectMomentum(Ptot/2, Ptot/2)

  return [op]

def P4_E1_A1_1_1(ops, Ptot):
  if Ptot.psq != 4:
    print("Invalid P^2 given")
    return

  if Ptot.x != 0:
    op1 = ops[2].projectMomentum(Ptot/2, Ptot/2)
    op2 = ops[3].projectMomentum(Ptot/2, Ptot/2)

  elif Ptot.y != 0:
    op1 = ops[1].projectMomentum(Ptot/2, Ptot/2)
    op2 = ops[3].projectMomentum(Ptot/2, Ptot/2)

  elif Ptot.z != 0:
    op1 = ops[1].projectMomentum(Ptot/2, Ptot/2)
    op2 = ops[2].projectMomentum(Ptot/2, Ptot/2)

  return [op1, op2]

def test_irrep(ops, irrep, op_name, use_generators):
  print("TESTING: {}".format(op_name))

  op_rep = OperatorRepresentation(*ops)
  t_irrep = op_rep.littleGroupContents(True, use_generators)

  if irrep != t_irrep:
    print("\tFAILED - expected {}, transforms as {}".format(irrep, t_irrep))
    return 0
  else:
    print("\tPASSED - expected {}, transforms as {}".format(irrep, t_irrep))
    return 1


g = Gamma()

C5p = Array(g.chargeConj * g.five * g.parityPlus)
C1p = Array(g.chargeConj * g.one * g.parityPlus)
C2p = Array(g.chargeConj * g.two * g.parityPlus)
C3p = Array(g.chargeConj * g.three * g.parityPlus)

i = DiracIdx('i')
j = DiracIdx('j')

def flavor_term(fl1_i, fl2_i):

  fl1_0_j = fl1_i * C5p[i,j]
  fl1_1_j = fl1_i * C1p[i,j]
  fl1_2_j = fl1_i * C2p[i,j]
  fl1_3_j = fl1_i * C3p[i,j]

  fl1_0_0_op = Operator(fl1_0_j.subs(j,0))
  fl1_0_1_op = Operator(fl1_0_j.subs(j,1))
  fl1_0_2_op = Operator(fl1_0_j.subs(j,2))
  fl1_0_3_op = Operator(fl1_0_j.subs(j,3))

  fl1_1_0_op = Operator(fl1_1_j.subs(j,0))
  fl1_1_1_op = Operator(fl1_1_j.subs(j,1))
  fl1_1_2_op = Operator(fl1_1_j.subs(j,2))
  fl1_1_3_op = Operator(fl1_1_j.subs(j,3))

  fl1_2_0_op = Operator(fl1_2_j.subs(j,0))
  fl1_2_1_op = Operator(fl1_2_j.subs(j,1))
  fl1_2_2_op = Operator(fl1_2_j.subs(j,2))
  fl1_2_3_op = Operator(fl1_2_j.subs(j,3))

  fl1_3_0_op = Operator(fl1_3_j.subs(j,0))
  fl1_3_1_op = Operator(fl1_3_j.subs(j,1))
  fl1_3_2_op = Operator(fl1_3_j.subs(j,2))
  fl1_3_3_op = Operator(fl1_3_j.subs(j,3))


  fl2_0_op = Operator(fl2_i.subs(i,0))
  fl2_1_op = Operator(fl2_i.subs(i,1))
  fl2_2_op = Operator(fl2_i.subs(i,2))
  fl2_3_op = Operator(fl2_i.subs(i,3))


  op_0 = fl1_0_0_op*fl2_0_op + fl1_0_1_op*fl2_1_op + fl1_0_2_op*fl2_2_op + fl1_0_3_op*fl2_3_op
  op_1 = fl1_1_0_op*fl2_0_op + fl1_1_1_op*fl2_1_op + fl1_1_2_op*fl2_2_op + fl1_1_3_op*fl2_3_op
  op_2 = fl1_2_0_op*fl2_0_op + fl1_2_1_op*fl2_1_op + fl1_2_2_op*fl2_2_op + fl1_2_3_op*fl2_3_op
  op_3 = fl1_3_0_op*fl2_0_op + fl1_3_1_op*fl2_1_op + fl1_3_2_op*fl2_2_op + fl1_3_3_op*fl2_3_op

  return [op_0, op_1, op_2, op_3]


def Baryon(q1, q2, q3, alpha):
  aB = ColorIdx('aB')
  bB = ColorIdx('bB')
  cB = ColorIdx('cB')

  iB = DiracIdx('iB')
  jB = DiracIdx('jB')

  baryon = Eijk(aB, bB, cB) * q2[aB,iB] * C5p[iB,jB] * q3[bB,jB] * q1[cB,alpha]
  return baryon


if __name__ == "__main__":
  main()
