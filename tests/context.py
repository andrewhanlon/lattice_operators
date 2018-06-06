import os
import sys

this_file = os.path.abspath(__file__)
operators_dir = os.path.join(os.path.dirname(this_file), "..")
operators_dir = os.path.normpath(operators_dir)
sys.path.insert(0, operators_dir)

import operators
