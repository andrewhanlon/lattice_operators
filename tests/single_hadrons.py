from sympy import KroneckerDelta, Array, sqrt

from context import operators

from operators.operators import QuarkField, AntiQuarkField, ColorIdx, DiracIdx
from operators.tensors import Gamma

g = Gamma()

u = QuarkField.create('u')
ubar = AntiQuarkField.create('u')
d = QuarkField.create('d')
dbar = AntiQuarkField.create('d')
s = QuarkField.create('s')
sbar = AntiQuarkField.create('s')

a = ColorIdx('a')
b = ColorIdx('b')
i = DiracIdx('i')
j = DiracIdx('j')

pion_p = KroneckerDelta(a, b)*dbar[a,i]*Array(g.five)[i,j]*u[b,j]
pion_0 = 1/sqrt(2)*KroneckerDelta(a, b)*(dbar[a,i]*Array(g.five)[i,j]*d[b,j] - ubar[a,i]*Array(g.five)[i,j]*u[b,j])
pion_m = -KroneckerDelta(a, b)*ubar[a,i]*Array(g.five)[i,j]*d[b,j]

kaon_p = KroneckerDelta(a, b)*sbar[a,i]*Array(g.five)[i,j]*u[b,j]
kaon_0 = KroneckerDelta(a, b)*sbar[a,i]*Array(g.five)[i,j]*d[b,j]
