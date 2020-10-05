from sympy import KroneckerDelta, Eijk, Array, sqrt

from context import operators

from operators.operators import QuarkField, AntiQuarkField, ColorIdx, DiracIdx
from operators.tensors import Gamma

g = Gamma()

u = QuarkField.create('u')
ubar = AntiQuarkField.create('u')
d = QuarkField.create('d')
dbar = AntiQuarkField.create('d')
up = QuarkField.create('up')
upbar = AntiQuarkField.create('up')
dp = QuarkField.create('dp')
dpbar = AntiQuarkField.create('dp')
s = QuarkField.create('s')
sbar = AntiQuarkField.create('s')

a = ColorIdx('a')
b = ColorIdx('b')
c = ColorIdx('c')
i = DiracIdx('i')
j = DiracIdx('j')
k = DiracIdx('k')

pion_p = KroneckerDelta(a, b)*dbar[a,i]*Array(g.five)[i,j]*u[b,j]
pion_0 = 1/sqrt(2)*KroneckerDelta(a, b)*(dbar[a,i]*Array(g.five)[i,j]*d[b,j] - ubar[a,i]*Array(g.five)[i,j]*u[b,j])
pion_m = -KroneckerDelta(a, b)*ubar[a,i]*Array(g.five)[i,j]*d[b,j]

pionp_p = KroneckerDelta(a, b)*dpbar[a,i]*Array(g.five)[i,j]*up[b,j]
pionp_0 = 1/sqrt(2)*KroneckerDelta(a, b)*(dpbar[a,i]*Array(g.five)[i,j]*dp[b,j] - upbar[a,i]*Array(g.five)[i,j]*up[b,j])
pionp_m = -KroneckerDelta(a, b)*upbar[a,i]*Array(g.five)[i,j]*dp[b,j]

kaon_p = KroneckerDelta(a, b)*sbar[a,i]*Array(g.five)[i,j]*u[b,j]
kaon_0 = KroneckerDelta(a, b)*sbar[a,i]*Array(g.five)[i,j]*d[b,j]

Cg5 = Array(g.chargeConj*g.five)

#nucleon_1 = Eijk(a,b,c)*u[a,0]*u[b,j]*Cg5[j,k]*d[c,k]
#nucleon_2 = Eijk(a,b,c)*u[a,1]*u[b,j]*Cg5[j,k]*d[c,k]

nucleon_p_1 = Eijk(a,b,c)*(0.5*u[a,0]*u[b,1]*d[c,0] - u[a,0]*u[b,0]*d[c,1] + 0.5*u[a,1]*u[b,0]*d[c,0])
nucleon_p_2 = Eijk(a,b,c)*(-0.5*u[a,1]*u[b,0]*d[c,1] + u[a,1]*u[b,1]*d[c,0] - 0.5*u[a,0]*u[b,1]*d[c,1])

nucleon_p = [nucleon_p_1, nucleon_p_2]

nucleon_m_1 = -Eijk(a,b,c)*(0.5*d[a,0]*d[b,1]*u[c,0] - d[a,0]*d[b,0]*u[c,1] + 0.5*d[a,1]*d[b,0]*u[c,0])
nucleon_m_2 = -Eijk(a,b,c)*(-0.5*d[a,1]*d[b,0]*u[c,1] + d[a,1]*d[b,1]*u[c,0] - 0.5*d[a,0]*d[b,1]*u[c,1])

nucleon_m = [nucleon_m_1, nucleon_m_2]

