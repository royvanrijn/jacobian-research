"""Exact finite-group verification of the six-point ramified-3 calculation.

All operations are in integers modulo 3; there is no numerical tolerance.
"""
from itertools import product
import json
from pathlib import Path

P = 3

def radd(x,y):return tuple((a+b)%P for a,b in zip(x,y))
def rneg(x):return tuple(-a%P for a in x)
def rmul(x,y):
    a,b=x;c,d=y
    return ((a*c-b*d)%P,(a*d+b*c-b*d)%P)
def rbar(x):return ((x[0]-x[1])%P,-x[1]%P)
def rnorm(x):return (x[0]**2-x[0]*x[1]+x[1]**2)%P
def qmul(x,y):
    z,w=x[:2],x[2:];u,v=y[:2],y[2:]
    return radd(rmul(z,u),rneg(rmul(w,rbar(v))))+radd(rmul(z,v),rmul(w,rbar(u)))
def qnorm(x):return (rnorm(x[:2])+rnorm(x[2:]))%P
def canonical(x):
    t=next(a for a in x if a)
    return tuple(a*pow(t,-1,P)%P for a in x)
def mul(x,y):return canonical(qmul(x,y))
def inverse(x):
    t=pow(qnorm(x),-1,P)
    return tuple(a*t%P for a in rbar(x[:2])+rneg(x[2:]))

I=(1,0,0,0);W=(0,1,0,0);J=(0,0,1,0)
G={canonical(x) for x in product(range(P),repeat=4) if qnorm(x)}
W2=qmul(W,W)
Gamma={canonical(x) for x in [I,W,W2,J,qmul(W,J),qmul(W2,J)]}
assert len(G)==36 and len(Gamma)==6
def coset(x):return min(mul(g,x) for g in Gamma)
Y=sorted({coset(x) for x in G});index={x:i for i,x in enumerate(Y)}
assert len(Y)==6
chi=[1 if qnorm(x)==1 else -1 for x in Y]
assert chi.count(1)==chi.count(-1)==3
alpha=(1,0,1,0)
beta=qmul(qmul(alpha,W),inverse(alpha))
assert beta==(1,0,1,2)
assert beta not in Gamma and canonical(qmul(beta,beta)) not in Gamma
assert len({coset(qmul(alpha,w)) for w in [I,W,W2]})==3
assert {x for x in G if x[:2]==(0,0) or x[2:]==(0,0)}==Gamma
generators=[(1,0)+w[:2] for w in [I,W,W2]]
T=[[0]*6 for _ in range(6)]
for i,g in enumerate(Y):
    for a in generators:T[i][index[coset(mul(a,g))]]+=1
assert all(T[i][j]==int(chi[i]!=chi[j]) for i in range(6) for j in range(6))
result={
    'group_size':len(G),'unit_image_size':len(Gamma),'quotient_size':len(Y),
    'quotient_representatives':Y,'norm_characters':chi,'T2_matrix':T,
    'T2_is_complete_bipartite_exactly':True,
    'gamma_transitive_on_opposite_norm_class':True,
    'bad_ramified_residue_classes_equal_gamma':True,
    'alpha_omega_alpha_inverse_mod3':beta,
    'all_arithmetic':'exact integers modulo 3',
}
print(json.dumps(result,indent=2))
Path(__file__).with_name('mixed_spectral_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
