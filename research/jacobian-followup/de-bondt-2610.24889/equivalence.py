"""Exact monomial conjugacy between the repository K12 and paper section 4.5."""
from __future__ import annotations
import json
import sys
import sympy as s
from verify import X, ROOT, coefficient_rank


def main():
    oldscripts=ROOT.parents[1]/'archive/non-elliptic/scripts'
    sys.path.insert(0,str(oldscripts))
    from audit_macfarlane_g20_dimension_reduction import build_maps
    old=build_maps();oldx=old['x'];z=s.symbols('z1:13')
    sub=dict(zip(oldx[:12],z));sub[oldx[12]]=-z[1]**2
    transformed=[s.expand(f.subs(sub,simultaneous=True)) for f in old['F13']]
    transformed[3]=s.expand(transformed[3]-transformed[7]**2)
    K=transformed[:12]
    a,b,c=X[:3];v=1+2*a*b
    f=[a-3*a*a*b-a**3*c,3*v*v*a*c+(9*v*v-6*v-2)*b,
       v**3*c+(6*v*v+2*v)*b*b,
       X[3]+a*a*b*b,X[4]+a*a*b,X[5]+a*b*c+3*b*b,
       X[6]+a*a,X[7]+a*b,X[8]+a*c,
       X[9]+3*b*c-2*b*X[5],X[10]+a*c-a*X[5]-c*X[4],
       X[11]+7*b*b-2*c*X[3]-a*X[9]]
    J=[f[0]+f[6]*f[8],f[1]-12*f[4]*f[5]-12*f[7]*f[10],
       f[2]-8*f[3]*f[5]-4*f[4]*f[9]-4*f[7]*f[11],f[3]-f[7]**2]+f[4:]
    J=list(map(s.expand,J))
    permutation=[0,1,2,3,5,4,9,7,10,11,8,6]
    scales=[1,2,2,4,4,2,4,2,12,4,2,s.Rational(-1,2)]
    Tx=[scale*X[j] for j,scale in zip(permutation,scales)]
    TJ=[scale*J[j] for j,scale in zip(permutation,scales)]
    assert all(s.expand(k.subs(dict(zip(z,Tx)),simultaneous=True)-tj)==0 for k,tj in zip(K,TJ))
    T=s.Matrix(Tx).jacobian(X[:12])
    assert T.det()!=0
    result={'identity':'K(Tx)=T J(x); J is section 4.5, not the five-cubic section 4.6 map',
            'T_permutation_1_based':[i+1 for i in permutation],
            'T_scales':[str(c) for c in scales],
            'T_determinant':str(T.det()),'polynomial_identities':12}
    (ROOT/'equivalence.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
