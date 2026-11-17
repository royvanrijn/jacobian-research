#!/usr/bin/env python3
"""Reproduce the complete9-parameter target-locus inputs;retained logs suffice for replay.
This preparation script makes no solver calls. The saved model.sing and certify.sing
are the exact Singular queries used for discovery and seven membership certificates.
"""
from pathlib import Path
import sympy as s
import json
B=Path(__file__).resolve().parent;O=B/'raw13-target-field/target-locus'
def main():
 O.mkdir(exist_ok=True);X=s.symbols('x y z');x,y,z=X;U=s.symbols('u0:3');V=s.symbols('v0:3');W=s.symbols('w0:3');vs=W+V+U;Q=x*x;S=y**3+27*x*z*z-9*y*z;T=y*y-4*x*S;R=[X[i]+U[i]*Q+V[i]*S+W[i]*T for i in range(3)];J=s.Poly(s.expand(s.det(s.Matrix(R).jacobian(X))-1),*X);eqs=[co for e,co in J.terms() if co]
 payload={'source_variables':list(map(str,X)),'parameters':list(map(str,vs)),'target_basis':[str(Q),str(S),str(T)],'target_map':list(map(str,R)),'equations':[[{'exponent':list(e),'coefficient':str(co)} for e,co in s.Poly(f,*vs).terms()] for f in eqs],'scope':'Complete tangent-identity target maps in the seven-dimensional degree8 function-field intersection;constant-Jacobian locus,all coefficient strata.'}
 (O/'input.json').write_text(json.dumps(payload,indent=2)+'\n');prefix='ring r=0,('+','.join(map(str,vs))+'),dp;\nshort=0;\nideal I='+','.join(str(s.expand(f)).replace('**','^') for f in eqs)+';\n'
 (O/'model.sing').write_text(prefix+'ideal G=std(I);\nprint("BASIS_BEGIN");print(G);print("BASIS_END");\nquit;\n')
 (O/'certify.sing').write_text(prefix+'ideal K=u0,v0,v1,v2,w0,w1,w2;\nmatrix C=lift(I,K);\nint i,j;\nfor(j=1;j<=size(K);j++){print("CERT_BEGIN|"+string(K[j]));for(i=1;i<=size(I);i++){print(string(i)+"|"+string(C[i,j]));}print("CERT_END");}\nquit;\n')
 print('PREPARED27 determinant generators and exact retained certificate query;no solver launched')
if __name__=='__main__':main()
