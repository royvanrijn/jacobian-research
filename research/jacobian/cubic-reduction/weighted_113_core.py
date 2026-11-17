#!/usr/bin/env python3
"""Small remaining quadratic-graph core after the exact highest-v case split."""
from pathlib import Path
import hashlib,json
import sympy as s
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
u,beta,P,T=s.symbols('u beta P T');A=-3-beta;alpha=2*beta+4;G0=A*u-beta;G1=alpha*u+beta;H1=-(G1/beta)**2;H0=1+(2*alpha/beta-1)*u+(alpha**2/beta**2-2*alpha/beta+P)*u**2+T*u**3
E0=s.expand(G0*s.diff(H0,u)-3*A*H0-1);E1=s.expand(2*G1*s.diff(H0,u)-3*alpha*H0+G0*s.diff(H1,u)-4*A*H1);V=(T,P,beta);eqs=[s.Poly(s.fraction(s.cancel(c))[0],*V).monic().as_expr() for f in [E0,E1] for e,c in s.Poly(f,u).terms() if c];out=HERE/'weighted-113/core';out.mkdir(exist_ok=True);payload=dict(source_sha256=sha(Path(__file__)),variables=list(map(str,V)),equations=[[[list(e),str(c)] for e,c in s.Poly(f,*V).terms() if c] for f in eqs],scope='Necessary remaining core atK=H2=0, G1=alpha*u+beta, beta!=0; source normalizations forceA=-3-beta,alpha=2beta+4 and the displayedH0,H1. Denominators are beta powers only; exact unit excludes every remaining point. All branch/normalization derivations need independent replay and written coverage proof.');(out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');print('CORE',len(eqs))
