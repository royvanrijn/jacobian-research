#!/usr/bin/env python3
"""Exact cover/shear, whole-fibre and omission controls; not formalization."""
from pathlib import Path
from itertools import product,combinations
from math import gcd
import hashlib,json,platform,time
import sympy as s
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
checks=[]
def ck(label,value):
    assert value,label
    checks.append(label)
def dim(gb,n):
    leading=[p.LM(order=gb.order).exponents for p in gb.polys]
    if (0,)*n in leading:return 0
    bounds=[]
    for i in range(n):
        v=[m[i] for m in leading if m[i] and all(a==0 for j,a in enumerate(m) if j!=i)]
        assert v,'positive dimensional quotient'
        bounds.append(min(v))
    return sum(not any(all(a>=b for a,b in zip(v,m)) for m in leading) for v in product(*(range(x) for x in bounds)))
def main():
    started=time.monotonic();u,v,w,t,e,U,x,y,z=s.symbols('u v w t e U x y z')
    fibres=[];etale=[];examples=[]
    for k in range(2,6):
        phi=s.cyclotomic_poly(k,e)
        def red(expr):return s.Poly(s.expand(expr),e,domain=s.QQ.poly_ring(u,v,w,t,U,x,y,z)).rem(s.Poly(phi,e,domain=s.QQ.poly_ring(u,v,w,t,U,x,y,z))).as_expr().expand()
        R=red(s.div((1+(e-1)*t)**k-1,t*(t-1),t)[0]);j=[w,v*R.subs(t,-u*v),1-(e-1)*u*v]
        G=u+u*u*v-w**k;B=x**k*y-z**k+1
        source=s.groebner([G,phi],u,v,w,e,domain=s.QQ)
        # Q=0 is the old surjective sample; Q=1 and Q_A move its omitted lines.
        for label,Q in [('zero',s.Integer(0)),('one',s.Integer(1)),('two_omissions',s.cancel((e*(1-U)*(1-U/2)-e)/U))]:
            hh=e+U*Q;bb=red(s.div(hh**k-1,U,U)[0]);shift=x**k*Q.subs(U,x**k)
            alpha=[x,y+s.div((z+shift)**k-z**k,x**k,x)[0],z+shift]
            ck(f'k{k} Q{label} polynomial shear preserves cover',red(B.subs({x:alpha[0],y:alpha[1],z:alpha[2]},simultaneous=True)-B)==0)
            minus=[x,y+s.div((z-shift)**k-z**k,x**k,x)[0],z-shift]
            ck(f'k{k} Q{label} polynomial shear inverse',all(red(a.subs({x:minus[0],y:minus[1],z:minus[2]},simultaneous=True)-b)==0 for a,b in zip(alpha,[x,y,z])))
            Aj=[red(q.subs({x:j[0],y:j[1],z:j[2]},simultaneous=True)) for q in alpha]
            F=[Aj[0]**k,Aj[1],Aj[0]*Aj[2]]
            ck(f'k{k} Q{label} entire source maps into target',source.reduce(s.expand(F[0]+F[0]**2*F[1]-F[2]**k))[1]==0)
            ck(f'k{k} Q{label} actual cyclic source-target equivariance',all(red(F[i].subs(w,e*w)-(e*F[i] if i==2 else F[i]))==0 for i in range(3)))
            ss=u+u*u*v;z0=1-(e-1)*u*v
            vpoly=s.Poly(F[1],w);assert all(m[0]%k==0 for m,c in vpoly.terms())
            barV=sum(c*ss**(m[0]//k) for m,c in vpoly.terms())
            jac=s.Matrix([ss,barV]).jacobian([u,v]).det()
            ck(f'k{k} Q{label} coarse quotient Jacobian identity',red(jac-k*(1-e)*(z0+ss*Q.subs(U,ss))**(k-1))==0)

            ck(f'k{k} Q{label} b polynomial identity',red(U*bb-(hh**k-1))==0)
            # Every target below satisfies the target equation over the exact cyclotomic field.
            b0=bb.subs(U,0)
            cases=[('vertical rank one',(0,b0,0),1),('vertical rank two',(0,b0+1,0),2)]
            aa=2**k;ha=hh.subs(U,aa);ba=bb.subs(U,aa)
            cases.append(('partial omitted curve',(aa,ba,2*ha),k-1))
            # W=0 off the prescribed omission roots retains all k sheets.
            cases.append(('complete cover fibre',(aa,-s.Rational(1,aa),0),k))
            roots=[]
            if label=='one':roots=[-e]
            elif label=='two_omissions':roots=[s.Integer(1),s.Integer(2)]
            for root in roots:
                # -1/root is reduced using the cyclotomic relation before forming ideals.
                inv=s.invert(s.Poly(root,e),s.Poly(phi,e)).as_expr()
                cases.append(('empty prescribed point',(root,-inv,0),0))
            for name,point,expected in cases:
                ck(f'k{k} Q{label} {name} target lies on surface',red(point[0]+point[0]**2*point[1]-point[2]**k)==0)
                gb=s.groebner([G,phi]+[s.expand(F[i]-point[i]) for i in range(3)],u,v,w,e,order='grevlex',domain=s.QQ)
                actual=dim(gb,4);degree=int(s.degree(phi,e))
                if name=='complete cover fibre':
                    q0=red(Q.subs(U,aa));h0=red(hh.subs(U,aa));invh=s.invert(s.Poly(h0,e),s.Poly(phi,e)).as_expr();inve=s.invert(s.Poly(e-1,e),s.Poly(phi,e)).as_expr()
                    uc=red(aa*(e-1)*invh);vc=red((1+aa*q0)*h0*inve**2/aa)
                    ck(f'k{k} Q{label} full marked fibre has fixed u,v and w^k=c^k',all(gb.reduce(s.expand(q))[1]==0 for q in [u-uc,v-vc,w**k-aa]))

                ck(f'k{k} Q{label} {name} complete quotient',actual==expected*degree)
                fibres.append({'degree':k,'Q':label,'stratum':name,'target':[str(a) for a in point],'whole_geometric_rank':expected,'Q_quotient_dimension':actual,'cyclotomic_degree':degree})
            if label=='one' and k in [2,3]:
                matrix=s.Matrix([[s.diff(G,a) for a in [u,v,w]]]).col_join(s.Matrix(F).jacobian([u,v,w]))
                minors=[s.expand(matrix[list(rows),:].det()) for rows in combinations(range(4),3)]
                ck(f'k{k} Q1 whole source etaleness ideal',list(s.groebner([G,phi]+minors,u,v,w,e))==[1]);etale.append(k)
            examples.append({'degree':k,'Q':label,'omission_count':len(roots)})
    # Arbitrarily many rational omissions are constructed by an all-m formula;
    # these controls check whole empty-fibre ideals, not found rational roots.
    for m in range(1,5):
        hh=-s.prod(1-U/a for a in range(1,m+1));Q=s.div(hh+1,U,U)[0]
        ck(f'prescribed m{m} polynomial interpolation identity',s.expand(-1+U*Q-hh)==0)
        z0=1+2*u*v;q=Q.subs(U,w*w)
        F=[w*w,4*v+2*z0*q+w*w*q*q,w*(z0+w*w*q)]
        for a in range(1,m+1):
            gb=s.groebner([u+u*u*v-w*w,F[0]-a,F[1]+s.Rational(1,a),F[2]],u,v,w)
            ck(f'prescribed m{m} actual empty fibre at a{a}',list(gb)==[1])
        ck(f'prescribed m{m} distinct roots and omission components',s.Poly(hh,U).degree()==m and s.gcd(s.Poly(hh,U),s.Poly(s.diff(hh,U),U)).degree()==0)
    # All nonzero quadratic torsor parameters at once over Q[a,a^{-1}].
    a,r=s.symbols('a r');z0=1+2*u*v
    F=[r*z*z,4*v+2*z0+r*z*z,z*(z0+r*z*z)]
    ideal=[a*(u+u*u*v)-z*z,a*r-1,F[0]-4,F[1]+s.Rational(1,4),F[2]]
    whole=s.groebner(ideal,u,v,z,r,a,domain=s.QQ)
    expected=s.groebner([u+s.Rational(8,3),v-s.Rational(15,16),z*z-4*a,a*r-1],u,v,z,r,a,domain=s.QQ)
    ck('whole quadratic twisted family has exactly the torsor marked fibre',whole==expected)
    ck('whole cyclic fixed scheme has two disjoint components',list(s.groebner([u,1+u*v],u,v))==[1])
    # Exact intersection of the two moved boundary lines in the quadratic cover.
    gb=s.groebner([x*x*y-z*z+1,z+1-x*x,y+2-x*x,z-1+x*x],x,y,z)
    ck('quadratic moved boundary intersection has two points',dim(gb,3)==2)
    ck('quadratic moved boundary maps to the unique omitted point',all(gb.reduce(q)[1]==0 for q in [x*x-1,y+1,x*z]))
    for n in range(2,66,2):
        perm=tuple(i+1 if i%2==0 else i-1 for i in range(n))
        ck(f'even symmetric rank{n} explicit 2-derangement',all(perm[i]!=i and perm[perm[i]]==i for i in range(n)))
    for k in [2,3]:
        for ell in [2,3,5,7]:
            ck(f'cyclic H1 torsion k{k} mod{ell}',(s.Matrix([[k]]).rank(iszerofunc=lambda a:a%ell==0)==1)==(gcd(k,ell)==1))
    manifest=json.loads((HERE/'sources.json').read_text())
    for src in manifest['local_sources']:ck('local source pin '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    for src in manifest['external_sources']:ck('primary source pin '+src['path'],hashlib.sha256((HERE/src['path']).read_bytes()).hexdigest()==src['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,'complete_fibre_controls':fibres,'deformation_cases':examples,'global_etaleness_certificates':etale,'python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-started,3),'sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']},'assurance':'Written all-characteristic unit/Picard Galois lifting, characteristic-zero cohomological consequences, pluricanonical primitive-monodromy stable atomicity, canonical cyclic torsor descent and all-k cover image proofs. Exact bounded controls certify polynomial shears, complete fibre ideals including empty fibres, two global etaleness ideals and finite omission/component examples. They do not formalize purity, compactification, etale cohomology, arbitrary-m topology or global equivariant gluing. Construction credited to Dubouloz-Palka. No independent whole-proof replay, formal verification, external review or novelty certification.'}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({'result':'pass','checks':len(checks),'complete_fibres':len(fibres),'global_etaleness':etale,'seconds':receipt['runtime_seconds']}))
if __name__=='__main__':main()
