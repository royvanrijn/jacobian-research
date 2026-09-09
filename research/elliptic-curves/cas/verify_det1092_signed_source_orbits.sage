#!/usr/bin/env sage-python
"""Rebuild signed source maps and verify saved Hensel/Gauss root certificates.

No factorization, root-finding, lattice-reduction call, or elliptic point
search. The certificate inequalities themselves prove completeness. <=25s.
"""
import hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,matrix,vector,gcd,lcm,Zmod
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_signed_source_orbits_v1';PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
PENCIL=ART/'det1092_norm8_seed_cover_v2/generic.json'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
protocol=read(OUT/'protocol.json');d=read(PARENT);p=read(PENCIL)
for path,value in protocol['inputs'].items():assert sha(ROOT/path)==value
R=PolynomialRing(QQ,'u');u=R.gen();K=R.fraction_field()
def dec(a):return K(R(a['numerator']))/R(a['denominator'])
aa=[dec(x) for x in d['a_invariants']];a1,a2,a3,a4,a6=aa;b2=a1*a1+4*a2
h=R(p['pole_h']);shift=R(p['shift']);cx=K(R(p['nx']))/(h*h);cy=K(R(p['ny']))/(h**3)
G=matrix(ZZ,d['generic_height_gram']);w=vector(ZZ,[0]*14+[1,-1,0]);zstar=QQ(protocol['carrier_label'])
counts={};root_checks=0;files=[]
for i,row in enumerate(d['basis_weierstrass_coordinates']):
    x,y=map(dec,row);assert y*y+a1*x*y+a3*y==x**3+a2*x*x+a4*x+a6
    X=x+b2/12;Y=y+(a1*x+a3)/2
    for sign in [-1,1]:
        path=OUT/('source-%02d-%s.json'%(i,'minus' if sign<0 else 'plus'));saved=read(path);files.append(path)
        slope=(sign*Y+cy)/(X-cx);z=(h*slope+shift)/(h*h)
        assert z==dec(saved['z_map']);degree=max(z.numerator().degree(),z.denominator().degree())
        assert degree==G[i,i]-sign*(G*w)[i]==saved['degree']
        status=saved['status'];counts[status]=counts.get(status,0)+1
        if degree<=1:
            assert status=='INHERITED_ALTERNATE_SECTION' if degree==1 else z!=zstar
            continue
        f=R(z.numerator()-zstar*z.denominator());f*=lcm([q.denominator() for q in f.list()]);f/=gcd([ZZ(q) for q in f.list()]);f=R(f)
        if f.leading_coefficient()<0:f=-f
        assert f==R(saved['primitive_polynomial']) and f.degree()==degree and f[0]
        proof=saved['root_proof'];prime=ZZ(proof['prime']);assert prime in protocol['limits']['prime_pool'] and prime.is_prime()
        Fp=PolynomialRing(Zmod(prime),'u');fp=Fp(f)
        assert fp.degree()==degree
        roots=[int(a) for a in range(int(prime)) if fp(a)==0]
        assert roots==proof['all_roots_mod_p'] and all(fp.derivative()(a)!=0 for a in roots)
        A=abs(ZZ(f[0]));B=abs(ZZ(f.leading_coefficient()));H=max(A,B)
        assert [str(A),str(B),str(H)]==[proof[k] for k in ['numerator_bound','denominator_bound','common_bound']]
        assert [r['initial_root'] for r in proof['lifts']]==roots
        for lift in proof['lifts']:
            M=ZZ(lift['modulus']);r=ZZ(lift['lift']);assert M==prime**M.valuation(prime) and M>8*H*H
            assert r%prime==lift['initial_root'] and PolynomialRing(Zmod(M),'u')(f)(r)==0
            a,b=[vector(ZZ,v) for v in lift['reduced_basis']]
            assert abs(matrix(ZZ,[a,b]).det())==M and all((v[0]-r*v[1])%M==0 for v in [a,b])
            assert a*a<=b*b and 2*abs(a*b)<=a*a
            if lift['decision']=='NO_BOUNDED_NUMERATOR_DENOMINATOR_VECTOR':assert a*a>2*H*H
            else:
                assert lift['decision']=='ONLY_POSSIBLE_RATIO_IS_NOT_A_ROOT' and a[1]
                candidate=QQ(a[0])/a[1];assert candidate==QQ(lift['candidate']) and f(candidate)!=0
            root_checks+=1
        assert status=='NO_RATIONAL_SOURCE_POINT' and not saved['rational_preimages'] and not saved['infinity_preimage']
assert counts=={'NO_RATIONAL_SOURCE_POINT':32,'INHERITED_ALTERNATE_SECTION':2}
result={'status':'PASS_INDEPENDENT_INFINITE_TRANSLATION_ORBIT_OBSTRUCTION',
    'classification':'new verified application of fibre invariance and existing Picard-image obstruction',
    'counts':counts,'Hensel_Gauss_certificates':root_checks,
    'inputs':{str(path.relative_to(ROOT)):sha(path) for path in [PARENT,PENCIL,OUT/'protocol.json',*files,Path(__file__)]},
    'deduction':'The first seed and its hyperelliptic conjugate cannot lie on the orbit of any signed displayed generic-section curve or O under the entire group of generic alternate-fibration translations and inversion.',
    'boundary':'Not all generic MW section curves or alternating-fibration automorphisms are covered. No other representative of the original residual direction is excluded.'}
path=OUT/'independent-replay.json'
if path.exists():assert read(path)==result
else:
    with path.open('x') as stream:json.dump(result,stream,indent=2,sort_keys=True);stream.write('\n')
print(result['status'],'root certificates',root_checks,flush=True)
